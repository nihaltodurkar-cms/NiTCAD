// Energy-balance transport, solved (Unit 52; DECISIONS.md ED1-ED13): equilibrium (EC5), the
// uniform-field limit at 1e3 and 1e4 V/cm (EC3), the drift-diffusion limit (EC4), current
// conservation (EC8) and velocity overshoot in Aste and Vahldieck's ballistic diode (EC10).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

constexpr double q = base::q_C, k_over_q = base::k_B_J_per_K / base::q_C;

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

// A 1D device on [0, length] with doping N(x) (positive donors, negative acceptors) and ohmic
// contacts at both ends.
template <class Doping>
device::Device line(double length_cm, int nodes, Doping doping,
                    physics::Semiconductor material = physics::silicon()) {
    mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, length_cm, nodes));
    const auto n = static_cast<std::size_t>(nodes);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double N = doping(m.points()[i][0]);
        (N > 0.0 ? donors : acceptors)[i] = std::abs(N);
    }
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", std::move(material)}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"left", device::ContactKind::ohmic, {0}},
                                                  {"right", device::ContactKind::ohmic,
                                                   {static_cast<mesh::NodeId>(n - 1)}}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

solve::BiasOptions options(bool energy, bool field = false) {
    solve::BiasOptions o;
    o.models.energy_balance = energy;
    o.models.field_mobility = field;
    o.step_control.max_cuts = 20;
    return o;
}

// The right contact's bias at each point, the left at 0.
std::vector<std::vector<double>> ramp(std::initializer_list<double> v) {
    std::vector<std::vector<double>> p;
    for (const double x : v) p.push_back({0.0, x});
    return p;
}

results::Sweep sweep(const device::Device& d, std::initializer_list<double> biases,
                     const solve::BiasOptions& o) {
    const auto points = ramp(biases);
    auto s = solve::sweep_bias(d, points, o);
    if (!s) FAIL(s.error().message);
    if (s->stopped) FAIL(s->stopped->message);
    REQUIRE(s->points.size() == points.size());
    return std::move(*s);
}

// EC8 as ruled 2026-10-29: the terminal currents sum to zero within 1e-10 of the largest or within
// the state's current resolution, whichever is larger. The sum is the Newton residual's, so plain
// drift-diffusion meets no tighter bound (measured on EC4's diode at -0.4 V: 1.2e-10 A/cm2 of
// 4.4e-4, against a resolution of 2.3e-7).
void conserved(const results::BiasPoint& p) {
    CAPTURE(p.terminal_current, p.terminal_current_resolution, p.bias_V);
    double sum = 0.0, largest = 0.0, resolution = 0.0;
    for (std::size_t c = 0; c < p.terminal_current.size(); ++c) {
        sum += p.terminal_current[c];
        largest = std::max(largest, std::abs(p.terminal_current[c]));
        resolution = std::max(resolution, p.terminal_current_resolution[c]);
    }
    REQUIRE(std::abs(sum) <= std::max(1e-10 * largest, resolution));
}

}  // namespace

TEST_CASE("energy balance: equilibrium is drift-diffusion's with the lattice temperature (EC5)") {
    const device::Device d = line(1e-4, 101, [](double x) { return x < 0.5e-4 ? -1e17 : 2e17; });
    const results::Sweep eb = sweep(d, {0.0}, options(true));
    const results::Sweep dd = sweep(d, {0.0}, options(false));
    const results::NodeFields& f = eb.points[0].fields;
    const results::NodeFields& g = dd.points[0].fields;
    for (std::size_t i = 0; i < f.n_cm3.size(); ++i) {
        CAPTURE(i);
        REQUIRE(std::abs(f.electron_temperature_K[i] - 300.0) <= 1e-12 * 300.0);
        REQUIRE(std::abs(f.hole_temperature_K[i] - 300.0) <= 1e-12 * 300.0);
        REQUIRE(std::abs(f.potential_V[i] - g.potential_V[i]) <=
                1e-12 * std::max(1.0, std::abs(g.potential_V[i])));
        REQUIRE(std::abs(f.n_cm3[i] - g.n_cm3[i]) <= 1e-12 * g.n_cm3[i]);
        REQUIRE(std::abs(f.p_cm3[i] - g.p_cm3[i]) <= 1e-12 * g.p_cm3[i]);
    }
    REQUIRE(g.electron_temperature_K.empty());
}

TEST_CASE("energy balance: the steady energy budget closes (EC7)") {
    // EC7 as ruled 2026-10-29: the delivered power equals the relaxation power, the recombined
    // pairs' gap and thermal energy and the contact terms (energy flux into the contacts, the
    // contact edges' half Joule heat, minus the Peltier term), within 1e-9 of the largest term.
    // A forward-biased pn diode (p 1e17 below 0.5 um, n 2e17 above, 1 um) in 1D, 2D and 3D; the
    // state is rebuilt from the solution's fields.
    for (const int D : {1, 2, 3}) {
        CAPTURE(D);
        const auto x = uniform(0.0, 1e-4, D == 3 ? 41 : 81);
        const auto y = uniform(0.0, 0.2e-4, 4);
        mesh::Mesh m = D == 1   ? *mesh::make_tensor_grid(x)
                       : D == 2 ? *mesh::make_tensor_grid(x, y)
                                : *mesh::make_tensor_grid(x, y, y);
        const std::size_t n = m.node_count();
        std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (m.points()[i][0] < 0.5e-4) {
                acceptors[i] = 1e17;
            } else {
                donors[i] = 2e17;
            }
        }
        const auto left = m.find_boundary("x_min")->nodes, right = m.find_boundary("x_max")->nodes;
        auto made = device::Device::create({.mesh = std::move(m),
                                            .temperature_K = 300.0,
                                            .regions = {{"silicon", physics::silicon()}},
                                            .node_region = std::vector<device::RegionId>(n, 0),
                                            .donors = std::move(donors),
                                            .acceptors = std::move(acceptors),
                                            .contacts = {{"anode", device::ContactKind::ohmic, left},
                                                         {"cathode", device::ContactKind::ohmic, right}}});
        REQUIRE(made.has_value());
        const device::Device& d = *made;
        const results::Sweep s = sweep(d, {-0.3, -0.6, -0.8}, options(true));
        const auto scaling = *assemble::make_scaling(d);
        auto dd = assemble::DriftDiffusion::create(d, scaling, {.energy_balance = true});
        REQUIRE(dd.has_value());
        for (const results::BiasPoint& b : s.points) {
            REQUIRE(dd->set_bias(b.bias_V).has_value());
            std::vector<double> state(dd->unknowns());
            const results::NodeFields& f = b.fields;
            for (std::size_t i = 0; i < n; ++i) {
                state[5 * i] = f.potential_V[i] / scaling.V_T;
                state[5 * i + 1] = f.n_cm3[i] / scaling.Ns;
                state[5 * i + 2] = f.p_cm3[i] / scaling.Ns;
                state[5 * i + 3] = (f.electron_temperature_K[i] - 300.0) / 300.0;
                state[5 * i + 4] = (f.hole_temperature_K[i] - 300.0) / 300.0;
            }
            const auto e = dd->energy_powers(state);
            const double out = e.relaxation + e.recombination + e.contact_flux + e.contact_joule -
                               e.peltier;
            const double scale = std::max({std::abs(e.input), std::abs(e.relaxation),
                                           std::abs(e.recombination), std::abs(e.contact_flux),
                                           std::abs(e.contact_joule), std::abs(e.peltier)});
            std::printf("EC7 %dD %.1f V: input %.6e relax %.6e rec %.6e cflux %.6e cjoule %.6e "
                        "peltier %.6e residual %.3e; miss %.3e, miss with residual %.3e\n",
                        D, b.bias_V[1], e.input, e.relaxation, e.recombination, e.contact_flux,
                        e.contact_joule, e.peltier, e.residual, (e.input - out) / scale,
                        (e.input - out - e.residual) / scale);
            // The residual's share: at most the power a continuity imbalance of the state's current
            // resolution carries at the largest band-edge energy (both scaled, k T_L units).
            double edge = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                edge = std::max({edge, std::abs(b.bands.conduction_eV[i]),
                                 std::abs(b.bands.valence_eV[i])});
            }
            const double resolution = b.terminal_current_resolution[0] /
                                      (scaling.J0 * std::pow(scaling.L_D, D - 1));
            const double allowed = std::max(1e-9 * scale, resolution * edge / scaling.V_T);
            std::printf("EC7 %dD %.1f V: residual %.3e of allowed %.3e\n", D, b.bias_V[1],
                        e.residual, allowed);
            CAPTURE(b.bias_V, e.input, out, e.residual, scale, allowed);
            REQUIRE(std::abs(e.input - out - e.residual) <= 1e-9 * scale);
            REQUIRE(std::abs(e.residual) <= allowed);
        }
    }
}

TEST_CASE("energy balance: a uniform field gives Canali's velocity and the local temperature (EC3)") {
    // Resistors 3 um long (86 electron and 75 hole relaxation lengths v_sat tau_w), 10 nm nodes;
    // the middle third is uniform at 1e3 and 1e4 V/cm. At 1e5 V/cm it is not (ruled 2026-10-29):
    // in saturation the cold carriers entering at the contact are not screened within any
    // practical length (1e17: E from 2.3e4 to 1.6e5 V/cm across 3 um, 6% across the middle third
    // of 30 um), so 1e5 V/cm is the discrete bulk check of tests/assemble/energy_balance_test.cpp.
    // mu0 is drift-diffusion's own (constant mobility) at the same bias, so the comparison does
    // not restate the doping mobility models.
    const double L = 3e-4;
    for (const double N : {1e16, 1e17, 1e18}) {
        for (const bool electrons : {true, false}) {
            const device::Device d =
                line(L, 301, [&](double) { return electrons ? N : -N; });
            const physics::CanaliParameters sat = physics::saturation(
                physics::silicon(), electrons ? physics::Carrier::electron
                                              : physics::Carrier::hole);
            const double tau_w = electrons ? physics::silicon_energy_relaxation.tau_n_s
                                           : physics::silicon_energy_relaxation.tau_p_s;
            // 1e3 and 1e4 V/cm, with the points between as continuation.
            const auto biases = {0.1, 0.3, 1.0, 3.0};
            const results::Sweep eb = sweep(d, biases, options(true));
            const results::Sweep dd = sweep(d, biases, options(false));
            for (const std::size_t k : {std::size_t{1}, std::size_t{3}}) {
                const results::BiasPoint& p = eb.points[k];
                const results::BiasPoint& r = dd.points[k];
                conserved(p);
                const double J = std::abs(p.terminal_current[0]), Jr = std::abs(r.terminal_current[0]);
                double worst_v = 0.0, worst_T = 0.0;
                for (std::size_t i = 100; i < 200; ++i) {
                    const auto& f = p.fields;
                    const double E = std::abs(f.potential_V[i + 1] - f.potential_V[i]) / 1e-6;  // 10 nm nodes
                    const double c = electrons ? 0.5 * (f.n_cm3[i] + f.n_cm3[i + 1])
                                               : 0.5 * (f.p_cm3[i] + f.p_cm3[i + 1]);
                    const double v = J / (q * c);
                    const auto& g = r.fields;
                    const double Er = std::abs(g.potential_V[i + 1] - g.potential_V[i]) / 1e-6;
                    const double cr = electrons ? 0.5 * (g.n_cm3[i] + g.n_cm3[i + 1])
                                                : 0.5 * (g.p_cm3[i] + g.p_cm3[i + 1]);
                    const double mu0 = Jr / (q * cr * Er);
                    const double expected = physics::canali_mobility(mu0, E, sat).mobility * E;
                    const auto& T = electrons ? f.electron_temperature_K : f.hole_temperature_K;
                    const double Te = 0.5 * (T[i] + T[i + 1]);
                    const double T_expected = 300.0 + 2.0 / 3.0 * tau_w * v * E / k_over_q;
                    worst_v = std::max(worst_v, std::abs(v - expected) / expected);
                    worst_T = std::max(worst_T, std::abs(Te - T_expected) / T_expected);
                }
                UNSCOPED_INFO("N " << N << (electrons ? " n" : " p") << " point " << k
                                   << ": velocity " << worst_v << ", temperature " << worst_T);
                CAPTURE(N, electrons, k, worst_v, worst_T);
                REQUIRE(worst_v <= 1e-3);
                REQUIRE(worst_T <= 1e-3);
            }
        }
    }
}

TEST_CASE("energy balance: the drift-diffusion limit (EC4, tau_w -> 0 and v_sat -> infinity)") {
    // EC4 as ruled 2026-10-29: with tau_w -> 0 alone the carriers stay in the local balance, so
    // the mobility tends to Canali's mu(E), not mu0; constant-mobility drift-diffusion is the limit
    // with v_sat -> infinity too. The tau_w-only comparison is the next test case, printed only.
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.energy_relaxation = {.tau_n_s = 1e-18, .tau_p_s = 1e-18};
    p.electron_saturation.v_sat_cm_s = 1e30;
    p.hole_saturation.v_sat_cm_s = 1e30;
    const auto doping = [](double x) { return x < 0.5e-4 ? -1e17 : 2e17; };
    const device::Device eb_device = line(1e-4, 101, doping, *physics::Semiconductor::create(p));
    const device::Device dd_device = line(1e-4, 101, doping);
    // Forward (the n side negative): the reverse currents lie below the resolution (EC8's note).
    const auto biases = {-0.2, -0.4, -0.6, -0.7, -0.8};
    const results::Sweep eb = sweep(eb_device, biases, options(true));
    const results::Sweep dd = sweep(dd_device, biases, options(false));
    for (std::size_t k = 0; k < eb.points.size(); ++k) {
        conserved(eb.points[k]);
        const double a = eb.points[k].terminal_current[1], b = dd.points[k].terminal_current[1];
        UNSCOPED_INFO("point " << k << ": relative difference " << std::abs(a - b) / std::abs(b));
        CAPTURE(k, a, b);
        REQUIRE(std::abs(a - b) <= 1e-6 * std::abs(b));
    }
}

TEST_CASE("energy balance: tau_w -> 0 alone, against drift-diffusion with Canali (EC4, measured)") {
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.energy_relaxation = {.tau_n_s = 1e-18, .tau_p_s = 1e-18};
    const auto doping = [](double x) { return x < 0.5e-4 ? -1e17 : 2e17; };
    const device::Device eb_device = line(1e-4, 101, doping, *physics::Semiconductor::create(p));
    const device::Device dd_device = line(1e-4, 101, doping);
    const auto biases = {-1.0, 0.2, 0.4, 0.6, 0.7, 1.0};
    const results::Sweep eb = sweep(eb_device, biases, options(true));
    const results::Sweep dd = sweep(dd_device, biases, options(false, true));
    for (std::size_t k = 0; k < eb.points.size(); ++k) {
        const double a = eb.points[k].terminal_current[1], b = dd.points[k].terminal_current[1];
        std::printf("EC4 tau_w -> 0 vs DD + Canali, point %zu: %.3e (EB %.6e, DD %.6e A/cm2)\n", k,
                    std::abs(a - b) / std::abs(b), a, b);
    }
}

TEST_CASE("energy balance: velocity overshoot in the ballistic diode (EC10)") {
    // Aste and Vahldieck, arXiv physics/0312021: n+ 1e18 / n 2e15 / n+ 1e18, 0.1 um each, 1 V;
    // mu0(N) = 80 + 1350 / (1 + (N / 1.12e17)^0.72), v_s = 1.03e7 cm/s, beta = 2. This unit's
    // constant tau_w (silicon's), c = 5/2. 1 nm nodes.
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.electron_mobility.mu_min = 80.0;
    p.electron_mobility.mu_max = 1430.0;
    p.electron_mobility.N_ref = 1.12e17;
    p.electron_mobility.alpha = 0.72;
    p.electron_saturation = {.v_sat_cm_s = 1.03e7, .beta = 2.0};
    p.energy_relaxation = physics::silicon_energy_relaxation;
    const physics::Semiconductor si = *physics::Semiconductor::create(p);
    const auto doping = [](double x) { return x < 0.1e-4 || x > 0.2e-4 ? 1e18 : 2e15; };
    const device::Device d = line(0.3e-4, 301, doping, si);
    const auto biases = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const results::Sweep eb = sweep(d, biases, options(true));
    const results::Sweep dd = sweep(d, biases, options(false, true));
    const auto peak = [](const results::BiasPoint& b) {
        const double J = std::abs(b.terminal_current[0]);
        double best = 0.0;
        std::size_t at = 0;
        for (std::size_t i = 0; i < b.fields.n_cm3.size(); ++i) {
            const double v = J / (q * b.fields.n_cm3[i]);
            if (v > best) {
                best = v;
                at = i;
            }
        }
        return std::pair{best, at};
    };
    const auto [v_eb, at_eb] = peak(eb.points.back());
    const auto [v_dd, at_dd] = peak(dd.points.back());
    const double T_max = *std::ranges::max_element(eb.points.back().fields.electron_temperature_K);
    std::printf("EC10: peak velocity %.4f v_s at x = %.0f nm (Aste-Vahldieck EBM 1.87 v_s); "
                "DD + Canali %.6f v_s at %.0f nm; peak T_n %.0f K; J %.4e A/cm2 (DD %.4e)\n",
                v_eb / 1.03e7, static_cast<double>(at_eb), v_dd / 1.03e7,
                static_cast<double>(at_dd), T_max, eb.points.back().terminal_current[1],
                dd.points.back().terminal_current[1]);
    conserved(eb.points.back());
    REQUIRE(v_eb > 1.3 * 1.03e7);
    REQUIRE(at_eb >= 170);  // in the n region, within 30 nm of the drain junction at 200 nm
    REQUIRE(at_eb <= 200);
    // Drift-diffusion's clause as ruled 2026-10-29: its J / (q n) carries diffusion as well as
    // drift (1.34 v_s at the drain junction, printed above), so the bound is on its drift velocity
    // mu_C(mu0_e, E_e) E_e on every edge, and energy balance's peak exceeds its J / (q n) peak.
    const results::BiasPoint& r = dd.points.back();
    const physics::CanaliParameters& sat = physics::saturation(si, physics::Carrier::electron);
    for (std::size_t i = 0; i + 1 < r.fields.potential_V.size(); ++i) {
        const auto mu = [&](std::size_t k) {
            return physics::caughey_thomas_mobility(si, physics::Carrier::electron,
                                                    doping(1e-7 * static_cast<double>(k)), 300.0);
        };
        const double mu0 = 2.0 * mu(i) * mu(i + 1) / (mu(i) + mu(i + 1));
        const double E = std::abs(r.fields.potential_V[i + 1] - r.fields.potential_V[i]) / 1e-7;
        CAPTURE(i, E);
        REQUIRE(physics::canali_mobility(mu0, E, sat).mobility * E <= 1.03e7 * (1.0 + 1e-12));
    }
    REQUIRE(v_eb > v_dd);
}

namespace {

// The legacy MOSFET of electrothermal_test.cpp (n+ source and drain, N_A 1e17, 5 nm oxide, n+
// polysilicon gate), gate length Lg, 5 nm nodes at the junctions. Contacts: source, drain, body,
// gate.
device::Device mosfet(double Lg, physics::Semiconductor material = physics::silicon()) {
    constexpr double Lsd = 3e-5, depth = 2e-5;
    const double L = 2 * Lsd + Lg;
    const auto x = legacy_graded_mesh(L, {Lsd, Lsd + Lg}, 5e-7, 5e-6, 1.15);
    const auto y = legacy_graded_mesh(depth, 0.0, depth / 400.0, depth / 20.0, 1.15);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    const auto sd = [](double px, double py, double edge, bool source) {
        const double s = source ? edge - px : px - edge;
        return 1e19 * std::exp(-(py * py) / (2.0 * 5e-6 * 5e-6)) * 0.5 *
               std::erfc(-s / (std::sqrt(2.0) * 1e-6));
    };
    std::vector<double> donors(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        donors[i] = sd(p[0], p[1], Lsd, true) + sd(p[0], p[1], Lsd + Lg, false);
    }
    std::vector<mesh::NodeId> source, drain, gate;
    for (const mesh::NodeId v : m.find_boundary("y_min")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        (px <= Lsd ? source : px >= Lsd + Lg ? drain : gate).push_back(v);
    }
    auto body = m.find_boundary("y_max")->nodes;
    device::Contact g{"gate", device::ContactKind::gate, std::move(gate)};
    g.gate = {.boundary = "y_min", .oxide_thickness_cm = 5e-7};
    auto dev = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", std::move(material)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::vector<double>(n, 1e17),
         .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                      {"drain", device::ContactKind::ohmic, std::move(drain)},
                      {"body", device::ContactKind::ohmic, std::move(body)},
                      std::move(g)}});
    if (!dev) FAIL(dev.error().message);
    return std::move(*dev);
}

// The gate up to 1.5 V, then the drain up to 1.5 V.
std::vector<std::vector<double>> mosfet_points() {
    std::vector<std::vector<double>> points;
    for (int k = 1; k <= 6; ++k) points.push_back({0.0, 0.0, 0.0, 0.25 * k});
    for (int k = 1; k <= 15; ++k) points.push_back({0.0, 0.1 * k, 0.0, 1.5});
    return points;
}

results::Sweep mosfet_sweep(const device::Device& d, const solve::BiasOptions& o) {
    const auto t0 = std::chrono::steady_clock::now();
    auto s = solve::sweep_bias(d, mosfet_points(), o);
    if (!s) FAIL(s.error().message);
    if (s->stopped) FAIL(s->stopped->message);
    std::size_t iterations = 0, cuts = 0;
    for (const results::BiasPoint& p : s->points) {
        iterations += p.convergence.iterations.size();
        cuts += static_cast<std::size_t>(p.bias_substeps);
    }
    const std::chrono::duration<double> t = std::chrono::steady_clock::now() - t0;
    std::printf("  MOSFET sweep (%zu nodes, %s): %.1f s, %zu Newton iterations, %zu substeps\n",
                d.mesh().node_count(), o.models.energy_balance ? "energy balance" : "drift-diffusion",
                t.count(), iterations, cuts);
    std::fflush(stdout);
    return std::move(*s);
}

}  // namespace

TEST_CASE("energy balance: the MOSFET's drift-diffusion limit (EC4)", "[.mosfet]") {
    // EC4's MOSFET half (the double limit as ruled 2026-10-29) at L = 0.6 um, over the sweep. With
    // v_sat = 1e30 the channel's drift is unclamped and its first-order tau_w heating is about
    // 4e-4 of T_L at 1e-18: measured 1.23e-6 there, 1.23e-8 at 1e-20, so tau_w is 1e-20 here
    // (ruled 2026-10-29).
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.energy_relaxation = {.tau_n_s = 1e-20, .tau_p_s = 1e-20};
    p.electron_saturation.v_sat_cm_s = 1e30;
    p.hole_saturation.v_sat_cm_s = 1e30;
    const results::Sweep eb = mosfet_sweep(mosfet(6e-5, *physics::Semiconductor::create(p)),
                                           options(true));
    const results::Sweep dd = mosfet_sweep(mosfet(6e-5), options(false));
    double worst = 0.0;
    std::size_t skipped = 0;
    for (std::size_t k = 0; k < eb.points.size(); ++k) {
        conserved(eb.points[k]);
        for (const std::size_t c : {std::size_t{0}, std::size_t{1}}) {  // source, drain
            const double a = eb.points[k].terminal_current[c], b = dd.points[k].terminal_current[c];
            // A current below drift-diffusion's resolution (the gate ramp's) is not resolved and
            // is skipped (ruled 2026-10-29).
            if (std::abs(b) <= dd.points[k].terminal_current_resolution[c]) {
                ++skipped;
                continue;
            }
            CAPTURE(k, c, a, b);
            worst = std::max(worst, std::abs(a - b) / std::abs(b));
            REQUIRE(std::abs(a - b) <= 1e-6 * std::abs(b));
        }
    }
    std::printf("EC4 MOSFET: worst relative difference %.3e, %zu of %zu currents unresolved (skipped)\n",
                worst, skipped, 2 * eb.points.size());
}

TEST_CASE("energy balance: MOSFET channel length (EC11)", "[.mosfet]") {
    // At L = 2 um and 0.1 um, Vg = 1.5 V, Vd = 1.5 V (saturation): at 2 um energy balance and
    // drift-diffusion with Canali agree within 5%; at 0.1 um the energy-balance current is larger.
    // Both reported.
    double ratio[2] = {};
    for (const double Lg : {2e-4, 1e-5}) {
        CAPTURE(Lg);
        const device::Device d = mosfet(Lg);
        const results::Sweep eb = mosfet_sweep(d, options(true));
        const results::Sweep dd = mosfet_sweep(d, options(false, true));
        conserved(eb.points.back());
        const double a = eb.points.back().terminal_current[1];
        const double b = dd.points.back().terminal_current[1];
        ratio[Lg > 1e-4 ? 0 : 1] = a / b;
        std::printf("EC11: L = %.1f um: I_d energy balance %.6e, drift-diffusion + Canali %.6e "
                    "A/cm (ratio %.4f)\n",
                    Lg * 1e4, a, b, a / b);
    }
    REQUIRE(std::abs(ratio[0] - 1.0) <= 0.05);
    REQUIRE(ratio[1] > 1.0);
}
