// Small-signal (AC) runs (ARCHITECTURE.md section 11, Unit 22): input errors; conservation and
// gauge invariance of the admittance; reciprocity in equilibrium; the low-frequency limits (the
// DC conductance and the quasi-static capacitance); the analytic RC and dielectric-relaxation
// responses of a MOS capacitor (meshed oxide and lumped gate); the small-signal fields; 2D and 3D
// extrusions; cancellation, progress and the run record; the MOS capacitor in inversion, a gated
// diode, the conductance method of interface traps, junction and long-diode admittances. The
// comparison with transient sine runs (11 s in Release) carries the hidden tag [.transient] and
// runs in the Release configuration only (ctest solve_transient).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <initializer_list>
#include <numbers>
#include <stop_token>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/cv.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "NiTCAD/solve/waveform.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;
using base::ErrorCode;
using Complex = std::complex<double>;

namespace {

constexpr double two_pi = 2.0 * std::numbers::pi;

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

// The Unit 9 diode: legacy graded mesh, 1e17 / 1e17, junction at 1 um; anode (contact 0) on x_min.
device::Device diode(mesh::Mesh m) {
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

std::vector<double> diode_axis() { return legacy_graded_mesh(2e-4, {1e-4}, 1e-8, 1e-6); }

// The Unit 21 1D p-type MOS capacitor: a metal gate (4.1 eV) on t_ox of oxide over L of silicon
// (N_A), the substrate contact on x_max; meshed (oxide region, electrode on x_min) or lumped
// (Unit 12 gate). The silicon mesh is graded from 1 nm at the surface.
device::Device moscap(double t_ox, double L, double N_A, bool meshed,
                      std::vector<device::Interface> interfaces = {},
                      physics::SemiconductorParameters silicon = physics::silicon_parameters) {
    std::vector<double> si = legacy_graded_mesh(L, {0.0}, 1e-7, L / 200.0);
    std::vector<double> x;
    if (meshed) {
        const double h = si[1];
        for (int k = 0; k <= 20; ++k) x.push_back(-t_ox + (t_ox - 0.5 * h) * k / 20);
        for (const double v : si) x.push_back(v + 0.5 * h);
    } else {
        x = si;
    }
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> acceptors(n, N_A);
    if (meshed) {
        for (std::size_t i = 0; i < n; ++i) {
            if (m.points()[i][0] < 0.0) {
                region[i] = 1;
                acceptors[i] = 0.0;
            }
        }
    }
    auto top = m.find_boundary("x_min")->nodes;
    auto substrate = m.find_boundary("x_max")->nodes;
    const auto kind = meshed ? device::ContactKind::electrode : device::ContactKind::gate;
    device::Contact gate{"gate", kind, std::move(top)};
    if (meshed) {
        gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.1};
    } else {
        gate.gate = {.boundary = "x_min",
                     .oxide_thickness_cm = t_ox,
                     .electrode = device::GateElectrode::metal,
                     .work_function_eV = 4.1};
    }
    std::vector<device::Region> regions{{"silicon", *physics::Semiconductor::create(silicon)}};
    if (meshed) regions.push_back({"oxide", physics::silicon_dioxide()});
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = std::move(regions),
                                    .node_region = std::move(region),
                                    .donors = std::vector<double>(n, 0.0),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {std::move(gate),
                                                 {"substrate", device::ContactKind::ohmic,
                                                  std::move(substrate)}},
                                    .interfaces = std::move(interfaces)});
}

// The Unit 21 assembler's 2D MOS structure: oxide on [-5 nm, 0) with the electrode on part of
// x_min, p-silicon (1e17) on (0, 120 nm] with an n+ pocket (1e19) and the substrate on x_max, a
// source contact on the pocket's y_max face (its first node on an interface edge); a fixed
// charge, two trap levels, a trap band and surface recombination at the interface.
device::Device mos2d() {
    const auto x = *mesh::straddle_interface(uniform(-5e-7, 1.2e-5, 30), 0.0, 1e-7);
    const auto t = uniform(0.0, 3e-6, 4);
    mesh::Mesh m = *mesh::make_tensor_grid(x, t);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        const bool oxide = p[0] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (oxide) continue;
        acceptors[i] = 1e17;
        if (p[1] > 2e-6 && p[0] < 3e-6) donors[i] = 1e19;
    }
    std::vector<mesh::NodeId> electrode, source;
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        if (m.points()[static_cast<std::size_t>(v)][1] < 1.5e-6) electrode.push_back(v);
    }
    for (const mesh::NodeId v : m.find_boundary("y_max")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        if (px > 0.0 && px < 3e-6) source.push_back(v);
    }
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = 3e11;
    f.traps.levels = {{.type = physics::TrapType::donor, .density_cm2 = 2e11, .energy_eV = -0.2,
                       .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 3e-16},
                      {.type = physics::TrapType::acceptor, .density_cm2 = 1e11,
                       .energy_eV = 0.25, .sigma_n_cm2 = 2e-16, .sigma_p_cm2 = 1e-15}};
    f.traps.bands = {{.type = physics::TrapType::acceptor, .density_cm2_eV = 5e11,
                      .energy_low_eV = -0.1, .energy_high_eV = 0.4}};
    f.recombination_velocity_n_cm_s = 1e3;
    f.recombination_velocity_p_cm_s = 4e2;
    device::Contact gate{"gate", device::ContactKind::electrode, std::move(electrode)};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)},
                      {"source", device::ContactKind::ohmic, std::move(source)}},
         .interfaces = {std::move(f)}});
}

std::vector<double> decades(double from, double to, int per_decade) {
    std::vector<double> f;
    const int n = static_cast<int>(std::lround(std::log10(to / from) * per_decade));
    for (int k = 0; k <= n; ++k) {
        f.push_back(from * std::pow(10.0, static_cast<double>(k) / per_decade));
    }
    return f;
}

results::SmallSignal run(const device::Device& d, std::vector<std::vector<double>> points,
                         std::vector<double> frequencies, solve::BiasOptions steady = {}) {
    const solve::SmallSignalOptions o{.steady = std::move(steady),
                                      .frequencies_Hz = std::move(frequencies)};
    auto r = solve::solve_small_signal(d, points, o);
    REQUIRE(r.has_value());
    if (r->stopped) {
        CAPTURE(r->stopped->message);
        FAIL("the small-signal run stopped");
    }
    return std::move(*r);
}

// The series circuit of the gate: C (the oxide and the semiconductor's surface layer, the
// quasi-static capacitance) in series with the substrate, R = L / (q mu_p N_A) shunted by its
// geometric capacitance C_g = eps_si / L (Unit 21's RC response, in the frequency domain).
Complex rc_admittance(double omega, double C, double R, double Cg) {
    const Complex j{0.0, 1.0};
    return 1.0 / (1.0 / (j * omega * C) + R / (1.0 + j * omega * R * Cg));
}

}  // namespace

TEST_CASE("small signal: input errors") {
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const std::vector<std::vector<double>> points{{0.3, 0.0}};
    const auto error = [&](solve::SmallSignalOptions o) {
        const auto r = solve::solve_small_signal(d, points, o);
        REQUIRE_FALSE(r.has_value());
        return r.error();
    };
    REQUIRE(error({}).code == ErrorCode::invalid_input);
    REQUIRE(error({.frequencies_Hz = {1.0, -1.0}}).context->index == 1);
    REQUIRE(error({.frequencies_Hz = {std::nan("")}}).code == ErrorCode::invalid_input);
    REQUIRE(error({.frequencies_Hz = {1.0, 2.0}, .field_frequencies_Hz = {3.0}}).code ==
            ErrorCode::invalid_input);
    const std::vector<std::vector<double>> bad{{0.3}};
    const auto r = solve::solve_small_signal(d, bad, {.frequencies_Hz = {1.0}});
    REQUIRE(r.error().code == ErrorCode::invalid_input);
}

TEST_CASE("small signal: MOS capacitor RC response, meshed oxide and lumped gate") {
    // 10 nm of oxide over 100 um of 1e16 p-silicon, accumulated at -2 V (Unit 21's fixture):
    // Y = 1 / (1 / (i omega C) + R / (1 + i omega R C_g)), tau = R (C + C_g) = 5.0 ns, from
    // omega tau = 1e-4 to 100: measured 6.7e-4 (meshed and lumped; Unit 21 measured the step
    // response 0.09% short, the accumulation layer conducting better than the bulk). At
    // omega tau = 1e-4 Im Y / omega is the quasi-static capacitance to 6e-9.
    for (const bool meshed : {true, false}) {
        const double L = 1e-2, N_A = 1e16;
        const device::Device d = moscap(1e-6, L, N_A, meshed);
        const double V0 = -2.0, dV = 1e-3;
        const auto q_at = [&](double V) {
            return solve::solve_bias(d, std::vector<double>{V, 0.0})->gate_charge[0];
        };
        const double C = (q_at(V0 + dV) - q_at(V0 - dV)) / (2.0 * dV);
        const physics::Semiconductor si = physics::silicon();
        const double mu = physics::caughey_thomas_mobility(si, physics::Carrier::hole, N_A, 300.0);
        const double R = L / (base::q_C * mu * N_A);
        const double Cg = si.parameters().eps_r * base::eps0_F_per_cm / L;
        const double tau = R * (C + Cg);
        std::vector<double> f;
        for (const double wt : {1e-4, 1e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0}) {
            f.push_back(wt / (two_pi * tau));
        }
        const auto r = run(d, {{V0, 0.0}}, f);
        double worst = 0.0;
        for (std::size_t k = 0; k < f.size(); ++k) {
            const Complex y = r.admittance(0, k, 0, 0);
            const Complex expected = rc_admittance(two_pi * f[k], C, R, Cg);
            worst = std::max(worst, std::abs(y / expected - 1.0));
            CAPTURE(meshed, f[k] * two_pi * tau, y, expected);
        }
        CAPTURE(meshed, worst, r.capacitance(0, 0, 0, 0) / C - 1.0);
        REQUIRE(worst < 3e-3);
        // The low-frequency limit is the quasi-static capacitance.
        REQUIRE(std::abs(r.capacitance(0, 0, 0, 0) / C - 1.0) < 1e-5);
    }
}

TEST_CASE("small signal: dielectric relaxation behind a thick oxide") {
    // 100 um of oxide over 10 um of 1e16 p-silicon (Unit 21's fixture): the corner at
    // omega = 1 / (R (C + C_g)), 1.03 times the dielectric relaxation rate: measured 1.3e-4.
    const double L = 1e-3, N_A = 1e16;
    const device::Device d = moscap(1e-2, L, N_A, true);
    const double V0 = -1.0, dV = 1e-3;
    const auto q_at = [&](double V) {
        return solve::solve_bias(d, std::vector<double>{V, 0.0})->gate_charge[0];
    };
    const double C = (q_at(V0 + dV) - q_at(V0 - dV)) / (2.0 * dV);
    const physics::Semiconductor si = physics::silicon();
    const double mu = physics::caughey_thomas_mobility(si, physics::Carrier::hole, N_A, 300.0);
    const double R = L / (base::q_C * mu * N_A);
    const double Cg = si.parameters().eps_r * base::eps0_F_per_cm / L;
    const double tau = R * (C + Cg);
    std::vector<double> f;
    for (const double wt : {0.01, 0.1, 1.0, 10.0}) f.push_back(wt / (two_pi * tau));
    const auto r = run(d, {{V0, 0.0}}, f);
    double worst = 0.0;
    for (std::size_t k = 0; k < f.size(); ++k) {
        const Complex expected = rc_admittance(two_pi * f[k], C, R, Cg);
        worst = std::max(worst, std::abs(r.admittance(0, k, 0, 0) / expected - 1.0));
    }
    CAPTURE(worst);
    REQUIRE(worst < 1e-3);
}

TEST_CASE("small signal: the admittance conserves current and is gauge invariant") {
    // The 2D MOS structure with traps, an electrode and two ohmic contacts (one on an interface
    // edge): every column (sum over the contacts' currents) sums to zero within its reported
    // resolution, and every row (a common bias shift) likewise, at every frequency. The DC
    // conductances here are some 1e-12 S/cm, below the resolution of the (psi, n, p) state (a
    // tenth of it to several times it), so the sums are measured against the resolution, not
    // against the largest entry; where the displacement current dominates they vanish to 3e-15 of
    // it.
    const device::Device d = mos2d();
    const std::vector<double> f{0.0, 1e2, 1e5, 1e8, 1e11, 1e13};
    for (const std::vector<double>& bias :
         {std::vector<double>{-0.4, 0.0, 0.1}, std::vector<double>{0.0, 0.0, 0.0}}) {
        const auto r = run(d, {bias}, f);
        double column = 0.0, row = 0.0, resolved = 0.0;
        for (std::size_t k = 0; k < f.size(); ++k) {
            const std::vector<double>& res = r.points[0].resolution[k];
            double largest = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    largest = std::max(largest, std::abs(r.admittance(0, k, i, j)));
                }
            }
            double res_all = 0.0;
            for (const double v : res) res_all = std::max(res_all, v);
            for (std::size_t j = 0; j < 3; ++j) {
                Complex sum_i, sum_j;
                for (std::size_t i = 0; i < 3; ++i) {
                    sum_i += r.admittance(0, k, i, j);
                    sum_j += r.admittance(0, k, j, i);
                }
                column = std::max(column, std::abs(sum_i) / res[j]);
                row = std::max(row, std::abs(sum_j) / (3.0 * res_all));
                if (f[k] >= 1e11) resolved = std::max(resolved, std::abs(sum_i) / largest);
            }
        }
        CAPTURE(bias, column, row, resolved);
        REQUIRE(column <= 1.0);  // measured at most 0.43
        REQUIRE(row <= 1.0);
        REQUIRE(resolved < 1e-13);
    }
}

TEST_CASE("small signal: in equilibrium the admittance is symmetric") {
    // Reciprocity of a passive device in equilibrium (Onsager): Y_ij = Y_ji within the
    // resolution of the two entries (measured at most 0.06 of it). Out of equilibrium it fails
    // by up to 114 times the resolution here: a biased device is not reciprocal.
    const device::Device d = mos2d();
    const std::vector<double> f{0.0, 1e2, 1e5, 1e8, 1e11, 1e13};
    double equilibrium = 0.0, biased = 0.0;
    for (const bool zero : {true, false}) {
        const auto r = run(d, {zero ? std::vector<double>{0.0, 0.0, 0.0}
                                    : std::vector<double>{-0.4, 0.0, 0.1}},
                           f);
        for (std::size_t k = 0; k < f.size(); ++k) {
            const std::vector<double>& res = r.points[0].resolution[k];
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < i; ++j) {
                    const double miss =
                        std::abs(r.admittance(0, k, i, j) - r.admittance(0, k, j, i)) /
                        (res[i] + res[j]);
                    (zero ? equilibrium : biased) = std::max(zero ? equilibrium : biased, miss);
                }
            }
        }
    }
    CAPTURE(equilibrium, biased);
    REQUIRE(equilibrium <= 1.0);
    REQUIRE(biased > 10.0);
}

TEST_CASE("small signal: at zero frequency the admittance is the DC conductance") {
    // The diode's Y(0) against the fourth-order central difference of the DC current
    // (dV = 3e-4 V): measured 4.1e-8 at 0.5 V and 2.6e-10 at 0.7 V. Lower down the DC current
    // itself is resolved only to about 1e-4 of its derivative (at 0.3 V), and at 0 V both the DC
    // difference and Y(0) lie below their resolutions.
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    for (const double V : {0.5, 0.7}) {
        const auto r = run(d, {{V, 0.0}}, {0.0, 1.0});
        const double h = 3e-4;
        const auto I = [&](double v) {
            return solve::solve_bias(d, std::vector<double>{v, 0.0})->terminal_current[0];
        };
        const double G =
            (8.0 * (I(V + h) - I(V - h)) - (I(V + 2.0 * h) - I(V - 2.0 * h))) / (12.0 * h);
        const Complex y0 = r.admittance(0, 0, 0, 0);
        CAPTURE(V, y0.real() / G - 1.0, r.admittance(0, 1, 0, 0).real() / G - 1.0);
        REQUIRE(y0.imag() == 0.0);
        REQUIRE(std::abs(y0.real() / G - 1.0) < 1e-6);
        REQUIRE(std::abs(r.admittance(0, 1, 0, 0).real() / G - 1.0) < 1e-6);
        REQUIRE(r.points[0].resolution[0][0] < 1e-6 * y0.real());
    }
    // At 0 V the conductance is unresolved: Y(0) lies within its resolution of zero.
    const auto r = run(d, {{0.0, 0.0}}, {0.0});
    REQUIRE(std::abs(r.admittance(0, 0, 0, 0)) <= r.points[0].resolution[0][0]);
}

TEST_CASE("small signal: the fields of a driven contact") {
    // The potential response is 1 V/V on the driven contact's nodes and 0 on the others'; the
    // densities do not move on an ohmic contact.
    const device::Device d = mos2d();
    const solve::SmallSignalOptions o{.frequencies_Hz = {1e3, 1e9}, .field_frequencies_Hz = {1e9}};
    const auto r = *solve::solve_small_signal(d, std::vector<std::vector<double>>{{-0.4, 0.0, 0.1}},
                                              o);
    REQUIRE(!r.stopped);
    REQUIRE(r.points[0].fields.size() == 3);
    for (const auto& f : r.points[0].fields) {
        REQUIRE(f.frequency_Hz == 1e9);
        for (std::size_t c = 0; c < 3; ++c) {
            for (const mesh::NodeId v : d.contacts()[c].nodes) {
                const auto i = static_cast<std::size_t>(v);
                REQUIRE(std::abs(f.potential[i] - (c == f.contact ? 1.0 : 0.0)) < 1e-12);
                if (d.contacts()[c].kind == device::ContactKind::ohmic) {
                    REQUIRE(f.n_cm3[i] == Complex{});
                    REQUIRE(f.p_cm3[i] == Complex{});
                }
            }
        }
    }
}

TEST_CASE("small signal: y-uniform 2D and 3D devices reproduce 1D") {
    const auto x = diode_axis();
    const auto y = uniform(0.0, 2e-5, 3);
    const std::vector<double> f{0.0, 1e6, 1e9, 1e11};
    const auto r1 = run(diode(*mesh::make_tensor_grid(x)), {{0.5, 0.0}}, f);
    const auto r2 = run(diode(*mesh::make_tensor_grid(x, y)), {{0.5, 0.0}}, f);
    const auto r3 = run(diode(*mesh::make_tensor_grid(x, y, y)), {{0.5, 0.0}}, f);
    // Per unit width (2e-5 cm) and area (4e-10 cm^2), within the resolutions of the two runs.
    double worst = 0.0;
    for (std::size_t k = 0; k < f.size(); ++k) {
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                const Complex y1 = r1.admittance(0, k, i, j);
                const double res1 = r1.points[0].resolution[k][j];
                const double d2 = std::abs(r2.admittance(0, k, i, j) / 2e-5 - y1) /
                                  (res1 + r2.points[0].resolution[k][j] / 2e-5);
                const double d3 = std::abs(r3.admittance(0, k, i, j) / 4e-10 - y1) /
                                  (res1 + r3.points[0].resolution[k][j] / 4e-10);
                worst = std::max({worst, d2, d3});
            }
        }
    }
    CAPTURE(worst);
    REQUIRE(worst <= 1.0);
}

TEST_CASE("small signal: cancellation keeps the completed points; progress is monotonic") {
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const std::vector<std::vector<double>> points{{0.1, 0.0}, {0.2, 0.0}, {0.3, 0.0}};
    const solve::SmallSignalOptions o{.frequencies_Hz = {1.0, 1e3, 1e6}};
    std::stop_source stop;
    std::vector<solve::Progress> events;
    solve::RunControl control{stop.get_token(), [&](const solve::Progress& p) {
                                  events.push_back(p);
                                  if (p.phase == solve::Phase::small_signal && p.point == 1 &&
                                      p.iteration == 2) {
                                      stop.request_stop();
                                  }
                              }};
    const auto r = solve::solve_small_signal(d, points, o, nullptr, control);
    REQUIRE(r.has_value());
    REQUIRE(r->stopped.has_value());
    REQUIRE(r->stopped->code == ErrorCode::cancelled);
    REQUIRE(r->points.size() == 1);
    REQUIRE(r->points[0].admittance.size() == 3);
    for (std::size_t e = 1; e < events.size(); ++e) {
        const auto& a = events[e - 1];
        const auto& b = events[e];
        const auto key = [](const solve::Progress& p) {
            return std::tuple{static_cast<int>(p.phase), p.point, p.iteration};
        };
        REQUIRE(key(a) < key(b));
    }
    REQUIRE(events.back().phase == solve::Phase::small_signal);
    REQUIRE(events.back().frequency_Hz == 1e3);
}

TEST_CASE("small signal: the run record") {
    const device::Device d = moscap(1e-6, 1e-3, 1e16, false);
    const std::vector<std::vector<double>> points{{-1.0, 0.0}};
    solve::SmallSignalOptions a{.frequencies_Hz = {1e3}};
    solve::SmallSignalOptions b = a;
    b.frequencies_Hz = {1e4};
    solve::SmallSignalOptions c = a;
    c.steady.equations = solve::Equations::equilibrium_poisson;
    solve::SmallSignalOptions e = c;
    e.steady.models.srh = false;  // read by the small-signal system even with a quasi-static point
    const auto id = [&](const solve::SmallSignalOptions& o) {
        return solve::make_run_record(d, o, points).input_identity;
    };
    REQUIRE(id(a) != id(b));
    REQUIRE(id(a) != id(c));
    REQUIRE(id(c) != id(e));
    REQUIRE(id(a) == solve::solve_small_signal(d, points, a)->run.input_identity);
    const auto settings = solve::make_run_record(d, c, points).settings;
    const std::pair<std::string, double> quasi_static{"small_signal.operating_point_equations", 1.0};
    REQUIRE(std::ranges::find(settings, quasi_static) != settings.end());
}


namespace {

// The high-frequency capacitance of a 1D lumped-gate MOS capacitor at an equilibrium state,
// computed here on the device's own mesh by the box method: the holes follow the potential
// (dp = p dpsi / V_T), the electrons keep their number with a uniform quasi-Fermi level
// (dn = n (dpsi - dphi) / V_T, sum V dn = 0), so the linearized Poisson equation
//     sum eps A / h (dpsi_b - dpsi_a) + q V (dp - dn) + C_ox (dV_G - dpsi_0) = 0
// with dpsi = 0 on the substrate contact is two tridiagonal solves (dV_G = 1 and dphi = 1) and the
// constraint. Returns C = C_ox (1 - dpsi_0) [F/cm^2].
double frozen_minority_capacitance(const device::Device& d, const results::BiasPoint& p,
                                   double C_ox) {
    const mesh::Mesh& m = d.mesh();
    const std::size_t n = m.node_count();
    const double eps = physics::silicon_parameters.eps_r * base::eps0_F_per_cm;
    const double V_T = base::thermal_voltage(d.temperature_K());
    std::vector<double> diag(n, 0.0), lower(n, 0.0), upper(n, 0.0);
    for (const mesh::Edge& e : m.edges()) {
        const auto a = static_cast<std::size_t>(e.first), b = static_cast<std::size_t>(e.second);
        REQUIRE(b == a + 1);
        const double g = eps * e.coupling_area / e.length;
        diag[a] -= g;
        diag[b] -= g;
        upper[a] = g;
        lower[b] = g;
    }
    for (std::size_t i = 0; i < n; ++i) {
        diag[i] -= base::q_C * m.volumes()[i] * (p.fields.p_cm3[i] + p.fields.n_cm3[i]) / V_T;
    }
    diag[0] -= C_ox;
    const auto solve = [&](double gate, double phi) {
        std::vector<double> rhs(n), c(n), x(n);
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] = -base::q_C * m.volumes()[i] * p.fields.n_cm3[i] * phi / V_T;
        }
        rhs[0] -= C_ox * gate;
        // The substrate node: dpsi = 0.
        std::vector<double> dg = diag, lo = lower, up = upper;
        dg[n - 1] = 1.0;
        lo[n - 1] = 0.0;
        rhs[n - 1] = 0.0;
        // Thomas algorithm.
        c[0] = up[0] / dg[0];
        rhs[0] /= dg[0];
        for (std::size_t i = 1; i < n; ++i) {
            const double den = dg[i] - lo[i] * c[i - 1];
            c[i] = up[i] / den;
            rhs[i] = (rhs[i] - lo[i] * rhs[i - 1]) / den;
        }
        x[n - 1] = rhs[n - 1];
        for (std::size_t i = n - 1; i-- > 0;) x[i] = rhs[i] - c[i] * x[i + 1];
        return x;
    };
    const std::vector<double> u = solve(1.0, 0.0), w = solve(0.0, 1.0);
    double nu = 0.0, nw = 0.0;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const double vn = m.volumes()[i] * p.fields.n_cm3[i];
        nu += vn * u[i];
        nw += vn * (w[i] - 1.0);
    }
    const double phi = -nu / nw;
    return C_ox * (1.0 - (u[0] + phi * w[0]));
}

// A 2D gated diode: a lumped gate (10 nm, metal 4.1 eV) over x in [1, 4] um of 1e16 p-silicon,
// an n+ region (1e19) on x < 1 um, y < 0.2 um with the source contact on its top, the body
// contact on the bottom (2 um deep).
device::Device gated_diode() {
    std::vector<double> x;
    for (int k = 0; k <= 40; ++k) x.push_back(4e-4 * k / 40);
    const auto y = legacy_graded_mesh(2e-4, {0.0}, 1e-7, 1e-5);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if (p[0] < 1e-4 && p[1] < 2e-5) donors[i] = 1e19;
    }
    std::vector<mesh::NodeId> source, gate;
    for (const mesh::NodeId v : m.find_boundary("y_min")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        if (px <= 0.5e-4) source.push_back(v);
        if (px >= 1e-4) gate.push_back(v);
    }
    auto body = m.find_boundary("y_max")->nodes;
    device::Contact g{"gate", device::ContactKind::gate, std::move(gate)};
    g.gate = {.boundary = "y_min",
              .oxide_thickness_cm = 1e-6,
              .electrode = device::GateElectrode::metal,
              .work_function_eV = 4.1};
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::vector<double>(n, 1e16),
         .contacts = {std::move(g), {"source", device::ContactKind::ohmic, std::move(source)},
                      {"body", device::ContactKind::ohmic, std::move(body)}}});
}

// The single-level trap fixture of the conductance method: 10 nm of oxide (meshed) on 2 um of
// 1e16 p-silicon at -0.7 V (depletion near flat band, the surface holes some 3e13 cm^-3 and the
// electrons negligible), band-gap narrowing off so that n_ie is the intrinsic density. The level
// (acceptor, N_t = 1e10 cm^-2, sigma 1e-15 cm^2, v_th 1e7 cm/s) sits where its hole emission
// rate equals the surface holes' capture rate (p1 = p_s, f near 1/2). With it: the occupancy f
// (from the trapped charge), tau = 1 / (c_p (p_I + p1)) = f / (c_p p1) and
// C_it = q N_t f (1 - f) / V_T.
struct TrapFixture {
    device::Device device;
    solve::BiasOptions options;
    double V0, C_ox, tau, C_it;
};

TrapFixture trap_fixture() {
    const double t_ox = 1e-6, L = 2e-4, N_A = 1e16, V0 = -0.7, T = 300.0;
    const double V_T = base::thermal_voltage(T);
    solve::BiasOptions o;
    o.models.bgn = false;
    const device::Device plain = moscap(t_ox, L, N_A, true);
    const auto p0 = *solve::solve_bias(plain, std::vector<double>{V0, 0.0}, o);
    std::size_t s = 0;  // the first silicon node
    while (plain.mesh().points()[s][0] < 0.0) ++s;
    const double n_i = physics::intrinsic_density(physics::silicon(), T);
    const double E_t = -V_T * std::log(p0.fields.p_cm3[s] / n_i);  // p1 = p_s
    const double N_t = 1e10, sigma = 1e-15, v_th = 1e7;
    device::Interface f{"silicon", "oxide"};
    f.traps.levels = {{.type = physics::TrapType::acceptor, .density_cm2 = N_t, .energy_eV = E_t,
                       .sigma_n_cm2 = sigma, .sigma_p_cm2 = sigma}};
    device::Device d = moscap(t_ox, L, N_A, true, {f});
    const auto dc = *solve::solve_bias(d, std::vector<double>{V0, 0.0}, o);
    const double occupied = -dc.interface_trap_charge[0] / (base::q_C * N_t);
    const double p1 = n_i * std::exp(-E_t / V_T);
    return {std::move(d), o, V0, 3.9 * base::eps0_F_per_cm / t_ox,
            occupied / (sigma * v_th * p1),
            base::q_C * N_t * occupied * (1.0 - occupied) / V_T};
}

// The semiconductor admittance behind the oxide: 1 / (1 / Y - 1 / (i omega C_ox)).
Complex behind_oxide(Complex y, double omega, double C_ox) {
    return 1.0 / (1.0 / y - 1.0 / (Complex{0.0, omega * C_ox}));
}

// The first Fourier component of the total current of contact c over the last period (its last
// `steps` fixed steps) of a transient run, by the trapezoid rule on the time points: Y A for a
// sine A sin(omega t).
Complex fundamental(const results::Transient& t, std::size_t c, double f, std::size_t steps) {
    const double T = 1.0 / f, w = two_pi * f;
    double a = 0.0, b = 0.0;
    for (std::size_t k = t.points.size() - steps; k < t.points.size(); ++k) {
        const double t0 = t.points[k - 1].time_s, t1 = t.points[k].time_s;
        const double i0 = t.points[k - 1].terminal_current[c];
        const double i1 = t.points[k].terminal_current[c];
        a += 0.5 * (t1 - t0) * (i0 * std::sin(w * t0) + i1 * std::sin(w * t1));
        b += 0.5 * (t1 - t0) * (i0 * std::cos(w * t0) + i1 * std::cos(w * t1));
    }
    return {2.0 * a / T, 2.0 * b / T};
}

// The admittance of a sine of A = 1 mV on contact 0 at frequency f, from transient runs of eight
// periods in fixed BDF2 steps (`steps` per period), one per step count.
std::vector<Complex> transient_admittance(const device::Device& d, double V0, double f,
                                          const solve::BiasOptions& steady,
                                          std::initializer_list<int> steps) {
    const double A = 1e-3, P = 1.0 / f;
    std::vector<Complex> y;
    for (const int n : steps) {
        const std::vector<solve::Waveform> w{*solve::Waveform::sine(V0, A, f),
                                             solve::Waveform::constant(0.0)};
        const solve::TransientOptions o{.steady = steady,
                                        .t_end_s = 8.0 * P,
                                        .dt_initial_s = P / n,
                                        .adaptive = false};
        const auto t = solve::solve_transient(d, w, o);
        REQUIRE(t.has_value());
        REQUIRE(!t->stopped);
        y.push_back(fundamental(*t, 0, f, static_cast<std::size_t>(n)) / A);
    }
    return y;
}

}  // namespace

TEST_CASE("small signal: a MOS capacitor in inversion, high frequency and its limit") {
    // 10 nm of oxide on 2 um of 1e16 p-silicon, metal gate 4.1 eV, at 1.5 V (strong inversion).
    // The operating point is the quasi-static one (drift-diffusion cannot reach it, Unit 12); the
    // small-signal system is drift-diffusion. At 1 MHz the minority carriers, supplied only by
    // generation and by diffusion from the substrate, cannot follow, and C is the frozen-minority
    // capacitance computed above: measured 3.1e-9 (accumulation-side holes and the inversion
    // layer's redistribution are in both). Below about 1 Hz the response of the isolated
    // inversion layer is the minority current across the depleted region, which the (psi, n, p)
    // state does not resolve (as the DC current there, Unit 12): the run reports it, |Y| below
    // its resolution. The resolved low-frequency curve needs a source of minority carriers (the
    // gated diode below).
    const double t_ox = 1e-6, L = 2e-4, N_A = 1e16, V0 = 1.5;
    const double C_ox = 3.9 * base::eps0_F_per_cm / t_ox;
    solve::BiasOptions qs;
    qs.equations = solve::Equations::equilibrium_poisson;
    const device::Device d = moscap(t_ox, L, N_A, false);
    const auto r = run(d, {{V0, 0.0}}, {1e-4, 1e6}, qs);
    const double C_hf = frozen_minority_capacitance(d, r.points[0].dc, C_ox);
    CAPTURE(C_hf, r.capacitance(0, 1, 0, 0) / C_hf - 1.0);
    REQUIRE(std::abs(r.capacitance(0, 1, 0, 0) / C_hf - 1.0) < 1e-6);
    REQUIRE(std::abs(r.admittance(0, 1, 0, 0)) > 1e3 * r.points[0].resolution[1][0]);
    REQUIRE(std::abs(r.admittance(0, 0, 0, 0)) < r.points[0].resolution[0][0]);
}

TEST_CASE("small signal: a gated diode gives the low-frequency capacitance") {
    // With an n+ source beside the gate the inversion layer is fed through the channel, and from
    // 10 Hz to 100 kHz the gate capacitance is the quasi-static dQ/dV (measured to below 1e-6,
    // the difference's own error) in accumulation, depletion and inversion. It is resolved from
    // 1 kHz up; at 10 Hz the resolution, a bound dominated by the conduction terms, is above the
    // capacitive current, which still matches. The operating points are quasi-static (every ohmic
    // contact at 0 V).
    const device::Device d = gated_diode();
    solve::BiasOptions qs;
    qs.equations = solve::Equations::equilibrium_poisson;
    const std::vector<double> f{10.0, 1e3, 1e5};
    for (const double V0 : {-1.5, -0.5, 1.5}) {
        const auto q_at = [&](double V) {
            return solve::solve_bias(d, std::vector<double>{V, 0.0, 0.0}, qs)->gate_charge[0];
        };
        const double h = 1e-3;
        const double C_qs =
            (8.0 * (q_at(V0 + h) - q_at(V0 - h)) - (q_at(V0 + 2.0 * h) - q_at(V0 - 2.0 * h))) /
            (12.0 * h);
        const auto r = run(d, {{V0, 0.0, 0.0}}, f, qs);
        for (std::size_t k = 0; k < f.size(); ++k) {
            CAPTURE(V0, f[k], r.capacitance(0, k, 0, 0) / C_qs - 1.0);
            REQUIRE(std::abs(r.capacitance(0, k, 0, 0) / C_qs - 1.0) < 1e-5);
            if (f[k] >= 1e3) {
                REQUIRE(std::abs(r.admittance(0, k, 0, 0)) > 10.0 * r.points[0].resolution[k][0]);
            }
        }
    }
}

TEST_CASE("small signal: an interface trap level, the conductance method") {
    // Nicollian and Goetzberger: the semiconductor admittance behind the oxide has
    // G_p / omega = C_it omega tau / (1 + omega^2 tau^2), peaking at omega tau = 1 at C_it / 2,
    // and its capacitance falls by C_it from low to high frequency. Measured from omega tau = 0.01
    // to 3 within 2.5e-3 of the peak value (0.09% at the peak); above that the substrate series
    // resistance adds omega^2 R C^2, which the method corrects separately. The capacitance step:
    // 2.1e-4 of C_it.
    const TrapFixture t = trap_fixture();
    std::vector<double> freq;
    const std::vector<double> wt{0.01, 0.1, 0.3, 1.0, 3.0, 100.0};
    for (const double v : wt) freq.push_back(v / (two_pi * t.tau));
    const auto r = run(t.device, {{t.V0, 0.0}}, freq, t.options);
    double worst = 0.0;
    std::vector<double> Gp, Cs;
    for (std::size_t k = 0; k < freq.size(); ++k) {
        const double w = two_pi * freq[k];
        const Complex ys = behind_oxide(r.admittance(0, k, 0, 0), w, t.C_ox);
        Gp.push_back(ys.real() / w);
        Cs.push_back(ys.imag() / w);
        if (wt[k] <= 3.0) {
            const double expected = t.C_it * wt[k] / (1.0 + wt[k] * wt[k]);
            worst = std::max(worst, std::abs(Gp.back() - expected) / (0.5 * t.C_it));
        }
    }
    CAPTURE(t.tau, t.C_it, Gp, worst);
    REQUIRE(worst < 5e-3);
    REQUIRE(Gp[3] > Gp[2]);  // the peak at omega tau = 1
    REQUIRE(Gp[3] > Gp[4]);
    REQUIRE(std::abs((Cs.front() - Cs.back()) / t.C_it - 1.0) < 1e-3);
}

TEST_CASE("small signal: analysis::conductance_peak finds the trap level (Unit 24)") {
    // The same level, a quarter decade apart from omega tau = 0.1 to 3.2 and off-centre by 0.1
    // decade: the extractor's peak against C_it / 2 at omega tau = 1 (the synthetic
    // tests/analysis/cv_test.cpp case measured 0.31% and 0.41% on this sampling).
    const TrapFixture t = trap_fixture();
    std::vector<double> freq;
    for (int k = -4; k <= 2; ++k) freq.push_back(std::pow(10.0, 0.25 * k + 0.1) / (two_pi * t.tau));
    const auto r = run(t.device, {{t.V0, 0.0}}, freq, t.options);
    const auto peak = analysis::conductance_peak(r, 0, 0, t.C_ox);
    REQUIRE(peak.has_value());
    const double value = peak->conductance_over_omega.value / (0.5 * t.C_it) - 1.0;
    const double where = peak->frequency.value * two_pi * t.tau - 1.0;
    std::printf("trap level: conductance peak %.3e off, frequency %.3e off\n", value, where);
    REQUIRE(std::abs(value) < 5e-3);
    REQUIRE(std::abs(where) < 1e-2);
}

TEST_CASE("small signal: transient sine responses converge to the admittance", "[.transient]") {
    // A 1 mV sine through Unit 21's transient run, eight periods in fixed BDF2 steps, the
    // fundamental of the last period: it converges to Y at second order (ratio 4 per halving).
    // The accumulated MOS-C at omega tau = 1 (RC): 2.3e-4, 5.8e-5, 1.5e-5 at 200, 400, 800 steps
    // per period. The trap level of the conductance method at omega tau = 1 (the traps'
    // dynamics in time against their elimination in the frequency domain): 3.1e-4, 7.8e-5,
    // 2.1e-5.
    {
        const double L = 1e-2, N_A = 1e16, V0 = -2.0;
        const device::Device d = moscap(1e-6, L, N_A, false);
        const auto q_at = [&](double V) {
            return solve::solve_bias(d, std::vector<double>{V, 0.0})->gate_charge[0];
        };
        const double C = (q_at(V0 + 1e-3) - q_at(V0 - 1e-3)) / 2e-3;
        const physics::Semiconductor si = physics::silicon();
        const double mu = physics::caughey_thomas_mobility(si, physics::Carrier::hole, N_A, 300.0);
        const double tau = L / (base::q_C * mu * N_A) *
                           (C + si.parameters().eps_r * base::eps0_F_per_cm / L);
        const double f = 1.0 / (two_pi * tau);
        const Complex ya = run(d, {{V0, 0.0}}, {f}).admittance(0, 0, 0, 0);
        const auto yt = transient_admittance(d, V0, f, {}, {200, 400, 800});
        const double e[3] = {std::abs(yt[0] / ya - 1.0), std::abs(yt[1] / ya - 1.0),
                             std::abs(yt[2] / ya - 1.0)};
        CAPTURE(e[0], e[1], e[2]);
        REQUIRE(e[0] / e[1] > 3.5);
        REQUIRE(e[1] / e[2] > 3.5);
        REQUIRE(e[2] < 3e-5);
    }
    {
        const TrapFixture t = trap_fixture();
        const double f = 1.0 / (two_pi * t.tau);
        const Complex ya = run(t.device, {{t.V0, 0.0}}, {f}, t.options).admittance(0, 0, 0, 0);
        const auto yt = transient_admittance(t.device, t.V0, f, t.options, {200, 400, 800});
        const double e[3] = {std::abs(yt[0] / ya - 1.0), std::abs(yt[1] / ya - 1.0),
                             std::abs(yt[2] / ya - 1.0)};
        CAPTURE(e[0], e[1], e[2]);
        REQUIRE(e[0] / e[1] > 3.5);
        REQUIRE(e[1] / e[2] > 3.5);
        REQUIRE(e[2] < 4e-5);
    }
}

namespace {

// An abrupt p+n junction (N_A on x < 1 um, N_D beyond, 6 um): anode (contact 0) on x_min.
device::Device p_plus_n(double NA, double ND) {
    const auto x = legacy_graded_mesh(6e-4, {1e-4}, 1e-8, 2e-6);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 1e-4) {
            acceptors[i] = NA;
        } else {
            donors[i] = ND;
        }
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

}  // namespace

TEST_CASE("small signal: a reverse-biased junction, depletion capacitance") {
    // An abrupt p+n junction (1e19 / 1e16) at -2 V: C = sqrt(q eps N / (2 (V_bi - V - 2 V_T))),
    // 1/N = 1/N_A + 1/N_D (the depletion approximation with the majority-carrier tails, Sze):
    // measured 0.88% above it at 1 kHz and 1 MHz alike (the approximation's own error).
    const double NA = 1e19, ND = 1e16;
    const device::Device d = p_plus_n(NA, ND);
    solve::BiasOptions o;
    o.models.bgn = false;
    const auto r = run(d, {{-2.0, 0.0}}, {1e3, 1e6}, o);
    const double V_T = base::thermal_voltage(300.0);
    const double n_i = physics::intrinsic_density(physics::silicon(), 300.0);
    const double Vbi = V_T * std::log(NA * ND / (n_i * n_i));
    const double N = 1.0 / (1.0 / NA + 1.0 / ND);
    const double eps = physics::silicon_parameters.eps_r * base::eps0_F_per_cm;
    const double C = std::sqrt(base::q_C * eps * N / (2.0 * (Vbi + 2.0 - 2.0 * V_T)));
    for (std::size_t k = 0; k < 2; ++k) {
        CAPTURE(k, r.capacitance(0, k, 0, 0) / C - 1.0);
        REQUIRE(std::abs(r.capacitance(0, k, 0, 0) / C - 1.0) < 1.5e-2);
    }
    REQUIRE(std::abs(r.capacitance(0, 0, 0, 0) / r.capacitance(0, 1, 0, 0) - 1.0) < 1e-6);
}

TEST_CASE("small signal: analysis::doping_profile of a reverse-biased junction (Unit 24)") {
    // The same junction from -0.5 to -4 V at 1 kHz: d(1/C^2)/dV gives the n side's doping at the
    // depletion edge, w = eps / C. The majority tails shift 1/C^2 by a constant to first order and
    // leave a second-order error of order (L_D / w)^2 (L_D = 41 nm, w = 10 to 19 L_D here): the
    // profile reads below N and approaches it with depth (measured 1.9% low at -0.5 V, 0.5% at -4 V).
    const double NA = 1e19, ND = 1e16;
    solve::BiasOptions o;
    o.models.bgn = false;
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 14; ++k) points.push_back({-0.5 - 0.25 * k, 0.0});
    const auto r = run(p_plus_n(NA, ND), points, {1e3}, o);
    const auto C = analysis::capacitance_curve(r, 0, 0, 0, 0);
    REQUIRE(C.has_value());
    const double eps = physics::silicon_parameters.eps_r * base::eps0_F_per_cm;
    const auto p = analysis::doping_profile(*C, eps);
    REQUIRE(p.has_value());
    const double N = 1.0 / (1.0 / NA + 1.0 / ND);
    double worst = 0.0;
    for (const double v : p->doping_cm3) worst = std::max(worst, std::abs(v / N - 1.0));
    std::printf("junction profile: doping within %.3e of N over depth %.3e to %.3e cm\n", worst,
                p->depth_cm.front(), p->depth_cm.back());
    CAPTURE(p->doping_cm3);
    REQUIRE(worst < 2.5e-2);
    REQUIRE(std::abs(p->doping_cm3.back() / N - 1.0) < 1e-2);
    for (std::size_t k = 1; k < p->doping_cm3.size(); ++k) {
        REQUIRE(p->doping_cm3[k] < N);
        REQUIRE(p->doping_cm3[k] > p->doping_cm3[k - 1]);
    }
    REQUIRE(p->depth_cm.back() > p->depth_cm.front());
}

TEST_CASE("small signal: a long diode, the diffusion admittance") {
    // A p+n diode (1e19 / 1e16, band-gap narrowing and Auger off), lifetimes 1 us everywhere,
    // 300 um of n side (9.4 diffusion lengths) at 0.5 V (low injection): the hole diffusion
    // admittance is G_0 sqrt(1 + i omega tau) (Shockley), so Re Y / Re Y(0) =
    // Re sqrt(1 + i omega tau), the junction capacitance being imaginary. Measured from
    // omega tau = 0.01 to 3 within 3.2e-3 (the electrons injected into the p+ side and the
    // recombination in the depleted region); above that the n side series resistance with the
    // junction capacitance adds omega^2 R C_j^2.
    physics::SemiconductorParameters si = physics::silicon_parameters;
    si.lifetime.tau_n0 = 1e-6;
    si.lifetime.tau_p0 = 1e-6;
    si.lifetime.N_ref = 1e30;  // the same lifetime at every doping
    const double W = 3e-2;
    const auto x = legacy_graded_mesh(1e-4 + W, {1e-4}, 1e-8, W / 400.0);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 1e-4) {
            acceptors[i] = 1e19;
        } else {
            donors[i] = 1e16;
        }
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(si)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    solve::BiasOptions o;
    o.models.bgn = false;
    o.models.auger = false;
    const double tau = 1e-6;
    const std::vector<double> wt{0.01, 0.1, 0.3, 1.0, 3.0};
    std::vector<double> f{0.0};
    for (const double v : wt) f.push_back(v / (two_pi * tau));
    const auto r = run(d, {{0.5, 0.0}}, f, o);
    const double G0 = r.conductance(0, 0, 0, 0);
    double worst = 0.0;
    for (std::size_t k = 0; k < wt.size(); ++k) {
        const double expected = std::sqrt(Complex{1.0, wt[k]}).real();
        worst = std::max(worst, std::abs(r.conductance(0, k + 1, 0, 0) / G0 / expected - 1.0));
    }
    CAPTURE(worst);
    REQUIRE(worst < 5e-3);
}
