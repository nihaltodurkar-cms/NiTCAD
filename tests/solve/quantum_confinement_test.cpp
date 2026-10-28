// MLDA quantum confinement on MOS capacitors (Unit 47; DECISIONS.md QC5, QC7-QC10): equilibrium
// with the model on carries no current, n p = n_ie^2 f_n f_p and nothing recombines; 2D = 1D; and,
// Release only (hidden tag [.mlda]): the inversion layer of a 1D MOS-C with a meshed 5 nm oxide
// (the density peak inside the silicon, lower inversion charge and capacitance, the threshold
// shift, accumulation), its convergence with the surface spacing, and a lumped gate against the
// meshed oxide.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/quantum_confinement.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "NiTCAD/solve/waveform.hpp"

using namespace NiTCAD;

namespace {

constexpr double t_ox = 5e-7;  // [cm]
constexpr double N_A = 1e17;   // [cm^-3]

assemble::PhysicsModels models(bool on, bool fermi_dirac = false) {
    assemble::PhysicsModels m;
    m.fermi_dirac = fermi_dirac;
    if (on) m.quantum_confinement = assemble::QuantumConfinement::mlda;
    return m;
}

// The silicon axis from the interface (x = 0): spacing h to 3 nm, then growing by 1.1 to 5 nm
// spacing, to 300 nm.
std::vector<double> silicon_axis(double h) {
    std::vector<double> x;
    double v = 0.0;
    for (; v < 3e-7 - 1e-12; v += h) x.push_back(v);
    double step = h;
    while (v < 3e-5) {
        x.push_back(v);
        step = std::min(step * 1.1, 5e-7);
        v += step;
    }
    x.push_back(3e-5);
    return x;
}

// A p-silicon MOS-C, the substrate contact on x_max. Meshed: a 5 nm SiO2 region (0.25 nm spacing)
// with an n+ poly electrode on x_min, the interface at x = 0 between nodes at -+h/2, so the first
// silicon node lies h/2 inside. Lumped: an n+ poly gate of the same oxide on x_min, the first node
// on the interface. With `width`, a y-uniform 2D device of that width (4 nodes across).
// `recombining` puts traps and surface recombination on the meshed interface.
device::Device mosc(double h, bool lumped = false, double width = 0.0, bool recombining = false) {
    std::vector<double> x;
    if (lumped) {
        x = silicon_axis(h);
    } else {
        for (int k = 0; k <= 20; ++k) x.push_back(-t_ox + k * (t_ox / 20.0));
        for (const double v : silicon_axis(h)) {
            if (v > 0.0) x.push_back(v);
        }
        x = *mesh::straddle_interface(x, 0.0, h);
    }
    std::vector<double> t;
    for (int k = 0; k < 4; ++k) t.push_back(width * k / 3.0);
    mesh::Mesh m = width > 0.0 ? *mesh::make_tensor_grid(x, t) : *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 1);
    std::vector<double> acceptors(n, N_A);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 0.0) {
            region[i] = 0;
            acceptors[i] = 0.0;
        }
    }
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Contact gate{"gate", lumped ? device::ContactKind::gate : device::ContactKind::electrode,
                         m.find_boundary("x_min")->nodes};
    if (lumped) {
        gate.gate = {.boundary = "x_min", .oxide_thickness_cm = t_ox,
                     .oxide_relative_permittivity = physics::silicon_dioxide().parameters().eps_r};
    } else {
        gate.electrode = {.kind = device::GateElectrode::n_poly};
    }
    std::vector<device::Interface> interfaces;
    if (recombining) {
        device::Interface f{"oxide", "silicon"};
        f.traps.levels = {{.type = physics::TrapType::acceptor, .density_cm2 = 1e11,
                           .energy_eV = 0.15, .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 1e-15}};
        f.recombination_velocity_n_cm_s = 1e3;
        f.recombination_velocity_p_cm_s = 1e3;
        interfaces.push_back(std::move(f));
    }
    std::vector<device::Region> regions;
    regions.push_back({"oxide", physics::silicon_dioxide()});
    regions.push_back({"silicon", physics::silicon()});
    if (lumped) regions.erase(regions.begin());
    if (lumped) std::fill(region.begin(), region.end(), 0);
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

// The quasi-static sweep of V_G over `vg` (substrate at 0).
std::vector<results::BiasPoint> quasi_static(const device::Device& d, bool on,
                                             std::span<const double> vg) {
    std::vector<std::vector<double>> points;
    for (const double v : vg) points.push_back({v, 0.0});
    solve::BiasOptions o;
    o.models = models(on);
    o.equations = solve::Equations::equilibrium_poisson;
    auto s = solve::sweep_bias(d, points, o);
    if (!s) FAIL(s.error().message);
    REQUIRE_FALSE(s->stopped.has_value());
    REQUIRE(s->points.size() == vg.size());
    return std::move(s->points);
}

// The electron sheet charge q int n dx over the silicon [C/cm^2].
double inversion_charge(const device::Device& d, const results::BiasPoint& p) {
    double q = 0.0;
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        if (!d.is_insulator(static_cast<mesh::NodeId>(i))) q += p.fields.n_cm3[i] * d.mesh().volumes()[i];
    }
    return base::q_C * q;
}

std::vector<double> range(double a, double b, double step) {
    std::vector<double> v;
    for (int k = 0; a + k * step <= b + 1e-9; ++k) v.push_back(a + k * step);
    return v;
}

// V_FB of the classical device: the gate charge's zero crossing.
double flat_band(const device::Device& d) {
    const auto vg = range(-2.0, 1.0, 0.05);
    const auto pts = quasi_static(d, false, vg);
    double vfb = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 1; k < pts.size() && std::isnan(vfb); ++k) {
        const double a = pts[k - 1].gate_charge[0], b = pts[k].gate_charge[0];
        if ((a <= 0.0) != (b <= 0.0)) vfb = vg[k - 1] + (vg[k] - vg[k - 1]) * a / (a - b);
    }
    INFO("the flat band lies in -2 .. 1 V");
    REQUIRE(std::isfinite(vfb));
    return vfb;
}

}  // namespace

TEST_CASE("quantum confinement: equilibrium with the model on (QC5)") {
    for (const bool fd : {false, true}) {
        const device::Device d = mosc(0.1e-7, false, 0.0, true);
        const auto c = *assemble::mlda_factors(d, models(true, fd));
        const physics::Semiconductor si = physics::silicon();
        for (const double vg : {-2.0, 1.5}) {  // accumulation and inversion
            solve::BiasOptions o;
            o.models = models(true, fd);
            // From the quasi-static state at the same bias: drift-diffusion from the equilibrium
            // guess does not reach inversion with only a substrate contact (Unit 12, 6.4).
            solve::BiasOptions qs = o;
            qs.equations = solve::Equations::equilibrium_poisson;
            const std::vector<double> bias{vg, 0.0};
            const auto start = solve::solve_bias(d, bias, qs);
            if (!start) FAIL(start.error().message);
            const auto p = solve::solve_bias(d, bias, o, &start->fields);
            if (!p) FAIL(p.error().message);
            CAPTURE(fd, vg, p->terminal_current[1], p->terminal_current_resolution[1]);
            REQUIRE(std::abs(p->terminal_current[1]) <= p->terminal_current_resolution[1]);
            const auto& f = p->fields;
            for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
                if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
                const double n = f.n_cm3[i], pp = f.p_cm3[i];
                const double tau = physics::scharfetter_lifetime(si, physics::Carrier::electron, N_A);
                const double scale = n * pp / (tau * (n + pp));
                CAPTURE(i, n, pp, p->rates.srh[i], scale);
                REQUIRE(std::abs(p->rates.srh[i]) <= 1e-9 * scale);
                const double h = d.mesh().volumes()[i];
                REQUIRE(std::abs(p->rates.interface_n[i]) <= 1e-9 * 1e3 * n * pp / (n + pp) / h);
                if (fd) continue;
                // n p = n_ie^2 f_n f_p (Boltzmann; n_ie with band-gap narrowing at N_A), on the
                // equilibrium state. The drift-diffusion state's own product carries the rounding
                // of its minority density, up to 3.3e-12 here at p = 1e-4 cm^-3 (7.5e-13 without
                // the model; reported by the next case).
                const double nie = physics::effective_intrinsic_density(si, N_A, 300.0);
                const double np = start->fields.n_cm3[i] * start->fields.p_cm3[i];
                REQUIRE(std::abs(np - nie * nie * c.f_n[i] * c.f_p[i]) <=
                        1e-12 * nie * nie * c.f_n[i] * c.f_p[i]);
            }
        }
    }
}

// Reported, not gated: the worst relative n p error of the drift-diffusion and quasi-static states
// with and without the model.
TEST_CASE("quantum confinement: n p at equilibrium with and without the model (QC5)", "[.mlda]") {
    const device::Device d = mosc(0.1e-7, false, 0.0, true);
    const physics::Semiconductor si = physics::silicon();
    const double nie = physics::effective_intrinsic_density(si, N_A, 300.0);
    for (const bool on : {false, true}) {
        const auto c = *assemble::mlda_factors(d, models(true));
        for (const double vg : {-2.0, 1.5}) {
            solve::BiasOptions o;
            o.models = models(on);
            o.newton.tol_update = 1e-13;
            solve::BiasOptions qs = o;
            qs.equations = solve::Equations::equilibrium_poisson;
            const std::vector<double> bias{vg, 0.0};
            const auto start = solve::solve_bias(d, bias, qs);
            if (!start) FAIL(start.error().message);
            const auto p = solve::solve_bias(d, bias, o, &start->fields);
            if (!p) FAIL(p.error().message);
            double worst = 0.0, worst_qs = 0.0;
            std::size_t at = 0;
            for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
                if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
                const double target = nie * nie * (on ? c.f_n[i] * c.f_p[i] : 1.0);
                const double e = std::abs(p->fields.n_cm3[i] * p->fields.p_cm3[i] / target - 1.0);
                if (e > worst) {
                    worst = e;
                    at = i;
                }
                worst_qs = std::max(
                    worst_qs,
                    std::abs(start->fields.n_cm3[i] * start->fields.p_cm3[i] / target - 1.0));
            }
            std::printf("QC5 n p: model %s, V_G %.1f V: drift-diffusion worst %.2e (node %zu, p %.3e), "
                        "quasi-static %.2e\n",
                        on ? "on" : "off", vg, worst, at, p->fields.p_cm3[at], worst_qs);
        }
    }
}

TEST_CASE("quantum confinement: a y-uniform 2D MOS-C is the 1D one (QC10)") {
    const auto vg = range(-1.5, 1.5, 0.5);
    const device::Device one = mosc(0.1e-7);
    const double width = 1e-5;
    const device::Device two = mosc(0.1e-7, false, width);
    const auto a = quasi_static(one, true, vg);
    const auto b = quasi_static(two, true, vg);
    for (std::size_t k = 0; k < vg.size(); ++k) {
        CAPTURE(vg[k], a[k].gate_charge[0], b[k].gate_charge[0] / width);
        REQUIRE(std::abs(b[k].gate_charge[0] / width - a[k].gate_charge[0]) <=
                1e-12 * std::abs(a[k].gate_charge[0]));
    }
}

TEST_CASE("quantum confinement: a long transient step ends at the steady state (QC10)") {
    // The MOS-C's gate stepped from 0 into accumulation (-1.5 V): the majority holes follow within
    // the dielectric relaxation time, so by 1 us the run is at the DC state of -1.5 V.
    const device::Device d = mosc(0.1e-7);
    solve::TransientOptions o{.t_end_s = 1e-6, .dt_initial_s = 1e-16};
    o.steady.models = models(true);
    const std::vector<solve::Waveform> wf{*solve::Waveform::step(0.0, -1.5, 0.0),
                                          solve::Waveform::constant(0.0)};
    const auto run = solve::solve_transient(d, wf, o);
    if (!run) FAIL(run.error().message);
    REQUIRE_FALSE(run->stopped.has_value());
    solve::BiasOptions b;
    b.models = models(true);
    const auto dc = solve::solve_bias(d, std::vector<double>{-1.5, 0.0}, b);
    if (!dc) FAIL(dc.error().message);
    const auto& f = run->snapshots.back().fields;
    for (std::size_t i = 0; i < f.potential_V.size(); ++i) {
        CAPTURE(i);
        REQUIRE(std::abs(f.potential_V[i] - dc->fields.potential_V[i]) <= 1e-9);
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        REQUIRE(std::abs(f.p_cm3[i] / dc->fields.p_cm3[i] - 1.0) <= 1e-7);
    }
}

TEST_CASE("quantum confinement: the inversion layer of a MOS-C (QC7)", "[.mlda]") {
    const device::Device d = mosc(0.05e-7);
    const double vfb = flat_band(d);
    const auto vg = range(vfb - 2.0, vfb + 2.5, 0.05);
    const auto classical = quasi_static(d, false, vg);
    const auto quantum = quasi_static(d, true, vg);
    std::printf("QC7: V_FB = %.4f V\n", vfb);
    // Strong inversion, V_G - V_FB = 2 V: the density peak and the charge.
    const auto at = [&](double v) {
        return static_cast<std::size_t>(std::lround((v - vg.front()) / 0.05));
    };
    const std::size_t k2 = at(vfb + 2.0);
    const auto peak = [&](const results::BiasPoint& p) {
        double best = -1.0, where = 0.0;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
            if (p.fields.n_cm3[i] > best) {
                best = p.fields.n_cm3[i];
                where = d.mesh().points()[i][0];
            }
        }
        return where;
    };
    const double z_classical = peak(classical[k2]), z_quantum = peak(quantum[k2]);
    const double q_classical = inversion_charge(d, classical[k2]);
    const double q_quantum = inversion_charge(d, quantum[k2]);
    std::printf("QC7: peak %.3f nm (classical %.3f nm); Q_inv %.4e (classical %.4e) C/cm^2\n",
                z_quantum * 1e7, z_classical * 1e7, q_quantum, q_classical);
    CHECK(z_classical < 0.1e-7);
    CHECK(z_quantum >= 0.3e-7);
    CHECK(z_quantum <= 2e-7);
    CHECK(q_quantum < q_classical);
    // The quasi-static gate capacitance at V_G - V_FB = 2 V.
    const auto capacitance = [&](const std::vector<results::BiasPoint>& p) {
        return (p[k2 + 1].gate_charge[0] - p[k2 - 1].gate_charge[0]) / 0.1;
    };
    const double c_classical = capacitance(classical), c_quantum = capacitance(quantum);
    const double c_ox = physics::silicon_dioxide().parameters().eps_r *
                        base::eps0_F_per_cm / t_ox;
    std::printf("QC7: C/C_ox %.4f (classical %.4f)\n", c_quantum / c_ox, c_classical / c_ox);
    CHECK(c_quantum < c_classical);
    // The threshold at Q_inv = 1e-7 C/cm^2, interpolated in ln Q.
    const auto threshold = [&](const std::vector<results::BiasPoint>& p) {
        double v = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t k = 1; k < p.size() && std::isnan(v); ++k) {
            const double a = inversion_charge(d, p[k - 1]), b = inversion_charge(d, p[k]);
            if (a < 1e-7 && b >= 1e-7) v = vg[k - 1] + 0.05 * std::log(1e-7 / a) / std::log(b / a);
        }
        INFO("the sweep crosses the threshold");
        REQUIRE(std::isfinite(v));
        return v;
    };
    const double shift = threshold(quantum) - threshold(classical);
    std::printf("QC7: threshold shift %.1f mV\n", shift * 1e3);
    CHECK(shift >= 0.010);
    CHECK(shift <= 0.150);
    // Accumulation, V_G - V_FB = -2 V: less hole charge, so less gate charge.
    const std::size_t ka = at(vfb - 2.0);
    std::printf("QC7: accumulation Q_G %.4e (classical %.4e) C/cm^2\n", quantum[ka].gate_charge[0],
                classical[ka].gate_charge[0]);
    CHECK(quantum[ka].gate_charge[0] < 0.0);
    CHECK(std::abs(quantum[ka].gate_charge[0]) < std::abs(classical[ka].gate_charge[0]));
}

TEST_CASE("quantum confinement: convergence with the surface spacing (QC8)", "[.mlda]") {
    const double vfb = flat_band(mosc(0.05e-7));
    std::vector<double> q;
    for (const double h : {0.2e-7, 0.1e-7, 0.05e-7}) {
        const device::Device d = mosc(h);
        const std::vector<double> vg{vfb + 2.0};
        q.push_back(inversion_charge(d, quasi_static(d, true, vg)[0]));
        std::printf("QC8: h = %.2f nm, Q_inv = %.6e C/cm^2\n", h * 1e7, q.back());
    }
    const double d1 = std::abs(q[1] - q[0]), d2 = std::abs(q[2] - q[1]);
    CAPTURE(q[0], q[1], q[2], d1, d2);
    CHECK(d2 <= 0.01 * q[2]);
    CHECK(d2 < d1);
}

TEST_CASE("quantum confinement: a lumped gate against the meshed oxide (QC9)", "[.mlda]") {
    const device::Device meshed = mosc(0.05e-7);
    const device::Device lumped = mosc(0.05e-7, true);
    const double vfb = flat_band(meshed);
    const std::vector<double> vg{vfb + 2.0};
    const double qm = inversion_charge(meshed, quasi_static(meshed, true, vg)[0]);
    const double ql = inversion_charge(lumped, quasi_static(lumped, true, vg)[0]);
    std::printf("QC9: Q_inv meshed %.5e, lumped %.5e C/cm^2 (ratio %.4f)\n", qm, ql, ql / qm);
    CHECK(std::abs(ql - qm) <= 0.05 * qm);
}
