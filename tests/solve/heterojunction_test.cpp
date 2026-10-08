// Heterojunctions in the device solves (ARCHITECTURE.md section 11, Unit 15; legacy
// tests/test_m33_interface.py: G1 detailed balance per carrier, G2 the affinity step's effect and
// direction, the isotype barrier, S2 thermionic emission; the exact equilibrium of an abrupt
// heterojunction; composition with Fermi-Dirac statistics and with a gate).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/fermi_dirac.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "legacy_graded_mesh.hpp"
#include "../physics/references.hpp"

using namespace NiTCAD;

namespace {

constexpr double T = 300.0;

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// The legacy M33 mesh: graded_mesh(2e-4, [1e-4], 1e-8, 1e-6, 1.12).
std::vector<double> m33_axis() { return legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6, 1.12); }

constexpr auto thermionic = device::InterfaceTransport::thermionic_emission;

// `left` on x < 1 um with net doping C_left, `right` beyond with C_right (positive: donors), ohmic
// contacts on x_min and x_max, the interface declared with `transport`. With `transverse`, a
// y-uniform 2D device of that width.
device::Device junction(const physics::SemiconductorParameters& left,
                        const physics::SemiconductorParameters& right, double C_left,
                        double C_right, const std::vector<double>& transverse = {},
                        device::InterfaceTransport transport =
                            device::InterfaceTransport::drift_diffusion) {
    mesh::Mesh m = transverse.empty() ? *mesh::make_tensor_grid(m33_axis())
                                      : *mesh::make_tensor_grid(m33_axis(), transverse);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const bool is_left = m.points()[i][0] < 1e-4;
        region[i] = is_left ? 0 : 1;
        const double C = is_left ? C_left : C_right;
        (C > 0.0 ? donors[i] : acceptors[i]) = std::abs(C);
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"left", *physics::Semiconductor::create(left)},
                     {"right", *physics::Semiconductor::create(right)}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .interfaces = {{"left", "right", transport}}});
}

physics::SemiconductorParameters silicon_with(double chi_eV, double Eg0_eV = 1.17) {
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.electron_affinity_eV = chi_eV;
    p.Eg0_eV = Eg0_eV;
    return p;
}

// The legacy M33 models: band-gap narrowing off, everything else default.
constexpr assemble::PhysicsModels m33{.bgn = false};

double current(const device::Device& d, double bias, const assemble::PhysicsModels& models = m33) {
    const std::vector<double> v{bias, 0.0};
    const auto s = solve::solve_bias(d, v, {.models = models});
    REQUIRE(s.has_value());
    return s->terminal_current[0];
}

// Index of the last node left of the interface.
std::size_t interface_node(const std::vector<double>& x) {
    std::size_t k = 0;
    while (x[k + 1] < 1e-4) ++k;
    return k;
}

// The exact Boltzmann equilibrium of an abrupt heterojunction (no band-gap narrowing). On each side
// the first integral of Poisson's equation gives D^2 = 2 q V_T eps (n_b (e^u - 1 - u) +
// p_b (e^-u - 1 + u)), u = (phi - phi_b) / V_T, with n_b, p_b the side's neutral densities; D is
// continuous at the interface and phi_b2 - phi_b1 follows from a common Fermi level:
// E_vac = E_F - kT ln(n_b / Nc) + chi on each side, phi = -E_vac / q. Returns phi(0) - phi_b1
// and D(0).
struct Interface {
    double phi0;  // interface potential relative to the left bulk [V]
    double D;     // |D| at the interface [C/cm^2]
    double vbi;   // phi_b2 - phi_b1 [V]
};

Interface exact_interface(const physics::SemiconductorParameters& a,
                          const physics::SemiconductorParameters& b, double Ca, double Cb) {
    const double VT = base::thermal_voltage(T);
    struct Side {
        double eps, nb, pb, Nc, chi;
    };
    const auto side = [&](const physics::SemiconductorParameters& p, double C) {
        const physics::Semiconductor m = *physics::Semiconductor::create(p);
        const auto e = physics::boltzmann_neutral_equilibrium(C, physics::intrinsic_density(m, T));
        return Side{p.eps_r * base::eps0_F_per_cm, e.n, e.p, physics::conduction_band_dos(m, T),
                    p.electron_affinity_eV};
    };
    const Side l = side(a, Ca), r = side(b, Cb);
    const double vbi = (VT * std::log(r.nb / r.Nc) - r.chi) - (VT * std::log(l.nb / l.Nc) - l.chi);
    const auto D2 = [&](const Side& s, double u) {
        return 2.0 * base::q_C * VT * s.eps *
               (s.nb * (std::expm1(u) - u) + s.pb * (std::expm1(-u) + u));
    };
    // phi0 between the bulks; D_left^2 rises from 0 and D_right^2 falls to 0 across that range.
    double lo = std::min(0.0, vbi), hi = std::max(0.0, vbi);
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double g = D2(l, mid / VT) - D2(r, (mid - vbi) / VT);
        ((g < 0.0) == (vbi > 0.0) ? lo : hi) = mid;
    }
    const double phi0 = 0.5 * (lo + hi);
    return {phi0, std::sqrt(D2(l, phi0 / VT)), vbi};
}

}  // namespace

TEST_CASE("heterojunction: the solved interface matches the exact first integral") {
    // Abrupt junctions of GaAs and Al0.3Ga0.7As, which differ in every band parameter (chi by
    // 0.255 eV, Eg by 0.374 eV, Nc, Nv and eps_r): an isotype n-N junction and a p-N junction, plus
    // a pure affinity-and-permittivity step in silicon. Boltzmann, no band-gap narrowing. The
    // interface lies at the midpoint of the edge joining the regions; its potential is the
    // permittivity-weighted mean of the edge's ends (exact for the two linear half-edges), and D is
    // the harmonic-mean permittivity times the edge field.
    const physics::SemiconductorParameters gaas = physics::gallium_arsenide_parameters;
    const physics::SemiconductorParameters algaas = *physics::algaas_parameters(0.3);
    physics::SemiconductorParameters si2 = silicon_with(4.30);
    si2.eps_r = 13.0;
    struct Case {
        const char* name;
        physics::SemiconductorParameters a, b;
        double Ca, Cb;
    };
    const Case cases[] = {
        {"n-GaAs / N-AlGaAs", gaas, algaas, 1e17, 1e17},
        {"p-GaAs / N-AlGaAs", gaas, algaas, -1e17, 1e17},
        {"N-AlGaAs / n-GaAs", algaas, gaas, 1e17, 3e16},
        {"n-Si / n-Si(chi + 0.25, eps 13)", physics::silicon_parameters, si2, 1e17, 1e17},
    };
    const auto x = m33_axis();
    const std::size_t k = interface_node(x);
    for (const Case& c : cases) {
        CAPTURE(c.name);
        const device::Device d = junction(c.a, c.b, c.Ca, c.Cb);
        const auto eq = solve::solve_equilibrium(d, {.models = m33});
        REQUIRE(eq.has_value());
        const auto& phi = eq->fields.potential_V;
        const Interface exact = exact_interface(c.a, c.b, c.Ca, c.Cb);
        const double ea = c.a.eps_r, eb = c.b.eps_r;
        const double phi0 = (ea * phi[k] + eb * phi[k + 1]) / (ea + eb) - phi.front();
        const double D = base::eps0_F_per_cm * (2.0 * ea * eb / (ea + eb)) *
                         std::abs(phi[k + 1] - phi[k]) / (x[k + 1] - x[k]);
        UNSCOPED_INFO(c.name << ": V_bi " << phi.back() - phi.front() << " (exact " << exact.vbi
                             << "), phi(0) " << phi0 << " (exact " << exact.phi0 << "), D " << D
                             << " (exact " << exact.D << ")");
        REQUIRE(std::abs((phi.back() - phi.front()) - exact.vbi) <= 1e-12 * std::abs(exact.vbi));
        REQUIRE(std::abs(phi0 - exact.phi0) <= 2e-4 * std::abs(exact.vbi));
        REQUIRE(close(D, exact.D, 2e-3));
    }
}

TEST_CASE("heterojunction: equilibrium carries no current, per carrier (legacy G1)") {
    // The legacy fixture: p-Si / n-Si with the right side's chi shifted, and a pure gap step at
    // fixed chi (electrons see no barrier, holes the whole step). The fluxes at the equilibrium
    // state are rounding, each held to 64 eps of its one-sided terms (as Unit 14); thermionic
    // emission on the interface edge too.
    struct Case {
        double chi, Eg0;
    };
    const auto x = m33_axis();
    const double VT = base::thermal_voltage(T);
    for (const Case& c : {Case{3.75, 1.17}, Case{3.95, 1.17}, Case{4.05, 1.17}, Case{4.15, 1.17},
                          Case{4.35, 1.17}, Case{4.05, 0.97}, Case{4.05, 1.37}}) {
        for (const bool te : {false, true}) {
            CAPTURE(c.chi, c.Eg0, te);
            const device::Device d =
                junction(physics::silicon_parameters, silicon_with(c.chi, c.Eg0), -1e17, 1e17, {},
                         te ? thermionic : device::InterfaceTransport::drift_diffusion);
            const assemble::PhysicsModels models = m33;
            const auto eq = solve::solve_equilibrium(d, {.models = models});
            REQUIRE(eq.has_value());
            const auto scaling = *assemble::make_scaling(d);
            const auto system = *assemble::DriftDiffusion::create(d, scaling, models);
            std::vector<double> psi(x.size());
            for (std::size_t i = 0; i < x.size(); ++i) {
                psi[i] = eq->fields.potential_V[i] / scaling.V_T;
            }
            const auto currents = system.edge_currents(system.state_from_potential(psi));
            const auto& f = eq->fields;
            for (std::size_t e = 0; e + 1 < x.size(); ++e) {
                CAPTURE(e);
                const double drive = 3.0 + std::abs(psi[e + 1] - psi[e]) + 1.0 / VT;
                const double floor = 64.0 * 2.2204460492503131e-16 * base::q_C * 1360.0 * VT /
                                     (x[e + 1] - x[e]) * drive;
                REQUIRE(std::abs(currents[e].first * scaling.J0) <=
                        floor * (f.n_cm3[e] + f.n_cm3[e + 1]));
                REQUIRE(std::abs(currents[e].second * scaling.J0) <=
                        floor * (f.p_cm3[e] + f.p_cm3[e + 1]));
            }
            // The legacy gate on its own published currents: zero bias through solve_bias.
            const std::vector<double> zero{0.0, 0.0};
            const auto s = solve::solve_bias(d, zero, {.models = models});
            REQUIRE(s.has_value());
            for (std::size_t e = 0; e < s->edge_current_n.size(); ++e) {
                REQUIRE(std::abs(s->edge_current_n[e]) < 1e-7);
                REQUIRE(std::abs(s->edge_current_p[e]) < 1e-7);
            }
        }
    }
}

TEST_CASE("heterojunction: the affinity step moves the solution, in the physical direction (G2)") {
    // Legacy G2 on its own fixture (p-Si / n-Si, right chi varied, 0.4 V): raising the n side's chi
    // lowers its E_c, so the step electrons climb into the p side rises and the forward current
    // falls, monotonically; the legacy measured 2.749e-4 to 2.650e-4 A/cm^2 from 3.85 to 4.25 eV.
    // NiTCAD's band shift without band-gap narrowing is the legacy affinity gauge exactly, so the
    // values are reproduced.
    const double chis[] = {3.85, 3.95, 4.05, 4.15, 4.25};
    std::vector<double> J;
    for (const double chi : chis) {
        J.push_back(current(junction(physics::silicon_parameters, silicon_with(chi), -1e17, 1e17),
                            0.4));
    }
    UNSCOPED_INFO("J(0.4 V): " << J[0] << " ... " << J[4]);
    for (std::size_t k = 0; k + 1 < J.size(); ++k) REQUIRE(J[k] > J[k + 1]);
    REQUIRE(close(J[0], 2.749e-4, 2e-3));
    REQUIRE(close(J[4], 2.650e-4, 2e-3));
    // And psi moves: at equilibrium the potential differs by more than 1 mV somewhere.
    const auto a = solve::solve_equilibrium(
        junction(physics::silicon_parameters, silicon_with(4.05), -1e17, 1e17), {.models = m33});
    const auto b = solve::solve_equilibrium(
        junction(physics::silicon_parameters, silicon_with(4.25), -1e17, 1e17), {.models = m33});
    double moved = 0.0;
    for (std::size_t i = 0; i < a->fields.potential_V.size(); ++i) {
        moved = std::max(moved, std::abs(a->fields.potential_V[i] - b->fields.potential_V[i]));
    }
    REQUIRE(moved > 0.19);  // the right side's whole potential moves by the 0.2 eV step
}

TEST_CASE("heterojunction: an isotype step is a barrier of either sign (legacy G2)") {
    // n-Si / n-Si(chi), 1e17, 0.1 V: the majority current is largest at no step and falls for
    // either sign. The legacy measured 6.42e3 A/cm^2 flat, 4.34e3 at -0.20 eV, 3.11e3 at +0.20 eV.
    const auto run = [&](double chi) {
        return current(junction(physics::silicon_parameters, silicon_with(chi), 1e17, 1e17), 0.1);
    };
    const double flat = run(4.05), lo = run(3.85), hi = run(4.25);
    UNSCOPED_INFO("J flat " << flat << ", -0.2 eV " << lo << ", +0.2 eV " << hi);
    REQUIRE(close(flat, 6.42e3, 2e-3));
    REQUIRE(close(lo, 4.34e3, 2e-3));
    REQUIRE(close(hi, 3.11e3, 2e-3));
}

TEST_CASE("heterojunction: thermionic emission limits the current, more so behind a higher step") {
    // Legacy S2 fixture: n-Si / n-Si(chi), 1e17, 0.1 V, the interface declared thermionic. At the
    // real emission velocity the interface cuts the current relative to drift-diffusion, and the
    // cut deepens with the step.
    const auto ratio = [&](double chi) {
        const device::Device d = junction(physics::silicon_parameters, silicon_with(chi), 1e17,
                                          1e17);
        const device::Device t = junction(physics::silicon_parameters, silicon_with(chi), 1e17,
                                          1e17, {}, thermionic);
        return current(t, 0.1) / current(d, 0.1);
    };
    const double r0 = ratio(4.05), r1 = ratio(4.20), r2 = ratio(4.35);
    UNSCOPED_INFO("J(TE) / J(DD): " << r0 << ", " << r1 << ", " << r2);
    REQUIRE(r0 < 1.0);
    REQUIRE(r1 < r0);
    REQUIRE(r2 < r1);
    REQUIRE(r2 > 0.0);
}

TEST_CASE("heterojunction: a declared interface without a step is the emission resistance") {
    // Two silicon regions with identical parameters, declared thermionic: the interface is a pure
    // velocity limit, J = q v (n1 - n2), i.e. a series resistance V_T / (q v n) at small bias
    // (v the emission velocity, n = 1e17). In a uniform resistor the drift-diffusion edge it
    // replaces is a 1e-8 cm slice (resistance 1e-3 of it). So V / J(TE) - V / J(DD) = R_int.
    const device::Device d = junction(physics::silicon_parameters, physics::silicon_parameters,
                                      1e17, 1e17);
    const device::Device t = junction(physics::silicon_parameters, physics::silicon_parameters,
                                      1e17, 1e17, {}, thermionic);
    const double V = 0.01, VT = base::thermal_voltage(T);
    const double v = physics::emission_velocity_cm_s(physics::conduction_band_dos(
                                                         physics::silicon(), T),
                                                     T);
    const double R_int = VT / (base::q_C * v * 1e17);
    const double added = V / current(t, V) - V / current(d, V);
    UNSCOPED_INFO("added resistance " << added << " Ohm cm^2, V_T / (q v n) = " << R_int);
    REQUIRE(close(added, R_int, 5e-3));
}

TEST_CASE("heterojunction: thermionic emission becomes drift-diffusion as the velocity grows") {
    // The legacy G-5 limit, now through the materials' Richardson constants: an emitter that
    // fast limits nothing, so the current tends to drift-diffusion's (the legacy measured 1.0032
    // at its largest factor; an ideal interface slightly beats the drift-diffusion edge it
    // replaces), and the ratio rises monotonically.
    std::vector<double> ratios;
    for (const double A : {0.0, 1e4, 1e6, 1e8}) {
        physics::SemiconductorParameters a = physics::silicon_parameters, b = silicon_with(4.35);
        a.richardson = b.richardson = {.electron = A, .hole = A};
        const device::Device d = junction(a, b, 1e17, 1e17);
        const device::Device t = junction(a, b, 1e17, 1e17, {}, thermionic);
        ratios.push_back(current(t, 0.1) / current(d, 0.1));
    }
    UNSCOPED_INFO("J(TE) / J(DD) for A* = DOS, 1e4, 1e6, 1e8: " << ratios[0] << ", " << ratios[1]
                                                                 << ", " << ratios[2] << ", "
                                                                 << ratios[3]);
    REQUIRE(ratios[0] < ratios[1]);
    REQUIRE(ratios[1] < ratios[2]);
    REQUIRE(ratios[2] <= ratios[3] + 1e-9);
    REQUIRE(std::abs(ratios[3] - 1.0) < 0.01);
}

TEST_CASE("heterojunction: composes with Fermi-Dirac statistics and band-gap narrowing") {
    // The legacy refused band_offset="affinity" with fd: NiTCAD's shift only moves eta, so the
    // composition is exact. A degenerate n+ GaAs (1e19, eta_c > 0) against p-Al0.3Ga0.7As, with
    // thermionic emission: flat Fermi level at equilibrium and a converging forward sweep. The
    // AlGaAs minority electrons are 2e-11 cm^-3, 2e-30 of the largest density: the sweep converges
    // because the Newton measure floors densities at 1e-20 of the largest (6.2, Unit 15).
    const device::Device d = junction(physics::gallium_arsenide_parameters,
                                      *physics::algaas_parameters(0.3), 1e19, -1e17, {},
                                      thermionic);
    const assemble::PhysicsModels models{.fermi_dirac = true};
    const auto eq = solve::solve_equilibrium(d, {.models = models});
    REQUIRE(eq.has_value());
    // The electron quasi-Fermi level, from each node's own statistics, is flat to the solve.
    const auto x = m33_axis();
    const auto scaling = *assemble::make_scaling(d);
    double spread = 0.0;
    double first = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const physics::Semiconductor& m = d.material(static_cast<mesh::NodeId>(i));
        const double Nc = physics::conduction_band_dos(m, T);
        const double eta_c = physics::inverse_fermi_half(eq->fields.n_cm3[i] / Nc);
        // E_F = E_c + kT eta_c, E_c = E_vac - chi - dEg / 2 (band-gap narrowing is shared
        // between the bands), E_vac = -q phi up to a constant [eV].
        const double dEg =
            physics::bandgap_narrowing_eV(m, d.total_impurity(static_cast<mesh::NodeId>(i)));
        const double ef = scaling.V_T * eta_c - m.parameters().electron_affinity_eV - 0.5 * dEg -
                          eq->fields.potential_V[i];
        if (i == 0) first = ef;
        spread = std::max(spread, std::abs(ef - first));
    }
    UNSCOPED_INFO("electron Fermi level spread " << spread << " V");
    REQUIRE(spread < 1e-9);
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 6; ++k) points.push_back({-0.1 * k, 0.0});
    const auto sweep = solve::sweep_bias(d, points, {.models = models});
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
}

TEST_CASE("heterojunction: a gate on the shifted material sees its own flat band") {
    // A p-Si MOS capacitor whose substrate end (node 0, the reference material) is silicon with
    // chi + 0.3 eV: the gate (x_max) sits on plain silicon, where the band shift is -0.3 eV / V_T.
    // At V_G = Phi + V_T eta_bulk (the flat-band voltage of the silicon under the gate) the gate
    // charge is zero; without the shift in the electrode potential it would be C_ox 0.3 V.
    const auto x = m33_axis();
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) region[i] = m.points()[i][0] < 0.3e-4 ? 0 : 1;
    device::Contact gate{"gate", device::ContactKind::gate, m.find_boundary("x_max")->nodes};
    gate.gate = {.boundary = "x_max", .oxide_thickness_cm = 5e-7};
    auto substrate = m.find_boundary("x_min")->nodes;
    const physics::Semiconductor si = physics::silicon();
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"shifted", *physics::Semiconductor::create(silicon_with(4.35))},
                     {"silicon", si}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::vector<double>(n, 1e17),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
    const double VT = base::thermal_voltage(T);
    const double eta_bulk = physics::boltzmann_neutral_equilibrium(
                                -1e17, physics::intrinsic_density(si, T))
                                .eta;
    const double vfb =
        assemble::gate_intrinsic_offset_V(d.contacts()[0].gate, si, T) + VT * eta_bulk;
    const std::vector<std::vector<double>> points{{vfb, 0.0}, {vfb + 0.3, 0.0}};
    solve::BiasOptions o{.models = m33, .equations = solve::Equations::equilibrium_poisson};
    const auto sweep = solve::sweep_bias(d, points, o);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    const double cox = 3.9 * base::eps0_F_per_cm / 5e-7;
    const double q_fb = sweep->points[0].gate_charge[0];
    UNSCOPED_INFO("gate charge at flat band " << q_fb << " C/cm^2, C_ox 0.3 V = " << cox * 0.3);
    REQUIRE(std::abs(q_fb) < 1e-6 * cox * VT);
    REQUIRE(std::abs(sweep->points[1].gate_charge[0]) > 0.1 * cox * 0.3);
}

TEST_CASE("heterojunction: a y-uniform 2D device reproduces 1D") {
    const physics::SemiconductorParameters right = *physics::algaas_parameters(0.3);
    const std::vector<double> y{0.0, 0.5e-4, 1e-4};
    const auto d1 =
        junction(physics::gallium_arsenide_parameters, right, -1e17, 1e17, {}, thermionic);
    const auto d2 =
        junction(physics::gallium_arsenide_parameters, right, -1e17, 1e17, y, thermionic);
    const assemble::PhysicsModels models{};
    const double j1 = current(d1, 1.2, models);
    const double j2 = current(d2, 1.2, models);
    UNSCOPED_INFO("1D " << j1 << " A/cm^2, 2D " << j2 / 1e-4 << " A/cm^2 per width");
    REQUIRE(close(j2 / 1e-4, j1, 1e-9));
}

// Unit 15 fixes: graded compositions, incomplete ionization, radiative recombination, the band
// diagram, the current resolution, a non-planar interface.

namespace {

// n-type 1e17 GaAs graded to Al0.3Ga0.7As over [0.9, 1.1] um in `steps` regions (no interface
// declared), then Al0.3Ga0.7As to 2 um.
device::Device graded(int steps) {
    mesh::Mesh m = *mesh::make_tensor_grid(m33_axis());
    const std::size_t n = m.node_count();
    // Each node takes the composition of its step; a region is made only for the steps that hold
    // a node (a fine staircase has more steps than the mesh has nodes in the grade).
    std::vector<device::Region> regions;
    std::vector<int> region_of_step(static_cast<std::size_t>(steps) + 1, -1);
    std::vector<device::RegionId> region(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double f = std::clamp((m.points()[i][0] - 0.9e-4) / 0.2e-4, 0.0, 1.0);
        const auto k = static_cast<std::size_t>(std::lround(steps * f));
        if (region_of_step[k] < 0) {
            region_of_step[k] = static_cast<int>(regions.size());
            regions.push_back({"x" + std::to_string(k),
                               *physics::Semiconductor::create(*physics::algaas_parameters(
                                   0.3 * static_cast<double>(k) / steps))});
        }
        region[i] = region_of_step[k];
    }
    auto a = m.find_boundary("x_min")->nodes;
    auto c = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = std::move(regions),
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 1e17),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"a", device::ContactKind::ohmic, std::move(a)},
                      {"c", device::ContactKind::ohmic, std::move(c)}}});
}

// A uniform device of one material with donors and acceptors on 21 nodes over 1 um.
device::Device uniform(const physics::SemiconductorParameters& p, double ND, double NA,
                       double temperature) {
    std::vector<double> x;
    for (int i = 0; i <= 20; ++i) x.push_back(0.05e-4 * i);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    auto a = m.find_boundary("x_min")->nodes;
    auto c = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = temperature,
         .regions = {{"bulk", *physics::Semiconductor::create(p)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, ND),
         .acceptors = std::vector<double>(n, NA),
         .contacts = {{"a", device::ContactKind::ohmic, std::move(a)},
                      {"c", device::ContactKind::ohmic, std::move(c)}}});
}

}  // namespace

TEST_CASE("heterojunction: a graded composition converges with its staircase") {
    // Undeclared steps are drift-diffusion, so refining the staircase converges to the graded
    // device (the weakness Unit 15 fixed: with emission on every step a 10-step grade cut the
    // current to 0.554 and the answer depended on the step count).
    const double j10 = current(graded(10), 0.1), j20 = current(graded(20), 0.1),
                 j40 = current(graded(40), 0.1);
    UNSCOPED_INFO("J(0.1 V) with 10, 20, 40 steps: " << j10 << ", " << j20 << ", " << j40);
    REQUIRE(std::abs(j20 / j40 - 1.0) < std::abs(j10 / j40 - 1.0));
    REQUIRE(std::abs(j20 / j40 - 1.0) < 2e-3);
}

TEST_CASE("heterojunction: incomplete ionization in the device solves") {
    // Boron in silicon at 1e16, Fermi-Dirac, no narrowing: the solved bulk's ionized fraction
    // (p - n) / N_A equals the neutral root bisected in double-double arithmetic
    // (../physics/references.hpp) at 77 and 300 K; 4H-SiC nitrogen and aluminium at 1e17, 300 K.
    struct Case {
        physics::SemiconductorParameters p;
        double ND, NA, T;
        bool fd;
    };
    const Case cases[] = {
        {physics::silicon_parameters, 0.0, 1e16, 77.0, true},
        {physics::silicon_parameters, 0.0, 1e16, 300.0, true},
        {physics::silicon_carbide_4h_parameters, 1e17, 0.0, 300.0, false},
        {physics::silicon_carbide_4h_parameters, 0.0, 1e17, 300.0, false},
    };
    for (const Case& c : cases) {
        CAPTURE(c.T, c.ND, c.NA);
        const device::Device d = uniform(c.p, c.ND, c.NA, c.T);
        const assemble::PhysicsModels models{
            .srh = false, .bgn = false, .fermi_dirac = c.fd, .incomplete_ionization = true};
        const auto eq = solve::solve_equilibrium(
            d, {.newton = {.max_iterations = 200, .tol_update = 1e-12}, .models = models});
        REQUIRE(eq.has_value());
        const double N = c.ND > 0.0 ? c.ND : c.NA;
        const double carriers = c.ND > 0.0 ? eq->fields.n_cm3[10] - eq->fields.p_cm3[10]
                                           : eq->fields.p_cm3[10] - eq->fields.n_cm3[10];
        const physics::Semiconductor m = *physics::Semiconductor::create(c.p);
        const double ni = physics::intrinsic_density(m, c.T);
        const double VT = base::thermal_voltage(c.T);
        const reference::Neutral r = reference::neutral_equilibrium(
            {.donors = c.ND,
             .acceptors = c.NA,
             .donor_kT = c.p.ionization.donor_eV / VT,
             .acceptor_kT = c.p.ionization.acceptor_eV / VT,
             .donor_degeneracy = c.p.ionization.donor_degeneracy,
             .acceptor_degeneracy = c.p.ionization.acceptor_degeneracy},
            ni, std::log(physics::conduction_band_dos(m, c.T) / ni),
            std::log(physics::valence_band_dos(m, c.T) / ni), c.fd);
        const double fraction =
            ((c.ND > 0.0 ? r.n - r.p : r.p - r.n) / reference::DD(N)).value();
        CAPTURE(carriers / N, fraction);
        REQUIRE(close(carriers / N, fraction, 1e-10));
    }
}

TEST_CASE("heterojunction: a 4H-SiC p-n diode with incomplete ionization") {
    // Forward sweep to 3 V with aluminium 10% ionized: it converges, the band diagram's Fermi
    // levels sit at the contacts' biases, and the reported current resolution separates the
    // unresolved low-bias currents from the resolved ones.
    const physics::SemiconductorParameters sic = physics::silicon_carbide_4h_parameters;
    const device::Device d = junction(sic, sic, -1e17, 1e17);
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 6; ++k) points.push_back({0.5 * k, 0.0});
    const auto sweep = solve::sweep_bias(d, points, {.models = {.incomplete_ionization = true}});
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    for (const results::BiasPoint& p : sweep->points) {
        CAPTURE(p.bias_V[0]);
        // Contact nodes: hole Fermi level at the anode is -V, electron level at the cathode 0.
        REQUIRE(std::abs(p.bands.hole_fermi_eV.front() + p.bias_V[0]) < 1e-9);
        REQUIRE(std::abs(p.bands.electron_fermi_eV.back()) < 1e-9);
        UNSCOPED_INFO("V " << p.bias_V[0] << ": J " << p.terminal_current[0] << ", resolution "
                           << p.terminal_current_resolution[0]);
    }
    REQUIRE(std::abs(sweep->points[2].terminal_current[0]) <
            sweep->points[2].terminal_current_resolution[0]);  // 1 V: unresolved
    REQUIRE(sweep->points[6].terminal_current[0] >
            1e3 * sweep->points[6].terminal_current_resolution[0]);  // 3 V: resolved
}

TEST_CASE("heterojunction: radiative recombination in GaAs, none in silicon") {
    // A long-base GaAs p-n diode (200 um a side, 1e17 / 1e17) with radiative recombination alone
    // (SRH and Auger off): in low injection a minority carrier recombines with the lifetime
    // tau = 1 / (B N) (R = B (n p - n_i^2) = B N dn), and with both sides many diffusion lengths
    // long (L_n about 13 um, L_p about 2.4 um) the current is the long-base ideal-diode current
    //     J = q n_i^2 (sqrt(D_n / tau_n) / N_A + sqrt(D_p / tau_p) / N_D) (e^(V / V_T) - 1),
    // D = mu V_T with each side's Caughey-Thomas mobility. Silicon (B = 0) is bit-identical with
    // the model on or off.
    const physics::SemiconductorParameters gaas = physics::gallium_arsenide_parameters;
    const double N = 1e17, side = 2e-2;
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2.0 * side, side, 1e-7, 2e-5));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < side ? acceptors : donors)[i] = N;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"gaas", *physics::Semiconductor::create(gaas)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    const physics::Semiconductor g = *physics::Semiconductor::create(gaas);
    const double VT = base::thermal_voltage(T), ni = physics::intrinsic_density(g, T);
    const double Dn = physics::caughey_thomas_mobility(g, physics::Carrier::electron, N, T) * VT;
    const double Dp = physics::caughey_thomas_mobility(g, physics::Carrier::hole, N, T) * VT;
    const double tau = 1.0 / (gaas.radiative_cm3_s * N);
    CAPTURE(std::sqrt(Dn * tau), std::sqrt(Dp * tau));
    REQUIRE(side > 8.0 * std::sqrt(Dn * tau));  // long base on both sides
    const assemble::PhysicsModels radiative_only{.srh = false, .auger = false};
    double worst = 0.0;
    for (const double V : {0.6, 0.8, 1.0}) {
        const double J = current(d, V, radiative_only);
        const double ideal = base::q_C * ni * ni *
                             (std::sqrt(Dn / tau) / N + std::sqrt(Dp / tau) / N) *
                             std::expm1(V / VT);
        UNSCOPED_INFO("V " << V << ": J " << J << ", long-base ideal " << ideal);
        worst = std::max(worst, std::abs(J / ideal - 1.0));
    }
    CAPTURE(worst);
    REQUIRE(worst < 0.02);
    // Without any recombination the minority carriers diffuse to the contacts instead: the
    // short-base current q n_i^2 (D_n + D_p) / (N W) (e^(V / V_T) - 1), W = side, about L / W of
    // the radiative one.
    const double none = current(d, 0.8, {.srh = false, .auger = false, .radiative = false});
    const double short_base = base::q_C * ni * ni * (Dn + Dp) / (N * side) * std::expm1(0.8 / VT);
    CAPTURE(none, short_base);
    REQUIRE(std::abs(none / short_base - 1.0) < 0.02);
    const device::Device s =
        junction(physics::silicon_parameters, physics::silicon_parameters, -1e17, 1e17);
    REQUIRE(current(s, 0.5, {}) == current(s, 0.5, {.radiative = false}));
}

TEST_CASE("heterojunction: the band diagram at equilibrium and under bias") {
    // GaAs | Al0.3Ga0.7As p-n: at equilibrium both Fermi levels are 0 everywhere and the band step
    // at the interface is the affinity step plus the potential drop across the edge; at 1.2 V the
    // contacts hold the Fermi levels at -1.2 and 0 eV.
    const physics::SemiconductorParameters gaas = physics::gallium_arsenide_parameters;
    const physics::SemiconductorParameters algaas = *physics::algaas_parameters(0.3);
    const device::Device d = junction(gaas, algaas, -1e17, 1e17);
    const auto eq = solve::solve_equilibrium(d, {.models = m33});
    REQUIRE(eq.has_value());
    const auto x = m33_axis();
    const std::size_t k = interface_node(x);
    double worst = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        worst = std::max({worst, std::abs(eq->bands.electron_fermi_eV[i]),
                          std::abs(eq->bands.hole_fermi_eV[i])});
    }
    REQUIRE(worst < 1e-9);
    const auto& phi = eq->fields.potential_V;
    const double step = eq->bands.conduction_eV[k + 1] - eq->bands.conduction_eV[k];
    REQUIRE(std::abs(step - (-(algaas.electron_affinity_eV - gaas.electron_affinity_eV) -
                             (phi[k + 1] - phi[k]))) < 1e-12);
    const std::vector<double> v{1.2, 0.0};
    const auto b = solve::solve_bias(d, v, {.models = m33});
    REQUIRE(b.has_value());
    REQUIRE(std::abs(b->bands.hole_fermi_eV.front() + 1.2) < 1e-9);
    REQUIRE(std::abs(b->bands.electron_fermi_eV.back()) < 1e-9);
    // The gap is Eg everywhere (no narrowing at 1e17 in these sets).
    REQUIRE(std::abs(b->bands.conduction_eV[5] - b->bands.valence_eV[5] -
                     physics::band_gap_eV(*physics::Semiconductor::create(gaas), T)) < 1e-12);
}

TEST_CASE("heterojunction: the current resolution bounds the edge-to-edge spread") {
    // Unit 11's finding: in a 1e19 / 1e18 silicon diode the total current varies edge to edge by
    // rounding. The reported resolution covers that spread and is within 1e4 of it.
    const physics::SemiconductorParameters si = physics::silicon_parameters;
    const device::Device d = junction(si, si, -1e19, 1e18);
    const std::vector<double> v{0.5, 0.0};
    const auto p = solve::solve_bias(d, v);
    REQUIRE(p.has_value());
    double lo = 1e300, hi = -1e300;
    for (std::size_t e = 0; e < p->edge_current_n.size(); ++e) {
        const double J = p->edge_current_n[e] + p->edge_current_p[e];
        lo = std::min(lo, J);
        hi = std::max(hi, J);
    }
    const double spread = hi - lo, resolution = p->terminal_current_resolution[0];
    UNSCOPED_INFO("spread " << spread << " A/cm^2, resolution " << resolution);
    REQUIRE(spread <= resolution);
    REQUIRE(resolution < 1e4 * spread);
}

TEST_CASE("heterojunction: a non-planar interface in 2D solves with a flat Fermi level") {
    // GaAs on x < 1 um for y < 0.5 um and on x < 1.2 um above (a stepped interface), p-n at
    // x = 1 um, the interface declared thermionic: flat Fermi levels at equilibrium and a
    // converging forward sweep.
    const auto x = m33_axis();
    const std::vector<double> y{0.0, 0.25e-4, 0.5e-4, 0.75e-4, 1e-4};
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& pt = m.points()[i];
        region[i] = pt[0] < (pt[1] < 0.5e-4 ? 1e-4 : 1.2e-4) ? 0 : 1;
        (pt[0] < 1e-4 ? acceptors[i] : donors[i]) = 1e17;
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"gaas",
                      *physics::Semiconductor::create(physics::gallium_arsenide_parameters)},
                     {"algaas", *physics::Semiconductor::create(*physics::algaas_parameters(0.3))}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .interfaces = {{"gaas", "algaas", thermionic}}});
    const auto eq = solve::solve_equilibrium(d);
    REQUIRE(eq.has_value());
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        worst = std::max({worst, std::abs(eq->bands.electron_fermi_eV[i]),
                          std::abs(eq->bands.hole_fermi_eV[i])});
    }
    REQUIRE(worst < 1e-9);
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 6; ++k) points.push_back({0.2 * k, 0.0});
    const auto sweep = solve::sweep_bias(d, points);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    const auto& last = sweep->points.back();
    REQUIRE(last.terminal_current[0] > 0.0);
    REQUIRE(std::abs(last.terminal_current[0] + last.terminal_current[1]) <
            1e-6 * last.terminal_current[0]);
}

TEST_CASE("heterojunction: the run identity follows the interface transport") {
    // Drift-diffusion reads each interface's transport; the quasi-static sweep has no continuity
    // equations and ignores it.
    const physics::SemiconductorParameters gaas = physics::gallium_arsenide_parameters;
    const physics::SemiconductorParameters algaas = *physics::algaas_parameters(0.3);
    const device::Device dd = junction(gaas, algaas, -1e17, 1e17);
    const device::Device te =
        junction(gaas, algaas, -1e17, 1e17, {}, device::InterfaceTransport::thermionic_emission);
    const std::vector<std::vector<double>> points{{0.0, 0.0}};
    REQUIRE(solve::make_run_record(dd, {}, points).input_identity !=
            solve::make_run_record(te, {}, points).input_identity);
    solve::BiasOptions quasi_static;
    quasi_static.equations = solve::Equations::equilibrium_poisson;
    REQUIRE(solve::make_run_record(dd, quasi_static, points).input_identity ==
            solve::make_run_record(te, quasi_static, points).input_identity);
}
