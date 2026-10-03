// Heterojunctions in the device solves (ARCHITECTURE.md section 11, Unit 15; legacy
// tests/test_m33_interface.py: G1 detailed balance per carrier, G2 the affinity step's effect and
// direction, the isotype barrier, S2 thermionic emission; the exact equilibrium of an abrupt
// heterojunction; composition with Fermi-Dirac statistics and with a gate).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
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
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

constexpr double T = 300.0;

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// The legacy M33 mesh: graded_mesh(2e-4, [1e-4], 1e-8, 1e-6, 1.12).
std::vector<double> m33_axis() { return legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6, 1.12); }

// `left` on x < 1 um with net doping C_left, `right` beyond with C_right (positive: donors), ohmic
// contacts on x_min and x_max. With `transverse`, a y-uniform 2D device of that width.
device::Device junction(const physics::SemiconductorParameters& left,
                        const physics::SemiconductorParameters& right, double C_left,
                        double C_right, const std::vector<double>& transverse = {}) {
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
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
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
                junction(physics::silicon_parameters, silicon_with(c.chi, c.Eg0), -1e17, 1e17);
            assemble::PhysicsModels models = m33;
            models.thermionic_emission = te;
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
    // Legacy S2 fixture: n-Si / n-Si(chi), 1e17, 0.1 V. At the real emission velocity the interface
    // cuts the current relative to drift-diffusion, and the cut deepens with the step. (At no step
    // the two materials are one, so there is no interface edge and the ratio is exactly 1; the
    // legacy, which keyed interfaces on the material object, applied emission there too.)
    const auto ratio = [&](double chi) {
        const device::Device d = junction(physics::silicon_parameters, silicon_with(chi), 1e17,
                                          1e17);
        assemble::PhysicsModels te = m33;
        te.thermionic_emission = true;
        return current(d, 0.1, te) / current(d, 0.1);
    };
    const double r0 = ratio(4.05), r1 = ratio(4.20), r2 = ratio(4.35);
    UNSCOPED_INFO("J(TE) / J(DD): " << r0 << ", " << r1 << ", " << r2);
    REQUIRE(r0 == 1.0);
    REQUIRE(r1 < 1.0);
    REQUIRE(r2 < r1);
    REQUIRE(r2 > 0.0);
}

TEST_CASE("heterojunction: composes with Fermi-Dirac statistics and band-gap narrowing") {
    // The legacy refused band_offset="affinity" with fd: NiTCAD's shift only moves eta, so the
    // composition is exact. A degenerate n+ GaAs (1e19, eta_c > 0) against p-Al0.3Ga0.7As, with
    // thermionic emission: flat Fermi level at equilibrium and a converging forward sweep. The
    // AlGaAs minority electrons are 2e-11 cm^-3, 2e-30 of the largest density: the sweep converges
    // because the Newton measure floors densities at 1e-20 of the largest (6.2, Unit 15).
    const device::Device d = junction(physics::gallium_arsenide_parameters,
                                      *physics::algaas_parameters(0.3), 1e19, -1e17);
    const assemble::PhysicsModels models{.fermi_dirac = true, .thermionic_emission = true};
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
    const auto d1 = junction(physics::gallium_arsenide_parameters, right, -1e17, 1e17);
    const auto d2 = junction(physics::gallium_arsenide_parameters, right, -1e17, 1e17, y);
    assemble::PhysicsModels models{.thermionic_emission = true};
    const double j1 = current(d1, 1.2, models);
    const double j2 = current(d2, 1.2, models);
    UNSCOPED_INFO("1D " << j1 << " A/cm^2, 2D " << j2 / 1e-4 << " A/cm^2 per width");
    REQUIRE(close(j2 / 1e-4, j1, 1e-9));
}

