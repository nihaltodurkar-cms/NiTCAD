// Fermi-Dirac statistics in the device solves (ARCHITECTURE.md section 11, Unit 14; legacy
// tests/test_m13_solver.py: the scheme's detailed balance, G4 neutrality and built-in potential,
// G4(c) generalized mass action, G6(b) Boltzmann equivalence, G7(a) degenerate bulk, G7(d) the
// degenerate MOS C_max direction).
//
// References: Fermi-Dirac neutral roots bisected in double-double arithmetic (silicon,
// 300 K, no band-gap narrowing, this code's n_i, Nc and Nv; ../physics/references.hpp).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
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
#include "legacy_moscap.hpp"
#include "../physics/references.hpp"

using namespace NiTCAD;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

constexpr double T = 300.0;

// N_A on x < 1 um, N_D beyond, on the legacy junction mesh; anode on x_min, cathode on x_max.
device::Device diode(double NA, double ND) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
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
         .temperature_K = T,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

// Uniform N_D on the given mesh, contacts on x_min and x_max.
device::Device uniform(mesh::Mesh m, double ND) {
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, ND),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}}});
}

std::vector<double> uniform_axis(double length, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = length * i / (nodes - 1);
    return a;
}

double Nc() { return physics::conduction_band_dos(physics::silicon(), T); }
double Nv() { return physics::valence_band_dos(physics::silicon(), T); }
double ni() { return physics::intrinsic_density(physics::silicon(), T); }

// The double-double Fermi-Dirac neutral root at net doping C (../physics/references.hpp), with
// this code's n_i, Nc and Nv.
reference::Neutral fd_root(double C) {
    reference::Dopants d;
    (C > 0.0 ? d.donors : d.acceptors) = std::abs(C);
    return reference::neutral_equilibrium(d, ni(), std::log(Nc() / ni()), std::log(Nv() / ni()),
                                          true);
}

// The reduced Fermi energies of a density: eta = F_{1/2}^-1(n / N).
double eta_of(double density, double N) { return physics::inverse_fermi_half(density / N); }

}  // namespace

TEST_CASE("fermi-dirac: a uniform degenerate bulk holds the neutral root (legacy G7(a))") {
    // N_D = 1e20, no band-gap narrowing: n = N_D + p everywhere and the Fermi level 2.43 kT above
    // the conduction band edge (double-double root); the minority density is n_i^2 gamma_n / N_D, a
    // third of the Boltzmann value.
    const assemble::PhysicsModels models{.srh = false, .bgn = false, .fermi_dirac = true};
    const reference::Neutral root = fd_root(1e20);
    const double p_ref = root.p.value();
    for (const int D : {1, 2}) {
        CAPTURE(D);
        const auto x = uniform_axis(1e-4, 21);
        const device::Device d =
            D == 1 ? uniform(*mesh::make_tensor_grid(x), 1e20)
                   : uniform(*mesh::make_tensor_grid(x, uniform_axis(2e-5, 3)), 1e20);
        const auto eq = solve::solve_equilibrium(d, {.models = models});
        REQUIRE(eq.has_value());
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            CAPTURE(i);
            REQUIRE(close(eq->fields.n_cm3[i], 1e20, 1e-13));
            REQUIRE(close(eq->fields.p_cm3[i], p_ref, 1e-11));
        }
        REQUIRE(std::abs(eta_of(eq->fields.n_cm3[0], Nc()) -
                         (root.eta - std::log(Nc() / ni()))) <= 1e-8);
        const auto boltzmann =
            solve::solve_equilibrium(d, {.models = {.srh = false, .bgn = false}});
        REQUIRE(close(boltzmann->fields.p_cm3[0], ni() * ni() / 1e20, 1e-12));
        REQUIRE(eq->fields.p_cm3[0] < 0.4 * boltzmann->fields.p_cm3[0]);
    }
}

TEST_CASE("fermi-dirac: built-in potential of a degenerate junction (legacy G4(d))") {
    // p+ 1e20 / n 1e17, no band-gap narrowing: V_bi = V_T (eta(1e17) - eta(-1e20)) from the
    // double-double roots, 28 mV above the Boltzmann value because the degenerate side needs a
    // deeper Fermi level for the same hole density. The legacy gate was 1e-3 V; both ends are Dirichlet
    // contact values, so the solve reproduces them to rounding.
    const device::Device d = diode(1e20, 1e17);
    const auto fd = solve::solve_equilibrium(d, {.models = {.bgn = false, .fermi_dirac = true}});
    const auto bz = solve::solve_equilibrium(d, {.models = {.bgn = false}});
    REQUIRE(fd.has_value());
    REQUIRE(bz.has_value());
    const double VT = base::thermal_voltage(T);
    const double eta_n = fd_root(1e17).eta, eta_p = -fd_root(-1e20).eta;
    const double vbi_ref = VT * (eta_n + eta_p);
    const auto& phi = fd->fields.potential_V;
    const double vbi = phi.back() - phi.front();
    REQUIRE(close(vbi, vbi_ref, 1e-13));
    const double shift = vbi - (bz->fields.potential_V.back() - bz->fields.potential_V.front());
    UNSCOPED_INFO("V_bi " << vbi << " V, above Boltzmann by " << shift << " V");
    const double expected = VT * ((eta_n - std::asinh(1e17 / (2.0 * ni()))) +
                                  (eta_p - std::asinh(1e20 / (2.0 * ni()))));
    REQUIRE(close(shift, expected, 1e-11));
    // The 1e20 side is degenerate (legacy G3 companion: eta_p > 2).
    REQUIRE(eta_of(fd->fields.p_cm3.front(), Nv()) > 2.0);
}

TEST_CASE("fermi-dirac: generalized mass action at every node of a solved junction (G4(c))") {
    // n p e^(eta_n + eta_p) / (F(eta_n) F(eta_p)) = n_i^2, with each eta from its own density:
    // at equilibrium eta_n + eta_p = -E_g / kT, so this is Nc Nv e^(-E_g/kT). No band-gap
    // narrowing (n_i is then one number).
    for (const auto& [NA, ND] : {std::pair{1e17, 1e17}, std::pair{1e20, 1e17}}) {
        CAPTURE(NA, ND);
        const device::Device d = diode(NA, ND);
        const auto eq = solve::solve_equilibrium(
            d, {.newton = {.tol_update = 1e-12}, .models = {.bgn = false, .fermi_dirac = true}});
        REQUIRE(eq.has_value());
        double worst = 0.0;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            const double n = eq->fields.n_cm3[i], p = eq->fields.p_cm3[i];
            const double en = eta_of(n, Nc()), ep = eta_of(p, Nv());
            const double lhs = (n * p) * std::exp(en + ep) /
                               (physics::fermi_half(en).value * physics::fermi_half(ep).value);
            worst = std::max(worst, std::abs(lhs / (ni() * ni()) - 1.0));
        }
        UNSCOPED_INFO("worst mass-action error " << worst);
        REQUIRE(worst <= 1e-10);
    }
}

TEST_CASE("fermi-dirac: nondegenerate devices reproduce Boltzmann statistics (legacy G6(b))") {
    // 1e16 / 1e16: the majority carriers deviate from Boltzmann by delta = e^eta / 2^(3/2) =
    // 1.2e-4 (eta = -8 below the band edge). Gates as the legacy: densities within 3 delta (two
    // degeneracy factors), currents at 0.5 V within 20 delta.
    const device::Device d = diode(1e16, 1e16);
    const double delta = std::exp(eta_of(1e16, Nc())) / std::pow(2.0, 1.5);
    const std::vector<double> zero{0.0, 0.0}, forward{0.5, 0.0};
    const solve::BiasOptions fd{.models = {.fermi_dirac = true}};
    const auto bz0 = solve::solve_bias(d, zero), fd0 = solve::solve_bias(d, zero, fd);
    REQUIRE(bz0.has_value());
    REQUIRE(fd0.has_value());
    double worst = 0.0;
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        worst = std::max({worst, std::abs(fd0->fields.n_cm3[i] / bz0->fields.n_cm3[i] - 1.0),
                          std::abs(fd0->fields.p_cm3[i] / bz0->fields.p_cm3[i] - 1.0)});
    }
    UNSCOPED_INFO("delta " << delta << ", worst density deviation " << worst);
    REQUIRE(worst <= 3.0 * delta);
    const auto bz = solve::solve_bias(d, forward), f = solve::solve_bias(d, forward, fd);
    REQUIRE(bz.has_value());
    REQUIRE(f.has_value());
    double worst_current = 0.0;
    for (std::size_t k = 0; k < bz->edge_current_n.size(); ++k) {
        const double jb = bz->edge_current_n[k] + bz->edge_current_p[k];
        const double jf = f->edge_current_n[k] + f->edge_current_p[k];
        worst_current = std::max(worst_current, std::abs(jf / jb - 1.0));
    }
    UNSCOPED_INFO("worst current deviation " << worst_current);
    REQUIRE(worst_current <= 20.0 * delta);
}

TEST_CASE("fermi-dirac: equilibrium carries no current across a degenerate junction") {
    // The scheme gate (legacy test_sg_scheme_detailed_balance_degenerate_step): the fluxes at the
    // equilibrium state. 1e20 / 1e17 with band-gap narrowing on, so n_ie and both degeneracy
    // factors vary along the device. Every edge is zero to rounding, each held to 64 eps of its
    // one-sided terms, bounded by q mu_max V_T / h (n_a + n_b) (2 + |delta|), delta the driving
    // term with its degeneracy part (B(-x) = B(x) + x).
    const device::Device d = diode(1e20, 1e17);
    const assemble::PhysicsModels models{.fermi_dirac = true};
    const auto eq = solve::solve_equilibrium(d, {.models = models});
    REQUIRE(eq.has_value());
    const auto scaling = *assemble::make_scaling(d);
    const auto system = *assemble::DriftDiffusion::create(d, scaling, models);
    std::vector<double> psi(d.mesh().node_count());
    for (std::size_t i = 0; i < psi.size(); ++i) psi[i] = eq->fields.potential_V[i] / scaling.V_T;
    const std::vector<double> state = system.state_from_potential(psi);
    const auto currents = system.edge_currents(state);
    const double VT = base::thermal_voltage(T);
    const auto x = legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6);
    const physics::Semiconductor si = physics::silicon();
    const auto& f = eq->fields;
    const auto log_gamma = [&](std::size_t i, bool electrons) {
        const double N_total = d.total_impurity(static_cast<mesh::NodeId>(i));
        const double nie = physics::effective_intrinsic_density(si, N_total, T);
        const double density = electrons ? f.n_cm3[i] : f.p_cm3[i];
        return physics::fermi_dirac_degeneracy(nie, std::log((electrons ? Nc() : Nv()) / nie),
                                               density)
            .log_gamma;
    };
    double worst = 0.0;
    for (std::size_t k = 0; k + 1 < x.size(); ++k) {
        const double dpsi = std::abs(psi[k + 1] - psi[k]);
        const double dn = dpsi + 1.0 + std::abs(log_gamma(k + 1, true) - log_gamma(k, true));
        const double dp = dpsi + 1.0 + std::abs(log_gamma(k + 1, false) - log_gamma(k, false));
        const double prefactor = base::q_C * 1360.0 * VT / (x[k + 1] - x[k]);
        const double floor = 64.0 * 2.2204460492503131e-16 * prefactor;
        const double bound_n = floor * (2.0 + dn) * (f.n_cm3[k] + f.n_cm3[k + 1]);
        const double bound_p = floor * (2.0 + dp) * (f.p_cm3[k] + f.p_cm3[k + 1]);
        const double jn = currents[k].first * scaling.J0, jp = currents[k].second * scaling.J0;
        CAPTURE(k);
        REQUIRE(std::abs(jn) <= bound_n);
        REQUIRE(std::abs(jp) <= bound_p);
        worst = std::max({worst, std::abs(jn) / bound_n, std::abs(jp) / bound_p});
    }
    UNSCOPED_INFO("worst edge current over its rounding bound " << worst);
    // The bias solve at zero bias starts there and stops after one correction: the equilibrium
    // potential leaves a few ulp of charge per node (a potential of 22 V_T resolves e^psi only to
    // 22 eps), which that step removes.
    const std::vector<double> zero{0.0, 0.0};
    const auto s = solve::solve_bias(d, zero, {.models = models});
    REQUIRE(s.has_value());
    REQUIRE(s->convergence.iterations.size() == 1);
    REQUIRE(s->convergence.iterations.front().update < 1e-9);
}

TEST_CASE("fermi-dirac: forward sweep of a degenerate p+ n diode") {
    // p+ 1e20 / n 1e17, all default models (band-gap narrowing, SRH, Auger) on, 0 to 0.8 V.
    // Fermi-Dirac statistics lower the minority electron density of the degenerate p+ side by its
    // hole degeneracy factor (n_0 = n_ie^2 gamma_n gamma_p / p_0, gamma_n = 1), so the electron
    // current injected into it falls by gamma_p = 0.335 at every low-injection bias, while the
    // hole current injected into the nondegenerate n side stays within that side's own degeneracy
    // correction (e^eta / 2^(3/2) = 1.2e-3). The diode current falls by 1.7-2.4%, the electron
    // share of it. Measured at 0.2-0.7 V: electron ratio 0.33538-0.33545, hole ratio
    // 0.99933-0.99936.
    const device::Device d = diode(1e20, 1e17);
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 8; ++k) points.push_back({0.1 * k, 0.0});
    const auto fd = solve::sweep_bias(d, points, {.models = {.fermi_dirac = true}});
    const auto bz = solve::sweep_bias(d, points);
    REQUIRE(fd.has_value());
    REQUIRE(bz.has_value());
    REQUIRE_FALSE(fd->stopped.has_value());
    REQUIRE_FALSE(bz->stopped.has_value());
    const physics::Semiconductor si = physics::silicon();
    const double nie = physics::effective_intrinsic_density(si, 1e20, T);
    const double gamma_p =
        std::exp(physics::fermi_dirac_degeneracy(nie, std::log(Nv() / nie), 1e20).log_gamma);
    const double delta_n = std::exp(eta_of(1e17, Nc())) / std::pow(2.0, 1.5);
    for (std::size_t k = 1; k < points.size(); ++k) {
        CAPTURE(points[k][0]);
        const auto& f = fd->points[k];
        const auto& b = bz->points[k];
        REQUIRE(f.terminal_current[0] > fd->points[k - 1].terminal_current[0]);
        REQUIRE(f.convergence.iterations.size() <= 8);
        // Below 0.2 V the current (1e-12 A/cm^2 at 0.1 V) is rounding noise in a 1e20 device
        // (ARCHITECTURE.md 5, Unit 11); above 0.7 V high injection sets in.
        if (k < 2 || k > 7) continue;
        const double electrons = f.edge_current_n.front() / b.edge_current_n.front();
        const double holes = f.edge_current_p.back() / b.edge_current_p.back();
        UNSCOPED_INFO("V " << points[k][0] << ": electron ratio " << electrons << " (gamma_p "
                           << gamma_p << "), hole ratio " << holes << ", diode ratio "
                           << f.terminal_current[0] / b.terminal_current[0]);
        REQUIRE(std::abs(electrons / gamma_p - 1.0) <= 2e-3);
        REQUIRE(std::abs(holes - 1.0) <= delta_n);
        REQUIRE(f.terminal_current[0] < b.terminal_current[0]);
    }
    // The current extracted at both contacts agrees (continuity), as with Boltzmann.
    const auto& last = fd->points.back();
    REQUIRE(std::abs(last.terminal_current[0] + last.terminal_current[1]) <=
            1e-8 * std::abs(last.terminal_current[0]));
}

TEST_CASE("fermi-dirac: the degenerate MOS C_max sits below the Boltzmann value (legacy G7(d))") {
    // Legacy fixture: N_A = 1e18, t_ox = 5 nm, n+ polysilicon gate, quasi-static C-V from -2.5 to
    // 3 V. With Fermi-Dirac statistics the accumulation and inversion layers' density grows like
    // eta^(3/2) instead of e^eta, so their differential capacitance is finite and C_max falls below
    // the classical value, which is within 3% of C_ox. Legacy gate: a reduction of 2-30%.
    mesh::Mesh m = *mesh::make_tensor_grid(LegacyMOSCapacitor(-1e18, 5e-7).x);
    const std::size_t n = m.node_count();
    auto gate_nodes = m.find_boundary("x_min")->nodes;
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Contact gate{"gate", device::ContactKind::gate, std::move(gate_nodes)};
    gate.gate = {.boundary = "x_min", .oxide_thickness_cm = 5e-7};
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::vector<double>(n, 1e18),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
    std::vector<double> Vg;
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 110; ++k) {
        Vg.push_back(-2.5 + 0.05 * k);
        points.push_back({Vg.back(), 0.0});
    }
    const auto c_max = [&](bool fermi_dirac) {
        solve::BiasOptions o{.models = {.fermi_dirac = fermi_dirac},
                             .equations = solve::Equations::equilibrium_poisson};
        const auto sweep = solve::sweep_bias(d, points, o);
        REQUIRE(sweep.has_value());
        REQUIRE_FALSE(sweep->stopped.has_value());
        std::vector<double> Q;
        for (const auto& p : sweep->points) Q.push_back(p.gate_charge[0]);
        const std::vector<double> C = numpy_gradient(Q, Vg);
        return *std::max_element(C.begin(), C.end());
    };
    const double cox = legacy_eps_ox_r * base::eps0_F_per_cm / 5e-7;
    const double classical = c_max(false), degenerate = c_max(true);
    const double reduction = 1.0 - degenerate / classical;
    UNSCOPED_INFO("C_max / C_ox: Boltzmann " << classical / cox << ", Fermi-Dirac "
                                             << degenerate / cox << ", reduction " << reduction);
    REQUIRE(classical >= 0.97 * cox);
    REQUIRE(reduction >= 0.02);
    REQUIRE(reduction <= 0.30);
}

TEST_CASE("fermi-dirac: the run identity includes the statistics") {
    const device::Device d = diode(1e17, 1e17);
    const std::vector<std::vector<double>> points{{0.0, 0.0}};
    for (const auto equations :
         {solve::Equations::drift_diffusion, solve::Equations::equilibrium_poisson}) {
        solve::BiasOptions o{.equations = equations};
        const auto base = solve::make_run_record(d, o, points);
        o.models.fermi_dirac = true;
        const auto record = solve::make_run_record(d, o, points);
        REQUIRE(record.input_identity != base.input_identity);
        bool found = false;
        for (const auto& [name, value] : record.settings) {
            if (name == "models.fermi_dirac") found = value == 1.0;
        }
        REQUIRE(found);
    }
}

