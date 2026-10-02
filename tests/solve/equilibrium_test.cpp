// Equilibrium solve of a silicon pn diode (ARCHITECTURE.md section 11, Unit 8 gates: built-in
// potential within 2e-3 V; bulk neutrality; convergence; non-convergence returns an error).
//
// The legacy built-in-potential gate (test_device1d_native_gates.py::test_g1_built_in_potential)
// takes psi[-1] - psi[0]; both ends are Dirichlet contacts, so it checks the contact values, not
// the solve. The solve itself is checked by the peak field against the depletion approximation,
// global and bulk neutrality, quadratic convergence and the 2D/3D reduction to 1D.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/newton.hpp"

using namespace NiTCAD;
using base::ErrorCode;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// The legacy fixture's mesh, graded_mesh(2e-4, [1e-4], h_min=1e-8, h_max=1e-6): spacing 1e-8 at
// the junction, growing by 1.2 per node up to 1e-6.
std::vector<double> junction_axis(double h_min = 1e-8, double h_max = 1e-6) {
    std::vector<double> d{0.0};
    double h = h_min;
    while (d.back() + h < 1e-4) {
        d.push_back(d.back() + h);
        h = std::min(h * 1.2, h_max);
    }
    if (1e-4 - d.back() < 0.5 * h) d.pop_back();
    d.push_back(1e-4);
    std::vector<double> x;
    for (auto it = d.rbegin(); it != d.rend(); ++it) x.push_back(1e-4 - *it);
    for (std::size_t k = 1; k < d.size(); ++k) x.push_back(1e-4 + d[k]);
    return x;
}

std::vector<double> uniform_axis(double length, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = length * i / (nodes - 1);
    return a;
}

device::Device diode(mesh::Mesh m, double NA = 1e17, double ND = 1e17, double T = 300.0) {
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

double built_in(double NA, double ND, double T) {
    const double ni = physics::intrinsic_density(physics::silicon(), T);
    return base::thermal_voltage(T) * std::log(NA * ND / (ni * ni));
}

}  // namespace

TEST_CASE("equilibrium: built-in potential within 2e-3 V (legacy gate)") {
    for (const double T : {300.0, 400.0}) {
        const auto x = junction_axis();
        const auto d = diode(*mesh::make_tensor_grid(x), 1e17, 1e17, T);
        const auto s = solve::solve_equilibrium(d);
        REQUIRE(s.has_value());
        const double vbi = s->potential_V.back() - s->potential_V.front();
        CAPTURE(T, vbi);
        REQUIRE(std::abs(vbi - built_in(1e17, 1e17, T)) < 2e-3);
    }
}

TEST_CASE("equilibrium: converges quadratically, with one analysis") {
    const auto d = diode(*mesh::make_tensor_grid(junction_axis()));
    const auto scaling = *assemble::make_scaling(d);
    const auto system = *assemble::EquilibriumPoisson::create(d, scaling);
    auto solver = *linalg::LinearSolver::create({});
    std::vector<double> psi = system.charge_neutral_potential();
    const auto report = solve::newton_solve(system, psi, {}, solver);
    REQUIRE(report.has_value());
    UNSCOPED_INFO("iterations " << report->iterations);
    REQUIRE(report->iterations <= 30);
    REQUIRE(solver.analyses() == 1);
    // The last three corrections shrink quadratically.
    const auto& u = report->updates;
    REQUIRE(u.size() >= 3);
    for (std::size_t k = u.size() - 2; k < u.size(); ++k) {
        CAPTURE(k, u[k - 1], u[k]);
        REQUIRE(u[k] <= 10.0 * u[k - 1] * u[k - 1]);
    }
    std::vector<double> f(psi.size());
    system.residual(psi, f);
    double worst = 0.0;
    for (const double v : f) worst = std::max(worst, std::abs(v));
    REQUIRE(worst <= 1e-10);
}

TEST_CASE("equilibrium: bulk and global neutrality") {
    const auto x = junction_axis();
    const auto d = diode(*mesh::make_tensor_grid(x));
    const auto s = *solve::solve_equilibrium(d);
    double net = 0.0, depleted = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double C = d.net_doping(static_cast<mesh::NodeId>(i));
        const double rho = s.p_cm3[i] - s.n_cm3[i] + C;  // charge density / q
        // More than 0.4 um from the junction (the depletion region is about 0.15 um wide).
        if (std::abs(x[i] - 1e-4) > 0.4e-4) {
            CAPTURE(x[i]);
            REQUIRE(std::abs(rho) <= 1e-6 * std::abs(C));
        }
        const double v = d.mesh().volumes()[i];
        net += v * rho;
        depleted += v * std::abs(rho);
    }
    // The field vanishes at both neutral contacts, so the total charge is zero (Gauss).
    REQUIRE(std::abs(net) <= 1e-8 * depleted);
    // n p = n_i^2 everywhere (carriers are slaved to the potential at equilibrium).
    const double ni = s.scaling.n_i;
    for (std::size_t i = 0; i < x.size(); ++i) {
        REQUIRE(close(s.n_cm3[i] * s.p_cm3[i], ni * ni, 1e-13));
    }
}

TEST_CASE("equilibrium: peak field matches the depletion approximation") {
    // Symmetric abrupt junctions, depletion approximation with the 2 V_T correction for the
    // majority-carrier tails: E_max = sqrt(2 q (V_bi - 2 V_T) N_A N_D / (eps (N_A + N_D))).
    // Measured: 0.01% (1e16), 0.10% (1e17), 0.17% (1e18); without the correction about 3%.
    const double VT = base::thermal_voltage(300.0);
    const double eps = 11.7 * base::eps0_F_per_cm;
    const auto peak_field = [](const std::vector<double>& x, const std::vector<double>& phi) {
        double peak = 0.0;
        for (std::size_t k = 0; k + 1 < x.size(); ++k) {
            peak = std::max(peak, std::abs(phi[k + 1] - phi[k]) / (x[k + 1] - x[k]));
        }
        return peak;
    };
    const auto depletion = [&](double NA, double ND, double v) {
        return std::sqrt(2.0 * base::q_C * v * NA * ND / (eps * (NA + ND)));
    };
    const auto x = junction_axis();
    for (const double N : {1e16, 1e17, 1e18}) {
        const auto s = *solve::solve_equilibrium(diode(*mesh::make_tensor_grid(x), N, N));
        const double peak = peak_field(x, s.potential_V);
        const double vbi = built_in(N, N, 300.0);
        CAPTURE(N, peak);
        REQUIRE(close(peak, depletion(N, N, vbi - 2.0 * VT), 0.005));
        REQUIRE_FALSE(close(peak, depletion(N, N, vbi), 0.02));  // the correction matters
    }
    // Asymmetric junctions are outside the approximation: holes from the 1e18 side spill a few
    // Debye lengths (4 nm) into the 1e16 side, about 1e11 cm^-2 against a depletion charge of
    // 3e11 cm^-2, and raise the field at the junction well above it (measured +39%).
    const auto s = *solve::solve_equilibrium(diode(*mesh::make_tensor_grid(x), 1e18, 1e16));
    const double approx = depletion(1e18, 1e16, built_in(1e18, 1e16, 300.0) - 2.0 * VT);
    REQUIRE(peak_field(x, s.potential_V) > 1.2 * approx);
}

TEST_CASE("equilibrium: a y-uniform 2D and 3D device reproduce 1D") {
    const auto x = junction_axis(5e-8, 2e-6);
    const auto s1 = *solve::solve_equilibrium(diode(*mesh::make_tensor_grid(x)));
    const auto y = uniform_axis(1e-4, 4);
    const auto z = uniform_axis(5e-5, 3);
    const auto s2 = *solve::solve_equilibrium(diode(*mesh::make_tensor_grid(x, y)));
    const auto s3 = *solve::solve_equilibrium(diode(*mesh::make_tensor_grid(x, y, z)));
    double worst2 = 0.0, worst3 = 0.0;
    for (std::size_t k = 0; k < s2.potential_V.size(); ++k) {
        worst2 = std::max(worst2, std::abs(s2.potential_V[k] - s1.potential_V[k % x.size()]));
    }
    for (std::size_t k = 0; k < s3.potential_V.size(); ++k) {
        worst3 = std::max(worst3, std::abs(s3.potential_V[k] - s1.potential_V[k % x.size()]));
    }
    UNSCOPED_INFO("2D " << worst2 << " V, 3D " << worst3 << " V");
    REQUIRE(worst2 <= 1e-9);  // legacy 2D gate: y-independent to atol 1e-9
    REQUIRE(worst3 <= 1e-9);
}

TEST_CASE("equilibrium: an undoped device has zero potential") {
    const auto s =
        *solve::solve_equilibrium(diode(*mesh::make_tensor_grid(junction_axis()), 0.0, 0.0));
    for (const double v : s.potential_V) REQUIRE(v == 0.0);
    REQUIRE(s.newton.iterations == 1);
}

TEST_CASE("equilibrium: non-convergence and bad options return errors") {
    const auto d = diode(*mesh::make_tensor_grid(junction_axis()));
    const auto few = solve::solve_equilibrium(d, {.newton = {.max_iterations = 2}});
    REQUIRE_FALSE(few.has_value());
    REQUIRE(few.error().code == ErrorCode::non_convergence);
    REQUIRE(few.error().context->index == 2);
    REQUIRE(solve::solve_equilibrium(d, {.Ns_override = -1.0}).error().code ==
            ErrorCode::invalid_input);
    REQUIRE(solve::solve_equilibrium(d, {.linear = {.threads = 4}}).error().code ==
            ErrorCode::invalid_input);
}

