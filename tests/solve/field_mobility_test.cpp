// Field-dependent (Canali) mobility at device level (ARCHITECTURE.md section 11, Unit 13).
//
// Gates:
// - a uniformly doped resistor with ohmic contacts has the exact discrete solution of a linear
//   potential and uniform densities, so its current is J = q (n mu_n(E) + p mu_p(E)) E with
//   E = V / L to rounding: velocity saturation for electrons (beta = 2) and holes (beta = 1), the
//   low-field limit, and y-uniform 2D and 3D reproducing 1D;
// - Newton keeps its quadratic convergence (the Jacobian includes the mobility's field derivative);
// - the legacy J(0.5 V) diode: the effect of field mobility measured and current continuity kept;
// - the legacy MOSFET (mosfet_test.cpp, Release only, [.mosfet]): the drain current at high drain
//   bias falls below the constant-mobility one, and hardly at low drain bias.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

constexpr double L = 1e-5;  // 100 nm resistor
constexpr double N = 1e16;

std::vector<double> uniform_axis(double length, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = length * i / (nodes - 1);
    return a;
}

// A uniformly doped resistor, contact 0 on x_min (biased), contact 1 on x_max (grounded).
device::Device resistor(mesh::Mesh m, bool n_type) {
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, n_type ? N : 0.0),
         .acceptors = std::vector<double>(n, n_type ? 0.0 : N),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}}});
}

solve::BiasOptions field_mobility() {
    solve::BiasOptions o;
    o.models.field_mobility = true;
    return o;
}

std::vector<std::vector<double>> ramp(double step, int count) {
    std::vector<std::vector<double>> points;
    for (int k = 1; k <= count; ++k) points.push_back({step * k, 0.0});
    return points;
}

results::Sweep sweep(const device::Device& d, const std::vector<std::vector<double>>& points,
                     const solve::BiasOptions& options) {
    auto s = solve::sweep_bias(d, points, options);
    REQUIRE(s.has_value());
    if (s->stopped) {
        CAPTURE(s->stopped->message, s->points.size());
        FAIL("the sweep stopped");
    }
    return std::move(*s);
}

// q (n mu_n(E) + p mu_p(E)) E [A/cm^2] of the uniform resistor at field E.
double resistor_current(bool n_type, double E, bool with_field) {
    const physics::Semiconductor si = physics::silicon();
    const double nie = physics::effective_intrinsic_density(si, N, 300.0);
    const auto eq = physics::boltzmann_neutral_equilibrium(n_type ? N : -N, nie);
    double J = 0.0;
    for (const auto [carrier, density] :
         {std::pair{physics::Carrier::electron, eq.n}, std::pair{physics::Carrier::hole, eq.p}}) {
        const double mu0 = physics::caughey_thomas_mobility(si, carrier, N, 300.0);
        const double mu = with_field
                              ? physics::canali_mobility(mu0, E, physics::saturation(si, carrier))
                                    .mobility
                              : mu0;
        J += base::q_C * density * mu * E;
    }
    return J;
}

}  // namespace

TEST_CASE("field mobility: a uniform resistor saturates exactly as Canali") {
    // 0.1 V steps to 10 V: fields 1e4 to 1e6 V/cm.
    const auto points = ramp(0.1, 100);
    for (const bool n_type : {true, false}) {
        CAPTURE(n_type);
        const auto d = resistor(*mesh::make_tensor_grid(uniform_axis(L, 41)), n_type);
        const auto s = sweep(d, points, field_mobility());
        double worst = 0.0;
        int most_iterations = 0;
        for (const auto& p : s.points) {
            const double E = p.bias_V[0] / L;
            const double J = p.terminal_current[0];
            worst = std::max(worst, std::abs(J / resistor_current(n_type, E, true) - 1.0));
            REQUIRE(std::abs(p.terminal_current[1] + J) <= 1e-12 * J);  // Kirchhoff
            most_iterations = std::max(most_iterations,
                                       static_cast<int>(p.convergence.iterations.size()));
        }
        // Velocity saturation at 1e6 V/cm: electrons (beta = 2, x = mu0 E / v_sat = 117) within
        // 0.1% of v_sat; holes (beta = 1, x = 48.5) at x / (1 + x) = 98.0% of it.
        const physics::Semiconductor si = physics::silicon();
        const auto& sat = physics::saturation(
            si, n_type ? physics::Carrier::electron : physics::Carrier::hole);
        const double v = s.points.back().terminal_current[0] / (base::q_C * N);
        CAPTURE(worst, most_iterations, v / sat.v_sat_cm_s);
        REQUIRE(worst < 1e-10);
        REQUIRE(v < sat.v_sat_cm_s);
        REQUIRE(v > (n_type ? 0.999 : 0.979) * sat.v_sat_cm_s);
        // With constant mobility the resistor is linear in psi (2 iterations a point); the Canali
        // factor makes it nonlinear: measured at most 5 (electrons) and 7 (holes, beta = 1, at
        // the lowest fields), falling to 3.
        REQUIRE(most_iterations <= (n_type ? 5 : 7));
    }
}

TEST_CASE("field mobility: off, the resistor is ohmic; on, the ratio is the Canali factor") {
    const auto d = resistor(*mesh::make_tensor_grid(uniform_axis(L, 41)), true);
    const auto points = ramp(0.1, 10);
    const auto on = sweep(d, points, field_mobility());
    const auto off = sweep(d, points, {});
    for (std::size_t k = 0; k < points.size(); ++k) {
        const double E = points[k][0] / L;
        CAPTURE(E);
        const double J_on = on.points[k].terminal_current[0];
        const double J_off = off.points[k].terminal_current[0];
        REQUIRE(std::abs(J_off / resistor_current(true, E, false) - 1.0) < 1e-10);
        REQUIRE(std::abs(J_on / J_off - resistor_current(true, E, true) /
                                            resistor_current(true, E, false)) < 1e-10);
        REQUIRE(J_on < J_off);
    }
}

TEST_CASE("field mobility: a y-uniform 2D and 3D resistor reproduces 1D") {
    const auto x = uniform_axis(L, 21);
    const std::vector<double> y{0.0, 2e-6, 5e-6}, z{0.0, 1e-6, 4e-6};
    const auto points = ramp(0.25, 12);  // to 3e5 V/cm
    const auto one = sweep(resistor(*mesh::make_tensor_grid(x), true), points, field_mobility());
    const auto two = sweep(resistor(*mesh::make_tensor_grid(x, y), true), points, field_mobility());
    const auto three =
        sweep(resistor(*mesh::make_tensor_grid(x, y, z), true), points, field_mobility());
    for (std::size_t k = 0; k < points.size(); ++k) {
        const double J = one.points[k].terminal_current[0];
        REQUIRE(std::abs(two.points[k].terminal_current[0] / (J * y.back()) - 1.0) < 1e-12);
        REQUIRE(std::abs(three.points[k].terminal_current[0] / (J * y.back() * z.back()) - 1.0) <
                1e-12);
    }
}

TEST_CASE("field mobility: Newton stays quadratic with the exact Jacobian") {
    // One 0.5 V step on the resistor (from 2 to 2.5 V, fields 2e5 to 2.5e5 V/cm): the last
    // corrections shrink quadratically.
    const auto d = resistor(*mesh::make_tensor_grid(uniform_axis(L, 41)), true);
    const std::vector<std::vector<double>> points{{2.0, 0.0}, {2.5, 0.0}};
    const auto s = sweep(d, ramp(0.1, 20), field_mobility());
    auto step = solve::solve_bias(d, points[1], field_mobility(), &s.points.back().fields);
    REQUIRE(step.has_value());
    const auto& it = step->convergence.iterations;
    REQUIRE(it.size() >= 3);
    const double e1 = it[it.size() - 3].update, e2 = it[it.size() - 2].update;
    CAPTURE(it.size(), e1, e2, it.back().update);
    REQUIRE(e2 <= 10.0 * e1 * e1);
}

TEST_CASE("field mobility: the legacy J(0.5 V) diode") {
    // Legacy test_device1d_native_gates _diode fixture: 1e17 / 1e17, graded mesh, 250 nodes.
    const auto x = legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const auto d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    std::vector<std::vector<double>> points;
    for (int k = 1; k <= 5; ++k) points.push_back({0.1 * k, 0.0});
    const auto on = sweep(d, points, field_mobility());
    const auto off = sweep(d, points, {});
    const double J_on = on.points.back().terminal_current[0];
    const double J_off = off.points.back().terminal_current[0];
    double spread = 0.0;
    const auto& p = on.points.back();
    for (std::size_t k = 0; k < p.edge_current_n.size(); ++k) {
        spread = std::max(spread,
                          std::abs(p.edge_current_n[k] + p.edge_current_p[k] - J_on) / J_on);
    }
    CAPTURE(J_on, J_off, J_on / J_off - 1.0, spread);
    REQUIRE(std::abs(J_off - 1.280e-2) / 1.280e-2 < 1e-2);  // the Unit 9 gate, unchanged
    REQUIRE(J_on < J_off);
    REQUIRE(J_on / J_off > 0.9);  // a short-base diode at 0.5 V is diffusion-limited
    REQUIRE(spread < 1e-6);       // total current the same on every edge (legacy gate)
}

TEST_CASE("field mobility: the switch is in the run record, except in the quasi-static sweep") {
    const auto d = resistor(*mesh::make_tensor_grid(uniform_axis(L, 11)), true);
    const std::vector<std::vector<double>> points{{0.0, 0.0}};
    const auto on = solve::make_run_record(d, field_mobility(), points);
    const auto off = solve::make_run_record(d, {}, points);
    REQUIRE(on.input_identity != off.input_identity);
    bool found = false;
    for (const auto& [name, value] : on.settings) {
        if (name == "models.field_mobility") found = value == 1.0;
    }
    REQUIRE(found);
    solve::BiasOptions qs_on = field_mobility(), qs_off;
    qs_on.equations = qs_off.equations = solve::Equations::equilibrium_poisson;
    REQUIRE(solve::make_run_record(d, qs_on, points).input_identity ==
            solve::make_run_record(d, qs_off, points).input_identity);
}
