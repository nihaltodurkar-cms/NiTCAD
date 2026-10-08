// Bias solve of a silicon pn diode: the first end-to-end gates (ARCHITECTURE.md section 11,
// Unit 9):
// J(0.5 V) = 1.280e-2 A/cm^2 within 1%; the ideal-diode law; current continuity; mesh independence;
// uniform 2D/3D reproduce 1D; reverse bias (and the 6.10 pivot-ratio gate).
//
// Fixtures follow the legacy tests:
// - tests/test_device1d_native_gates.py::_diode: graded_mesh(2e-4, [1e-4], h_min=1e-8, h_max=1e-6)
//   (ratio 1.15, 250 nodes), 1e17 / 1e17, Caughey-Thomas mobility, Scharfetter SRH. The legacy also
//   enables Auger and band-gap narrowing; BGN is zero at 1e17 (below bgn_N0 = 1.3e17) and Auger is
//   not implemented (its effect is measured and recorded in ARCHITECTURE.md 6.2, Unit 9).
// - tests/test_validation.py::_diode: the same with ratio 1.12 (262 nodes) and BGN off.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <numeric>
#include <utility>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/assemble/bernoulli.hpp"
#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/newton.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;
using base::ErrorCode;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

std::vector<double> gate_mesh() { return legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6); }
std::vector<double> validation_mesh() { return legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6, 1.12); }

std::vector<double> uniform_axis(double length, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = length * i / (nodes - 1);
    return a;
}

// Anode (contact 0) on x_min, the p side; cathode (contact 1) on x_max.
device::Device diode(mesh::Mesh m, double NA = 1e17, double ND = 1e17) {
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

device::Device diode_1d(const std::vector<double>& x) { return diode(*mesh::make_tensor_grid(x)); }

// Anode current density of a sweep over the given anode voltages (cathode grounded), each point
// started from the previous one.
std::vector<double> sweep(const device::Device& d, const std::vector<double>& volts,
                          const solve::BiasOptions& options = {}) {
    std::vector<double> J;
    std::optional<results::NodeFields> state;
    for (const double V : volts) {
        const std::vector<double> bias{V, 0.0};
        auto s = solve::solve_bias(d, bias, options, state ? &*state : nullptr);
        REQUIRE(s.has_value());
        J.push_back(s->terminal_current[0]);
        state = std::move(s->fields);
    }
    return J;
}

}  // namespace

TEST_CASE("bias: the ported legacy graded_mesh follows the legacy specification") {
    // The legacy pytcad/mesh.py graded_mesh puts ceil(int dx / s(x)) cells under the spacing
    // target s(x) = min(h_max, h_min + (ratio - 1) |x - x_focus|), then limits neighbouring cells
    // to `ratio` and rescales them to span [0, L]. Checked here against that specification: the
    // cell count from the integral in closed form (ln((h_min + g d) / h_min) / g up to the cap,
    // then d / h_max), the end points, monotonicity and the gradient limit. (The node positions
    // are not compared with the legacy's own output, which would need the legacy run; the Unit 9
    // J(0.5 V) gate below exercises them against the legacy current.)
    struct Case {
        double L, focus, h_min, h_max, ratio;
        std::size_t nodes;  // the legacy fixtures' counts, also from the closed form below
    };
    for (const Case& c : {Case{2e-4, 1e-4, 1e-8, 1e-6, 1.15, 250},
                          Case{2e-4, 1e-4, 1e-8, 1e-6, 1.12, 262},
                          Case{2e-4, 1e-4, 5e-9, 5e-7, 1.15, 450}}) {
        CAPTURE(c.ratio, c.h_min);
        const double g = c.ratio - 1.0;
        const auto side = [&](double D) {  // int_0^D dd / min(h_max, h_min + g d)
            const double knee = (c.h_max - c.h_min) / g;
            if (D <= knee) return std::log((c.h_min + g * D) / c.h_min) / g;
            return std::log(c.h_max / c.h_min) / g + (D - knee) / c.h_max;
        };
        const double cells = std::ceil(side(c.focus) + side(c.L - c.focus));
        REQUIRE(static_cast<std::size_t>(cells) + 1 == c.nodes);
        const auto x = legacy_graded_mesh(c.L, c.focus, c.h_min, c.h_max, c.ratio);
        REQUIRE(x.size() == c.nodes);
        REQUIRE(x.front() == 0.0);
        REQUIRE(x.back() == c.L);
        for (std::size_t k = 1; k < x.size(); ++k) REQUIRE(x[k] > x[k - 1]);
        for (std::size_t k = 2; k < x.size(); ++k) {
            const double a = x[k - 1] - x[k - 2], b = x[k] - x[k - 1];
            REQUIRE(std::max(a, b) / std::min(a, b) <= c.ratio * (1.0 + 1e-9));
        }
    }
}

TEST_CASE("bias: J(0.5 V) = 1.280e-2 A/cm^2 within 1% (section 10 gate)") {
    const auto d = diode_1d(gate_mesh());
    const std::vector<double> bias{0.5, 0.0};
    const auto s = solve::solve_bias(d, bias);
    REQUIRE(s.has_value());
    const double J = s->terminal_current[0];
    UNSCOPED_INFO("J(0.5 V) = " << J << " A/cm^2, " << s->convergence.iterations.size()
                                << " iterations");
    REQUIRE(std::abs(J - 1.280e-2) / 1.280e-2 < 1e-2);
    // Kirchhoff, and continuity: the total current is the same on every edge. Both hold to the
    // Newton tolerance (relative updates below 1e-8); the gate is the legacy spread, 1e-6.
    REQUIRE(std::abs(s->terminal_current[0] + s->terminal_current[1]) <= 1e-6 * J);
    double spread = 0.0;
    for (std::size_t k = 0; k < s->edge_current_n.size(); ++k) {
        spread = std::max(spread, std::abs(s->edge_current_n[k] + s->edge_current_p[k] - J) / J);
    }
    UNSCOPED_INFO("Kirchhoff " << (s->terminal_current[0] + s->terminal_current[1]) / J
                               << ", spread " << spread);
    REQUIRE(spread <= 1e-6);
}

TEST_CASE("bias: ideal-diode law (legacy test_ideal_diode_law)") {
    const auto d = diode_1d(validation_mesh());
    std::vector<double> V;
    for (int k = 0; k <= 12; ++k) V.push_back(0.05 * k);
    const auto J = sweep(d, V);

    // Short-base saturation current: J_s = q n_i^2 (D_p / (W_n N_D) + D_n / (W_p N_A)), with the
    // neutral widths 1 um - W/2 and W the zero-bias depletion width.
    const double VT = base::thermal_voltage(300.0);
    const double ni = physics::intrinsic_density(physics::silicon(), 300.0);
    const double eps = 11.7 * base::eps0_F_per_cm;
    const double Vbi = VT * std::log(1e34 / (ni * ni));
    const double W = std::sqrt(2.0 * eps * Vbi / base::q_C * (1.0 / 1e17 + 1.0 / 1e17));
    const double Wn = 1e-4 - W / 2.0, Wp = Wn;
    const physics::Semiconductor si = physics::silicon();
    const double Dn =
        physics::caughey_thomas_mobility(si, physics::Carrier::electron, 1e17, 300.0) * VT;
    const double Dp =
        physics::caughey_thomas_mobility(si, physics::Carrier::hole, 1e17, 300.0) * VT;
    const double Js = base::q_C * ni * ni * (Dp / (Wn * 1e17) + Dn / (Wp * 1e17));
    const double ratio = J[10] / (Js * std::exp(V[10] / VT));  // V = 0.5
    UNSCOPED_INFO("J/J_ideal(0.5 V) = " << ratio);
    REQUIRE((ratio > 0.85 && ratio < 1.15));
    // Ideality within 0.02 of 1 for V >= 0.3 V.
    for (std::size_t k = 6; k + 1 < V.size(); ++k) {
        const double n = (V[k + 1] - V[k]) / (VT * std::log(J[k + 1] / J[k]));
        CAPTURE(V[k], n);
        REQUIRE(std::abs(n - 1.0) < 0.02);
    }
}

TEST_CASE("bias: ideality from a fit over 0.3 to 0.7 V (legacy g3)") {
    const auto d = diode_1d(gate_mesh());
    std::vector<double> V;
    for (int k = 0; k <= 14; ++k) V.push_back(0.05 * k);
    const auto J = sweep(d, V);
    // Least-squares slope of ln J against V over 0.3 ... 0.7 V.
    double sv = 0, sy = 0, svv = 0, svy = 0, count = 0;
    for (std::size_t k = 6; k < V.size(); ++k) {
        const double y = std::log(J[k]);
        sv += V[k];
        sy += y;
        svv += V[k] * V[k];
        svy += V[k] * y;
        count += 1;
    }
    const double slope = (count * svy - sv * sy) / (count * svv - sv * sv);
    const double ideality = 1.0 / (slope * base::thermal_voltage(300.0));
    UNSCOPED_INFO("ideality " << ideality);
    REQUIRE(std::abs(ideality - 1.0) < 0.05);
}

TEST_CASE("bias: the analysis extractors read the diode's ideality (Unit 24)") {
    // analysis::diode_fit on the sweep's curve is the least-squares fit above (legacy g3, V16:
    // ideality 1.004); the local ideality stays within 0.02 of 1 from 0.3 to 0.6 V (legacy
    // test_ideal_diode_law); the currents are far above their resolution.
    const auto d = diode_1d(gate_mesh());
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 14; ++k) points.push_back({0.05 * k, 0.0});
    const auto s = solve::sweep_bias(d, points);
    REQUIRE(s.has_value());
    REQUIRE(!s->stopped);
    const auto curve = analysis::current_curve(*s, 0, 0);
    REQUIRE(curve.has_value());
    REQUIRE(curve->resolution.size() == points.size());
    double sv = 0, sy = 0, svv = 0, svy = 0, count = 0;
    for (std::size_t k = 6; k < points.size(); ++k) {
        const double y = std::log(curve->y[k]);
        sv += curve->x[k];
        sy += y;
        svv += curve->x[k] * curve->x[k];
        svy += curve->x[k] * y;
        count += 1;
    }
    const double by_hand =
        (count * svv - sv * sv) / ((count * svy - sv * sy) * base::thermal_voltage(300.0));
    const auto fit = analysis::diode_fit(*curve, 0.3, 0.7, 300.0);
    REQUIRE(fit.has_value());
    const auto local = analysis::local_ideality(*analysis::slice(*curve, 0.3, 0.6), 300.0);
    const auto rs = analysis::series_resistance(*curve, 0.3, 0.7, 300.0);
    REQUIRE(rs.has_value());
    std::printf("diode: ideality %.6f (by hand %.6f) I_s %.4e, R_s fit: n %.6f R_s %.4e\n",
                fit->ideality.value, by_hand, fit->saturation_current.value, rs->ideality.value,
                rs->resistance.value);
    REQUIRE(fit->ideality.window.first == 6);
    REQUIRE(fit->ideality.window.last == 14);
    REQUIRE(!fit->ideality.below_resolution);
    REQUIRE(std::abs(fit->ideality.value / by_hand - 1.0) < 1e-12);
    REQUIRE(std::abs(fit->ideality.value - 1.0) < 0.05);
    for (const double n : *local) REQUIRE(std::abs(n - 1.0) < 0.02);
}

TEST_CASE("bias: a 2x finer mesh moves J(0.5 V) by less than 3% (legacy g4)") {
    const std::vector<double> bias{0.5, 0.0};
    const double coarse = solve::solve_bias(diode_1d(gate_mesh()), bias)->terminal_current[0];
    const auto fine_mesh = legacy_graded_mesh(2e-4, 1e-4, 5e-9, 5e-7);
    const double fine = solve::solve_bias(diode_1d(fine_mesh), bias)->terminal_current[0];
    UNSCOPED_INFO("coarse " << coarse << ", fine " << fine);
    REQUIRE(std::abs(fine - coarse) / std::abs(coarse) < 0.03);
}

TEST_CASE("bias: y-uniform 2D and 3D devices reproduce the 1D current") {
    const auto x = gate_mesh();
    const auto y = uniform_axis(1e-4, 3);
    const auto z = uniform_axis(5e-5, 3);
    const std::vector<double> bias{0.5, 0.0};
    const double J1 = solve::solve_bias(diode_1d(x), bias)->terminal_current[0];
    const auto s2 = solve::solve_bias(diode(*mesh::make_tensor_grid(x, y)), bias);
    const auto s3 = solve::solve_bias(diode(*mesh::make_tensor_grid(x, y, z)), bias);
    REQUIRE(s2.has_value());
    REQUIRE(s3.has_value());
    // 2D current is per unit depth (A/cm), over a 1 um wide device; 3D is in A over 1 um x 0.5 um.
    const double J2 = s2->terminal_current[0] / 1e-4;
    const double J3 = s3->terminal_current[0] / (1e-4 * 5e-5);
    UNSCOPED_INFO("J1 " << J1 << ", 2D " << J2 / J1 - 1 << ", 3D " << J3 / J1 - 1);
    REQUIRE(close(J2, J1, 1e-9));
    REQUIRE(close(J3, J1, 1e-9));
}

TEST_CASE("bias: reverse leakage is generation-limited, and its Jacobians are not flagged") {
    // Legacy test_reverse_saturation: J < 0 and the exponent of -J against (0.83 - V) lies in
    // (0.5, 1.5). Here the linear solver keeps its default pivot-ratio check (6.10, Unit 9 gate:
    // one-sided Scharfetter-Gummel links of order e^-40 must not be taken for a floating region).
    const auto d = diode_1d(validation_mesh());
    const std::vector<double> V{-0.5, -2.0, -8.0};
    std::vector<double> J;
    double smallest = 1.0;
    std::optional<results::NodeFields> state;
    for (const double v : V) {
        const std::vector<double> bias{v, 0.0};
        auto s = solve::solve_bias(d, bias, {}, state ? &*state : nullptr);
        REQUIRE(s.has_value());
        J.push_back(s->terminal_current[0]);
        smallest = std::min(smallest, s->convergence.smallest_pivot_ratio.value_or(1.0));
        state = std::move(s->fields);
    }
    UNSCOPED_INFO("J " << J[0] << ", " << J[1] << ", " << J[2] << "; smallest pivot ratio "
                       << smallest);
    for (const double j : J) REQUIRE(j < 0.0);
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t k = 0; k < V.size(); ++k) {
        const double lx = std::log(0.83 - V[k]), ly = std::log(-J[k]);
        sx += lx;
        sy += ly;
        sxx += lx * lx;
        sxy += lx * ly;
    }
    const double exponent = (3.0 * sxy - sx * sy) / (3.0 * sxx - sx * sx);
    UNSCOPED_INFO("exponent " << exponent);
    REQUIRE((exponent > 0.5 && exponent < 1.5));
    REQUIRE(smallest > linalg::SolverConfig{}.min_pivot_ratio);
}

TEST_CASE("bias: zero bias is equilibrium") {
    const auto d = diode_1d(gate_mesh());
    const std::vector<double> bias{0.0, 0.0};
    const auto s = solve::solve_bias(d, bias);
    REQUIRE(s.has_value());
    const auto eq = solve::solve_equilibrium(d);
    for (std::size_t i = 0; i < eq->fields.potential_V.size(); ++i) {
        REQUIRE(std::abs(s->fields.potential_V[i] - eq->fields.potential_V[i]) <= 1e-12);
    }
    for (const double I : s->terminal_current) REQUIRE(std::abs(I) <= 1e-10);
    REQUIRE(s->convergence.iterations.size() == 1);
}

TEST_CASE("bias: without recombination the electron current is the same on every edge") {
    const auto d = diode_1d(gate_mesh());
    const std::vector<double> bias{0.4, 0.0};
    const auto s = solve::solve_bias(d, bias, {.models = {.srh = false, .auger = false}});
    REQUIRE(s.has_value());
    // Exact up to rounding. In the neutral n region the electron flux is a difference of one-sided
    // terms n2 B(delta) and n1 B(-delta) up to 2.7e9 times larger (measured: deviation 4.9e-6 of
    // Jn there), so each edge is held to 64 eps of its own one-sided terms (q mu_max V_T bounds the
    // edge factor).
    const auto x = gate_mesh();
    const double VT = base::thermal_voltage(300.0);
    const double Jn = s->edge_current_n.front();  // p side: minority electrons, no cancellation
    for (std::size_t k = 0; k < s->edge_current_n.size(); ++k) {
        const double delta = (s->fields.potential_V[k + 1] - s->fields.potential_V[k]) / VT;
        const double one_sided = base::q_C * 1360.0 * VT / (x[k + 1] - x[k]) *
                                 (s->fields.n_cm3[k + 1] * assemble::bernoulli(delta) +
                                  s->fields.n_cm3[k] * assemble::bernoulli(-delta));
        CAPTURE(k);
        REQUIRE(std::abs(s->edge_current_n[k] - Jn) <= 64.0 * 2.2204460492503131e-16 * one_sided +
                                                           1e-9 * std::abs(Jn));
    }
    // With SRH the electron current changes along the device while the total does not.
    const auto with_srh = solve::solve_bias(d, bias);
    REQUIRE_FALSE(close(with_srh->edge_current_n.front(), with_srh->edge_current_n.back(), 1e-3));
}

TEST_CASE("bias: Newton analyzes once, and converges to a quadratic finish") {
    const auto d = diode_1d(gate_mesh());
    const auto scaling = *assemble::make_scaling(d);
    auto system = *assemble::DriftDiffusion::create(d, scaling);
    const auto eq = *solve::solve_equilibrium(d);
    std::vector<double> psi(eq.fields.potential_V.size());
    for (std::size_t i = 0; i < psi.size(); ++i) psi[i] = eq.fields.potential_V[i] / scaling.V_T;
    const std::vector<double> bias{0.5, 0.0};
    REQUIRE(system.set_bias(bias).has_value());
    auto x = system.state_from_potential(psi);
    auto solver = *linalg::LinearSolver::create({});
    results::ConvergenceRecord record;
    REQUIRE(solve::newton_solve(system, x, {}, solver, record).has_value());
    REQUIRE(solver.analyses() == 1);
    std::vector<double> u;
    for (const auto& r : record.iterations) u.push_back(r.update);
    REQUIRE(u.back() < 1e-8);
    REQUIRE(u[u.size() - 1] <= 10.0 * u[u.size() - 2] * u[u.size() - 2]);
}

TEST_CASE("bias: invalid input and non-convergence return errors") {
    const auto d = diode_1d(gate_mesh());
    const std::vector<double> one{0.5};
    REQUIRE(solve::solve_bias(d, one).error().code == ErrorCode::invalid_input);
    const std::vector<double> bad{std::nan(""), 0.0};
    REQUIRE(solve::solve_bias(d, bad).error().code == ErrorCode::invalid_input);
    const std::vector<double> bias{0.5, 0.0};
    results::NodeFields empty;
    REQUIRE(solve::solve_bias(d, bias, {}, &empty).error().code == ErrorCode::invalid_input);
    const auto few = solve::solve_bias(d, bias, {.newton = {.max_iterations = 2}});
    REQUIRE_FALSE(few.has_value());
    REQUIRE(few.error().code == ErrorCode::non_convergence);
}

