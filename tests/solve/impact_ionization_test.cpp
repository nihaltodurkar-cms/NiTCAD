// Impact ionization and breakdown (ARCHITECTURE.md section 11, Unit 19): the local model's current
// balance (the simulated reverse current of a one-sided junction against the ionization integral
// on its own field profile); avalanche breakdown at the ionization-integral condition (the
// primary check) and against the analytic one-sided breakdown voltage; the trace's points are
// steady states, its input errors, its stops and its progress; the run record; an open-base
// transistor's snapback; a curved junction against a planar one; y-uniform extrusions; transient
// and small-signal runs with the generation.
// The traces take about a minute each in Release: the breakdown gates, and the extrusion and
// transient runs with multiplication (minutes in Debug), carry the hidden tag
// [.breakdown] and run in the Release configuration only (ctest solve_breakdown).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <stop_token>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/control.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "NiTCAD/solve/waveform.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

// The legacy M15 one-sided abrupt junction: 6 um, the light p side (N) on x < 3 um, n+ 1e19 beyond;
// anode (contact 0) on x_min; legacy graded mesh, 10 nm at the junction. With tau, both lifetimes
// tau at every doping.
device::Device one_sided(double N, double tau = 0.0) {
    physics::SemiconductorParameters si = physics::silicon_parameters;
    if (tau > 0.0) si.lifetime = {.tau_n0 = tau, .tau_p0 = tau, .N_ref = 1e30};
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(6e-4, {3e-4}, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 3e-4) {
            acceptors[i] = N;
        } else {
            donors[i] = 1e19;
        }
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(si)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

const physics::ImpactIonizationCoefficients& electrons() {
    static const physics::Semiconductor si = physics::silicon();
    return physics::impact_ionization(si, physics::Carrier::electron);
}
const physics::ImpactIonizationCoefficients& holes() {
    static const physics::Semiconductor si = physics::silicon();
    return physics::impact_ionization(si, physics::Carrier::hole);
}

// On a 1D state, from x_min (where the electrons enter the depleted p side): phi at the nodes,
// phi(x) = integral_0^x (alpha_n - alpha_p), and the electron ionization integral
// I_n = integral alpha_n e^-phi, alpha at each edge's field, by the midpoint rule per edge.
// I_n = 1 is avalanche breakdown.
struct Ionization {
    std::vector<double> phi;
    double integral;
};

Ionization ionization(const device::Device& d, const std::vector<double>& potential) {
    Ionization r{std::vector<double>(potential.size(), 0.0), 0.0};
    for (std::size_t i = 0; i + 1 < potential.size(); ++i) {
        const double h = d.mesh().points()[i + 1][0] - d.mesh().points()[i][0];
        const double E = std::abs(potential[i + 1] - potential[i]) / h;
        const double an = physics::impact_ionization_coefficient(electrons(), 1.0, E).alpha;
        const double ap = physics::impact_ionization_coefficient(holes(), 1.0, E).alpha;
        r.integral += an * std::exp(-(r.phi[i] + 0.5 * (an - ap) * h)) * h;
        r.phi[i + 1] = r.phi[i] + (an - ap) * h;
    }
    return r;
}

// The current a 1D state carries by the local model's balance: integrating the continuity of the
// electron current, dJ_n/dx = alpha_n J_n + alpha_p J_p + q g, with J = J_n + J_p constant,
//     J (1 - I_n) = J_n0 + J_p0 e^-phi(W) + q integral g e^-phi,
// g the net SRH generation (both lifetimes tau, midgap level, n_ie the intrinsic density: band-gap
// narrowing and Auger off), J_n0 and J_p0 the electron current entering at x_min and the hole
// current entering at x_max, every current counted towards -x (a reverse current).
double balance_current(const device::Device& d, const results::BiasPoint& p, double tau) {
    const Ionization ion = ionization(d, p.fields.potential_V);
    const double n_i = physics::intrinsic_density(physics::silicon(), 300.0);
    double source = 0.0;
    for (std::size_t i = 1; i + 1 < ion.phi.size(); ++i) {
        const double n = p.fields.n_cm3[i], h = p.fields.p_cm3[i];
        const double U = (n * h - n_i * n_i) / (tau * (n + n_i) + tau * (h + n_i));
        source += d.mesh().volumes()[i] * (-U) * std::exp(-ion.phi[i]);
    }
    const double Jn0 = -p.edge_current_n.front(), Jp0 = -p.edge_current_p.back();
    return (Jn0 + Jp0 * std::exp(-ion.phi.back()) + base::q_C * source) / (1.0 - ion.integral);
}

// The analytic one-sided breakdown voltage (the legacy analysis layer): the depletion
// approximation's triangular field E(x) = q N x / eps on [0, W], W = sqrt(2 eps (V_bi + V) / q N),
// V_bi = V_T ln(N N_D / n_i^2), and the electron ionization integral from the depletion edge
// (where the electrons enter) reaching 1; by bisection, the integral by the midpoint rule on 20000
// cells.
double analytic_breakdown(double N, double N_D = 1e19) {
    const double eps = physics::silicon_parameters.eps_r * base::eps0_F_per_cm;
    const double V_T = base::thermal_voltage(300.0);
    const double n_i = physics::intrinsic_density(physics::silicon(), 300.0);
    const double Vbi = V_T * std::log(N * N_D / (n_i * n_i));
    const auto integral = [&](double V) {
        const double W = std::sqrt(2.0 * eps * (Vbi + V) / (base::q_C * N));
        const int cells = 20000;
        const double h = W / cells;
        double phi = 0.0, I = 0.0;
        for (int k = 0; k < cells; ++k) {
            const double E = base::q_C * N * (k + 0.5) * h / eps;
            const double an = physics::impact_ionization_coefficient(electrons(), 1.0, E).alpha;
            const double ap = physics::impact_ionization_coefficient(holes(), 1.0, E).alpha;
            I += an * std::exp(-(phi + 0.5 * (an - ap) * h)) * h;
            phi += (an - ap) * h;
        }
        return I;
    };
    double lo = 1.0, hi = 500.0;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        (integral(mid) < 1.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

solve::TraceOptions reverse_trace(double end_V, double current_limit, double tau) {
    solve::TraceOptions o;
    o.steady.models.bgn = false;
    if (tau > 0.0) o.steady.models.auger = false;
    o.steady.models.impact_ionization = true;
    o.start_V = {0.0, 0.0};
    o.end_V = end_V;
    o.step_V = 0.5;
    o.max_step_V = 2.0;
    o.current_limit = current_limit;
    return o;
}

}  // namespace

TEST_CASE("impact ionization: the reverse current follows the local model's balance",
          "[.breakdown]") {
    // Lifetimes of 10 ps put the generation current (1e-3 to 1e-1 A/cm^2) well above the state's
    // current resolution (1e-4), so the multiplied current is resolved at every bias. The
    // simulated current against J (1 - I_n) = J_n0 + J_p0 e^-phi(W) + q integral g e^-phi on its
    // own field and densities: measured below 1e-6 at I_n below 0.25, 2.6e-4 at I_n = 0.94
    // (multiplication 18) and 0.97% at I_n = 0.998 (multiplication 500), where 1 / (1 - I_n)
    // amplifies the discretization error of the integral.
    const double tau = 1e-11;
    const device::Device d = one_sided(1e16, tau);
    const auto t = solve::trace_bias(d, reverse_trace(-80.0, 1.0, tau));
    REQUIRE(t.has_value());
    REQUIRE(!t->stopped);
    double worst_low = 0.0, worst = 0.0, deepest = 0.0;
    for (std::size_t k = 1; k < t->points.size(); ++k) {
        const results::BiasPoint& p = t->points[k];
        const double J = -p.terminal_current[0];
        const double I = ionization(d, p.fields.potential_V).integral;
        const double miss = std::abs(J / balance_current(d, p, tau) - 1.0);
        REQUIRE(J > 3.0 * p.terminal_current_resolution[0]);
        deepest = std::max(deepest, I);
        if (I <= 0.95) {
            worst_low = std::max(worst_low, miss);
        } else {
            worst = std::max(worst, miss);
        }
    }
    CAPTURE(t->points.size(), worst_low, worst, deepest);
    REQUIRE(deepest > 0.99);
    REQUIRE(worst_low < 5e-4);
    REQUIRE(worst < 2e-2);
}

TEST_CASE("impact ionization: avalanche breakdown at the ionization-integral condition",
          "[.breakdown]") {
    // The primary check: breakdown is where the electron ionization integral on the simulated
    // field profile reaches 1, and there the simulated current runs away. One-sided junctions
    // (N = 1e16 and 2e16 under n+ 1e19) with 10 ps lifetimes, so the current is resolved all the
    // way through breakdown (the lifetime sets the leakage, not the field), traced past the point
    // where the integral crosses 1: the bias of the crossing (interpolated between the traced
    // points) against the analytic one-sided breakdown voltage of the depletion approximation
    // (measured 55.34 against 55.26 V and 35.95 against 35.86 V); there the current is more than
    // 10 times its value at 90% of that bias. The legacy fixture itself (N = 1e16, the legacy
    // lifetimes, its leakage four decades below the state's current resolution) traced until the
    // current reaches 1e-4 A/cm^2, above the resolution: the integral there is 1 (measured
    // 0.99965) and the bias within 1% of the analytic value. Its trace crawls through the last
    // 0.1 V before breakdown, where the current is not resolved (ARCHITECTURE.md 6.2, Unit 19), so
    // its corrector is capped at 25 iterations.
    for (const auto& [N, end] : {std::pair{1e16, -80.0}, std::pair{2e16, -55.0}}) {
        const double tau = 1e-11;
        const device::Device d = one_sided(N, tau);
        const auto t = solve::trace_bias(d, reverse_trace(end, 10.0, tau));
        REQUIRE(t.has_value());
        REQUIRE(!t->stopped);
        double V_bd = 0.0, J_bd = 0.0;
        double previous_I = 0.0, previous_V = 0.0;
        for (const auto& p : t->points) {
            const double I = ionization(d, p.fields.potential_V).integral;
            if (I >= 1.0 && previous_I < 1.0) {
                const double f = (1.0 - previous_I) / (I - previous_I);
                V_bd = previous_V + f * (-p.bias_V[0] - previous_V);
                J_bd = -p.terminal_current[0];
            }
            previous_I = I;
            previous_V = -p.bias_V[0];
        }
        double J_before = 0.0;  // at 90% of the breakdown voltage
        for (const auto& p : t->points) {
            if (-p.bias_V[0] <= 0.9 * V_bd) J_before = -p.terminal_current[0];
        }
        const double analytic = analytic_breakdown(N);
        CAPTURE(N, V_bd, analytic, V_bd / analytic - 1.0, J_bd, J_before, t->points.size());
        REQUIRE(V_bd > 0.0);
        REQUIRE(std::abs(V_bd / analytic - 1.0) < 1e-2);
        REQUIRE(J_bd > 10.0 * J_before);
    }
    {
        const device::Device d = one_sided(1e16);
        auto o = reverse_trace(-80.0, 1e-4, 0.0);
        o.steady.newton.max_iterations = 25;
        const auto t = solve::trace_bias(d, o);
        REQUIRE(t.has_value());
        REQUIRE(!t->stopped);
        const results::BiasPoint& last = t->points.back();
        const double V_bd = -last.bias_V[0];
        const double I = ionization(d, last.fields.potential_V).integral;
        const double analytic = analytic_breakdown(1e16);
        double below = 0.0;  // the current at 90% of the breakdown voltage
        for (const auto& p : t->points) {
            if (-p.bias_V[0] <= 0.9 * V_bd) below = -p.terminal_current[0];
        }
        CAPTURE(V_bd, I, analytic, below, t->points.size());
        REQUIRE(-last.terminal_current[0] >= 1e-4);
        REQUIRE(-last.terminal_current[0] > 3.0 * last.terminal_current_resolution[0]);
        REQUIRE(std::abs(I - 1.0) < 2e-3);
        REQUIRE(std::abs(V_bd / analytic - 1.0) < 1e-2);
        REQUIRE(below < 1e-2 * 1e-4);
    }
}

namespace {

// A symmetric 1e17 diode, 2 um, anode (contact 0) on x_min.
device::Device forward_diode(
    const physics::SemiconductorParameters& si = physics::silicon_parameters) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, {1e-4}, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(si)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

solve::TraceOptions forward_trace() {
    solve::TraceOptions o;
    o.start_V = {0.0, 0.0};
    o.end_V = 0.6;
    o.step_V = 0.05;
    o.max_step_V = 0.2;
    return o;
}

}  // namespace

TEST_CASE("impact ionization: the trace's points are steady states") {
    // Forward bias without ionization (a monotonic branch): each traced point, re-solved by
    // sweep_bias at its own bias from its own state, is that state (currents to 1e-9); the biases
    // rise monotonically and the last lands exactly on end_V.
    const device::Device d = forward_diode();
    const auto t = solve::trace_bias(d, forward_trace());
    REQUIRE(t.has_value());
    REQUIRE(!t->stopped);
    REQUIRE(t->points.back().bias_V[0] == 0.6);
    double worst = 0.0;
    for (std::size_t k = 1; k < t->points.size(); ++k) {
        const results::BiasPoint& p = t->points[k];
        REQUIRE(p.bias_V[0] > t->points[k - 1].bias_V[0]);
        const auto again = solve::solve_bias(d, p.bias_V, {}, &p.fields);
        REQUIRE(again.has_value());
        worst = std::max(worst, std::abs(again->terminal_current[0] - p.terminal_current[0]) /
                                    p.terminal_current_resolution[0]);
    }
    CAPTURE(worst, t->points.size());
    REQUIRE(worst <= 1.0);
}

TEST_CASE("impact ionization: trace input errors") {
    const device::Device d = forward_diode();
    const auto error = [&](const solve::TraceOptions& o) {
        const auto r = solve::trace_bias(d, o);
        REQUIRE_FALSE(r.has_value());
        return r.error().code;
    };
    const auto with = [](auto change) {
        solve::TraceOptions o = forward_trace();
        change(o);
        return o;
    };
    using base::ErrorCode;
    REQUIRE(error(with([](auto& o) {
                o.steady.equations = solve::Equations::equilibrium_poisson;
            })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.contact = 2; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.start_V = {0.0}; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.start_V = {std::nan(""), 0.0}; })) ==
            ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.end_V = 0.0; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.end_V = std::nan(""); })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.step_V = 0.0; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.step_V = 0.5; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.min_step_V = 0.1; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.current_limit = -1.0; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.max_points = 0; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) { o.density_floor_cm3 = -1.0; })) == ErrorCode::invalid_input);
    REQUIRE(error(with([](auto& o) {
                o.steady.models.impact_ionization = true;
                o.steady.models.impact_current_resolution = -1.0;
            })) == ErrorCode::invalid_input);
}

TEST_CASE("impact ionization: a trace stops at its current limit, point count, or on request") {
    const device::Device d = forward_diode();
    {
        solve::TraceOptions o = forward_trace();
        o.current_limit = 0.1;  // A/cm^2, reached before 0.6 V (0.61 A/cm^2)
        const auto t = solve::trace_bias(d, o);
        REQUIRE(t.has_value());
        REQUIRE(!t->stopped);
        REQUIRE(t->points.back().terminal_current[0] >= 0.1);
        REQUIRE(t->points[t->points.size() - 2].terminal_current[0] < 0.1);
        REQUIRE(t->points.back().bias_V[0] < 0.6);
    }
    {
        solve::TraceOptions o = forward_trace();
        o.max_points = 3;
        const auto t = solve::trace_bias(d, o);
        REQUIRE(t.has_value());
        REQUIRE(!t->stopped);
        REQUIRE(t->points.size() == 3);
    }
    // Cancelled during the third step: the starting point and two steps are kept; progress is the
    // starting point (point 0 of 1), then each step attempt (point = attempt + 1, count 0).
    std::stop_source stop;
    std::vector<solve::Progress> events;
    solve::RunControl control{stop.get_token(), [&](const solve::Progress& p) {
                                  events.push_back(p);
                                  if (p.point == 3 && p.iteration == 1) stop.request_stop();
                              }};
    const auto t = solve::trace_bias(d, forward_trace(), nullptr, control);
    REQUIRE(t.has_value());
    REQUIRE(t->stopped.has_value());
    REQUIRE(t->stopped->code == base::ErrorCode::cancelled);
    REQUIRE(t->points.size() == 3);
    const auto key = [](const solve::Progress& p) {
        return std::tuple{static_cast<int>(p.phase), p.point, p.iteration};
    };
    for (std::size_t e = 1; e < events.size(); ++e) REQUIRE(key(events[e - 1]) < key(events[e]));
    std::size_t bias_events = 0;
    for (const solve::Progress& p : events) {
        if (p.phase != solve::Phase::bias) continue;
        ++bias_events;
        REQUIRE(p.point_count == (p.point == 0 ? 1u : 0u));
    }
    REQUIRE(bias_events > 0);
    REQUIRE(events.back().point == 3);
}

TEST_CASE("impact ionization: the run record") {
    // The model's settings and coefficients enter a run's identity only when it is on, so the runs
    // without it keep theirs; the trace's own options enter it too.
    const device::Device d = forward_diode();
    const std::vector<std::vector<double>> points{{0.1, 0.0}};
    solve::BiasOptions off;
    solve::BiasOptions on = off;
    on.models.impact_ionization = true;
    solve::BiasOptions resolved = on;
    resolved.models.impact_current_resolution = 0.1;
    const auto record = [&](const solve::BiasOptions& o) {
        return solve::make_run_record(d, o, points);
    };
    const auto has = [](const results::RunRecord& r, const std::string& name) {
        return std::ranges::any_of(r.settings, [&](const auto& s) { return s.first == name; });
    };
    REQUIRE(!has(record(off), "models.impact_ionization"));
    REQUIRE(!has(record(off), "models.impact_current_resolution"));
    REQUIRE(has(record(on), "models.impact_ionization"));
    REQUIRE(has(record(on), "models.impact_current_resolution"));
    REQUIRE(record(off).input_identity != record(on).input_identity);
    REQUIRE(record(on).input_identity != record(resolved).input_identity);
    // The coefficients: a device with another hole set differs only with the model on.
    physics::SemiconductorParameters si = physics::silicon_parameters;
    si.impact_ionization.hole.A_high_per_cm *= 2.0;
    si.impact_ionization.phonon_energy_eV = 0.05;
    const device::Device d2 = forward_diode(si);
    REQUIRE(solve::make_run_record(d2, off, points).input_identity == record(off).input_identity);
    REQUIRE(solve::make_run_record(d2, on, points).input_identity != record(on).input_identity);
    // The trace.
    solve::TraceOptions a = forward_trace();
    solve::TraceOptions b = a;
    b.max_step_V = 0.3;
    solve::TraceOptions c = a;
    c.current_limit = 1.0;
    const auto id = [&](const solve::TraceOptions& o) {
        return solve::make_run_record(d, o).input_identity;
    };
    REQUIRE(id(a) != id(b));
    REQUIRE(id(a) != id(c));
    REQUIRE(id(a) == solve::trace_bias(d, a)->run.input_identity);
    REQUIRE(has(solve::make_run_record(d, c), "trace.current_limit"));
}

namespace {

// An open-base n+ p n- n+ transistor (1D, two contacts) at T: emitter n+ 1e19 on [0, 1) um, base
// p 1e17 on [1, 1.5) um, collector n- 1e16 to 5 um, n+ 1e19 to 6 um; emitter on x_min, collector on
// x_max; lifetimes 1e-7 s (Scharfetter, N_ref 5e16).
device::Device open_base(double T) {
    physics::SemiconductorParameters si = physics::silicon_parameters;
    si.lifetime = {.tau_n0 = 1e-7, .tau_p0 = 1e-7, .N_ref = 5e16};
    const double xb = 1e-4, xc = 1.5e-4;
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(6e-4, {xb, xc, 5e-4}, 2e-8, 2e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = m.points()[i][0];
        if (x < xb) {
            donors[i] = 1e19;
        } else if (x < xc) {
            acceptors[i] = 1e17;
        } else if (x < 5e-4) {
            donors[i] = 1e16;
        } else {
            donors[i] = 1e19;
        }
    }
    auto e = m.find_boundary("x_min")->nodes;
    auto c = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"silicon", *physics::Semiconductor::create(si)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"emitter", device::ContactKind::ohmic, std::move(e)},
                      {"collector", device::ContactKind::ohmic, std::move(c)}}});
}

}  // namespace

TEST_CASE("impact ionization: an open-base transistor snaps back", "[.breakdown]") {
    // At 350 K (the leakage, 1e-5 A/cm^2 and more, above the state's resolution) the collector
    // current of the open-base transistor rises with the bias until avalanche multiplication
    // turns the base current around (BV_CEO, measured 17.27 V at 1.8 mA/cm^2); then the current
    // keeps rising while the bias falls: the curve folds back, and a bias sweep, which can only
    // step the bias, stays on the low-current branch (measured 2.2e-4 A/cm^2 at 9.92 V, where the
    // trace carries 23 A/cm^2). On the high-current branch the base current vanishes as the holes
    // the collector generates replace it: alpha_T M_n = 1, with alpha_T the electron fraction of
    // the current entering the collector at the metallurgical junction and M_n = 1 / (1 - I_n),
    // I_n the electron ionization integral over the collector (measured 0.9993 at 35 mA/cm^2,
    // 1.0002 at 0.58 A/cm^2; below 1 at the fold by the collector leakage's share). Every traced
    // point, re-solved at its own bias from its own state, is that state within the current
    // resolution.
    const double T = 350.0, xc = 1.5e-4;
    const device::Device d = open_base(T);
    solve::TraceOptions o;
    o.steady.models.bgn = false;
    o.steady.models.impact_ionization = true;
    o.contact = 1;
    o.start_V = {0.0, 0.0};
    o.end_V = 80.0;
    o.step_V = 0.5;
    o.max_step_V = 2.0;
    o.current_limit = 20.0;
    const auto t = solve::trace_bias(d, o);
    REQUIRE(t.has_value());
    REQUIRE(!t->stopped);
    std::size_t top = 0;
    for (std::size_t k = 0; k < t->points.size(); ++k) {
        if (t->points[k].bias_V[1] > t->points[top].bias_V[1]) top = k;
    }
    const results::BiasPoint& fold = t->points[top];
    const results::BiasPoint& last = t->points.back();
    CAPTURE(fold.bias_V[1], fold.terminal_current[1], last.bias_V[1], last.terminal_current[1]);
    REQUIRE(top > 0);
    REQUIRE(top + 1 < t->points.size());
    REQUIRE(last.terminal_current[1] >= 20.0);
    REQUIRE(last.bias_V[1] < 0.7 * fold.bias_V[1]);
    REQUIRE(last.terminal_current[1] > 1e3 * fold.terminal_current[1]);
    // Past the fold the bias falls monotonically while the current rises.
    for (std::size_t k = top + 1; k < t->points.size(); ++k) {
        REQUIRE(t->points[k].terminal_current[1] > t->points[k - 1].terminal_current[1]);
        REQUIRE(t->points[k].bias_V[1] < t->points[k - 1].bias_V[1]);
    }
    // The sustaining condition, and the points are steady states.
    const double gamma = physics::impact_ionization_temperature_factor(0.063, T);
    double worst_resolve = 0.0, high = 0.0;
    for (std::size_t k = top; k < t->points.size(); ++k) {
        const results::BiasPoint& p = t->points[k];
        double phi = 0.0, In = 0.0;
        std::size_t at = 0;  // the first edge in the collector
        for (std::size_t i = 0; i + 1 < p.fields.potential_V.size(); ++i) {
            const double x0 = d.mesh().points()[i][0], x1 = d.mesh().points()[i + 1][0];
            if (x1 <= xc) {
                at = i + 1;
                continue;
            }
            const double h = x1 - x0;
            const double E = std::abs(p.fields.potential_V[i + 1] - p.fields.potential_V[i]) / h;
            const double an = physics::impact_ionization_coefficient(electrons(), gamma, E).alpha;
            const double ap = physics::impact_ionization_coefficient(holes(), gamma, E).alpha;
            In += an * std::exp(-(phi + 0.5 * (an - ap) * h)) * h;
            phi += (an - ap) * h;
        }
        const double J = p.terminal_current[1];
        const double alpha_T = -p.edge_current_n[at] / J;
        if (J >= 0.1) high = std::max(high, std::abs(alpha_T / (1.0 - In) - 1.0));
        const auto again = solve::solve_bias(d, p.bias_V, o.steady, &p.fields);
        REQUIRE(again.has_value());
        worst_resolve = std::max(worst_resolve, std::abs(again->terminal_current[1] - J) /
                                                    p.terminal_current_resolution[1]);
    }
    std::vector<std::vector<double>> ramp;
    for (double v = 0.5; v < last.bias_V[1]; v += 0.5) ramp.push_back({0.0, v});
    ramp.push_back(last.bias_V);
    const auto sweep = solve::sweep_bias(d, ramp, o.steady);
    REQUIRE(sweep.has_value());
    REQUIRE(!sweep->stopped);
    const double low = sweep->points.back().terminal_current[1];
    CAPTURE(high, worst_resolve, low);
    REQUIRE(high < 1e-3);
    REQUIRE(worst_resolve <= 1.0);
    REQUIRE(low < 1e-2 * last.terminal_current[1]);
}

namespace {

// A one-sided junction in 2D: p (N, lifetimes tau) with n+ 1e19 on x >= 3 um, either across the
// whole height (planar) or only for y < 1 um (its edge a cylindrical junction); anode on x_min,
// cathode on the n+ region's part of x_max. 6 um by 4 um.
device::Device junction_2d(bool curved, double N, double tau, double h_min = 2e-7) {
    physics::SemiconductorParameters si = physics::silicon_parameters;
    si.lifetime = {.tau_n0 = tau, .tau_p0 = tau, .N_ref = 1e30};
    const auto x = legacy_graded_mesh(6e-4, {3e-4}, h_min, 2e-5);
    const auto y = legacy_graded_mesh(4e-4, {1e-4}, h_min, 2e-5);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if (p[0] >= 3e-4 && (!curved || p[1] <= 1e-4)) {
            donors[i] = 1e19;
        } else {
            acceptors[i] = N;
        }
    }
    std::vector<mesh::NodeId> cathode;
    for (const mesh::NodeId v : m.find_boundary("x_max")->nodes) {
        if (!curved || m.points()[static_cast<std::size_t>(v)][1] <= 1e-4) cathode.push_back(v);
    }
    auto anode = m.find_boundary("x_min")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(si)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

}  // namespace

TEST_CASE("impact ionization: a curved junction breaks down before a planar one",
          "[.breakdown]") {
    // 2D, 6 um by 4 um of p 1e16 (lifetimes 10 ps, so the generation current is resolved), with an
    // n+ 1e19 region on x >= 3 um: across the whole height (planar, y-uniform) or only for
    // y <= 1 um (its edge a cylindrical junction of radius about 1 um); anode on x_min, cathode on
    // the n+ region's part of x_max. The field crowds at the curved edge. Swept to -27 V (1 V
    // steps, 0.25 V from -25 V), the current's growth over its value at -10 V: measured 5.7 times
    // (curved) against 2.2 (planar; the local model's multiplication at -27 V is 1.3). Traced, the
    // curved junction's current reaches 100 times its -10 V value at 27.8 V and the planar
    // junction's at 55.4 V (the 1D breakdown voltage is 55.3); those traces take minutes in 2D
    // (ARCHITECTURE.md 6.2, Unit 19), so the gate uses the sweep.
    const auto grown = [](bool curved) {
        const device::Device d = junction_2d(curved, 1e16, 1e-11);
        solve::BiasOptions o;
        o.models.bgn = false;
        o.models.auger = false;
        o.models.impact_ionization = true;
        std::vector<std::vector<double>> points;
        for (int v = 1; v <= 25; ++v) points.push_back({-static_cast<double>(v), 0.0});
        for (int k = 1; k <= 8; ++k) points.push_back({-25.0 - 0.25 * k, 0.0});
        const auto s = solve::sweep_bias(d, points, o);
        REQUIRE(s.has_value());
        REQUIRE(!s->stopped);
        return s->points.back().terminal_current[0] / s->points[9].terminal_current[0];
    };
    const double planar = grown(false), curved = grown(true);
    CAPTURE(planar, curved, curved / planar);
    REQUIRE(planar < 2.5);
    REQUIRE(curved > 2.0 * planar);
}

namespace {

// The 10 ps one-sided junction (N = 1e16) on a tensor grid with the legacy x axis and, for D > 1,
// y (and z) uniform over 0.2 um in 2 nodes.
device::Device one_sided_extruded(int D) {
    physics::SemiconductorParameters si = physics::silicon_parameters;
    si.lifetime = {.tau_n0 = 1e-11, .tau_p0 = 1e-11, .N_ref = 1e30};
    const auto x = legacy_graded_mesh(6e-4, {3e-4}, 1e-8, 1e-6);
    const std::vector<double> t{0.0, 2e-5};
    mesh::Mesh m = D == 1   ? *mesh::make_tensor_grid(x)
                   : D == 2 ? *mesh::make_tensor_grid(x, t)
                            : *mesh::make_tensor_grid(x, t, t);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 3e-4) {
            acceptors[i] = 1e16;
        } else {
            donors[i] = 1e19;
        }
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(si)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

solve::BiasOptions multiplied() {
    solve::BiasOptions o;
    o.models.bgn = false;
    o.models.auger = false;
    o.models.impact_ionization = true;
    return o;
}

std::vector<std::vector<double>> ramp_to(double V, double step) {
    std::vector<std::vector<double>> points;
    for (double v = -step; v > V; v -= step) points.push_back({v, 0.0});
    points.push_back({V, 0.0});
    return points;
}

}  // namespace

TEST_CASE("impact ionization: y-uniform 2D and 3D reproduce 1D", "[.breakdown]") {
    // At -45 V (multiplication about 4): the current per unit width (0.2 um) and area reproduces
    // the 1D current density; the transverse edges carry no current and no field, so the
    // reconstruction adds nothing transverse (measured 0 in 2D and 2.2e-16 in 3D).
    const auto points = ramp_to(-45.0, 9.0);
    const auto one = solve::sweep_bias(one_sided_extruded(1), points, multiplied());
    const auto two = solve::sweep_bias(one_sided_extruded(2), points, multiplied());
    const auto three = solve::sweep_bias(one_sided_extruded(3), points, multiplied());
    REQUIRE(one.has_value());
    REQUIRE(two.has_value());
    REQUIRE(three.has_value());
    REQUIRE(!one->stopped);
    REQUIRE(!two->stopped);
    REQUIRE(!three->stopped);
    const double J1 = one->points.back().terminal_current[0];
    const double J2 = two->points.back().terminal_current[0] / 2e-5;
    const double J3 = three->points.back().terminal_current[0] / 4e-10;
    CAPTURE(J1, J2 / J1 - 1.0, J3 / J1 - 1.0);
    REQUIRE(std::abs(J2 / J1 - 1.0) < 1e-9);
    REQUIRE(std::abs(J3 / J1 - 1.0) < 1e-9);
}

TEST_CASE("impact ionization: transient and small-signal runs carry the generation",
          "[.breakdown]") {
    // The 10 ps junction at -45 V (multiplication about 4): a transient run at that constant bias
    // stays at the steady current (Unit 21's steps include the generation), and the small-signal
    // conductance at f = 0 is the derivative of the multiplied DC current (Unit 22; fourth-order
    // differences, h = 0.1 V). Measured: drift 6.1e-11; Y against the differences 4.7e-6, which
    // is the differences' own scatter (about 3e-6 between h = 0.1 and 0.0125 V, each difference
    // carrying the bias solves' convergence error).
    const device::Device d = one_sided_extruded(1);
    const solve::BiasOptions o = multiplied();
    const auto sweep = solve::sweep_bias(d, ramp_to(-45.0, 5.0), o);
    REQUIRE(sweep.has_value());
    REQUIRE(!sweep->stopped);
    const results::BiasPoint& dc = sweep->points.back();
    // Transient.
    const std::vector<solve::Waveform> w{solve::Waveform::constant(-45.0),
                                         solve::Waveform::constant(0.0)};
    const auto t = solve::solve_transient(
        d, w, {.steady = o, .t_end_s = 1e-9, .output_times_s = {}}, &dc.fields);
    REQUIRE(t.has_value());
    REQUIRE(!t->stopped);
    double drift = 0.0;
    for (const auto& p : t->points) {
        drift = std::max(drift, std::abs(p.terminal_current[0] / dc.terminal_current[0] - 1.0));
    }
    // Small signal.
    const double h = 0.1;
    const auto I = [&](double V) {
        return solve::solve_bias(d, std::vector<double>{V, 0.0}, o, &dc.fields)
            ->terminal_current[0];
    };
    const double G = (8.0 * (I(-45.0 + h) - I(-45.0 - h)) - (I(-45.0 + 2 * h) - I(-45.0 - 2 * h))) /
                     (12.0 * h);
    const auto ac = solve::solve_small_signal(d, std::vector<std::vector<double>>{{-45.0, 0.0}},
                                              {.steady = o, .frequencies_Hz = {0.0}}, &dc.fields);
    REQUIRE(ac.has_value());
    REQUIRE(!ac->stopped);
    const double Y = ac->conductance(0, 0, 0, 0);
    CAPTURE(drift, G, Y, Y / G - 1.0);
    REQUIRE(drift < 1e-9);
    REQUIRE(std::abs(Y / G - 1.0) < 2e-5);
}
