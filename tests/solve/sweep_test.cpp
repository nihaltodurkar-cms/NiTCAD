// Sweeps, cancellation, progress and the run record (ARCHITECTURE.md section 11, Unit 10 gates:
// cancellation returns partial results; progress is monotonic).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <latch>
#include <optional>
#include <semaphore>
#include <stop_token>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/control.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;
using base::ErrorCode;
using solve::Phase;
using solve::Progress;

namespace {

device::Device diode(double donor_scale = 1.0) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 1e-4) {
            acceptors[i] = 1e17;
        } else {
            donors[i] = 1e17 * donor_scale;
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

// Anode at 0, 0.1, ..., 0.5 V; cathode grounded.
std::vector<std::vector<double>> forward_points() {
    std::vector<std::vector<double>> p;
    for (int k = 0; k <= 5; ++k) p.push_back({0.1 * k, 0.0});
    return p;
}

bool same_point(const results::BiasPoint& a, const results::BiasPoint& b) {
    return a.bias_V == b.bias_V && a.fields.potential_V == b.fields.potential_V &&
           a.fields.n_cm3 == b.fields.n_cm3 && a.terminal_current == b.terminal_current;
}

// Equal to a relative 1e-9: a chain of one-point solves passes each state through physical units
// (psi V_T / V_T), which moves the starting point by rounding.
bool close_point(const results::BiasPoint& a, const results::BiasPoint& b) {
    const auto close = [](const std::vector<double>& u, const std::vector<double>& v) {
        if (u.size() != v.size()) return false;
        for (std::size_t i = 0; i < u.size(); ++i) {
            if (std::abs(u[i] - v[i]) > 1e-9 * std::max(std::abs(u[i]), std::abs(v[i]))) {
                return false;
            }
        }
        return true;
    };
    return a.bias_V == b.bias_V && close(a.fields.potential_V, b.fields.potential_V) &&
           close(a.fields.n_cm3, b.fields.n_cm3) && close(a.terminal_current, b.terminal_current);
}

auto key(const Progress& e) { return std::tuple(static_cast<int>(e.phase), e.point, e.iteration); }

}  // namespace

TEST_CASE("sweep: completed run, matching one-point solves, with its record") {
    const auto d = diode();
    const auto points = forward_points();
    const auto sweep = solve::sweep_bias(d, points);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    REQUIRE_FALSE(sweep->unfinished.has_value());
    REQUIRE(sweep->points.size() == points.size());
    // The same chain of one-point solves, each started from the previous point's fields.
    std::optional<results::NodeFields> state;
    for (std::size_t k = 0; k < points.size(); ++k) {
        auto one = solve::solve_bias(d, points[k], {}, state ? &*state : nullptr);
        REQUIRE(one.has_value());
        CAPTURE(k);
        REQUIRE(close_point(*one, sweep->points[k]));
        REQUIRE(sweep->points[k].convergence.converged);
        state = std::move(one->fields);
    }
    REQUIRE(std::abs(sweep->points[5].terminal_current[0] - 1.280e-2) / 1.280e-2 < 1e-2);
    REQUIRE(sweep->run.input_identity != 0);
    bool found = false;
    for (const auto& [name, value] : sweep->run.settings) {
        if (name == "newton.tol_update") found = value == 1e-8;
    }
    REQUIRE(found);
}

TEST_CASE("sweep: progress is monotonic and matches the convergence records") {
    const auto d = diode();
    const auto points = forward_points();
    std::vector<Progress> events;
    const solve::RunControl control{.stop = {},
                                    .progress = [&](const Progress& e) { events.push_back(e); }};
    const auto sweep = solve::sweep_bias(d, points, {}, nullptr, control);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    REQUIRE_FALSE(events.empty());
    REQUIRE(events.front().phase == Phase::equilibrium);
    REQUIRE(events.front().iteration == 1);
    for (std::size_t k = 1; k < events.size(); ++k) {
        CAPTURE(k);
        REQUIRE(key(events[k - 1]) < key(events[k]));  // strictly increasing
        // Within a point the iterations count up by one; a new point starts at 1, after the
        // previous point's converged event.
        if (events[k].phase == events[k - 1].phase && events[k].point == events[k - 1].point) {
            REQUIRE(events[k].iteration == events[k - 1].iteration + 1);
            REQUIRE_FALSE(events[k - 1].converged);
        } else {
            REQUIRE(events[k].iteration == 1);
            REQUIRE(events[k - 1].converged);
        }
    }
    REQUIRE(events.back().converged);
    // The bias events are exactly the points' convergence records.
    std::size_t at = 0;
    while (events[at].phase == Phase::equilibrium) ++at;
    for (std::size_t p = 0; p < points.size(); ++p) {
        for (const results::IterationRecord& r : sweep->points[p].convergence.iterations) {
            const Progress& e = events[at++];
            REQUIRE(e.point == p);
            REQUIRE(e.point_count == points.size());
            REQUIRE(e.iteration == r.iteration);
            REQUIRE(e.update == r.update);
            REQUIRE(e.residual == r.residual);
        }
    }
    REQUIRE(at == events.size());
}

TEST_CASE("cancel: between bias points, the completed points are kept") {
    const auto d = diode();
    const auto points = forward_points();
    const auto full = *solve::sweep_bias(d, points);
    std::stop_source source;
    const auto stop_after_point_1 = [&](const Progress& e) {
        if (e.phase == Phase::bias && e.point == 1 && e.converged) source.request_stop();
    };
    const auto sweep = solve::sweep_bias(
        d, points, {}, nullptr, {.stop = source.get_token(), .progress = stop_after_point_1});
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->stopped.has_value());
    REQUIRE(sweep->stopped->code == ErrorCode::cancelled);
    REQUIRE(sweep->stopped->context->index == 2);  // stopped before point 2
    REQUIRE_FALSE(sweep->unfinished.has_value());
    REQUIRE(sweep->points.size() == 2);
    REQUIRE(same_point(sweep->points[0], full.points[0]));
    REQUIRE(same_point(sweep->points[1], full.points[1]));
    REQUIRE(sweep->run.input_identity == full.run.input_identity);
}

TEST_CASE("cancel: inside a point, at the next Newton iteration") {
    const auto d = diode();
    const auto points = forward_points();
    std::stop_source source;
    const auto stop_in_point_3 = [&](const Progress& e) {
        if (e.phase == Phase::bias && e.point == 3 && e.iteration == 2) source.request_stop();
    };
    const auto sweep = solve::sweep_bias(
        d, points, {}, nullptr, {.stop = source.get_token(), .progress = stop_in_point_3});
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->stopped->code == ErrorCode::cancelled);
    REQUIRE(sweep->stopped->context->index == 3);  // before Newton iteration 3
    REQUIRE(sweep->points.size() == 3);
    REQUIRE(sweep->unfinished.has_value());
    REQUIRE(sweep->unfinished->iterations.size() == 2);
    REQUIRE_FALSE(sweep->unfinished->converged);
}

TEST_CASE("cancel: before the start nothing is solved") {
    const auto d = diode();
    std::stop_source source;
    source.request_stop();
    int events = 0;
    const solve::RunControl control{.stop = source.get_token(),
                                    .progress = [&](const Progress&) { ++events; }};
    const auto sweep = solve::sweep_bias(d, forward_points(), {}, nullptr, control);
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->stopped->code == ErrorCode::cancelled);
    REQUIRE(sweep->points.empty());
    REQUIRE(events == 0);
    REQUIRE(solve::solve_equilibrium(d, {}, control).error().code == ErrorCode::cancelled);
    const std::vector<double> bias{0.5, 0.0};
    REQUIRE(solve::solve_bias(d, bias, {}, nullptr, control).error().code == ErrorCode::cancelled);
}

TEST_CASE("cancel: from another thread (6.9: the solve runs on a worker)") {
    // The worker pauses in its first unconverged bias event (point 1, iteration 1; the 0 V point
    // converges in one iteration) until the main thread has requested the stop, so the outcome does
    // not depend on timing: the next safe point is Newton iteration 2 of point 1.
    const auto d = diode();
    const auto points = forward_points();
    std::stop_source source;
    std::latch reached(1);
    std::binary_semaphore resume(0);
    std::optional<std::expected<results::Sweep, base::Error>> outcome;
    {
        std::jthread worker([&] {
            bool first = true;
            const auto pause_once = [&](const Progress& e) {
                if (e.phase == Phase::bias && !e.converged && first) {
                    first = false;
                    reached.count_down();
                    resume.acquire();
                }
            };
            outcome = solve::sweep_bias(d, points, {}, nullptr,
                                        {.stop = source.get_token(), .progress = pause_once});
        });
        reached.wait();
        source.request_stop();
        resume.release();
    }  // joins
    REQUIRE(outcome.has_value());
    REQUIRE(outcome->has_value());
    const results::Sweep& sweep = **outcome;
    REQUIRE(sweep.stopped->code == ErrorCode::cancelled);
    REQUIRE(sweep.stopped->context->index == 2);  // before Newton iteration 2
    REQUIRE(sweep.points.size() == 1);
    REQUIRE(sweep.unfinished->iterations.size() == 1);
}

TEST_CASE("sweep: a point that does not converge stops the run, keeping the earlier points") {
    // With at most 8 iterations the equilibrium (7) and the 0 V point (1) converge, but the jump
    // to 0.8 V does not.
    const auto d = diode();
    const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.8, 0.0}, {0.9, 0.0}};
    const auto sweep = solve::sweep_bias(d, points, {.newton = {.max_iterations = 8}});
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->stopped->code == ErrorCode::non_convergence);
    REQUIRE(sweep->points.size() == 1);
    REQUIRE(sweep->unfinished->iterations.size() == 8);
}

TEST_CASE("sweep: invalid input is an error before anything is solved") {
    const auto d = diode();
    int events = 0;
    const solve::RunControl control{.stop = {}, .progress = [&](const Progress&) { ++events; }};
    const std::vector<std::vector<double>> none;
    REQUIRE(solve::sweep_bias(d, none, {}, nullptr, control).error().code ==
            ErrorCode::invalid_input);
    const std::vector<std::vector<double>> short_point{{0.1, 0.0}, {0.2}};
    const auto e = solve::sweep_bias(d, short_point, {}, nullptr, control);
    REQUIRE(e.error().code == ErrorCode::invalid_input);
    REQUIRE(e.error().context->index == 1);
    const std::vector<std::vector<double>> nan_point{{std::nan(""), 0.0}};
    REQUIRE(solve::sweep_bias(d, nan_point, {}, nullptr, control).error().code ==
            ErrorCode::invalid_input);
    REQUIRE(events == 0);
}

TEST_CASE("run record: the identity follows every input") {
    const auto d = diode();
    const auto points = forward_points();
    const solve::BiasOptions options;
    const auto base_id = solve::make_run_record(d, options, points).input_identity;
    REQUIRE(solve::make_run_record(diode(), options, points).input_identity == base_id);
    REQUIRE(solve::make_run_record(diode(1.01), options, points).input_identity != base_id);
    solve::BiasOptions tighter = options;
    tighter.newton.tol_update = 1e-9;
    REQUIRE(solve::make_run_record(d, tighter, points).input_identity != base_id);
    auto other_points = points;
    other_points.back()[0] = 0.55;
    REQUIRE(solve::make_run_record(d, options, other_points).input_identity != base_id);
    const auto eq = *solve::solve_equilibrium(d);
    REQUIRE(solve::make_run_record(d, options, points, &eq.fields).input_identity != base_id);
    solve::BiasOptions scaled = options;
    scaled.Ns_override = 1e16;
    const auto record = solve::make_run_record(d, scaled, points);
    bool found = false;
    for (const auto& [name, value] : record.settings) {
        if (name == "scaling.Ns_override") found = value == 1e16;
    }
    REQUIRE(found);
}
