// Bias step cutting (Unit 28, solve::StepControl): a point that does not converge from the previous
// state is reached through shorter steps, the result is the point's own steady state, and with the
// control off (or never needed) nothing changes.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stop_token>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/control.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;
using base::ErrorCode;

namespace {

// The Unit 10 sweep fixture: the legacy 1e17/1e17 silicon diode.
device::Device diode() {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
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

double relative(const std::vector<double>& u, const std::vector<double>& v) {
    REQUIRE(u.size() == v.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        worst = std::max(worst, std::abs(u[i] - v[i]) /
                                    std::max({std::abs(u[i]), std::abs(v[i]), 1e-300}));
    }
    return worst;
}

}  // namespace

TEST_CASE("step control: a step too long for Newton is cut and reaches the point's steady state") {
    // 8 iterations do not take the equilibrium state to 0.8 V (the Unit 10 failing-point fixture).
    const auto d = diode();
    const std::vector<std::vector<double>> jump{{0.8, 0.0}};
    solve::BiasOptions options{.newton = {.max_iterations = 8, .tol_update = 1e-12}};
    const auto plain = solve::sweep_bias(d, jump, options);
    REQUIRE(plain.has_value());
    REQUIRE(plain->stopped->code == ErrorCode::non_convergence);

    options.step_control.max_cuts = 10;
    const auto cut = solve::sweep_bias(d, jump, options);
    REQUIRE(cut.has_value());
    REQUIRE_FALSE(cut->stopped.has_value());
    REQUIRE(cut->points.size() == 1);
    const results::BiasPoint& p = cut->points.front();
    REQUIRE(p.bias_V == jump.front());  // the requested bias, exactly
    REQUIRE(p.bias_substeps > 0);
    REQUIRE(p.convergence.converged);

    // The same steady state as a fine sweep reaching 0.8 V with no cutting.
    std::vector<std::vector<double>> fine;
    for (int k = 1; k <= 16; ++k) fine.push_back({0.05 * k, 0.0});
    const auto reference =
        solve::sweep_bias(d, fine, {.newton = {.tol_update = 1e-12}});
    REQUIRE(reference.has_value());
    REQUIRE_FALSE(reference->stopped.has_value());
    const results::BiasPoint& r = reference->points.back();
    REQUIRE(r.bias_V == p.bias_V);
    CHECK(relative(p.fields.potential_V, r.fields.potential_V) < 1e-9);
    CHECK(relative(p.fields.n_cm3, r.fields.n_cm3) < 1e-9);
    CHECK(relative(p.fields.p_cm3, r.fields.p_cm3) < 1e-9);
    CHECK(relative(p.terminal_current, r.terminal_current) < 1e-9);
}

TEST_CASE("step control: off, or never needed, changes nothing") {
    const auto d = diode();
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 5; ++k) points.push_back({0.1 * k, 0.0});
    const auto off = solve::sweep_bias(d, points);
    solve::BiasOptions on;
    on.step_control.max_cuts = 20;
    const auto with = solve::sweep_bias(d, points, on);
    REQUIRE(off.has_value());
    REQUIRE(with.has_value());
    REQUIRE(off->points.size() == with->points.size());
    for (std::size_t k = 0; k < points.size(); ++k) {
        const auto& a = off->points[k];
        const auto& b = with->points[k];
        CHECK(b.bias_substeps == 0);
        CHECK(a.fields.potential_V == b.fields.potential_V);
        CHECK(a.fields.n_cm3 == b.fields.n_cm3);
        CHECK(a.fields.p_cm3 == b.fields.p_cm3);
        CHECK(a.terminal_current == b.terminal_current);
        CHECK(a.convergence.iterations.size() == b.convergence.iterations.size());
    }
    // The identity of a run without it is the one it always had; with it, it differs.
    CHECK(off->run.input_identity != with->run.input_identity);
    const auto named = [](const results::RunRecord& r, const std::string& name) {
        return std::ranges::any_of(r.settings, [&](const auto& s) { return s.first == name; });
    };
    CHECK_FALSE(named(off->run, "step_control.max_cuts"));
    CHECK(named(with->run, "step_control.max_cuts"));
}

TEST_CASE("step control: a step that fails at the smallest cut stops the run, keeping points") {
    const auto d = diode();
    const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.8, 0.0}};
    solve::BiasOptions options{.newton = {.max_iterations = 8}};
    options.step_control.max_cuts = 1;  // steps down to 0.4 V
    const auto sweep = solve::sweep_bias(d, points, options);
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->stopped.has_value());
    CHECK(sweep->stopped->code == ErrorCode::non_convergence);
    CHECK(sweep->stopped->message.find("2^-1") != std::string::npos);
    CHECK(sweep->points.size() == 1);
    REQUIRE(sweep->unfinished.has_value());
    CHECK(sweep->unfinished->iterations.size() == 8);
}

TEST_CASE("step control: progress stays strictly increasing across the cut attempts") {
    const auto d = diode();
    const std::vector<std::vector<double>> points{{0.2, 0.0}, {0.8, 0.0}};
    std::vector<solve::Progress> events;
    const solve::RunControl control{.stop = {},
                                    .progress = [&](const solve::Progress& e) {
                                        events.push_back(e);
                                    }};
    solve::BiasOptions options{.newton = {.max_iterations = 8}};
    options.step_control.max_cuts = 10;
    const auto sweep = solve::sweep_bias(d, points, options, nullptr, control);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    REQUIRE(sweep->points.back().bias_substeps > 0);
    const auto key = [](const solve::Progress& e) {
        return std::tuple(static_cast<int>(e.phase), e.point, e.iteration);
    };
    for (std::size_t i = 1; i < events.size(); ++i) CHECK(key(events[i - 1]) < key(events[i]));
    // The last event closes the last point.
    CHECK(events.back().converged);
    CHECK(events.back().point == 1);
}

TEST_CASE("step control: cancellation during the cut attempts keeps the completed points") {
    const auto d = diode();
    const std::vector<std::vector<double>> points{{0.2, 0.0}, {0.8, 0.0}};
    std::stop_source source;
    const solve::RunControl control{
        .stop = source.get_token(), .progress = [&](const solve::Progress& e) {
            // Past the 8 iterations of the failed full step: inside the cut attempts.
            if (e.phase == solve::Phase::bias && e.point == 1 && e.iteration > 10) {
                source.request_stop();
            }
        }};
    solve::BiasOptions options{.newton = {.max_iterations = 8}};
    options.step_control.max_cuts = 10;
    const auto sweep = solve::sweep_bias(d, points, options, nullptr, control);
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->stopped.has_value());
    CHECK(sweep->stopped->code == ErrorCode::cancelled);
    CHECK(sweep->points.size() == 1);
}

TEST_CASE("step control: a given initial state's first point is not cut") {
    const auto d = diode();
    const std::vector<std::vector<double>> zero{{0.0, 0.0}};
    const auto start = solve::sweep_bias(d, zero);
    REQUIRE(start.has_value());
    const std::vector<std::vector<double>> jump{{0.8, 0.0}};
    solve::BiasOptions options{.newton = {.max_iterations = 8}};
    options.step_control.max_cuts = 10;
    const auto sweep = solve::sweep_bias(d, jump, options, &start->points.front().fields);
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->stopped.has_value());
    CHECK(sweep->stopped->code == ErrorCode::non_convergence);
}

TEST_CASE("step control: invalid options are errors before anything is solved") {
    const auto d = diode();
    const std::vector<std::vector<double>> points{{0.1, 0.0}};
    for (const solve::StepControl s :
         {solve::StepControl{.max_cuts = -1}, solve::StepControl{.max_cuts = 53},
          solve::StepControl{.max_cuts = 4, .growth = 0.5},
          solve::StepControl{.max_cuts = 4, .max_iterations = 0},
          solve::StepControl{.max_cuts = 4, .growth = std::nan("")}}) {
        solve::BiasOptions options;
        options.step_control = s;
        const auto e = solve::sweep_bias(d, points, options);
        REQUIRE_FALSE(e.has_value());
        CHECK(e.error().code == ErrorCode::invalid_input);
    }
}
