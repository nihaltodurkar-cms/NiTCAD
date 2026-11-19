// Trace resume (Unit 61, RD5; gate RG5, solver part): a trace stopped part-way and continued from
// its resume state equals the uninterrupted trace bit for bit from its last accepted point on: on
// the diode, and on the electrothermal runaway fold (Unit 36's case) at points before, at and past
// the fold. A stop lands either right after an accepted point or after an attempt the trace then
// rejects; both resume from the last accepted point. Another run's state, a finished one, one with
// an initial state are refused.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cstddef>
#include <cstdint>
#include <stop_token>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
bool same(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!same(a[i], b[i])) return false;
    return true;
}

void same_point(const results::BiasPoint& a, const results::BiasPoint& b) {
    CHECK(same(a.bias_V, b.bias_V));
    CHECK(same(a.terminal_current, b.terminal_current));
    CHECK(same(a.fields.potential_V, b.fields.potential_V));
    CHECK(same(a.fields.n_cm3, b.fields.n_cm3));
    CHECK(same(a.fields.p_cm3, b.fields.p_cm3));
    CHECK(same(a.fields.temperature_K, b.fields.temperature_K));
    CHECK(a.convergence.iterations.size() == b.convergence.iterations.size());
}

// The Unit 9 diode (anode, contact 0, on x_min).
device::Device diode() {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, {1e-4}, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

// Unit 36's runaway fold: a 1 um pn bar, kappa 1e3 times silicon's, R_th 2 K cm^2/W to 300 K at
// each end (electrothermal_test.cpp).
device::Device fold_bar() {
    std::vector<double> x(41);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = 1e-4 * static_cast<double>(i) / 40.0;
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n), acceptors(n);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 0.5e-4 ? acceptors : donors)[i] = 1e17;
    physics::SemiconductorParameters material = physics::silicon_parameters;
    material.thermal.conductivity_W_cmK *= 1e3;
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(material)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}},
         .thermal_contacts = {{"left", "x_min", device::ThermalContactKind::resistance, 300.0, 2.0},
                              {"right", "x_max", device::ThermalContactKind::resistance, 300.0, 2.0}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

// Converged corrector attempts of a trace (the starting point's sweep excluded).
std::size_t converged_attempts(const device::Device& d, const solve::TraceOptions& o) {
    std::size_t count = 0;
    solve::RunControl control;
    control.progress = [&](const solve::Progress& p) {
        if (p.phase == solve::Phase::bias && p.point_count == 0 && p.converged) ++count;
    };
    REQUIRE(solve::trace_bias(d, o, nullptr, control).has_value());
    return count;
}

// The trace stopped at the j-th converged attempt (j = 0: right after the starting point).
results::Sweep stopped_at(const device::Device& d, const solve::TraceOptions& o, std::size_t j) {
    std::stop_source source;
    std::size_t count = 0;
    solve::RunControl control{.stop = source.get_token()};
    control.progress = [&](const solve::Progress& p) {
        if (p.phase != solve::Phase::bias || !p.converged) return;
        if (p.point_count == 0) ++count;
        if (count == j) source.request_stop();
    };
    auto r = solve::trace_bias(d, o, nullptr, control);
    REQUIRE(r.has_value());
    REQUIRE(r->stopped.has_value());
    REQUIRE(r->resume.has_value());
    return std::move(*r);
}

// Stops `o` at each of the attempts `at`, continues and compares with the uninterrupted trace.
void joins(const device::Device& d, const solve::TraceOptions& o, const std::vector<std::size_t>& at) {
    const auto full = solve::trace_bias(d, o);
    REQUIRE(full.has_value());
    REQUIRE(full->resume.has_value());
    CHECK(full->resume->finished == !full->stopped.has_value());
    for (const std::size_t j : at) {
        INFO("stopped at converged attempt " << j);
        const auto first = stopped_at(d, o, j);
        const std::size_t k = first.resume->points - 1;
        REQUIRE(k + 1 == first.points.size());
        same_point(first.resume->point, full->points[k]);
        const auto rest = solve::trace_bias(d, o, nullptr, {}, &*first.resume);
        if (!rest) FAIL(rest.error().message);
        CHECK(rest->run.input_identity == full->run.input_identity);
        REQUIRE(rest->points.size() == full->points.size() - k);
        for (std::size_t i = 0; i < rest->points.size(); ++i) {
            INFO("point " << k + i);
            same_point(rest->points[i], full->points[k + i]);
        }
        CHECK(rest->stopped.has_value() == full->stopped.has_value());
        if (rest->stopped) CHECK(rest->stopped->message == full->stopped->message);
        REQUIRE(rest->resume.has_value());
        CHECK(rest->resume->points == full->resume->points);
        CHECK(rest->resume->attempt == full->resume->attempt);
        CHECK(same(rest->resume->step, full->resume->step));
    }
}

solve::TraceOptions diode_trace() {
    solve::TraceOptions o;
    o.contact = 0;
    o.start_V = {0.0, 0.0};
    o.end_V = 0.9;
    o.step_V = 0.05;
    o.max_step_V = 0.1;
    return o;
}

}  // namespace

TEST_CASE("trace resume: the diode continues bit for bit (RG5)", "[solve][resume]") {
    const device::Device d = diode();
    const solve::TraceOptions o = diode_trace();
    const std::size_t attempts = converged_attempts(d, o);
    REQUIRE(attempts > 6);
    joins(d, o, {0, 1, 2, attempts / 2, attempts - 1});
}

TEST_CASE("trace resume: the electrothermal fold continues bit for bit, before, at and past it (RG5)",
          "[solve][resume]") {
    const device::Device d = fold_bar();
    solve::TraceOptions o;
    o.steady.models.electrothermal = true;
    o.steady.newton.tol_update = 1e-12;
    o.contact = 0;
    o.start_V = {0.5, 0.0};
    o.end_V = 1.0;
    o.step_V = 0.02;
    o.max_step_V = 0.02;
    const auto full = solve::trace_bias(d, o);
    REQUIRE(full.has_value());
    std::size_t top = 0;
    for (std::size_t k = 0; k < full->points.size(); ++k)
        if (full->points[k].bias_V[0] > full->points[top].bias_V[0]) top = k;
    REQUIRE(top > 0);
    REQUIRE(top + 2 < full->points.size());
    // Converged attempts up to and past the fold's point (attempts >= accepted points).
    const std::size_t attempts = converged_attempts(d, o);
    REQUIRE(attempts >= full->points.size() - 1);
    joins(d, o, {1, top - 1, top, top + 1, attempts - 1});
}

TEST_CASE("trace resume: refusals (RG5)", "[solve][resume]") {
    const device::Device d = diode();
    const solve::TraceOptions o = diode_trace();
    const auto part = stopped_at(d, o, 3);
    const auto refused = [&](const solve::TraceOptions& options, const results::TraceResume& r,
                             const std::string& what, const results::NodeFields* initial = nullptr) {
        const auto t = solve::trace_bias(d, options, initial, {}, &r);
        REQUIRE_FALSE(t.has_value());
        INFO(t.error().message);
        CHECK(t.error().message.find(what) != std::string::npos);
    };
    solve::TraceOptions other = o;
    other.end_V = 0.8;
    refused(other, *part.resume, "belongs to another run");
    refused(o, *part.resume, "takes no initial state", &part.points.front().fields);
    auto finished = solve::trace_bias(d, o);
    REQUIRE(finished.has_value());
    REQUIRE(finished->resume->finished);
    refused(o, *finished->resume, "nothing is left");
    results::TraceResume cut = *part.resume;
    cut.tangent.pop_back();
    refused(o, cut, "does not fit");
}
