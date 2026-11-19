// `initial = <checkpoint>` (Unit 61, RD2; gate RG3): a transient, small-signal, trace and sweep run
// from a checkpoint of its own starting state equals the run without it to 1e-10 relative, with no
// more Newton iterations at the first state; refused with restart and on another mesh.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_device;
using deck_test::must;
using deck_test::plan_of;
using deck_test::scratch;
namespace fs = std::filesystem;

namespace {

// |a - b| <= 1e-10 max(|a|, |b|) + 1e-10 of the field's largest magnitude (its own scale).
void close(std::span<const double> a, std::span<const double> b, const char* what) {
    INFO(what);
    REQUIRE(a.size() == b.size());
    double scale = 0.0;
    for (const double x : a) scale = std::max(scale, std::abs(x));
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double allowed = 1e-10 * std::max(std::abs(a[i]), std::abs(b[i])) + 1e-10 * scale;
        INFO(i << ": " << a[i] << " against " << b[i]);
        CHECK(std::abs(a[i] - b[i]) <= allowed);
    }
}

void close(const results::NodeFields& a, const results::NodeFields& b) {
    close(a.potential_V, b.potential_V, "potential");
    close(a.n_cm3, b.n_cm3, "n");
    close(a.p_cm3, b.p_cm3, "p");
}

void close(const results::BiasPoint& a, const results::BiasPoint& b) {
    close(a.fields, b.fields);
    close(a.terminal_current, b.terminal_current, "terminal current");
}

// The run of `deck` in `dir`, and the Newton iterations of its first state: every progress event
// before the first one of a later bias point or a later time.
struct Counted {
    deck::RunReport report;
    int first_iterations = 0;
};

Counted run_counted(const std::string& text, const fs::path& dir) {
    const deck::Plan plan = must(text, dir);
    Counted out;
    bool past_first = false;
    solve::RunControl control;
    control.progress = [&](const solve::Progress& p) {
        if (p.point > 0 || p.time_s > 0.0 || p.frequency_Hz > 0.0) past_first = true;
        if (!past_first) ++out.first_iterations;
    };
    auto r = deck::run(plan, control);
    if (!r) FAIL(r.error().message);
    out.report = std::move(*r);
    return out;
}

// The checkpoint of the steady state at anode = 0.4 V, written by a one-point sweep.
fs::path start_checkpoint(const fs::path& dir) {
    const auto r = deck::run(must(diode_device + "[sweep]\ncontact anode = 0.4\n[output]\ncheckpoint = start.nitckpt\n", dir));
    if (!r) FAIL(r.error().message);
    return dir / "start.nitckpt";
}

}  // namespace

TEST_CASE("deck initial: each analysis from a checkpoint of its start equals the run without it (RG3)",
          "[deck][initial]") {
    const fs::path dir = scratch("initial_runs");
    start_checkpoint(dir);
    const std::string with = "initial = start.nitckpt\n";

    SECTION("sweep") {
        const std::string body = "[sweep]\ncontact anode = 0.4 : 0.6 : 0.1\n";
        const auto a = run_counted(diode_device + body, dir);
        const auto b = run_counted(diode_device + body + with, dir);
        REQUIRE(a.report.sweep.has_value());
        REQUIRE(b.report.sweep.has_value());
        REQUIRE(a.report.sweep->points.size() == b.report.sweep->points.size());
        for (std::size_t k = 0; k < a.report.sweep->points.size(); ++k) close(a.report.sweep->points[k], b.report.sweep->points[k]);
        INFO(a.first_iterations << " against " << b.first_iterations);
        CHECK(b.first_iterations <= a.first_iterations);
        CHECK(b.first_iterations > 0);
    }
    SECTION("small signal") {
        const std::string body = "[small_signal]\ncontact anode = 0.4, 0.5\nfrequencies = 1MHz, 1GHz\n";
        const auto a = run_counted(diode_device + body, dir);
        const auto b = run_counted(diode_device + body + with, dir);
        REQUIRE(a.report.small_signal.has_value());
        REQUIRE(b.report.small_signal.has_value());
        const auto& pa = a.report.small_signal->points;
        const auto& pb = b.report.small_signal->points;
        REQUIRE(pa.size() == pb.size());
        for (std::size_t k = 0; k < pa.size(); ++k) {
            close(pa[k].dc, pb[k].dc);
            REQUIRE(pa[k].admittance.size() == pb[k].admittance.size());
            for (std::size_t q = 0; q < pa[k].admittance.size(); ++q) {
                const auto& ya = pa[k].admittance[q];
                const auto& yb = pb[k].admittance[q];
                REQUIRE(ya.size() == yb.size());
                double scale = 0.0;
                for (const auto& y : ya) scale = std::max(scale, std::abs(y));
                for (std::size_t i = 0; i < ya.size(); ++i) CHECK(std::abs(ya[i] - yb[i]) <= 1e-10 * scale);
            }
        }
        CHECK(b.first_iterations <= a.first_iterations);
    }
    SECTION("trace") {
        const std::string body = "[trace]\ncontact = anode\nstart anode = 0.4\nend = 0.6\n";
        const auto a = run_counted(diode_device + body, dir);
        const auto b = run_counted(diode_device + body + with, dir);
        REQUIRE(a.report.sweep.has_value());
        REQUIRE(b.report.sweep.has_value());
        REQUIRE(a.report.sweep->points.size() == b.report.sweep->points.size());
        for (std::size_t k = 0; k < a.report.sweep->points.size(); ++k) close(a.report.sweep->points[k], b.report.sweep->points[k]);
        CHECK(b.first_iterations <= a.first_iterations);
    }
    SECTION("transient") {
        // Fixed steps, so the comparison is of the same time points.
        const std::string body = "[transient]\ncontact anode = step(0.4, 0.5, 100ps)\nt_end = 300ps\n"
                                 "dt_initial = 10ps\nadaptive = off\noutput_times = 100ps, 200ps\n";
        const auto a = run_counted(diode_device + body, dir);
        const auto b = run_counted(diode_device + body + with, dir);
        REQUIRE(a.report.transient.has_value());
        REQUIRE(b.report.transient.has_value());
        const auto& ta = *a.report.transient;
        const auto& tb = *b.report.transient;
        REQUIRE(ta.points.size() == tb.points.size());
        for (std::size_t k = 0; k < ta.points.size(); ++k) {
            CHECK(ta.points[k].time_s == tb.points[k].time_s);
            close(ta.points[k].terminal_current, tb.points[k].terminal_current, "terminal current");
        }
        REQUIRE(ta.snapshots.size() == tb.snapshots.size());
        for (std::size_t k = 0; k < ta.snapshots.size(); ++k) close(ta.snapshots[k].fields, tb.snapshots[k].fields);
        CHECK(b.first_iterations <= a.first_iterations);
    }
}

TEST_CASE("deck initial: refused with restart, on another mesh, and in [adapt] (RG3)", "[deck][initial]") {
    const fs::path dir = scratch("initial_refusals");
    start_checkpoint(dir);
    const std::string d = diode_device;
    const auto lines = static_cast<int>(std::count(d.begin(), d.end(), '\n'));
    const auto at = [lines](int k) { return "t.deck:" + std::to_string(lines + 1 + k) + ": "; };

    const auto both = plan_of(d + "[sweep]\ncontact anode = 0.4 : 0.6 : 0.1\nrestart = start.nitckpt\ninitial = start.nitckpt\n", dir);
    REQUIRE_FALSE(both.has_value());
    INFO(both.error().message);
    CHECK(both.error().message.find(at(3) + "initial and restart together") != std::string::npos);

    // A checkpoint of a 101-node diode, given to the 201-node one.
    std::string coarse = d;
    coarse.replace(coarse.find("linspace(0, 2um, 201)"), 21, "linspace(0, 2um, 101)");
    const auto c = deck::run(must(coarse + "[sweep]\ncontact anode = 0.4\n[output]\ncheckpoint = coarse.nitckpt\n", dir));
    REQUIRE(c.has_value());
    for (const std::string section : {"[sweep]\ncontact anode = 0.4\n",
                                      "[small_signal]\ncontact anode = 0.4\nfrequencies = 1e6\n",
                                      "[trace]\ncontact = anode\nstart anode = 0.4\nend = 0.6\n",
                                      "[transient]\ncontact anode = constant(0.4)\nt_end = 1e-10\n"}) {
        INFO(section);
        const auto r = deck::run(must(d + section + "initial = coarse.nitckpt\n", dir));
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == base::ErrorCode::invalid_input);
    }

    const auto missing = deck::run(must(d + "[sweep]\ncontact anode = 0.4\ninitial = none.nitckpt\n", dir));
    CHECK_FALSE(missing.has_value());
}
