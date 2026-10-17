// Running [transient], [small_signal] and [trace] decks (Unit 43; DECISIONS.md AN5, gates AG2 and
// AG4): each run equals the C++ call with the options the deck states, on the hand-built diode, bit
// for bit and with the same run identity; the files it writes load equal to its results, and its
// table has the stated columns and one row per step, (point, frequency) or point.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <complex>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_device;
using deck_test::hand_diode;
using deck_test::must;
using deck_test::scratch;
namespace fs = std::filesystem;

namespace {

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
    }
    return true;
}

bool same_bits(const std::vector<std::complex<double>>& a, const std::vector<std::complex<double>>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i].real()) != std::bit_cast<std::uint64_t>(b[i].real()) ||
            std::bit_cast<std::uint64_t>(a[i].imag()) != std::bit_cast<std::uint64_t>(b[i].imag())) {
            return false;
        }
    }
    return true;
}

std::string text_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

long lines_of(const std::string& s) { return static_cast<long>(std::count(s.begin(), s.end(), '\n')); }

}  // namespace

TEST_CASE("deck run: a [transient] deck equals solve_transient; its file and table", "[deck]") {
    const fs::path dir = scratch("run_transient");
    const deck::Plan plan = must(diode_device + R"(
[transient]
contact anode = ramp(0, 0.3, 0, 2e-10)
t_end = 4e-10
dt_initial = 1e-10
adaptive = off
integrator = backward_euler
output_times = 2e-10

[output]
transient = out.nittran
table = tr.csv
)", dir);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->transient.has_value());
    const results::Transient& t = *report->transient;
    REQUIRE_FALSE(t.stopped.has_value());

    solve::TransientOptions o;
    o.t_end_s = 4e-10;
    o.dt_initial_s = 1e-10;
    o.adaptive = false;
    o.integrator = solve::Integrator::backward_euler;
    o.output_times_s = {2e-10};
    const std::vector<solve::Waveform> w{*solve::Waveform::ramp(0, 0.3, 0, 2e-10),
                                         solve::Waveform::constant(0.0)};
    const auto hand = solve::solve_transient(hand_diode(), w, o);
    REQUIRE(hand.has_value());
    CHECK(t.run.input_identity == hand->run.input_identity);
    REQUIRE(t.points.size() == hand->points.size());
    REQUIRE(t.points.size() == 5);  // the start and four fixed steps
    for (std::size_t k = 0; k < t.points.size(); ++k) {
        CHECK(t.points[k].time_s == hand->points[k].time_s);
        CHECK(same_bits(t.points[k].terminal_current, hand->points[k].terminal_current));
        CHECK(same_bits(t.points[k].displacement_current, hand->points[k].displacement_current));
    }
    REQUIRE(t.snapshots.size() == hand->snapshots.size());
    CHECK(same_bits(t.snapshots.back().fields.n_cm3, hand->snapshots.back().fields.n_cm3));

    CHECK(report->written == std::vector<fs::path>{dir / "out.nittran", dir / "tr.csv"});
    const auto loaded = io::load_transient(dir / "out.nittran");
    REQUIRE(loaded.has_value());
    CHECK(loaded->run.input_identity == t.run.input_identity);
    CHECK(loaded->run.provenance == t.run.provenance);
    CHECK(same_bits(loaded->points.back().terminal_current, t.points.back().terminal_current));
    CHECK(same_bits(loaded->snapshots.back().fields.potential_V, t.snapshots.back().fields.potential_V));

    const std::string table = text_of(dir / "tr.csv");
    CHECK(table.starts_with("time s,step s,order,anode V,cathode V,anode A,cathode A,anode conduction A,"
                            "cathode conduction A,anode displacement A,cathode displacement A\n0,0,0,0,0,"));
    CHECK(lines_of(table) == 1 + 5);
    // The last step, at t_end (the solver's sum of steps), 17 digits.
    CHECK(table.find("\n4.0000000000000001e-10,1.0000000000000002e-10,1,0.29999999999999999,0,") !=
          std::string::npos);
    CHECK_FALSE(fs::exists(dir / "tr.csv.tmp"));
}

TEST_CASE("deck run: a [small_signal] deck equals solve_small_signal; its files and table", "[deck]") {
    const fs::path dir = scratch("run_small_signal");
    const deck::Plan plan = must(diode_device + R"(
[small_signal]
contact anode = -0.2, 0
frequencies = logspace(1e3, 1e9, 3)
field_frequencies = 1e6

[output]
small_signal = out.nitac
checkpoint = out.nitckpt
table = ac.csv
)", dir);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->small_signal.has_value());
    const results::SmallSignal& s = *report->small_signal;
    REQUIRE_FALSE(s.stopped.has_value());

    const std::vector<std::vector<double>> points{{-0.2, 0.0}, {0.0, 0.0}};
    const auto hand = solve::solve_small_signal(
        hand_diode(), points, {.frequencies_Hz = plan.small_signal.frequencies_Hz, .field_frequencies_Hz = {1e6}});
    REQUIRE(hand.has_value());
    CHECK(plan.small_signal.frequencies_Hz == std::vector<double>{1e3, 1e6, 1e9});
    CHECK(s.run.input_identity == hand->run.input_identity);
    REQUIRE(s.points.size() == 2);
    for (std::size_t k = 0; k < 2; ++k) {
        for (std::size_t f = 0; f < 3; ++f) CHECK(same_bits(s.points[k].admittance[f], hand->points[k].admittance[f]));
        REQUIRE(s.points[k].fields.size() == 2);  // at 1 MHz, one per driven contact
        for (std::size_t f = 0; f < 2; ++f) {
            CHECK(s.points[k].fields[f].contact == hand->points[k].fields[f].contact);
            CHECK(same_bits(s.points[k].fields[f].potential, hand->points[k].fields[f].potential));
        }
    }

    CHECK(report->written == std::vector<fs::path>{dir / "out.nitac", dir / "ac.csv", dir / "out.nitckpt"});
    const auto loaded = io::load_small_signal(dir / "out.nitac");
    REQUIRE(loaded.has_value());
    CHECK(loaded->run.input_identity == s.run.input_identity);
    CHECK(loaded->run.provenance == s.run.provenance);
    CHECK(same_bits(loaded->points[1].admittance[2], s.points[1].admittance[2]));
    const auto c = io::load_checkpoint(dir / "out.nitckpt", plan.device.mesh());
    REQUIRE(c.has_value());
    CHECK(c->bias_V == std::vector<double>{0.0, 0.0});
    CHECK(c->run_identity == s.run.input_identity);
    CHECK(same_bits(c->fields.potential_V, s.points[1].dc.fields.potential_V));

    const std::string table = text_of(dir / "ac.csv");
    CHECK(table.starts_with("point,frequency Hz,anode V,cathode V,G anode,anode S,G anode,cathode S,"
                            "G cathode,anode S,G cathode,cathode S,C anode,anode F,C anode,cathode F,"
                            "C cathode,anode F,C cathode,cathode F\n0,1000,-0.20000000000000001,0,"));
    CHECK(lines_of(table) == 1 + 2 * 3);
    // The row of point 1 at 1 GHz, rebuilt from the run: G then C, row-major, 17 digits.
    const auto g17 = [](double v) {
        char b[32];
        std::snprintf(b, sizeof b, "%.17g", v);
        return std::string(b);
    };
    std::string row = "\n1,1000000000,0,0";
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) row += "," + g17(s.conductance(1, 2, i, j));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) row += "," + g17(s.capacitance(1, 2, i, j));
    }
    CHECK(table.find(row + "\n") != std::string::npos);
}

TEST_CASE("deck run: a [trace] deck equals trace_bias; its files and table", "[deck]") {
    const fs::path dir = scratch("run_trace");
    // A solver option too, so the trace is seen to run with the plan's options as its steady ones.
    const deck::Plan plan = must(diode_device + R"(
[solver]
tol_update = 1e-12
[trace]
contact = anode
end = 0.3
step = 0.05
max_step = 0.1

[output]
sweep = out.nitsweep
checkpoint = out.nitckpt
table = tr.csv
)", dir);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    const results::Sweep& s = *report->sweep;
    REQUIRE_FALSE(s.stopped.has_value());

    solve::TraceOptions o;
    o.contact = 0;
    o.start_V = {0.0, 0.0};
    o.end_V = 0.3;
    o.step_V = 0.05;
    o.max_step_V = 0.1;
    o.steady.newton.tol_update = 1e-12;
    const auto hand = solve::trace_bias(hand_diode(), o);
    REQUIRE(hand.has_value());
    CHECK(s.run.input_identity == hand->run.input_identity);
    REQUIRE(s.points.size() == hand->points.size());
    REQUIRE(s.points.size() > 2);
    for (std::size_t k = 0; k < s.points.size(); ++k) {
        CHECK(same_bits(s.points[k].bias_V, hand->points[k].bias_V));
        CHECK(same_bits(s.points[k].terminal_current, hand->points[k].terminal_current));
    }

    CHECK(report->written == std::vector<fs::path>{dir / "out.nitsweep", dir / "tr.csv", dir / "out.nitckpt"});
    const auto loaded = io::load_sweep(dir / "out.nitsweep");
    REQUIRE(loaded.has_value());
    CHECK(loaded->run.input_identity == s.run.input_identity);
    CHECK(same_bits(loaded->points.back().terminal_current, s.points.back().terminal_current));
    const auto c = io::load_checkpoint(dir / "out.nitckpt", plan.device.mesh());
    REQUIRE(c.has_value());
    CHECK(same_bits(c->bias_V, s.points.back().bias_V));

    const std::string table = text_of(dir / "tr.csv");
    CHECK(table.starts_with("point,anode V,cathode V,anode A,cathode A,substeps\n0,0,0,"));
    CHECK(lines_of(table) == 1 + static_cast<long>(s.points.size()));
}
