// The [transient], [small_signal] and [trace] sections (Unit 43; DECISIONS.md AN1-AN6, gate AG3):
// every key reaches its option, absent keys give the C++ defaults, logspace and the waveforms are
// read exactly, and every mistake names its line.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_device;
using deck_test::must;
using deck_test::plan_of;

namespace {

std::string failure(const std::string& text) {
    const auto p = plan_of(text);
    REQUIRE_FALSE(p.has_value());
    return p.error().message;
}

void expect(const std::string& message, const std::string& what) {
    INFO(message);
    CHECK(message.find(what) != std::string::npos);
}

// The diode with an isothermal sink, for the thermal keys.
const std::string thermal_diode = diode_device + R"(
[thermal_contact sink]
patch = x_max
kind = isothermal
temperature = 300
)";

}  // namespace

TEST_CASE("deck analyses: [transient] keys reach TransientOptions", "[deck]") {
    const auto defaults = must(diode_device + "[transient]\nt_end = 1e-9\n");
    REQUIRE(defaults.analysis == deck::Analysis::transient);
    const solve::TransientOptions cpp;
    const solve::TransientOptions& d = defaults.transient;
    CHECK(d.t_end_s == 1e-9);
    CHECK(d.integrator == cpp.integrator);
    CHECK(d.dt_initial_s == cpp.dt_initial_s);
    CHECK(d.dt_min_s == cpp.dt_min_s);
    CHECK(d.dt_max_s == cpp.dt_max_s);
    CHECK(d.adaptive == cpp.adaptive);
    CHECK(d.rtol == cpp.rtol);
    CHECK(d.density_ref_cm3 == cpp.density_ref_cm3);
    CHECK(d.output_times_s.empty());
    CHECK(d.fields_every_step == cpp.fields_every_step);
    REQUIRE(defaults.waveforms.size() == 2);
    for (const auto& w : defaults.waveforms) CHECK(w.parameters() == solve::Waveform::constant(0.0).parameters());
    CHECK(defaults.options.thermal_bias_K.empty());

    const auto p = must(thermal_diode + R"(
[transient]
contact anode = pulse(0, 0.6, 1e-10, 2e-10, 1e-11, 2e-11)
contact cathode = sine(0, 0.01, 1e9, 1e-10, 0.5)
t_end = 1e-9
dt_initial = 1e-12
dt_min = 1e-16
dt_max = 1e-10
integrator = backward_euler
adaptive = off
rtol = 1e-4
density_ref = 1e12
output_times = 2e-10, 5e-10
fields_every_step = on
thermal sink = 310
)");
    const solve::TransientOptions& o = p.transient;
    CHECK(o.t_end_s == 1e-9);
    CHECK(o.dt_initial_s == 1e-12);
    CHECK(o.dt_min_s == 1e-16);
    CHECK(o.dt_max_s == 1e-10);
    CHECK(o.integrator == solve::Integrator::backward_euler);
    CHECK_FALSE(o.adaptive);
    CHECK(o.rtol == 1e-4);
    CHECK(o.density_ref_cm3 == 1e12);
    CHECK(o.output_times_s == std::vector<double>{2e-10, 5e-10});
    CHECK(o.fields_every_step);
    CHECK(p.waveforms[0].parameters() ==
          solve::Waveform::pulse(0, 0.6, 1e-10, 2e-10, 1e-11, 2e-11)->parameters());
    CHECK(p.waveforms[1].kind() == solve::Waveform::Kind::sine);
    CHECK(p.waveforms[1].parameters() == solve::Waveform::sine(0, 0.01, 1e9, 1e-10, 0.5)->parameters());
    CHECK(p.options.thermal_bias_K == std::vector<std::vector<double>>{{310.0}});
}

TEST_CASE("deck analyses: every waveform form", "[deck]") {
    const auto wave = [](const std::string& w) {
        const auto p = must(diode_device + "[transient]\nt_end = 1e-9\ncontact anode = " + w + "\n");
        return p.waveforms[0].parameters();
    };
    CHECK(wave("constant(0.3)") == solve::Waveform::constant(0.3).parameters());
    CHECK(wave("step(0, 0.5, 1e-10)") == solve::Waveform::step(0, 0.5, 1e-10)->parameters());
    CHECK(wave("ramp(0, 0.5, 1e-10, 3e-10)") == solve::Waveform::ramp(0, 0.5, 1e-10, 3e-10)->parameters());
    CHECK(wave("pulse(0, 1, 1e-10, 2e-10)") == solve::Waveform::pulse(0, 1, 1e-10, 2e-10)->parameters());
    CHECK(wave("pwl(0 0, 1e-10 0.2, 1e-10 0.4, 5e-10 0.1)") ==
          solve::Waveform::piecewise_linear({{0, 0}, {1e-10, 0.2}, {1e-10, 0.4}, {5e-10, 0.1}})->parameters());
    CHECK(wave("sine(0.1, 0.02, 1e9)") == solve::Waveform::sine(0.1, 0.02, 1e9)->parameters());
    CHECK(wave(" step( 0 , 0.5 , 1e-10 ) ") == solve::Waveform::step(0, 0.5, 1e-10)->parameters());
}

TEST_CASE("deck analyses: [small_signal] keys and logspace", "[deck]") {
    const auto p = must(thermal_diode + R"(
[small_signal]
contact anode = 0 : 0.2 : 0.1
thermal sink = 300, 305, 310
frequencies = logspace(1e3, 1e9, 7)
field_frequencies = 1e3, 1e9
)");
    REQUIRE(p.analysis == deck::Analysis::small_signal);
    CHECK(p.points == std::vector<std::vector<double>>{{0.0, 0.0}, {0.1, 0.0}, {0.2, 0.0}});
    CHECK(p.options.thermal_bias_K == std::vector<std::vector<double>>{{300.0}, {305.0}, {310.0}});
    const auto& f = p.small_signal.frequencies_Hz;
    REQUIRE(f.size() == 7);
    CHECK(f.front() == 1e3);  // exact, not 10^3 rounded
    CHECK(f.back() == 1e9);
    for (std::size_t k = 1; k + 1 < f.size(); ++k) {
        CHECK(std::abs(f[k] / std::pow(10.0, 3.0 + static_cast<double>(k)) - 1.0) < 1e-14);
    }
    CHECK(p.small_signal.field_frequencies_Hz == std::vector<double>{1e3, 1e9});
    CHECK(p.restart.empty());

    // Without field frequencies, none; one point at 0 V when no contact is named.
    const auto q = must(diode_device + "[small_signal]\nfrequencies = 1e6\n");
    CHECK(q.points == std::vector<std::vector<double>>{{0.0, 0.0}});
    CHECK(q.small_signal.field_frequencies_Hz.empty());
    CHECK(q.small_signal.frequencies_Hz == std::vector<double>{1e6});

    // The bounds are exact even where 10^log10(x) does not round-trip.
    for (const auto& [a, b] : {std::pair{0.3, 7e5}, std::pair{1.7e-3, 0.9}, std::pair{3e3, 1.1e7},
                               std::pair{0.7, 0.07}}) {
        const auto r = must(diode_device + "[small_signal]\nfrequencies = logspace(" + std::to_string(a) +
                            ", " + std::to_string(b) + ", 5)\n");
        CAPTURE(a, b);
        CHECK(r.small_signal.frequencies_Hz.front() == std::stod(std::to_string(a)));
        CHECK(r.small_signal.frequencies_Hz.back() == std::stod(std::to_string(b)));
    }

    // logspace also in a [sweep] series.
    const auto s = must(diode_device + "[sweep]\ncontact anode = logspace(0.01, 0.1, 2)\n");
    CHECK(s.points == std::vector<std::vector<double>>{{0.01, 0.0}, {0.1, 0.0}});
}

TEST_CASE("deck analyses: [trace] keys reach TraceOptions", "[deck]") {
    const auto defaults = must(diode_device + "[trace]\ncontact = anode\nend = 0.5\n");
    REQUIRE(defaults.analysis == deck::Analysis::trace);
    const solve::TraceOptions cpp;
    const solve::TraceOptions& d = defaults.trace;
    CHECK(d.contact == 0);
    CHECK(d.end_V == 0.5);
    CHECK(d.start_V == std::vector<double>{0.0, 0.0});
    CHECK(d.step_V == cpp.step_V);
    CHECK(d.max_step_V == cpp.max_step_V);
    CHECK(d.min_step_V == cpp.min_step_V);
    CHECK(d.current_limit == cpp.current_limit);
    CHECK(d.max_points == cpp.max_points);
    CHECK(d.density_floor_cm3 == cpp.density_floor_cm3);

    const auto p = must(thermal_diode + R"(
[trace]
contact = cathode
end = -5
start anode = 0.1
start cathode = -0.2
step = 0.01
max_step = 0.5
min_step = 1e-6
current_limit = 10
max_points = 300
density_floor = 1e8
thermal sink = 320
)");
    const solve::TraceOptions& o = p.trace;
    CHECK(o.contact == 1);
    CHECK(o.end_V == -5.0);
    CHECK(o.start_V == std::vector<double>{0.1, -0.2});
    CHECK(o.step_V == 0.01);
    CHECK(o.max_step_V == 0.5);
    CHECK(o.min_step_V == 1e-6);
    CHECK(o.current_limit == 10.0);
    CHECK(o.max_points == 300);
    CHECK(o.density_floor_cm3 == 1e8);
    CHECK(p.options.thermal_bias_K == std::vector<std::vector<double>>{{320.0}});
}

TEST_CASE("deck analyses: outputs each analysis writes", "[deck]") {
    const auto outputs = [](const std::string& analysis, const std::string& output) {
        return diode_device + analysis + "[output]\n" + output + "\n";
    };
    const std::string transient = "[transient]\nt_end = 1e-9\n";
    const std::string small_signal = "[small_signal]\nfrequencies = 1e6\n";
    const std::string trace = "[trace]\ncontact = anode\nend = 0.5\n";

    const auto t = must(outputs(transient, "transient = t.nittran\ntable = t.csv\ndevice = d.nitdev\nmesh = m.nitmesh"));
    CHECK(t.outputs.transient.filename() == "t.nittran");
    CHECK(t.outputs.table.filename() == "t.csv");
    const auto s = must(outputs(small_signal, "small_signal = s.nitac\ncheckpoint = c.nitckpt\ntable = s.csv"));
    CHECK(s.outputs.small_signal.filename() == "s.nitac");
    CHECK(s.outputs.checkpoint.filename() == "c.nitckpt");
    const auto r = must(outputs(trace, "sweep = r.nitsweep\ncheckpoint = c.nitckpt\ntable = r.csv"));
    CHECK(r.outputs.sweep.filename() == "r.nitsweep");

    expect(failure(outputs(transient, "checkpoint = c.nitckpt")), "a transient run writes no checkpoint");
    expect(failure(outputs(transient, "sweep = s.nitsweep")), "a transient run writes no sweep");
    expect(failure(outputs(transient, "small_signal = s.nitac")), "a transient run writes no small_signal");
    expect(failure(outputs(small_signal, "sweep = s.nitsweep")), "a small-signal run writes no sweep");
    expect(failure(outputs(small_signal, "transient = t.nittran")), "a small-signal run writes no transient");
    expect(failure(outputs(trace, "transient = t.nittran")), "a trace run writes no transient");
    expect(failure(outputs(trace, "small_signal = s.nitac")), "a trace run writes no small_signal");
    expect(failure(outputs("[sweep]\ncontact anode = 0\n", "transient = t.nittran")),
           "a sweep run writes no transient");
    expect(failure(outputs("[equilibrium]\n", "small_signal = s.nitac")),
           "an equilibrium run writes no small_signal");
}

TEST_CASE("deck analyses: every mistake names its line", "[deck]") {
    // "t.deck:<line>: " for the k-th line after diode_device (k = 0: the section header).
    const std::string d = diode_device;
    const auto lines = static_cast<int>(std::count(d.begin(), d.end(), '\n'));
    const auto at = [lines](int k) { return "t.deck:" + std::to_string(lines + 1 + k) + ": "; };
    expect(failure(d + "[transient]\ncontact anode = step(0, 1, 1e-10)\n"), at(0) + "[transient] needs t_end");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact anode = wobble(1)\n"),
           at(2) + "unknown waveform 'wobble'");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact anode = 0.5\n"), at(2) + "'0.5' is not a waveform");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact anode = step(0, 1)\n"),
           at(2) + "step takes (v0, v1, t)");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact anode = pulse(0, 1, 0, 1e-9, 1e-11)\n"),
           at(2) + "pulse takes both rise and fall, or neither");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact anode = ramp(0, 1, 2e-10, 1e-10)\n"), at(2) + "");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact anode = pwl(0 0, 1e-10)\n"),
           at(2) + "pwl takes corners 't v'");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact anode = sine(0, 1, -5)\n"), at(2) + "");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact anode = constant(1)\ncontact anode = constant(2)\n"),
           at(3) + "a second waveform for 'anode'");
    expect(failure(d + "[transient]\nt_end = 1e-9\ncontact gate = constant(1)\n"), at(2) + "no contact named 'gate'");
    expect(failure(d + "[transient]\nt_end = 1e-9\nintegrator = rk4\n"), at(2) + "integrator is bdf2 or backward_euler");
    expect(failure(d + "[transient]\nt_end = 1e-9\nadaptive = maybe\n"), at(2) + "'maybe' is not on or off");
    expect(failure(d + "[transient]\nt_end = 1e-9\nthermal sink = 300\n"), at(2) + "no thermal contact named 'sink'");
    expect(failure(d + "[transient]\nt_end = 1e-9\nrestart = c.nitckpt\n"),
           at(2) + "restart is only for [sweep] (DK9), not [transient]");
    expect(failure(d + "[transient]\nt_end = 1e-9\nsteps = 10\n"), at(2) + "unknown key 'steps' in [transient]");

    expect(failure(d + "[small_signal]\ncontact anode = 0\n"), at(0) + "[small_signal] needs frequencies");
    expect(failure(d + "[small_signal]\nfrequencies = logspace(0, 1e6, 7)\n"),
           at(1) + "logspace needs positive bounds");
    expect(failure(d + "[small_signal]\nfrequencies = logspace(1, 1e6, 1)\n"),
           at(1) + "logspace needs a count from 2");
    expect(failure(d + "[small_signal]\nfrequencies = logspace(1, 1e6)\n"),
           at(1) + "logspace takes (first, last, count)");
    expect(failure(d + "[small_signal]\nfrequencies = 1e6\nrestart = c.nitckpt\n"),
           at(2) + "restart is only for [sweep] (DK9), not [small_signal]");
    expect(failure(d + "[small_signal]\nfrequencies = 1e6\nvoltage anode = 1\n"),
           at(2) + "unknown key 'voltage anode' in [small_signal]");

    expect(failure(d + "[trace]\nend = 1\n"), at(0) + "[trace] needs contact");
    expect(failure(d + "[trace]\ncontact = anode\n"), at(0) + "[trace] needs end");
    expect(failure(d + "[trace]\ncontact = gate\nend = 1\n"), at(1) + "no contact named 'gate'");
    expect(failure(d + "[trace]\ncontact = anode\nend = 1\nstart anode = 0\nstart anode = 1\n"),
           at(4) + "a second start for 'anode'");
    expect(failure(d + "[trace]\ncontact = anode\nend = 1\nmax_points = 0\n"), at(3) + "");
    expect(failure(d + "[trace]\ncontact = anode\nend = 1\nrestart = c.nitckpt\n"),
           at(3) + "restart is only for [sweep] (DK9), not [trace]");

    expect(failure(d + "[transient]\nt_end = 1e-9\n[trace]\ncontact = anode\nend = 1\n"),
           "the deck needs exactly one analysis");
    expect(failure(d + "[sweep]\ncontact anode = 0\n[small_signal]\nfrequencies = 1e6\n"),
           "the deck needs exactly one analysis");
}
