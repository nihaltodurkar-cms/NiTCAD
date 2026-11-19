// Time and frequency units in the deck (Unit 61, RD1; gate RG2): a unit shifts the decimal exponent,
// so every suffixed value is the double its bare exponent form is, for every key and waveform
// argument; a unit of the wrong kind, or an unknown one, names its line.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "fixtures.hpp"
#include "values.hpp"

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

}  // namespace

TEST_CASE("deck units: each unit is its exponent, bit for bit (RG2)", "[deck][units]") {
    using deck::detail::hertz;
    using deck::detail::number;
    using deck::detail::seconds;
    const std::vector<std::pair<std::string, std::string>> times = {
        {"2.5ns", "2.5e-9"}, {"1s", "1"},       {"3ms", "3e-3"},     {"0.1us", "0.1e-6"},
        {"7ps", "7e-12"},    {"1.25fs", "1.25e-15"}, {"2.5e-3ns", "2.5e-12"}, {"2.5e+3ps", "2.5e-9"},
        {"-4ns", "-4e-9"},   {"+4ns", "4e-9"},  {"3 ns", "3e-9"},    {"0.3ns", "0.3e-9"},
    };
    for (const auto& [suffixed, bare] : times) {
        INFO(suffixed);
        const auto a = seconds(suffixed);
        const auto b = number(bare);
        REQUIRE(a.has_value());
        REQUIRE(b.has_value());
        CHECK(*a == *b);
    }
    // The unit is read as the literal would be.
    CHECK(*seconds("0.3ns") == 0.3e-9);
    const std::vector<std::pair<std::string, std::string>> frequencies = {
        {"10GHz", "1e10"}, {"1Hz", "1"}, {"2.5kHz", "2.5e3"}, {"0.7MHz", "0.7e6"}, {"3THz", "3e12"}, {"1e-3GHz", "1e6"},
    };
    for (const auto& [suffixed, bare] : frequencies) {
        INFO(suffixed);
        const auto a = hertz(suffixed);
        REQUIRE(a.has_value());
        CHECK(*a == *number(bare));
    }
    CHECK(*seconds("4e-9") == 4e-9);  // bare numbers unchanged
    CHECK(*hertz("1e6") == 1e6);
}

TEST_CASE("deck units: every time and frequency key and argument takes a unit (RG2)", "[deck][units]") {
    const std::string bare = diode_device + R"([transient]
contact anode = pulse(0, 0.5, 1e-10, 2e-10, 1e-11, 2e-11)
t_end = 1e-9
dt_initial = 1e-13
dt_min = 1e-16
dt_max = 1e-10
output_times = 2e-10, 5e-10
)";
    const std::string with_units = diode_device + R"([transient]
contact anode = pulse(0, 0.5, 100ps, 200ps, 10ps, 20ps)
t_end = 1ns
dt_initial = 0.1ps
dt_min = 0.1fs
dt_max = 100ps
output_times = 200ps, 0.5ns
)";
    const auto a = must(bare);
    const auto b = must(with_units);
    CHECK(a.transient.t_end_s == b.transient.t_end_s);
    CHECK(a.transient.dt_initial_s == b.transient.dt_initial_s);
    CHECK(a.transient.dt_min_s == b.transient.dt_min_s);
    CHECK(a.transient.dt_max_s == b.transient.dt_max_s);
    CHECK(a.transient.output_times_s == b.transient.output_times_s);
    REQUIRE(a.waveforms.size() == b.waveforms.size());
    for (std::size_t k = 0; k < a.waveforms.size(); ++k) CHECK(a.waveforms[k].parameters() == b.waveforms[k].parameters());
    CHECK(solve::make_run_record(a.device, a.transient, a.waveforms).input_identity ==
          solve::make_run_record(b.device, b.transient, b.waveforms).input_identity);

    const std::vector<std::pair<std::string, std::string>> waveforms = {
        {"step(0, 1, 2e-9)", "step(0, 1, 2ns)"},
        {"ramp(0, 1, 1e-9, 3e-9)", "ramp(0, 1, 1ns, 3ns)"},
        {"pwl(0 0, 1e-9 0.5, 2e-9 0.5)", "pwl(0s 0, 1ns 0.5, 2ns 0.5)"},
        {"sine(0, 0.01, 1e9, 1e-10, 0.5)", "sine(0, 0.01, 1GHz, 100ps, 0.5)"},
    };
    for (const auto& [plain, united] : waveforms) {
        INFO(united);
        const auto p = must(diode_device + "[transient]\nt_end = 4e-9\ncontact anode = " + plain + "\n");
        const auto q = must(diode_device + "[transient]\nt_end = 4ns\ncontact anode = " + united + "\n");
        CHECK(p.waveforms[1].parameters() == q.waveforms[1].parameters());
    }

    const auto s = must(diode_device + "[small_signal]\ncontact anode = 0.5\nfrequencies = logspace(1e3, 1e9, 7)\n"
                                       "field_frequencies = 1e6\n");
    const auto t = must(diode_device + "[small_signal]\ncontact anode = 0.5\nfrequencies = logspace(1kHz, 1GHz, 7)\n"
                                       "field_frequencies = 1MHz\n");
    CHECK(s.small_signal.frequencies_Hz == t.small_signal.frequencies_Hz);
    CHECK(s.small_signal.field_frequencies_Hz == t.small_signal.field_frequencies_Hz);
    const auto u = must(diode_device + "[small_signal]\ncontact anode = 0.5\nfrequencies = 1kHz, 10kHz, 1MHz\n");
    CHECK(u.small_signal.frequencies_Hz == std::vector<double>{1e3, 1e4, 1e6});
}

TEST_CASE("deck units: wrong and unknown units name their line (RG2)", "[deck][units]") {
    // "t.deck:<line>: " for the k-th line after diode_device (k = 0: the section header).
    const std::string d = diode_device;
    const auto lines = static_cast<int>(std::count(d.begin(), d.end(), '\n'));
    const auto at = [lines](int k) { return "t.deck:" + std::to_string(lines + 1 + k) + ": "; };
    const std::string tr = d + "[transient]\n";
    expect(failure(tr + "t_end = 1MHz\n"), at(1) + "'1MHz' is a frequency, not a time");
    expect(failure(tr + "t_end = 1nm\n"), at(1) + "'1nm': unknown unit 'nm' for a time");
    expect(failure(tr + "t_end = 1e-9\ndt_max = 2ks\n"), at(2) + "'2ks': unknown unit 'ks'");
    expect(failure(tr + "t_end = 1e-9\noutput_times = 1ns, 2GHz\n"), at(2) + "'2GHz' is a frequency, not a time");
    expect(failure(tr + "t_end = 1e-9\ncontact anode = step(0, 1ns, 1ns)\n"), at(2) + "'1ns' has a unit; a voltage is in V");
    expect(failure(tr + "t_end = 1e-9\ncontact anode = sine(0, 1, 1ns)\n"), at(2) + "'1ns' is a time, not a frequency");
    expect(failure(tr + "t_end = 1e-9\ncontact anode = sine(0, 1, 1GHz, 1GHz)\n"),
           at(2) + "'1GHz' is a frequency, not a time");
    expect(failure(tr + "t_end = 1e-9\ncontact anode = pwl(0 0, 1ns 1V)\n"), at(2) + "'1V' is not a number");
    expect(failure(tr + "t_end = 1e-9\ncontact anode = pwl(0 0, 1ns 1ns)\n"), at(2) + "'1ns' has a unit; a voltage is in V");
    expect(failure(d + "[small_signal]\ncontact anode = 0.5\nfrequencies = 1ns\n"),
           at(2) + "'1ns' is a time, not a frequency");
    expect(failure(tr + "t_end = ns\n"), at(1) + "'ns' is not a number");
    expect(failure(tr + "t_end = 1e5e2ns\n"), at(1) + "'1e5e2ns' is not a number");
}
