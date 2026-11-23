// [heavy_ion] decks (Unit 65; DECISIONS.md HD8, HD9; gates HC9, HC10): a deck run equals the C++
// call bit for bit with the same run identity; the table's collected-charge columns equal
// analysis::collected_charge; the generation field and the deposit round-trip through the transient
// file and plot; refusals and missing keys; the example deck.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/physics/heavy_ion.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_device;
using deck_test::hand_diode;
using deck_test::must;
using deck_test::plan_of;
using deck_test::scratch;
namespace fs = std::filesystem;

namespace {

bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
bool same(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(), [](double x, double y) { return same(x, y); });
}

std::string text_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

const std::string strike = R"(
[heavy_ion strike]
start = 0.9um
direction = 1
length = 0.2um
let = 0.001
area = 1um2
t0 = 20ps
sigma = 2ps
)";

const std::string transient = R"(
[transient]
contact anode = constant(-2)
t_end = 60ps
)";

// The first error message of the deck `text`, or "" when it builds.
std::string refusal(const std::string& text) {
    const auto p = plan_of(text);
    return p ? std::string() : p.error().message;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream s(line);
    for (std::string item; std::getline(s, item, ',');) out.push_back(item);
    return out;
}

}  // namespace

TEST_CASE("heavy ion deck: equals solve_transient; table, file and plot (HC9, HC10)", "[deck]") {
    const fs::path dir = scratch("heavy_ion");
    const deck::Plan plan = must(diode_device + transient + strike + R"(
[output]
transient = out.nittran
table = tr.csv
)", dir);
    REQUIRE(plan.transient.heavy_ions.size() == 1);
    const physics::HeavyIonTrack& t = plan.transient.heavy_ions[0];
    CHECK(t.start_cm[0] == 0.9 * 1e-4);
    CHECK(t.length_cm == 0.2 * 1e-4);
    CHECK(t.let_C_per_cm == 0.001 * 1e-12 / 1e-4);
    CHECK(t.area_cm2 == 1e-8);
    CHECK(t.t0_s == 20e-12);
    CHECK(t.sigma_s == 2e-12);

    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->transient.has_value());
    const results::Transient& r = *report->transient;
    REQUIRE_FALSE(r.stopped.has_value());

    solve::TransientOptions o;
    o.t_end_s = 60e-12;
    o.heavy_ions = {physics::HeavyIonTrack{.start_cm = {0.9 * 1e-4, 0.0, 0.0}, .direction = {1.0, 0.0, 0.0},
                                           .length_cm = 0.2 * 1e-4, .let_C_per_cm = 0.001 * 1e-12 / 1e-4,
                                           .area_cm2 = 1e-8, .t0_s = 20e-12, .sigma_s = 2e-12}};
    const std::vector<solve::Waveform> w{solve::Waveform::constant(-2.0), solve::Waveform::constant(0.0)};
    const auto hand = solve::solve_transient(hand_diode(), w, o);
    if (!hand) FAIL(hand.error().message);
    CHECK(r.run.input_identity == hand->run.input_identity);
    REQUIRE(r.points.size() == hand->points.size());
    bool equal = true;
    for (std::size_t k = 0; k < r.points.size(); ++k) {
        equal = equal && same(r.points[k].time_s, hand->points[k].time_s) &&
                same(r.points[k].terminal_current, hand->points[k].terminal_current) &&
                r.points[k].heavy_ion_charge.has_value() &&
                same(*r.points[k].heavy_ion_charge, *hand->points[k].heavy_ion_charge);
    }
    CHECK(equal);
    // The deposit: 0 at the start, the whole pulse by t_end (t0 + 20 sigma).
    CHECK(*r.points.front().heavy_ion_charge == 0.0);
    CHECK(*r.points.back().heavy_ion_charge > 0.0);

    // The table's collected-charge columns are analysis::collected_charge's, the last the deposit.
    const std::string table = text_of(dir / "tr.csv");
    std::stringstream lines(table);
    std::string header;
    std::getline(lines, header);
    CHECK(header.ends_with(",anode collected C,cathode collected C,deposited C"));
    const auto q_cathode = *analysis::collected_charge(r, 1);
    std::size_t row = 0;
    double worst = 0.0;
    for (std::string line; std::getline(lines, line); ++row) {
        const auto cells = split(line);
        REQUIRE(cells.size() >= 3);
        worst = std::max(worst, std::abs(std::stod(cells[cells.size() - 2]) - q_cathode[row]));
        CHECK(std::stod(cells.back()) == *r.points[row].heavy_ion_charge);  // 17 digits round-trip
    }
    CHECK(row == r.points.size());
    CHECK(worst == 0.0);
    std::printf("HC9: collected by the cathode %.6e of %.6e deposited (C/cm^2)\n", q_cathode.back(),
                *r.points.back().heavy_ion_charge);

    // The file: the generation field and the deposit round-trip; the field plots.
    const auto loaded = io::load_transient(dir / "out.nittran");
    if (!loaded) FAIL(loaded.error().message);
    REQUIRE(loaded->snapshots.size() == r.snapshots.size());
    CHECK(same(loaded->snapshots.back().rates.heavy_ion, r.snapshots.back().rates.heavy_ion));
    CHECK(r.snapshots.back().rates.heavy_ion.size() == plan.device.mesh().node_count());
    CHECK(r.snapshots.back().rates.net_n.size() == plan.device.mesh().node_count());  // rates on (default)
    CHECK(same(*loaded->points.back().heavy_ion_charge, *r.points.back().heavy_ion_charge));
    // At t_end (t0 + 20 sigma) the rate is 0 to double; a point at the peak has it.
    results::BiasPoint peak;
    peak.fields = r.snapshots.back().fields;
    peak.bands = r.snapshots.back().bands;
    peak.rates.heavy_ion.assign(plan.device.mesh().node_count(), 0.0);
    {
        solve::TransientOptions o2 = o;
        o2.output_times_s = {20e-12};
        const auto at_peak = solve::solve_transient(hand_diode(), w, o2);
        REQUIRE(at_peak.has_value());
        REQUIRE(at_peak->snapshots.size() == 3);
        peak.rates = at_peak->snapshots[1].rates;
        const auto g = render::field_values(plan.device, peak, render::Field::heavy_ion_generation);
        if (!g) FAIL(g.error().message);
        const double expected = physics::heavy_ion_pairs(*physics::check_heavy_ion_track(o.heavy_ions[0], 1), 1,
                                                         {1e-4, 0.0, 0.0}) *
                                physics::heavy_ion_time_profile(o.heavy_ions[0], 20e-12);
        CHECK(*std::ranges::max_element(*g) == expected);
        CHECK(*std::ranges::min_element(*g) == 0.0);
    }
    CHECK(render::find_field("Gion") == render::Field::heavy_ion_generation);
    results::BiasPoint none;
    none.fields = r.snapshots.back().fields;
    none.bands = r.snapshots.back().bands;
    CHECK_FALSE(render::field_values(plan.device, none, render::Field::heavy_ion_generation).has_value());
}

TEST_CASE("heavy ion deck: refusals and missing keys (HC9, HD6, HD9)", "[deck]") {
    CHECK(refusal(diode_device + transient + strike).empty());
    CHECK(refusal(diode_device + "\n[sweep]\ncontact anode = 0\n" + strike).find("needs a [transient]") !=
          std::string::npos);
    for (const char* key : {"start", "direction", "length", "let", "area", "t0", "sigma"}) {
        INFO(key);
        std::string s = strike;
        const auto at = s.find(std::string("\n") + key + " ");
        REQUIRE(at != std::string::npos);
        s.erase(at + 1, s.find('\n', at + 1) - at);
        CHECK_FALSE(refusal(diode_device + transient + s).empty());
    }
    CHECK(refusal(diode_device + transient + strike + "let_mev = 1\n").find("not both") != std::string::npos);
    CHECK(refusal(diode_device + transient + strike + "depth = 1um\n").find("depth") != std::string::npos);
    CHECK(refusal(diode_device + transient + strike + "colour = red\n").find("unknown key") != std::string::npos);
    std::string early = strike;
    early.replace(early.find("t0 = 20ps"), 9, "t0 = 7ps");
    CHECK(refusal(diode_device + transient + early).find("4 sigma") != std::string::npos);
    std::string mev = strike;
    mev.replace(mev.find("let = 0.001"), 11, "let_mev = 0.1");
    const auto p = plan_of(diode_device + transient + mev);
    REQUIRE(p.has_value());
    CHECK(p->transient.heavy_ions[0].let_C_per_cm == physics::let_from_mev_cm2_mg(0.1));
    // Two tracks.
    std::string second = strike;
    second.replace(second.find("strike"), 6, "second");
    const auto two = plan_of(diode_device + transient + strike + second);
    REQUIRE(two.has_value());
    CHECK(two->transient.heavy_ions.size() == 2);
}

TEST_CASE("heavy ion deck: the example deck runs (HC10)", "[deck]") {
    // A copy in a scratch directory: the run writes its outputs beside the deck.
    const fs::path deck = scratch("heavy_ion_example") / "seu_diode_2d.deck";
    fs::copy_file(fs::path(NITCAD_EXAMPLES_DIR) / "seu_diode_2d.deck", deck);
    const auto d = deck::read(deck);
    if (!d) FAIL(d.error().message);
    const auto plan = deck::build(*d);
    if (!plan) FAIL(plan.error().message);
    const auto start = std::chrono::steady_clock::now();
    const auto report = deck::run(*plan);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (!report) FAIL(report.error().message);
    REQUIRE(report->transient.has_value());
    const results::Transient& r = *report->transient;
    REQUIRE_FALSE(r.stopped.has_value());
    const auto q = *analysis::collected_charge(r, 1);
    std::printf("HC10: example deck %.1f s, %zu points; collected %.4e of %.4e deposited (C/cm); %zu warnings\n",
                seconds, r.points.size(), q.back(), *r.points.back().heavy_ion_charge, r.run.warnings.size());
    CHECK(r.run.warnings.empty());
    CHECK(seconds < 120.0);
}
