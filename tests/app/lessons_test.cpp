// Unit 26p (GU6, GUC4): the shipped MOSFET guide lessons parse, and each one's claim holds in its own results.
// [lessons] is fast (parsing); [.lessons] runs every lesson's sweep or runs through nitcad (Release only, ctest
// app_lessons), each claim checked against the figures the guide would report.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "fixtures.hpp"
#include "NiTCAD/app/dimensions.hpp"
#include "NiTCAD/app/lesson.hpp"
#include "NiTCAD/app/sweep_results.hpp"
#include "NiTCAD/deck/dimension_sweep.hpp"

using namespace NiTCAD::app;
namespace fs = std::filesystem;

namespace {

const fs::path lessons_dir = fs::path(NITCAD_EXAMPLES_DIR) / "lessons";

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

// One figure per value: a sweep of `name` over `list` of the deck `text`, run by nitcad in a scratch folder.
std::vector<double> sweep_figure(const std::string& scratch, const std::string& text, const std::string& name,
                                 const std::string& list, const std::string& figure) {
    const fs::path dir = app_test::scratch(scratch);
    const fs::path deck = app_test::write(dir, "lesson.deck", text);
    const auto rows = dims::read_dimensions(text);
    REQUIRE(rows);
    const auto row = std::ranges::find(*rows, name, &dims::Dimension::name);
    REQUIRE(row != rows->end());
    const auto values = NiTCAD::deck::sweep_values(list, row->length);
    REQUIRE(values);
    // A sweep solved once per run of the tests, however many of its figures are checked.
    static std::map<std::string, SweepFigures> solved;
    const std::string key = text + "\n" + name + "=" + list;
    if (!solved.contains(key)) {
        std::string joined;
        for (const auto& v : *values) joined += (joined.empty() ? "" : ",") + v;
        const auto cli = app_test::nitcad("sweep " + app_test::quoted(deck) + " " + name + "=" + joined);
        INFO(cli.output);
        REQUIRE(cli.status == 0);
        const auto s = load_sweep(text, deck, name, *values);
        REQUIRE(s);
        solved[key] = sweep_figures(*s);
    }
    const SweepFigures& f = solved[key];
    const auto col = std::ranges::find(f.names, figure);
    REQUIRE(col != f.names.end());
    REQUIRE(f.runs.size() == values->size());  // every run converged and loaded
    std::vector<double> out;
    for (const auto& r : f.rows) out.push_back(r[static_cast<std::size_t>(col - f.names.begin())]);
    for (std::size_t k = 0; k < out.size(); ++k) WARN(scratch << " " << figure << " [" << (*values)[k] << "] = " << out[k]);
    return out;
}

// The lesson's first sweep step: its parameter, values and the figure asked for.
struct Swept {
    std::string text, name, list;
};
Swept lesson_sweep(const std::string& file) {
    const std::string text = read_file(lessons_dir / file);
    const auto l = parse_lesson(text);
    REQUIRE(l);
    for (const LessonStep& s : l->steps) {
        if (s.action.kind == LessonAction::Kind::sweep) return {text, s.action.name, s.action.value};
    }
    FAIL_CHECK("no sweep step in " << file);
    return {};
}

std::vector<double> lesson_figure(const std::string& file, const std::string& figure) {
    const Swept s = lesson_sweep(file);
    return sweep_figure(file, s.text, s.name, s.list, figure);
}

// A run lesson's figure, with and without its edit: the deck run as a one-value sweep of its drain bias.
std::pair<double, double> edited_figure(const std::string& file, const std::string& section, const std::string& key,
                                        const std::string& value, const std::string& figure) {
    const std::string text = read_file(lessons_dir / file);
    const auto edited = set_entry(text, section, key, value);
    REQUIRE(edited);
    const auto vd = [&](const std::string& t) {
        const auto start = t.find("\nVd = ");
        REQUIRE(start != std::string::npos);
        const auto from = start + 6;
        return t.substr(from, t.find(' ', from) - from);
    };
    const double before = sweep_figure(file + "_before", text, "Vd", vd(text), figure).at(0);
    const double after = sweep_figure(file + "_after", *edited, "Vd", vd(text), figure).at(0);
    return {before, after};
}

bool rising(const std::vector<double>& v) {
    for (std::size_t k = 1; k < v.size(); ++k) {
        if (!(v[k] > v[k - 1])) return false;
    }
    return v.size() > 1;
}

bool falling(const std::vector<double>& v) {
    for (std::size_t k = 1; k < v.size(); ++k) {
        if (!(v[k] < v[k - 1])) return false;
    }
    return v.size() > 1;
}

}  // namespace

TEST_CASE("lessons: every shipped lesson parses, thirteen in three tracks (GU6)", "[lessons]") {
    const auto ls = list_lessons(lessons_dir);
    REQUIRE(ls.size() == 13);
    int basics = 0, intermediate = 0, advanced = 0;
    for (std::size_t k = 0; k < ls.size(); ++k) {
        CAPTURE(ls[k].title);
        CHECK(ls[k].order == static_cast<int>(k) + 1);
        CHECK(ls[k].steps.size() >= 3);
        basics += ls[k].track == LessonTrack::basics;
        intermediate += ls[k].track == LessonTrack::intermediate;
        advanced += ls[k].track == LessonTrack::advanced;
        for (const LessonStep& s : ls[k].steps) {
            CHECK_FALSE(s.prose.empty());
            for (const std::string& r : s.report) {
                CHECK((r == "Threshold (max gm)" || r == "Subthreshold swing" || r == "Max transconductance" ||
                       r == "On/off ratio" || r == "Current at the last bias"));
            }
        }
    }
    CHECK(basics == 4);
    CHECK(intermediate == 4);
    CHECK(advanced == 5);
}

TEST_CASE("lessons 1: the transfer curve gives a threshold, a swing and an on/off ratio", "[.lessons]") {
    const std::string text = read_file(lessons_dir / "01_gate_and_channel.deck");
    for (const char* f : {"Threshold (max gm)", "Subthreshold swing", "On/off ratio"}) {
        CHECK(std::isfinite(sweep_figure("lesson01", text, "Vd", "0.05", f).at(0)));
    }
}

TEST_CASE("lessons 2: V_th and the swing rise with Na", "[.lessons]") {
    CHECK(rising(lesson_figure("02_body_doping.deck", "Threshold (max gm)")));
    CHECK(rising(lesson_figure("02_body_doping.deck", "Subthreshold swing")));
}

TEST_CASE("lessons 3: V_th and the swing fall as the oxide thins", "[.lessons]") {
    CHECK(falling(lesson_figure("03_oxide_thickness.deck", "Threshold (max gm)")));
    CHECK(falling(lesson_figure("03_oxide_thickness.deck", "Subthreshold swing")));
}

TEST_CASE("lessons 4: the saturation current rises with Vg", "[.lessons]") {
    CHECK(rising(lesson_figure("04_drain_bias.deck", "Current at the last bias")));
}

TEST_CASE("lessons 5: the swing stays above 60 mV/dec and rises with Na", "[.lessons]") {
    const auto ss = lesson_figure("05_subthreshold_swing.deck", "Subthreshold swing");
    for (const double v : ss) CHECK(v >= 60.0);
    CHECK(rising(ss));
}

TEST_CASE("lessons 6: V_th follows the work function, about 0.3 V a step; the swing does not move", "[.lessons]") {
    const auto vth = lesson_figure("06_work_function.deck", "Threshold (max gm)");
    REQUIRE(vth.size() == 3);
    for (std::size_t k = 1; k < 3; ++k) {
        CHECK(vth[k] - vth[k - 1] > 0.25);
        CHECK(vth[k] - vth[k - 1] < 0.35);
    }
    const auto ss = lesson_figure("06_work_function.deck", "Subthreshold swing");
    CHECK(std::abs(ss.back() - ss.front()) < 3.0);
}

TEST_CASE("lessons 7: the peak transconductance rises with the source/drain doping", "[.lessons]") {
    CHECK(rising(lesson_figure("07_source_drain_doping.deck", "Max transconductance")));
}

TEST_CASE("lessons 8: the surface mobility lowers g_m and hardly moves V_th", "[.lessons]") {
    const auto gm = edited_figure("08_surface_mobility.deck", "models", "surface_mobility", "lombardi", "Max transconductance");
    CHECK(gm.second < gm.first);
    const auto vth = edited_figure("08_surface_mobility.deck", "models", "surface_mobility", "lombardi", "Threshold (max gm)");
    CHECK(std::abs(vth.second - vth.first) < 0.1);
}

TEST_CASE("lessons 9: V_th rolls off as Lg shrinks", "[.lessons]") {
    CHECK(falling(lesson_figure("09_short_channel.deck", "Threshold (max gm)")));
}

TEST_CASE("lessons 10: V_th falls as the drain bias rises (DIBL)", "[.lessons]") {
    CHECK(falling(lesson_figure("10_dibl.deck", "Threshold (max gm)")));
}

TEST_CASE("lessons 11: the short gate carries more current", "[.lessons]") {
    const auto id = lesson_figure("11_output_conductance.deck", "Current at the last bias");
    REQUIRE(id.size() == 2);
    CHECK(id[1] > id[0]);
}

TEST_CASE("lessons 12: velocity saturation lowers g_m", "[.lessons]") {
    const auto gm = edited_figure("12_velocity_saturation.deck", "models", "field_mobility", "off", "Max transconductance");
    CHECK(gm.first < gm.second);
}

TEST_CASE("lessons 13: the two finer meshes agree more closely than the coarse one", "[.lessons]") {
    const auto vth = lesson_figure("13_mesh_and_models.deck", "Threshold (max gm)");
    REQUIRE(vth.size() == 3);
    CHECK(std::abs(vth[2] - vth[1]) < std::abs(vth[1] - vth[0]));
    CHECK(std::abs(vth[2] - vth[1]) < 0.01);
}
