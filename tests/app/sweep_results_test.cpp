// Unit 26n slice 3 (DT7): a dimension sweep's results: the runs' folders loaded, the overlay, Unit 24's figures
// against the swept parameter and their CSV.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "NiTCAD/app/sweep_results.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

std::string read(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

// The example PN diode swept over Lp by `nitcad sweep`, in a scratch folder: the deck's text and path.
struct Swept {
    std::string text;
    fs::path path;
};

Swept diode_swept(const std::string& name, const std::string& values) {
    const fs::path dir = app_test::scratch(name);
    Swept s{read(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "pn_diode_1d.deck"), {}};
    s.path = app_test::write(dir, "pn.deck", s.text);
    const auto r = app_test::nitcad("sweep " + app_test::quoted(s.path) + " Lp=" + values);
    INFO(r.output);
    REQUIRE(r.status == 0);
    return s;
}

std::string mosfet_with_wf() {
    std::string t = app_test::mosfet2d;
    t.replace(t.find("work_function = 4.1"), 19, "work_function = $wf");
    return "[parameters]\nwf = 4.1  # Gate work function\n\n" + t;
}

}  // namespace

TEST_CASE("sweep results: each value's folder, loaded where it has run", "[app][sweep_results]") {
    const Swept s = diode_swept("sweep_results_load", "1um,1.5um");
    const auto r = app::load_sweep(s.text, s.path, "Lp", {"1um", "1.5um", "2um"});
    REQUIRE(r);
    CHECK(r->name == "Lp");
    CHECK(r->length);
    CHECK(r->directory == s.path.parent_path() / "pn_sweep_Lp");
    REQUIRE(r->runs.size() == 3);
    CHECK(r->runs[0].value == "1um");
    CHECK(r->runs[0].numeric == 1e-4);
    CHECK(r->runs[0].folder == r->directory / "Lp=1um");
    REQUIRE(r->runs[0].workspace);
    REQUIRE(r->runs[1].workspace);
    CHECK(r->runs[0].workspace->has_results());
    CHECK_FALSE(r->runs[2].workspace);
    CHECK(r->runs[2].error == "not run");

    CHECK_FALSE(app::load_sweep(s.text, s.path, "Wfin", {"1nm"}));
    CHECK_FALSE(app::load_sweep(s.text, s.path, "Lp", {"1"}));
    CHECK_FALSE(app::load_sweep("[parameters\n", s.path, "Lp", {"1um"}));
}

TEST_CASE("sweep results: labels and the shown value", "[app][sweep_results]") {
    app::SweepResults s;
    s.name = "Lg";
    s.length = true;
    app::SweepRun a{.value = "20nm", .numeric = 20e-7, .folder = {}, .workspace = std::nullopt, .error = {}};
    app::SweepRun b{.value = "1.5um", .numeric = 1.5e-4, .folder = {}, .workspace = std::nullopt, .error = {}};
    CHECK(app::run_label(s, a) == "Lg = 20 nm");
    CHECK(app::run_label(s, b) == "Lg = 1.5 \xC2\xB5m");
    CHECK(app::shown_value(s, b) == 1500.0);
    CHECK(app::value_title(s) == "Lg (nm)");
    app::SweepResults n;
    n.name = "Nd";
    app::SweepRun c{.value = "1e18", .numeric = 1e18, .folder = {}, .workspace = std::nullopt, .error = {}};
    CHECK(app::run_label(n, c) == "Nd = 1e18");
    CHECK(app::shown_value(n, c) == 1e18);
    CHECK(app::value_title(n) == "Nd");
}

TEST_CASE("sweep results: the overlay holds one curve per loaded run, viridis in sweep order", "[app][sweep_results]") {
    const Swept s = diode_swept("sweep_results_overlay", "1um,1.5um,2um");
    const auto r = app::load_sweep(s.text, s.path, "Lp", {"1um", "1.5um", "2um", "3um"});
    REQUIRE(r);
    const auto plot = app::sweep_overlay(*r);
    REQUIRE(plot);
    REQUIRE(plot->series.size() == 3);
    CHECK(plot->series[0].name == "Lp = 1 \xC2\xB5m");
    CHECK(plot->series[2].name == "Lp = 2 \xC2\xB5m");
    CHECK(plot->x.label == "V(anode) (V)");
    const auto lum = [](render::Rgb c) { return 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b; };
    CHECK(lum(plot->series[0].color) < lum(plot->series[1].color));
    CHECK(lum(plot->series[1].color) < lum(plot->series[2].color));
    // The curves are the runs' own: a longer p side carries less current at the same bias.
    CHECK(std::abs(plot->series[0].y.back()) > std::abs(plot->series[2].y.back()));
    CHECK(plot->series[0].x == plot->series[2].x);

    app::SweepResults none;
    none.name = r->name;
    none.length = true;
    for (const auto& run : r->runs) none.runs.push_back({.value = run.value, .numeric = run.numeric, .folder = run.folder, .workspace = std::nullopt, .error = "not run"});
    CHECK_FALSE(app::sweep_overlay(none));
}

TEST_CASE("sweep results: the figures against the parameter, and their CSV", "[app][sweep_results]") {
    const Swept s = diode_swept("sweep_results_figures", "1um,1.5um");
    const auto r = app::load_sweep(s.text, s.path, "Lp", {"1um", "1.5um", "2um"});
    REQUIRE(r);
    const app::SweepFigures f = app::sweep_figures(*r);
    REQUIRE(f.names.size() == 1);
    CHECK(f.names[0] == "Current at the last bias");
    CHECK(f.units[0] == "A/cm\xC2\xB2");
    REQUIRE(f.rows.size() == 2);  // the unrun value has no row
    CHECK(f.runs == std::vector<std::size_t>{0, 1});
    const auto own = app::extracted_figures(*r->runs[1].workspace);
    REQUIRE(own.size() == 1);
    CHECK(f.rows[1][0] == own[0].value);

    const std::string csv = app::figures_csv(*r, f);
    const std::string header = "Lp (nm),Current at the last bias (A/cm\xC2\xB2)\n";
    REQUIRE(csv.starts_with(header));
    CHECK(csv.substr(header.size()).starts_with("1000,"));
    CHECK(std::ranges::count(csv, '\n') == 3);
    CHECK(app::figures_csv_path(*r) == r->directory / "figures.csv");

    const auto plot = app::figure_plot(*r, f, 0);
    REQUIRE(plot);
    REQUIRE(plot->series.size() == 1);
    CHECK(plot->series[0].x == std::vector<double>{1000.0, 1500.0});
    CHECK(plot->series[0].markers);
    CHECK(plot->x.label == "Lp (nm)");
    CHECK_FALSE(app::figure_plot(*r, f, 1));
}

TEST_CASE("sweep results: a transfer sweep's thresholds follow the swept work function", "[app][sweep_results]") {
    const fs::path dir = app_test::scratch("sweep_results_mosfet");
    const std::string text = mosfet_with_wf();
    const fs::path path = app_test::write(dir, "m.deck", text);
    const auto cli = app_test::nitcad("sweep " + app_test::quoted(path) + " wf=4.1,4.3");
    INFO(cli.output);
    REQUIRE(cli.status == 0);
    const auto r = app::load_sweep(text, path, "wf", {"4.1", "4.3"});
    REQUIRE(r);
    CHECK_FALSE(r->length);
    const app::SweepFigures f = app::sweep_figures(*r);
    const auto vth = std::ranges::find(f.names, std::string("Threshold (max gm)"));
    REQUIRE(vth != f.names.end());
    const auto j = static_cast<std::size_t>(vth - f.names.begin());
    REQUIRE(f.rows.size() == 2);
    CHECK(f.rows[1][j] > f.rows[0][j] + 0.1);  // 0.2 eV more work function, a higher threshold
    const auto plot = app::sweep_overlay(*r, true);
    REQUIRE(plot);
    CHECK(plot->series.size() == 2);
    CHECK(plot->x.label == "V(gate) (V)");
    CHECK(app::figures_csv(*r, f).starts_with("wf,"));
}

TEST_CASE("sweep results: the dialog's values, a bare number taking the row's unit", "[app][sweep_results]") {
    const std::string deck = "[parameters]\nLg = 20nm\nLsd = 0.1um\nx = $Lg + $Lsd\nNd = 1e18\n";
    CHECK(app::typed_sweep_values(deck, "Lg", "20, 30nm, 45") == std::vector<std::string>{"20nm", "30nm", "45nm"});
    CHECK(app::typed_sweep_values(deck, "Lsd", "0.1 : 0.3 : 0.1") == std::vector<std::string>{"0.1um", "0.2um", "0.3um"});
    CHECK(app::typed_sweep_values(deck, "Nd", "1e17, 1e18") == std::vector<std::string>{"1e17", "1e18"});
    CHECK_FALSE(app::typed_sweep_values(deck, "Nd", "1nm"));
    CHECK_FALSE(app::typed_sweep_values(deck, "Lg", "abc"));
    CHECK_FALSE(app::typed_sweep_values(deck, "Lg", ""));
    CHECK_FALSE(app::typed_sweep_values(deck, "x", "1nm"));  // a derived row is not swept
    CHECK_FALSE(app::typed_sweep_values(deck, "Wfin", "1nm"));
    const auto e = app::typed_sweep_values(deck, "Lg", "20, abc");
    REQUIRE_FALSE(e);
    CHECK(e.error().find("abc") != std::string::npos);
}
