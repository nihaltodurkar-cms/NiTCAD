// Unit 26n slice 3 (DT6, DC5): parameter overrides (`--set`), dimension sweeps and their folders.

#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/dimension_sweep.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/io/files.hpp"
#include "fixtures.hpp"

namespace fs = std::filesystem;
using namespace NiTCAD;
using deck_test::scratch;

namespace {

bool has(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

const std::string params_deck = "[parameters]\n"
                                "# group: Main\n"
                                "Lg = 20nm  # Gate length\n"
                                "Nd = 1e18\n"
                                "x1 = $Lg + 5nm\n"
                                "\n"
                                "[mesh]\n"
                                "type = tensor\n"
                                "x = linspace(0, $x1, 11)\n";

deck::Deck parsed(const std::string& text) {
    auto d = deck::parse(text, "t.deck", "C:/somewhere");
    REQUIRE(d);
    return *d;
}

std::string error_of(const std::string& text, std::vector<deck::Override> o) {
    const auto r = deck::with_overrides(parsed(text), o);
    REQUIRE_FALSE(r.has_value());
    return r.error().message;
}

std::string bytes_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

void write(const fs::path& p, const std::string& s) {
    std::ofstream f(p, std::ios::binary);
    f << s;
}

struct Result {
    int status = 0;
    std::string out;
};

Result nitcad(const std::string& arguments) {
    const std::string command = "\"\"" NITCAD_CLI_PATH "\" " + arguments + " 2>&1\"";
    FILE* pipe = _popen(command.c_str(), "r");
    REQUIRE(pipe != nullptr);
    std::string out;
    char buffer[512];
    while (std::fgets(buffer, sizeof buffer, pipe)) out += buffer;
    return {_pclose(pipe), out};
}

// A .nitsweep with its provenance (build, start, wall time) cleared and written again: its bytes.
std::string sweep_bytes(const fs::path& file, const fs::path& out) {
    auto s = io::load_sweep(file);
    if (!s) FAIL(file.string() << ": " << s.error().message);
    s->run.provenance = {};
    REQUIRE(io::save_sweep(out, *s).has_value());
    return bytes_of(out);
}

}  // namespace

TEST_CASE("overrides: parse_override splits name=value", "[deck][dimension_sweep]") {
    const auto o = deck::parse_override("Lg=30nm");
    REQUIRE(o);
    CHECK(o->name == "Lg");
    CHECK(o->value == "30nm");
    const auto spaced = deck::parse_override(" Lg = 30nm ");
    REQUIRE(spaced);
    CHECK(spaced->name == "Lg");
    CHECK(spaced->value == "30nm");
    CHECK_FALSE(deck::parse_override("Lg"));
    CHECK_FALSE(deck::parse_override("=30nm"));
    CHECK_FALSE(deck::parse_override("Lg="));
}

TEST_CASE("overrides: the deck as if the value were written on its line", "[deck][dimension_sweep]") {
    const auto d = deck::with_overrides(parsed(params_deck), {{.name = "Lg", .value = "30nm"}});
    REQUIRE(d);
    std::string expected = params_deck;
    expected.replace(expected.find("Lg = 20nm"), 9, "Lg = 30nm");
    CHECK(d->text == expected);
    CHECK(d->source == "t.deck");
    CHECK(d->directory == fs::path("C:/somewhere"));
    const auto p = deck::parameters(*d);
    REQUIRE(p);
    CHECK((*p)[2].value == 35e-7);  // x1 follows Lg

    // CRLF lines keep their CR; two overrides at once; a number for a number.
    std::string crlf;
    for (const char c : params_deck) {
        if (c == '\n') crlf += '\r';
        crlf += c;
    }
    const auto both = deck::with_overrides(parsed(crlf), {{.name = "Nd", .value = "2e18"}, {.name = "Lg", .value = "1um"}});
    REQUIRE(both);
    std::string want = crlf;
    want.replace(want.find("Lg = 20nm"), 9, "Lg = 1um");
    want.replace(want.find("Nd = 1e18"), 9, "Nd = 2e18");
    CHECK(both->text == want);

    // No override: the same deck text.
    const auto none = deck::with_overrides(parsed(params_deck), {});
    REQUIRE(none);
    CHECK(none->text == params_deck);
}

TEST_CASE("overrides: an unknown name or a wrong kind names the parameter", "[deck][dimension_sweep]") {
    const std::string unknown = error_of(params_deck, {{.name = "Wfin", .value = "8nm"}});
    CHECK(has(unknown, "Wfin"));
    CHECK(has(unknown, "no parameter"));
    const std::string number_for_length = error_of(params_deck, {{.name = "Lg", .value = "30"}});
    CHECK(has(number_for_length, "Lg"));
    CHECK(has(number_for_length, "length"));
    const std::string length_for_number = error_of(params_deck, {{.name = "Nd", .value = "5nm"}});
    CHECK(has(length_for_number, "Nd"));
    CHECK(has(length_for_number, "number"));
    const std::string junk = error_of(params_deck, {{.name = "Lg", .value = "abc"}});
    CHECK(has(junk, "Lg"));
    const std::string no_section = error_of("[mesh]\ntype = tensor\nx = linspace(0, 1um, 3)\n", {{.name = "Lg", .value = "1nm"}});
    CHECK(has(no_section, "Lg"));
    // The override is applied to the parameter named, not another.
    const auto d = deck::with_overrides(parsed(params_deck), {{.name = "Nd", .value = "3e18"}});
    REQUIRE(d);
    const auto p = deck::parameters(*d);
    REQUIRE(p);
    CHECK((*p)[0].value == 20e-7);
    CHECK((*p)[1].value == 3e18);
}

TEST_CASE("sweep values: a list as written, a range in its first value's unit", "[deck][dimension_sweep]") {
    const auto list = deck::sweep_values("20nm, 30nm,45nm", true);
    REQUIRE(list);
    CHECK(*list == std::vector<std::string>{"20nm", "30nm", "45nm"});
    const auto range = deck::sweep_values("20nm : 30nm : 5nm", true);
    REQUIRE(range);
    CHECK(*range == std::vector<std::string>{"20nm", "25nm", "30nm"});
    const auto um = deck::sweep_values("0.1um : 0.3um : 0.1um", true);
    REQUIRE(um);
    CHECK(*um == std::vector<std::string>{"0.1um", "0.2um", "0.3um"});
    const auto mixed = deck::sweep_values("10nm : 0.03um : 10nm", true);
    REQUIRE(mixed);
    CHECK(*mixed == std::vector<std::string>{"10nm", "20nm", "30nm"});
    const auto numbers = deck::sweep_values("1e17, 1e18", false);
    REQUIRE(numbers);
    CHECK(*numbers == std::vector<std::string>{"1e17", "1e18"});
    const auto nrange = deck::sweep_values("1 : 3 : 1", false);
    REQUIRE(nrange);
    CHECK(*nrange == std::vector<std::string>{"1", "2", "3"});

    CHECK_FALSE(deck::sweep_values("20, 30", true));         // a length needs its unit
    CHECK_FALSE(deck::sweep_values("20nm", false));          // a number has none
    CHECK_FALSE(deck::sweep_values("20nm : 10nm : 5nm", true));  // the step points away
    CHECK_FALSE(deck::sweep_values("20nm : 30nm", true));
    CHECK_FALSE(deck::sweep_values("", true));
    CHECK_FALSE(deck::sweep_values("20nm, 20nm", true));     // a value twice would share a folder
    CHECK_FALSE(deck::sweep_values("$Lsd", true));           // literals only
}

TEST_CASE("sweep folders: beside the deck, one per value", "[deck][dimension_sweep]") {
    const fs::path dir = deck::sweep_directory("C:/w/mosfet_2d.deck", "Lg");
    CHECK(dir == fs::path("C:/w/mosfet_2d_sweep_Lg"));
    CHECK(deck::sweep_point_directory(dir, "Lg", "20nm") == fs::path("C:/w/mosfet_2d_sweep_Lg/Lg=20nm"));
    CHECK(deck::sweep_point_directory(dir, "Lg", "30nm") != deck::sweep_point_directory(dir, "Lg", "20nm"));
}

TEST_CASE("sweep outputs: moved into the value's folder", "[deck][dimension_sweep]") {
    deck::Outputs p;
    p.sweep = "C:/w/out.nitsweep";
    p.device = "C:/w/sub/out.nitdev";
    p.table = "C:/elsewhere/t.csv";
    p.adapt_meshes = "C:/w/passes";
    const deck::Outputs q = deck::with_outputs_in(p, "C:/w", "C:/w/d_sweep_Lg/Lg=20nm");
    CHECK(q.sweep == fs::path("C:/w/d_sweep_Lg/Lg=20nm/out.nitsweep"));
    CHECK(q.device == fs::path("C:/w/d_sweep_Lg/Lg=20nm/sub/out.nitdev"));
    CHECK(q.table == fs::path("C:/w/d_sweep_Lg/Lg=20nm/t.csv"));
    CHECK(q.adapt_meshes == fs::path("C:/w/d_sweep_Lg/Lg=20nm/passes"));
    CHECK(q.mesh.empty());  // not written stays not written
}

TEST_CASE("DC5: a sweep's runs are the runs of copies with the value written in", "[deck][dimension_sweep]") {
    const fs::path dir = scratch("dimension_sweep_dc5");
    const fs::path deck_path = dir / "pn.deck";
    std::string text = bytes_of(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "pn_diode_1d.deck");
    REQUIRE(has(text, "Lp = 1.2um"));
    write(deck_path, text);

    const Result r = nitcad("sweep \"" + deck_path.string() + "\" Lp=1um,1.5um");
    INFO(r.out);
    REQUIRE(r.status == 0);
    CHECK(has(r.out, "Lp = 1um"));
    CHECK(has(r.out, "Lp = 1.5um"));

    for (const std::string v : {"1um", "1.5um"}) {
        const fs::path copy_dir = dir / ("copy_" + v);
        fs::create_directories(copy_dir);
        std::string copy = text;
        copy.replace(copy.find("Lp = 1.2um"), 10, "Lp = " + v);
        write(copy_dir / "pn.deck", copy);
        const Result c = nitcad("run \"" + (copy_dir / "pn.deck").string() + "\"");
        INFO(c.out);
        REQUIRE(c.status == 0);
        const fs::path point = dir / "pn_sweep_Lp" / ("Lp=" + v);
        for (const std::string f : {"pn_diode_1d.nitdev", "pn_diode_1d.csv"}) {
            INFO(f);
            REQUIRE(fs::exists(point / f));
            CHECK(bytes_of(point / f) == bytes_of(copy_dir / f));
        }
        CHECK(sweep_bytes(point / "pn_diode_1d.nitsweep", dir / "a.nitsweep") ==
              sweep_bytes(copy_dir / "pn_diode_1d.nitsweep", dir / "b.nitsweep"));
    }
    // The deck itself is not edited and its own outputs are not written.
    CHECK(bytes_of(deck_path) == text);
    CHECK_FALSE(fs::exists(dir / "pn_diode_1d.nitsweep"));
}

TEST_CASE("DC5: run --set writes the deck's outputs at the value given", "[deck][dimension_sweep]") {
    const fs::path dir = scratch("dimension_sweep_set");
    std::string text = bytes_of(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "pn_diode_1d.deck");
    write(dir / "pn.deck", text);
    const Result r = nitcad("run \"" + (dir / "pn.deck").string() + "\" --set Lp=1um --set Na=2e16");
    INFO(r.out);
    REQUIRE(r.status == 0);
    const fs::path copy_dir = dir / "copy";
    fs::create_directories(copy_dir);
    std::string copy = text;
    copy.replace(copy.find("Lp = 1.2um"), 10, "Lp = 1um");
    copy.replace(copy.find("Na = 1e16"), 9, "Na = 2e16");
    write(copy_dir / "pn.deck", copy);
    REQUIRE(nitcad("run \"" + (copy_dir / "pn.deck").string() + "\"").status == 0);
    CHECK(bytes_of(dir / "pn_diode_1d.nitdev") == bytes_of(copy_dir / "pn_diode_1d.nitdev"));
    CHECK(bytes_of(dir / "pn_diode_1d.csv") == bytes_of(copy_dir / "pn_diode_1d.csv"));
    CHECK(bytes_of(dir / "pn.deck") == text);
}

TEST_CASE("DC5: --set and sweep refuse an unknown name or a wrong kind as usage errors", "[deck][dimension_sweep]") {
    const fs::path dir = scratch("dimension_sweep_usage");
    write(dir / "pn.deck", bytes_of(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "pn_diode_1d.deck"));
    const std::string deck = "\"" + (dir / "pn.deck").string() + "\"";
    for (const std::string& args : {"run " + deck + " --set Wfin=8nm", "run " + deck + " --set Lp=5",
                                   "check " + deck + " --set Na=1nm", "sweep " + deck + " Wfin=1nm,2nm",
                                   "sweep " + deck + " Lp=1,2", "sweep " + deck + " Na=1nm"}) {
        const Result r = nitcad(args);
        INFO(args << "\n" << r.out);
        CHECK(r.status == 64);
        const bool named = has(r.out, "Wfin") || has(r.out, "Lp") || has(r.out, "Na");
        CHECK(named);
    }
    CHECK(nitcad("sweep " + deck).status == 64);
    CHECK(nitcad("run " + deck + " --set").status == 64);
    CHECK_FALSE(fs::exists(dir / "pn_sweep_Wfin"));
    CHECK_FALSE(fs::exists(dir / "pn_sweep_Lp"));
}
