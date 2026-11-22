// Unit 26n slice 0 (DT1): the dimensions model — a deck's [parameters] as labelled rows, a row's value written by a
// one-line minimal edit, the mouse wheel's steps.
#include "NiTCAD/app/dimensions.hpp"

#include "NiTCAD/app/geometry.hpp"
#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/deck/syntax.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace NiTCAD::app;

namespace {

const std::string mos2d = R"([parameters]
# group: Gate
Lg = 20nm  # Gate length
Lsd = 50nm
x_d = $Lsd + $Lg
Nsd = 1e20  # S/D doping

[mesh]
type = tensor
x = linspace(0, $Lsd, 6) + linspace($Lsd, $x_d, 5) + linspace($x_d, $x_d + $Lsd, 6)
y = linspace(0, 40nm, 9) + linspace(40nm, 42nm, 3)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0 $x_d + $Lsd  40nm 42nm

[doping]
profile = acceptor constant 1e16
profile = donor constant $Nsd x=0:$Lsd
profile = donor constant $Nsd x=$x_d:$x_d + $Lsd

[contact source]
kind = ohmic
patch = x_min

[contact drain]
kind = ohmic
patch = x_max

[contact gate]
kind = electrode
patch = y_max
box = $Lsd $x_d  42nm 42nm
electrode = metal
work_function = 4.1

[equilibrium]
)";

std::vector<std::string> lines_of(const std::string& s) {
    std::vector<std::string> v;
    std::size_t at = 0;
    while (true) {
        const auto end = s.find('\n', at);
        v.push_back(s.substr(at, end == std::string::npos ? std::string::npos : end - at));
        if (end == std::string::npos) break;
        at = end + 1;
    }
    return v;
}

// The lines of b that differ from a's (same line count expected).
std::vector<std::string> changed(const std::string& a, const std::string& b) {
    const auto x = lines_of(a), y = lines_of(b);
    REQUIRE(x.size() == y.size());
    std::vector<std::string> out;
    for (std::size_t k = 0; k < x.size(); ++k) {
        if (x[k] != y[k]) out.push_back(y[k]);
    }
    return out;
}

double number_of(const std::string& expr) { return std::stod(expr); }

}  // namespace

TEST_CASE("dimensions: a deck's parameters read as labelled, grouped rows in deck order", "[app][dimensions]") {
    const auto rows = dims::read_dimensions(mos2d);
    REQUIRE(rows);
    REQUIRE(rows->size() == 4);
    const dims::Dimension& lg = (*rows)[0];
    CHECK(lg.name == "Lg");
    CHECK(lg.label == "Gate length");
    CHECK(lg.group == "Gate");
    CHECK(lg.expr == "20nm");
    CHECK(lg.length);
    CHECK(std::abs(lg.value - 2e-6) < 1e-18);
    CHECK(lg.value_text == "20 nm");
    CHECK(lg.line == 3);
    CHECK((*rows)[1].label == "Lsd");  // no comment: the name
    CHECK((*rows)[1].group == "Gate");
    const dims::Dimension& xd = (*rows)[2];
    CHECK(xd.expr == "$Lsd + $Lg");
    CHECK(std::abs(xd.value - 7e-6) < 1e-18);
    CHECK(xd.value_text == "70 nm");
    const dims::Dimension& nsd = (*rows)[3];
    CHECK(nsd.label == "S/D doping");
    CHECK_FALSE(nsd.length);
    CHECK(nsd.value == 1e20);
    CHECK(nsd.value_text == "1e+20");
}

TEST_CASE("dimensions: a deck without [parameters] has no rows; one that does not parse is an error", "[app][dimensions]") {
    const auto none = dims::read_dimensions("[mesh]\ntype = tensor\nx = linspace(0, 1um, 11)\n\n[region si]\nmaterial = silicon\n");
    REQUIRE(none);
    CHECK(none->empty());
    CHECK_FALSE(dims::read_dimensions("[parameters]\nLg = \n"));
}

TEST_CASE("dimensions: a typed value changes only its parameter's line, its comment kept", "[app][dimensions]") {
    SECTION("a length") {
        const auto out = dims::set_dimension(mos2d, "Lg", "30nm");
        REQUIRE(out);
        CHECK(changed(mos2d, *out) == std::vector<std::string>{"Lg = 30nm  # Gate length"});
    }
    SECTION("a bare number in a length row takes the row's unit") {
        const auto out = dims::set_dimension(mos2d, "Lg", "30");
        REQUIRE(out);
        CHECK(changed(mos2d, *out) == std::vector<std::string>{"Lg = 30nm  # Gate length"});
    }
    SECTION("a derived row's expression replaced as typed") {
        const auto out = dims::set_dimension(mos2d, "x_d", "$Lsd + 2nm");
        REQUIRE(out);
        CHECK(changed(mos2d, *out) == std::vector<std::string>{"x_d = $Lsd + 2nm"});
    }
    SECTION("a number") {
        const auto out = dims::set_dimension(mos2d, "Nsd", "5e19");
        REQUIRE(out);
        CHECK(changed(mos2d, *out) == std::vector<std::string>{"Nsd = 5e19  # S/D doping"});
    }
    SECTION("the same value: the text unchanged") {
        const auto out = dims::set_dimension(mos2d, "Lg", "20nm");
        REQUIRE(out);
        CHECK(*out == mos2d);
    }
}

TEST_CASE("dimensions: a refused value names the parameter and writes nothing", "[app][dimensions]") {
    for (const char* typed : {"abc", "$nope", "20 nm nm", ""}) {
        INFO(typed);
        const auto out = dims::set_dimension(mos2d, "Lg", typed);
        REQUIRE_FALSE(out);
        CHECK(out.error().find("Lg") != std::string::npos);
    }
    CHECK_FALSE(dims::set_dimension(mos2d, "nope", "1nm"));
}

TEST_CASE("dimensions: the wheel steps a number by a tenth of a decade, a length on a 1-2-5 ladder", "[app][dimensions]") {
    const auto rows = dims::read_dimensions(mos2d);
    REQUIRE(rows);
    dims::Dimension nsd = (*rows)[3];
    const std::string up = dims::step_dimension(nsd, 1);
    CHECK(std::abs(number_of(up) / (1e20 * std::pow(10.0, 0.1)) - 1.0) < 1e-12);
    nsd.expr = up;
    nsd.value = number_of(up);
    CHECK(std::abs(number_of(dims::step_dimension(nsd, -1)) / 1e20 - 1.0) < 1e-12);

    dims::Dimension lg = (*rows)[0];
    const auto length = [&](const std::string& expr) {
        lg.expr = expr;
        lg.value = static_cast<double>(*geo::parse_length(expr)) * geo::quantum_cm;
        return lg;
    };
    CHECK(dims::step_dimension(length("20nm"), 1) == "22nm");
    CHECK(dims::step_dimension(length("22nm"), -1) == "20nm");
    CHECK(dims::step_dimension(length("1nm"), 1) == "1.1nm");
    CHECK(dims::step_dimension(length("0.1nm"), -1) == "0.1nm");  // never below one lattice step
    CHECK(dims::step_dimension(length("2um"), 1) == "2.2um");     // the row's unit kept
    CHECK(dims::step_dimension(length("20nm"), 3) == "26nm");

    // A derived row is not stepped (its expression would be lost).
    CHECK(dims::step_dimension((*rows)[2], 1) == "$Lsd + $Lg");
}

TEST_CASE("dimensions: the 3D modeller keeps a parameter line's comment when it rewrites the line", "[app][dimensions][parameters]") {
    const std::string deck = R"([parameters]
x0 = 30nm
Lg = 20nm  # Gate length

[mesh]
type = generate3d
layered = on
domain = 0 100nm 0 60nm 0 14nm
patch body = 0 0 0  100nm 60nm 0
max_size = 0.05um
interface_resolution = 2nm

[region si]
material = silicon

[region gate]
material = sio2
box = $x0 $x0 + $Lg 40nm 50nm 0 14nm

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";
    auto m = solid::read_model(deck);
    REQUIRE(m);
    REQUIRE(solid::set_parameter(*m, "Lg", "25nm"));
    const auto out = solid::write_model(deck, *m);
    REQUIRE(out);
    CHECK(changed(deck, *out) == std::vector<std::string>{"Lg = 25nm  # Gate length"});
}

// --- Task 2: dimension spans (DT2) ----------------------------------------------------------------------------------

namespace {

bool near(double a, double b) { return std::abs(a - b) <= 1e-12 + 1e-9 * std::abs(b); }

const std::string pn1d = R"([parameters]
Ln = 0.8um  # n side
Lp = 1.2um  # p side

[mesh]
type = tensor
x = linspace(0, $Ln, 81) + linspace($Ln, $Ln + $Lp, 121)

[region silicon]
material = silicon

[doping]
profile = donor constant 1e18 x=0:$Ln
profile = acceptor constant 1e16

[contact cathode]
kind = ohmic
patch = x_min

[contact anode]
kind = ohmic
patch = x_max

[equilibrium]
)";

const std::string mos3d = R"([parameters]
tox = 2nm  # Oxide
L = 100nm

[mesh]
type = tensor
x = linspace(0, $L, 11)
y = linspace(0, 40nm, 9) + linspace(40nm, 40nm + $tox, 3)
z = linspace(0, 30nm, 4)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0 $L  40nm 40nm + $tox  0 30nm

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = y_min

[equilibrium]
)";

}  // namespace

TEST_CASE("dimension spans: a length's line is the shortest extent it changes with one end fixed", "[app][dimensions]") {
    SECTION("2D: the gate length on the gate contact") {
        const auto s = dims::dimension_spans(mos2d, "Lg");
        REQUIRE(s.size() == 1);
        CHECK(s[0].axis == 0);
        CHECK(near(s[0].lo, 5e-6));
        CHECK(near(s[0].hi, 7e-6));
        CHECK(s[0].source == "contact gate");
        CHECK(near(s[0].at[0], 5e-6));
        CHECK(near(s[0].at[1], 4.2e-6));
    }
    SECTION("2D: a derived edge, the gate's far side") {
        const auto s = dims::dimension_spans(mos2d, "x_d");
        REQUIRE(s.size() == 1);
        CHECK(s[0].source == "contact gate");
        CHECK(near(s[0].hi, 7e-6));
    }
    SECTION("a number has no line") { CHECK(dims::dimension_spans(mos2d, "Nsd").empty()); }
    SECTION("an unknown name has none") { CHECK(dims::dimension_spans(mos2d, "nope").empty()); }
    SECTION("1D: the p side from the mesh's second piece") {
        const auto s = dims::dimension_spans(pn1d, "Lp");
        REQUIRE(s.size() == 1);
        CHECK(s[0].axis == 0);
        CHECK(near(s[0].lo, 0.8e-4));
        CHECK(near(s[0].hi, 2.0e-4));
        CHECK(s[0].source == "mesh x");
    }
    SECTION("1D: the n side from the doping window (before the mesh at a tie)") {
        const auto s = dims::dimension_spans(pn1d, "Ln");
        REQUIRE(s.size() == 1);
        CHECK(near(s[0].hi, 0.8e-4));
        CHECK(s[0].source == "doping 1");
    }
    SECTION("3D: the oxide thickness along y on the oxide box") {
        const auto s = dims::dimension_spans(mos3d, "tox");
        REQUIRE(s.size() == 1);
        CHECK(s[0].axis == 1);
        CHECK(near(s[0].lo, 4e-6));
        CHECK(near(s[0].hi, 4.2e-6));
        CHECK(s[0].source == "region oxide");
        CHECK(near(s[0].at[2], 0.0));
    }
}
