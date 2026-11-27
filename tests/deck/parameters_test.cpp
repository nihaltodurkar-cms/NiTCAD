// Unit 26m slice 3 (MT8, MC7): [parameters] and `$name` references.

#include <catch2/catch_test_macros.hpp>

#include <format>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/io/files.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_deck;
using deck_test::must;
using deck_test::plan_of;

namespace {

bool has(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

std::string failure(const std::string& text) {
    const auto p = plan_of(text);
    REQUIRE_FALSE(p.has_value());
    return p.error().message;
}

std::vector<deck::Parameter> params_of(const std::string& text) {
    const auto d = deck::parse(text, "t.deck");
    REQUIRE(d);
    const auto p = deck::parameters(*d);
    REQUIRE(p);
    return *p;
}

std::string replaced(const std::string& value, const std::string& params) {
    const auto p = params_of("[parameters]\n" + params);
    const auto r = deck::substitute(value, p);
    REQUIRE(r);
    return *r;
}

// Two plans the same bits: the run identity (mesh, device, options, points) and the mesh digest.
void same_bits(const deck::Plan& a, const deck::Plan& b) {
    CHECK(solve::make_run_record(a.device, a.options, a.points).input_identity ==
          solve::make_run_record(b.device, b.options, b.points).input_identity);
    CHECK(io::mesh_digest(a.device.mesh()) == io::mesh_digest(b.device.mesh()));
}

}  // namespace

TEST_CASE("parameters: a deck with references builds to the bits of the values written out (MC7)", "[deck][parameters]") {
    SECTION("the diode: [mesh], [doping], a number parameter") {
        const std::string refs = R"(
[parameters]
L = 2um
n = 201
N = 1e17
mid = 1um

[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, $L, $n)

[region silicon]
material = silicon

[doping]
profile = acceptor constant $N x=-1:$mid
profile = donor constant $N x=$mid:3um

[contact anode]
kind = ohmic
patch = x_min

[contact cathode]
kind = ohmic
patch = x_max

[sweep]
contact anode = 0 : 0.2 : 0.1
)";
        same_bits(must(refs), must(diode_deck));
    }
    SECTION("a layered generate2d deck: solids, the domain, patches, sums and tied faces") {
        const std::string tail = R"(
[doping]
profile = acceptor constant 1e16

[contact source]
kind = ohmic
patch = source

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";
        const std::string written = R"(
[mesh]
type = generate2d
layered = on
domain = 0 2um 0 1um
patch source = 0 0, 0 0.6um
patch body = 0.2um 0, 0.8um 0
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon

[region n]
material = silicon
polygon = 1um 0, 2um 0, 2um 0.5um, 2um 0.8um, 1um 0.8um

[region oxide]
material = sio2
box = 0 2um 0.8um 1um
)" + tail;
        // half + half is 1um to the bit (0.5 * 1e-4 is exact, so is their sum); x0 + w the same.
        const std::string refs = R"(
[parameters]
W = 2um
H = 1um
half = 0.5um
x0 = $half
w = $half
tox = 0.2um
top = $H - $tox
s = 20nm
h = 0.25um

[mesh]
type = generate2d
layered = on
domain = 0 $W 0 $H
patch source = 0 0, 0 0.6um
patch body = 0.2um 0, 0.8um 0
max_size = $h
interface_resolution = $s

[region p]
material = silicon

[region n]
material = silicon
polygon = $half + $half 0, $W 0, $W 0.5um, $W $top, $x0 + $w $top

[region oxide]
material = sio2
box = 0 $W $top $H
)" + tail;
        // top = 1um - 0.2um is not 0.8um to the bit: compare with the value the substitution writes out.
        const auto p = params_of(refs);
        const deck::Parameter* top = nullptr;
        for (const auto& q : p) {
            if (q.name == "top") top = &q;
        }
        REQUIRE(top != nullptr);
        CHECK(top->length);
        CHECK(top->text == std::format("{}cm", 1e-4 - 0.2 * 1e-4));
        std::string out = written;
        const auto swap = [&](const std::string& from, const std::string& to) {
            const auto at = out.find(from);
            REQUIRE(at != std::string::npos);
            out.replace(at, from.size(), to);
        };
        swap("box = 0 2um 0.8um 1um", "box = 0 2um " + top->text + " 1um");
        swap("2um 0.8um, 1um 0.8um", "2um " + top->text + ", 1um " + top->text);
        same_bits(must(refs), must(out));
    }
}

TEST_CASE("parameters: substitution rules (MQ1 (a))", "[deck][parameters]") {
    CHECK(replaced("box = 0 $L", "L = 20nm") == "box = 0 20nm");                       // a lone reference: as written
    CHECK(replaced("x=$a:$b", "a = 1um\nb = 3um") == "x=1um:3um");
    CHECK(replaced("linspace(0, $L, $n)", "L = 2um\nn = 9") == "linspace(0, 2um, 9)");
    CHECK(replaced("$a + $b", "a = 1um\nb = 1um") == std::format("{}cm", 2e-4));       // a sum: in cm
    CHECK(replaced("$a-5nm", "a = 1um") == std::format("{}cm", 1e-4 - 5 * 1e-7));      // no spaces: an operator
    CHECK(replaced("5nm + $a", "a = 1um") == std::format("{}cm", 5 * 1e-7 + 1e-4));    // a literal first
    CHECK(replaced("$a -5nm", "a = 1um") == "1um -5nm");                               // space before only: two values
    CHECK(replaced("$a - 0", "a = 1um") == std::format("{}cm", 1e-4));                 // unitless in a sum: cm (DK5)
    CHECK(replaced("linspace(0, 2um, 9) + linspace(2um, $L, 9)", "L = 3um") == "linspace(0, 2um, 9) + linspace(2um, 3um, 9)");
    CHECK(replaced("no references 1um + 2um", "a = 1um") == "no references 1um + 2um");
    // A parameter of parameters, followed in any order.
    const auto p = params_of("[parameters]\nb = $a + 1um\na = 1um\n");
    REQUIRE(p.size() == 2);
    CHECK(p[0].name == "a");
    CHECK(p[1].name == "b");
    CHECK(p[1].value == 2e-4);
    CHECK(p[1].length);
    CHECK(params_of("[parameters]\nn = 9\n")[0].length == false);
    // A sum of lengths alone: its value written out, a length.
    const auto lit = params_of("[parameters]\ns = 10nm + 20nm\n");
    CHECK(lit[0].length);
    CHECK(lit[0].text == std::format("{}cm", 10 * 1e-7 + 20 * 1e-7));
    CHECK(lit[0].value == 10 * 1e-7 + 20 * 1e-7);
}

TEST_CASE("parameters: errors on the referring line (MC7)", "[deck][parameters]") {
    const std::string body = "\n[device]\ntemperature = 300\n\n[mesh]\ntype = tensor\nx = linspace(0, 2um, 201)\n";
    // undefined, on the [mesh] line (line 7 of body)
    CHECK(has(failure(body + "y = $nope\n"), "t.deck:8: undefined parameter '$nope'"));
    // a cycle, on the parameter's line
    const auto cycle = failure("[parameters]\na = $b\nb = $a\n" + body);
    CHECK(has(cycle, "t.deck:"));
    CHECK(has(cycle, "a cycle of parameters: $a -> $b -> $a"));
    CHECK(has(failure("[parameters]\na = $a\n" + body), "t.deck:2: a cycle of parameters: $a -> $a"));
    // a number in a sum
    CHECK(has(failure("[parameters]\nn = 9\n" + body + "y = $n + 1um\n"), "t.deck:10: '$n' is a number (9), not a length"));
    // neither a length nor a number
    CHECK(has(failure("[parameters]\nm = silicon\n" + body), "t.deck:2: parameter 'm' = 'silicon' is neither"));
    // a bad name, twice, a second section, a name on it
    CHECK(has(failure("[parameters]\n2x = 1um\n" + body), "t.deck:2: malformed key"));  // the parser's own
    CHECK(has(failure("[parameters]\na.b = 1um\n" + body), "t.deck:2: 'a.b' is not a parameter name"));
    CHECK(has(failure("[parameters]\na = 1um\na = 2um\n" + body), "t.deck:3: parameter 'a' given twice (first on line 2)"));
    CHECK(has(failure("[parameters]\na = 1um\n[parameters]\nb = 1um\n" + body), "t.deck:3: a second [parameters] section"));
    CHECK(has(failure("[parameters p]\na = 1um\n" + body), "t.deck:1: [parameters] takes no name"));
    // outside [mesh], [region], [contact], [doping] and [sweep]
    CHECK(has(failure("[parameters]\nT = 300\n\n[device]\ntemperature = $T\n"), "t.deck:5: a parameter reference outside [mesh], [region], [contact], [doping] and [sweep]"));
    // a '$' not starting a reference
    CHECK(has(failure(body + "y = 1um$a\n"), "t.deck:8: "));
    CHECK(has(failure(body + "y = $ 1um\n"), "a '$' not starting a reference"));
}

// Unit 26p (GU1): [sweep] takes $name, so a bias can be a parameter.
namespace {
std::string sweep_value(const std::string& text, const std::string& key) {
    const auto d = deck::parse(text, "t.deck");
    REQUIRE(d);
    const auto s = deck::substituted(*d);
    REQUIRE(s);
    for (const auto& sec : s->sections) {
        if (sec.kind != "sweep") continue;
        for (const auto& e : sec.entries) {
            if (e.key == key) return e.value;
        }
    }
    FAIL("no such [sweep] entry");
    return {};
}
}  // namespace

TEST_CASE("parameters: a number parameter substitutes in [sweep] (GU1)", "[deck][parameters]") {
    CHECK(sweep_value("[parameters]\nVg = 0.8\n\n[sweep]\ncontact gate = $Vg\n", "contact gate") == "0.8");
}

TEST_CASE("parameters: a range in [sweep] takes parameters (GU1)", "[deck][parameters]") {
    CHECK(sweep_value("[parameters]\nVd = 1.2\n\n[sweep]\ncontact drain = 0 : $Vd : 0.05\n", "contact drain") ==
          "0 : 1.2 : 0.05");
}

TEST_CASE("parameters: a missing parameter in [sweep] is refused on its line (GU1)", "[deck][parameters]") {
    const auto d = deck::parse("[parameters]\nVd = 1.2\n\n[sweep]\ncontact drain = $nope\n", "t.deck");
    REQUIRE(d);
    const auto s = deck::substituted(*d);
    REQUIRE_FALSE(s);
    CHECK(has(s.error().message, "t.deck:5: undefined parameter '$nope'"));
}
