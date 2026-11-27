// Unit 68 slice 3 (MAD10; DECISIONS.md amendment 7): the III-V names in decks: `ingaas` with no
// fraction is the legacy set, unrecorded and bit-identical (MAC1); with one, the alloy; a material
// without SRH lifetimes is refused while srh is on; every new set solves a p-n junction.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/physics/iii_v.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;
using deck_test::plan_of;

namespace {

std::string diode(const std::string& material, const std::string& models = "", const std::string& run = "[equilibrium]\n") {
    return R"([mesh]
type = tensor
x = linspace(0, 2um, 201)

[region body]
material = )" + material + R"(

[doping]
profile = acceptor constant 1e17 x=-1:1um
profile = donor constant 1e17 x=1um:3um

[contact anode]
kind = ohmic
patch = x_min

[contact cathode]
kind = ohmic
patch = x_max

[models]
)" + models + "\n" + run;
}

const physics::SemiconductorParameters& used(const deck::Plan& p) {
    return std::get<physics::Semiconductor>(p.device.regions()[0].material).parameters();
}

}  // namespace

TEST_CASE("III-V decks: ingaas is the legacy set without a fraction, the alloy with one", "[deck][unit68][iii_v]") {
    const deck::Plan legacy = must(diode("ingaas"));
    CHECK(used(legacy) == physics::indium_gallium_arsenide_parameters);
    CHECK(legacy.materials.empty());  // not recorded: its result files do not change (MAC1)
    const deck::Plan alloy = must(diode("ingaas 0.53", "srh = off"));
    CHECK(used(alloy) == *physics::ingaas_parameters(0.53));
    REQUIRE(alloy.materials.size() == 1);
    CHECK(alloy.materials[0].material == "ingaas");
    CHECK(alloy.materials[0].composition == std::vector<double>{0.53});
    const auto two = plan_of(diode("ingaas 0.5 0.5"));
    REQUIRE_FALSE(two.has_value());
    CHECK(two.error().message.find("ingaas takes the In fraction") != std::string::npos);
    const auto out = plan_of(diode("inalas 1.2", "srh = off"));
    REQUIRE_FALSE(out.has_value());
    CHECK(out.error().message.find("[0, 1]") != std::string::npos);
    // A graded region names the alloy entry.
    const std::string graded = R"([mesh]
type = tensor
x = linspace(0, 1um, 11)
[region body]
material = ingaas
composition = linear 0.4 0.6 x 0 1um
[contact a]
kind = ohmic
patch = x_min
[models]
srh = off
[equilibrium]
)";
    CHECK(plan_of(graded).has_value());
}

TEST_CASE("III-V decks: no SRH lifetimes is refused while srh is on", "[deck][unit68][iii_v]") {
    for (const char* m : {"insb", "ingaas 0.53", "inalas 0.52"}) {
        CAPTURE(m);
        const std::string sweep = "[sweep]\ncontact anode = 0 : 0.1 : 0.1\n";  // equilibrium needs no SRH
        const auto r = deck::run(must(diode(m, "", sweep)));
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("region 'body' has no SRH lifetimes") != std::string::npos);
        CHECK(deck::run(must(diode(m, "srh = off", sweep))).has_value());
    }
}

TEST_CASE("III-V decks: every new set solves a p-n junction at equilibrium", "[deck][unit68][iii_v]") {
    for (const char* m : {"gaas_2001", "inp", "inas", "gap", "alas", "gasb", "insb", "ingaas 0.53", "inalas 0.52",
                          "inalas 0.2"}) {
        CAPTURE(m);
        const bool lifetimes = std::string(m) != "insb" && std::string(m).find(' ') == std::string::npos;
        const auto r = deck::run(must(diode(m, lifetimes ? "" : "srh = off")));
        if (!r) FAIL(r.error().message);
        for (const double n : r->equilibrium->fields.n_cm3) REQUIRE((std::isfinite(n) && n > 0.0));
        for (const double p : r->equilibrium->fields.p_cm3) REQUIRE((std::isfinite(p) && p > 0.0));
    }
}
