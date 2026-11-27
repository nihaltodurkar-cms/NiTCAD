// Unit 68 slice 2 (MAD9; DECISIONS.md amendment 6): nitride polarization charge in decks.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <string>

#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/physics/polarization.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;
using deck_test::plan_of;

namespace {

// A 1D Al0.25Ga0.75N (25 nm) on GaN (1 um) stack, the surface at x_max.
std::string stack(const std::string& models, const std::string& barrier_keys = "") {
    return R"([mesh]
type = tensor
x = linspace(0, 0.9um, 46) + linspace(0.9um, 1.025um, 126)

[region buffer]
material = gan

[region barrier]
material = algan 0.25
box = 1um 1.025um
)" + barrier_keys + R"(
[doping]
profile = donor constant 1e16

[contact body]
kind = ohmic
patch = x_min

[contact top]
kind = ohmic
patch = x_max

[models]
)" + models + R"(
[equilibrium]
)";
}

std::string failure(const std::string& text) {
    const auto p = plan_of(text);
    REQUIRE_FALSE(p.has_value());
    return p.error().message;
}

double P_cm2(const char* m, double x, double a) {
    return physics::polarization_sheet_density_cm2(physics::nitride_polarization(*physics::nitride_lattice(m, x), a));
}

// The bound charge [cm^-2] on the nodes with x in [lo, hi].
double charge(const deck::Plan& p, double lo, double hi) {
    const auto b = p.device.bound_charge();
    const auto pts = p.device.mesh().points();
    const auto vol = p.device.mesh().volumes();
    double q = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        if (pts[i][0] >= lo && pts[i][0] <= hi) q += b[i] * vol[i];
    }
    return q;
}

}  // namespace

TEST_CASE("MAD9: Al0.25Ga0.75N pseudomorphic on GaN, Ambacher's polarization", "[deck][unit68][polarization]") {
    const physics::NitrideLattice g = *physics::nitride_lattice("gan", 0.0);
    const physics::NitrideLattice a = *physics::nitride_lattice("algan", 0.25);
    CHECK(std::abs(a.a0_A - 3.16975) < 1e-12);
    CHECK(std::abs(a.P_sp - (-0.0440625)) < 1e-12);
    const double sigma = physics::nitride_polarization(g, g.a0_A) - physics::nitride_polarization(a, g.a0_A);
    CHECK(std::abs(sigma - 0.019361) < 1e-6);  // by hand, amendment 6 values
    CHECK(std::abs(physics::polarization_sheet_density_cm2(sigma) / 1.2084e13 - 1.0) < 1e-3);
    CHECK_FALSE(physics::nitride_lattice("silicon", 0.0).has_value());
}

TEST_CASE("MAD9: the stack's bound charge, interface and free surfaces", "[deck][unit68][polarization]") {
    const deck::Plan off = must(stack(""));
    CHECK(off.device.bound_charge().empty());
    CHECK(must(stack("polarization = off\n")).device.bound_charge().empty());

    const deck::Plan on = must(stack("polarization = on\n"));
    const double a_gan = physics::nitride_lattice("gan", 0.0)->a0_A;
    const double P_gan = P_cm2("gan", 0, a_gan), P_bar = P_cm2("algan", 0.25, a_gan);
    const double sigma = P_gan - P_bar;
    std::printf("MAD9 Al0.25Ga0.75N/GaN: sigma %.4e cm^-2\n", sigma);
    CHECK(std::abs(charge(on, 0.99e-4, 1.01e-4) / sigma - 1.0) < 1e-12);  // the interface
    CHECK(std::abs(charge(on, 1.0249e-4, 1.026e-4) / P_bar - 1.0) < 1e-12);  // the surface: P.n, n = +c
    CHECK(std::abs(charge(on, -1.0, 1e-7) / -P_gan - 1.0) < 1e-12);       // the bottom: n = -c
    CHECK(std::abs(charge(on, -1.0, 1.0)) < 1e-9 * std::abs(sigma));       // bound charge is neutral

    const deck::Plan half = must(stack("polarization = on\npolarization_scale = 0.5\n"));
    CHECK(std::abs(charge(half, 0.99e-4, 1.01e-4) / (0.5 * sigma) - 1.0) < 1e-12);
    const deck::Plan flipped = must(stack("polarization = on\nc_axis = -x\n"));
    CHECK(std::abs(charge(flipped, 0.99e-4, 1.01e-4) / -sigma - 1.0) < 1e-12);
    // Fully relaxed: the spontaneous part alone.
    const double a_bar = physics::nitride_lattice("algan", 0.25)->a0_A;
    const deck::Plan relaxed = must(stack("polarization = on\n", "relaxation = 1\n"));
    CHECK(std::abs(charge(relaxed, 0.99e-4, 1.01e-4) / (P_gan - P_cm2("algan", 0.25, a_bar)) - 1.0) < 1e-12);
}

TEST_CASE("MAD9: the 2DEG forms at equilibrium; off equals absent", "[deck][unit68][polarization]") {
    const auto electrons_near = [](const std::string& text) {
        const deck::Plan p = must(text);
        const auto r = deck::run(p);
        if (!r) FAIL(r.error().message);
        REQUIRE(r->equilibrium.has_value());
        const auto pts = p.device.mesh().points();
        const auto vol = p.device.mesh().volumes();
        double ns = 0.0;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            if (pts[i][0] > 0.95e-4 && pts[i][0] < 1e-4) ns += r->equilibrium->fields.n_cm3[i] * vol[i];
        }
        return ns;
    };
    const double ns_on = electrons_near(stack("polarization = on\n"));
    const double ns_off = electrons_near(stack(""));
    std::printf("MAD9 electrons in the 50 nm below the interface: on %.3e, off %.3e cm^-2\n", ns_on, ns_off);
    CHECK(ns_on > 5e12);
    CHECK(ns_off < 0.1 * ns_on);  // the band offset alone
    CHECK(electrons_near(stack("polarization = off\n")) == ns_off);
}

TEST_CASE("MAD9: polarization refusals", "[deck][unit68][polarization]") {
    const auto has = [](const std::string& m, const std::string& what) {
        INFO(m);
        CHECK(m.find(what) != std::string::npos);
    };
    has(failure(stack("polarization_scale = 0.5\n")), "'polarization_scale' needs polarization = on");
    has(failure(stack("c_axis = x\n")), "'c_axis' needs polarization = on");
    has(failure(stack("", "relaxation = 0.5\n")), "'relaxation' needs [models] polarization = on");
    has(failure(stack("polarization = on\npolarization_scale = 2\n")), "polarization_scale must lie in 0..1");
    has(failure(stack("polarization = on\nc_axis = y\n")), "c_axis is x, y or z");
    has(failure(stack("polarization = on\n", "relaxation = 1.5\n")), "relaxation must lie in 0..1");
    std::string si = stack("polarization = on\n");
    si.replace(si.find("material = gan"), 14, "material = silicon");
    si.replace(si.find("material = algan 0.25"), 21, "material = silicon");
    has(failure(si), "no region is a nitride");
}
