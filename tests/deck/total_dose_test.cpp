// Total-dose oxide charge in decks (Unit 69, DECISIONS.md TDC3-TDC5): a gate's and an interface's
// dose build the fixed charge and trap bands written out by hand, and the dose keys' refusals.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/physics/total_dose.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;
using deck_test::plan_of;

namespace {

std::string failure(const std::string& text) {
    const auto p = plan_of(text);
    REQUIRE_FALSE(p.has_value());
    return p.error().message;
}

bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
bool same(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::ranges::equal(a, b, [](double x, double y) { return same(x, y); });
}

std::string g17(double v) {
    char s[64];
    std::snprintf(s, sizeof s, "%.17g", v);
    return s;
}

// A 2D compact-gate MOSFET (cross_run_learning_test's family member L_g 400 nm, t_ox 5 nm).
std::string mosfet(const std::string& gate_keys) {
    return R"([mesh]
type = tensor
x = linspace(0, 1000nm, 21)
y = linspace(0, 200nm, 11)

[region silicon]
material = silicon

[doping]
profile = acceptor constant 1e17
profile = donor gaussian 1e19 axis=y center=0 width=50nm x=0:300nm lateral=14.142135624nm
profile = donor gaussian 1e19 axis=y center=0 width=50nm x=700nm:1000nm lateral=14.142135624nm

[contact source]
kind = ohmic
patch = y_min
box = 0 300.5nm  0 0

[contact drain]
kind = ohmic
patch = y_min
box = 699.5nm 1000nm  0 0

[contact body]
kind = ohmic
patch = y_max

[contact gate]
kind = gate
patch = y_min
box = 301nm 699nm  0 0
)" + gate_keys + R"(
[sweep]
contact drain = 0.05
contact gate = 0 : 1 : 0.5
)";
}

// tests/decks/oxide_traps_2d.deck's device with its [interface] keys replaced.
std::string oxide(const std::string& interface_keys) {
    return R"([mesh]
type = tensor
x = linspace(0, 2um, 41)
y = linspace(0, 0.9um, 19) + linspace(0.9um, 0.95um, 6) + linspace(0.95um, 1um, 3)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0 2um  0.95um 1um

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = y_min

[contact gate]
kind = electrode
patch = y_max
box = 0.5um 1.5um  1um 1um
electrode = metal
work_function = 4.1

[interface silicon oxide]
)" + interface_keys + R"(
[equilibrium]
)";
}


const device::Interface& only_interface(const deck::Plan& p) {
    REQUIRE(p.device.interfaces().size() == 1);
    return p.device.interfaces()[0];
}

}  // namespace

TEST_CASE("TDC3: a gate's dose is its fixed charge", "[total_dose]") {
    const double n_ot = physics::oxide_trapped_charge_cm2(1e5, 1e6, 5 * 1e-7, 0.3);  // 5nm as the deck reads it
    // dose before oxide_thickness: the gate's t_ox is used wherever it is written.
    const deck::Plan a = must(mosfet("fixed_charge = 1e10\ndose = 1e5\ndose_field = 1e6\nhole_trapping = 0.3\n"
                                     "oxide_thickness = 5nm\n"));
    const deck::Plan b = must(mosfet("oxide_thickness = 5nm\nfixed_charge = " + g17(1e10 + n_ot) + "\n"));
    const auto* ga = a.device.find_contact("gate");
    const auto* gb = b.device.find_contact("gate");
    REQUIRE(ga != nullptr);
    REQUIRE(gb != nullptr);
    CHECK(ga->gate.fixed_charge_cm2 > 1e10);
    CHECK(same(ga->gate.fixed_charge_cm2, gb->gate.fixed_charge_cm2));

    const auto ra = deck::run(a);
    const auto rb = deck::run(b);
    REQUIRE(ra.has_value());
    REQUIRE(rb.has_value());
    REQUIRE(ra->sweep->points.size() == 3);
    REQUIRE(rb->sweep->points.size() == 3);
    for (std::size_t k = 0; k < 3; ++k) {
        const auto& pa = ra->sweep->points[k];
        const auto& pb = rb->sweep->points[k];
        CHECK(same(pa.fields.potential_V, pb.fields.potential_V));
        CHECK(same(pa.terminal_current, pb.terminal_current));
    }
}

TEST_CASE("TDC4: an interface's dose is its fixed charge and two trap bands", "[total_dose]") {
    const std::string user_band = "trap_band = donor 1e11 -0.4 -0.1 1e-15 2e-15\n";
    const double n_ot = physics::oxide_trapped_charge_cm2(1e5, 1e6, 50 * 1e-7, 0.3);  // 50nm as the deck reads it
    const deck::Plan a = must(oxide("fixed_charge = 5e10\n" + user_band +
                                    "dose = 1e5\ndose_field = 1e6\nhole_trapping = 0.3\n"
                                    "dose_oxide_thickness = 50nm\ninterface_trap_yield = 1e6\n"));
    const deck::Plan b = must(oxide("fixed_charge = " + g17(5e10 + n_ot) + "\n" + user_band +
                                    "trap_band = donor 1e11 -0.5 0\ntrap_band = acceptor 1e11 0 0.5\n"));
    const auto& fa = only_interface(a);
    const auto& fb = only_interface(b);
    CHECK(same(fa.fixed_charge_cm2, fb.fixed_charge_cm2));
    REQUIRE(fa.traps.bands.size() == 3);
    CHECK(fa.traps.bands[0].sigma_p_cm2 == 2e-15);  // the user's band first
    CHECK(fa.traps == fb.traps);
}

TEST_CASE("TDC4: hole_trapping 0 leaves the fixed charge", "[total_dose]") {
    const auto& f = only_interface(
        must(oxide("fixed_charge = 5e10\ndose = 1e5\nhole_trapping = 0\ndose_oxide_thickness = 50nm\n")));
    CHECK(f.fixed_charge_cm2 == 5e10);
    CHECK(f.traps.bands.empty());
}

TEST_CASE("TDC1: without a dose the fixed charge is as written, -0 too", "[total_dose]") {
    CHECK(std::signbit(only_interface(must(oxide("fixed_charge = -0\n"))).fixed_charge_cm2));
    const auto* g = must(mosfet("oxide_thickness = 5nm\nfixed_charge = -0\n")).device.find_contact("gate");
    REQUIRE(g != nullptr);
    CHECK(std::signbit(g->gate.fixed_charge_cm2));
}

TEST_CASE("TDC5: dose refusals", "[total_dose]") {
    // The [interface] section is line 27; its keys start at line 28.
    const auto iface = [](const std::string& keys) { return failure(oxide(keys)); };
    const auto has = [](const std::string& m, const std::string& line, const std::string& what) {
        INFO(m);
        CHECK(m.find("t.deck:" + line + ":") != std::string::npos);
        CHECK(m.find(what) != std::string::npos);
    };
    const std::string ok = "hole_trapping = 0.3\ndose_oxide_thickness = 50nm\n";
    has(iface("dose = 0\n" + ok), "28", "dose must be finite and positive");
    has(iface("dose = -1\n" + ok), "28", "dose must be finite and positive");
    has(iface("dose_field = 1e6\n"), "28", "'dose_field' needs 'dose'");
    has(iface("interface_trap_yield = 1e6\n"), "28", "'interface_trap_yield' needs 'dose'");
    has(iface("dose = 1e5\ndose_oxide_thickness = 50nm\n"), "28", "dose needs 'hole_trapping'");
    has(iface("dose = 1e5\nhole_trapping = 0.3\n"), "28", "dose needs 'dose_oxide_thickness'");
    has(iface("dose = 1e5\nhole_trapping = 1.5\ndose_oxide_thickness = 50nm\n"), "29", "hole_trapping must lie in 0..1");
    has(iface("dose = 1e5\nhole_trapping = -0.1\ndose_oxide_thickness = 50nm\n"), "29",
        "hole_trapping must lie in 0..1");
    has(iface("dose = 1e5\n" + ok + "dose_field = -1\n"), "31", "dose_field must be finite and >= 0");
    has(iface("dose = 1e5\n" + ok + "interface_trap_yield = -1\n"), "31", "interface_trap_yield must be finite and >= 0");

    // The gate section is line 28; its keys start at line 32.
    const std::string gate = "oxide_thickness = 5nm\ndose = 1e5\nhole_trapping = 0.3\n";
    has(failure(mosfet(gate + "interface_trap_yield = 1e6\n")), "35", "'interface_trap_yield' is for an [interface]");
    has(failure(mosfet(gate + "dose_oxide_thickness = 5nm\n")), "35", "'dose_oxide_thickness' is for an [interface]");
    has(failure(mosfet("dose = 1e5\nhole_trapping = 0.3\n")), "32", "dose needs the gate's oxide_thickness");

    // On an ohmic contact the dose keys are refused with the gate keys.
    std::string ohmic = mosfet("oxide_thickness = 5nm\n");
    ohmic.replace(ohmic.find("patch = y_max\n"), 14, "patch = y_max\ndose = 1e5\n");
    has(failure(ohmic), "24", "the dose keys are for a gate");
}
