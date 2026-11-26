// Hurkx trap-assisted tunnelling in decks (Unit 67, DECISIONS.md TAC8, TAC9, TAC10): the diode device checks,
// the [models] values and masses, the errors naming the values and the refused combination, the
// run record and identity, and a stopped transient with the model on continued bit for bit.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <stop_token>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_deck;
using deck_test::diode_device;
using deck_test::must;
using deck_test::plan_of;
using deck_test::scratch;
namespace fs = std::filesystem;

namespace {

std::string failure(const std::string& text) {
    const auto p = plan_of(text);
    REQUIRE_FALSE(p.has_value());
    return p.error().message;
}

double setting(const results::RunRecord& r, const std::string& key) {
    for (const auto& [k, v] : r.settings) {
        if (k == key) return v;
    }
    return -1.0;
}

bool named(const results::RunRecord& r, const std::string& key) {
    return std::ranges::any_of(r.settings, [&](const auto& s) { return s.first == key; });
}

bool same(double a, double b) {
    return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}
bool same(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::ranges::equal(a, b, [](double x, double y) { return same(x, y); });
}

const std::string transient_body = R"([transient]
contact anode = pulse(0, -2, 100ps, 200ps, 50ps, 50ps)
t_end = 600ps
dt_initial = 0.1ps
)";

const std::string tat_models = "[models]\ntrap_assisted_tunnelling = hurkx\n";

deck::RunReport run_deck(const std::string& text, const fs::path& dir, std::size_t stop_after = 0) {
    const deck::Plan plan = must(text, dir);
    std::stop_source source;
    solve::RunControl control{.stop = source.get_token()};
    std::size_t accepted = 0;
    double last_time = 0.0;
    if (stop_after > 0) {
        control.progress = [&](const solve::Progress& p) {
            if (p.phase != solve::Phase::transient || !p.converged) return;
            if (p.time_s > last_time) {
                last_time = p.time_s;
                if (++accepted == stop_after) source.request_stop();
            }
        };
    }
    auto r = deck::run(plan, control);
    if (!r) FAIL(r.error().message);
    REQUIRE(r->transient.has_value());
    return std::move(*r);
}

}  // namespace

TEST_CASE("deck build: trap_assisted_tunnelling and its masses (Unit 67, TAC9)", "[deck]") {
    const auto off = must(diode_deck).options.models;
    CHECK(off.trap_assisted_tunnelling == assemble::TrapTunnelling::off);
    CHECK(off.tat_mass_n == 0.25);
    CHECK(off.tat_mass_p == 0.25);
    const auto on = must(diode_deck + tat_models + "tat_mass_n = 0.3\ntat_mass_p = 0.4\n").options.models;
    CHECK(on.trap_assisted_tunnelling == assemble::TrapTunnelling::hurkx);
    CHECK(on.tat_mass_n == 0.3);
    CHECK(on.tat_mass_p == 0.4);
    CHECK(must(diode_deck + "[models]\ntrap_assisted_tunnelling = off\n").options.models.trap_assisted_tunnelling ==
          assemble::TrapTunnelling::off);
    CHECK(failure(diode_deck + "[models]\ntrap_assisted_tunnelling = schenk\n")
              .find("trap_assisted_tunnelling is off or hurkx") != std::string::npos);
}

TEST_CASE("deck run: trap_assisted_tunnelling without srh is refused by name (Unit 67, TAC9)", "[deck]") {
    const deck::Plan p = must(diode_deck + tat_models + "srh = off\n");
    const auto r = deck::run(p);
    REQUIRE_FALSE(r.has_value());
    CAPTURE(r.error().message);
    CHECK(r.error().message.find("srh") != std::string::npos);
}

TEST_CASE("deck run: trap_assisted_tunnelling in the run record and identity (Unit 67, TAC10)", "[deck]") {
    deck::Plan off = must(diode_deck);
    deck::Plan on = must(diode_deck + tat_models + "tat_mass_p = 0.5\n");
    off.outputs = {};
    on.outputs = {};
    const auto a = deck::run(off);
    const auto b = deck::run(on);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    const auto& ra = a->sweep->run;
    const auto& rb = b->sweep->run;
    for (const char* key : {"models.trap_assisted_tunnelling", "models.tat_mass_n", "models.tat_mass_p"}) {
        CAPTURE(key);
        CHECK_FALSE(named(ra, key));
        CHECK(named(rb, key));
    }
    CHECK(setting(rb, "models.trap_assisted_tunnelling") == 1.0);
    CHECK(setting(rb, "models.tat_mass_n") == 0.25);
    CHECK(setting(rb, "models.tat_mass_p") == 0.5);
    CHECK(ra.input_identity != rb.input_identity);
}

TEST_CASE("deck transient restart: with trap-assisted tunnelling, bit for bit (Unit 67, TAC8)",
          "[deck][restart]") {
    const fs::path dir = scratch("tat_transient_restart");
    const std::string deck = diode_device + tat_models + transient_body;
    const auto full = run_deck(deck, dir);
    REQUIRE_FALSE(full.transient->stopped.has_value());
    const auto& f = *full.transient;
    const auto first = run_deck(deck + "[output]\ncheckpoint = part.nitckpt\n", dir, 20);
    REQUIRE(first.transient->stopped.has_value());
    const std::size_t k = first.transient->points.size() - 1;
    const auto rest = run_deck(deck + "restart = part.nitckpt\n", dir);
    REQUIRE_FALSE(rest.transient->stopped.has_value());
    const auto& r = *rest.transient;
    CHECK(r.run.input_identity == f.run.input_identity);
    REQUIRE(r.points.size() == f.points.size() - k);
    for (std::size_t i = 0; i < r.points.size(); ++i) {
        INFO("point " << k + i);
        CHECK(same(r.points[i].time_s, f.points[k + i].time_s));
        CHECK(same(r.points[i].terminal_current, f.points[k + i].terminal_current));
        CHECK(same(r.points[i].contact_charge, f.points[k + i].contact_charge));
    }
}

namespace {

// TAC8 device checks: the deck tests' 1D diode, plain SRH and with Hurkx TAT. The reverse current
// is read as the weighted terminal current (Unit 32): the contact-edge current's rounding
// (BiasPoint::terminal_current_resolution, about 4e-7 A/cm^2 here) is far above the
// generation current (about 5e-9 A/cm^2).
struct Diode {
    deck::Plan plan;
    solve::BiasOptions srh, tat;
};

Diode tac8_diode() {
    Diode d{must(diode_deck)};
    d.srh = d.plan.options;
    d.srh.weighted_current = true;
    d.tat = d.srh;
    d.tat.models.trap_assisted_tunnelling = assemble::TrapTunnelling::hurkx;
    return d;
}

std::vector<std::vector<double>> ramp(double to, double step) {
    std::vector<std::vector<double>> v;
    const int n = static_cast<int>(std::lround(std::abs(to) / step));
    for (int k = 0; k <= n; ++k) v.push_back({std::copysign(step * k, to), 0.0});
    return v;
}

}  // namespace

TEST_CASE("Unit 67 TAC8: the reverse-biased diode leaks more with TAT, the ratio rising", "[deck]") {
    const Diode d = tac8_diode();
    const auto pts = ramp(-5.0, 0.25);
    const auto s = solve::sweep_bias(d.plan.device, pts, d.srh);
    const auto t = solve::sweep_bias(d.plan.device, pts, d.tat);
    REQUIRE(s.has_value());
    REQUIRE(t.has_value());
    REQUIRE_FALSE(s->stopped);
    REQUIRE_FALSE(t->stopped);
    REQUIRE(t->points.size() == pts.size());
    double last = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const auto &ps = s->points[i], &pt = t->points[i];
        const double Is = ps.weighted_current[0], It = pt.weighted_current[0];
        const std::size_t ns = ps.convergence.iterations.size(), nt = pt.convergence.iterations.size();
        CAPTURE(pts[i][0], Is, It, ns, nt);
        CHECK(pt.convergence.converged);
        CHECK(nt <= ns + 2);
        if (i == 0) continue;  // 0 V: both currents are zero to rounding
        // Each current is resolved: its weighted resolution far below it.
        CHECK(ps.weighted_current_resolution[0] < 1e-6 * std::abs(Is));
        CHECK(pt.weighted_current_resolution[0] < 1e-6 * std::abs(It));
        CHECK(std::abs(It) >= std::abs(Is));
        const double ratio = It / Is;
        CHECK(ratio > last);
        last = ratio;
        if (i % 4 == 0)
            std::printf("TAC8 V %+.2f  J srh %+.5e  tat %+.5e A/cm^2  ratio %.4f  Newton %zu / %zu\n",
                        pts[i][0], Is, It, ratio, ns, nt);
    }
}

TEST_CASE("Unit 67 TAC8: the forward-biased diode, recorded", "[deck]") {
    const Diode d = tac8_diode();
    const auto pts = ramp(0.7, 0.1);
    const auto s = solve::sweep_bias(d.plan.device, pts, d.srh);
    const auto t = solve::sweep_bias(d.plan.device, pts, d.tat);
    REQUIRE(s.has_value());
    REQUIRE(t.has_value());
    REQUIRE_FALSE(t->stopped);
    for (std::size_t i = 3; i < pts.size(); ++i) {
        const double Is = s->points[i].weighted_current[0], It = t->points[i].weighted_current[0];
        CHECK(t->points[i].convergence.converged);
        std::printf("TAC8 forward V %.1f  J srh %+.5e  tat %+.5e A/cm^2 (%+.3f%%)\n", pts[i][0], Is, It,
                    100.0 * (It / Is - 1.0));
    }
}

TEST_CASE("Unit 67 TAC8: small-signal conductance at -2 V against the DC slope", "[deck]") {
    const Diode d = tac8_diode();
    const auto sweep = solve::sweep_bias(d.plan.device, ramp(-2.0, 0.25), d.tat);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped);
    const results::BiasPoint& dc = sweep->points.back();
    const double h = 0.05;
    const auto I = [&](double V) {
        return solve::solve_bias(d.plan.device, std::vector<double>{V, 0.0}, d.tat, &dc.fields)
            ->weighted_current[0];
    };
    const double G =
        (8.0 * (I(-2.0 + h) - I(-2.0 - h)) - (I(-2.0 + 2 * h) - I(-2.0 - 2 * h))) / (12.0 * h);
    solve::BiasOptions steady = d.tat;
    steady.weighted_current = false;
    const auto ac = solve::solve_small_signal(d.plan.device, std::vector<std::vector<double>>{{-2.0, 0.0}},
                                              {.steady = steady, .frequencies_Hz = {0.0}}, &dc.fields);
    REQUIRE(ac.has_value());
    REQUIRE(!ac->stopped);
    const double Y = ac->conductance(0, 0, 0, 0);
    std::printf("TAC8 small signal at -2 V: G(DC) %+.6e  Y %+.6e S/cm^2 (%+.2e)  Y resolution %.3e\n", G, Y,
                Y / G - 1.0, ac->points[0].resolution[0][0]);
    // Owner ruling (2026-11-25): the admittance has no weighted counterpart, so at this conductance
    // it is compared within its own reported resolution, not to 1e-4 relative.
    CHECK(std::abs(Y - G) <= ac->points[0].resolution[0][0]);
}
