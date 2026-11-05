// TVB-32, the consistency suite on real solves of the three TVR-1 decks (DECISIONS Unit 55 section 11.6, TVR1-Q5 and
// the TVB-32 rulings of 11.13): each deck is built and run as `nitcad run` would (outputs off), and every check
// whose inputs the run supplies is applied at every operating point; the outcome of each is held to the expectation
// frozen in 11.13 before any run (pass, or NOT APPLICABLE with the named predicate). Any other outcome is a finding.
// C6, C9, C10 and C11 need inputs no deck produces (a transient, a DC stencil, a limit series, a symmetric pair):
// listed, not run (their real-solve evidence is in their own tests). Per check the worst residual/tolerance ratio is
// printed.
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace c = NiTCAD::verify::consistency;
namespace v = NiTCAD::verify;
namespace results = NiTCAD::results;
using v::Outcome;

namespace {

std::filesystem::path source(const char* file) { return std::filesystem::path(NITCAD_SOURCE_DIR) / file; }

// The frozen expectation for a check gated by C1's predicates: pass at exactly zero bias with one thermal
// temperature, else NOT APPLICABLE naming "all biases 0".
bool at_zero(const results::BiasPoint& p) {
    if (std::any_of(p.bias_V.begin(), p.bias_V.end(), [](double b) { return b != 0.0; })) return false;
    return p.thermal_bias_K.empty() ||
           std::all_of(p.thermal_bias_K.begin(), p.thermal_bias_K.end(),
                       [&](double t) { return t == p.thermal_bias_K.front(); });
}

// Per check: how many results had each outcome, and the worst judged ratio.
struct Tally {
    std::map<std::string, int> outcomes;
    double worst = 0.0;
    std::string where;
};
std::map<std::string, Tally> tally;

const char* name_of(Outcome o) {
    switch (o) {
        case Outcome::pass: return "pass";
        case Outcome::fail: return "FAIL";
        case Outcome::not_applicable: return "not applicable";
        case Outcome::not_evaluated: return "NOT EVALUATED";
        default: return "other";
    }
}

// Holds a result to its expectation and records it.
void expect(const c::CheckResult& r, Outcome want, const char* predicate, const std::string& at) {
    auto& t = tally[r.id];
    ++t.outcomes[name_of(r.comparison.outcome)];
    for (const auto& p : r.parts) {
        if (p.outcome != Outcome::pass || !(p.tolerance > 0.0)) continue;
        const double q = p.residual / p.tolerance;
        if (q > t.worst) {
            t.worst = q;
            t.where = at + ", " + p.name + " at " + p.where;
        }
    }
    CAPTURE(r.id, at, r.comparison.reason);
    for (const auto& p : r.parts) {
        CAPTURE(p.name, name_of(p.outcome), p.residual, p.tolerance, p.where, p.reason);
        if (p.outcome == Outcome::fail || p.outcome == Outcome::not_evaluated) FAIL_CHECK("part " << p.name);
    }
    CHECK(r.comparison.outcome == want);
    if (want == Outcome::not_applicable) CHECK(r.comparison.reason.find(predicate) != std::string::npos);
}

void report(const char* deck) {
    for (const auto& [id, t] : tally) {
        std::string counts;
        for (const auto& [o, n] : t.outcomes) counts += (counts.empty() ? "" : ", ") + std::to_string(n) + " " + o;
        std::printf("TVB-32 %s %s: %s; worst %.3g of the bound (%s)\n", deck, id.c_str(), counts.c_str(), t.worst,
                    t.where.c_str());
    }
    tally.clear();
}

// C1-C5 and C12 at one operating point, against the frozen expectations.
void dc_checks(const results::BiasPoint& p, const NiTCAD::mesh::Mesh& m, const c::DeviceData& d,
               const std::string& at) {
    const Outcome zero = at_zero(p) ? Outcome::pass : Outcome::not_applicable;
    expect(c::check_c1(p, m, d), zero, "all biases 0", at);
    expect(c::check_c2(p, m, d), zero, "all biases 0", at);
    expect(c::check_c3(p, m, d), Outcome::pass, "", at);
    expect(c::check_c4(p, m, d), Outcome::pass, "", at);
    expect(c::check_c5(p, m, d), Outcome::not_applicable, "electrothermal on", at);
    expect(c::check_c12(p.fields, &p, m, d), Outcome::pass, "", at);
}

std::string bias_label(const results::BiasPoint& p) {
    std::string s = "bias (";
    for (std::size_t i = 0; i < p.bias_V.size(); ++i) {
        char b[32];
        std::snprintf(b, sizeof b, "%s%.6g", i == 0 ? "" : ", ", p.bias_V[i]);
        s += b;
    }
    return s + ") V";
}

struct Built {
    NiTCAD::deck::Plan plan;
    c::DeviceData data;
};

Built built(const char* file) {
    auto deck = NiTCAD::deck::read(source(file));
    REQUIRE(deck.has_value());
    auto plan = NiTCAD::deck::build(*deck);
    REQUIRE(plan.has_value());
    plan->outputs = {};  // nothing written
    auto data = NiTCAD::verify::benchmark::device_data(plan->device, plan->options.models);
    return {std::move(*plan), std::move(data)};
}

void equilibrium_checks(const Built& b) {
    NiTCAD::solve::EquilibriumOptions eo;
    eo.models = b.plan.options.models;
    const auto eq = NiTCAD::solve::solve_equilibrium(b.plan.device, eo);
    REQUIRE(eq.has_value());
    expect(c::check_c1(*eq, b.plan.device.mesh(), b.data), Outcome::pass, "", "equilibrium");
    expect(c::check_c3(*eq, b.plan.device.mesh(), b.data), Outcome::pass, "", "equilibrium");
}

void small_signal_checks(const results::SmallSignal& run, const NiTCAD::mesh::Mesh& m, const c::DeviceData& d,
                         const std::string& prefix) {
    for (std::size_t k = 0; k < run.points.size(); ++k) {
        const auto& p = run.points[k];
        const std::string at = prefix + bias_label(p.dc);
        dc_checks(p.dc, m, d, at);
        const Outcome zero = at_zero(p.dc) ? Outcome::pass : Outcome::not_applicable;
        for (const auto& r : {c::check_c7(run, k, m, d), c::check_c8(run, k, m, d)}) {
            expect(r, zero, "all biases 0", at);
            if (zero == Outcome::pass && d.contacts.size() == 2) {
                CHECK(std::any_of(r.diagnostics.begin(), r.diagnostics.end(),
                                  [](const c::Diagnostic& g) { return g.name == "two contacts"; }));
            }
        }
    }
}

}  // namespace

TEST_CASE("TVB-32: the consistency suite on the DC sweeps of the pn and Schottky diode decks", "[verify][consistency]") {
    for (const char* file : {"examples/pn_diode_1d.deck", "decks/schottky_diode_1d.deck"}) {
        CAPTURE(file);
        const Built b = built(file);
        const auto report_run = NiTCAD::deck::run(b.plan);
        REQUIRE(report_run.has_value());
        REQUIRE(report_run->sweep.has_value());
        const auto& sweep = *report_run->sweep;
        CAPTURE(sweep.stopped.has_value());
        REQUIRE(!sweep.stopped.has_value());
        REQUIRE(sweep.points.size() == b.plan.points.size());
        equilibrium_checks(b);
        for (const auto& p : sweep.points) dc_checks(p, b.plan.device.mesh(), b.data, bias_label(p));
        report(file);
    }
}

TEST_CASE("TVB-32: the consistency suite on the MOS-C deck's C-V run and a zero-bias AC point", "[verify][consistency]") {
    const char* file = "examples/moscap_cv.deck";
    const Built b = built(file);
    const auto report_run = NiTCAD::deck::run(b.plan);
    REQUIRE(report_run.has_value());
    REQUIRE(report_run->small_signal.has_value());
    const auto& run = *report_run->small_signal;
    REQUIRE(!run.stopped.has_value());
    REQUIRE(run.points.size() == b.plan.points.size());
    equilibrium_checks(b);
    small_signal_checks(run, b.plan.device.mesh(), b.data, "");

    // The zero-bias AC point added through the API (TVR1-Q5), at the deck's own frequencies; the deck not edited.
    auto o = b.plan.small_signal;
    o.steady = b.plan.options;
    const std::vector<std::vector<double>> zero{std::vector<double>(b.plan.device.contacts().size(), 0.0)};
    const auto added = NiTCAD::solve::solve_small_signal(b.plan.device, zero, o);
    REQUIRE(added.has_value());
    REQUIRE(!added->stopped.has_value());
    small_signal_checks(*added, b.plan.device.mesh(), b.data, "added AC point, ");
    report(file);
}
