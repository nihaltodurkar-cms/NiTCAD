// TVB-32, the consistency suite on real solves of the three TVR-1 decks (DECISIONS Unit 55 section 11.6, TVR1-Q5 and
// the TVB-32 rulings of 11.13). The runs are verify::benchmark::tvb32_entries (src/verify/benchmark/tvb32.cpp, moved
// there by RC-Q1 with the decks embedded, RC-Q5); this test holds every entry to the expectation frozen in 11.13
// before any run: C1, C2, C7, C8 pass at a point whose biases are exactly 0 (one thermal temperature), else NOT
// APPLICABLE naming "all biases 0"; C3, C4, C12 pass; C5 NOT APPLICABLE ("electrothermal on"); C1 and C3 pass on the
// equilibrium; with two contacts a zero-bias C7 or C8 carries its "two contacts" diagnostic. Any other outcome, and
// any failed or NOT EVALUATED part, is a finding. Per deck and check the worst residual/tolerance ratio is printed.
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/benchmark/tvb32.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdio>
#include <map>
#include <string>

namespace c = NiTCAD::verify::consistency;
namespace bench = NiTCAD::verify::benchmark;
using NiTCAD::verify::Outcome;

namespace {

const char* name_of(Outcome o) {
    switch (o) {
        case Outcome::pass: return "pass";
        case Outcome::fail: return "FAIL";
        case Outcome::not_applicable: return "not applicable";
        case Outcome::not_evaluated: return "NOT EVALUATED";
        default: return "other";
    }
}

struct Tally {
    std::map<std::string, int> outcomes;
    double worst = 0.0;
    std::string where;
};

}  // namespace

TEST_CASE("TVB-32: the consistency suite on the three decks (frozen outcomes)", "[verify][consistency]") {
    const bench::Tvb32Run run = bench::tvb32_entries();
    CAPTURE(run.run_error.value_or(""));
    REQUIRE_FALSE(run.run_error.has_value());
    std::map<std::string, Tally> tally;  // "deck check"
    std::map<std::string, int> per_deck;
    for (const bench::Tvb32Entry& e : run.entries) {
        const c::CheckResult& r = e.result;
        ++per_deck[e.deck];
        const std::string& id = r.id;
        Outcome want = Outcome::pass;
        const char* predicate = "";
        if (e.at != "equilibrium" && (id == "C1" || id == "C2" || id == "C7" || id == "C8")) {
            want = e.zero_bias ? Outcome::pass : Outcome::not_applicable;
            predicate = "all biases 0";
        } else if (id == "C5") {
            want = Outcome::not_applicable;
            predicate = "electrothermal on";
        }
        Tally& t = tally[e.deck + " " + id];
        ++t.outcomes[name_of(r.comparison.outcome)];
        for (const auto& p : r.parts) {
            if (p.outcome != Outcome::pass || !(p.tolerance > 0.0)) continue;
            const double q = p.residual / p.tolerance;
            if (q > t.worst) {
                t.worst = q;
                t.where = e.at + ", " + p.name + " at " + p.where;
            }
        }
        CAPTURE(e.deck, e.at, id, r.comparison.reason);
        for (const auto& p : r.parts) {
            CAPTURE(p.name, name_of(p.outcome), p.residual, p.tolerance, p.where, p.reason);
            if (p.outcome == Outcome::fail || p.outcome == Outcome::not_evaluated) FAIL_CHECK("part " << p.name);
        }
        CHECK(r.comparison.outcome == want);
        if (want == Outcome::not_applicable) CHECK(r.comparison.reason.find(predicate) != std::string::npos);
        if ((id == "C7" || id == "C8") && want == Outcome::pass && e.contacts == 2) {
            CHECK(std::any_of(r.diagnostics.begin(), r.diagnostics.end(),
                              [](const c::Diagnostic& g) { return g.name == "two contacts"; }));
        }
    }
    // pn: 15 points x 6 + 2 equilibrium; Schottky: 29 x 6 + 2; MOS-C: 20 AC points (the C-V run and the added one)
    // x 8 + 2 (the counts of the as-built table before the move).
    CHECK(per_deck["pn_diode_1d"] == 92);
    CHECK(per_deck["schottky_diode_1d"] == 176);
    CHECK(per_deck["moscap_cv"] == 162);
    for (const auto& [key, t] : tally) {
        std::string counts;
        for (const auto& [o, n] : t.outcomes) counts += (counts.empty() ? "" : ", ") + std::to_string(n) + " " + o;
        std::printf("TVB-32 %s: %s; worst %.3g of the bound (%s)\n", key.c_str(), counts.c_str(), t.worst,
                    t.where.c_str());
    }
}
