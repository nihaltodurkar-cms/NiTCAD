// TVB-31.01-31.05, the band set and narrowing (DECISIONS Unit 55 15.3, 15.9; TVR-3 step 3). Frozen expected outcomes,
// written before the first run:
//   - no run error, no fail and no not_evaluated in any entry;
//   - findings only where 15.9 and the register put them: the two grade-L table rows of 31.01 (eps_r, electron
//     affinity: no source read), in 31.04 n_i at 9.31 K on both sets (known finding (vi), MI-Q12, ruled after the
//     first run found it; the frozen criterion is unchanged), and in 31.05 every in-range comparison against the published Slotboom form (the
//     documented departure F1, known finding (F1));
//   - every other in-range comparison passes; out-of-range points, and n_i where it is not a normal double, are
//     not_applicable;
//   - each comparison's outcome is well formed (verify::validate).
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <string_view>

namespace v = NiTCAD::verify;
namespace bench = NiTCAD::verify::benchmark;
using v::Outcome;

namespace {

bool starts_with(std::string_view s, std::string_view p) { return s.substr(0, p.size()) == p; }

bool expected_finding(const bench::BenchmarkResult& r, const bench::ComparisonRecord& c) {
    if (r.id == "TVB-31.01") {
        return c.quantity == "parameter legacy_silicon.eps_r" ||
               c.quantity == "parameter legacy_silicon.electron_affinity_eV";
    }
    if (r.id == "TVB-31.04") {
        // Known finding (vi) (15.13, MI-Q12, ruled after the first run): n_i where e^(-Eg / 2kT) is subnormal.
        return (c.quantity == "legacy: n_i" || c.quantity == "default: n_i") && c.parameter > 9.0 && c.parameter < 9.7;
    }
    if (r.id == "TVB-31.05") return starts_with(c.quantity, "published form: ");
    return false;
}

struct Counts {
    int pass = 0, finding = 0, not_applicable = 0, other = 0;
};

Counts check(const bench::BenchmarkResult& r) {
    INFO(r.id);
    CHECK_FALSE(r.run_error.has_value());
    CHECK(r.notes.empty());  // in particular no "known finding ... no longer reproduces"
    Counts n;
    for (const auto& c : r.comparisons) {
        INFO(c.quantity << " at " << c.parameter << ": " << v::to_string(c.comparison.outcome) << " ("
                        << c.comparison.reason << "), sim " << c.q_sim << ", ref " << c.q_ref << ", metric "
                        << c.metric);
        CHECK(v::validate(c.comparison).empty());
        switch (c.comparison.outcome) {
            case Outcome::pass:
                CHECK_FALSE((expected_finding(r, c) && r.id != "TVB-31.04"));  // (vi): its region may pass at its edge
                ++n.pass;
                break;
            case Outcome::finding:
                CHECK(expected_finding(r, c));
                ++n.finding;
                break;
            case Outcome::not_applicable: ++n.not_applicable; break;
            default:
                ++n.other;
                FAIL_CHECK("fail or not_evaluated");
                break;
        }
    }
    return n;
}

}  // namespace

TEST_CASE("TVB-31.01: legacy band set (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_01();
    const Counts n = check(r);
    CHECK(n.finding == 2);
    // 64 points a decade over 1-500 K (173 temperatures with 300 K), six quantities each, all passing.
    CHECK(n.pass >= 6 * 173 + 5);  // + the five PB table rows
    CHECK(n.not_applicable == 6 * 2);
}

TEST_CASE("TVB-31.02: Passler gap (frozen outcomes)", "[verify][tvb31]") {
    const Counts n = check(bench::run_tvb31_02());
    CHECK(n.finding == 0);
    CHECK(n.pass >= 2 * 166 + 4);
    CHECK(n.not_applicable == 2 * 2);
}

TEST_CASE("TVB-31.03: Couderc DOS mass fits (frozen outcomes)", "[verify][tvb31]") {
    const Counts n = check(bench::run_tvb31_03());
    CHECK(n.finding == 0);
    CHECK(n.pass >= 8 * 166 + 8);
    CHECK(n.not_applicable == 8 * 2);
}

TEST_CASE("TVB-31.04: n_i, intrinsic level, d ln n_i / dT (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_04();
    const Counts n = check(r);
    CHECK(n.finding == 2);
    CHECK(n.pass >= 3 * 166);
    // n_i is not a normal double at the lowest temperatures only.
    for (const auto& c : r.comparisons) {
        if (c.comparison.outcome == Outcome::not_applicable && c.quantity.ends_with("n_i")) {
            CHECK((c.parameter < 12.0 || c.parameter > 400.0));
        }
    }
}

TEST_CASE("TVB-31.05: Slotboom narrowing (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_05();
    const Counts n = check(r);
    // 7 decades at 64 a decade (449 points), N0 and its two neighbours, N0 (1 + 10^-j) for j = 1-12; three quantities,
    // two forms.
    CHECK(n.pass >= 3 * 464 + 2);  // + the two PB table rows
    CHECK(n.finding >= 3 * 464);
    CHECK(n.not_applicable == 2 * 3 * 3);
}
