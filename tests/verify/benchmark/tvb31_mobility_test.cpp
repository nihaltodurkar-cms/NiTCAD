// TVB-31.06-31.10, mobility (DECISIONS Unit 55 15.3, 15.9; TVR-3 step 4). Frozen expected outcomes, written before
// the first run:
//   - no run error, fail or not_evaluated;
//   - findings only: 31.06's eight grade-L table rows (F3); 31.07's published-form slopes (known finding (F6)) and
//     its two tag-C T_exponent rows (pinned: "regression pin, not evidence"); 31.08's two Canali beta rows, outside
//     the law's printed digits (F2); 31.10's tag-C regularization row, and its mu0 partials for a >= 20 (known
//     finding (vii), MI-Q14, ruled after the first run found it);
//   - every other in-range comparison passes (31.08's two v_sat rows within the laws' printed digits).
#include "tvb31_check.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>

using tvb31_test::check;
using tvb31_test::starts_with;
namespace bench = NiTCAD::verify::benchmark;

TEST_CASE("TVB-31.06: Caughey-Thomas (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto n = check(bench::run_tvb31_06(bench::Tier::full), [](const bench::ComparisonRecord& c) {
        return starts_with(c.quantity, "parameter ") && !c.quantity.ends_with("T_exponent");
    });
    CHECK(n.finding == 8);
    CHECK(n.pass >= 2 * 2 * (514 + 4 * 26) + 2);
    CHECK(n.not_applicable == 2 * 2 * 3);
}

TEST_CASE("TVB-31.07: Masetti (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto n = check(bench::run_tvb31_07(bench::Tier::full), [](const bench::ComparisonRecord& c) {
        return starts_with(c.quantity, "published form: ") || c.quantity.ends_with("T_exponent");
    });
    CHECK(n.finding == 2 * 514 + 2);
    CHECK(n.pass >= 2 * (2 * 514 + 2 * 2 * 26) + 17);
    CHECK(n.not_applicable == 2 * (3 + 2 * 2 * 2));
}

TEST_CASE("TVB-31.08: Canali (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto n = check(bench::run_tvb31_08(bench::Tier::full), [](const bench::ComparisonRecord& c) {
        return c.quantity.ends_with("_saturation.beta");
    });
    CHECK(n.finding == 2);
    CHECK(n.pass >= 2 * 2 * 2 * 514 + 2);
    CHECK(n.not_applicable == 2 * 2 * 2);
}

TEST_CASE("TVB-31.09: Lombardi (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto n = check(bench::run_tvb31_09(bench::Tier::full), [](const bench::ComparisonRecord&) { return false; });
    CHECK(n.finding == 0);
    CHECK(n.pass >= 2 * 4 * (2 * 449 + 26) + 8);
    CHECK(n.not_applicable == 2 * 4 * 2 * 2);
}

TEST_CASE("TVB-31.10: energy-balance mobility (frozen outcomes)", "[verify][.tvb31_full]") {
    // Known finding (vii) (15.14, MI-Q14, ruled after the first run found it): the mu0 partials for a >= 20; the region
    // may pass where the cancellation happens to stay within the bound.
    const auto n = check(
        bench::run_tvb31_10(bench::Tier::full),
        [](const bench::ComparisonRecord& c) {
            return c.quantity == "parameter energy_balance_regularization" ||
                   ((c.quantity.ends_with(" dmu/dmu0") || c.quantity.ends_with(" dE_eff/dmu0")) && c.parameter >= 20.0);
        },
        true);
    CHECK(n.finding >= 1 + 6 * 100);
    CHECK(n.pass >= 3 * 6 * (11 + 1153) + 2 - n.finding);
    CHECK(n.not_applicable == 3 * 6 * 2);
}
