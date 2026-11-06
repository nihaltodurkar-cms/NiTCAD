// TVB-31.30-31.31, the material tables (DECISIONS Unit 55 15.3, 15.9, 15.19; TVR-3 step 9). Frozen expected outcomes,
// written before the first run:
//   - no run error, fail or not_evaluated (no uncovered field, MI-Q25);
//   - 31.30: 7 rows; one pass (the PA gap 9.0), six findings (eps_r, chi_ox, kappa and rho c L; both masses C);
//   - 31.31: 86 rows; three passes (the PA radiative coefficients of InGaAs, 4H-SiC and AlGaAs), every other row a
//     finding (Ge's PA radiative row a register mismatch); 7 AlGaAs laws at 10 compositions, all findings (PC) but
//     N_c at x = 0.45, beyond its printed range x < 0.41 (not_applicable).
#include "tvb31_check.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>

using tvb31_test::check;
using tvb31_test::starts_with;
namespace bench = NiTCAD::verify::benchmark;

namespace {
bool pa_match(const bench::ComparisonRecord& c) {
    return c.quantity == "parameter silicon_dioxide_tunnelling.band_gap_eV" ||
           c.quantity == "parameter indium_gallium_arsenide.radiative_cm3_s" ||
           c.quantity == "parameter silicon_carbide_4h.radiative_cm3_s" ||
           c.quantity == "parameter algaas.radiative_cm3_s";
}
}  // namespace

TEST_CASE("TVB-31.30: SiO2 data (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_30(bench::Tier::fast);
    const auto n = check(r, [](const bench::ComparisonRecord& c) { return !pa_match(c); });
    CHECK(r.comparisons.size() == 7u);
    CHECK(n.pass == 1);
    CHECK(n.finding == 6);
}

TEST_CASE("TVB-31.31: non-silicon material sets (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_31(bench::Tier::fast);
    const auto n = check(r, [](const bench::ComparisonRecord& c) {
        return starts_with(c.quantity, "composition ") || (starts_with(c.quantity, "parameter ") && !pa_match(c));
    });
    CHECK(r.comparisons.size() == 86u + 7 * 10);
    CHECK(n.pass == 3);
    CHECK(n.not_applicable == 1);
    CHECK(n.finding == 83 + 69);
    // Added after the first run (a tightening, TVM2): a PC or L row is a finding whether or not it matches, so the
    // counts above cannot see a changed value. Pin the split of 15.19's as-built table: 21 PC rows matching the
    // printed value, 25 PC rows and Ge's PA radiative row not matching, 36 L rows; 4 law points within their digits.
    int pc_met = 0, pc_missed = 0, l_rows = 0, pa_mismatch = 0, law_met = 0;
    for (const auto& c : r.comparisons) {
        const std::string& why = c.comparison.reason;
        const bool met = why.size() >= 4 && why.compare(why.size() - 4, 4, " met") == 0;
        if (starts_with(c.quantity, "composition ")) {
            if (met) ++law_met;
        } else if (starts_with(why, "provenance insufficient (PC)")) {
            ++(met ? pc_met : pc_missed);
        } else if (starts_with(why, "provenance insufficient (L")) {
            ++l_rows;
        } else if (starts_with(why, "register mismatch")) {
            ++pa_mismatch;
        }
    }
    CHECK(pc_met == 21);
    CHECK(pc_missed == 25);
    CHECK(pa_mismatch == 1);
    CHECK(l_rows == 36);
    CHECK(law_met == 4);
}
