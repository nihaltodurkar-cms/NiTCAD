// TVB-31.17-31.19, incomplete ionization and the neutral equilibria (DECISIONS Unit 55 15.3, 15.9; TVR-3 step 6).
// Frozen expected outcomes, written before the first run:
//   - no run error, fail or not_evaluated;
//   - findings only: 31.17's four PC table rows (the legacy 45 meV levels and degeneracies 2 and 4); and, ruled
//     after the first run found it (MI-Q16), 31.19's Fermi-Dirac complete-ionization eta where |eta| is small (x);
//   - every other in-range comparison passes (31.18's ten PA Table I rows bit-equal to the printed decimals);
//   - not_applicable: eta +-50 and depth 0.35 eV (31.17), N 1e22 (31.18), net doping +-1e22 (31.19).
// Grid sizes: 31.17 has 2905 points per (g, b); 31.18 451 per dopant; 31.19 2693 complete-ionization roots per
// statistics and 1435 ionized roots per (statistics, levels).
#include "tvb31_check.hpp"

#include <catch2/catch_test_macros.hpp>

using tvb31_test::check;
using tvb31_test::starts_with;
namespace bench = NiTCAD::verify::benchmark;

TEST_CASE("TVB-31.17: incomplete ionization (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_17(bench::Tier::full);
    const auto n = check(r, [](const bench::ComparisonRecord& c) { return starts_with(c.quantity, "parameter "); });
    CHECK(r.comparisons.size() == 6u * 2905 * 2 + 4);
    CHECK(n.finding == 4);
    CHECK(n.not_applicable == 6 * 15 * 2);
    CHECK(n.pass == static_cast<int>(r.comparisons.size()) - 4 - n.not_applicable);
}

TEST_CASE("TVB-31.18: Altermatt-Schenk levels (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_18(bench::Tier::full);
    const auto n = check(r, [](const bench::ComparisonRecord&) { return false; });
    CHECK(r.comparisons.size() == 2u * 451 * 2 + 10);
    CHECK(n.not_applicable == 4);
    CHECK(n.pass == static_cast<int>(r.comparisons.size()) - 4);
}

TEST_CASE("TVB-31.19: neutral equilibria (frozen outcomes)", "[verify][.tvb31_full]") {
    // Known finding (x) (15.16, MI-Q16, ruled after the first run found it): the Fermi-Dirac complete-ionization eta
    // within the residual's rounding floor, where |eta| is small; it passes elsewhere.
    const auto r = bench::run_tvb31_19(bench::Tier::full);
    const auto n = check(
        r, [](const bench::ComparisonRecord& c) { return starts_with(c.quantity, "eta (Fermi-Dirac, complete)"); },
        true);
    CHECK(r.comparisons.size() == 3u * (2 * 2693 + 4 * 1435));
    CHECK(n.finding > 0);
    CHECK(n.not_applicable == 2 * 2 * 3);
    CHECK(n.pass + n.finding == static_cast<int>(r.comparisons.size()) - n.not_applicable);
}
