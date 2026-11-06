// TVB-31.20-31.22, impact ionization and band-to-band generation (DECISIONS Unit 55 15.3, 15.9; TVR-3 step 7).
// Frozen expected outcomes, written before the first run:
//   - no run error, fail or not_evaluated;
//   - findings only: 31.20's seven PC coefficient rows (F4) and its tag-C gamma(300 K) pin ("regression pin, not
//     evidence"); 31.21's two grade-L Kane rows (F5); and, ruled after the first run found them (MI-Q18), 31.22's
//     short segments (xi) and its path-rate dG/dI_ik as z -> 0 (xii);
//   - every other in-range comparison passes: 31.20's PB phonon-energy row, alpha and d alpha / dE on both branches,
//     the documented cut gamma B / E > 700 by 15.9's cut rule, gamma and d gamma / dT; 31.21's G and dG/dF with the cut
//     B / F > 700; 31.22's kappa, both antiderivatives, the segment integrals and partials, the path rate and its
//     four partials, eq (8) (outside (xi) and (xii));
//   - not_applicable: T 150 and 600 K (31.20); E 1e8 V/cm (31.20, 31.21); delta -0.5 and 1.5, and F 1e6 V/m where the
//     exponent exceeds 700 (31.22).
// Grid sizes: 31.20 has 30 temperatures and 225 fields per (carrier, T); 31.21 222 fields; 31.22 1057 deltas, 186
// segments and 194 fields per band, and 672 path-rate points.
#include "tvb31_check.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>

using tvb31_test::check;
using tvb31_test::starts_with;
namespace bench = NiTCAD::verify::benchmark;

TEST_CASE("TVB-31.20: van Overstraeten-de Man (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_20(bench::Tier::full);
    const auto n = check(r, [](const bench::ComparisonRecord& c) {
        return starts_with(c.quantity, "parameter ") && !c.quantity.ends_with("phonon_energy_eV");
    });
    CHECK(r.comparisons.size() == 30u * 2 + 2u * 3 * 225 * 2 + 9);
    CHECK(n.finding == 8);
    CHECK(n.not_applicable == 2 * 2 + 2 * 3 * 2);
    CHECK(n.pass == static_cast<int>(r.comparisons.size()) - 8 - n.not_applicable);
}

TEST_CASE("TVB-31.21: Kane (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_21(bench::Tier::full);
    const auto n = check(r, [](const bench::ComparisonRecord& c) { return starts_with(c.quantity, "parameter "); });
    CHECK(r.comparisons.size() == 222u * 2 + 2);
    CHECK(n.finding == 2);
    CHECK(n.not_applicable == 2);
    CHECK(n.pass == static_cast<int>(r.comparisons.size()) - 4);
}

TEST_CASE("TVB-31.22: direct-gap WKB (frozen outcomes)", "[verify][.tvb31_full]") {
    // Known findings (xi) and (xii) (15.17, MI-Q18, ruled after the first run found them): wkb_segment's short
    // segments (the antiderivative difference's rounding floor; below D = 1e-10 the midpoint rule's truncation) and
    // wkb_path_rate's dG/dI_ik as z -> 0 (its rounding floor). Each region also passes where the miss is small.
    const auto r = bench::run_tvb31_22(bench::Tier::full);
    const auto n = check(
        r,
        [](const bench::ComparisonRecord& c) {
            return starts_with(c.quantity, "path dG/dI_ik") ||
                   (c.quantity.find(" [") != std::string::npos && !starts_with(c.quantity, "path "));
        },
        true);
    CHECK(r.comparisons.size() == 2u * (1057 * 3 + 186 * 6 + 194) + 672 * 5);
    CHECK(n.finding > 0);
    CHECK(n.not_applicable == 2 * (2 * 3 + 1));
    CHECK(n.pass + n.finding == static_cast<int>(r.comparisons.size()) - n.not_applicable);
}
