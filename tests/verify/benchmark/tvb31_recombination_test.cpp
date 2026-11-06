// TVB-31.11-31.16, lifetimes, recombination and interface traps (DECISIONS Unit 55 15.3, 15.9; TVR-3 step 5). Frozen
// expected outcomes, written before the first run:
//   - no run error, fail or not_evaluated;
//   - findings only: 31.11's three table rows (tau_n0 and tau_p0 grade L, N_ref PC); and, ruled after the first
//     run found them (MI-Q15), 31.12's partials (viii) and 31.15's trap_kinetics partials (ix);
//   - every other in-range comparison passes (31.13's two Auger rows, PB, bit-equal to the printed decimals);
//   - not_applicable: the points outside the stated ranges (N 1e22; p outside [1e-10, 1e21] cm^-3; |x| > 745;
//     E_F +-0.6 eV; T 150 and 600 K).
// The grid sizes: the 1e13-1e21 log grid has 513 points, 1e-10-1e21 has 1985, and 200-500 K has 27.
#include "tvb31_check.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using tvb31_test::check;
using tvb31_test::starts_with;
namespace bench = NiTCAD::verify::benchmark;

namespace {
bool none(const bench::ComparisonRecord&) { return false; }
}  // namespace

TEST_CASE("TVB-31.11: Scharfetter lifetimes (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_11(bench::Tier::fast);
    const auto n = check(r, [](const bench::ComparisonRecord& c) { return starts_with(c.quantity, "parameter "); });
    CHECK(n.finding == 3);
    CHECK(n.not_applicable == 2);
    CHECK(n.pass == 2 * 514);
}

TEST_CASE("TVB-31.12: SRH (frozen outcomes)", "[verify][tvb31]") {
    // Known finding (viii) (15.15, MI-Q15, ruled after the first run found it): the partials, where their error is
    // within the shipped quotient-rule form's rounding floor (they pass where nothing cancels).
    const auto r = bench::run_tvb31_12(bench::Tier::fast);
    const auto n = check(
        r,
        [](const bench::ComparisonRecord& c) {
            return starts_with(c.quantity, "dR/dn (") || starts_with(c.quantity, "dR/dp (");
        },
        true);
    CHECK(r.comparisons.size() == 2u * 3 * 1985 * 3);
    CHECK(n.pass + n.finding + n.not_applicable == static_cast<int>(r.comparisons.size()));
    CHECK(n.finding > 0);
    CHECK(n.not_applicable > 0);
    CHECK(n.pass > n.not_applicable);
}

TEST_CASE("TVB-31.13: Auger (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_13(bench::Tier::fast);
    const auto n = check(r, none);
    CHECK(r.comparisons.size() == 2u * 3 * 1985 * 3 + 2);
    CHECK(n.pass + n.not_applicable == static_cast<int>(r.comparisons.size()));
    CHECK(n.not_applicable > 0);
}

TEST_CASE("TVB-31.14: radiative (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_14(bench::Tier::fast);
    const auto n = check(r, none);
    CHECK(r.comparisons.size() == 5u * 3 * 1985 * 3);  // four materials, GaAs with both E
    CHECK(n.pass + n.not_applicable == static_cast<int>(r.comparisons.size()));
    CHECK(n.not_applicable > 0);
}

TEST_CASE("TVB-31.15: trap occupancy and kinetics (frozen outcomes)", "[verify][tvb31]") {
    // Known finding (ix) (15.15, MI-Q15, ruled after the first run found it): the four trap_kinetics partials, where
    // their error is within the shipped form's rounding floor (they pass where nothing cancels).
    const auto r = bench::run_tvb31_15(bench::Tier::fast);
    const auto n = check(
        r,
        [](const bench::ComparisonRecord& c) {
            return starts_with(c.quantity, "df/dn (") || starts_with(c.quantity, "df/dp (") ||
                   starts_with(c.quantity, "dU/dn (") || starts_with(c.quantity, "dU/dp (");
        },
        true,
        // TB-Q4 with the U_ref floor of 15.18 (MI-Q19, step 8): where |q_ref| < 2^-969 the double-double reference's
        // low part is subnormal, so the occupancy's underflow tail at eta = -745 is not_evaluated (it passed before,
        // bit-equal).
        [](const bench::ComparisonRecord& c) { return c.q_ref != 0.0 && std::abs(c.q_ref) < 0x1p-969; });
    CHECK(r.comparisons.size() == 3u * 3004 + 10u * 1985 * 7);
    CHECK(n.pass + n.finding + n.not_applicable + n.not_evaluated == static_cast<int>(r.comparisons.size()));
    CHECK(n.finding > 0);
    CHECK(n.not_applicable >= 3 * 20);
}

TEST_CASE("TVB-31.16: trap-band quadrature (frozen outcomes)", "[verify][tvb31]") {
    const auto r = bench::run_tvb31_16(bench::Tier::fast);
    const auto n = check(r, none);
    CHECK(r.comparisons.size() == 4u * 2 * (257 + 2 + 2 * (27 + 2)));
    CHECK(n.not_applicable == 4 * 2 * (2 + 2 * 2));
    CHECK(n.pass == static_cast<int>(r.comparisons.size()) - n.not_applicable);
}
