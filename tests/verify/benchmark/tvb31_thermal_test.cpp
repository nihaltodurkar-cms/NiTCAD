// TVB-31.26-31.29, emission, MLDA, conductivity and the thermopower factor (DECISIONS Unit 55 15.3, 15.9; TVR-3 step
// 8). Frozen expected outcomes, written before the first run:
//   - no run error or fail;
//   - findings only: 31.27's two PC quantization masses; 31.28's three register rows (1.48 PC, the exponent's
//     mismatch 1.33 against the printed 1.3, the heat capacity L); and 31.29's dh/dx for x < -1 within its rounding
//     floor, (xiii) (MI-Q21, ruled before the first run from the shipped form);
//   - not_evaluated only in 31.29, at F_-1, F_0, F_1 and dh/dx where |q_ref| < 2^-969 (MI-Q19: U_ref floored at the
//     subnormal spacing exceeds 1e-3 of the bound);
//   - every other in-range comparison passes;
//   - not_applicable: 31.26 T 150 and 600 K, and on the Couderc set T above 400 K; 31.27 z / lambda above 10 and a /
//     lambda 1e-9 and 1e4; 31.28 T 150 and 700 K; 31.29 x = +-760.
// Grid sizes: 31.26 has 3 x 385 (N, T) points and 30 temperatures per (set, form, carrier), the legacy
// density-of-states form adding the Schottky overload; 31.27 56 wavelengths, 769 z and 710 a; 31.28 35 temperatures;
// 31.29 10104 x.
#include "tvb31_check.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using tvb31_test::check;
using tvb31_test::never;
using tvb31_test::starts_with;
namespace bench = NiTCAD::verify::benchmark;

TEST_CASE("TVB-31.26: emission velocity (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_26(bench::Tier::full);
    const auto n = check(r, never);
    CHECK(r.comparisons.size() == 3u * 385 + 2 * (90 + 60) + 2 * (60 + 60));
    CHECK(n.not_applicable == 2 * (6 + 4) + 2 * (18 + 18));
    CHECK(n.pass == static_cast<int>(r.comparisons.size()) - n.not_applicable);
}

TEST_CASE("TVB-31.27: MLDA (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_27(bench::Tier::full);
    const auto n = check(r, [](const bench::ComparisonRecord& c) { return starts_with(c.quantity, "parameter "); });
    CHECK(r.comparisons.size() == 56u + 769 + 710 + 2);
    CHECK(n.finding == 2);
    CHECK(n.not_applicable == 128 + 2);
    CHECK(n.pass == static_cast<int>(r.comparisons.size()) - 2 - n.not_applicable);
}

TEST_CASE("TVB-31.28: lattice conductivity (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_28(bench::Tier::full);
    const auto n = check(r, [](const bench::ComparisonRecord& c) { return starts_with(c.quantity, "parameter "); });
    CHECK(r.comparisons.size() == 35u * 2 + 3);
    CHECK(n.finding == 3);
    CHECK(n.not_applicable == 4);
    CHECK(n.pass == static_cast<int>(r.comparisons.size()) - 3 - 4);
}

TEST_CASE("TVB-31.29: Fermi orders and h (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_29(bench::Tier::full);
    const auto n = check(
        r, [](const bench::ComparisonRecord& c) { return c.quantity == "dh/dx" && c.parameter < -1.0; }, true,
        [](const bench::ComparisonRecord& c) {
            return (c.quantity == "F_-1" || c.quantity == "F_0" || c.quantity == "F_1" || c.quantity == "dh/dx") &&
                   c.q_ref != 0.0 && std::abs(c.q_ref) < 0x1p-969;
        });
    CHECK(r.comparisons.size() == 10104u * 5);
    CHECK(n.finding > 0);
    CHECK(n.not_evaluated > 0);
    CHECK(n.not_applicable == 2 * 5);
    CHECK(n.pass + n.finding + n.not_evaluated == static_cast<int>(r.comparisons.size()) - n.not_applicable);
}
