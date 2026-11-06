// TVB-31.23-31.25, the WKB action and the Tsu-Esaki supply and current (DECISIONS Unit 55 15.3, 15.9; TVR-3 step 8).
// Frozen expected outcomes, written before the first run:
//   - no run error, fail or not_evaluated in 31.23 and 31.24 (no finding in 31.23): every comparison passes
//     (31.23's segment action and partials, the scale, S, dS/dE and TC of the four barriers; 31.24's N and its three
//     partials) but
//     31.24's points with a level more than 60 kT away (not_applicable); and, ruled after the first run found it
//     (MI-Q22), 31.24's dN/dE where both levels are many kT above E (xiv);
//   - 31.25 has no comparisons: 33 diagnostic notes (32 points and the worst), each "diagnostic (not judged, MI-Q8)";
//     8 more of a kinked barrier were added after the first run (no criterion: diagnostics only).
// Grid sizes: 31.23 has 304 (a, b) pairs, the scale, and 295, 295, 295 and 301 energies; 31.24 5 x 1025 points.
#include "tvb31_check.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

#include "NiTCAD/base/constants.hpp"

using tvb31_test::check;
using tvb31_test::never;
using tvb31_test::starts_with;
namespace bench = NiTCAD::verify::benchmark;

TEST_CASE("TVB-31.23: WKB action (frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_23(bench::Tier::full);
    const auto n = check(r, never);
    CHECK(r.comparisons.size() == 304u * 3 + 1 + (3 * 295 + 301) * 3);
    CHECK(n.pass == static_cast<int>(r.comparisons.size()));
}

TEST_CASE("TVB-31.24: Tsu-Esaki supply (frozen outcomes)", "[verify][.tvb31_full]") {
    // Known finding (xiv) (15.18, MI-Q22, ruled after the first run found it): dN/dE within the rounding floor of
    // sigma(x2) - sigma(x1) where both levels are many kT above E; it passes elsewhere.
    const auto r = bench::run_tvb31_24(bench::Tier::full);
    const auto n =
        check(r, [](const bench::ComparisonRecord& c) { return starts_with(c.quantity, "dN/dE"); }, true);
    CHECK(r.comparisons.size() == 5u * 1025 * 4);
    // The not_applicable points, by the grid's own definition.
    const double kT = NiTCAD::base::k_B_eV_per_K * 300.0, F1 = 0.2;
    int na = 0;
    for (const double dF : {0.0, 1e-9, 1e-3, 0.1, 1.0}) {
        for (int i = -512; i <= 512; ++i) {
            const double E = F1 - (i / 8.0) * kT;
            if (std::abs((F1 - E) / kT) > 60.0 || std::abs((F1 - dF - E) / kT) > 60.0) na += 4;
        }
    }
    CHECK(n.not_applicable == na);
    CHECK(n.finding > 0);
    CHECK(n.pass + n.finding == static_cast<int>(r.comparisons.size()) - na);
}

TEST_CASE("TVB-31.25: Tsu-Esaki current (diagnostic, frozen outcomes)", "[verify][.tvb31_full]") {
    const auto r = bench::run_tvb31_25(bench::Tier::full);
    CHECK_FALSE(r.run_error.has_value());
    CHECK(r.comparisons.empty());
    REQUIRE(r.notes.size() == 41u);  // 33 as frozen, and 8 of the kinked barrier added after the first run
    for (const std::string& note : r.notes) {
        INFO(note);
        CHECK(starts_with(note, "diagnostic (not judged, MI-Q8): "));
    }
}
