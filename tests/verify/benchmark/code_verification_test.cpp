// TVR-2 code-verification benchmarks (DECISIONS Unit 55 14.4, criteria 14.7, frozen before code): their frozen
// expected outcomes. TVB-11 (uniform resistor) and TVB-12 (Scharfetter-Gummel at constant field) are D-class:
// every comparison passes, to 64 eps N^2.
#include <catch2/catch_test_macros.hpp>

#include <cstdio>

#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

namespace bench = NiTCAD::verify::benchmark;
using NiTCAD::verify::Outcome;

namespace {

void require_all_pass(const bench::BenchmarkResult& r, std::size_t expected_comparisons) {
    INFO(r.id);
    REQUIRE_FALSE(r.run_error.has_value());
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    REQUIRE(r.comparisons.size() == expected_comparisons);
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.q_sim, c.q_ref, c.metric, c.criterion, c.comparison.reason);
        std::printf("%s %s: %.3e (%s)\n", r.id.c_str(), c.quantity.c_str(), c.metric, c.criterion.c_str());
        CHECK(c.comparison.outcome == Outcome::pass);
    }
    CHECK(bench::overall(r) == Outcome::pass);
}

}  // namespace

TEST_CASE("TVB-11: the uniform resistor, to rounding (frozen outcome: pass)", "[tvr2]") {
    require_all_pass(bench::run_tvb11(), 6);  // current and psi at 1D N = 32, 128 and 2D 32 x 32
}

TEST_CASE("TVB-12: Scharfetter-Gummel exact at constant field (frozen outcome: pass)", "[tvr2]") {
    require_all_pass(bench::run_tvb12(), 8);  // n and p at Peclet 0.01, 1, 10, 40
}
