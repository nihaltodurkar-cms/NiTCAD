// TVR-4 device-theory benchmarks (DECISIONS Unit 55 16.4; criteria 16.8 and 16.19, frozen before code): their
// frozen expected outcomes.
// TVB-06 (dielectric stack, insulator-only device, 16.22 DT-Q19 A) is D-class: C from both admittances and the
// field ratio pass to 64 eps N^2 at every level; fast runs 1D N = 32, 128 (16.7), full adds 2D 32 x 32 and runs in
// the Release-only ctest verify_benchmark.tvr4_full with TVB-05.
#include <catch2/catch_test_macros.hpp>

#include <cstdio>

#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

namespace bench = NiTCAD::verify::benchmark;
using NiTCAD::verify::Outcome;

namespace {

void require_all_pass(const bench::BenchmarkResult& r, std::size_t expected_comparisons) {
    INFO(r.id);
    if (r.run_error) FAIL(*r.run_error);
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

TEST_CASE("TVB-06: the dielectric stack, fast tier (frozen outcome: pass)", "[tvr4]") {
    require_all_pass(bench::run_tvb06(bench::Tier::fast), 6);  // C twice and the field ratio, N = 32 and 128
}

TEST_CASE("TVB-06: the dielectric stack, full tier (frozen outcome: pass)", "[.tvr4_full]") {
    require_all_pass(bench::run_tvb06(bench::Tier::full), 9);  // + 2D 32 x 32
}

// TVB-05 (coaxial capacitor, 2D generated, M = 32 ... 256 sectors) is E-class and full tier only (16.7): its order
// and finest accuracy pass. Release-only ctest verify_benchmark.tvr4_full (hidden tag), as tvb31_full.
TEST_CASE("TVB-05: the coaxial capacitor, full tier (frozen outcome: pass)", "[.tvr4_full]") {
    require_all_pass(bench::run_tvb05(bench::Tier::full), 2);  // the order of C and C on the finest level
}

TEST_CASE("TVB-05: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb05(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}
