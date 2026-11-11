// TVR-5a benchmarks, classical interfaces and generation (DECISIONS Unit 55 17.4; criteria 17.8 and 17.14, frozen
// before code): their frozen expected outcomes. Everything is in the full tier (17.7) and runs in the Release-only
// ctest verify_benchmark.tvr5a_full.
// TVB-08 (equilibrium heterojunction): expectation written before the first run: every comparison passes (three
// junctions, each 4 D comparisons of the three equilibria, the psi order, and 5 finest-level accuracies).
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdio>
#include <string>

#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

namespace bench = NiTCAD::verify::benchmark;
using NiTCAD::verify::Outcome;

namespace {

void require_all_pass(const bench::BenchmarkResult& r, std::size_t expected_comparisons) {
    INFO(r.id);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("%s note: %s\n", r.id.c_str(), n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("%s [%s] %s: %.4g (%s) %s\n", r.id.c_str(), std::string(to_string(c.comparison.outcome)).c_str(),
                    c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == expected_comparisons);
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.q_sim, c.q_ref, c.metric, c.criterion, c.comparison.reason);
        CHECK(c.comparison.outcome == Outcome::pass);
    }
    CHECK(bench::overall(r) == Outcome::pass);
}

}  // namespace

TEST_CASE("TVB-08: the equilibrium heterojunction, full tier (frozen outcome: pass)", "[.tvr5a_full]") {
    require_all_pass(bench::run_tvb08(bench::Tier::full), 30);
}

TEST_CASE("TVB-08: nothing in the fast tier", "[tvr5a]") {
    const auto r = bench::run_tvb08(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}
