// TVR-4 device-theory benchmarks (DECISIONS Unit 55 16.4; criteria 16.8 and 16.19, frozen before code): their
// frozen expected outcomes.
// TVB-06 (dielectric stack, insulator-only device, 16.22 DT-Q19 A) is D-class: C from both admittances and the
// field ratio pass to 64 eps N^2 at every level; fast runs 1D N = 32, 128 (16.7), full adds 2D 32 x 32 and runs in
// the Release-only ctest verify_benchmark.tvr4_full with TVB-05.
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
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

// TVB-05 (coaxial capacitor, 2D generated, M = 32 ... 256 sectors) is E-class and full tier only (16.7). Release-only
// ctest verify_benchmark.tvr4_full (hidden tag), as tvb31_full. Recorded as found (16.24, DT-Q20 A): the finest
// accuracy passes; the order is not_evaluated, its successive orders 2.0545, 2.0825, 2.1846 spreading 0.102 > 0.1
// (the M = 32 pair), with p-hat 2.054 in the band. Neither the criteria nor the ladder changed (TVM2).
// The (h, tol) factorial (16.2 (p), DT-Q23 A) needs M = 512, about 4 minutes a solve, so it runs only in the nightly
// tier (verify_full.tvb05, local; 16.30). In full the finest accuracy is therefore held at not_evaluated (4.3).
TEST_CASE("TVB-05: the coaxial capacitor, full tier (frozen outcome: order and accuracy not_evaluated)",
          "[.tvr4_full]") {
    const auto r = bench::run_tvb05(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    REQUIRE(r.comparisons.size() == 2);
    const auto& order = r.comparisons[0];
    const auto& finest = r.comparisons[1];
    CAPTURE(order.metric, order.comparison.reason, finest.metric);
    std::printf("TVB-05 %s: %.4f (%s)\n", order.quantity.c_str(), order.metric, order.comparison.reason.c_str());
    std::printf("TVB-05 %s: %.3e (%s)\n", finest.quantity.c_str(), finest.metric, finest.criterion.c_str());
    CHECK(order.quantity == "C order");
    CHECK(order.comparison.outcome == Outcome::not_evaluated);
    CHECK(order.comparison.reason == "not in the asymptotic range");
    CHECK(order.metric >= 1.8);  // p-hat itself is in the band
    CHECK(order.metric <= 2.3);
    CHECK(finest.metric <= 1e-3);  // the accuracy itself (3.7e-5), not judged without the factorial
    CHECK(finest.comparison.outcome == Outcome::not_evaluated);
    CHECK(finest.comparison.reason == "the (h, tol) factorial was not evaluated");
    CHECK(bench::overall(r) == Outcome::not_evaluated);
}

TEST_CASE("TVB-05: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb05(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-15 (minority diffusion with SRH, the TVQ3 hook's G) is E (linear model, low injection) + A (injection ladder),
// full tier only (16.7): the order and finest accuracy of Delta n pass, and the departure's order in G and its
// G -> 0 limit pass (16.8, 16.19).
TEST_CASE("TVB-15: minority diffusion, full tier (frozen outcome: pass)", "[.tvr4_full]") {
    require_all_pass(bench::run_tvb15(bench::Tier::full), 4);
}

TEST_CASE("TVB-15: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb15(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-17 (depletion capacitance against the QE model, DT-Q4's predicate) is E, full tier only (16.7). Recorded as
// found (16.28, DT-Q22 A): at each of V = -0.5, -1, -2 V the finest accuracy of C_j passes, and the order is
// not_evaluated (spread 0.13-0.23 from the N = 64 pair, h about the Debye length) with p-hat in the band. Neither
// criteria, ladder nor construction changed (TVM2).
TEST_CASE("TVB-17: depletion capacitance, full tier (frozen outcome: orders not_evaluated, accuracy pass)",
          "[.tvr4_full]") {
    const auto r = bench::run_tvb17(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    REQUIRE(r.comparisons.size() == 6);
    for (std::size_t k = 0; k < 6; k += 2) {
        const auto& order = r.comparisons[k];
        const auto& finest = r.comparisons[k + 1];
        CAPTURE(order.quantity, order.metric, order.comparison.reason, finest.quantity, finest.metric);
        std::printf("%s: %.4f (%s)\n", order.quantity.c_str(), order.metric, order.comparison.reason.c_str());
        std::printf("%s: %.3e (%s)\n", finest.quantity.c_str(), finest.metric, finest.criterion.c_str());
        CHECK(order.comparison.outcome == Outcome::not_evaluated);
        CHECK(order.comparison.reason == "not in the asymptotic range");
        CHECK(order.metric >= 1.8);
        CHECK(order.metric <= 2.3);
        CHECK(finest.comparison.outcome == Outcome::pass);
    }
    CHECK(bench::overall(r) == Outcome::not_evaluated);
}

TEST_CASE("TVB-17: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb17(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-16 (ideal diode, three configurations; A, 16.19 frozen; DT-Q24 A, DT-Q25 A; construction 16.31) is full tier
// only (16.7). Frozen outcome: every J and derivative comparison passes.
TEST_CASE("TVB-16: the ideal diode, full tier (frozen outcome: pass)", "[.tvr4_full]") {
    const auto r = bench::run_tvb16(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("TVB-16 note: %s\n", n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("TVB-16 [%s] %s: %.4g (%s) %s\n", std::string(to_string(c.comparison.outcome)).c_str(),
                    c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == 41);
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.metric, c.criterion, c.comparison.reason);
        CHECK(c.comparison.outcome == Outcome::pass);
    }
}

TEST_CASE("TVB-16: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb16(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}
