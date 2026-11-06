// TVR-4 device-theory benchmarks (DECISIONS Unit 55 16.4; criteria 16.8 and 16.19, frozen before code): their
// frozen expected outcomes.
// TVB-06 (dielectric stack, insulator-only device, 16.22 DT-Q19 A) is D-class: C from both admittances and the
// field ratio pass to 64 eps N^2 at every level; fast runs 1D N = 32, 128 (16.7), full adds 2D 32 x 32 and runs in
// the Release-only ctest verify_benchmark.tvr4_full with TVB-05.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

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
// only (16.7). Frozen outcome: as found (16.35), 13 pass and 28 not_evaluated.
TEST_CASE("TVB-16: the ideal diode, full tier (frozen outcome: as found, 16.37)", "[.tvr4_full]") {
    const auto r = bench::run_tvb16(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("TVB-16 note: %s\n", n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("TVB-16 [%s] %s: %.4g (%s) %s\n", std::string(to_string(c.comparison.outcome)).c_str(),
                    c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == 41);
    // Frozen as found after DT-Q32 A (DECISIONS Unit 55 16.37): these pass; the long base at 0.5 and 0.6 V are known
    // finding (x) (DT-Q27 A); every other comparison is not_evaluated (U_disc not estimated: outside the asymptotic
    // range and the finest-triple order outside 2 +- 0.5 or oscillating; or a factorial state not reached).
    // The passes' bounds as found, to 1e-6: they pin U_num's wiring (U_disc, the factorial, eps).
    const std::vector<std::string> findings{"J (ii) SRH, long base, V = 0.5 V", "J (ii) SRH, long base, V = 0.6 V",
                                            "ideality (ii) SRH, long base, V = 0.5 .. 0.6 V"};
    const std::vector<std::pair<std::string, double>> passing{
        {"J (i) no G-R, short base, V = 0.2 V", 0.012969149188851229},
        {"J (i) no G-R, short base, V = 0.3 V", 0.012900968947520106},
        {"J (i) no G-R, short base, V = 0.4 V", 0.01284897206434447},
        {"J (i) no G-R, short base, V = 0.5 V", 0.012984518235497101},
        {"J (i) no G-R, short base, V = 0.6 V", 0.02062468841401274},
        {"ideality (i) no G-R, short base, V = 0.2 .. 0.3 V", 0.015638331230131353},
        {"ideality (i) no G-R, short base, V = 0.3 .. 0.4 V", 0.015576529070217962},
        {"ideality (i) no G-R, short base, V = 0.4 .. 0.5 V", 0.015640656838250177},
        {"ideality (i) no G-R, short base, V = 0.5 .. 0.6 V", 0.02048675753687821},
        {"J (ii) SRH, short base, V = 0.5 V", 0.012982402586626456},
        {"J (ii) SRH, short base, V = 0.6 V", 0.02062463318097297},
        {"J (ii) SRH, short base, V = -1 V", 0.012000282544006202},
        {"J (ii) SRH, short base, V = -2 V", 0.011961397873277196},
        {"ideality (ii) SRH, short base, V = 0.5 .. 0.6 V", 0.020490327411534195},
        {"J (ii) SRH, long base, V = 0.4 V", 0.0003116398018498728},
        {"J (ii) SRH, long base, V = -0.5 V", 0.0008452570169339834},
        {"J (ii) SRH, long base, V = -1 V", 0.0004328477970241768},
        {"J (ii) SRH, long base, V = -2 V", 0.00022225934829510192},
        {"reverse exponent (ii) SRH, long base, V = -0.5 .. -1 V", 0.0027973850856005764}};
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.metric, c.criterion, c.comparison.reason);
        const auto it = std::find_if(passing.begin(), passing.end(), [&](const auto& p) { return p.first == c.quantity; });
        const bool pass = it != passing.end();
        const bool finding = std::find(findings.begin(), findings.end(), c.quantity) != findings.end();
        CHECK(c.comparison.outcome == (pass ? Outcome::pass : finding ? Outcome::finding : Outcome::not_evaluated));
        if (pass) {
            const double bound = std::stod(c.criterion.substr(c.criterion.rfind("= ") + 2));
            CHECK(std::abs(bound - it->second) <= 1e-6 * it->second);
        }
    }
}

TEST_CASE("TVB-16: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb16(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-18 (the diode's admittance at forward bias, 5 um bases; A, 16.19 frozen; DT-Q29-Q31 A; construction 16.34) is
// full tier only (16.7). Frozen outcome, set before its first run: every comparison passes (63: |Y|, Re Y, arg Y at
// 3 biases x 7 frequencies).
TEST_CASE("TVB-18: the diode admittance, full tier (frozen outcome: pass)", "[.tvr4_full]") {
    const auto r = bench::run_tvb18(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("TVB-18 note: %s\n", n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("TVB-18 [%s] %s: %.4g (%s) %s\n", std::string(to_string(c.comparison.outcome)).c_str(),
                    c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == 63);
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.metric, c.criterion, c.comparison.reason);
        CHECK(c.comparison.outcome == Outcome::pass);
    }
}

TEST_CASE("TVB-18: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb18(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}
