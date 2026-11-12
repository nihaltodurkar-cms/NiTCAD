// TVR-5a benchmarks, classical interfaces and generation (DECISIONS Unit 55 17.4; criteria 17.8 and 17.14, frozen
// before code): their frozen outcomes. Everything is in the full tier (17.7) and runs in the Release-only ctest
// verify_benchmark.tvr5a_full.
//
// TVB-08 (equilibrium heterojunction). The expectation written before the first run was that all 30 comparisons
// pass. It was not met. Recorded as found (17.17, IG-Q19 A), with no criterion, ladder or mesh changed (TVM2):
//   - the three equilibria agree at every level of every junction (12 pass);
//   - p-N: the order and all five finest-level accuracies pass;
//   - n-N: the order is not_evaluated (p-hat 1.98, spread 0.109 above 0.1); psi, the interface level and the left
//     band bending miss at N = 512 (2.5e-4, 2.3e-4 against 1e-4; 1.16e-3 against 1e-3); D and the right bending pass;
//   - chi-eps: the order passes; psi and the interface level miss (2.3e-4, 2.1e-4); D and both bendings pass.
// The misses are second-order discretization error on the frozen ladder (N = 1024, run by the factorial, gives
// 6.2e-5 and 5.8e-5), not a defect found in the solver.
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

namespace bench = NiTCAD::verify::benchmark;
using NiTCAD::verify::Outcome;

namespace {

// A frozen outcome recorded as found: per comparison, in order, its quantity, outcome and a window on its metric
// (0: not pinned).
struct Found {
    std::string quantity;
    Outcome outcome;
    double metric_above = 0.0, metric_below = 0.0;
};

void require_as_found(const bench::BenchmarkResult& r, const std::vector<Found>& expected) {
    INFO(r.id);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("%s note: %s\n", r.id.c_str(), n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("%s [%s] %s: %.4g (%s) %s\n", r.id.c_str(), std::string(to_string(c.comparison.outcome)).c_str(),
                    c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const auto& c = r.comparisons[i];
        CAPTURE(c.quantity, c.metric, c.criterion, c.comparison.reason);
        CHECK(c.quantity == expected[i].quantity);
        CHECK(c.comparison.outcome == expected[i].outcome);
        if (expected[i].metric_above > 0.0) CHECK(c.metric > expected[i].metric_above);
        if (expected[i].metric_below > 0.0) CHECK(c.metric < expected[i].metric_below);
    }
}

}  // namespace

TEST_CASE("TVB-08: the equilibrium heterojunction, full tier (frozen outcome: as found, 17.17)", "[.tvr5a_full]") {
    const auto pass = Outcome::pass, fail = Outcome::fail, ne = Outcome::not_evaluated;
    std::vector<Found> expected;
    // Per junction: the order's outcome and window, then psi, interface level, D, left and right bending.
    struct Junction {
        const char* name;
        Found order, psi, level, field, left, right;
    };
    const Junction junctions[] = {
        {"n-N", {"", ne, 1.95, 2.0}, {"", fail, 2.2e-4, 2.7e-4}, {"", fail, 2.0e-4, 2.5e-4}, {"", pass, 1.4e-4, 1.8e-4},
         {"", fail, 1.05e-3, 1.3e-3}, {"", pass, 2.5e-4, 3.2e-4}},
        {"p-N", {"", pass, 1.99, 2.01}, {"", pass, 6.0e-6, 7.4e-6}, {"", pass, 1.2e-8, 1.6e-8},
         {"", pass, 3.6e-6, 4.5e-6}, {"", pass, 2.5e-8, 3.3e-8}, {"", pass, 2.4e-8, 3.1e-8}},
        {"chi-eps", {"", pass, 1.96, 2.0}, {"", fail, 2.1e-4, 2.5e-4}, {"", fail, 1.85e-4, 2.3e-4},
         {"", pass, 1.4e-4, 1.8e-4}, {"", pass, 2.4e-4, 3.0e-4}, {"", pass, 8.2e-4, 9.9e-4}},
    };
    for (const Junction& j : junctions) {
        const std::string at = j.name;
        for (const int N : {64, 128, 256, 512}) {
            // The agreement is bounded by the Newton term; its ratio to the bound stays below 0.1.
            expected.push_back({"equilibria agree, " + at + ", N = " + std::to_string(N), pass, 0.0, 0.1});
        }
        const auto add = [&](const std::string& q, const Found& f) {
            expected.push_back({q, f.outcome, f.metric_above, f.metric_below});
        };
        add("psi Linf order, " + at, j.order);
        add("psi Linf finest (N = 512), " + at, j.psi);
        add("interface level finest (N = 512), " + at, j.level);
        add("fields (D), both sides finest (N = 512), " + at, j.field);
        add("band bending, left finest (N = 512), " + at, j.left);
        add("band bending, right finest (N = 512), " + at, j.right);
    }
    const auto r = bench::run_tvb08(bench::Tier::full);
    require_as_found(r, expected);
    CHECK(bench::overall(r) == fail);
}

// TVB-22 (interface trap level), as rebuilt once under IG-Q22 B (17.19: N_t 1e11, sigma_p 1e-18, biases V_FB - 0.20
// and - 0.25 V). The expectation written before the rerun was that all 40 comparisons pass (two trap types x two
// biases, each two orders, Y at seven frequencies and the peak's place). Recorded as found (17.20):
//   - every comparison is applicable (P1 <= 8e-9, P2 <= 5.6e-6, P3 <= 4e-10);
//   - 30 pass: Y at w tau >= 0.3, the G_p / w order and the peak's place, in all four runs, and Y at w tau = 0.1 in
//     the two runs where the resolution bound lies under a tenth of the tolerance;
//   - 10 are NOT EVALUATED: Y at w tau = 1e-2 and the C_p order (taken there) in all four runs, and Y at
//     w tau = 0.1 in two; the admittance's resolution bound (5e-8 to 9e-8 S/cm^2) is above a tenth of the tolerance.
//     At w tau = 0.1 the bound is within a factor 2 of that threshold in every run, so those rows are held to
//     "pass or not_evaluated".
// Every measured error is pinned: |Y - Y_ref| is 4.5e-4 to 2.8e-3 of the tolerance, the orders 1.89 to 2.07 and
// 2.00, the peak's place within 3e-6.
// The first construction (17.18: 9 pass, 11 not_evaluated, 20 not_applicable) is in the history at eb73e39.
TEST_CASE("TVB-22: the trap level's admittance, full tier (frozen outcome: as found, 17.20)", "[.tvr5a_full]") {
    const auto r = bench::run_tvb22(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("%s note: %s\n", r.id.c_str(), n.c_str());
    REQUIRE(r.comparisons.size() == 40);
    const auto pass = Outcome::pass, ne = Outcome::not_evaluated;
    const char* wt[] = {"0.01", "0.1", "0.3", "1", "3", "10", "100"};
    std::size_t i = 0;
    int passes = 0;
    for (const char* type : {"acceptor", "donor"}) {
        for (const char* bias : {"-0.2", "-0.25"}) {
            const std::string at = std::string(type) + ", V_FB " + bias + " V";
            // The row at index i: its name, the outcomes allowed, and a window on the metric.
            const auto row = [&](const std::string& quantity, bool borderline, Outcome judged, double above,
                                 double below) {
                const auto& c = r.comparisons[i++];
                std::printf("%s [%s] %s: %.4g (%s) %s\n", r.id.c_str(),
                            std::string(to_string(c.comparison.outcome)).c_str(), c.quantity.c_str(), c.metric,
                            c.criterion.c_str(), c.comparison.reason.c_str());
                CAPTURE(c.quantity, c.metric, c.comparison.reason);
                CHECK(c.quantity == quantity);
                if (borderline) {
                    CHECK((c.comparison.outcome == pass || c.comparison.outcome == ne));
                } else {
                    CHECK(c.comparison.outcome == judged);
                }
                if (c.comparison.outcome == ne) CHECK(c.comparison.reason.find("admittance resolution") == 0);
                if (c.comparison.outcome == pass) ++passes;
                CHECK(c.metric > above);
                CHECK(c.metric < below);
            };
            row("G_p / w order (w tau = 1), " + at, false, pass, 1.85, 2.1);
            row("C_p order (w tau = 1e-2), " + at, false, ne, 1.99, 2.02);
            for (std::size_t k = 0; k < 7; ++k) {
                // The metric is |Y - Y_ref| / (w C_T) over the tolerance.
                row(std::string("Y finest (N = 512), w tau = ") + wt[k] + ", " + at, k == 1, k == 0 ? ne : pass, 2e-4,
                    4e-3);
            }
            row("G_p / w peak place finest (N = 512), " + at, false, pass, 4e-7, 4e-6);
        }
    }
    // The w tau = 1 error of each run (rows 5, 15, 25, 35), within a quarter of its value as found.
    const double at_one[] = {7.197e-4, 1.833e-3, 1.091e-3, 7.429e-4};
    for (std::size_t j = 0; j < 4; ++j) {
        const double m = r.comparisons[10 * j + 5].metric;
        CAPTURE(j, m);
        CHECK(m > 0.75 * at_one[j]);
        CHECK(m < 1.25 * at_one[j]);
    }
    CHECK(passes >= 28);
    CHECK(passes <= 32);
    CHECK(bench::overall(r) == ne);
}

TEST_CASE("TVB-22: nothing in the fast tier", "[tvr5a]") {
    const auto r = bench::run_tvb22(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-23 (avalanche multiplication), as rebuilt once under IG-Q24 B (17.22: M from the generation's share; V(M_c) by
// a bias root). The expectation before the rerun: all 26 pass. The first construction's outcome as found (17.21: 9
// pass, 7 fail, 10 not_evaluated) is in the history at 94e7256.
TEST_CASE("TVB-23: avalanche multiplication, full tier (expectation before the rerun, 17.22)", "[.tvr5a_full]") {
    const auto r = bench::run_tvb23(bench::Tier::full);
    std::vector<Found> expected;
    for (const char* carrier : {"M_n", "M_p"}) {
        for (const char* I : {"0.1", "0.3", "0.5", "0.7", "0.9"}) {
            // The bias in the name is the reference's; matched by prefix below.
            expected.push_back({std::string("1 - 1/M order, ") + carrier + ", I_n " + I, Outcome::pass});
            expected.push_back({std::string("1 - 1/M finest (N = 1024), ") + carrier + ", I_n " + I, Outcome::pass});
        }
    }
    for (const char* carrier : {"M_n", "M_p"}) {
        for (const char* M : {"10", "100", "1000"}) {
            expected.push_back({std::string("V at ") + carrier + " = " + M + " (A), N = 1024", Outcome::pass});
        }
    }
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("%s note: %s\n", r.id.c_str(), n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("%s [%s] %s: %.4g (%s) %s\n", r.id.c_str(), std::string(to_string(c.comparison.outcome)).c_str(),
                    c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const auto& c = r.comparisons[i];
        CAPTURE(c.quantity, c.metric, c.comparison.reason);
        CHECK(c.quantity.starts_with(expected[i].quantity));
        CHECK(c.comparison.outcome == expected[i].outcome);
    }
}

TEST_CASE("TVB-23: nothing in the fast tier", "[tvr5a]") {
    const auto r = bench::run_tvb23(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

TEST_CASE("TVB-08: nothing in the fast tier", "[tvr5a]") {
    const auto r = bench::run_tvb08(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}
