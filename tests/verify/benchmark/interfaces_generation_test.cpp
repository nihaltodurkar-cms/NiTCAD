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

// TVB-22 (interface trap level). The expectation written before the first run was that all 40 comparisons pass
// (two trap types x two biases, each two orders, Y at seven frequencies and the peak's place). It was not met.
// Recorded as found (17.18, IG-Q21 A), nothing changed after the run:
//   - at V_FB - 0.3 V all 20 comparisons are NOT APPLICABLE: P1 (the minority carriers' share of the level's
//     exchange) is 3.1e-6 and 2.7e-6 against the frozen limit 1e-6;
//   - at V_FB - 0.2 V the comparisons at w tau <= 0.3, and the C_p order (taken at w tau = 1e-2), are NOT EVALUATED:
//     the admittance's resolution bound (about 1e-7 S/cm^2) is above a tenth of the tolerance there; at w tau = 1
//     the bound is within 10 % of that threshold on either side, so those rows (Y, the G_p / w order, the peak) are
//     held to "pass or not_evaluated"; at w tau >= 3 every comparison passes.
// Every measured error is far inside its criterion wherever it is formed (pinned below): |Y - Y_ref| at most 1.5 %
// of the tolerance at -0.2 V and 5.3 % at -0.3 V, orders 1.99 to 2.07, the peak's place within 3e-6.
TEST_CASE("TVB-22: the trap level's admittance, full tier (frozen outcome: as found, 17.18)", "[.tvr5a_full]") {
    const auto r = bench::run_tvb22(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("%s note: %s\n", r.id.c_str(), n.c_str());
    REQUIRE(r.comparisons.size() == 40);
    const auto pass = Outcome::pass, ne = Outcome::not_evaluated, na = Outcome::not_applicable;
    const char* wt[] = {"0.01", "0.1", "0.3", "1", "3", "10", "100"};
    std::size_t i = 0;
    int passes = 0;
    for (const char* type : {"acceptor", "donor"}) {
        for (const char* bias : {"-0.2", "-0.3"}) {
            const std::string at = std::string(type) + ", V_FB " + bias + " V";
            const bool shallow = std::string(bias) == "-0.2";
            // The row at index i: its name, the outcomes allowed, and a window on the metric.
            const auto row = [&](const std::string& quantity, bool borderline, Outcome judged, double above,
                                 double below) {
                const auto& c = r.comparisons[i++];
                std::printf("%s [%s] %s: %.4g (%s) %s\n", r.id.c_str(),
                            std::string(to_string(c.comparison.outcome)).c_str(), c.quantity.c_str(), c.metric,
                            c.criterion.c_str(), c.comparison.reason.c_str());
                CAPTURE(c.quantity, c.metric, c.comparison.reason);
                CHECK(c.quantity == quantity);
                if (!shallow) {
                    CHECK(c.comparison.outcome == na);
                    CHECK(c.comparison.reason.find("P1 ") == 0);
                } else if (borderline) {
                    CHECK((c.comparison.outcome == pass || c.comparison.outcome == ne));
                } else {
                    CHECK(c.comparison.outcome == judged);
                }
                if (c.comparison.outcome == pass) ++passes;
                CHECK(c.metric > above);
                CHECK(c.metric < below);
            };
            row("G_p / w order (w tau = 1), " + at, true, pass, 1.95, 2.1);
            row("C_p order (w tau = 1e-2), " + at, false, ne, 1.99, 2.01);
            for (std::size_t k = 0; k < 7; ++k) {
                // The metric is |Y - Y_ref| / (w C_T) over the tolerance.
                row(std::string("Y finest (N = 512), w tau = ") + wt[k] + ", " + at, k == 3, k < 3 ? ne : pass,
                    shallow ? 0.010 : 0.040, shallow ? 0.017 : 0.058);
            }
            row("G_p / w peak place finest (N = 512), " + at, true, pass, 1e-6, 5e-6);
        }
    }
    CHECK(passes >= 6);  // the six rows at w tau >= 3 of the two -0.2 V runs
    CHECK(passes <= 12);
    CHECK(bench::overall(r) == ne);
}

TEST_CASE("TVB-22: nothing in the fast tier", "[tvr5a]") {
    const auto r = bench::run_tvb22(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

TEST_CASE("TVB-08: nothing in the fast tier", "[tvr5a]") {
    const auto r = bench::run_tvb08(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}
