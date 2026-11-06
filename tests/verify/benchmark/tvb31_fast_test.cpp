// TVB-31.06-31.31 in the fast tier (DECISIONS Unit 55 15.20, MI-Q9): one point a decade and every 16th point of the
// linear grids, with every switch point and neighbour kept, so each fast point is a point of the full grid, whose
// outcomes are frozen in the [tvb31_full] tests (Release). The rules here are the per-comparison rules those tests
// hold, without their counts:
//   - no run error and no fail; every comparison well formed;
//   - not_evaluated only by TB-Q4 (a reason starting "U_ref ");
//   - notes only 31.25's diagnostics, or a known-finding staleness note (a fast grid may miss a finding's region;
//     staleness is judged in the full tier);
//   - at least one comparison in every entry but 31.25 (diagnostics only).
// TVB-31.01-31.05 take no tier: their frozen tests run in both builds.
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <string_view>

namespace bench = NiTCAD::verify::benchmark;
namespace v = NiTCAD::verify;

namespace {
bool starts_with(std::string_view s, std::string_view p) { return s.substr(0, p.size()) == p; }
bool ends_with(std::string_view s, std::string_view p) {
    return s.size() >= p.size() && s.substr(s.size() - p.size()) == p;
}
}  // namespace

TEST_CASE("TVB-31.06-31.31: fast tier", "[verify][tvb31]") {
    int ran = 0;
    for (const bench::Benchmark& b : bench::benchmarks()) {
        const std::string_view id = b.id;
        if (!starts_with(id, "TVB-31.") || id < std::string_view("TVB-31.06")) continue;
        INFO(id);
        const bench::BenchmarkResult r = b.run(bench::Tier::fast);
        ++ran;
        CHECK_FALSE(r.run_error.has_value());
        if (id != "TVB-31.25") CHECK_FALSE(r.comparisons.empty());
        for (const std::string& note : r.notes) {
            INFO(note);
            CHECK((starts_with(note, "diagnostic (not judged, MI-Q8): ") ||
                   (starts_with(note, "known finding ") && ends_with(note, "no longer reproduces: remove the entry"))));
        }
        for (const auto& c : r.comparisons) {
            INFO(c.quantity << " at " << c.parameter << ": " << v::to_string(c.comparison.outcome) << " ("
                            << c.comparison.reason << ")");
            CHECK(v::validate(c.comparison).empty());
            CHECK(c.comparison.outcome != v::Outcome::fail);
            if (c.comparison.outcome == v::Outcome::not_evaluated) CHECK(starts_with(c.comparison.reason, "U_ref "));
        }
    }
    CHECK(ran == 26);
}
