// Shared outcome check of the TVB-31 entry tests (Unit 55, TVR-3): no run error and no stale-finding note; every
// comparison well formed; a pass is never an expected finding (unless `pass_allowed` says the region may pass at its
// edge); a finding is always expected; fail is refused, and not_evaluated unless the caller allows it (TB-Q4). Returns
// the counts.
#pragma once

#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string_view>

namespace tvb31_test {

namespace v = NiTCAD::verify;
namespace bench = NiTCAD::verify::benchmark;

struct Counts {
    int pass = 0, finding = 0, not_applicable = 0, not_evaluated = 0, other = 0;
};

inline bool starts_with(std::string_view s, std::string_view p) { return s.substr(0, p.size()) == p; }

inline bool never(const bench::ComparisonRecord&) { return false; }

// `allowed_not_evaluated` names the comparisons that may be not_evaluated (TB-Q4; none by default).
template <class Expected, class AllowedNE = bool (*)(const bench::ComparisonRecord&)>
Counts check(const bench::BenchmarkResult& r, Expected&& expected_finding, bool pass_allowed = false,
             AllowedNE allowed_not_evaluated = never) {
    INFO(r.id);
    CHECK_FALSE(r.run_error.has_value());
    CHECK(r.notes.empty());
    Counts n;
    for (const auto& c : r.comparisons) {
        INFO(c.quantity << " at " << c.parameter << ": " << v::to_string(c.comparison.outcome) << " ("
                        << c.comparison.reason << "), sim " << c.q_sim << ", ref " << c.q_ref << ", metric "
                        << c.metric);
        CHECK(v::validate(c.comparison).empty());
        switch (c.comparison.outcome) {
            case v::Outcome::pass:
                if (!pass_allowed) CHECK_FALSE(expected_finding(c));
                ++n.pass;
                break;
            case v::Outcome::finding:
                CHECK(expected_finding(c));
                ++n.finding;
                break;
            case v::Outcome::not_applicable: ++n.not_applicable; break;
            case v::Outcome::not_evaluated:
                if (allowed_not_evaluated(c)) {
                    ++n.not_evaluated;
                    break;
                }
                ++n.other;
                FAIL_CHECK("not_evaluated where none is expected");
                break;
            default:
                ++n.other;
                FAIL_CHECK("fail");
                break;
        }
    }
    return n;
}

}  // namespace tvb31_test
