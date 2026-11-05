// The registry, known findings, summary and report of TVR-1 (DECISIONS Unit 55 11.7; RC-Q1, RC-Q2, RC-Q4, RC-Q6).
#include "NiTCAD/io/verify_report.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

namespace bench = NiTCAD::verify::benchmark;
namespace io = NiTCAD::io;
using NiTCAD::verify::Outcome;

namespace {

bench::ComparisonRecord record(const char* quantity, double parameter, Outcome o) {
    bench::ComparisonRecord c;
    c.quantity = quantity;
    c.parameter = parameter;
    c.comparison.outcome = o;
    c.comparison.reason = o == Outcome::pass ? "" : "criterion missed";
    return c;
}

bool above_ten(double p) { return p > 10.0; }

constexpr bench::KnownFinding table[] = {{"TVB-X", "q", &above_ten, "(z)", "a ruling"}};

}  // namespace

TEST_CASE("registry: the shipped benchmarks, each starting in the fast tier") {
    const auto b = bench::benchmarks();
    REQUIRE(b.size() == 8);  // TVR-1 plus TVR-2's TVB-04, TVB-07, TVB-11 and TVB-12 (DECISIONS Unit 55 14)
    const char* ids[] = {"TVB-01", "TVB-02", "TVB-03", "TVB-04", "TVB-07", "TVB-11", "TVB-12", "TVB-32"};
    for (std::size_t i = 0; i < b.size(); ++i) {
        CHECK(b[i].id == ids[i]);
        CHECK(b[i].tier == bench::Tier::fast);
        CHECK(b[i].run != nullptr);
    }
    // (iii) and (i) fixed in Unit 56 (FQ4, FQ1); (ii); TVB-04 (iv) in 2D and 3D (CV-Q12)
    REQUIRE(bench::known_findings().size() == 3);
}

TEST_CASE("registry: a known finding turns a fail inside its region into finding, and only there") {
    bench::BenchmarkResult r;
    r.id = "TVB-X";
    r.comparisons = {record("q", 20.0, Outcome::fail), record("q", 5.0, Outcome::fail),
                     record("other", 20.0, Outcome::fail), record("q", 30.0, Outcome::pass)};
    bench::apply_known_findings(r, table);
    CHECK(r.comparisons[0].comparison.outcome == Outcome::finding);
    CHECK(r.comparisons[0].comparison.reason == "known finding (z) (a ruling): criterion missed");
    CHECK(r.comparisons[1].comparison.outcome == Outcome::fail);  // outside the region
    CHECK(r.comparisons[2].comparison.outcome == Outcome::fail);  // another quantity
    CHECK(r.comparisons[3].comparison.outcome == Outcome::pass);
    CHECK(r.notes.empty());
    CHECK(bench::overall(r) == Outcome::fail);

    // Another benchmark's result is untouched.
    bench::BenchmarkResult other = r;
    other.id = "TVB-Y";
    other.comparisons[1].comparison.outcome = Outcome::fail;
    bench::apply_known_findings(other, table);
    CHECK(other.comparisons[1].comparison.outcome == Outcome::fail);
}

TEST_CASE("registry: a known finding that no longer reproduces adds the note") {
    bench::BenchmarkResult r;
    r.id = "TVB-X";
    r.comparisons = {record("q", 20.0, Outcome::pass), record("q", 25.0, Outcome::pass)};
    bench::apply_known_findings(r, table);
    REQUIRE(r.notes.size() == 1);
    CHECK(r.notes[0] == "known finding (z) no longer reproduces: remove the entry");
    CHECK(bench::overall(r) == Outcome::pass);
}

TEST_CASE("registry: the overall outcome and the summary") {
    bench::BenchmarkResult r;
    r.id = "TVB-X";
    CHECK(bench::overall(r) == Outcome::not_applicable);  // nothing judged: "unverified"
    r.comparisons = {record("a", 0, Outcome::not_applicable)};
    CHECK(bench::overall(r) == Outcome::not_applicable);
    r.comparisons.push_back(record("b", 0, Outcome::pass));
    CHECK(bench::overall(r) == Outcome::pass);
    r.comparisons.push_back(record("c", 0, Outcome::finding));
    CHECK(bench::overall(r) == Outcome::finding);
    r.comparisons.push_back(record("d", 0, Outcome::not_evaluated));
    CHECK(bench::overall(r) == Outcome::not_evaluated);
    r.comparisons.push_back(record("e", 0, Outcome::fail));
    CHECK(bench::overall(r) == Outcome::fail);
    bench::BenchmarkResult broken;
    broken.run_error = "the deck did not run";
    CHECK(bench::overall(broken) == Outcome::fail);

    bench::BenchmarkResult na;
    na.comparisons = {record("a", 0, Outcome::not_applicable)};
    const std::vector<bench::BenchmarkResult> all{r, broken, na};
    const bench::Summary s = bench::summarize(all);
    CHECK(s.fail == 2);
    CHECK(s.run_errors == 1);
    CHECK(s.unverified == 1);
    CHECK(bench::summary_line(s, all.size()) ==
          "verify: 3 benchmarks; 0 pass, 0 finding, 2 fail, 1 not applicable (1 unverified), 0 not evaluated, "
          "1 run errors");
}

TEST_CASE("registry: the CSV quotes fields with commas and quotes, and prints shortest round-trip numbers") {
    io::VerifyReport rep;
    io::VerifyBenchmark b;
    b.id = "TVB-X";
    b.comparisons.push_back({"q, at 1", "V", 0.1, 700.0, 1.0 / 3.0, 0.0, "not estimated", 2.5e-16,
                             "say \"x\"", "pass", ""});
    rep.benchmarks.push_back(b);
    CHECK(bench::report_csv(rep) ==
          "benchmark,quantity,unit,parameter,q_sim,q_ref,u_ref,u_num,metric,criterion,outcome,reason\n"
          "TVB-X,\"q, at 1\",V,0.1,700,0.3333333333333333,0,not estimated,2.5e-16,\"say \"\"x\"\"\",pass,\n");
}

TEST_CASE("registry: the report round-trips through io, and two runs are equal except the clock (VR7)") {
    auto run = [] {
        const std::vector<bench::BenchmarkResult> results{bench::run_tvb01(), bench::run_tvb03()};
        NiTCAD::results::Provenance p;
        p.version = "test";
        p.started_utc = "2026-11-05T00:00:00.000Z";
        p.wall_s = 1.25;
        return bench::make_report(results, p, 1, "fast");
    };
    const io::VerifyReport a = run();
    io::VerifyReport b = run();
    b.provenance.started_utc = "2026-11-05T01:02:03.456Z";
    b.provenance.wall_s = 9.5;
    CHECK(a != b);
    b.provenance.started_utc = a.provenance.started_utc;
    b.provenance.wall_s = a.provenance.wall_s;
    CHECK(a == b);
    REQUIRE(a.benchmarks.size() == 2);
    CHECK(a.benchmarks[0].outcome == "pass");  // TVB-01: finding (i) fixed in Unit 56 (FQ1)
    CHECK(a.benchmarks[1].outcome == "pass");  // TVB-03: finding (iii) fixed in Unit 56 (FQ4)
    CHECK(a.reference_evidence.size() == 2);

    const auto path = std::filesystem::temp_directory_path() / "nitcad_verify_report_test.nvr";
    REQUIRE(io::save_verify_report(path, a).has_value());
    const auto loaded = io::load_verify_report(path);
    REQUIRE(loaded.has_value());
    CHECK(*loaded == a);
    std::filesystem::remove(path);
}
