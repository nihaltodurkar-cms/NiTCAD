// The registry, known findings, summary and report of TVR-1 (DECISIONS Unit 55 11.7; RC-Q1, RC-Q2, RC-Q4, RC-Q6).
#include "NiTCAD/io/verify_report.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
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

TEST_CASE("registry: the shipped benchmarks and the tier each starts in") {
    const auto b = bench::benchmarks();
    // TVR-1, TVR-2's TVB-04, 07, 09, 11, 12 and 14 (DECISIONS Unit 55 14), TVR-3's TVB-31.01-31, TVR-4's TVB-05, 06, 15-18 (16).
    REQUIRE(b.size() == 47);
    const char* ids[] = {"TVB-01",    "TVB-02",    "TVB-03",    "TVB-04",    "TVB-05",    "TVB-06",    "TVB-07",    "TVB-09",
                         "TVB-11",    "TVB-12",    "TVB-14",    "TVB-15",    "TVB-16",    "TVB-17",    "TVB-18",    "TVB-31.01",
                         "TVB-31.02", "TVB-31.03", "TVB-31.04", "TVB-31.05", "TVB-31.06", "TVB-31.07", "TVB-31.08", "TVB-31.09",
                         "TVB-31.10", "TVB-31.11", "TVB-31.12", "TVB-31.13", "TVB-31.14", "TVB-31.15",
                         "TVB-31.16", "TVB-31.17", "TVB-31.18", "TVB-31.19", "TVB-31.20", "TVB-31.21",
                         "TVB-31.22", "TVB-31.23", "TVB-31.24", "TVB-31.25", "TVB-31.26", "TVB-31.27",
                         "TVB-31.28", "TVB-31.29", "TVB-31.30", "TVB-31.31", "TVB-32"};
    REQUIRE(std::size(ids) == b.size());
    for (std::size_t i = 0; i < b.size(); ++i) {
        CHECK(b[i].id == ids[i]);
        // TVB-05, 15, 16, 17 and 18 start in full (DECISIONS Unit 55 16.7, DT-Q9); every other benchmark in fast.
        const bool full = b[i].id == "TVB-05" || b[i].id == "TVB-15" || b[i].id == "TVB-16" ||
                          b[i].id == "TVB-17" || b[i].id == "TVB-18";
        CHECK(b[i].tier == (full ? bench::Tier::full : bench::Tier::fast));
        CHECK(b[i].run != nullptr);
    }
    // (iii) and (i) fixed in Unit 56 (FQ4, FQ1); (ii); TVB-04 (iv) in 2D and 3D (CV-Q12); TVB-14 (iv), 6 (CV-Q20),
    // and (v), 6 (CV-Q21); TVB-31.04 (vi), 2 (MI-Q12); TVB-31.05 (F1) and TVB-31.07 (F6) (15.4, MI-Q4); TVB-31.10 (vii), 6 (MI-Q14); TVB-16 (x), 2 (DT-Q27 A).
    // TVB-31.12 (viii) and TVB-31.15 (ix) (MI-Q15) are classified inside their entries, by the rounding floor of the
    // shipped form, not by a region of the parameter.
    REQUIRE(bench::known_findings().size() == 27);
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

// The (h, tol) factorial helper (DECISIONS 4.1 (b), (c); Unit 55 16.2 (p), DT-Q23 A).
TEST_CASE("registry: the (h, tol) factorial and its block on an E pass", "[registry]") {
    using namespace NiTCAD::verify;
    // Additive errors: e = a(h) + b(tol) has no interaction.
    const std::vector<double> parts{1e-4, 2e-6};
    const auto add = benchmark::factorial(1e-3 + 1e-7, 1e-3 + 1e-9, 2.5e-4 + 1e-7, 2.5e-4 + 1e-9, parts);
    REQUIRE(add.ok);
    CHECK(std::abs(add.eps_int) <= 1e-18);
    CHECK_FALSE(add.interacting);
    CHECK(std::abs(add.u - (1e-4 + 2e-6 + std::abs(add.eps_int))) <= 1e-19);
    // An interaction: the tolerance effect grows on the finer mesh. eps_int = -(4e-5 - 1e-5),
    // U = 3e-5 + 1e-4 + 2e-6.
    const auto mix = benchmark::factorial(1e-3, 1e-3 - 1e-5, 2.5e-4, 2.5e-4 - 4e-5, parts);
    REQUIRE(mix.ok);
    CHECK(std::abs(mix.eps_int - (-3e-5)) <= 1e-18);
    CHECK(std::abs(mix.u - 1.32e-4) <= 1e-18);
    CHECK(mix.interacting);
    CHECK(mix.note.find("interacting errors") != std::string::npos);
    // The block: a pass becomes not_evaluated; a fail is kept; nothing changes without an interaction.
    benchmark::ComparisonRecord c;
    c.comparison = Comparison{Outcome::pass, "", {}};
    benchmark::block_if_interacting(c, add);
    CHECK(c.comparison.outcome == Outcome::pass);
    benchmark::block_if_interacting(c, mix);
    CHECK(c.comparison.outcome == Outcome::not_evaluated);
    CHECK(c.comparison.reason.find("interacting errors") != std::string::npos);
    c.comparison = Comparison{Outcome::fail, "x", {}};
    benchmark::block_if_interacting(c, mix);
    CHECK(c.comparison.outcome == Outcome::fail);
    c.comparison = Comparison{Outcome::not_applicable, "predicate", {}};
    benchmark::block_if_interacting(c, mix);
    CHECK(c.comparison.outcome == Outcome::not_applicable);
    CHECK(c.comparison.reason == "predicate");
    // Not finite: not ok, and a pass is not kept.
    const auto bad = benchmark::factorial(std::nan(""), 1.0, 1.0, 1.0, parts);
    CHECK_FALSE(bad.ok);
    c.comparison = Comparison{Outcome::pass, "", {}};
    benchmark::block_if_interacting(c, bad);
    CHECK(c.comparison.outcome == Outcome::not_evaluated);
    CHECK_FALSE(benchmark::factorial(1.0, 1.0, 1.0, 1.0, std::vector<double>{INFINITY}).ok);
}

// The shared U_disc rule (DECISIONS 4.2, 4.3; Unit 55 DT-Q26 A, DT-Q28 B, DT-Q32 A) on exact sequences f(h) = c + a h^p.
TEST_CASE("registry: U_disc, the in-range GCI and the two-grid fallback", "[registry]") {
    const double c = 1.0, a = 1e-3, h0 = 0.01;
    const std::array<double, 4> h{h0, 2 * h0, 4 * h0, 8 * h0};
    const auto quad = [&](double x) { return c + a * x * x; };
    // f(h/2) giving the finest triple (h/2, h, 2h) the observed order q: f0 - f(h/2) = (f1 - f0) / 2^q.
    const auto half_for = [&](double q) { return quad(h0) - 3.0 * a * h0 * h0 / std::pow(2.0, q); };
    const double tol = 1e-8 * a * h0 * h0;  // f ~ 1 + 1e-7: the differences keep ~8 digits
    // Order 2 everywhere: GCI = 1.25 |f1 - f0| / 3 = 1.25 a h0^2; two-grid = 3 |f(h/2) - f0| 4/3 = 3 a h0^2.
    {
        const std::array<double, 4> f{quad(h[0]), quad(h[1]), quad(h[2]), quad(h[3])};
        const auto d = bench::u_disc(h, f, h0 / 2, quad(h0 / 2));
        REQUIRE(d.u.has_value());
        CHECK(std::abs(d.order - 2.0) <= 1e-6);
        CHECK(std::abs(d.order_coarse - 2.0) <= 1e-6);
        CHECK(std::abs(d.order_fine - 2.0) <= 1e-6);
        CHECK(std::abs(*d.u - 3.0 * a * h0 * h0) <= tol);
        CHECK(d.how == "two-grid (> GCI)");
        // The finest triple at order 1 (f(h/2) off the h^2 line): only the GCI is kept.
        const auto g = bench::u_disc(h, f, h0 / 2, half_for(1.0));
        REQUIRE(g.u.has_value());
        CHECK(std::abs(g.order_fine - 1.0) <= 1e-6);
        CHECK(std::abs(*g.u - 1.25 * a * h0 * h0) <= tol);
        CHECK(g.how.rfind("GCI; two-grid: finest-triple order 1", 0) == 0);
        // At order 3: not kept either.
        CHECK(bench::u_disc(h, f, h0 / 2, half_for(3.0)).how.rfind("GCI; two-grid: finest-triple order 3", 0) == 0);
        // f(h/2) past f(h)'s side of the limit: the finest triple oscillates.
        const auto o = bench::u_disc(h, f, h0 / 2, quad(h0) + a * h0 * h0);
        CHECK(o.how == "GCI; two-grid: the finest triple oscillates or has no positive order");
    }
    // Order 2 on the finest triples, a coarsest level off the h^2 line: out of range; the two-grid estimate stands in.
    {
        const std::array<double, 4> f{quad(h[0]), quad(h[1]), quad(h[2]), c + 2.0 * a * h[3] * h[3]};
        const auto e = bench::u_disc(h, f, h0 / 2, quad(h0 / 2));
        REQUIRE(e.u.has_value());
        CHECK(std::abs(*e.u - 3.0 * a * h0 * h0) <= tol);
        CHECK(e.how.rfind("two-grid; outside the asymptotic range (orders 2", 0) == 0);
        // With the finest triple at order 1 as well: nothing.
        const auto d = bench::u_disc(h, f, h0 / 2, half_for(1.0));
        CHECK_FALSE(d.u.has_value());
        CHECK(d.how.rfind("outside the asymptotic range (orders 2", 0) == 0);
        CHECK(d.how.find("; two-grid: finest-triple order 1") != std::string::npos);
    }
    // Orders 2 and about 2.3 (the coarsest level 30 % further off): outside 0.1, though within 1.
    {
        const std::array<double, 4> f{quad(h[0]), quad(h[1]), quad(h[2]), c + 1.3 * a * h[3] * h[3]};
        const auto d = bench::u_disc(h, f, h0 / 2, half_for(1.0));
        CHECK(d.order_coarse - d.order > 0.1);
        CHECK(d.order_coarse - d.order < 1.0);
        CHECK_FALSE(d.u.has_value());
        CHECK(d.how.rfind("outside the asymptotic range", 0) == 0);
    }
    // The band's edges: orders 1.45 and 2.55 not kept, 1.55 and 2.45 kept, at 4 |f(h/2) - f0|.
    {
        const std::array<double, 4> f{quad(h[0]), quad(h[1]), quad(h[2]), quad(h[1])};  // no three-grid estimate
        CHECK_FALSE(bench::u_disc(h, f, h0 / 2, half_for(1.45)).u.has_value());
        CHECK_FALSE(bench::u_disc(h, f, h0 / 2, half_for(2.55)).u.has_value());
        for (const double q : {1.55, 2.45}) {
            const auto d = bench::u_disc(h, f, h0 / 2, half_for(q));
            REQUIRE(d.u.has_value());
            CHECK(std::abs(*d.u - 4.0 * 3.0 * a * h0 * h0 / std::pow(2.0, q)) <= tol);
        }
    }
    // The coarser triple oscillates: no order there.
    {
        const std::array<double, 4> f{quad(h[0]), quad(h[1]), quad(h[2]), quad(h[1])};
        const auto d = bench::u_disc(h, f, h0 / 2, half_for(1.0));
        CHECK_FALSE(d.u.has_value());
        CHECK(d.how.rfind("asymptotic check: no order on the coarser triple", 0) == 0);
    }
    // The finest triple of the ladder oscillates.
    {
        const std::array<double, 4> f{1.0, 1.1, 1.05, 1.2};
        const auto d = bench::u_disc(h, f, h0 / 2, 1.0);
        CHECK_FALSE(d.u.has_value());
        CHECK(d.how.rfind("GCI: oscillatory or a non-positive order", 0) == 0);
    }
}
