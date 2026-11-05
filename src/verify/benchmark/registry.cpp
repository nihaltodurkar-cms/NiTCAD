// The benchmark registry, known findings, summary and report (Unit 55, TVR-1; registry.hpp).
#include "NiTCAD/verify/benchmark/registry.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <utility>

namespace NiTCAD::verify::benchmark {

namespace {

constexpr std::array<Benchmark, 4> shipped{{
    {"TVB-01", "Constants and thermal voltage against CODATA 2018", Tier::fast, &run_tvb01},
    {"TVB-02", "Fermi-Dirac integrals F_1/2, F_-1/2, F_3/2 and the inverse of F_1/2", Tier::fast, &run_tvb02},
    {"TVB-03", "The Bernoulli function B(x) = x / (e^x - 1)", Tier::fast, &run_tvb03},
    {"TVB-32", "The consistency suite C1-C12 on the three TVR-1 decks", Tier::fast, &run_tvb32},
}};

bool everywhere(double) { return true; }

constexpr std::array<KnownFinding, 1> findings{{
    {"TVB-02", "F_-1/2", &everywhere, "(ii)",
     "DECISIONS 11.13 TVB-02 finding (ii) widened (option A), accepted as a limitation (Unit 56, 13): F_-1/2 "
     "misses from eta = 13.76, a quadrature-rounding floor of the derivative-consistent table"},
}};

bool starts_with(std::string_view s, std::string_view prefix) { return s.substr(0, prefix.size()) == prefix; }

std::string reference_text(const ReferenceInfo& i) {
    return std::string(to_string(i.cls)) + " " + std::string(to_string(i.independence)) + " " +
           std::string(to_string(i.provenance)) + "; " + i.citation + "; " + i.version;
}

std::string csv_field(std::string_view s) {
    if (s.find_first_of(",\"\n") == std::string_view::npos) return std::string(s);
    std::string out = "\"";
    for (const char c : s) {
        if (c == '"') out += '"';
        out += c;
    }
    return out + "\"";
}

}  // namespace

std::span<const Benchmark> benchmarks() { return shipped; }
std::span<const KnownFinding> known_findings() { return findings; }

void apply_known_findings(BenchmarkResult& r, std::span<const KnownFinding> table) {
    for (const KnownFinding& k : table) {
        if (k.benchmark != r.id) continue;
        int in_region = 0, failing = 0;
        for (ComparisonRecord& c : r.comparisons) {
            if (!starts_with(c.quantity, k.quantity) || !k.in_region(c.parameter)) continue;
            ++in_region;
            if (c.comparison.outcome != Outcome::fail) continue;
            ++failing;
            c.comparison.outcome = Outcome::finding;
            c.comparison.reason = "known finding " + std::string(k.label) + " (" + std::string(k.ruling) +
                                  "): " + c.comparison.reason;
        }
        if (in_region > 0 && failing == 0) {
            r.notes.push_back("known finding " + std::string(k.label) + " no longer reproduces: remove the entry");
        }
    }
}

Outcome overall(const BenchmarkResult& r) {
    if (r.run_error) return Outcome::fail;
    bool ne = false, finding = false, pass = false;
    for (const ComparisonRecord& c : r.comparisons) {
        switch (c.comparison.outcome) {
            case Outcome::fail: return Outcome::fail;
            case Outcome::not_evaluated: ne = true; break;
            case Outcome::finding: finding = true; break;
            case Outcome::pass: pass = true; break;
            case Outcome::not_applicable: break;
        }
    }
    if (ne) return Outcome::not_evaluated;
    if (finding) return Outcome::finding;
    if (pass) return Outcome::pass;
    return Outcome::not_applicable;
}

Summary summarize(std::span<const BenchmarkResult> results) {
    Summary s;
    for (const BenchmarkResult& r : results) {
        if (r.run_error) ++s.run_errors;
        switch (overall(r)) {
            case Outcome::pass: ++s.pass; break;
            case Outcome::finding: ++s.finding; break;
            case Outcome::fail: ++s.fail; break;
            case Outcome::not_evaluated: ++s.not_evaluated; break;
            case Outcome::not_applicable:
                ++s.not_applicable;
                ++s.unverified;
                break;
        }
    }
    return s;
}

std::string summary_line(const Summary& s, std::size_t n) {
    return "verify: " + std::to_string(n) + " benchmarks; " + std::to_string(s.pass) + " pass, " +
           std::to_string(s.finding) + " finding, " + std::to_string(s.fail) + " fail, " +
           std::to_string(s.not_applicable) + " not applicable (" + std::to_string(s.unverified) + " unverified), " +
           std::to_string(s.not_evaluated) + " not evaluated, " + std::to_string(s.run_errors) + " run errors";
}

io::VerifyReport make_report(std::span<const BenchmarkResult> results, results::Provenance provenance,
                             std::uint32_t threads, std::string tier) {
    io::VerifyReport rep;
    rep.provenance = std::move(provenance);
    rep.threads = threads;
    rep.tier = std::move(tier);
    rep.reference_evidence = {
        "ctest verify.* (TC1-TC9 L0 tests, metrics VR4, contract VR3): build-time, not re-run",
        "ctest verify_reference_boundary, verify_consistency_boundary (VR2): build-time"};
    for (const BenchmarkResult& r : results) {
        io::VerifyBenchmark b;
        b.id = r.id;
        b.title = r.title;
        b.reference = reference_text(r.info);
        const Outcome o = overall(r);
        b.outcome = o == Outcome::not_applicable ? "unverified" : std::string(to_string(o));
        b.run_error = r.run_error.value_or("");
        b.notes = r.notes;
        for (const ComparisonRecord& c : r.comparisons) {
            b.comparisons.push_back({c.quantity, c.unit, c.parameter, c.q_sim, c.q_ref, c.u_ref, "not estimated",
                                     c.metric, c.criterion, std::string(to_string(c.comparison.outcome)),
                                     c.comparison.reason});
        }
        rep.benchmarks.push_back(std::move(b));
    }
    return rep;
}

std::string report_csv(const io::VerifyReport& rep) {
    std::string s = "benchmark,quantity,unit,parameter,q_sim,q_ref,u_ref,u_num,metric,criterion,outcome,reason\n";
    for (const io::VerifyBenchmark& b : rep.benchmarks) {
        for (const io::VerifyComparison& c : b.comparisons) {
            s += csv_field(b.id) + "," + csv_field(c.quantity) + "," + csv_field(c.unit) + "," +
                 shortest(c.parameter) + "," + shortest(c.q_sim) + "," + shortest(c.q_ref) + "," + shortest(c.u_ref) +
                 "," + csv_field(c.u_num) + "," + shortest(c.metric) + "," + csv_field(c.criterion) + "," +
                 csv_field(c.outcome) + "," + csv_field(c.reason) + "\n";
        }
    }
    return s;
}

std::string shortest(double x) {
    char buf[64];
    const auto r = std::to_chars(buf, buf + sizeof buf, x);
    return std::string(buf, r.ptr);
}

}  // namespace NiTCAD::verify::benchmark
