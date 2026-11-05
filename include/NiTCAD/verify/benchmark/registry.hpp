// The benchmark registry, known findings, summary and report of TVR-1 (Unit 55; DECISIONS 11.7, 11.8, owner rulings
// RC-Q1, RC-Q2, RC-Q4, RC-Q6 in 11.13).
//
// A benchmark is a function that runs the code under test against its independent reference and returns every
// comparison judged through verify::judge (the contract first). The tests call the same functions and hold the
// outcomes to the frozen expectations; `nitcad verify` runs them and writes the report. The reference's own L0 tests
// (TC1-TC9, metrics, contract) are ctest's, named in the report as build-time evidence, not re-run.
//
// Known findings (RC-Q2): a fail inside a region an owner ruling has recorded becomes `finding` with the ruling cited,
// so it is counted without failing the run (TVD-23); any other fail stays `fail`. A known finding whose region holds
// no fail any more leaves its comparisons as they are and adds the note "known finding (x) no longer reproduces:
// remove the entry".
#pragma once

#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/io/verify_report.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

namespace NiTCAD::verify::benchmark {

enum class Tier { fast, full, nightly };

struct ComparisonRecord {
    std::string quantity;  // e.g. "F_-1/2", "C3 [pn_diode_1d] bias (0.1) V"
    std::string unit;
    double parameter = 0.0;  // the swept argument (eta, x, T), or 0
    double q_sim = 0.0;
    double q_ref = 0.0;
    double u_ref = 0.0;
    double metric = 0.0;  // the judged value (a relative error, a residual / tolerance ratio, ...)
    std::string criterion;
    Comparison comparison;
};

struct BenchmarkResult {
    std::string id;
    std::string title;
    ReferenceInfo info;
    std::vector<ComparisonRecord> comparisons;
    std::vector<std::string> notes;
    std::optional<std::string> run_error;  // the benchmark could not run (a solve or a deck failed)
};

struct Benchmark {
    std::string_view id;
    std::string_view title;
    Tier tier;
    BenchmarkResult (*run)();
};

// The shipped benchmarks of TVR-1: TVB-01, TVB-02, TVB-03, TVB-32 (all fast).
[[nodiscard]] std::span<const Benchmark> benchmarks();

// The benchmark functions, for the tests.
[[nodiscard]] BenchmarkResult run_tvb01();
[[nodiscard]] BenchmarkResult run_tvb02();
[[nodiscard]] BenchmarkResult run_tvb03();
[[nodiscard]] BenchmarkResult run_tvb11();  // TVR-2 (DECISIONS Unit 55 14)
[[nodiscard]] BenchmarkResult run_tvb12();
[[nodiscard]] BenchmarkResult run_tvb32();

struct KnownFinding {
    std::string_view benchmark;
    std::string_view quantity;           // a comparison's quantity starts with this
    bool (*in_region)(double parameter);  // the ruled region
    std::string_view label;              // "(i)", "(ii)", "(iii)"
    std::string_view ruling;             // where it is recorded
};

// Findings (i)-(iii) of DECISIONS 11.13.
[[nodiscard]] std::span<const KnownFinding> known_findings();

// Applies a table of known findings to a result (fail -> finding inside a region; the stale note).
void apply_known_findings(BenchmarkResult& r, std::span<const KnownFinding> table);

// The overall outcome: fail if any comparison fails or the run failed; else not_evaluated if any is; else finding if
// any is; else pass if any passes; else not_applicable (printed "unverified").
[[nodiscard]] Outcome overall(const BenchmarkResult& r);

struct Summary {
    int pass = 0, finding = 0, fail = 0, not_applicable = 0, not_evaluated = 0, unverified = 0, run_errors = 0;
};
[[nodiscard]] Summary summarize(std::span<const BenchmarkResult> results);
// "verify: 4 benchmarks; 2 pass, 2 finding, 0 fail, ..." (one line, no newline).
[[nodiscard]] std::string summary_line(const Summary& s, std::size_t benchmarks);

// The report (RC-Q4) and its CSV summary: one header line, then one line per comparison.
[[nodiscard]] io::VerifyReport make_report(std::span<const BenchmarkResult> results, results::Provenance provenance,
                                           std::uint32_t threads, std::string tier);
[[nodiscard]] std::string report_csv(const io::VerifyReport& report);

// The shortest round-trip form of a double (std::to_chars), as the contract's reasons print numbers.
[[nodiscard]] std::string shortest(double x);

}  // namespace NiTCAD::verify::benchmark
