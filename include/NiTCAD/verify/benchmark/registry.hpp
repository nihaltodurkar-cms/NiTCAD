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
    Tier tier;  // the lowest tier it runs in: `nitcad verify --tier t` runs every benchmark with tier <= t
    BenchmarkResult (*run)(Tier);  // at the requested tier (a ladder may grow with it; CV-Q9)
};

// The shipped benchmarks: TVR-1 (TVB-01, 02, 03, 32), TVR-2 (TVB-04, 07, 09, 11, 12, 14) and TVR-3 (TVB-31.xx), in
// id order.
[[nodiscard]] std::span<const Benchmark> benchmarks();

// The MKL PARDISO threads of the benchmark solves that use it (TVB-04's 3D levels; DECISIONS Unit 55 14.9), set by
// `nitcad verify --threads`; 1 by default. Those solves need the MKL runtime loaded (`--mkl`).
void set_pardiso_threads(int threads);
[[nodiscard]] int pardiso_threads();

// The benchmark functions, for the tests.
[[nodiscard]] BenchmarkResult run_tvb01();
[[nodiscard]] BenchmarkResult run_tvb02();
[[nodiscard]] BenchmarkResult run_tvb03();
[[nodiscard]] BenchmarkResult run_tvb04(Tier tier);  // fast: 1D; full and nightly: + 2D and 3D tensor
[[nodiscard]] BenchmarkResult run_tvb07(Tier tier);  // fast: N = 64, 128 (V_bi only); full: + 256, 512 and the E gates
[[nodiscard]] BenchmarkResult run_tvb09(Tier tier);  // fast: one bias (V_FB + 1.5 V); full: the V_G grid
[[nodiscard]] BenchmarkResult run_tvb14(Tier tier);  // fast: 1D tensor; full: + 2D/3D tensor, 2D/3D generated, swept
[[nodiscard]] BenchmarkResult run_tvb11();  // TVR-2 (DECISIONS Unit 55 14)
[[nodiscard]] BenchmarkResult run_tvb12();
// TVR-3 (DECISIONS Unit 55 15): TVB-31, the model-implementation catalogue.
[[nodiscard]] BenchmarkResult run_tvb31_01();  // legacy band set
[[nodiscard]] BenchmarkResult run_tvb31_02();  // Passler gap
[[nodiscard]] BenchmarkResult run_tvb31_03();  // Couderc DOS mass fits
[[nodiscard]] BenchmarkResult run_tvb31_04();  // n_i, intrinsic level, d ln n_i / dT
[[nodiscard]] BenchmarkResult run_tvb31_05();  // Slotboom narrowing and n_ie
[[nodiscard]] BenchmarkResult run_tvb31_06(Tier tier);  // Caughey-Thomas (tiers: grid density, MI-Q9)
[[nodiscard]] BenchmarkResult run_tvb31_07(Tier tier);  // Masetti
[[nodiscard]] BenchmarkResult run_tvb31_08(Tier tier);  // Canali
[[nodiscard]] BenchmarkResult run_tvb31_09(Tier tier);  // Lombardi
[[nodiscard]] BenchmarkResult run_tvb31_10(Tier tier);  // energy-balance mobility
[[nodiscard]] BenchmarkResult run_tvb31_11(Tier tier);  // Scharfetter lifetimes
[[nodiscard]] BenchmarkResult run_tvb31_12(Tier tier);  // SRH
[[nodiscard]] BenchmarkResult run_tvb31_13(Tier tier);  // Auger
[[nodiscard]] BenchmarkResult run_tvb31_14(Tier tier);  // radiative
[[nodiscard]] BenchmarkResult run_tvb31_15(Tier tier);  // trap occupancy and kinetics
[[nodiscard]] BenchmarkResult run_tvb31_16(Tier tier);  // trap-band quadrature
[[nodiscard]] BenchmarkResult run_tvb31_17(Tier tier);  // incomplete ionization
[[nodiscard]] BenchmarkResult run_tvb31_18(Tier tier);  // Altermatt-Schenk levels
[[nodiscard]] BenchmarkResult run_tvb31_19(Tier tier);  // neutral equilibria
[[nodiscard]] BenchmarkResult run_tvb32();

struct KnownFinding {
    std::string_view benchmark;
    std::string_view quantity;           // a comparison's quantity starts with this
    bool (*in_region)(double parameter);  // the ruled region
    std::string_view label;              // "(i)", "(ii)", "(iii)"
    std::string_view ruling;             // where it is recorded
};

// Findings (i)-(iii) of DECISIONS 11.13, (iv)-(v) of 14.9 and 14.11, and the documented departures of 15.4 (F1).
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
