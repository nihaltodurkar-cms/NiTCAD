// The verification report on disk (Unit 55, TVR-1, DECISIONS 11.7; owner ruling RC-Q3, the one io addition of TVR-1).
//
// Plain data, so io knows nothing of the verify layer (which fills it): per benchmark, the reference's identity and
// grades, each compared quantity with q_sim, q_ref, U_ref, the U_num components ("not estimated" in TVR-1), the
// metric, the criterion, the outcome and its reason. Written in io's container (magic, format version, content kind
// verify_report, length, keyed payload, FNV-1a 64), atomically, doubles by bit pattern, so a loaded report equals the
// saved one field for field. Two runs on one machine and thread count give equal reports except
// provenance.started_utc and provenance.wall_s (RC-Q4).
#pragma once

#include <cstdint>
#include <expected>
#include <filesystem>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/run.hpp"

namespace NiTCAD::io {

struct VerifyComparison {
    std::string quantity;
    std::string unit;
    double parameter = 0.0;  // the swept argument (eta, x, T), or 0
    double q_sim = 0.0;
    double q_ref = 0.0;
    double u_ref = 0.0;
    std::string u_num;  // the U_num components, or "not estimated"
    double metric = 0.0;
    std::string criterion;
    std::string outcome;  // "pass", "finding", "fail", "not_applicable", "not_evaluated"
    std::string reason;
    bool operator==(const VerifyComparison&) const = default;
};

struct VerifyBenchmark {
    std::string id;
    std::string title;
    std::string reference;  // identity and grades: class, independence, provenance, citation, version
    std::string outcome;    // the benchmark's overall outcome (or "unverified": every comparison not_applicable)
    std::string run_error;  // empty unless the benchmark could not run
    std::vector<std::string> notes;
    std::vector<VerifyComparison> comparisons;
    bool operator==(const VerifyBenchmark&) const = default;
};

struct VerifyReport {
    results::Provenance provenance;
    std::uint32_t threads = 0;  // the thread count of the run (RC-Q4)
    std::string tier;
    std::vector<std::string> reference_evidence;  // the build-time L0 evidence named, not re-run (RC-Q1)
    std::vector<VerifyBenchmark> benchmarks;
    bool operator==(const VerifyReport&) const = default;
};

[[nodiscard]] std::expected<void, base::Error> save_verify_report(const std::filesystem::path& path,
                                                                  const VerifyReport& report);
[[nodiscard]] std::expected<VerifyReport, base::Error> load_verify_report(const std::filesystem::path& path);

}  // namespace NiTCAD::io
