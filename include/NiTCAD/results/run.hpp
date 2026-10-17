// A run: its record and the bias points it produced (ARCHITECTURE.md 6.6, 6.9): plain data.
//
// A run that is cancelled or fails part-way keeps every completed point; `stopped` then says why
// (cancelled, non_convergence, ...), and `unfinished` holds the convergence history of the point
// that was being solved, if one was. A run that finished has neither.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/convergence.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::results {

// Which build, on which platform and when, produced a run (Unit 42; DECISIONS.md PV1-PV3). Not
// part of the input identity: two runs of one input on different builds share an identity and
// differ here. Empty (strings empty, numbers 0) for a run read from a file written before Unit 42.
// The linear backend and its thread count are settings (RunRecord::settings), not repeated here.
struct Provenance {
    // The build.
    std::string version;      // NiTCAD's (the CMake project version)
    std::string revision;     // the git commit built; "unknown" outside a git checkout
    bool dirty = false;       // the working tree had uncommitted changes when it was built
    std::string build_type;   // "Release", "Debug"
    std::string compiler;     // e.g. "MSVC 19.51.36260.0"
    std::string eigen_version;
    // The run.
    std::vector<std::int64_t> eigen_cache_bytes;  // L1, L2, L3 Eigen blocks by (fixed, Unit 41)
    std::string mkl_version;  // of the loaded MKL runtime; empty when none is loaded
    std::string os;           // e.g. "Windows 10.0.26300"
    std::string cpu;          // the processor's brand string
    std::uint32_t logical_cpus = 0;
    std::string started_utc;  // ISO 8601, e.g. "2026-10-17T09:12:03.123Z"
    double wall_s = 0.0;      // from the run's start to its result; 0 for a record made alone
};

struct RunRecord {
    // A 64-bit FNV-1a digest of the device description and the solve options. Equal inputs give
    // equal digests on the same build; it identifies a run, it is not a checksum to compare across
    // machines (6.8).
    std::uint64_t input_identity = 0;
    // The options in effect, as (name, value) pairs, e.g. ("newton.tol_update", 1e-8).
    std::vector<std::pair<std::string, double>> settings;
    Provenance provenance;
};

struct Sweep {
    RunRecord run;
    std::vector<BiasPoint> points;
    std::optional<base::Error> stopped;
    std::optional<ConvergenceRecord> unfinished;
};

}  // namespace NiTCAD::results
