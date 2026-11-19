// A run: its record and the bias points it produced (ARCHITECTURE.md 6.6, 6.9): plain data.
//
// A run that is cancelled or fails part-way keeps every completed point; `stopped` then says why
// (cancelled, non_convergence, ...), and `unfinished` holds the convergence history of the point
// that was being solved, if one was. A run that finished has neither.
#pragma once

#include <cstddef>
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

    bool operator==(const Provenance&) const = default;
};

struct RunRecord {
    // A 64-bit FNV-1a digest of the device description and the solve options. Equal inputs give
    // equal digests on the same build; it identifies a run, it is not a checksum to compare across
    // machines (6.8).
    std::uint64_t input_identity = 0;
    // The options in effect, as (name, value) pairs, e.g. ("newton.tol_update", 1e-8).
    std::vector<std::pair<std::string, double>> settings;
    Provenance provenance;
    // The device mesh's mesh::graph_digest (Unit 45, RS3): pairs a result with its mesh exactly.
    // Not part of input_identity (which covers the mesh already); 0 when unknown (a file written
    // before Unit 45).
    std::uint64_t mesh_digest = 0;
    // For a mesh generated from a polyhedral complex (Unit 53, FD14): the deck's mesh type (e.g.
    // "extrude3d") and mesh::generation_digest of the complex and options. Empty and 0 otherwise;
    // not part of input_identity.
    std::string mesh_type;
    std::uint64_t mesh_generation_digest = 0;
};

// A trace's resume state (Unit 61, RD5): what solve::trace_bias needs to continue after its last
// accepted point exactly as the uninterrupted trace would, every double bit for bit. Taken at the
// acceptance of `point`; a trace stopped later (cancelled, or a failed attempt) resumes from it and
// replays those attempts.
struct TraceResume {
    std::vector<double> state;    // the scaled state, lambda (bias / V_T) last
    std::vector<double> tangent;  // the unit tangent there, in the arc metric
    double step = 0.0;            // the next arc step ds
    std::size_t entry = 0;        // the entry kept for the arc row (the pivot)
    std::size_t attempt = 0;      // attempts so far (progress and error indices)
    std::size_t points = 0;       // accepted points so far, the starting point included
    double current_floor_A = 0.0; // I0 of the arc metric, from the starting point
    BiasPoint point;              // the last accepted point (a restarted trace's first)
    // The trace's identity without its initial state (as TransientResume::identity).
    std::uint64_t identity = 0;
    // The trace ended here (landed on its end, reached its current limit or max_points): nothing
    // is left to continue.
    bool finished = false;
};

struct Sweep {
    RunRecord run;
    std::vector<BiasPoint> points;
    std::optional<base::Error> stopped;
    std::optional<ConvergenceRecord> unfinished;
    // A trace's resume state at its last accepted point (Unit 61, RD5); empty for a sweep.
    std::optional<TraceResume> resume;
};

}  // namespace NiTCAD::results
