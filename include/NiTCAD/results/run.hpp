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

struct RunRecord {
    // A 64-bit FNV-1a digest of the device description and the solve options. Equal inputs give
    // equal digests on the same build; it identifies a run, it is not a checksum to compare across
    // machines (6.8).
    std::uint64_t input_identity = 0;
    // The options in effect, as (name, value) pairs, e.g. ("newton.tol_update", 1e-8).
    std::vector<std::pair<std::string, double>> settings;
};

struct Sweep {
    RunRecord run;
    std::vector<BiasPoint> points;
    std::optional<base::Error> stopped;
    std::optional<ConvergenceRecord> unfinished;
};

}  // namespace NiTCAD::results
