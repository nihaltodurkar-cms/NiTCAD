// Cancellation and progress for the solve layer (ARCHITECTURE.md 6.9).
//
// Cancellation: the stop token is checked at safe points only, before every Newton iteration and
// between bias points. A factorization or linear solve in progress is not interrupted; that latency
// is a stated limit. A cancelled solve returns ErrorCode::cancelled; a cancelled sweep keeps its
// completed points (results/run.hpp).
// Progress: the callback is called on the solving thread after every Newton iteration. Events are
// strictly increasing in (phase, point, iteration). The callback must not throw and must not call
// back into the solve; an application forwards events to its UI thread itself (by posting a
// message), the solver never calls UI code.
#pragma once

#include <cstddef>
#include <functional>
#include <stop_token>

namespace NiTCAD::solve {

enum class Phase : int { equilibrium = 0, bias = 1 };

struct Progress {
    Phase phase;
    std::size_t point;        // bias point index (0 in the equilibrium phase)
    std::size_t point_count;  // bias points in the run (1 in the equilibrium phase)
    int iteration;            // 1-based within the point
    double update;            // largest scaled correction of this iteration
    double residual;          // largest |F| before the correction
    bool converged;           // true on the iteration that converged
};

struct RunControl {
    std::stop_token stop;
    std::function<void(const Progress&)> progress;
};

}  // namespace NiTCAD::solve
