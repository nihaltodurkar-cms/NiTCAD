// Newton convergence history of one solve (ARCHITECTURE.md 6.6): plain data.
#pragma once

#include <optional>
#include <vector>

namespace NiTCAD::results {

struct IterationRecord {
    int iteration;    // 1-based
    double update;    // largest scaled correction (the convergence measure, before damping)
    double residual;  // largest |F| at the iterate the correction was computed from
};

struct ConvergenceRecord {
    std::vector<IterationRecord> iterations;
    bool converged = false;
    // Smallest smallest-to-largest pivot ratio over the factorizations, when the linear solver
    // reports one (ARCHITECTURE.md 6.10).
    std::optional<double> smallest_pivot_ratio;
};

}  // namespace NiTCAD::results
