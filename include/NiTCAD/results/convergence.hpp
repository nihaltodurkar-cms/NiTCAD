// Newton convergence history of one solve (ARCHITECTURE.md 6.6): plain data.
#pragma once

#include <cstddef>
#include <optional>
#include <vector>

namespace NiTCAD::results {

struct IterationRecord {
    int iteration;    // 1-based
    double update;    // largest scaled correction (the convergence measure, before damping)
    double residual;  // largest |F| at the iterate the correction was computed from
    // The fraction of the damped step taken (solve::LineSearch, Unit 30); 1 without a search.
    double step = 1.0;
};

struct ConvergenceRecord {
    std::vector<IterationRecord> iterations;
    bool converged = false;
    // Smallest smallest-to-largest pivot ratio over the factorizations, when the linear solver
    // reports one (ARCHITECTURE.md 6.10).
    std::optional<double> smallest_pivot_ratio;
    // The last factorization's matrix nonzeros and factor fill nnz(L + U), when the solver reports
    // it (Unit 17 ID5's cost model). Diagnostics of the run: not saved, not observed.
    std::size_t matrix_nonzeros = 0;
    std::optional<std::size_t> factor_fill;
};

}  // namespace NiTCAD::results
