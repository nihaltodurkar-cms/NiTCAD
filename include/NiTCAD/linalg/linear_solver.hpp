// Backend-neutral sparse linear solver (ARCHITECTURE.md 6.10, D3).
//
// Three separate steps, so the symbolic analysis can be reused across Newton iterations:
//   analyze(A)    pattern only: ordering and elimination structure;
//   factorize(A)  numeric factorization; if A's pattern differs from the analyzed one (or
//                 nothing is analyzed yet) it analyzes first, so a moving pattern stays correct;
//   solve(b, x)   solves with the last factorization and checks the result.
//
// Acceptance is decided here, not by the backend: a solve succeeds only if x is finite and
// ||A x - b||_2 / ||b||_2 <= max_relative_residual (the absolute residual when b = 0), even if
// the backend reported success.
//
// Failures are values (6.7): invalid_input (non-finite A or b, empty system, bad
// configuration), singular_system, inaccurate_solve (residual above the threshold or non-finite
// x), resource_exhausted (out of memory). After a failed analyze or factorize nothing is kept:
// the next factorize analyzes afresh. Wrong dimensions or overlapping b and x are programmer
// errors and go through NITCAD_EXPECTS.
//
// Backend types do not appear in this header.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace NiTCAD::linalg {

enum class SolverBackend : std::uint8_t {
    // Eigen SparseLU with COLAMD column ordering. Single-threaded.
    eigen_sparse_lu,
};

struct SolverConfig {
    SolverBackend backend = SolverBackend::eigen_sparse_lu;
    // Thread count (6.8: default 1). eigen_sparse_lu accepts only 1.
    int threads = 1;
    // Largest accepted relative residual. Legacy value for its direct solver (PARDISO): 1e-6.
    double max_relative_residual = 1e-6;
};

struct SolveReport {
    double relative_residual;
    // Iterations used, for iterative backends only; empty for a direct backend.
    std::optional<std::size_t> iterations;
};

namespace detail {
class SparseLuBackend;
}

class LinearSolver {
public:
    // Errors: invalid_input if threads is not supported by the backend, or
    // max_relative_residual is not a positive finite number.
    [[nodiscard]] static std::expected<LinearSolver, base::Error> create(const SolverConfig& config);

    LinearSolver(LinearSolver&&) noexcept;
    LinearSolver& operator=(LinearSolver&&) noexcept;
    LinearSolver(const LinearSolver&) = delete;
    LinearSolver& operator=(const LinearSolver&) = delete;
    ~LinearSolver();

    // Requires a square matrix.
    [[nodiscard]] std::expected<void, base::Error> analyze(const SparseMatrix& a);
    // Requires a square matrix.
    [[nodiscard]] std::expected<void, base::Error> factorize(const SparseMatrix& a);
    // Requires is_factorized(), b.size() == x.size() == the factorized dimension, and b and x not
    // overlapping. On failure x holds no meaningful result.
    [[nodiscard]] std::expected<SolveReport, base::Error> solve(std::span<const double> b,
                                                                std::span<double> x);

    [[nodiscard]] bool is_analyzed() const noexcept { return analyzed_; }
    [[nodiscard]] bool is_factorized() const noexcept { return factorized_; }
    // Successful analyses and factorizations since creation (to observe symbolic reuse).
    [[nodiscard]] std::size_t analyses() const noexcept { return analyses_; }
    [[nodiscard]] std::size_t factorizations() const noexcept { return factorizations_; }
    [[nodiscard]] const SolverConfig& config() const noexcept { return config_; }

private:
    explicit LinearSolver(const SolverConfig& config);

    [[nodiscard]] std::expected<void, base::Error> analyze_unchecked(const SparseMatrix& a);
    void reset() noexcept;

    SolverConfig config_;
    std::unique_ptr<detail::SparseLuBackend> backend_;
    SparseMatrix matrix_;            // pattern of the analysis; values of the factorization
    std::vector<double> residual_;   // work vector for the acceptance check
    bool analyzed_ = false;
    bool factorized_ = false;
    std::size_t analyses_ = 0;
    std::size_t factorizations_ = 0;
};

}  // namespace NiTCAD::linalg
