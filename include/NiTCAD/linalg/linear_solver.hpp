// Backend-neutral sparse linear solver (ARCHITECTURE.md 6.10, D3).
//
// Three separate steps, so the symbolic analysis can be reused across Newton iterations:
//   analyze(A)    pattern only: ordering and elimination structure;
//   factorize(A)  numeric factorization; if A's pattern differs from the analyzed one (or
//                 nothing is analyzed yet) it analyzes first, so a moving pattern stays correct;
//   solve(b, x)   solves with the last factorization and checks the result.
// With the Eigen backend the reusable part is only the column ordering and elimination tree:
// every factorization redoes its own symbolic work because row pivoting depends on the values
// (measured: the ordering is 2-6% of a 2D factorization and 0.1% of a 3D one). With PARDISO the
// whole analysis (matching, scaling, nested-dissection ordering, symbolic factorization) is done
// once per pattern, at the first factorization after it (the matching and scaling need values),
// and every later factorization of that pattern is numeric only (backend_counts()).
//
// What this layer adds around the backend, identically for every backend:
// - Equilibration (optional, on by default): rows, then columns, are scaled by powers of two so
//   their largest entry lies in [0.5, 1). Exact in floating point; it changes only pivoting.
// - Singularity checks: a row or column without a nonzero value, an exactly zero pivot, and a
//   smallest-to-largest pivot ratio below min_pivot_ratio are singular_system. The last catches a
//   floating region (a connected block without a Dirichlet row), which leaves one pivot at
//   rounding level instead of exactly zero. It is a heuristic: a valid region anchored only
//   weakly (see min_pivot_ratio) is indistinguishable from a floating one and is rejected too.
//   A non-finite pivot (overflow during elimination) is inaccurate_solve.
// - Acceptance: a solve succeeds only if x is finite and its componentwise backward error, on the
//   original unscaled system, is at most max_backward_error, even if the backend reported success.
//   The componentwise backward error is max_i |b - A x|_i / (|A| |x| + |b|)_i (Oettli-Prager),
//   with the Arioli-Demmel-Duff denominator for rows where that one is at rounding level. It is
//   independent of row scaling, unlike ||b - A x|| / ||b||, which is reported for information.
//   Up to max_refinement_steps steps of iterative refinement are tried before a solve is rejected;
//   refinement stops at the first step that does not reduce the backward error, and the better
//   iterate is kept.
//
// Failures are values (6.7): invalid_input (non-finite A or b, empty system, bad configuration),
// singular_system, inaccurate_solve (backward error above the limit or non-finite x),
// resource_exhausted (out of memory). After a failed analyze or factorize nothing is kept: the
// next factorize analyzes afresh. Wrong dimensions, overlapping b and x, and using a moved-from
// solver are programmer errors and go through NITCAD_EXPECTS.
//
// A LinearSolver is not safe to use from two threads at once; separate solvers are independent.
// Backend types do not appear in this header.
#pragma once

#include <complex>
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
    // Eigen SparseLU with COLAMD column ordering. Single-threaded. The default.
    eigen_sparse_lu,
    // Intel MKL PARDISO (Unit 18), opt-in: real and complex unsymmetric matrices, METIS nested
    // dissection, PARDISO's weighted matching and scaling, no iterative refinement of its own (this
    // layer's applies). It does not fail on a tiny pivot: it replaces a pivot below 1e-13 of the
    // largest by that size and counts it (FactorizationReport::perturbed_pivots). Requires the MKL
    // runtime to be loaded first (mkl_runtime.hpp). Results agree with eigen_sparse_lu within the
    // acceptance tolerances, not bit for bit, and may differ at rounding level between thread
    // counts (6.8).
    mkl_pardiso,
};

struct SolverConfig {
    SolverBackend backend = SolverBackend::eigen_sparse_lu;
    // Thread count (6.8: default 1). eigen_sparse_lu accepts only 1; mkl_pardiso 1 to 1024 (MKL
    // threads for this solver's calls only, set thread-locally around each one).
    int threads = 1;
    // Scale rows and columns by powers of two before factorizing.
    bool equilibrate = true;
    // Largest accepted componentwise backward error. A backward-stable solve gives a few
    // multiples of 1.1e-16; a breakdown gives 1e-3 or more.
    double max_backward_error = 1e-8;
    // Iterative-refinement steps tried when a solve is above max_backward_error.
    int max_refinement_steps = 1;
    // A factorization with smallest |pivot| / largest |pivot| below this is singular_system;
    // 0 disables the check, and is required when equilibrate is false (raw pivots of a badly
    // scaled matrix span many decades). Measured on equilibrated systems:
    // - floating regions (singular): 1e-16 to 2.3e-14, rising slowly with size;
    // - well-anchored valid systems (graded meshes up to 1e8, scaled coupled systems): 2.7e-8 up;
    // - weakly anchored valid systems overlap the singular range: a region tied to the rest only
    //   by a zeroth-order term delta (relative to its couplings g) gives about n delta / (2 g), and
    //   one joined by a single weak link w gives about w / 2. Such a region with n = 1000 is
    //   rejected for delta / g = 1e-14 or w = 1e-12 (an SRH-only floating body on a fine mesh can
    //   reach this). No threshold separates the two; topology checks belong to the device layer.
    // mkl_pardiso reports no pivot ratio: its factor's diagonal is not a singularity measure (it
    // pivots statically; a valid equilibrated system measured 1.8e-13 with no perturbed pivot).
    // Instead, when this is positive, a perturbed pivot (below 1e-13 of the largest) is
    // singular_system: a floating region gives one. Between 1e-13 and min_pivot_ratio PARDISO
    // cannot tell a floating region from a valid system; the acceptance check still applies.
    double min_pivot_ratio = 1e-11;
};

struct FactorizationReport {
    // Smallest |pivot| / largest |pivot| (of the equilibrated matrix when equilibrating), if the
    // backend reports pivots.
    std::optional<double> pivot_ratio;
    // Pivots the backend perturbed instead of failing (mkl_pardiso); empty for eigen_sparse_lu.
    std::optional<std::size_t> perturbed_pivots;
};

// The backend's own work since the solver was created: symbolic analyses (an ordering and symbolic
// factorization), numeric factorizations, and triangular solves (one per solve() call and one per
// refinement step). With mkl_pardiso these are its phases 11, 22 and 33, counted where they are
// called; with eigen_sparse_lu its analyzePattern, factorize and solve.
struct BackendCounts {
    std::size_t analyses = 0;
    std::size_t factorizations = 0;
    std::size_t solves = 0;
};

struct SolveReport {
    // Componentwise backward error on the original system (the acceptance measure).
    double backward_error;
    // ||b - A x||_inf / ||b||_inf, or ||b - A x||_inf when b = 0. For information only.
    double relative_residual;
    // Iterative-refinement steps kept (a step that did not reduce the backward error is undone).
    std::size_t refinement_steps;
    // Iterations used, for iterative backends only; empty for a direct backend.
    std::optional<std::size_t> iterations;
};

namespace detail {
template <MatrixScalar Scalar>
class Backend;
}

// The solver of A x = b for a real (LinearSolver) or complex (ComplexLinearSolver) matrix. The
// complex solver serves AC small-signal analysis, (J + i omega C) x = b (decision A5); everything
// above applies to it with |.| the complex modulus: the equilibration scales are powers of two from
// the moduli of the entries, the pivot ratio is of the pivots' moduli, and the backward error is
// Oettli-Prager's with moduli. A value is finite when both its parts are.
template <MatrixScalar Scalar>
class BasicLinearSolver {
public:
    using Matrix = BasicSparseMatrix<Scalar>;

    // Errors: invalid_input if the backend is unknown, or mkl_pardiso without the MKL runtime
    // loaded, threads is not supported by the backend, max_backward_error is not a positive finite
    // number, max_refinement_steps is negative, min_pivot_ratio is not in [0, 1), or
    // min_pivot_ratio is nonzero while equilibrate is false.
    [[nodiscard]] static std::expected<BasicLinearSolver, base::Error> create(
        const SolverConfig& config);

    // The moved-from solver is left neither analyzed nor factorized, and must not be used again.
    BasicLinearSolver(BasicLinearSolver&& other) noexcept;
    BasicLinearSolver& operator=(BasicLinearSolver&& other) noexcept;
    BasicLinearSolver(const BasicLinearSolver&) = delete;
    BasicLinearSolver& operator=(const BasicLinearSolver&) = delete;
    ~BasicLinearSolver();

    // Requires a square matrix.
    [[nodiscard]] std::expected<void, base::Error> analyze(const Matrix& a);
    // Requires a square matrix.
    [[nodiscard]] std::expected<FactorizationReport, base::Error> factorize(const Matrix& a);
    // Requires is_factorized(), b.size() == x.size() == the factorized dimension, and b and x not
    // overlapping. On failure x holds no meaningful result.
    [[nodiscard]] std::expected<SolveReport, base::Error> solve(std::span<const Scalar> b,
                                                                std::span<Scalar> x);

    [[nodiscard]] bool is_analyzed() const noexcept { return analyzed_; }
    [[nodiscard]] bool is_factorized() const noexcept { return factorized_; }
    // Successful analyses and factorizations since creation (to observe symbolic reuse).
    [[nodiscard]] std::size_t analyses() const noexcept { return analyses_; }
    [[nodiscard]] std::size_t factorizations() const noexcept { return factorizations_; }
    // The backend's own analyses, factorizations and triangular solves (BackendCounts).
    [[nodiscard]] const BackendCounts& backend_counts() const noexcept;
    [[nodiscard]] const SolverConfig& config() const noexcept { return config_; }

private:
    explicit BasicLinearSolver(const SolverConfig& config);

    [[nodiscard]] std::expected<void, base::Error> analyze_unchecked(const Matrix& a);
    [[nodiscard]] std::expected<FactorizationReport, base::Error> factorize_unchecked(
        const Matrix& a);
    // Solves the equilibrated system for right-hand side rhs into x (unscaled).
    void scaled_solve(std::span<const Scalar> rhs, std::span<Scalar> x);
    // Componentwise backward error of x on the original system; fills residual_ with b - A x.
    [[nodiscard]] double backward_error(std::span<const Scalar> b, std::span<const Scalar> x);
    void reset() noexcept;

    SolverConfig config_;
    std::unique_ptr<detail::Backend<Scalar>> backend_;
    Matrix matrix_;                     // pattern of the analysis; original values of the factorization
    std::vector<Scalar> scaled_values_; // values handed to the backend
    // Power-of-two equilibration scales. Always computed (the backward error uses the column
    // scales); applied to the factorization only when config_.equilibrate.
    std::vector<double> row_scale_;
    std::vector<double> col_scale_;
    std::vector<double> weighted_row_max_;  // max_j |a_ij| col_scale_j, for the backward error
    std::vector<Scalar> residual_;      // work vectors for solve()
    std::vector<Scalar> work_rhs_;
    std::vector<Scalar> work_x_;
    std::vector<Scalar> correction_;
    std::vector<Scalar> previous_x_;
    bool analyzed_ = false;
    bool factorized_ = false;
    std::size_t analyses_ = 0;
    std::size_t factorizations_ = 0;
};

// Both instantiations are compiled once, in linear_solver.cpp.
extern template class BasicLinearSolver<double>;
extern template class BasicLinearSolver<std::complex<double>>;

using LinearSolver = BasicLinearSolver<double>;
using ComplexLinearSolver = BasicLinearSolver<std::complex<double>>;

}  // namespace NiTCAD::linalg
