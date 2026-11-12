// The interface every direct backend implements behind BasicLinearSolver (6.10). Private to
// linalg: backend types never appear in a public header.
#pragma once

#include <cstddef>
#include <expected>
#include <optional>
#include <span>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace NiTCAD::linalg::detail {

// Smallest |pivot| over largest |pivot| of a factorization, and the original column of the
// smallest one when the backend can tell it.
struct PivotRatio {
    double ratio;
    std::optional<std::size_t> column;
};

struct Factorization {
    std::optional<PivotRatio> pivots;              // when the backend reports its pivots
    std::optional<std::size_t> perturbed_pivots;   // when the backend perturbs tiny pivots
};

template <MatrixScalar Scalar>
class Backend {
public:
    virtual ~Backend() = default;
    // Records A's pattern for the analysis. Requires a square, non-empty matrix.
    [[nodiscard]] virtual std::expected<void, base::Error> analyze(
        const BasicSparseMatrix<Scalar>& a) = 0;
    // Numeric factorization of values laid out in the CSR order of the analyzed pattern.
    [[nodiscard]] virtual std::expected<Factorization, base::Error> factorize(
        std::span<const Scalar> csr_values) = 0;
    // x = A^-1 b with the last factorization. b and x have the analyzed dimension.
    virtual void solve(std::span<const Scalar> b, std::span<Scalar> x) = 0;
    // Whether solve_transposed is available (Unit 17 AQ4: the direct backends; an incomplete
    // factorization is a preconditioner, not a transposable solve).
    [[nodiscard]] virtual bool can_transpose() const noexcept { return false; }
    // x = A^-T b with the last factorization (the transpose, not the conjugate transpose, for a
    // complex A). Requires can_transpose(); b and x have the analyzed dimension.
    virtual void solve_transposed(std::span<const Scalar> /*b*/, std::span<Scalar> /*x*/) {}

    // The backend's own work, counted where it is done (BackendCounts).
    [[nodiscard]] const BackendCounts& counts() const noexcept { return counts_; }

protected:
    BackendCounts counts_;
};

}  // namespace NiTCAD::linalg::detail
