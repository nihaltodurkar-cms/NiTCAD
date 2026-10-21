// Eigen IncompleteLUT as the preconditioner of the gmres_ilut backend (Unit 44, LS1). Private to
// linalg: this header includes Eigen and is never included by a public header (6.10, Q5).
//
// In the Backend interface "factorize" builds the incomplete factors and "solve" applies them,
// x = M^-1 b, an approximation of A^-1 b: BasicLinearSolver runs GMRES around it and never takes
// its result as a solution. Counts: one analysis per analyze(), one factorization per ILUT build,
// one solve per application of M^-1.
#pragma once

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include <complex>
#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "backend.hpp"

namespace NiTCAD::linalg::detail {

// Scalar is double or std::complex<double>; both are instantiated in eigen_ilut.cpp.
template <MatrixScalar Scalar>
class IlutBackend final : public Backend<Scalar> {
public:
    IlutBackend(double drop_tolerance, int fill_factor);

    // Copies A's pattern and computes ILUT's fill-reducing ordering (AMD on the pattern of
    // A + A^T). Requires a square, non-empty matrix.
    [[nodiscard]] std::expected<void, base::Error> analyze(
        const BasicSparseMatrix<Scalar>& a) override;
    // Builds the incomplete factors from values in the CSR order of the analyzed pattern. ILUT does
    // not pivot: a zero pivot is replaced by sqrt(drop tolerance) times its row's 2-norm (Eigen's
    // rule), and a factor that is not finite is inaccurate_solve. Reports no pivots: the
    // incomplete factor's diagonal is no singularity measure.
    [[nodiscard]] std::expected<Factorization, base::Error> factorize(
        std::span<const Scalar> csr_values) override;
    // x = M^-1 b with the last build.
    void solve(std::span<const Scalar> b, std::span<Scalar> x) override;

private:
    using Matrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>;
    // Eigen keeps L and U in one row-major matrix; matrixL() and matrixU() return copies, so the
    // factors are read in place for the finiteness check (Eigen 5.0.1, see eigen_ilut.cpp).
    struct Ilut : Eigen::IncompleteLUT<Scalar, Index> {
        [[nodiscard]] const Matrix& factors() const noexcept { return this->m_lu; }
    };

    double drop_tolerance_;
    int fill_factor_;
    Matrix csr_;
    std::optional<Ilut> ilut_;  // recreated by every analysis, so no state outlives a failure
    bool built_ = false;
};

}  // namespace NiTCAD::linalg::detail
