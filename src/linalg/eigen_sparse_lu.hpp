// Eigen SparseLU backend (D3). Private to linalg: this header includes Eigen and is never
// included by a public header, so Eigen stays out of the interface (6.10, Q5).
#pragma once

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <complex>
#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "backend.hpp"

namespace NiTCAD::linalg::detail {

// Scalar is double or std::complex<double>; both are instantiated in eigen_sparse_lu.cpp.
template <MatrixScalar Scalar>
class SparseLuBackend final : public Backend<Scalar> {
public:
    // Builds the column-compressed copy of A's pattern and runs the symbolic analysis.
    // Requires a square, non-empty matrix.
    [[nodiscard]] std::expected<void, base::Error> analyze(
        const BasicSparseMatrix<Scalar>& a) override;
    // Numeric factorization of values laid out in the CSR order of the analyzed pattern; its pivot
    // ratio and the original column of the smallest pivot.
    [[nodiscard]] std::expected<Factorization, base::Error> factorize(
        std::span<const Scalar> csr_values) override;
    // x = A^-1 b with the last factorization. b and x have the analyzed dimension.
    void solve(std::span<const Scalar> b, std::span<Scalar> x) override;

private:
    // COLAMD, not AMD on A^T + A: legacy measured AMD as pathological with SparseLU on the coupled
    // psi/n/p Jacobians (more than 15 s where COLAMD took about 45 ms, 11,640 unknowns).
    using Matrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, Index>;
    using Lu = Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<Index>>;

    Matrix csc_;
    std::vector<Index> csr_to_csc_;      // position in csc_ of each CSR entry
    std::vector<Index> factored_to_original_;  // original column of each factored column
    std::optional<Lu> lu_;               // recreated by every analysis, so no state outlives a failure
};

}  // namespace NiTCAD::linalg::detail
