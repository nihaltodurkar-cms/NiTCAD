// Eigen SparseLU backend (D3). Private to linalg: this header includes Eigen and is never
// included by a public header, so Eigen stays out of the interface (6.10, Q5).
#pragma once

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace NiTCAD::linalg::detail {

// Smallest |pivot| over largest |pivot| of a factorization, and the original column of the
// smallest one.
struct PivotRatio {
    double ratio;
    std::size_t column;
};

class SparseLuBackend {
public:
    // Builds the column-compressed copy of A's pattern and runs the symbolic analysis.
    // Requires a square, non-empty matrix.
    [[nodiscard]] std::expected<void, base::Error> analyze(const SparseMatrix& a);
    // Numeric factorization of values laid out in the CSR order of the analyzed pattern.
    [[nodiscard]] std::expected<PivotRatio, base::Error> factorize(std::span<const double> csr_values);
    // x = A^-1 b with the last factorization. b and x have the analyzed dimension.
    void solve(std::span<const double> b, std::span<double> x);

private:
    // COLAMD, not AMD on A^T + A: legacy measured AMD as pathological with SparseLU on the coupled
    // psi/n/p Jacobians (more than 15 s where COLAMD took about 45 ms, 11,640 unknowns).
    using Matrix = Eigen::SparseMatrix<double, Eigen::ColMajor, Index>;
    using Lu = Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<Index>>;

    Matrix csc_;
    std::vector<Index> csr_to_csc_;      // position in csc_ of each CSR entry
    std::vector<Index> factored_to_original_;  // original column of each factored column
    std::optional<Lu> lu_;               // recreated by every analysis, so no state outlives a failure
};

}  // namespace NiTCAD::linalg::detail
