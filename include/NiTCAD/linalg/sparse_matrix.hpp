// Sparse matrix storage owned by linalg (ARCHITECTURE.md 6.10, item 1).
//
// Compressed sparse row (CSR): row_offsets has rows() + 1 entries; the column indices of row r
// are col_indices[row_offsets[r] .. row_offsets[r + 1]), strictly increasing; values runs in
// parallel with col_indices. No backend type appears here: a backend converts this storage
// itself.
//
// The pattern is fixed at construction and the values may be rewritten in place, so a Newton
// loop can refill the same pattern and let the solver reuse its symbolic analysis.
//
// The scalar is double, or std::complex<double> for AC small-signal analysis, which solves
// (J + i omega C) x = b as the legacy code did (decision A5). SparseMatrix is the real matrix.
#pragma once

#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::linalg {

// 32-bit row and column indices: the index type of the Eigen backend and of LP64 PARDISO,
// so neither needs an index conversion per factorization. It bounds the number of nonzeros of
// A at 2^31 - 1, far beyond what a direct factorization can handle (decision A5).
using Index = std::int32_t;

template <class T>
concept MatrixScalar = std::same_as<T, double> || std::same_as<T, std::complex<double>>;

template <MatrixScalar Scalar>
struct BasicTriplet {
    Index row;
    Index col;
    Scalar value;
};

template <MatrixScalar Scalar>
class BasicSparseMatrix {
public:
    // The empty 0 x 0 matrix.
    BasicSparseMatrix() = default;

    // Builds the CSR matrix from (row, col, value) entries in any order.
    // - Duplicate (row, col) entries are summed, in input order, so the result is deterministic.
    // - Explicit zeros are kept in the pattern: an entry that happens to be zero in one Newton
    //   iteration must not change the pattern.
    // - Values are not checked here; non-finite values are rejected when factorized.
    // Errors: invalid_input for a negative dimension or an index out of range;
    // resource_exhausted if the number of nonzeros does not fit Index.
    [[nodiscard]] static std::expected<BasicSparseMatrix, base::Error> from_triplets(
        Index rows, Index cols, std::span<const BasicTriplet<Scalar>> entries);

    [[nodiscard]] Index rows() const noexcept { return rows_; }
    [[nodiscard]] Index cols() const noexcept { return cols_; }
    [[nodiscard]] std::size_t nonzeros() const noexcept { return values_.size(); }

    [[nodiscard]] std::span<const Index> row_offsets() const noexcept { return row_offsets_; }
    [[nodiscard]] std::span<const Index> col_indices() const noexcept { return col_indices_; }
    [[nodiscard]] std::span<const Scalar> values() const noexcept { return values_; }
    // Rewrites values in place; the pattern cannot change through this view.
    [[nodiscard]] std::span<Scalar> values() noexcept { return values_; }

    // True if both matrices have the same dimensions, row offsets and column indices.
    [[nodiscard]] bool has_same_pattern(const BasicSparseMatrix& other) const noexcept;

    // y = A x. Requires x.size() == cols(), y.size() == rows(), and x and y not overlapping.
    void multiply(std::span<const Scalar> x, std::span<Scalar> y) const;

private:
    Index rows_ = 0;
    Index cols_ = 0;
    std::vector<Index> row_offsets_ = std::vector<Index>(1, 0);
    std::vector<Index> col_indices_;
    std::vector<Scalar> values_;
};

// Both instantiations are compiled once, in sparse_matrix.cpp.
extern template class BasicSparseMatrix<double>;
extern template class BasicSparseMatrix<std::complex<double>>;

using Triplet = BasicTriplet<double>;
using SparseMatrix = BasicSparseMatrix<double>;
using ComplexTriplet = BasicTriplet<std::complex<double>>;
using ComplexSparseMatrix = BasicSparseMatrix<std::complex<double>>;

}  // namespace NiTCAD::linalg
