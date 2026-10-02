#include "NiTCAD/linalg/sparse_matrix.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <utility>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::linalg {

namespace {

base::Error invalid(std::string message, std::size_t entry) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = entry, .value = std::nullopt}};
}

template <class T>
bool overlaps(std::span<const T> a, std::span<const T> b) noexcept {
    if (a.empty() || b.empty()) {
        return false;
    }
    const std::less<const T*> before;
    return before(a.data(), b.data() + b.size()) && before(b.data(), a.data() + a.size());
}

}  // namespace

template <MatrixScalar Scalar>
std::expected<BasicSparseMatrix<Scalar>, base::Error> BasicSparseMatrix<Scalar>::from_triplets(
    Index rows, Index cols, std::span<const BasicTriplet<Scalar>> entries) {
    if (rows < 0 || cols < 0) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           "sparse matrix dimensions must be non-negative",
                                           std::nullopt});
    }
    for (std::size_t k = 0; k < entries.size(); ++k) {
        const BasicTriplet<Scalar>& t = entries[k];
        if (t.row < 0 || t.row >= rows) {
            return std::unexpected(invalid("sparse matrix entry row index out of range", k));
        }
        if (t.col < 0 || t.col >= cols) {
            return std::unexpected(invalid("sparse matrix entry column index out of range", k));
        }
    }

    const auto n_rows = static_cast<std::size_t>(rows);

    // Stable counting sort by row: within a row, entries keep their input order.
    std::vector<std::size_t> start(n_rows + 1, 0);
    for (const BasicTriplet<Scalar>& t : entries) {
        ++start[static_cast<std::size_t>(t.row) + 1];
    }
    for (std::size_t r = 0; r < n_rows; ++r) {
        start[r + 1] += start[r];
    }
    std::vector<std::size_t> order(entries.size());
    {
        std::vector<std::size_t> next(start.begin(), start.end() - 1);
        for (std::size_t k = 0; k < entries.size(); ++k) {
            order[next[static_cast<std::size_t>(entries[k].row)]++] = k;
        }
    }

    // Per row: stable sort by column, then sum duplicates in input order.
    BasicSparseMatrix m;
    m.rows_ = rows;
    m.cols_ = cols;
    m.row_offsets_.assign(n_rows + 1, 0);
    m.col_indices_.reserve(entries.size());
    m.values_.reserve(entries.size());
    constexpr auto max_nonzeros = static_cast<std::size_t>(std::numeric_limits<Index>::max());
    for (std::size_t r = 0; r < n_rows; ++r) {
        const auto first = order.begin() + static_cast<std::ptrdiff_t>(start[r]);
        const auto last = order.begin() + static_cast<std::ptrdiff_t>(start[r + 1]);
        std::stable_sort(first, last, [&](std::size_t a, std::size_t b) {
            return entries[a].col < entries[b].col;
        });
        for (auto it = first; it != last; ++it) {
            const BasicTriplet<Scalar>& t = entries[*it];
            if (it != first && entries[*(it - 1)].col == t.col) {
                m.values_.back() += t.value;
            } else {
                m.col_indices_.push_back(t.col);
                m.values_.push_back(t.value);
            }
        }
        if (m.values_.size() > max_nonzeros) {
            return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                               "sparse matrix has too many nonzeros for Index",
                                               std::nullopt});
        }
        m.row_offsets_[r + 1] = static_cast<Index>(m.values_.size());
    }
    m.col_indices_.shrink_to_fit();
    m.values_.shrink_to_fit();
    return m;
}

template <MatrixScalar Scalar>
bool BasicSparseMatrix<Scalar>::has_same_pattern(const BasicSparseMatrix& other) const noexcept {
    return rows_ == other.rows_ && cols_ == other.cols_ && row_offsets_ == other.row_offsets_ &&
           col_indices_ == other.col_indices_;
}

template <MatrixScalar Scalar>
void BasicSparseMatrix<Scalar>::multiply(std::span<const Scalar> x, std::span<Scalar> y) const {
    NITCAD_EXPECTS(x.size() == static_cast<std::size_t>(cols_));
    NITCAD_EXPECTS(y.size() == static_cast<std::size_t>(rows_));
    NITCAD_EXPECTS(!overlaps<Scalar>(x, y));
    for (std::size_t r = 0; r < y.size(); ++r) {
        Scalar sum{};
        const auto end = static_cast<std::size_t>(row_offsets_[r + 1]);
        for (auto k = static_cast<std::size_t>(row_offsets_[r]); k < end; ++k) {
            sum += values_[k] * x[static_cast<std::size_t>(col_indices_[k])];
        }
        y[r] = sum;
    }
}

template class BasicSparseMatrix<double>;
template class BasicSparseMatrix<std::complex<double>>;

}  // namespace NiTCAD::linalg
