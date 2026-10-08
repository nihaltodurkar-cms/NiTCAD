#include "eigen_sparse_lu.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "NiTCAD/base/contract.hpp"

// The failure handling below depends on Eigen 5.0.1 internals: the text of lastErrorMessage(),
// info() being left unset when SparseLU cannot allocate its working memory, the supernodal storage
// read for the pivots, and COLAMD failure being checked only by eigen_assert. A different Eigen
// must be re-checked against these before this assertion is changed.
static_assert(EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 5 && EIGEN_MINOR_VERSION == 0 &&
                  EIGEN_PATCH_VERSION == 1,
              "eigen_sparse_lu.cpp was verified against Eigen 5.0.1 only");

namespace NiTCAD::linalg::detail {

namespace {

// A pivot as an error context value: itself when real, its modulus when complex.
double reported(double v) noexcept { return v; }
double reported(std::complex<double> v) noexcept { return std::abs(v); }

}  // namespace

template <MatrixScalar Scalar>
std::expected<void, base::Error> SparseLuBackend<Scalar>::analyze(
    const BasicSparseMatrix<Scalar>& a) {
    NITCAD_EXPECTS(a.rows() == a.cols() && a.rows() > 0);
    const auto n = static_cast<std::size_t>(a.cols());
    const std::span<const Index> row_offsets = a.row_offsets();
    const std::span<const Index> col_indices = a.col_indices();

    // CSR -> CSC by counting sort on the column. Rows are visited in increasing order, so the
    // row indices within each column come out sorted, as Eigen's compressed storage requires.
    csc_.resize(a.rows(), a.cols());
    csc_.resizeNonZeros(static_cast<Eigen::Index>(a.nonzeros()));
    Index* const outer = csc_.outerIndexPtr();
    Index* const inner = csc_.innerIndexPtr();
    std::fill(outer, outer + n + 1, Index{0});
    for (const Index c : col_indices) {
        ++outer[static_cast<std::size_t>(c) + 1];
    }
    for (std::size_t c = 0; c < n; ++c) {
        outer[c + 1] += outer[c];
    }
    std::vector<Index> next(outer, outer + n);
    csr_to_csc_.resize(a.nonzeros());
    for (std::size_t r = 0; r < n; ++r) {
        const auto end = static_cast<std::size_t>(row_offsets[r + 1]);
        for (auto k = static_cast<std::size_t>(row_offsets[r]); k < end; ++k) {
            const Index position = next[static_cast<std::size_t>(col_indices[k])]++;
            inner[static_cast<std::size_t>(position)] = static_cast<Index>(r);
            csr_to_csc_[k] = position;
        }
    }
    std::fill(csc_.valuePtr(), csc_.valuePtr() + a.nonzeros(), Scalar{});

    // Eigen sizes COLAMD's workspace in Index (32-bit) arithmetic (Colamd::recommended), and when
    // COLAMD then fails, Ordering.h writes the permutation out of bounds before anything can
    // check it. Refuse sizes whose workspace does not fit, using the same formula in 64 bits.
    {
        namespace colamd = Eigen::internal::Colamd;
        const auto nnz = static_cast<std::int64_t>(a.nonzeros());
        const auto cols = static_cast<std::int64_t>(n);
        const auto words = [](std::int64_t count, std::size_t structure_bytes) {
            return (count + 1) * static_cast<std::int64_t>(structure_bytes) /
                   static_cast<std::int64_t>(sizeof(Index));
        };
        const std::int64_t workspace = 2 * nnz + words(cols, sizeof(colamd::ColStructure<Index>)) +
                                       words(cols, sizeof(colamd::RowStructure<Index>)) + cols +
                                       nnz / 5;
        if (workspace > std::numeric_limits<Index>::max()) {
            lu_.reset();
            return std::unexpected(base::Error{
                base::ErrorCode::resource_exhausted,
                "matrix too large for Eigen's 32-bit COLAMD workspace",
                base::ErrorContext{.index = std::nullopt, .value = static_cast<double>(workspace)}});
        }
    }

    lu_.emplace();
    lu_->analyzePattern(csc_);

    // Eigen checks COLAMD's result only with eigen_assert, which is compiled out in release
    // builds. With the size guard above COLAMD cannot run out of workspace; verify the column
    // ordering is a permutation anyway, as a last line of defence.
    const auto& perm = lu_->colsPermutation().indices();
    factored_to_original_.assign(n, -1);
    bool valid = static_cast<std::size_t>(perm.size()) == n;
    for (std::size_t i = 0; valid && i < n; ++i) {
        const Index j = perm(static_cast<Eigen::Index>(i));
        valid = j >= 0 && static_cast<std::size_t>(j) < n &&
                factored_to_original_[static_cast<std::size_t>(j)] < 0;
        if (valid) {
            factored_to_original_[static_cast<std::size_t>(j)] = static_cast<Index>(i);
        }
    }
    if (!valid) {
        lu_.reset();
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                           "COLAMD ordering failed", std::nullopt});
    }
    return {};
}

template <MatrixScalar Scalar>
std::expected<PivotRatio, base::Error> SparseLuBackend<Scalar>::factorize(
    std::span<const Scalar> csr_values) {
    NITCAD_EXPECTS(lu_.has_value());
    NITCAD_EXPECTS(csr_values.size() == csr_to_csc_.size());
    Scalar* const values = csc_.valuePtr();
    for (std::size_t k = 0; k < csr_values.size(); ++k) {
        values[static_cast<std::size_t>(csr_to_csc_[k])] = csr_values[k];
    }

    lu_->factorize(csc_);
    // Eigen 5.0.1 leaves info() unset when it cannot allocate its working memory, but always sets
    // the message on failure. lu_ is recreated by every analysis and the caller re-analyzes after
    // a failure, so an empty message means this factorization did not fail; info() is read only
    // then.
    const std::string message = lu_->lastErrorMessage();
    if (!message.empty() || lu_->info() != Eigen::Success) {
        if (std::string_view(message).starts_with("UNABLE TO")) {
            return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                               "Eigen SparseLU: " + message, std::nullopt});
        }
        // An exactly zero pivot. Eigen calls it "STRUCTURALLY SINGULAR ... ZERO COLUMN AT k",
        // also when the cause is numerical; k is 1-based in factored column order.
        std::optional<std::size_t> column;
        constexpr std::string_view marker = "ZERO COLUMN AT ";
        if (const auto at = message.find(marker); at != std::string::npos) {
            std::size_t k = 0;
            const char* first = message.data() + at + marker.size();
            const auto [ptr, ec] = std::from_chars(first, message.data() + message.size(), k);
            if (ec == std::errc{} && k >= 1 && k <= factored_to_original_.size()) {
                column = static_cast<std::size_t>(factored_to_original_[k - 1]);
            }
        }
        return std::unexpected(base::Error{base::ErrorCode::singular_system,
                                           "zero pivot in LU factorization",
                                           base::ErrorContext{.index = column, .value = 0.0}});
    }

    // The pivots are the diagonal of U, stored in the diagonal blocks of the supernodal L (read the
    // same way as Eigen's own logAbsDeterminant()).
    const auto& l = lu_->matrixL().m_mapL;
    using Supernodal = std::remove_cvref_t<decltype(l)>;
    double smallest = std::numeric_limits<double>::infinity();
    double largest = 0.0;
    std::size_t smallest_at = 0;
    for (Eigen::Index j = 0; j < l.cols(); ++j) {
        for (typename Supernodal::InnerIterator it(l, j); it; ++it) {
            if (it.index() == j) {
                const double pivot = std::abs(it.value());
                if (!std::isfinite(pivot)) {
                    // Overflow during elimination; a NaN would otherwise be ignored by the
                    // comparisons below and the factorization reported as usable.
                    return std::unexpected(base::Error{
                        base::ErrorCode::inaccurate_solve,
                        "non-finite pivot (overflow in LU factorization)",
                        base::ErrorContext{
                            .index = static_cast<std::size_t>(
                                factored_to_original_[static_cast<std::size_t>(j)]),
                            .value = reported(it.value())}});
                }
                if (pivot < smallest) {
                    smallest = pivot;
                    smallest_at = static_cast<std::size_t>(j);
                }
                largest = std::max(largest, pivot);
            }
        }
    }
    return PivotRatio{.ratio = largest > 0.0 ? smallest / largest : 0.0,
                      .column = static_cast<std::size_t>(factored_to_original_[smallest_at])};
}

template <MatrixScalar Scalar>
void SparseLuBackend<Scalar>::solve(std::span<const Scalar> b, std::span<Scalar> x) {
    NITCAD_EXPECTS(lu_.has_value());
    const auto n = static_cast<Eigen::Index>(csc_.rows());
    NITCAD_EXPECTS(b.size() == static_cast<std::size_t>(n) && x.size() == b.size());
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    const Eigen::Map<const Vector> bv(b.data(), n);
    Eigen::Map<Vector> xv(x.data(), n);
    xv = lu_->solve(bv);
}

}  // namespace NiTCAD::linalg::detail

namespace NiTCAD::linalg::detail {

template class SparseLuBackend<double>;
template class SparseLuBackend<std::complex<double>>;

}  // namespace NiTCAD::linalg::detail
