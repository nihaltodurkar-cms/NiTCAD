#include "eigen_sparse_lu.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::linalg::detail {

std::expected<void, base::Error> SparseLuBackend::analyze(const SparseMatrix& a) {
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
    std::fill(csc_.valuePtr(), csc_.valuePtr() + a.nonzeros(), 0.0);

    lu_.emplace();
    lu_->analyzePattern(csc_);
    return {};
}

std::expected<void, base::Error> SparseLuBackend::factorize(std::span<const double> csr_values) {
    NITCAD_EXPECTS(lu_.has_value());
    NITCAD_EXPECTS(csr_values.size() == csr_to_csc_.size());
    double* const values = csc_.valuePtr();
    for (std::size_t k = 0; k < csr_values.size(); ++k) {
        values[static_cast<std::size_t>(csr_to_csc_[k])] = csr_values[k];
    }

    lu_->factorize(csc_);
    // Eigen 5.0.1 leaves info() unset when it cannot allocate its working memory, but always sets
    // the message on failure. lu_ is recreated by every analysis and the caller re-analyzes after
    // a failure, so an empty message means this factorization did not fail; info() is read only
    // then.
    const std::string message = lu_->lastErrorMessage();
    if (message.empty() && lu_->info() == Eigen::Success) {
        return {};
    }
    const bool out_of_memory = std::string_view(message).starts_with("UNABLE TO");
    return std::unexpected(base::Error{
        out_of_memory ? base::ErrorCode::resource_exhausted : base::ErrorCode::singular_system,
        "Eigen SparseLU factorization failed: " + message, std::nullopt});
}

void SparseLuBackend::solve(std::span<const double> b, std::span<double> x) {
    NITCAD_EXPECTS(lu_.has_value());
    const auto n = static_cast<Eigen::Index>(csc_.rows());
    NITCAD_EXPECTS(b.size() == static_cast<std::size_t>(n) && x.size() == b.size());
    const Eigen::Map<const Eigen::VectorXd> bv(b.data(), n);
    Eigen::Map<Eigen::VectorXd> xv(x.data(), n);
    xv = lu_->solve(bv);
}

}  // namespace NiTCAD::linalg::detail
