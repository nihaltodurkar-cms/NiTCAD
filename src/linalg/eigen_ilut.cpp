#include "eigen_ilut.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>

#include "NiTCAD/base/contract.hpp"

// The zero-pivot rule, the drop rules and the use of info() below were read from Eigen 5.0.1's
// IncompleteLUT.h. A different Eigen must be re-checked against them before this assertion is
// changed.
static_assert(EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 5 && EIGEN_MINOR_VERSION == 0 &&
                  EIGEN_PATCH_VERSION == 1,
              "eigen_ilut.cpp was verified against Eigen 5.0.1 only");

namespace NiTCAD::linalg::detail {

namespace {

bool is_finite(double v) noexcept { return std::isfinite(v); }
bool is_finite(std::complex<double> v) noexcept {
    return std::isfinite(v.real()) && std::isfinite(v.imag());
}

}  // namespace

template <MatrixScalar Scalar>
IlutBackend<Scalar>::IlutBackend(double drop_tolerance, int fill_factor)
    : drop_tolerance_(drop_tolerance), fill_factor_(fill_factor) {}

template <MatrixScalar Scalar>
std::expected<void, base::Error> IlutBackend<Scalar>::analyze(const BasicSparseMatrix<Scalar>& a) {
    NITCAD_EXPECTS(a.rows() == a.cols() && a.rows() > 0);
    built_ = false;
    const auto n = static_cast<std::size_t>(a.rows());
    // The CSR pattern as it is: both are row-major with sorted columns.
    csr_.resize(a.rows(), a.cols());
    csr_.resizeNonZeros(static_cast<Eigen::Index>(a.nonzeros()));
    std::ranges::copy(a.row_offsets(), csr_.outerIndexPtr());
    std::ranges::copy(a.col_indices(), csr_.innerIndexPtr());
    std::fill(csr_.valuePtr(), csr_.valuePtr() + a.nonzeros(), Scalar{});
    NITCAD_EXPECTS(static_cast<std::size_t>(csr_.outerSize()) == n);

    ilut_.emplace();
    ilut_->setDroptol(drop_tolerance_);
    ilut_->setFillfactor(fill_factor_);
    ilut_->analyzePattern(csr_);
    ++this->counts_.analyses;
    return {};
}

template <MatrixScalar Scalar>
std::expected<Factorization, base::Error> IlutBackend<Scalar>::factorize(
    std::span<const Scalar> csr_values) {
    NITCAD_EXPECTS(ilut_.has_value());
    NITCAD_EXPECTS(csr_values.size() == static_cast<std::size_t>(csr_.nonZeros()));
    built_ = false;
    std::ranges::copy(csr_values, csr_.valuePtr());
    ilut_->factorize(csr_);
    ++this->counts_.factorizations;
    // NumericalIssue only for a row without values, which BasicLinearSolver rejects before.
    if (ilut_->info() != Eigen::Success) {
        return std::unexpected(base::Error{base::ErrorCode::inaccurate_solve,
                                           "ILUT preconditioner could not be built", std::nullopt});
    }
    // Without pivoting an elimination can overflow, or a zero pivot left by a zero drop tolerance
    // divides by zero: the preconditioner is then unusable.
    const Matrix& factors = ilut_->factors();
    for (Eigen::Index r = 0; r < factors.outerSize(); ++r) {
        for (typename Matrix::InnerIterator it(factors, r); it; ++it) {
            if (!is_finite(it.value())) {
                return std::unexpected(base::Error{
                    base::ErrorCode::inaccurate_solve,
                    "ILUT preconditioner has a non-finite entry",
                    base::ErrorContext{.index = static_cast<std::size_t>(r), .value = std::nullopt}});
            }
        }
    }
    built_ = true;
    return Factorization{.pivots = std::nullopt, .perturbed_pivots = std::nullopt};
}

template <MatrixScalar Scalar>
void IlutBackend<Scalar>::solve(std::span<const Scalar> b, std::span<Scalar> x) {
    NITCAD_EXPECTS(built_);
    const auto n = static_cast<Eigen::Index>(csr_.rows());
    NITCAD_EXPECTS(b.size() == static_cast<std::size_t>(n) && x.size() == b.size());
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    const Eigen::Map<const Vector> bv(b.data(), n);
    Eigen::Map<Vector> xv(x.data(), n);
    xv = ilut_->solve(bv);
    ++this->counts_.solves;
}

template class IlutBackend<double>;
template class IlutBackend<std::complex<double>>;

}  // namespace NiTCAD::linalg::detail
