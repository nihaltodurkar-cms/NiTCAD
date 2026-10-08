#include "mkl_pardiso.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::linalg::detail {

namespace {

constexpr MklInt one = 1;
constexpr MklInt quiet = 0;
constexpr MklInt phase_analyze = 11, phase_factorize = 22, phase_solve = 33, phase_release = -1;

// PARDISO's error codes (oneMKL developer reference) as NiTCAD errors.
base::Error pardiso_error(MklInt code, const char* step) {
    base::ErrorCode kind = base::ErrorCode::inaccurate_solve;
    std::string what = "internal error";
    switch (code) {
        case -1: kind = base::ErrorCode::invalid_input; what = "inconsistent input"; break;
        case -2: kind = base::ErrorCode::resource_exhausted; what = "not enough memory"; break;
        case -4: kind = base::ErrorCode::singular_system; what = "zero pivot"; break;
        case -7: kind = base::ErrorCode::singular_system; what = "diagonal matrix is singular"; break;
        case -8: kind = base::ErrorCode::resource_exhausted; what = "32-bit integer overflow"; break;
        case -9: kind = base::ErrorCode::resource_exhausted; what = "not enough memory"; break;
        case -12: kind = base::ErrorCode::resource_exhausted; what = "64-bit interface needed"; break;
        default: break;
    }
    return {kind, std::string("PARDISO ") + step + ": " + what,
            base::ErrorContext{.index = std::nullopt, .value = static_cast<double>(code)}};
}

}  // namespace

template <MatrixScalar Scalar>
PardisoBackend<Scalar>::PardisoBackend(int threads) : api_(mkl_api()), threads_(threads) {
    NITCAD_EXPECTS(api_ != nullptr);
    NITCAD_EXPECTS(threads >= 1);
}

template <MatrixScalar Scalar>
PardisoBackend<Scalar>::~PardisoBackend() {
    release();
}

template <MatrixScalar Scalar>
void PardisoBackend<Scalar>::release() noexcept {
    if (!held_) return;
    MklInt perm = 0, error = 0;
    api_->pardiso(pt_.data(), &one, &one, &matrix_type, &phase_release, &n_, values_.data(),
                  row_offsets_.data(), columns_.data(), &perm, &one, iparm_.data(), &quiet,
                  nullptr, nullptr, &error);
    held_ = false;
    analyzed_ = false;
    pt_.fill(nullptr);
}

template <MatrixScalar Scalar>
MklInt PardisoBackend<Scalar>::call(MklInt phase, Scalar* b, Scalar* x) {
    const int previous = api_->set_num_threads_local(threads_);
    MklInt perm = 0, error = 0;
    api_->pardiso(pt_.data(), &one, &one, &matrix_type, &phase, &n_, values_.data(),
                  row_offsets_.data(), columns_.data(), &perm, &one, iparm_.data(), &quiet, b, x,
                  &error);
    api_->set_num_threads_local(previous);
    // Phase ab runs PARDISO's steps a to b: 1 analysis, 2 numeric factorization, 3 solve.
    const MklInt first = phase / 10, last = phase % 10;
    if (first <= 1 && 1 <= last) ++this->counts_.analyses;
    if (first <= 2 && 2 <= last) ++this->counts_.factorizations;
    if (first <= 3 && 3 <= last) ++this->counts_.solves;
    return error;
}

template <MatrixScalar Scalar>
std::expected<void, base::Error> PardisoBackend<Scalar>::analyze(
    const BasicSparseMatrix<Scalar>& a) {
    NITCAD_EXPECTS(a.rows() == a.cols() && a.rows() > 0);
    release();
    n_ = static_cast<MklInt>(a.rows());
    row_offsets_.assign(a.row_offsets().begin(), a.row_offsets().end());
    columns_.assign(a.col_indices().begin(), a.col_indices().end());
    values_.assign(a.nonzeros(), Scalar{});
    rhs_.resize(static_cast<std::size_t>(n_));

    iparm_.fill(0);
    api_->pardisoinit(pt_.data(), &matrix_type, iparm_.data());
    iparm_[0] = 1;   // the settings below, not the defaults
    iparm_[1] = 2;   // METIS nested dissection
    iparm_[5] = 0;   // the solution into x, b unchanged
    iparm_[7] = 0;   // no iterative refinement (BasicLinearSolver refines on the original system)
    iparm_[9] = 13;  // perturb pivots below 1e-13 of the largest
    iparm_[10] = 1;  // scaling and
    iparm_[12] = 1;  // weighted matching (the unsymmetric defaults), computed in phase 11
    iparm_[26] = 0;  // no matrix checker: BasicLinearSolver checked the input
    iparm_[27] = 0;  // double precision
    iparm_[34] = 1;  // zero-based indices
    analyzed_ = false;
    return {};
}

template <MatrixScalar Scalar>
std::expected<Factorization, base::Error> PardisoBackend<Scalar>::factorize(
    std::span<const Scalar> csr_values) {
    NITCAD_EXPECTS(csr_values.size() == values_.size());
    std::ranges::copy(csr_values, values_.begin());
    if (!analyzed_) {
        held_ = true;
        if (const MklInt error = call(phase_analyze, nullptr, nullptr); error != 0) {
            release();
            return std::unexpected(pardiso_error(error, "analysis"));
        }
        analyzed_ = true;
    }
    if (const MklInt error = call(phase_factorize, nullptr, nullptr); error != 0) {
        // BasicLinearSolver re-analyzes after a failed factorization (release, then phase 11).
        return std::unexpected(pardiso_error(error, "factorization"));
    }
    // No pivot ratio: the factor's diagonal (pardiso_getdiag) is not a singularity measure here.
    // PARDISO pivots statically (weighted matching, then pivoting within supernodes), and a valid
    // equilibrated system (tests/linalg badly_scaled_coupled) has a diagonal ratio of 1.8e-13 with
    // no perturbed pivot, under the 1e-11 that Eigen's partial-pivoting ratio is calibrated for.
    // The perturbed-pivot count is the check (BasicLinearSolver).
    return Factorization{
        .pivots = std::nullopt,
        .perturbed_pivots = static_cast<std::size_t>(std::max<MklInt>(iparm_[13], 0))};
}

template <MatrixScalar Scalar>
void PardisoBackend<Scalar>::solve(std::span<const Scalar> b, std::span<Scalar> x) {
    NITCAD_EXPECTS(analyzed_);
    NITCAD_EXPECTS(b.size() == rhs_.size() && x.size() == b.size());
    std::ranges::copy(b, rhs_.begin());
    if (call(phase_solve, rhs_.data(), x.data()) != 0) {
        std::ranges::fill(x, Scalar{std::numeric_limits<double>::quiet_NaN()});
    }
}

template class PardisoBackend<double>;
template class PardisoBackend<std::complex<double>>;

}  // namespace NiTCAD::linalg::detail
