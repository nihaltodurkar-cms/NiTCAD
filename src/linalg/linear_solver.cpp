#include "NiTCAD/linalg/linear_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <new>
#include <string>
#include <utility>

#include "NiTCAD/base/contract.hpp"
#include "eigen_sparse_lu.hpp"

namespace NiTCAD::linalg {

namespace {

base::Error error(base::ErrorCode code, std::string message,
                  std::optional<std::size_t> index = std::nullopt,
                  std::optional<double> value = std::nullopt) {
    return {code, std::move(message), base::ErrorContext{.index = index, .value = value}};
}

// The message fits the small-string buffer, so reporting it does not allocate.
base::Error out_of_memory() {
    return {base::ErrorCode::resource_exhausted, "out of memory", std::nullopt};
}

// Index of the first non-finite value, or values.size() if all are finite.
std::size_t first_non_finite(std::span<const double> values) noexcept {
    const auto it = std::ranges::find_if(values, [](double v) { return !std::isfinite(v); });
    return static_cast<std::size_t>(it - values.begin());
}

bool overlaps(std::span<const double> a, std::span<const double> b) noexcept {
    if (a.empty() || b.empty()) {
        return false;
    }
    const std::less<const double*> before;
    return before(a.data(), b.data() + b.size()) && before(b.data(), a.data() + a.size());
}

}  // namespace

std::expected<LinearSolver, base::Error> LinearSolver::create(const SolverConfig& config) {
    switch (config.backend) {
        case SolverBackend::eigen_sparse_lu:
            if (config.threads != 1) {
                return std::unexpected(error(base::ErrorCode::invalid_input,
                                             "eigen_sparse_lu supports exactly 1 thread",
                                             std::nullopt, static_cast<double>(config.threads)));
            }
            break;
        default:
            return std::unexpected(
                error(base::ErrorCode::invalid_input, "unknown linear solver backend"));
    }
    if (!std::isfinite(config.max_relative_residual) || config.max_relative_residual <= 0.0) {
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "max_relative_residual must be positive and finite",
                                     std::nullopt, config.max_relative_residual));
    }
    return LinearSolver(config);
}

LinearSolver::LinearSolver(const SolverConfig& config)
    : config_(config), backend_(std::make_unique<detail::SparseLuBackend>()) {}

LinearSolver::LinearSolver(LinearSolver&&) noexcept = default;
LinearSolver& LinearSolver::operator=(LinearSolver&&) noexcept = default;
LinearSolver::~LinearSolver() = default;

void LinearSolver::reset() noexcept {
    analyzed_ = false;
    factorized_ = false;
}

std::expected<void, base::Error> LinearSolver::analyze(const SparseMatrix& a) {
    NITCAD_EXPECTS(a.rows() == a.cols());
    try {
        return analyze_unchecked(a);
    } catch (const std::bad_alloc&) {
        reset();
        return std::unexpected(out_of_memory());
    }
}

std::expected<void, base::Error> LinearSolver::analyze_unchecked(const SparseMatrix& a) {
    reset();
    if (a.rows() == 0) {
        return std::unexpected(error(base::ErrorCode::invalid_input, "empty linear system"));
    }
    if (auto analyzed = backend_->analyze(a); !analyzed) {
        return analyzed;
    }
    matrix_ = a;
    residual_.resize(static_cast<std::size_t>(a.rows()));
    analyzed_ = true;
    ++analyses_;
    return {};
}

std::expected<void, base::Error> LinearSolver::factorize(const SparseMatrix& a) {
    NITCAD_EXPECTS(a.rows() == a.cols());
    try {
        factorized_ = false;
        if (!analyzed_ || !a.has_same_pattern(matrix_)) {
            if (auto analyzed = analyze_unchecked(a); !analyzed) {
                return analyzed;
            }
        }
        const std::span<const double> values = a.values();
        if (const std::size_t k = first_non_finite(values); k < values.size()) {
            const std::span<const Index> offsets = a.row_offsets();
            const auto row = std::ranges::upper_bound(offsets, static_cast<Index>(k)) -
                             offsets.begin() - 1;
            return std::unexpected(error(base::ErrorCode::invalid_input,
                                         "matrix contains a non-finite value",
                                         static_cast<std::size_t>(row), values[k]));
        }
        if (auto factorized = backend_->factorize(values); !factorized) {
            reset();
            return factorized;
        }
        std::ranges::copy(values, matrix_.values().begin());
        factorized_ = true;
        ++factorizations_;
        return {};
    } catch (const std::bad_alloc&) {
        reset();
        return std::unexpected(out_of_memory());
    }
}

std::expected<SolveReport, base::Error> LinearSolver::solve(std::span<const double> b,
                                                            std::span<double> x) {
    NITCAD_EXPECTS(factorized_);
    NITCAD_EXPECTS(b.size() == static_cast<std::size_t>(matrix_.rows()) && x.size() == b.size());
    NITCAD_EXPECTS(!overlaps(b, x));

    if (const std::size_t k = first_non_finite(b); k < b.size()) {
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "right-hand side contains a non-finite value", k, b[k]));
    }
    try {
        backend_->solve(b, x);
    } catch (const std::bad_alloc&) {
        return std::unexpected(out_of_memory());
    }
    if (const std::size_t k = first_non_finite(x); k < x.size()) {
        return std::unexpected(error(base::ErrorCode::inaccurate_solve,
                                     "solution contains a non-finite value", k));
    }

    // Acceptance check, independent of what the backend reported.
    matrix_.multiply(x, residual_);
    double residual_sq = 0.0;
    double b_sq = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        const double r = residual_[i] - b[i];
        residual_sq += r * r;
        b_sq += b[i] * b[i];
    }
    const double relative = b_sq > 0.0 ? std::sqrt(residual_sq / b_sq) : std::sqrt(residual_sq);
    if (!std::isfinite(relative) || relative > config_.max_relative_residual) {
        return std::unexpected(error(base::ErrorCode::inaccurate_solve,
                                     "relative residual above max_relative_residual",
                                     std::nullopt, relative));
    }
    return SolveReport{.relative_residual = relative, .iterations = std::nullopt};
}

}  // namespace NiTCAD::linalg
