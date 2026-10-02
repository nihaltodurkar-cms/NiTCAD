#include "NiTCAD/linalg/linear_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <new>
#include <string>
#include <utility>

#include "NiTCAD/base/contract.hpp"
#include "eigen_sparse_lu.hpp"

namespace NiTCAD::linalg {

namespace {

constexpr double epsilon = std::numeric_limits<double>::epsilon();

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

// The power of two s with largest * s in [0.5, 1), kept finite for extreme magnitudes.
double power_of_two_scale(double largest) noexcept {
    int exponent = 0;
    std::frexp(largest, &exponent);
    return std::ldexp(1.0, std::clamp(-exponent, -1022, 1022));
}

double inf_norm(std::span<const double> v) noexcept {
    double m = 0.0;
    for (const double x : v) {
        m = std::max(m, std::abs(x));
    }
    return m;
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
    if (!std::isfinite(config.max_backward_error) || config.max_backward_error <= 0.0) {
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "max_backward_error must be positive and finite",
                                     std::nullopt, config.max_backward_error));
    }
    if (config.max_refinement_steps < 0) {
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "max_refinement_steps must not be negative", std::nullopt,
                                     static_cast<double>(config.max_refinement_steps)));
    }
    if (!(config.min_pivot_ratio >= 0.0 && config.min_pivot_ratio < 1.0)) {
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "min_pivot_ratio must be in [0, 1)", std::nullopt,
                                     config.min_pivot_ratio));
    }
    if (config.min_pivot_ratio > 0.0 && !config.equilibrate) {
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "min_pivot_ratio is calibrated on equilibrated matrices; "
                                     "set it to 0 when equilibrate is false",
                                     std::nullopt, config.min_pivot_ratio));
    }
    return LinearSolver(config);
}

LinearSolver::LinearSolver(const SolverConfig& config)
    : config_(config), backend_(std::make_unique<detail::SparseLuBackend>()) {}

LinearSolver::LinearSolver(LinearSolver&& other) noexcept
    : config_(other.config_),
      backend_(std::move(other.backend_)),
      matrix_(std::move(other.matrix_)),
      scaled_values_(std::move(other.scaled_values_)),
      row_scale_(std::move(other.row_scale_)),
      col_scale_(std::move(other.col_scale_)),
      weighted_row_max_(std::move(other.weighted_row_max_)),
      residual_(std::move(other.residual_)),
      work_rhs_(std::move(other.work_rhs_)),
      work_x_(std::move(other.work_x_)),
      correction_(std::move(other.correction_)),
      previous_x_(std::move(other.previous_x_)),
      analyzed_(std::exchange(other.analyzed_, false)),
      factorized_(std::exchange(other.factorized_, false)),
      analyses_(other.analyses_),
      factorizations_(other.factorizations_) {}

LinearSolver& LinearSolver::operator=(LinearSolver&& other) noexcept {
    if (this != &other) {
        config_ = other.config_;
        backend_ = std::move(other.backend_);
        matrix_ = std::move(other.matrix_);
        scaled_values_ = std::move(other.scaled_values_);
        row_scale_ = std::move(other.row_scale_);
        col_scale_ = std::move(other.col_scale_);
        weighted_row_max_ = std::move(other.weighted_row_max_);
        residual_ = std::move(other.residual_);
        work_rhs_ = std::move(other.work_rhs_);
        work_x_ = std::move(other.work_x_);
        correction_ = std::move(other.correction_);
        previous_x_ = std::move(other.previous_x_);
        analyzed_ = std::exchange(other.analyzed_, false);
        factorized_ = std::exchange(other.factorized_, false);
        analyses_ = other.analyses_;
        factorizations_ = other.factorizations_;
    }
    return *this;
}

LinearSolver::~LinearSolver() = default;

void LinearSolver::reset() noexcept {
    analyzed_ = false;
    factorized_ = false;
}

std::expected<void, base::Error> LinearSolver::analyze(const SparseMatrix& a) {
    NITCAD_EXPECTS(backend_ != nullptr);
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
    // A row or column without entries makes the matrix structurally singular.
    const auto n = static_cast<std::size_t>(a.rows());
    const std::span<const Index> offsets = a.row_offsets();
    for (std::size_t r = 0; r < n; ++r) {
        if (offsets[r] == offsets[r + 1]) {
            return std::unexpected(
                error(base::ErrorCode::singular_system, "matrix row has no entries", r));
        }
    }
    {
        std::vector<bool> used(n, false);
        for (const Index c : a.col_indices()) {
            used[static_cast<std::size_t>(c)] = true;
        }
        if (const auto it = std::ranges::find(used, false); it != used.end()) {
            return std::unexpected(error(base::ErrorCode::singular_system,
                                         "matrix column has no entries",
                                         static_cast<std::size_t>(it - used.begin())));
        }
    }
    if (auto analyzed = backend_->analyze(a); !analyzed) {
        return analyzed;
    }
    matrix_ = a;
    scaled_values_.resize(a.nonzeros());
    for (auto* v : {&row_scale_, &col_scale_, &weighted_row_max_, &residual_, &work_rhs_, &work_x_,
                    &correction_, &previous_x_}) {
        v->resize(n);
    }
    analyzed_ = true;
    ++analyses_;
    return {};
}

std::expected<FactorizationReport, base::Error> LinearSolver::factorize(const SparseMatrix& a) {
    NITCAD_EXPECTS(backend_ != nullptr);
    NITCAD_EXPECTS(a.rows() == a.cols());
    try {
        return factorize_unchecked(a);
    } catch (const std::bad_alloc&) {
        reset();
        return std::unexpected(out_of_memory());
    }
}

std::expected<FactorizationReport, base::Error> LinearSolver::factorize_unchecked(
    const SparseMatrix& a) {
    factorized_ = false;
    if (!analyzed_ || !a.has_same_pattern(matrix_)) {
        if (auto analyzed = analyze_unchecked(a); !analyzed) {
            return std::unexpected(std::move(analyzed.error()));
        }
    }
    const std::span<const double> values = a.values();
    const std::span<const Index> offsets = a.row_offsets();
    const std::span<const Index> columns = a.col_indices();
    if (const std::size_t k = first_non_finite(values); k < values.size()) {
        const auto row = std::ranges::upper_bound(offsets, static_cast<Index>(k)) -
                         offsets.begin() - 1;
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "matrix contains a non-finite value",
                                     static_cast<std::size_t>(row), values[k]));
    }

    // Row scales from the row maxima, then column scales from the column maxima of the
    // row-scaled matrix. Always computed: the backward error uses the column scales.
    const std::size_t n = row_scale_.size();
    for (std::size_t r = 0; r < n; ++r) {
        double largest = 0.0;
        for (auto k = static_cast<std::size_t>(offsets[r]);
             k < static_cast<std::size_t>(offsets[r + 1]); ++k) {
            largest = std::max(largest, std::abs(values[k]));
        }
        if (largest == 0.0) {
            reset();
            return std::unexpected(
                error(base::ErrorCode::singular_system, "matrix row is all zeros", r));
        }
        row_scale_[r] = power_of_two_scale(largest);
    }
    std::ranges::fill(col_scale_, 0.0);  // column maxima first
    for (std::size_t r = 0; r < n; ++r) {
        for (auto k = static_cast<std::size_t>(offsets[r]);
             k < static_cast<std::size_t>(offsets[r + 1]); ++k) {
            double& m = col_scale_[static_cast<std::size_t>(columns[k])];
            m = std::max(m, std::abs(values[k]) * row_scale_[r]);
        }
    }
    for (std::size_t c = 0; c < n; ++c) {
        if (col_scale_[c] == 0.0) {
            reset();
            return std::unexpected(
                error(base::ErrorCode::singular_system, "matrix column is all zeros", c));
        }
        col_scale_[c] = power_of_two_scale(col_scale_[c]);
    }
    for (std::size_t r = 0; r < n; ++r) {
        double weighted = 0.0;
        for (auto k = static_cast<std::size_t>(offsets[r]);
             k < static_cast<std::size_t>(offsets[r + 1]); ++k) {
            const double c = col_scale_[static_cast<std::size_t>(columns[k])];
            weighted = std::max(weighted, std::abs(values[k]) * c);
            scaled_values_[k] = config_.equilibrate ? values[k] * row_scale_[r] * c : values[k];
        }
        weighted_row_max_[r] = weighted;
    }

    const auto pivots = backend_->factorize(scaled_values_);
    if (!pivots) {
        reset();
        return std::unexpected(pivots.error());
    }
    if (pivots->ratio < config_.min_pivot_ratio) {
        reset();
        return std::unexpected(error(base::ErrorCode::singular_system,
                                     "pivot ratio below min_pivot_ratio: the matrix is singular "
                                     "to working precision (for example a floating region)",
                                     pivots->column, pivots->ratio));
    }
    std::ranges::copy(values, matrix_.values().begin());
    factorized_ = true;
    ++factorizations_;
    return FactorizationReport{.pivot_ratio = pivots->ratio};
}

void LinearSolver::scaled_solve(std::span<const double> rhs, std::span<double> x) {
    if (!config_.equilibrate) {
        backend_->solve(rhs, x);
        return;
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        work_rhs_[i] = rhs[i] * row_scale_[i];
    }
    backend_->solve(work_rhs_, work_x_);
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] = work_x_[i] * col_scale_[i];
    }
}

double LinearSolver::backward_error(std::span<const double> b, std::span<const double> x) {
    const std::span<const Index> offsets = matrix_.row_offsets();
    const std::span<const Index> columns = matrix_.col_indices();
    const std::span<const double> values = matrix_.values();
    // ||x|| in column-equilibrated units, so the fallback below is invariant to column scaling
    // (with plain ||x||_inf, columns scaled by 1e+-10 make it fire on every row and hide errors).
    double z_norm = 0.0;
    for (std::size_t j = 0; j < x.size(); ++j) {
        z_norm = std::max(z_norm, std::abs(x[j]) / col_scale_[j]);
    }
    const double rounding = 1000.0 * static_cast<double>(b.size()) * epsilon;
    double worst = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        double ax = 0.0;
        double abs_ax = 0.0;
        for (auto k = static_cast<std::size_t>(offsets[i]);
             k < static_cast<std::size_t>(offsets[i + 1]); ++k) {
            const double term = values[k] * x[static_cast<std::size_t>(columns[k])];
            ax += term;
            abs_ax += std::abs(term);
        }
        const double r = b[i] - ax;
        residual_[i] = r;
        // Oettli-Prager denominator, replaced by the Arioli-Demmel-Duff one where it is at
        // rounding level (for example a row with b_i = 0 whose x entries nearly cancel).
        const double bound = weighted_row_max_[i] * z_norm;
        double denominator = abs_ax + std::abs(b[i]);
        if (denominator <= rounding * (bound + std::abs(b[i]))) {
            denominator = abs_ax + bound;
        }
        const double row_error = denominator > 0.0 ? std::abs(r) / denominator
                                 : r == 0.0       ? 0.0
                                                  : std::numeric_limits<double>::infinity();
        // A NaN (from overflow) must not be lost by max.
        worst = std::isnan(row_error) ? row_error : std::max(worst, row_error);
        if (std::isnan(worst)) {
            break;
        }
    }
    return worst;
}

std::expected<SolveReport, base::Error> LinearSolver::solve(std::span<const double> b,
                                                            std::span<double> x) {
    NITCAD_EXPECTS(backend_ != nullptr);
    NITCAD_EXPECTS(factorized_);
    NITCAD_EXPECTS(b.size() == static_cast<std::size_t>(matrix_.rows()) && x.size() == b.size());
    NITCAD_EXPECTS(!overlaps(b, x));

    if (const std::size_t k = first_non_finite(b); k < b.size()) {
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "right-hand side contains a non-finite value", k, b[k]));
    }
    double omega = 0.0;
    std::size_t steps = 0;
    try {
        scaled_solve(b, x);
        if (const std::size_t k = first_non_finite(x); k < x.size()) {
            return std::unexpected(error(base::ErrorCode::inaccurate_solve,
                                         "solution contains a non-finite value", k));
        }
        omega = backward_error(b, x);
        for (int attempt = 0;
             !(omega <= config_.max_backward_error) && attempt < config_.max_refinement_steps;
             ++attempt) {
            scaled_solve(residual_, correction_);
            std::ranges::copy(x, previous_x_.begin());
            for (std::size_t i = 0; i < x.size(); ++i) {
                x[i] += correction_[i];
            }
            const double refined = first_non_finite(x) < x.size()
                                       ? std::numeric_limits<double>::quiet_NaN()
                                       : backward_error(b, x);
            if (!(refined < omega)) {
                // No progress: keep the better iterate, restore its residual, and stop.
                std::ranges::copy(previous_x_, x.begin());
                omega = backward_error(b, x);
                break;
            }
            omega = refined;
            ++steps;
        }
    } catch (const std::bad_alloc&) {
        return std::unexpected(out_of_memory());
    }
    if (!(omega <= config_.max_backward_error)) {
        return std::unexpected(error(base::ErrorCode::inaccurate_solve,
                                     "componentwise backward error above max_backward_error",
                                     std::nullopt, omega));
    }
    const double b_norm = inf_norm(b);
    const double r_norm = inf_norm(residual_);
    return SolveReport{.backward_error = omega,
                       .relative_residual = b_norm > 0.0 ? r_norm / b_norm : r_norm,
                       .refinement_steps = steps,
                       .iterations = std::nullopt};
}

}  // namespace NiTCAD::linalg
