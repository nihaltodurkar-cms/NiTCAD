// Damped Newton iteration in scaled variables (ARCHITECTURE.md 6.2; legacy device.py solver
// options and Device1D::solve_equilibrium_boltzmann).
//
// Each iteration evaluates F and J at x, solves J dx = -F with the given LinearSolver, clips every
// component of dx to +-max_update (the legacy max_dpsi damping cap) and adds it to x.
// Convergence (decided at Unit 8, the question left open in 6.1): the iteration has converged when
// the largest component of the FULL correction dx, before clipping, is below tol_update. The legacy
// measured the clipped correction; with a componentwise clip and tol_update < max_update the two
// tests agree, but judging the full correction does not rely on that. The legacy 1D loop declares a
// residual tolerance (1e-10) and never uses it; there is none here, and each iteration's residual
// is recorded instead.
// With options.line_search (Unit 30) the damped step is shortened by backtracking until the squared
// residual falls enough (LineSearch below); the fraction taken is IterationRecord::step.
// The Jacobian pattern is fixed by the system, so the solver analyzes once per run.
// A system may replace the measure and the update with its own hooks (assemble::DriftDiffusion
// does, for relative density updates):
//     double update_size(std::span<const double> x, std::span<const double> dx) const;
//     void apply_update(std::span<double> x, std::span<const double> dx, double max_update) const;
#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <expected>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <stop_token>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/results/convergence.hpp"

namespace NiTCAD::solve {

// Backtracking line search on the damped step (Unit 30). The damped update x1 (the system's own
// clipping and density bounds applied to the full correction) is tried first; if it fails the
// test below, the step is halved along the segment from x to x1, x + lambda (x1 - x), at most
// max_backtracks times. Every trial is a convex combination of x and x1, so it keeps every bound
// x1 keeps (positive densities, temperatures in range). A trial whose residual is not finite
// fails. If no trial passes, the best one is taken (never lambda = 0), so the iteration always
// moves. Convergence is still judged on the full correction, so a short step is never taken for
// convergence, and a converging iteration (full correction below tol_update) takes the full
// damped step without a search. lambda = 1, when it passes, is exactly the update without a
// search, and its F and J are reused by the next iteration. Off by default: every existing run is
// unchanged.
// What a trial must achieve:
// - residual: the Armijo test on m(x) = sum_i F_i(x)^2, m(trial) <= (1 - 2 alpha lambda) m(x), or
//   the same test on the row-scaled m_s(x) = sum_i (F_i(x) / s_i)^2, or every row at its rounding,
//   |F_i(trial)| <= 64 eps s_i (Unit 50, DECISIONS.md RD1). s_i = sum_j |J_ij x_j| at the iterate
//   (1 where that is 0), the size of the terms row i sums. Near a solution m is the rounding of
//   the largest rows, and a row many decades smaller that still needs a real correction (a
//   minority density) is invisible in it; m_s sees it. Far from a solution the damped step need
//   not decrease m_s (measured: Armijo on m_s alone stalls an equilibrium solve), so m stays.
//   Trials are ranked on m when none passes, so an iteration differs from the plain test's only
//   where a trial passes that it rejected;
// - natural (the default): Deuflhard's natural monotonicity test, affine-invariant and so blind to
//   the very different scales of the rows: the simplified correction at the trial,
//   dx_bar = J(x)^-1 F(trial) with the factorization already made (one back-substitution, no
//   factorization), must be smaller than the full correction dx, |dx_bar| <= (1 - lambda / 4) |dx|,
//   both measured as the convergence test measures a correction (the system's update_size, or the
//   largest component).
enum class LineSearchTest : int { residual = 0, natural = 1 };

struct LineSearch {
    bool enabled = false;
    int max_backtracks = 10;            // the smallest trial is 2^-max_backtracks of the step
    double sufficient_decrease = 1e-4;  // alpha of the residual test, in (0, 0.5)
    LineSearchTest test = LineSearchTest::natural;
};

struct NewtonOptions {
    int max_iterations = 100;
    double tol_update = 1e-8;  // on the largest scaled correction
    double max_update = 5.0;   // componentwise damping cap on the scaled correction
    LineSearch line_search;
};

// Called after every iteration with its record and whether it converged (see control.hpp for the
// rules a callback must follow).
using IterationObserver = std::function<void(const results::IterationRecord&, bool converged)>;

// A system Newton can drive: a fixed number of unknowns, a Jacobian pattern, and F and J at x.
template <class S>
concept NewtonSystem = requires(const S& s, std::span<const double> x, std::span<double> f,
                                linalg::SparseMatrix& j) {
    { s.unknowns() } -> std::convertible_to<std::size_t>;
    { s.make_jacobian() } -> std::same_as<linalg::SparseMatrix>;
    s.evaluate(x, f, j);
};

// The convergence history goes into `record`, which is reset first and filled whatever the outcome.
// Errors:
// - invalid_input: max_iterations < 1, or tol_update or max_update not finite and positive, or
//   line_search with max_backtracks outside 0..60, sufficient_decrease outside (0, 0.5) or an
//   unknown test;
// - cancelled: the stop token was set before an iteration (context: the iteration it would have
//   been); the token is checked before every iteration, the first included;
// - non_convergence: no convergence within max_iterations (context: the iteration count and the
//   last full correction), or a non-finite residual (context: the iteration);
// - singular_system, inaccurate_solve, resource_exhausted: from the linear solver, with the
//   iteration added to the message.
// x is updated in place; on error it holds the last iterate. Precondition (NITCAD_EXPECTS via the
// system): x has system.unknowns() entries.
template <NewtonSystem System>
[[nodiscard]] std::expected<void, base::Error> newton_solve(
    const System& system, std::span<double> x, const NewtonOptions& options,
    linalg::LinearSolver& solver, results::ConvergenceRecord& record,
    const std::stop_token& stop = {}, const IterationObserver& observe = {}) {
    record = {};
    const auto invalid = [](const char* message, double value) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input, message,
            base::ErrorContext{.index = std::nullopt, .value = value}});
    };
    if (options.max_iterations < 1) {
        return invalid("Newton max_iterations must be at least 1", options.max_iterations);
    }
    if (!(std::isfinite(options.tol_update) && options.tol_update > 0.0)) {
        return invalid("Newton tol_update must be finite and positive", options.tol_update);
    }
    if (!(std::isfinite(options.max_update) && options.max_update > 0.0)) {
        return invalid("Newton max_update must be finite and positive", options.max_update);
    }
    const LineSearch& search = options.line_search;
    if (search.max_backtracks < 0 || search.max_backtracks > 60) {
        return invalid("Newton line_search.max_backtracks must lie in 0..60",
                       search.max_backtracks);
    }
    if (!(search.sufficient_decrease > 0.0 && search.sufficient_decrease < 0.5)) {
        return invalid("Newton line_search.sufficient_decrease must lie in (0, 0.5)",
                       search.sufficient_decrease);
    }
    if (search.test != LineSearchTest::residual && search.test != LineSearchTest::natural) {
        return invalid("Newton line_search.test is unknown", static_cast<double>(search.test));
    }

    const std::size_t n = system.unknowns();
    linalg::SparseMatrix jacobian = system.make_jacobian();
    std::vector<double> f(n), rhs(n), dx(n);
    // The damped update, as the system applies it.
    const auto damp = [&](std::span<double> to, std::span<const double> correction) {
        if constexpr (requires { system.apply_update(to, correction, options.max_update); }) {
            system.apply_update(to, correction, options.max_update);
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                to[i] += std::clamp(correction[i], -options.max_update, options.max_update);
            }
        }
    };
    // The size of a correction, as the convergence test measures it.
    const auto measure = [&](std::span<const double> at, std::span<const double> d) {
        if constexpr (requires { system.update_size(at, d); }) {
            return system.update_size(at, d);
        } else {
            (void)at;  // systems without update_size measure the largest component
            double size = 0.0;
            for (std::size_t i = 0; i < n; ++i) size = std::max(size, std::abs(d[i]));
            return size;
        }
    };
    // Line search work (LineSearch): the iterate, the damped update, a trial's residual and its
    // simplified correction.
    std::vector<double> x0, x1, f_trial, dx_bar;
    std::vector<double> scale;  // the residual test's row scales s_i (LineSearch)
    if (search.enabled) {
        f_trial.resize(n);
        dx_bar.resize(n);
        if (search.test == LineSearchTest::residual) scale.resize(n);
    }
    const auto merit = [n](std::span<const double> r) {
        double m = 0.0;
        for (std::size_t i = 0; i < n; ++i) m += r[i] * r[i];
        return std::isfinite(m) ? m : std::numeric_limits<double>::infinity();
    };
    const auto scaled_merit = [n](std::span<const double> r, std::span<const double> s) {
        double m = 0.0;
        for (std::size_t i = 0; i < n; ++i) m += (r[i] / s[i]) * (r[i] / s[i]);
        return std::isfinite(m) ? m : std::numeric_limits<double>::infinity();
    };
    bool evaluated = false;  // f and jacobian already hold F and J at x (an accepted trial)
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        if (stop.stop_requested()) {
            return std::unexpected(base::Error{
                base::ErrorCode::cancelled,
                "cancelled before Newton iteration " + std::to_string(iteration),
                base::ErrorContext{.index = static_cast<std::size_t>(iteration),
                                   .value = std::nullopt}});
        }
        if (!evaluated) system.evaluate(x, f, jacobian);
        evaluated = false;
        double residual = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(f[i])) {
                return std::unexpected(base::Error{
                    base::ErrorCode::non_convergence,
                    "Newton residual is not finite at iteration " + std::to_string(iteration),
                    base::ErrorContext{.index = static_cast<std::size_t>(iteration),
                                       .value = std::nullopt}});
            }
            residual = std::max(residual, std::abs(f[i]));
            rhs[i] = -f[i];
        }
        auto factored = solver.factorize(jacobian);
        if (!factored) {
            base::Error e = std::move(factored.error());
            e.message = "Newton iteration " + std::to_string(iteration) + ": " + e.message;
            return std::unexpected(std::move(e));
        }
        if (const auto ratio = factored->pivot_ratio) {
            record.smallest_pivot_ratio =
                std::min(record.smallest_pivot_ratio.value_or(*ratio), *ratio);
        }
        auto solved = solver.solve(rhs, dx);
        if (!solved) {
            base::Error e = std::move(solved.error());
            e.message = "Newton iteration " + std::to_string(iteration) + ": " + e.message;
            return std::unexpected(std::move(e));
        }
        const std::span<const double> step = dx;
        const double update = measure(x, step);
        const bool converged = update < options.tol_update;
        double taken = 1.0;
        if (!search.enabled || converged) {
            damp(x, step);
        } else {
            // Backtrack along the damped step (LineSearch).
            x0.assign(x.begin(), x.end());
            damp(x, step);
            x1.assign(x.begin(), x.end());
            const double m0 = merit(f);
            const bool natural = search.test == LineSearchTest::natural;
            double m0_scaled = 0.0;
            if (!natural) {
                // s_i = sum_j |J_ij x_j| at the iterate (jacobian still holds J(x)).
                const auto offsets = jacobian.row_offsets();
                const auto columns = jacobian.col_indices();
                const auto values = jacobian.values();
                for (std::size_t r = 0; r < n; ++r) {
                    double s = 0.0;
                    for (auto k = static_cast<std::size_t>(offsets[r]);
                         k < static_cast<std::size_t>(offsets[r + 1]); ++k) {
                        s += std::abs(values[k] * x0[static_cast<std::size_t>(columns[k])]);
                    }
                    scale[r] = s > 0.0 && std::isfinite(s) ? s : 1.0;
                }
                m0_scaled = scaled_merit(f, scale);
            }
            // A trial's figure (smaller is better) and whether it passes the test.
            const auto judge = [&](double lambda) -> std::pair<double, bool> {
                if (natural) {
                    for (std::size_t i = 0; i < n; ++i) {
                        if (!std::isfinite(f_trial[i])) {
                            return {std::numeric_limits<double>::infinity(), false};
                        }
                        rhs[i] = -f_trial[i];
                    }
                    if (!solver.solve(rhs, dx_bar)) {
                        return {std::numeric_limits<double>::infinity(), false};
                    }
                    const double size = measure(x, dx_bar);
                    return {size, size <= (1.0 - 0.25 * lambda) * update};
                }
                const double m = merit(f_trial);
                const double k = 1.0 - 2.0 * search.sufficient_decrease * lambda;
                if (m <= k * m0) return {m, true};
                // The row-scaled merit, or every row at its rounding (Unit 50, RD1).
                if (scaled_merit(f_trial, scale) <= k * m0_scaled) return {m, true};
                for (std::size_t i = 0; i < n; ++i) {
                    if (!(std::abs(f_trial[i]) <= 64.0 * std::numeric_limits<double>::epsilon() * scale[i])) {
                        return {m, false};
                    }
                }
                return {m, true};
            };
            double lambda = 1.0, best_lambda = 1.0;
            double best = std::numeric_limits<double>::infinity();
            bool accepted = false;
            for (int b = 0; b <= search.max_backtracks; ++b, lambda *= 0.5) {
                if (b > 0) {
                    for (std::size_t i = 0; i < n; ++i) x[i] = x0[i] + lambda * (x1[i] - x0[i]);
                }
                system.evaluate(x, f_trial, jacobian);
                const auto [figure, passes] = judge(lambda);
                if (passes) {
                    accepted = true;
                    taken = lambda;
                    f.swap(f_trial);  // F and J at the new iterate: the next iteration's
                    evaluated = true;
                    break;
                }
                if (figure < best) {
                    best = figure;
                    best_lambda = lambda;
                }
            }
            if (!accepted) {
                // No trial decreased enough: the one of smallest merit (the smallest trial when
                // none had a finite residual), never a zero step.
                if (!std::isfinite(best)) best_lambda = std::ldexp(1.0, -search.max_backtracks);
                taken = best_lambda;
                if (best_lambda == 1.0) {
                    std::copy(x1.begin(), x1.end(), x.begin());
                } else {
                    for (std::size_t i = 0; i < n; ++i) {
                        x[i] = x0[i] + best_lambda * (x1[i] - x0[i]);
                    }
                }
            }
        }
        record.iterations.push_back({iteration, update, residual, taken});
        record.converged = converged;
        if (observe) observe(record.iterations.back(), converged);
        if (converged) return {};
    }
    return std::unexpected(base::Error{
        base::ErrorCode::non_convergence,
        "Newton did not converge in " + std::to_string(options.max_iterations) + " iterations",
        base::ErrorContext{.index = static_cast<std::size_t>(options.max_iterations),
                           .value = record.iterations.back().update}});
}

}  // namespace NiTCAD::solve
