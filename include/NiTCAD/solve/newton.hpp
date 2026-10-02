// Damped Newton iteration in scaled variables (ARCHITECTURE.md 6.2; legacy device.py solver
// options and Device1D::solve_equilibrium_boltzmann).
//
// Each iteration evaluates F and J at x, solves J dx = -F with the given LinearSolver, clips every
// component of dx to +-max_update (the legacy max_dpsi damping cap) and adds it to x.
// Convergence (decided at Unit 8, the question left open in 6.1): the iteration has converged when
// the largest component of the FULL correction dx, before clipping, is below tol_update. The legacy
// measured the clipped correction; with a componentwise clip and tol_update < max_update the two
// tests agree, but judging the full correction does not rely on that. The legacy 1D loop declares a
// residual tolerance (1e-10) and never uses it; there is none here, and the final residual is
// reported instead.
// The Jacobian pattern is fixed by the system, so the solver analyzes once per run.
#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <expected>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace NiTCAD::solve {

struct NewtonOptions {
    int max_iterations = 100;
    double tol_update = 1e-8;  // on the largest scaled correction
    double max_update = 5.0;   // componentwise damping cap on the scaled correction
};

struct NewtonReport {
    int iterations = 0;
    std::vector<double> updates;  // largest full correction of each iteration
    double final_residual = 0.0;  // largest |F| at the last evaluated x (before its correction)
};

// A system Newton can drive: a fixed number of unknowns, a Jacobian pattern, and F and J at x.
template <class S>
concept NewtonSystem = requires(const S& s, std::span<const double> x, std::span<double> f,
                                linalg::SparseMatrix& j) {
    { s.unknowns() } -> std::convertible_to<std::size_t>;
    { s.make_jacobian() } -> std::same_as<linalg::SparseMatrix>;
    s.evaluate(x, f, j);
};

// Errors:
// - invalid_input: max_iterations < 1, or tol_update or max_update not finite and positive;
// - non_convergence: no convergence within max_iterations (context: the iteration count and the
//   last full correction), or a non-finite residual (context: the iteration);
// - singular_system, inaccurate_solve, resource_exhausted: from the linear solver, with the
//   iteration added to the message.
// x is updated in place; on error it holds the last iterate. Precondition (NITCAD_EXPECTS via the
// system): x has system.unknowns() entries.
template <NewtonSystem System>
[[nodiscard]] std::expected<NewtonReport, base::Error> newton_solve(
    const System& system, std::span<double> x, const NewtonOptions& options,
    linalg::LinearSolver& solver) {
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

    const std::size_t n = system.unknowns();
    linalg::SparseMatrix jacobian = system.make_jacobian();
    std::vector<double> f(n), rhs(n), dx(n);
    NewtonReport report;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        system.evaluate(x, f, jacobian);
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
        report.final_residual = residual;
        auto factored = solver.factorize(jacobian);
        if (!factored) {
            base::Error e = std::move(factored.error());
            e.message = "Newton iteration " + std::to_string(iteration) + ": " + e.message;
            return std::unexpected(std::move(e));
        }
        auto solved = solver.solve(rhs, dx);
        if (!solved) {
            base::Error e = std::move(solved.error());
            e.message = "Newton iteration " + std::to_string(iteration) + ": " + e.message;
            return std::unexpected(std::move(e));
        }
        double update = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            update = std::max(update, std::abs(dx[i]));
            x[i] += std::clamp(dx[i], -options.max_update, options.max_update);
        }
        report.iterations = iteration;
        report.updates.push_back(update);
        if (update < options.tol_update) return report;
    }
    return std::unexpected(base::Error{
        base::ErrorCode::non_convergence,
        "Newton did not converge in " + std::to_string(options.max_iterations) + " iterations",
        base::ErrorContext{.index = static_cast<std::size_t>(options.max_iterations),
                           .value = report.updates.back()}});
}

}  // namespace NiTCAD::solve
