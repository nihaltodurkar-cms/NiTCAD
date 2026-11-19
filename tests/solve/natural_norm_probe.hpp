// Unit 59 slice 0 (DECISIONS.md ND1): the natural test's contraction along Newton's path, in
// today's measure and in each candidate norm. Test code only: it drives assemble::System through
// its public Newton hooks, as newton_solve does, and reproduces newton_solve's plain and
// natural-search iterations (checked bit for bit against it by the callers).
#pragma once

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/results/convergence.hpp"
#include "NiTCAD/solve/newton.hpp"

namespace natural_norm {

using NiTCAD::assemble::FieldId;
using NiTCAD::assemble::System;

// Which unknowns a norm measures: not a held one (its Jacobian row is its own diagonal alone: a
// contact's Dirichlet value, an electrode node, an insulator node's density), and the field of
// each (densities scaled against their field's largest value, every other field like psi).
struct Layout {
    std::vector<int> field;  // slot in the node block, -1 if held
    std::vector<bool> density;
    std::size_t block = 0;
    Layout(const System& s, const NiTCAD::linalg::SparseMatrix& j) {
        const auto& dofs = s.dofs();
        block = dofs.block();
        const std::size_t n = s.unknowns();
        field.assign(n, -1);
        density.assign(n, false);
        const auto offsets = j.row_offsets();
        const auto columns = j.col_indices();
        const auto values = j.values();
        for (std::size_t r = 0; r < n; ++r) {
            bool off_diagonal = false;
            for (auto k = static_cast<std::size_t>(offsets[r]);
                 k < static_cast<std::size_t>(offsets[r + 1]); ++k) {
                if (static_cast<std::size_t>(columns[k]) != r && values[k] != 0.0) {
                    off_diagonal = true;
                }
            }
            if (!off_diagonal) continue;
            field[r] = static_cast<int>(r % block);
            const FieldId id = dofs.field_of(r).id;
            density[r] = id == FieldId::electrons || id == FieldId::holes;
        }
    }
};

enum class Kind { weighted_rms, potential_only, asinh_rms };

// A candidate norm (ND3), measured at the reference state x_ref (the iterate) for both
// corrections. tau: the density threshold, relative to the field's largest value at x_ref.
struct Norm {
    std::string name;
    Kind kind;
    double tau;
    [[nodiscard]] double operator()(const Layout& l, std::span<const double> x_ref,
                                    std::span<const double> d) const {
        std::vector<double> largest(l.block, 0.0), sum(l.block, 0.0);
        std::vector<std::size_t> count(l.block, 0);
        for (std::size_t i = 0; i < d.size(); ++i) {
            if (l.field[i] >= 0 && l.density[i]) {
                auto& m = largest[static_cast<std::size_t>(l.field[i])];
                m = std::max(m, std::abs(x_ref[i]));
            }
        }
        for (std::size_t i = 0; i < d.size(); ++i) {
            if (l.field[i] < 0) continue;
            const auto f = static_cast<std::size_t>(l.field[i]);
            if (kind == Kind::potential_only && l.density[i]) continue;
            const double w = l.density[i] ? std::max(std::abs(x_ref[i]), tau * largest[f])
                                          : std::max(std::abs(x_ref[i]), 1.0);
            double u = d[i] / w;
            if (kind == Kind::asinh_rms && l.density[i]) u = std::asinh(u);
            sum[f] += u * u;
            ++count[f];
        }
        double total = 0.0;
        std::size_t fields = 0;
        for (std::size_t f = 0; f < l.block; ++f) {
            if (count[f] == 0) continue;
            total += sum[f] / static_cast<double>(count[f]);
            ++fields;
        }
        return fields == 0 ? 0.0 : std::sqrt(total / static_cast<double>(fields));
    }
};

inline std::vector<Norm> candidates() {
    std::vector<Norm> c;
    for (const double tau : {1e-4, 1e-8, 1e-12}) {
        char name[32];
        std::snprintf(name, sizeof name, "(a) rms %.0e", tau);
        c.push_back({name, Kind::weighted_rms, tau});
    }
    c.push_back({"(b) psi only", Kind::potential_only, 0.0});
    for (const double tau : {1e-4, 1e-8, 1e-12}) {
        char name[32];
        std::snprintf(name, sizeof name, "(c) asinh %.0e", tau);
        c.push_back({name, Kind::asinh_rms, tau});
    }
    return c;
}

// theta = |dx_bar| / |dx| at one plain iteration (lambda = 1), per measure: [0] today's (the
// simplified correction measured at the trial), [1] today's measure at the iterate, then the
// candidates in candidates() order.
struct Row {
    int iteration;
    double update;
    std::vector<double> theta;
};

struct Path {
    std::vector<Row> rows;
    bool converged = false;
    NiTCAD::results::ConvergenceRecord record;  // the iterations as newton_solve records them
};

inline void solve_or_fail(NiTCAD::linalg::LinearSolver& solver, std::span<const double> b,
                          std::span<double> x) {
    const auto r = solver.solve(b, x);
    REQUIRE(r.has_value());
}

// Plain Newton (no line search) from x, as newton_solve takes it, with theta in every measure.
inline Path plain_path(const System& s, std::vector<double> x,
                       const NiTCAD::solve::NewtonOptions& o,
                       const NiTCAD::linalg::SolverConfig& config = {}) {
    const std::size_t n = s.unknowns();
    auto solver = NiTCAD::linalg::LinearSolver::create(config);
    REQUIRE(solver.has_value());
    auto j = s.make_jacobian();
    std::vector<double> f(n), rhs(n), dx(n), f1(n), dx_bar(n);
    const auto norms = candidates();
    Path path;
    for (int it = 1; it <= o.max_iterations; ++it) {
        s.evaluate(x, f, j);
        const Layout layout(s, j);
        double residual = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            residual = std::max(residual, std::abs(f[i]));
            rhs[i] = -f[i];
        }
        REQUIRE(solver->factorize(j).has_value());
        solve_or_fail(*solver, rhs, dx);
        const double update = s.update_size(x, dx);
        const bool converged = update < o.tol_update;
        const std::vector<double> x0 = x;
        s.apply_update(x, dx, o.max_update);
        path.record.iterations.push_back({it, update, residual, 1.0});
        if (converged) {
            path.converged = true;
            break;
        }
        s.residual(x, f1);
        for (std::size_t i = 0; i < n; ++i) rhs[i] = -f1[i];
        solve_or_fail(*solver, rhs, dx_bar);
        Row row{it, update, {}};
        row.theta.push_back(s.update_size(x, dx_bar) / update);
        row.theta.push_back(s.update_size(x0, dx_bar) / update);
        for (const auto& norm : norms) {
            const double full = norm(layout, x0, dx);
            row.theta.push_back(full > 0.0 ? norm(layout, x0, dx_bar) / full : 0.0);
        }
        path.rows.push_back(std::move(row));
    }
    path.record.converged = path.converged;
    return path;
}

// newton_solve's natural line search with the measure m(x_ref, d) on both sides of the test,
// x_ref the iterate; with `today` the simplified correction is measured at the trial with
// update_size instead (Unit 30's test exactly). Returns whether it converged; x is the last
// iterate; record as newton_solve's.
inline bool natural_search(const System& s, std::vector<double>& x,
                           const NiTCAD::solve::NewtonOptions& o, const Norm* norm,
                           NiTCAD::results::ConvergenceRecord& record,
                           const NiTCAD::linalg::SolverConfig& config = {}) {
    record = {};
    const std::size_t n = s.unknowns();
    auto solver = NiTCAD::linalg::LinearSolver::create(config);
    REQUIRE(solver.has_value());
    auto j = s.make_jacobian();
    std::vector<double> f(n), rhs(n), dx(n), f_trial(n), dx_bar(n), x1(n);
    const int backtracks = o.line_search.max_backtracks;
    for (int it = 1; it <= o.max_iterations; ++it) {
        s.evaluate(x, f, j);
        const Layout layout(s, j);
        double residual = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(f[i])) return false;
            residual = std::max(residual, std::abs(f[i]));
            rhs[i] = -f[i];
        }
        if (!solver->factorize(j).has_value()) return false;
        if (!solver->solve(rhs, dx).has_value()) return false;
        const double update = s.update_size(x, dx);
        const bool converged = update < o.tol_update;
        double taken = 1.0;
        if (converged) {
            s.apply_update(x, dx, o.max_update);
        } else {
            const std::vector<double> x0 = x;
            s.apply_update(x, dx, o.max_update);
            x1 = x;
            const double full = norm != nullptr ? (*norm)(layout, x0, dx) : update;
            double lambda = 1.0, best_lambda = 1.0;
            double best = std::numeric_limits<double>::infinity();
            bool accepted = false;
            for (int b = 0; b <= backtracks; ++b, lambda *= 0.5) {
                if (b > 0) {
                    for (std::size_t i = 0; i < n; ++i) x[i] = x0[i] + lambda * (x1[i] - x0[i]);
                }
                s.evaluate(x, f_trial, j);
                double size = std::numeric_limits<double>::infinity();
                bool finite = true;
                for (std::size_t i = 0; i < n; ++i) {
                    if (!std::isfinite(f_trial[i])) finite = false;
                    rhs[i] = -f_trial[i];
                }
                if (finite && solver->solve(rhs, dx_bar).has_value()) {
                    size = norm != nullptr ? (*norm)(layout, x0, dx_bar) : s.update_size(x, dx_bar);
                }
                if (size <= (1.0 - 0.25 * lambda) * full) {
                    accepted = true;
                    taken = lambda;
                    break;
                }
                if (size < best) {
                    best = size;
                    best_lambda = lambda;
                }
            }
            if (!accepted) {
                if (!std::isfinite(best)) best_lambda = std::ldexp(1.0, -backtracks);
                taken = best_lambda;
                if (best_lambda == 1.0) {
                    x = x1;
                } else {
                    for (std::size_t i = 0; i < n; ++i) {
                        x[i] = x0[i] + best_lambda * (x1[i] - x0[i]);
                    }
                }
            }
        }
        record.iterations.push_back({it, update, residual, taken});
        record.converged = converged;
        if (converged) return true;
    }
    return false;
}

// The same iterations, update for update and step for step.
inline bool same_iterations(const NiTCAD::results::ConvergenceRecord& a,
                            const NiTCAD::results::ConvergenceRecord& b) {
    if (a.iterations.size() != b.iterations.size()) return false;
    for (std::size_t k = 0; k < a.iterations.size(); ++k) {
        if (a.iterations[k].update != b.iterations[k].update ||
            a.iterations[k].step != b.iterations[k].step) {
            return false;
        }
    }
    return true;
}

// The largest relative difference between two states (each component against the larger
// magnitude of the two).
inline double largest_relative(std::span<const double> u, std::span<const double> v) {
    REQUIRE(u.size() == v.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        worst = std::max(worst, std::abs(u[i] - v[i]) /
                                    std::max({std::abs(u[i]), std::abs(v[i]), 1e-300}));
    }
    return worst;
}

// The largest difference between two states as the candidate norms weigh a component: psi
// against max(|psi|, 1), a density against max(|n|, 1e-8 of its field's largest).
inline double largest_weighted(const Layout& l, std::span<const double> u,
                               std::span<const double> v) {
    std::vector<double> largest(l.block, 0.0);
    for (std::size_t i = 0; i < u.size(); ++i) {
        if (l.field[i] >= 0 && l.density[i]) {
            auto& m = largest[static_cast<std::size_t>(l.field[i])];
            m = std::max(m, std::abs(v[i]));
        }
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        if (l.field[i] < 0) continue;
        const double w = l.density[i]
                             ? std::max(std::abs(v[i]),
                                        1e-8 * largest[static_cast<std::size_t>(l.field[i])])
                             : std::max(std::abs(v[i]), 1.0);
        worst = std::max(worst, std::abs(u[i] - v[i]) / w);
    }
    return worst;
}

// Prints a path's theta table and, per measure, how many of plain Newton's steps pass the
// natural test at lambda = 1 (theta <= 3/4), how many contract (theta < 1), and the largest theta.
inline void print_path(const char* label, const Path& p) {
    const auto norms = candidates();
    std::vector<std::string> names{"today (trial)", "today (iterate)"};
    for (const auto& c : norms) names.push_back(c.name);
    std::printf("PATH %s: %s in %zu iterations\n", label, p.converged ? "converged" : "FAILED",
                p.record.iterations.size());
    std::printf("  it    update   ");
    for (const auto& nm : names) std::printf(" %15s", nm.c_str());
    std::printf("\n");
    for (const auto& r : p.rows) {
        std::printf("  %3d  %.2e ", r.iteration, r.update);
        for (const double t : r.theta) std::printf(" %15.3e", t);
        std::printf("\n");
    }
    std::printf("  SUMMARY %s (of %zu non-final steps)\n", label, p.rows.size());
    for (std::size_t m = 0; m < names.size(); ++m) {
        int pass = 0, contract = 0;
        double worst = 0.0;
        for (const auto& r : p.rows) {
            pass += r.theta[m] <= 0.75 ? 1 : 0;
            contract += r.theta[m] < 1.0 ? 1 : 0;
            worst = std::max(worst, r.theta[m]);
        }
        std::printf("    %-16s pass %3d  contract %3d  largest %.3e\n", names[m].c_str(), pass,
                    contract, worst);
    }
}

}  // namespace natural_norm
