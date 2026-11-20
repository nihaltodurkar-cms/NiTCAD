// Unit 63 slice 0 (MD1): GMRES on one field's diagonal block of a coupled Jacobian, preconditioned by
// the slice-2 AMG on the block as it is, or by the AMG on its symmetrization (SymmetrizedAmg, MD2).
// Test-only; shared by the linalg fixture probe and the deck layer's 3D diode probe.
#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "../../src/linalg/amg.hpp"
#include "../../src/linalg/gmres.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace block_probe {

using NiTCAD::linalg::Index;
using NiTCAD::linalg::SparseMatrix;
using NiTCAD::linalg::Triplet;
namespace detail = NiTCAD::linalg::detail;

// The diagonal block of the unknowns u with u % stride == field (node-interleaved layout).
inline SparseMatrix block(const SparseMatrix& a, std::size_t stride, std::size_t field) {
    const auto n = static_cast<std::size_t>(a.rows());
    std::vector<Triplet> t;
    for (std::size_t r = field; r < n; r += stride) {
        for (Index k = a.row_offsets()[r]; k < a.row_offsets()[r + 1]; ++k) {
            const auto c = static_cast<std::size_t>(a.col_indices()[k]);
            if (c % stride != field) continue;
            t.push_back({static_cast<Index>(r / stride), static_cast<Index>(c / stride), a.values()[k]});
        }
    }
    const auto m = static_cast<Index>((n - field + stride - 1) / stride);
    return SparseMatrix::from_triplets(m, m, t).value();
}

inline detail::Csr to_csr(const SparseMatrix& a) {
    detail::Csr c;
    c.rows = c.cols = static_cast<std::size_t>(a.rows());
    c.offsets.assign(a.row_offsets().begin(), a.row_offsets().end());
    c.columns.assign(a.col_indices().begin(), a.col_indices().end());
    c.values.assign(a.values().begin(), a.values().end());
    return c;
}

// Componentwise backward error of y; r = c - A y.
inline double backward_error(const SparseMatrix& a, std::span<const double> c, std::span<const double> y,
                             std::span<double> r) {
    double worst = 0.0;
    for (Index i = 0; i < a.rows(); ++i) {
        double ax = 0.0;
        double abs_ax = 0.0;
        for (Index k = a.row_offsets()[i]; k < a.row_offsets()[i + 1]; ++k) {
            const double term = a.values()[k] * y[static_cast<std::size_t>(a.col_indices()[k])];
            ax += term;
            abs_ax += std::abs(term);
        }
        const auto u = static_cast<std::size_t>(i);
        r[u] = c[u] - ax;
        const double denominator = abs_ax + std::abs(c[u]);
        worst = std::max(worst, denominator > 0.0 ? std::abs(r[u]) / denominator : 0.0);
    }
    return worst;
}

struct Outcome {
    bool built = false;
    bool converged = false;
    std::size_t iterations = 0;
    double setup_s = 0.0;
    double solve_s = 0.0;
    std::size_t levels = 0;
    std::size_t coarsest = 0;
};

inline detail::AmgSettings default_amg() {
    const NiTCAD::linalg::SolverConfig c{};
    return {.strength_threshold = c.amg_strength_threshold,
            .smoother_degree = c.amg_smoother_degree,
            .coarse_size = static_cast<std::size_t>(c.amg_coarse_size),
            .max_levels = c.amg_max_levels};
}

// GMRES(30), at most 1000 products, to a componentwise backward error of 1e-8, from y = 0.
template <class Preconditioner>
Outcome run(const SparseMatrix& a, std::span<const double> c, Preconditioner& m) {
    using clock = std::chrono::steady_clock;
    Outcome out;
    const auto t0 = clock::now();
    if (!m.analyze(a) || !m.factorize(a.values())) return out;
    out.built = true;
    const auto t1 = clock::now();
    std::vector<double> y(c.size());
    detail::Gmres<double> gmres;
    const auto result = gmres.solve(
        {.restart = 30, .max_iterations = 1000, .tolerance = 1e-8},
        [&](std::span<const double> v, std::span<double> w) { a.multiply(v, w); },
        [&](std::span<const double> v, std::span<double> w) { m.solve(v, w); },
        [&](std::span<const double> yy, std::span<double> r) { return backward_error(a, c, yy, r); },
        std::span<double>(y));
    const auto t2 = clock::now();
    out.converged = result.converged;
    out.iterations = result.iterations;
    out.setup_s = std::chrono::duration<double>(t1 - t0).count();
    out.solve_s = std::chrono::duration<double>(t2 - t1).count();
    return out;
}

inline std::vector<double> rhs(std::size_t n) {
    std::vector<double> b(n);
    for (std::size_t i = 0; i < n; ++i) b[i] = std::sin(0.37 * static_cast<double>(i)) + 1e-3;
    return b;
}

inline void print(const char* what, const Outcome& o) {
    if (!o.built) {
        std::printf("    %-22s not built\n", what);
        return;
    }
    std::printf("    %-22s %s%4zu its  setup %.3f s  solve %.3f s  levels %zu (coarsest %zu)\n", what,
                o.converged ? "" : "FAILED ", o.iterations, o.setup_s, o.solve_s, o.levels, o.coarsest);
}

// Both preconditioners on one block; returns {plain, symmetrized}.
inline std::pair<Outcome, Outcome> compare(const std::string& name, const SparseMatrix& b) {
    const std::vector<double> c = rhs(static_cast<std::size_t>(b.rows()));
    detail::AmgBackend plain(default_amg(), nullptr);
    Outcome p = run(b, c, plain);
    if (p.built) {
        const auto sizes = plain.level_sizes();
        p.levels = sizes.size();
        p.coarsest = sizes.back();
    }
    detail::SymmetrizedAmg sym(default_amg(), nullptr);
    Outcome s = run(b, c, sym);
    if (s.built) {
        const auto sizes = sym.amg().level_sizes();
        s.levels = sizes.size();
        s.coarsest = sizes.back();
    }
    std::printf("  %s: %lld unknowns, %zu nonzeros", name.c_str(), static_cast<long long>(b.rows()),
                b.nonzeros());
    if (s.built) {
        std::printf("  (components %zu, unpaired %zu, clamped %zu)", sym.symmetrization().components,
                    sym.symmetrization().unpaired, sym.symmetrization().clamped);
    }
    std::printf("\n");
    print("AMG on the block", p);
    print("AMG symmetrized", s);
    return {p, s};
}

}  // namespace block_probe
