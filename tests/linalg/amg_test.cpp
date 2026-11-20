// The amg preconditioner of the gmres backend (Unit 58 SD3, slice 2): its kernels (transpose,
// product, coarsening) against dense references, the Galerkin operator, the configuration checks,
// and GMRES-AMG through the public interface: accuracy against an independent backward error,
// iterations independent of the mesh size on Poisson problems, the near-null vector under
// equilibration, reuse of the coarsening, the fallback, and bit-identical results for every thread
// count (SG2).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "../../src/linalg/amg.hpp"
#include "NiTCAD/base/thread_pool.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "fixtures.hpp"

using namespace NiTCAD::linalg;
using NiTCAD::base::ErrorCode;
using namespace fixtures;
namespace amg = NiTCAD::linalg::detail::amg;
using NiTCAD::linalg::detail::Csr;

namespace {

constexpr SolverConfig amg_strict{.backend = SolverBackend::gmres_ilut,
                                  .fallback = LinearFallback::none,
                                  .preconditioner = Preconditioner::amg};

Csr to_csr(const SparseMatrix& a) {
    Csr c;
    c.rows = static_cast<std::size_t>(a.rows());
    c.cols = static_cast<std::size_t>(a.cols());
    c.offsets.assign(a.row_offsets().begin(), a.row_offsets().end());
    c.columns.assign(a.col_indices().begin(), a.col_indices().end());
    c.values.assign(a.values().begin(), a.values().end());
    return c;
}

std::vector<double> dense(const Csr& a) {
    std::vector<double> d(a.rows * a.cols, 0.0);
    for (std::size_t i = 0; i < a.rows; ++i) {
        for (auto k = static_cast<std::size_t>(a.offsets[i]); k < static_cast<std::size_t>(a.offsets[i + 1]); ++k) {
            d[i * a.cols + static_cast<std::size_t>(a.columns[k])] += a.values[k];
        }
    }
    return d;
}

bool sorted_rows(const Csr& a) {
    for (std::size_t i = 0; i < a.rows; ++i) {
        for (auto k = static_cast<std::size_t>(a.offsets[i]) + 1; k < static_cast<std::size_t>(a.offsets[i + 1]); ++k) {
            if (!(a.columns[k - 1] < a.columns[k])) return false;
        }
    }
    return true;
}

Csr random_csr(std::size_t rows, std::size_t cols, double density, std::uint64_t seed) {
    Csr c;
    c.rows = rows;
    c.cols = cols;
    c.offsets.assign(rows + 1, 0);
    std::uint64_t k = seed;
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            if (hash01(k++) < density) {
                c.columns.push_back(static_cast<Index>(j));
                c.values.push_back(hash01(k++) - 0.5);
            }
        }
        c.offsets[i + 1] = static_cast<Index>(c.columns.size());
    }
    return c;
}

// The 2D five-point Laplacian on an m x m grid whose boundary rows are identity (Dirichlet) rows
// coupled into the interior ones, as a device Jacobian's contacts are: nonsymmetric, with isolated
// rows. conv adds upwind convection along x.
SparseMatrix poisson_2d(Index m, double conv = 0.0) {
    std::vector<Triplet> t;
    for (Index j = 0; j < m; ++j) {
        for (Index i = 0; i < m; ++i) {
            const Index k = j * m + i;
            if (i == 0 || j == 0 || i + 1 == m || j + 1 == m) {
                t.push_back({k, k, 1.0});
                continue;
            }
            t.push_back({k, k, 4.0 + conv});
            t.push_back({k, k - 1, -1.0 - conv});
            t.push_back({k, k + 1, -1.0});
            t.push_back({k, k - m, -1.0});
            t.push_back({k, k + m, -1.0});
        }
    }
    return from(m * m, t);
}

// The 3D seven-point Laplacian with identity boundary rows.
SparseMatrix poisson_3d(Index m) {
    std::vector<Triplet> t;
    for (Index z = 0; z < m; ++z) {
        for (Index y = 0; y < m; ++y) {
            for (Index x = 0; x < m; ++x) {
                const Index k = (z * m + y) * m + x;
                if (x == 0 || y == 0 || z == 0 || x + 1 == m || y + 1 == m || z + 1 == m) {
                    t.push_back({k, k, 1.0});
                    continue;
                }
                t.push_back({k, k, 6.0});
                for (const Index d : {Index{1}, m, m * m}) {
                    t.push_back({k, k - d, -1.0});
                    t.push_back({k, k + d, -1.0});
                }
            }
        }
    }
    return from(m * m * m, t);
}

std::vector<double> rhs(std::size_t n) {
    std::vector<double> b(n);
    for (std::size_t i = 0; i < n; ++i) b[i] = std::sin(0.37 * static_cast<double>(i)) + 1e-3;
    return b;
}

// Componentwise backward error, written independently of the solver's.
double backward_error(const SparseMatrix& a, std::span<const double> b, std::span<const double> x) {
    double worst = 0.0;
    for (Index i = 0; i < a.rows(); ++i) {
        double ax = 0.0;
        double abs_ax = 0.0;
        for (Index k = a.row_offsets()[i]; k < a.row_offsets()[i + 1]; ++k) {
            const double term = a.values()[k] * x[static_cast<std::size_t>(a.col_indices()[k])];
            ax += term;
            abs_ax += std::abs(term);
        }
        const auto u = static_cast<std::size_t>(i);
        const double denominator = abs_ax + std::abs(b[u]);
        worst = std::max(worst, denominator > 0.0 ? std::abs(b[u] - ax) / denominator : 0.0);
    }
    return worst;
}

struct Solved {
    std::vector<double> x;
    SolveReport report;
    FactorizationReport factorization;
};

Solved solve_with(const SparseMatrix& a, std::span<const double> b, const SolverConfig& config) {
    auto solver = LinearSolver::create(config).value();
    const auto f = solver.factorize(a);
    INFO("factorize: " << (f ? std::string("ok") : f.error().message));
    REQUIRE(f.has_value());
    Solved s{std::vector<double>(b.size()), {}, *f};
    const auto r = solver.solve(b, s.x);
    INFO("solve: " << (r ? std::string("ok") : r.error().message));
    REQUIRE(r.has_value());
    s.report = *r;
    REQUIRE(backward_error(a, b, s.x) <= config.max_backward_error);
    return s;
}

}  // namespace

// --- kernels ---------------------------------------------------------------------------------

TEST_CASE("amg kernels: transpose and product agree with dense arithmetic") {
    const Csr a = random_csr(37, 23, 0.2, 11);
    const Csr b = random_csr(23, 41, 0.15, 97);
    const Csr at = amg::transpose(a);
    REQUIRE(at.rows == 23);
    REQUIRE(at.cols == 37);
    REQUIRE(sorted_rows(at));
    const std::vector<double> da = dense(a);
    const std::vector<double> dat = dense(at);
    for (std::size_t i = 0; i < 37; ++i) {
        for (std::size_t j = 0; j < 23; ++j) REQUIRE(dat[j * 37 + i] == da[i * 23 + j]);
    }
    const Csr c = amg::product(a, b);
    REQUIRE(c.rows == 37);
    REQUIRE(c.cols == 41);
    REQUIRE(sorted_rows(c));
    const std::vector<double> db = dense(b);
    const std::vector<double> dc = dense(c);
    for (std::size_t i = 0; i < 37; ++i) {
        for (std::size_t j = 0; j < 41; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < 23; ++k) s += da[i * 23 + k] * db[k * 41 + j];
            REQUIRE(std::abs(dc[i * 41 + j] - s) <= 1e-15);
        }
    }
    // y = A x, serial and on a pool of 3 threads (several row chunks): the same bits as dense.
    {
        const Csr big = random_csr(5000, 300, 0.02, 41);
        std::vector<double> x(300);
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = hash01(7 * i) - 0.5;
        std::vector<double> serial(5000), pooled(5000);
        amg::multiply(nullptr, big, x, serial);
        NiTCAD::base::ThreadPool pool(3);
        amg::multiply(&pool, big, x, pooled);
        REQUIRE(serial == pooled);
        for (std::size_t i = 0; i < big.rows; ++i) {
            double s = 0.0;
            for (auto k = static_cast<std::size_t>(big.offsets[i]); k < static_cast<std::size_t>(big.offsets[i + 1]); ++k) {
                s += big.values[k] * x[static_cast<std::size_t>(big.columns[k])];
            }
            REQUIRE(serial[i] == s);
        }
    }
    // The product's pattern is structural: a cancelling entry stays.
    Csr x{.rows = 1, .cols = 2, .offsets = {0, 2}, .columns = {0, 1}, .values = {1.0, -1.0}};
    Csr y{.rows = 2, .cols = 1, .offsets = {0, 1, 2}, .columns = {0, 0}, .values = {2.0, 2.0}};
    const Csr z = amg::product(x, y);
    REQUIRE(z.nonzeros() == 1);
    REQUIRE(z.values[0] == 0.0);
}

TEST_CASE("amg kernels: coarsening of a 1D chain, isolated Dirichlet rows, determinism") {
    // Interior chain 1..n-2 with Dirichlet rows at both ends.
    const Index n = 32;
    std::vector<Triplet> t{{0, 0, 1.0}, {n - 1, n - 1, 1.0}};
    for (Index i = 1; i + 1 < n; ++i) {
        t.push_back({i, i - 1, -1.0});
        t.push_back({i, i, 2.0});
        t.push_back({i, i + 1, -1.0});
    }
    const Csr a = to_csr(from(n, t));
    const amg::Coarsening c = amg::coarsen(a, 0.08);
    REQUIRE(c.aggregate[0] == -1);
    REQUIRE(c.aggregate[static_cast<std::size_t>(n - 1)] == -1);
    // Pass 1 in node order: {1, 2}, then {3, 4, 5}, {6, 7, 8}, ... (node 1's neighbours are 2 only).
    REQUIRE(c.aggregate[1] == 0);
    REQUIRE(c.aggregate[2] == 0);
    REQUIRE(c.aggregate[3] == 1);
    REQUIRE(c.aggregate[4] == 1);
    REQUIRE(c.aggregate[5] == 1);
    REQUIRE(c.aggregate[6] == 2);
    for (std::size_t i = 1; i + 1 < static_cast<std::size_t>(n); ++i) {
        REQUIRE(c.aggregate[i] >= 0);
        REQUIRE(static_cast<std::size_t>(c.aggregate[i]) < c.aggregates);
    }
    REQUIRE(c.aggregates == 10);  // {1, 2}, nine of three, node 30 joins its neighbour's
    // The Dirichlet couplings (columns 0 and n-1 in rows 1 and n-2) are not strong.
    REQUIRE(c.strong[static_cast<std::size_t>(a.offsets[1])] == 0);
    // Dirichlet rows that keep their off-diagonal entries as explicit zeros (as the device
    // Jacobians' contact rows do) are isolated too.
    {
        std::vector<Triplet> z = t;
        z.push_back({0, 1, 0.0});
        z.push_back({n - 1, n - 2, 0.0});
        const Csr az = to_csr(from(n, z));
        REQUIRE(az.nonzeros() == a.nonzeros() + 2);
        const amg::Coarsening cz = amg::coarsen(az, 0.08);
        REQUIRE(cz.aggregate == c.aggregate);
    }
    const amg::Coarsening again = amg::coarsen(a, 0.08);
    REQUIRE(again.aggregate == c.aggregate);
    REQUIRE(again.strong == c.strong);
}

TEST_CASE("amg kernels: weak couplings are not strong, and strength is symmetric") {
    // Row 0 couples to 1 strongly; row 2 couples to 1 only weakly in its own row, but row 1 couples
    // to 2 strongly, so 1-2 is strong both ways.
    const SparseMatrix m = from(3, {{0, 0, 2.0}, {0, 1, -1.0}, {1, 0, -1.0}, {1, 1, 2.0},
                                    {1, 2, -1.0}, {2, 1, -1e-3}, {2, 2, 2.0}});
    const Csr a = to_csr(m);
    const amg::Coarsening c = amg::coarsen(a, 0.08);
    REQUIRE(c.strong == std::vector<std::uint8_t>{0, 1, 1, 0, 1, 1, 0});
    const amg::Coarsening loose = amg::coarsen(a, 0.0);
    REQUIRE(loose.strong == c.strong);
    // A large theta leaves no strong link: singletons.
    const SparseMatrix weak = from(2, {{0, 0, 4.0}, {0, 1, -1.0}, {1, 0, -1.0}, {1, 1, 4.0}});
    REQUIRE(amg::coarsen(to_csr(weak), 0.5).aggregates == 2);
    REQUIRE(amg::coarsen(to_csr(weak), 0.2).aggregates == 1);
}

TEST_CASE("amg: the coarse operator is the Galerkin product R A P with R = P^T") {
    const SparseMatrix a = poisson_2d(24, 0.5);
    detail::AmgBackend backend({.strength_threshold = 0.08, .smoother_degree = 2, .coarse_size = 10,
                                .max_levels = 3},
                               nullptr);
    REQUIRE(backend.analyze(a).has_value());
    REQUIRE(backend.factorize(a.values()).has_value());
    const std::vector<std::size_t> sizes = backend.level_sizes();
    REQUIRE(sizes.size() == 3);
    REQUIRE(sizes[0] == 576);
    REQUIRE(sizes[1] < sizes[0] / 4);
    REQUIRE(sizes[2] < sizes[1]);
    // R is P^T, and each coarse operator is P^T A P, by dense arithmetic.
    for (std::size_t l = 0; l + 1 < sizes.size(); ++l) {
        CAPTURE(l);
        const std::vector<double> da = dense(backend.level_matrix(l));
        const std::vector<double> dp = dense(backend.prolongator(l));
        const std::vector<double> dr = dense(backend.restrictor(l));
        const std::vector<double> dc = dense(backend.level_matrix(l + 1));
        const std::size_t n = sizes[l];
        const std::size_t nc = sizes[l + 1];
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < nc; ++j) REQUIRE(dr[j * n + i] == dp[i * nc + j]);
        }
        std::vector<double> ap(n * nc, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < n; ++k) {
                if (da[i * n + k] == 0.0) continue;
                for (std::size_t j = 0; j < nc; ++j) ap[i * nc + j] += da[i * n + k] * dp[k * nc + j];
            }
        }
        double largest = 0.0;
        for (const double v : dc) largest = std::max(largest, std::abs(v));
        for (std::size_t i = 0; i < nc; ++i) {
            for (std::size_t j = 0; j < nc; ++j) {
                double s = 0.0;
                for (std::size_t k = 0; k < n; ++k) s += dp[k * nc + i] * ap[k * nc + j];
                REQUIRE(std::abs(dc[i * nc + j] - s) <= 1e-13 * largest);
            }
        }
    }
}

// --- configuration ---------------------------------------------------------------------------

TEST_CASE("amg: default settings and the configuration checks") {
    const SolverConfig config;
    REQUIRE(config.preconditioner == Preconditioner::ilut);
    REQUIRE(config.amg_strength_threshold == 0.08);
    REQUIRE(config.amg_smoother_degree == 2);
    REQUIRE(config.amg_coarse_size == 1000);
    REQUIRE(config.amg_max_levels == 10);
    REQUIRE(LinearSolver::create(amg_strict).has_value());
    const auto invalid = [](SolverConfig c) {
        const auto s = LinearSolver::create(c);
        return !s && s.error().code == ErrorCode::invalid_input;
    };
    const auto with = [](auto set) {
        SolverConfig c = amg_strict;
        set(c);
        return c;
    };
    REQUIRE(invalid(with([](SolverConfig& c) { c.backend = SolverBackend::eigen_sparse_lu; })));
    REQUIRE(invalid(with([](SolverConfig& c) { c.preconditioner = static_cast<Preconditioner>(7); })));
    for (const double theta : {-0.1, 1.0, std::nan("")}) {
        REQUIRE(invalid(with([&](SolverConfig& c) { c.amg_strength_threshold = theta; })));
    }
    for (const int d : {0, 11}) REQUIRE(invalid(with([&](SolverConfig& c) { c.amg_smoother_degree = d; })));
    for (const int s : {0, 100001}) REQUIRE(invalid(with([&](SolverConfig& c) { c.amg_coarse_size = s; })));
    for (const int l : {1, 51}) REQUIRE(invalid(with([&](SolverConfig& c) { c.amg_max_levels = l; })));
    // The amg settings are validated for every backend, as the gmres ones are.
    REQUIRE(invalid(SolverConfig{.amg_max_levels = 1}));
    // Real only (small signal stays on SparseLU).
    const auto complex_solver = ComplexLinearSolver::create(amg_strict);
    REQUIRE(!complex_solver);
    REQUIRE(complex_solver.error().code == ErrorCode::invalid_input);
    REQUIRE(ComplexLinearSolver::create(SolverConfig{.backend = SolverBackend::gmres_ilut}).has_value());
    // The bounds themselves are accepted.
    REQUIRE(LinearSolver::create(with([](SolverConfig& c) {
                c.amg_strength_threshold = 0.0;
                c.amg_smoother_degree = 10;
                c.amg_coarse_size = 1;
                c.amg_max_levels = 50;
            })).has_value());
}

// --- solves ----------------------------------------------------------------------------------

TEST_CASE("amg: GMRES-AMG solves Poisson problems in few iterations, independent of the size") {
    std::vector<std::size_t> its;
    for (const Index m : {33, 65, 129, 257}) {
        CAPTURE(m);
        const SparseMatrix a = poisson_2d(m);
        const std::vector<double> b = rhs(static_cast<std::size_t>(m * m));
        const Solved s = solve_with(a, b, amg_strict);
        REQUIRE(!s.report.fell_back);
        its.push_back(s.report.iterations.value_or(0));
        std::printf("amg 2D %d^2: %zu iterations, fill %zu\n", m, its.back(),
                    s.factorization.fill.value_or(0));
    }
    REQUIRE(its.back() <= 20);
    REQUIRE(static_cast<double>(its.back()) <= 1.5 * static_cast<double>(its.front()));
    std::vector<std::size_t> its3;
    for (const Index m : {16, 32, 48}) {
        CAPTURE(m);
        const SparseMatrix a = poisson_3d(m);
        const std::vector<double> b = rhs(static_cast<std::size_t>(m * m * m));
        const Solved s = solve_with(a, b, amg_strict);
        its3.push_back(s.report.iterations.value_or(0));
        std::printf("amg 3D %d^3: %zu iterations\n", m, its3.back());
    }
    REQUIRE(static_cast<double>(its3.back()) <= 1.5 * static_cast<double>(its3.front()));
}

TEST_CASE("amg: the near-null vector follows the equilibration") {
    // A heterogeneous Poisson problem (box-method conductances 2^-10 to 2^10 in 8 x 8 blocks): its
    // unknowns' near-null vector is the constant, but the equilibrated matrix's columns are scaled
    // by varying powers of two, so on it the near-null vector is 1 / C.
    const Index m = 129;
    const auto kappa = [&](Index i, Index j) {
        const std::uint64_t block = static_cast<std::uint64_t>((j / 8) * 64 + i / 8);
        return std::ldexp(1.0, static_cast<int>(hash01(block) * 21.0) - 10);
    };
    std::vector<Triplet> t;
    for (Index j = 0; j < m; ++j) {
        for (Index i = 0; i < m; ++i) {
            const Index k = j * m + i;
            if (i == 0 || j == 0 || i + 1 == m || j + 1 == m) {
                t.push_back({k, k, 1.0});
                continue;
            }
            double diagonal = 0.0;
            for (const auto [di, dj] : {std::pair{-1, 0}, std::pair{1, 0}, std::pair{0, -1}, std::pair{0, 1}}) {
                const double g = std::min(kappa(i, j), kappa(i + di, j + dj));
                t.push_back({k, (j + dj) * m + i + di, -g});
                diagonal += g;
            }
            t.push_back({k, k, diagonal});
        }
    }
    const SparseMatrix a = from(m * m, t);
    const std::vector<double> b = rhs(static_cast<std::size_t>(m * m));
    const Solved homogeneous = solve_with(poisson_2d(m), b, amg_strict);
    const Solved s = solve_with(a, b, amg_strict);
    std::printf("amg heterogeneous 2D %d^2: %zu iterations (homogeneous %zu)\n", m,
                s.report.iterations.value_or(0), homogeneous.report.iterations.value_or(0));
    REQUIRE(s.report.iterations.value_or(0) <= 2 * homogeneous.report.iterations.value_or(0));
}

TEST_CASE("amg: convection-diffusion and the shared nonsymmetric systems") {
    // Moderate upwind convection converges; strong convection (cell Peclet above 2) is the known hard
    // case for whole-system AMG (DECISIONS.md Unit 58, slice 2): it may need the fallback.
    {
        const SparseMatrix a = poisson_2d(129, 0.5);
        const std::vector<double> b = rhs(129 * 129);
        const Solved s = solve_with(a, b, amg_strict);
        std::printf("amg convection 0.5: %zu iterations\n", s.report.iterations.value_or(0));
        REQUIRE(s.report.iterations.value_or(0) <= 40);
    }
    {
        const SparseMatrix a = poisson_2d(129, 2.0);
        const std::vector<double> b = rhs(129 * 129);
        SolverConfig c = amg_strict;
        c.fallback = LinearFallback::sparse_lu;
        const Solved s = solve_with(a, b, c);
        std::printf("amg convection 2.0: %zu iterations, fell back %d\n",
                    s.report.iterations.value_or(0), s.report.fell_back ? 1 : 0);
    }
    const System ns = nonsymmetric(3000);
    const Solved s = solve_with(ns.a, ns.b, amg_strict);
    REQUIRE(relative_error(s.x, ns.x_exact) <= 1e-6);
    // The coupled, badly scaled system: solved within the backward error, by AMG or by the fallback.
    const System coupled = badly_scaled_coupled(20);
    SolverConfig with_fallback = amg_strict;
    with_fallback.fallback = LinearFallback::sparse_lu;
    const Solved c = solve_with(coupled.a, coupled.b, with_fallback);
    std::printf("amg coupled: %zu iterations, fell back %d\n", c.report.iterations.value_or(0),
                c.report.fell_back ? 1 : 0);
}

TEST_CASE("amg: smoother degrees and level limits all converge") {
    const SparseMatrix a = poisson_2d(65, 0.3);
    const std::vector<double> b = rhs(65 * 65);
    for (const int degree : {1, 2, 3, 5}) {
        SolverConfig c = amg_strict;
        c.amg_smoother_degree = degree;
        CAPTURE(degree);
        (void)solve_with(a, b, c);
    }
    for (const int levels : {2, 3}) {
        SolverConfig c = amg_strict;
        c.amg_max_levels = levels;
        c.amg_coarse_size = 1;
        CAPTURE(levels);
        (void)solve_with(a, b, c);
    }
    // A strongly diagonally dominant matrix above coarse_size has no strong link: the coarsening
    // stalls at once and the only level is smoothed, not factorized (no stored hierarchy).
    {
        std::vector<Triplet> t;
        const Index m2 = 80;
        for (Index k = 0; k < m2 * m2; ++k) {
            t.push_back({k, k, 40.0});
            if (k % m2 > 0) t.push_back({k, k - 1, -1.0});
            if (k % m2 + 1 < m2) t.push_back({k, k + 1, -1.0});
            if (k >= m2) t.push_back({k, k - m2, -1.0});
            if (k + m2 < m2 * m2) t.push_back({k, k + m2, -1.0});
        }
        const SparseMatrix dominant = from(m2 * m2, t);
        const Solved s = solve_with(dominant, rhs(static_cast<std::size_t>(m2 * m2)), amg_strict);
        REQUIRE(s.factorization.fill == 0);
        REQUIRE(s.report.iterations.value_or(0) <= 10);
        detail::AmgBackend backend({.strength_threshold = 0.08, .smoother_degree = 2,
                                    .coarse_size = 1000, .max_levels = 10},
                                   nullptr);
        REQUIRE(backend.analyze(dominant).has_value());
        REQUIRE(backend.factorize(dominant.values()).has_value());
        REQUIRE(backend.level_sizes() == std::vector<std::size_t>{6400});
    }
    // A matrix smaller than the coarse size is solved by the coarse LU alone: one iteration.
    SolverConfig whole = amg_strict;
    whole.amg_coarse_size = 65 * 65;
    REQUIRE(solve_with(a, b, whole).report.iterations == 1);
}

TEST_CASE("amg: one coarsening per pattern, one hierarchy per matrix, one V-cycle per application") {
    const SparseMatrix a = poisson_2d(65, 0.2);
    SparseMatrix a2 = poisson_2d(65, 0.6);
    const std::vector<double> b = rhs(65 * 65);
    auto solver = LinearSolver::create(amg_strict).value();
    REQUIRE(solver.factorize(a).has_value());
    std::vector<double> x(b.size());
    const auto r1 = solver.solve(b, x);
    REQUIRE(r1.has_value());
    REQUIRE(solver.backend_counts().analyses == 1);
    REQUIRE(solver.backend_counts().factorizations == 1);
    // One V-cycle per GMRES iteration and one per cycle's update.
    REQUIRE(solver.backend_counts().solves >= *r1->iterations + 1);
    REQUIRE(solver.factorize(a2).has_value());
    const auto r2 = solver.solve(b, x);
    REQUIRE(r2.has_value());
    REQUIRE(backward_error(a2, b, x) <= amg_strict.max_backward_error);
    REQUIRE(solver.backend_counts().analyses == 1);
    REQUIRE(solver.backend_counts().factorizations == 2);
    // A new pattern analyzes again.
    const SparseMatrix other = poisson_2d(33);
    REQUIRE(solver.factorize(other).has_value());
    REQUIRE(solver.backend_counts().analyses == 2);
    // Repeated solves are bit-identical.
    std::vector<double> y1(33 * 33);
    std::vector<double> y2(33 * 33);
    const std::vector<double> b2 = rhs(33 * 33);
    REQUIRE(solver.solve(b2, y1).has_value());
    REQUIRE(solver.solve(b2, y2).has_value());
    REQUIRE(y1 == y2);
}

TEST_CASE("amg: the same coarsening is kept for a later matrix of the same pattern") {
    const SparseMatrix a = poisson_2d(65, 0.2);
    SparseMatrix a2 = a;
    // Values that would coarsen differently on their own: every x coupling weak.
    for (Index i = 0; i < a2.rows(); ++i) {
        for (Index k = a2.row_offsets()[i]; k < a2.row_offsets()[i + 1]; ++k) {
            const Index j = a2.col_indices()[k];
            if (j == i - 1 || j == i + 1) a2.values()[static_cast<std::size_t>(k)] *= 1e-4;
        }
    }
    detail::AmgBackend fresh({.strength_threshold = 0.08, .smoother_degree = 2, .coarse_size = 50,
                              .max_levels = 10},
                             nullptr);
    REQUIRE(fresh.analyze(a2).has_value());
    REQUIRE(fresh.factorize(a2.values()).has_value());
    detail::AmgBackend reused({.strength_threshold = 0.08, .smoother_degree = 2, .coarse_size = 50,
                               .max_levels = 10},
                              nullptr);
    REQUIRE(reused.analyze(a).has_value());
    REQUIRE(reused.factorize(a.values()).has_value());
    const auto first = reused.level_sizes();
    REQUIRE(reused.factorize(a2.values()).has_value());
    REQUIRE(reused.level_sizes() == first);
    REQUIRE(fresh.level_sizes() != first);
}

TEST_CASE("amg: a coarsest level that cannot be factorized falls back, or fails in strict mode") {
    // Every unknown on the coarsest level, which is exactly singular (two equal rows).
    const SparseMatrix a = from(3, {{0, 0, 1.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 1.0}, {2, 2, 1.0}});
    auto strict = LinearSolver::create(amg_strict).value();
    const auto f = strict.factorize(a);
    REQUIRE(!f);
    REQUIRE(f.error().code == ErrorCode::inaccurate_solve);
    SolverConfig fallback = amg_strict;
    fallback.fallback = LinearFallback::sparse_lu;
    auto s = LinearSolver::create(fallback).value();
    const auto g = s.factorize(a);
    REQUIRE(!g);
    REQUIRE(g.error().code == ErrorCode::singular_system);  // the fallback's own verdict
}

TEST_CASE("amg: every thread count gives the one-thread bits (Unit 58 SG2)") {
    // 200^2 unknowns: several row chunks (2048) and vector chunks (8192) on the finest levels.
    const Index m = 200;
    const SparseMatrix a = poisson_2d(m, 0.7);
    const std::vector<double> b = rhs(static_cast<std::size_t>(m * m));
    SolverConfig config = amg_strict;
    config.restart = 4;  // several cycles
    const Solved one = solve_with(a, b, config);
    for (const int threads : {2, 3, 4, 6, 8}) {
        CAPTURE(threads);
        config.threads = threads;
        const Solved s = solve_with(a, b, config);
        REQUIRE(s.report.iterations == one.report.iterations);
        REQUIRE(std::bit_cast<std::uint64_t>(s.report.backward_error) ==
                std::bit_cast<std::uint64_t>(one.report.backward_error));
        std::size_t differ = 0;
        for (std::size_t i = 0; i < b.size(); ++i) {
            if (std::bit_cast<std::uint64_t>(s.x[i]) != std::bit_cast<std::uint64_t>(one.x[i])) ++differ;
        }
        REQUIRE(differ == 0);
    }
    std::printf("amg threads: %zu unknowns, %zu iterations, identical at 2, 3, 4, 6, 8 threads\n",
                b.size(), one.report.iterations.value_or(0));
}
