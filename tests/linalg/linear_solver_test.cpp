// Behavioural tests of the linear solver, written against the public interface only (Q5).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <thread>
#include <utility>
#include <vector>

#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

using namespace NiTCAD::linalg;
using NiTCAD::base::ErrorCode;

namespace {

SparseMatrix tridiagonal(Index n, double sub, double diag, double super) {
    std::vector<Triplet> t;
    for (Index i = 0; i < n; ++i) {
        if (i > 0) t.push_back({i, i - 1, sub});
        t.push_back({i, i, diag});
        if (i + 1 < n) t.push_back({i, i + 1, super});
    }
    return SparseMatrix::from_triplets(n, n, t).value();
}

SparseMatrix from(Index n, std::vector<Triplet> t) {
    return SparseMatrix::from_triplets(n, n, t).value();
}

LinearSolver make_solver(SolverConfig config = {}) { return LinearSolver::create(config).value(); }

double max_abs_diff(std::span<const double> a, std::span<const double> b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

double max_abs(std::span<const double> a) {
    double m = 0.0;
    for (const double v : a) m = std::max(m, std::abs(v));
    return m;
}

// Deterministic value in [0, 1) for integer k (no library RNG, so identical on every build).
double hash01(std::uint64_t k) {
    k = (k ^ (k >> 30)) * 0xbf58476d1ce4e5b9ULL;
    k = (k ^ (k >> 27)) * 0x94d049bb133111ebULL;
    k ^= k >> 31;
    return static_cast<double>(k >> 11) * 0x1.0p-53;
}

// A nonsymmetric, diagonally dominant system with a chosen solution; b = A x_exact.
struct System {
    SparseMatrix a;
    std::vector<double> x_exact;
    std::vector<double> b;
};

System nonsymmetric(Index n) {
    System s{tridiagonal(n, -1.3, 2.5, -0.7), {}, {}};
    for (Index i = 0; i < n; ++i) s.x_exact.push_back(1.0 + std::sin(0.01 * i));
    s.b.resize(static_cast<std::size_t>(n));
    s.a.multiply(s.x_exact, s.b);
    return s;
}

// 1D box-method Laplacian on n nodes with edge conductances g; Dirichlet rows at both ends add
// 1 to the diagonal. Without them the matrix is singular (pure Neumann, a floating region).
SparseMatrix laplace_1d(const std::vector<double>& g, bool dirichlet) {
    const auto n = static_cast<Index>(g.size()) + 1;
    std::vector<Triplet> t;
    for (Index e = 0; e + 1 < n; ++e) {
        const double w = g[static_cast<std::size_t>(e)];
        t.insert(t.end(), {{e, e, w}, {e + 1, e + 1, w}, {e, e + 1, -w}, {e + 1, e, -w}});
    }
    if (dirichlet) t.insert(t.end(), {{0, 0, 1.0}, {n - 1, n - 1, 1.0}});
    return from(n, std::move(t));
}

std::vector<double> random_conductances(std::size_t edges) {
    std::vector<double> g(edges);
    for (std::size_t e = 0; e < edges; ++e) g[e] = 1.0 / (0.5 + hash01(e));
    return g;
}

// Three unknowns per node on an m x m grid, neighbours coupled by dense 3x3 blocks (the shape of a
// psi/n/p Jacobian), then rows scaled by 10^[-20, 20] and columns by 10^[-10, 10], as carrier
// densities spanning many decades do. x_exact is chosen and b = A x_exact.
System badly_scaled_coupled(int m) {
    const int nodes = m * m;
    const Index n = 3 * nodes;
    std::vector<Triplet> t;
    std::uint64_t k = 0;
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            const Index p = j * m + i;
            const std::array<std::array<int, 2>, 4> nb{{{i - 1, j}, {i + 1, j}, {i, j - 1}, {i, j + 1}}};
            for (const auto& [a, c] : nb) {
                if (a < 0 || a >= m || c < 0 || c >= m) continue;
                const Index q = c * m + a;
                for (Index r = 0; r < 3; ++r) {
                    for (Index s = 0; s < 3; ++s) {
                        const double w = 0.1 + 0.9 * hash01(++k);
                        t.push_back({3 * p + r, 3 * q + s, -w});
                        t.push_back({3 * p + r, 3 * p + s, (r == s ? 1.6 : 0.3) * w});
                    }
                }
            }
            for (Index r = 0; r < 3; ++r) t.push_back({3 * p + r, 3 * p + r, 1.0});
        }
    }
    std::vector<double> row(static_cast<std::size_t>(n)), col(static_cast<std::size_t>(n));
    for (std::size_t i = 0; i < row.size(); ++i) {
        row[i] = std::pow(10.0, 40.0 * hash01(1'000'000 + i) - 20.0);
        col[i] = std::pow(10.0, 20.0 * hash01(2'000'000 + i) - 10.0);
    }
    for (Triplet& e : t) {
        e.value *= row[static_cast<std::size_t>(e.row)] * col[static_cast<std::size_t>(e.col)];
    }
    System s{from(n, std::move(t)), {}, std::vector<double>(static_cast<std::size_t>(n))};
    for (std::size_t i = 0; i < col.size(); ++i) {
        s.x_exact.push_back((1.0 + 0.5 * std::sin(0.1 * static_cast<double>(i))) / col[i]);
    }
    s.a.multiply(s.x_exact, s.b);
    return s;
}

double relative_error(std::span<const double> x, std::span<const double> exact) {
    return max_abs_diff(x, exact) / max_abs(exact);
}

}  // namespace

// --- configuration -------------------------------------------------------------------------

TEST_CASE("linear solver: default configuration") {
    const SolverConfig config;
    REQUIRE(config.backend == SolverBackend::eigen_sparse_lu);
    REQUIRE(config.threads == 1);
    REQUIRE(config.equilibrate);
    REQUIRE(config.max_backward_error == 1e-8);
    REQUIRE(config.max_refinement_steps == 1);
    REQUIRE(config.min_pivot_ratio == 1e-11);
    const LinearSolver solver = make_solver();
    REQUIRE_FALSE(solver.is_analyzed());
    REQUIRE_FALSE(solver.is_factorized());
    REQUIRE(solver.analyses() == 0);
    REQUIRE(solver.factorizations() == 0);
}

TEST_CASE("linear solver: invalid configuration is an error") {
    for (const int threads : {0, 2, -1}) {
        const auto s = LinearSolver::create({.threads = threads});
        REQUIRE(s.error().code == ErrorCode::invalid_input);
        REQUIRE(s.error().context->value == static_cast<double>(threads));
    }
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    for (const double limit : {0.0, -1e-6, nan, inf}) {
        REQUIRE(LinearSolver::create({.max_backward_error = limit}).error().code ==
                ErrorCode::invalid_input);
    }
    REQUIRE(LinearSolver::create({.max_refinement_steps = -1}).error().code ==
            ErrorCode::invalid_input);
    for (const double ratio : {-1e-12, 1.0, nan}) {
        REQUIRE(LinearSolver::create({.min_pivot_ratio = ratio}).error().code ==
                ErrorCode::invalid_input);
    }
    REQUIRE(LinearSolver::create({.min_pivot_ratio = 0.0}).has_value());  // check disabled
    REQUIRE(LinearSolver::create({.equilibrate = false}).error().code == ErrorCode::invalid_input);
    REQUIRE(LinearSolver::create({.equilibrate = false, .min_pivot_ratio = 0.0}).has_value());
    REQUIRE(LinearSolver::create({.backend = static_cast<SolverBackend>(7)}).error().code ==
            ErrorCode::invalid_input);
}

// --- known solutions -----------------------------------------------------------------------

TEST_CASE("linear solver: analytic 1D Laplacian, tridiag(-1, 2, -1) x = 1") {
    // Exact solution x_i = i (n + 1 - i) / 2 for i = 1..n.
    constexpr Index n = 1000;
    const SparseMatrix a = tridiagonal(n, -1.0, 2.0, -1.0);
    const std::vector<double> b(n, 1.0);
    std::vector<double> x(n);

    LinearSolver solver = make_solver();
    const auto factored = solver.factorize(a);
    REQUIRE(factored.has_value());
    REQUIRE(factored->pivot_ratio.has_value());
    REQUIRE(*factored->pivot_ratio > 0.1);  // tridiag(-1, 2, -1): pivots (k + 1) / k
    const auto report = solver.solve(b, x);
    REQUIRE(report.has_value());
    REQUIRE(report->backward_error <= 1e-15);  // a few units of rounding (1.1e-16)
    REQUIRE(report->refinement_steps == 0);
    REQUIRE_FALSE(report->iterations.has_value());  // direct backend

    double worst = 0.0;
    for (Index i = 1; i <= n; ++i) {
        const double exact = 0.5 * i * (n + 1 - i);
        worst = std::max(worst, std::abs(x[static_cast<std::size_t>(i - 1)] - exact) / exact);
    }
    REQUIRE(worst <= 1e-9);
}

TEST_CASE("linear solver: analytic 2D Laplacian eigenvector (5-point stencil)") {
    // On an m x m interior grid with h = 1 / (m + 1), u_ij = sin(pi i h) sin(pi j h) satisfies
    // A u = lambda u with lambda = 4 (1 - cos(pi h)), so A x = lambda u has the solution x = u.
    // The pattern is not banded, so this exercises the ordering and fill-in.
    constexpr Index m = 40;
    constexpr Index n = m * m;
    const double h = 1.0 / (m + 1);
    const double lambda = 4.0 * (1.0 - std::cos(std::numbers::pi * h));
    std::vector<Triplet> t;
    std::vector<double> u(n);
    for (Index j = 0; j < m; ++j) {
        for (Index i = 0; i < m; ++i) {
            const Index k = j * m + i;
            t.push_back({k, k, 4.0});
            if (i > 0) t.push_back({k, k - 1, -1.0});
            if (i + 1 < m) t.push_back({k, k + 1, -1.0});
            if (j > 0) t.push_back({k, k - m, -1.0});
            if (j + 1 < m) t.push_back({k, k + m, -1.0});
            u[static_cast<std::size_t>(k)] =
                std::sin(std::numbers::pi * (i + 1) * h) * std::sin(std::numbers::pi * (j + 1) * h);
        }
    }
    std::vector<double> b(n);
    std::ranges::transform(u, b.begin(), [&](double v) { return lambda * v; });
    std::vector<double> x(n);

    LinearSolver solver = make_solver();
    REQUIRE(solver.factorize(from(n, t)).has_value());
    REQUIRE(solver.solve(b, x).has_value());
    REQUIRE(max_abs_diff(x, u) <= 1e-10);
}

TEST_CASE("linear solver: nonsymmetric systems are solved, not their transposes") {
    // [2 1 0] [1]   [ 4]
    // [0 3 1] [2] = [ 9]      A^T x = b has a different solution.
    // [1 0 4] [3]   [13]
    const SparseMatrix a = from(3, {{0, 0, 2.0}, {0, 1, 1.0}, {1, 1, 3.0}, {1, 2, 1.0},
                                    {2, 0, 1.0}, {2, 2, 4.0}});
    const std::array<double, 3> b{4.0, 9.0, 13.0};
    std::array<double, 3> x{};
    LinearSolver solver = make_solver();
    REQUIRE(solver.factorize(a).has_value());
    REQUIRE(solver.solve(b, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 3>{1.0, 2.0, 3.0}) <= 1e-14);

    const System s = nonsymmetric(500);
    std::vector<double> y(500);
    REQUIRE(solver.factorize(s.a).has_value());
    REQUIRE(solver.solve(s.b, y).has_value());
    REQUIRE(max_abs_diff(y, s.x_exact) <= 1e-12);
}

TEST_CASE("linear solver: zero diagonal entries are handled by pivoting") {
    // A permutation and a saddle-point (Dirichlet-multiplier-like) system.
    const SparseMatrix swap = from(2, {{0, 1, 1.0}, {1, 0, 1.0}});
    std::array<double, 2> x{};
    LinearSolver solver = make_solver();
    REQUIRE(solver.factorize(swap).has_value());
    REQUIRE(solver.solve(std::array<double, 2>{2.0, 3.0}, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 2>{3.0, 2.0}) == 0.0);

    // [2 1] [1]   [3]
    // [1 0] [1] = [1]
    const SparseMatrix saddle = from(2, {{0, 0, 2.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 0.0}});
    REQUIRE(solver.factorize(saddle).has_value());
    REQUIRE(solver.solve(std::array<double, 2>{3.0, 1.0}, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 2>{1.0, 1.0}) <= 1e-15);
}

// --- acceptance ----------------------------------------------------------------------------

TEST_CASE("linear solver: a backward-stable but ill-conditioned solve is accepted") {
    // Regression for the review: the uniform 1D Laplacian with n = 600,000 (condition about
    // 1.5e11) has a backward-stable solution whose ||r||/||b|| exceeds 1e-6; the former gate
    // rejected it.
    constexpr Index n = 600'000;
    const SparseMatrix a = tridiagonal(n, -1.0, 2.0, -1.0);
    const std::vector<double> b(n, 1.0);
    std::vector<double> x(n);
    LinearSolver solver = make_solver({.max_refinement_steps = 0});
    REQUIRE(solver.factorize(a).has_value());
    const auto report = solver.solve(b, x);
    REQUIRE(report.has_value());
    REQUIRE(report->backward_error <= 1e-15);
    REQUIRE(report->relative_residual > 1e-6);
}

TEST_CASE("linear solver: badly scaled systems need equilibration, and the gate sees it") {
    // Regression for the review: rows scaled by up to 1e+-20 hide inaccurate small rows from
    // ||r||/||b||. Without equilibration the componentwise backward error exposes the bad solve;
    // with it the forward error is at rounding level.
    const System s = badly_scaled_coupled(30);
    std::vector<double> x(s.b.size());

    LinearSolver raw = make_solver({.equilibrate = false, .max_refinement_steps = 0, .min_pivot_ratio = 0.0});
    REQUIRE(raw.factorize(s.a).has_value());
    const auto rejected = raw.solve(s.b, x);
    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error().code == ErrorCode::inaccurate_solve);
    REQUIRE(*rejected.error().context->value > 1e-6);

    LinearSolver scaled = make_solver();
    REQUIRE(scaled.factorize(s.a).has_value());
    const auto ok = scaled.solve(s.b, x);
    REQUIRE(ok.has_value());
    // Measured 1.9e-11: equilibration bounds each row's largest entry, but entries within a row
    // still span decades, so this is above the few-epsilon of a well-scaled system.
    REQUIRE(ok->backward_error <= 1e-10);
    REQUIRE(relative_error(x, s.x_exact) <= 1e-10);
}

TEST_CASE("linear solver: iterative refinement recovers a solve above the limit") {
    const System s = badly_scaled_coupled(30);
    std::vector<double> x(s.b.size());
    LinearSolver solver = make_solver({.equilibrate = false, .max_refinement_steps = 3, .min_pivot_ratio = 0.0});
    REQUIRE(solver.factorize(s.a).has_value());
    const auto report = solver.solve(s.b, x);
    REQUIRE(report.has_value());
    REQUIRE(report->refinement_steps >= 1);
    REQUIRE(report->backward_error <= 1e-8);
}

TEST_CASE("linear solver: the backward-error gate rejects a solve above the limit") {
    const System s = nonsymmetric(500);
    std::vector<double> x(500);

    LinearSolver reference = make_solver();
    REQUIRE(reference.factorize(s.a).has_value());
    const auto ok = reference.solve(s.b, x);
    REQUIRE(ok.has_value());
    const double omega = ok->backward_error;
    REQUIRE(omega > 0.0);
    REQUIRE(omega <= 1e-15);

    // Same system, same build: the same backward error, now above a stricter limit.
    LinearSolver strict =
        make_solver({.max_backward_error = omega / 2.0, .max_refinement_steps = 0});
    REQUIRE(strict.factorize(s.a).has_value());
    const auto rejected = strict.solve(s.b, x);
    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error().code == ErrorCode::inaccurate_solve);
    REQUIRE(rejected.error().context->value == omega);
    REQUIRE(strict.is_factorized());  // the factorization itself is still valid
}

TEST_CASE("linear solver: huge right-hand sides do not overflow the acceptance check") {
    // The former check summed squares and overflowed above about 1e154.
    const System s = nonsymmetric(100);
    std::vector<double> b(s.b);
    for (double& v : b) v *= 1e200;
    std::vector<double> x(100);
    LinearSolver solver = make_solver();
    REQUIRE(solver.factorize(s.a).has_value());
    const auto report = solver.solve(b, x);
    REQUIRE(report.has_value());
    for (std::size_t i = 0; i < x.size(); ++i) {
        REQUIRE(std::abs(x[i] / 1e200 - s.x_exact[i]) <= 1e-13);
    }
}

TEST_CASE("linear solver: a zero right-hand side gives x = 0 with zero residual") {
    const System s = nonsymmetric(50);
    const std::vector<double> b(50, 0.0);
    std::vector<double> x(50, 1.0);
    LinearSolver solver = make_solver();
    REQUIRE(solver.factorize(s.a).has_value());
    const auto report = solver.solve(b, x);
    REQUIRE(report.has_value());
    REQUIRE(report->relative_residual == 0.0);
    REQUIRE(report->backward_error == 0.0);
    REQUIRE(std::ranges::all_of(x, [](double v) { return v == 0.0; }));
}

TEST_CASE("linear solver: repeated solves are bit-identical") {
    const System s = badly_scaled_coupled(20);
    std::vector<double> first(s.b.size()), second(s.b.size());
    LinearSolver a = make_solver();
    LinearSolver b = make_solver();
    REQUIRE(a.factorize(s.a).has_value());
    REQUIRE(b.factorize(s.a).has_value());
    REQUIRE(a.solve(s.b, first).has_value());
    REQUIRE(b.solve(s.b, second).has_value());
    for (std::size_t i = 0; i < first.size(); ++i) {
        REQUIRE(std::bit_cast<std::uint64_t>(first[i]) == std::bit_cast<std::uint64_t>(second[i]));
    }
}

// --- singular systems ----------------------------------------------------------------------

TEST_CASE("linear solver: a numerically singular matrix is singular_system and keeps nothing") {
    const SparseMatrix singular = from(2, {{0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 2.0}, {1, 1, 4.0}});
    LinearSolver solver = make_solver();
    const auto r = solver.factorize(singular);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::singular_system);
    REQUIRE_FALSE(solver.is_factorized());
    REQUIRE_FALSE(solver.is_analyzed());
    REQUIRE(solver.factorizations() == 0);

    // The same solver recovers on the next matrix, with a fresh analysis.
    const SparseMatrix regular = from(2, {{0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 2.0}, {1, 1, 5.0}});
    REQUIRE(solver.factorize(regular).has_value());
    REQUIRE(solver.analyses() == 2);
    std::array<double, 2> x{};
    REQUIRE(solver.solve(std::array<double, 2>{3.0, 7.0}, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 2>{1.0, 1.0}) <= 1e-15);
}

TEST_CASE("linear solver: an exactly zero pivot names its original column") {
    // Pure Neumann Laplacian with unit conductances: the last pivot is exactly zero.
    const SparseMatrix a = laplace_1d(std::vector<double>(99, 1.0), false);
    LinearSolver solver = make_solver();
    const auto r = solver.factorize(a);
    REQUIRE(r.error().code == ErrorCode::singular_system);
    REQUIRE(r.error().message == "zero pivot in LU factorization");
    REQUIRE(r.error().context->index.has_value());
    REQUIRE(*r.error().context->index < 100);
}

TEST_CASE("linear solver: a floating region is singular_system, not a solution") {
    // Regression for the review: with non-uniform conductances the singular pivot is at rounding
    // level, not zero, and the former code accepted a solution with an arbitrary offset whenever
    // the right-hand side was consistent (charge neutral).
    const SparseMatrix a = laplace_1d(random_conductances(999), false);
    LinearSolver solver = make_solver();
    const auto r = solver.factorize(a);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::singular_system);
    REQUIRE(*r.error().context->value < 1e-13);
    REQUIRE_FALSE(solver.is_factorized());

    // The same mesh with its Dirichlet contacts is accepted.
    const auto contacted = solver.factorize(laplace_1d(random_conductances(999), true));
    REQUIRE(contacted.has_value());
    REQUIRE(*contacted->pivot_ratio > 1e-8);
}

TEST_CASE("linear solver: one floating region among contacted ones is found") {
    // Nodes 0..499 are a chain with a contact at node 0; nodes 500..999 form a separate chain
    // with no contact. The reported column lies in the floating chain.
    std::vector<Triplet> t;
    const std::vector<double> g = random_conductances(1000);
    for (Index i = 0; i + 1 < 1000; ++i) {
        if (i == 499) continue;
        const double w = g[static_cast<std::size_t>(i)];
        t.insert(t.end(), {{i, i, w}, {i + 1, i + 1, w}, {i, i + 1, -w}, {i + 1, i, -w}});
    }
    t.push_back({0, 0, 1.0});
    LinearSolver solver = make_solver();
    const auto r = solver.factorize(from(1000, t));
    REQUIRE(r.error().code == ErrorCode::singular_system);
    REQUIRE(*r.error().context->index >= 500);
}

TEST_CASE("linear solver: min_pivot_ratio = 0 disables the floating-region check") {
    const SparseMatrix a = laplace_1d(random_conductances(999), false);
    LinearSolver solver = make_solver({.min_pivot_ratio = 0.0});
    const auto r = solver.factorize(a);
    REQUIRE(r.has_value());
    REQUIRE(*r->pivot_ratio < 1e-13);
}

TEST_CASE("linear solver: a structurally singular matrix is singular_system with its index") {
    LinearSolver solver = make_solver();
    // Column 1 has no entries.
    const auto col = solver.factorize(from(3, {{0, 0, 1.0}, {1, 0, 1.0}, {2, 2, 1.0}}));
    REQUIRE(col.error().code == ErrorCode::singular_system);
    REQUIRE(col.error().context->index == 1);
    REQUIRE_FALSE(solver.is_analyzed());
    REQUIRE(solver.analyses() == 0);
    // Row 2 has no entries.
    const auto row = solver.analyze(from(3, {{0, 0, 1.0}, {1, 1, 1.0}, {0, 2, 1.0}}));
    REQUIRE(row.error().code == ErrorCode::singular_system);
    REQUIRE(row.error().context->index == 2);
}

TEST_CASE("linear solver: a row or column of explicit zeros is singular_system") {
    LinearSolver solver = make_solver();
    const auto row = solver.factorize(from(3, {{0, 0, 1.0}, {1, 1, 0.0}, {1, 2, 0.0}, {2, 2, 1.0}}));
    REQUIRE(row.error().code == ErrorCode::singular_system);
    REQUIRE(row.error().context->index == 1);
    const auto col = solver.factorize(from(3, {{0, 0, 1.0}, {1, 1, 1.0}, {0, 2, 0.0}, {2, 0, 1.0}, {1, 2, 0.0}}));
    REQUIRE(col.error().code == ErrorCode::singular_system);
    REQUIRE(col.error().context->index == 2);
}

// --- reuse and lifecycle -------------------------------------------------------------------

TEST_CASE("linear solver: symbolic analysis is reused while the pattern is unchanged") {
    System s = nonsymmetric(200);
    LinearSolver solver = make_solver();
    REQUIRE(solver.analyze(s.a).has_value());
    REQUIRE(solver.is_analyzed());
    REQUIRE_FALSE(solver.is_factorized());

    const std::vector<double> base_values(s.a.values().begin(), s.a.values().end());
    std::vector<double> x(200);
    for (int k = 0; k < 4; ++k) {
        // New values on the same pattern: (1 + 0.1 k) A, so x = x_exact / (1 + 0.1 k).
        const double scale = 1.0 + 0.1 * k;
        std::ranges::transform(base_values, s.a.values().begin(),
                               [&](double v) { return scale * v; });
        REQUIRE(solver.factorize(s.a).has_value());
        REQUIRE(solver.solve(s.b, x).has_value());
        std::vector<double> expected(s.x_exact);
        for (double& v : expected) v /= scale;
        REQUIRE(max_abs_diff(x, expected) <= 1e-12);
    }
    REQUIRE(solver.analyses() == 1);
    REQUIRE(solver.factorizations() == 4);
}

TEST_CASE("linear solver: explicit zeros keep the pattern when a value becomes nonzero") {
    // An entry that is zero in one Newton iteration and nonzero in the next must not trigger a new
    // analysis, and must be used by the factorization.
    auto a = SparseMatrix::from_triplets(
                 2, 2, std::array<Triplet, 4>{{{0, 0, 2.0}, {0, 1, 0.0}, {1, 0, 0.0}, {1, 1, 4.0}}})
                 .value();
    LinearSolver solver = make_solver();
    std::array<double, 2> x{};
    REQUIRE(solver.factorize(a).has_value());
    REQUIRE(solver.solve(std::array<double, 2>{2.0, 4.0}, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 2>{1.0, 1.0}) == 0.0);

    a.values()[1] = 1.0;  // (0, 1)
    a.values()[2] = 1.0;  // (1, 0)
    REQUIRE(solver.factorize(a).has_value());
    REQUIRE(solver.solve(std::array<double, 2>{3.0, 5.0}, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 2>{1.0, 1.0}) <= 1e-15);
    REQUIRE(solver.analyses() == 1);
}

TEST_CASE("linear solver: factorize analyzes when needed and re-analyzes a changed pattern") {
    const SparseMatrix a = from(3, {{0, 0, 2.0}, {1, 1, 3.0}, {2, 2, 4.0}});
    const SparseMatrix b = from(3, {{0, 0, 2.0}, {0, 2, 1.0}, {1, 1, 3.0}, {2, 2, 4.0}});
    const std::array<double, 3> rhs{2.0, 3.0, 4.0};
    std::array<double, 3> x{};
    LinearSolver solver = make_solver();

    REQUIRE(solver.factorize(a).has_value());  // no explicit analyze
    REQUIRE(solver.analyses() == 1);

    REQUIRE(solver.factorize(b).has_value());  // superset pattern
    REQUIRE(solver.analyses() == 2);
    REQUIRE(solver.solve(rhs, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 3>{0.5, 1.0, 1.0}) <= 1e-15);

    REQUIRE(solver.factorize(a).has_value());  // any change re-analyzes, subsets included
    REQUIRE(solver.analyses() == 3);
    REQUIRE(solver.solve(rhs, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 3>{1.0, 1.0, 1.0}) <= 1e-15);
    REQUIRE(solver.factorizations() == 3);
}

TEST_CASE("linear solver: non-finite matrix values or right-hand sides are invalid_input") {
    System s = nonsymmetric(10);
    std::vector<double> x(10);
    LinearSolver solver = make_solver();

    // Entry 0 is row 0; each later row starts with its sub-diagonal entry, so entry 3k - 1 is
    // in row k. Entry 14 is in row 5.
    s.a.values()[14] = std::numeric_limits<double>::quiet_NaN();
    const auto bad_a = solver.factorize(s.a);
    REQUIRE(bad_a.error().code == ErrorCode::invalid_input);
    REQUIRE(bad_a.error().context->index == 5);
    REQUIRE(solver.is_analyzed());  // the pattern analysis stays valid
    REQUIRE_FALSE(solver.is_factorized());

    s.a.values()[14] = -1.3;
    REQUIRE(solver.factorize(s.a).has_value());
    REQUIRE(solver.analyses() == 1);

    s.b[7] = std::numeric_limits<double>::infinity();
    const auto bad_b = solver.solve(s.b, x);
    REQUIRE(bad_b.error().code == ErrorCode::invalid_input);
    REQUIRE(bad_b.error().context->index == 7);
}

TEST_CASE("linear solver: an empty system is invalid_input") {
    LinearSolver solver = make_solver();
    const SparseMatrix empty;
    REQUIRE(solver.analyze(empty).error().code == ErrorCode::invalid_input);
    REQUIRE(solver.factorize(empty).error().code == ErrorCode::invalid_input);
    REQUIRE_FALSE(solver.is_analyzed());
}

TEST_CASE("linear solver: a moved solver keeps its factorization; the source keeps nothing") {
    const System s = nonsymmetric(20);
    LinearSolver first = make_solver();
    REQUIRE(first.factorize(s.a).has_value());
    LinearSolver second = std::move(first);
    REQUIRE(second.is_factorized());
    // solve() requires is_factorized(), so the moved-from solver fails that precondition
    // instead of dereferencing a moved-out backend.
    REQUIRE_FALSE(first.is_factorized());  // NOLINT(bugprone-use-after-move)
    REQUIRE_FALSE(first.is_analyzed());    // NOLINT(bugprone-use-after-move)
    std::vector<double> x(20);
    REQUIRE(second.solve(s.b, x).has_value());
    REQUIRE(max_abs_diff(x, s.x_exact) <= 1e-13);

    LinearSolver third = make_solver();
    third = std::move(second);
    REQUIRE(third.is_factorized());
    REQUIRE_FALSE(second.is_factorized());  // NOLINT(bugprone-use-after-move)
    REQUIRE(third.solve(s.b, x).has_value());
}

// --- limits of the singularity checks (follow-up review of the hardening) -------------------

namespace {

// n-node chain with unit conductances, a zeroth-order term delta on every node, and (if w > 0)
// edge (n/2 - 1, n/2) weakened to w and a Dirichlet row at node 0.
SparseMatrix anchored_chain(Index n, double delta, double w) {
    std::vector<Triplet> t;
    for (Index i = 0; i + 1 < n; ++i) {
        const double g = (w > 0.0 && i == n / 2 - 1) ? w : 1.0;
        t.insert(t.end(), {{i, i, g}, {i + 1, i + 1, g}, {i, i + 1, -g}, {i + 1, i, -g}});
    }
    for (Index i = 0; i < n; ++i) t.push_back({i, i, delta});
    if (w > 0.0) t.push_back({0, 0, 1.0});
    return from(n, std::move(t));
}

}  // namespace

TEST_CASE("linear solver: a weakly anchored region near min_pivot_ratio (known false positive)") {
    // Both systems are nonsingular and solve with a backward error of about 1e-16, but the pivot
    // ratio cannot tell them from a floating region: about n delta / 2 for a zeroth-order anchor
    // and w / 2 for a weak link. These cases pin where the default threshold cuts.
    LinearSolver solver = make_solver();
    const auto zeroth_ok = solver.factorize(anchored_chain(1000, 1e-12, 0.0));
    REQUIRE(zeroth_ok.has_value());
    REQUIRE(*zeroth_ok->pivot_ratio > 4e-10);
    REQUIRE(*zeroth_ok->pivot_ratio < 6e-10);
    REQUIRE(solver.factorize(anchored_chain(1000, 1e-14, 0.0)).error().code ==
            ErrorCode::singular_system);

    const auto link_ok = solver.factorize(anchored_chain(1000, 0.0, 1e-9));
    REQUIRE(link_ok.has_value());
    REQUIRE(*link_ok->pivot_ratio > 4e-10);
    REQUIRE(*link_ok->pivot_ratio < 6e-10);
    REQUIRE(solver.factorize(anchored_chain(1000, 0.0, 1e-12)).error().code ==
            ErrorCode::singular_system);

    // With the check disabled the weak link is solved accurately.
    LinearSolver unchecked = make_solver({.min_pivot_ratio = 0.0});
    REQUIRE(unchecked.factorize(anchored_chain(1000, 0.0, 1e-12)).has_value());
    std::vector<double> b(1000, 1.0), x(1000);
    const auto report = unchecked.solve(b, x);
    REQUIRE(report.has_value());
    REQUIRE(report->backward_error <= 1e-15);
}

TEST_CASE("linear solver: the singular index points into the offending block") {
    // A contacted 996-node chain plus a 4-node block at columns [lo, lo + 4) that is exactly
    // singular (unit conductances: an exactly zero pivot) or floating (random conductances: a
    // rounding-level pivot). The reported original column must lie in the block wherever it is.
    for (const bool exact : {true, false}) {
        for (const Index lo : {3, 417, 995}) {
            std::vector<Triplet> t;
            Index previous = -1;
            for (Index i = 0; i < 1000; ++i) {
                if (i >= lo && i < lo + 4) continue;
                t.push_back({i, i, previous < 0 ? 3.0 : 2.0});
                if (previous >= 0) t.insert(t.end(), {{i, previous, -1.0}, {previous, i, -1.0}});
                previous = i;
            }
            const std::array<double, 3> g = exact ? std::array<double, 3>{1.0, 1.0, 1.0}
                                                  : std::array<double, 3>{0.731, 1.377, 0.519};
            for (Index e = 0; e < 3; ++e) {
                const Index a = lo + e;
                const double w = g[static_cast<std::size_t>(e)];
                t.insert(t.end(), {{a, a, w}, {a + 1, a + 1, w}, {a, a + 1, -w}, {a + 1, a, -w}});
            }
            LinearSolver solver = make_solver();
            const auto r = solver.factorize(from(1000, t));
            REQUIRE(r.error().code == ErrorCode::singular_system);
            const std::size_t index = *r.error().context->index;
            REQUIRE(index >= static_cast<std::size_t>(lo));
            REQUIRE(index < static_cast<std::size_t>(lo + 4));
        }
    }
}

TEST_CASE("linear solver: a pivot that overflows is inaccurate_solve at factorization") {
    // Without equilibration, eliminating column 0 gives u_11 = -1e308 - 0.9 * 1.7e308 = -inf.
    const SparseMatrix a = from(2, {{0, 0, 1.0}, {0, 1, 1.7e308}, {1, 0, 0.9}, {1, 1, -1e308}});
    LinearSolver solver = make_solver({.equilibrate = false, .min_pivot_ratio = 0.0});
    const auto r = solver.factorize(a);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::inaccurate_solve);
    REQUIRE_FALSE(solver.is_factorized());

    // Equilibrated, the same matrix is harmless.
    REQUIRE(make_solver().factorize(a).has_value());
}

TEST_CASE("linear solver: independent solvers on separate threads match a serial run") {
    const System s = badly_scaled_coupled(20);
    const auto run = [&s](std::vector<double>& x) {
        LinearSolver solver = make_solver();
        x.assign(s.b.size(), 0.0);
        return solver.factorize(s.a).has_value() && solver.solve(s.b, x).has_value();
    };
    std::vector<double> serial;
    REQUIRE(run(serial));
    std::vector<double> first, second;
    bool ok_first = false;
    bool ok_second = false;
    {
        std::jthread t1([&] { ok_first = run(first); });
        std::jthread t2([&] { ok_second = run(second); });
    }
    REQUIRE(ok_first);
    REQUIRE(ok_second);
    for (std::size_t i = 0; i < serial.size(); ++i) {
        REQUIRE(std::bit_cast<std::uint64_t>(first[i]) == std::bit_cast<std::uint64_t>(serial[i]));
        REQUIRE(std::bit_cast<std::uint64_t>(second[i]) == std::bit_cast<std::uint64_t>(serial[i]));
    }
}
