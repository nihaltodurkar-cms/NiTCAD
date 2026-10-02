// Behavioural tests of the linear solver, written against the public interface only (Q5).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
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

}  // namespace

TEST_CASE("linear solver: default configuration") {
    const SolverConfig config;
    REQUIRE(config.backend == SolverBackend::eigen_sparse_lu);
    REQUIRE(config.threads == 1);
    REQUIRE(config.max_relative_residual == 1e-6);
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
    for (const double limit : {0.0, -1e-6, std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()}) {
        REQUIRE(LinearSolver::create({.max_relative_residual = limit}).error().code ==
                ErrorCode::invalid_input);
    }
    REQUIRE(LinearSolver::create({.backend = static_cast<SolverBackend>(7)}).error().code ==
            ErrorCode::invalid_input);
}

TEST_CASE("linear solver: analytic 1D Laplacian, tridiag(-1, 2, -1) x = 1") {
    // Exact solution x_i = i (n + 1 - i) / 2 for i = 1..n.
    constexpr Index n = 1000;
    const SparseMatrix a = tridiagonal(n, -1.0, 2.0, -1.0);
    const std::vector<double> b(n, 1.0);
    std::vector<double> x(n);

    LinearSolver solver = make_solver();
    REQUIRE(solver.factorize(a).has_value());
    const auto report = solver.solve(b, x);
    REQUIRE(report.has_value());
    // Backward-stable LU: residual of order eps ||A|| ||x|| / ||b|| = 2.2e-16 * 4 * 2.9e6 / 31.6,
    // about 8e-11 here.
    REQUIRE(report->relative_residual <= 1e-10);
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

TEST_CASE("linear solver: a structurally singular matrix is singular_system") {
    // Column 1 has no entries.
    const SparseMatrix a = from(3, {{0, 0, 1.0}, {1, 0, 1.0}, {2, 2, 1.0}});
    LinearSolver solver = make_solver();
    REQUIRE(solver.factorize(a).error().code == ErrorCode::singular_system);
    REQUIRE_FALSE(solver.is_factorized());
}

TEST_CASE("linear solver: the residual check rejects a solve above the threshold") {
    const System s = nonsymmetric(500);
    std::vector<double> x(500);

    LinearSolver reference = make_solver();
    REQUIRE(reference.factorize(s.a).has_value());
    const auto ok = reference.solve(s.b, x);
    REQUIRE(ok.has_value());
    const double residual = ok->relative_residual;
    REQUIRE(residual > 0.0);
    REQUIRE(residual <= 1e-14);

    // Same system, same build: the same residual, now above a stricter threshold.
    LinearSolver strict = make_solver({.max_relative_residual = residual / 2.0});
    REQUIRE(strict.factorize(s.a).has_value());
    const auto rejected = strict.solve(s.b, x);
    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error().code == ErrorCode::inaccurate_solve);
    REQUIRE(rejected.error().context->value == residual);
    REQUIRE(strict.is_factorized());  // the factorization itself is still valid
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

TEST_CASE("linear solver: a zero right-hand side gives x = 0 with zero residual") {
    const System s = nonsymmetric(50);
    const std::vector<double> b(50, 0.0);
    std::vector<double> x(50, 1.0);
    LinearSolver solver = make_solver();
    REQUIRE(solver.factorize(s.a).has_value());
    const auto report = solver.solve(b, x);
    REQUIRE(report.has_value());
    REQUIRE(report->relative_residual == 0.0);
    REQUIRE(std::ranges::all_of(x, [](double v) { return v == 0.0; }));
}

TEST_CASE("linear solver: an empty system is invalid_input") {
    LinearSolver solver = make_solver();
    const SparseMatrix empty;
    REQUIRE(solver.analyze(empty).error().code == ErrorCode::invalid_input);
    REQUIRE(solver.factorize(empty).error().code == ErrorCode::invalid_input);
    REQUIRE_FALSE(solver.is_analyzed());
}

TEST_CASE("linear solver: a moved solver keeps its factorization") {
    const System s = nonsymmetric(20);
    LinearSolver first = make_solver();
    REQUIRE(first.factorize(s.a).has_value());
    LinearSolver second = std::move(first);
    REQUIRE(second.is_factorized());
    std::vector<double> x(20);
    REQUIRE(second.solve(s.b, x).has_value());
    REQUIRE(max_abs_diff(x, s.x_exact) <= 1e-13);
}
