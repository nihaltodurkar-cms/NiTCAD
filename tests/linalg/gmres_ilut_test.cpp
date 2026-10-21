// The gmres_ilut backend (ARCHITECTURE.md 6.10, Unit 44, LG2), through the public interface only:
// the configuration checks, agreement with Eigen SparseLU on the shared real and complex systems,
// the SparseLU fallback and the strict mode (LQ2), singular systems, ILUT build failures, reuse
// counts and reproducibility. Results are compared with Eigen's within tolerances, never bit for
// bit; every accepted solve is checked again with an independent backward error.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "fixtures.hpp"

using namespace NiTCAD::linalg;
using NiTCAD::base::ErrorCode;
using namespace fixtures;

namespace {

constexpr SolverConfig gmres{.backend = SolverBackend::gmres_ilut};
constexpr SolverConfig strict{.backend = SolverBackend::gmres_ilut,
                              .fallback = LinearFallback::none};

template <class Solver = LinearSolver>
Solver make(SolverConfig config = strict) {
    return Solver::create(config).value();
}

// Componentwise backward error of x (Oettli-Prager), written independently of the solver's.
template <MatrixScalar Scalar>
double backward_error(const BasicSparseMatrix<Scalar>& a, std::span<const Scalar> b,
                      std::span<const Scalar> x) {
    double worst = 0.0;
    for (Index i = 0; i < a.rows(); ++i) {
        Scalar ax{};
        double abs_ax = 0.0;
        for (Index k = a.row_offsets()[i]; k < a.row_offsets()[i + 1]; ++k) {
            const Scalar term = a.values()[k] * x[static_cast<std::size_t>(a.col_indices()[k])];
            ax += term;
            abs_ax += std::abs(term);
        }
        const auto u = static_cast<std::size_t>(i);
        const double denominator = abs_ax + std::abs(b[u]);
        const double row = denominator > 0.0 ? std::abs(b[u] - ax) / denominator : 0.0;
        worst = std::isnan(row) ? row : std::max(worst, row);
    }
    return worst;
}

template <MatrixScalar Scalar>
struct Solved {
    std::vector<Scalar> x;
    SolveReport report;
};

template <MatrixScalar Scalar>
Solved<Scalar> solve_with(const BasicSparseMatrix<Scalar>& a, std::span<const Scalar> b,
                          SolverConfig config) {
    auto solver = BasicLinearSolver<Scalar>::create(config).value();
    const auto f = solver.factorize(a);
    INFO("factorize: " << (f ? std::string("ok") : f.error().message));
    REQUIRE(f.has_value());
    Solved<Scalar> s{std::vector<Scalar>(b.size()), {}};
    const auto r = solver.solve(b, s.x);
    INFO("solve: " << (r ? std::string("ok") : r.error().message));
    REQUIRE(r.has_value());
    s.report = *r;
    return s;
}

SparseMatrix laplace_2d(Index m) {
    std::vector<Triplet> t;
    for (Index j = 0; j < m; ++j) {
        for (Index i = 0; i < m; ++i) {
            const Index k = j * m + i;
            t.push_back({k, k, 4.0});
            if (i > 0) t.push_back({k, k - 1, -1.0});
            if (i + 1 < m) t.push_back({k, k + 1, -1.0});
            if (j > 0) t.push_back({k, k - m, -1.0});
            if (j + 1 < m) t.push_back({k, k + m, -1.0});
        }
    }
    return from(m * m, t);
}

}  // namespace

// --- configuration -------------------------------------------------------------------------

TEST_CASE("gmres_ilut: default settings and the configuration checks") {
    const SolverConfig config;
    REQUIRE(config.restart == 30);
    REQUIRE(config.max_iterations == 1000);
    REQUIRE(config.ilut_drop_tolerance == 1e-4);
    REQUIRE(config.ilut_fill_factor == 10);
    REQUIRE(config.fallback == LinearFallback::sparse_lu);
    REQUIRE(make(gmres).config().backend == SolverBackend::gmres_ilut);

    for (const int threads : {0, 2, -1}) {
        const auto s = LinearSolver::create({.backend = SolverBackend::gmres_ilut, .threads = threads});
        REQUIRE(s.error().code == ErrorCode::invalid_input);
        REQUIRE(s.error().context->value == static_cast<double>(threads));
    }
    const auto invalid = [](SolverConfig c) {
        c.backend = SolverBackend::gmres_ilut;
        const auto s = LinearSolver::create(c);
        return !s.has_value() && s.error().code == ErrorCode::invalid_input;
    };
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    for (const int restart : {0, -1, 1001}) REQUIRE(invalid({.restart = restart}));
    for (const int limit : {0, -5}) REQUIRE(invalid({.max_iterations = limit}));
    for (const double drop : {-1e-6, 1.0, nan}) REQUIRE(invalid({.ilut_drop_tolerance = drop}));
    for (const int fill : {0, 1001}) REQUIRE(invalid({.ilut_fill_factor = fill}));
    REQUIRE(invalid({.fallback = static_cast<LinearFallback>(7)}));
    // The bounds themselves are valid.
    for (const SolverConfig c : {SolverConfig{.restart = 1}, SolverConfig{.restart = 1000},
                                 SolverConfig{.max_iterations = 1},
                                 SolverConfig{.ilut_drop_tolerance = 0.0},
                                 SolverConfig{.ilut_fill_factor = 1},
                                 SolverConfig{.ilut_fill_factor = 1000}}) {
        REQUIRE_FALSE(invalid(c));
    }
    // Checked for every backend, so a bad setting is caught before a backend is switched.
    REQUIRE(LinearSolver::create({.restart = 0}).error().code == ErrorCode::invalid_input);
}

TEST_CASE("gmres_ilut: the direct backends report no iterations and no fallback") {
    const System s = nonsymmetric(100);
    const Solved<double> e = solve_with<double>(s.a, s.b, {});
    REQUIRE_FALSE(e.report.iterations.has_value());
    REQUIRE_FALSE(e.report.fell_back);
}

// --- agreement with SparseLU ---------------------------------------------------------------

TEST_CASE("gmres_ilut: agrees with Eigen SparseLU on the shared systems") {
    // Strict mode: a fallback cannot stand in for GMRES here.
    const auto check = [](const char* name, const SparseMatrix& a, std::span<const double> b,
                          std::span<const double> exact, double tolerance) {
        const Solved<double> e = solve_with<double>(a, b, {});
        const Solved<double> g = solve_with<double>(a, b, strict);
        const double difference = max_abs_diff(g.x, e.x) / max_abs(e.x);
        std::printf("gmres_ilut %-22s: %3zu iterations, |x_g - x_e| / |x_e| %.2e, backward error "
                    "%.2e (Eigen %.2e), error vs exact %.2e\n",
                    name, *g.report.iterations, difference, g.report.backward_error,
                    e.report.backward_error, relative_error(g.x, exact));
        REQUIRE_FALSE(g.report.fell_back);
        REQUIRE(*g.report.iterations >= 1);
        REQUIRE(g.report.refinement_steps == 0);
        REQUIRE(g.report.backward_error <= 1e-8);
        REQUIRE(backward_error<double>(a, b, g.x) <= 1e-8);
        REQUIRE(difference <= tolerance);
        REQUIRE(relative_error(g.x, exact) <= tolerance);
    };
    {
        const System s = nonsymmetric(2000);
        check("nonsymmetric", s.a, s.b, s.x_exact, 1e-7);
    }
    {
        const SparseMatrix a = laplace_1d(random_conductances(999), true);
        std::vector<double> exact(1000), b(1000);
        for (std::size_t i = 0; i < exact.size(); ++i) exact[i] = 1.0 + hash01(70'000 + i);
        a.multiply(exact, b);
        check("contacted Laplacian", a, b, exact, 1e-6);
    }
    {
        const SparseMatrix a = laplace_2d(40);
        std::vector<double> exact(1600), b(1600);
        for (std::size_t i = 0; i < exact.size(); ++i) exact[i] = 1.0 + hash01(80'000 + i);
        a.multiply(exact, b);
        check("2D Laplacian", a, b, exact, 1e-6);
    }
    {
        // Rows over 40 decades and columns over 20: needs the equilibration GMRES runs inside.
        const System s = badly_scaled_coupled(20);
        check("badly scaled coupled", s.a, s.b, s.x_exact, 1e-6);
    }
}

TEST_CASE("gmres_ilut: complex systems") {
    {
        const Index n = 400;
        std::vector<ComplexTriplet> t;
        for (Index i = 0; i < n; ++i) {
            if (i > 0) t.push_back({i, i - 1, Complex{-1.3, 0.4}});
            t.push_back({i, i, Complex{2.5, 1.0 + 0.001 * i}});
            if (i + 1 < n) t.push_back({i, i + 1, Complex{-0.7, -0.2}});
        }
        const ComplexSparseMatrix a = complex_from(n, t);
        std::vector<Complex> exact, b(static_cast<std::size_t>(n));
        for (Index i = 0; i < n; ++i) exact.push_back({1.0 + std::sin(0.01 * i), std::cos(0.07 * i)});
        a.multiply(exact, b);
        const Solved<Complex> s = solve_with<Complex>(a, b, strict);
        REQUIRE_FALSE(s.report.fell_back);
        REQUIRE(s.report.backward_error <= 1e-8);
        REQUIRE(backward_error<Complex>(a, b, s.x) <= 1e-8);
        REQUIRE(max_abs_diff(s.x, exact) / max_abs(exact) < 1e-7);
    }
    // G + i omega C against Eigen over nine decades of frequency, reusing the analysis.
    {
        const SparseMatrix lap = laplace_1d(random_conductances(299), true);
        const auto n = static_cast<std::size_t>(lap.rows());
        std::vector<double> c(n);
        for (std::size_t i = 0; i < n; ++i) c[i] = 0.5 + hash01(5'000 + i);
        std::vector<Complex> b(n), xg(n), xe(n);
        for (std::size_t i = 0; i < n; ++i) b[i] = {hash01(i), hash01(10'000 + i) - 0.5};
        auto g = make<ComplexLinearSolver>();
        auto e = ComplexLinearSolver::create({}).value();
        double worst = 0.0;
        for (const double omega : {1e-6, 0.3, 30.0, 1e3}) {
            const ComplexSparseMatrix a = ac_matrix(lap, c, omega);
            REQUIRE(g.factorize(a).has_value());
            REQUIRE(e.factorize(a).has_value());
            const auto r = g.solve(b, xg);
            REQUIRE(r.has_value());
            REQUIRE_FALSE(r->fell_back);
            REQUIRE(backward_error<Complex>(a, b, xg) <= 1e-8);
            REQUIRE(e.solve(b, xe).has_value());
            worst = std::max(worst, max_abs_diff(xg, xe) / max_abs(xe));
        }
        std::printf("gmres_ilut complex G + i omega C: worst |x_g - x_e| / |x_e| %.2e\n", worst);
        REQUIRE(worst < 1e-4);
        REQUIRE(g.backend_counts().analyses == 1);
        REQUIRE(g.backend_counts().factorizations == 4);
    }
}

// --- failure: the fallback and the strict mode (LQ2) ---------------------------------------

TEST_CASE("gmres_ilut: an unconverged solve is inaccurate_solve in strict mode") {
    // A weak preconditioner (fill 1, large drop) and three iterations cannot solve the 2D
    // Laplacian.
    const SparseMatrix a = laplace_2d(30);
    const std::vector<double> b(900, 1.0);
    std::vector<double> x(900);
    LinearSolver solver = make({.backend = SolverBackend::gmres_ilut, .max_iterations = 3,
                                .ilut_drop_tolerance = 0.5, .ilut_fill_factor = 1,
                                .fallback = LinearFallback::none});
    REQUIRE(solver.factorize(a).has_value());
    const auto r = solver.solve(b, x);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::inaccurate_solve);
    REQUIRE(r.error().message.find("in 3 iterations") != std::string::npos);
    REQUIRE(*r.error().context->value > 1e-8);
    REQUIRE(solver.is_factorized());  // the preconditioner is still valid
}

TEST_CASE("gmres_ilut: an unconverged solve falls back to SparseLU and says so") {
    const SparseMatrix a = laplace_2d(30);
    const std::vector<double> b(900, 1.0);
    std::vector<double> x(900);
    LinearSolver solver = make({.backend = SolverBackend::gmres_ilut, .max_iterations = 3,
                                .ilut_drop_tolerance = 0.5, .ilut_fill_factor = 1});
    REQUIRE(solver.factorize(a).has_value());
    const auto first = solver.solve(b, x);
    REQUIRE(first.has_value());
    REQUIRE(first->fell_back);
    REQUIRE(*first->iterations == 3);  // the GMRES iterations spent before falling back
    REQUIRE(first->backward_error <= 1e-8);
    REQUIRE(backward_error<double>(a, b, x) <= 1e-8);
    const Solved<double> e = solve_with<double>(a, b, {});
    REQUIRE(max_abs_diff(x, e.x) / max_abs(e.x) <= 1e-13);  // the direct solve, refinement kept

    // The same factorization stays direct: no GMRES iterations are spent again.
    const std::size_t applications = solver.backend_counts().solves;
    const auto second = solver.solve(b, x);
    REQUIRE(second.has_value());
    REQUIRE(second->fell_back);
    REQUIRE(*second->iterations == 0);
    REQUIRE(solver.backend_counts().solves == applications);

    // A new factorization of the same pattern (no new analysis) tries GMRES again: with the
    // couplings 1e-4 of the diagonal it converges within the three iterations.
    SparseMatrix weak = a;
    for (Index r = 0; r < weak.rows(); ++r) {
        for (Index k = weak.row_offsets()[r]; k < weak.row_offsets()[r + 1]; ++k) {
            if (weak.col_indices()[k] != r) weak.values()[k] *= 1e-4;
        }
    }
    REQUIRE(solver.factorize(weak).has_value());
    REQUIRE(solver.analyses() == 1);
    const auto third = solver.solve(b, x);
    REQUIRE(third.has_value());
    REQUIRE_FALSE(third->fell_back);
    REQUIRE(*third->iterations <= 3);
}

TEST_CASE("gmres_ilut: a floating region is not detected (as with min_pivot_ratio = 0)") {
    // Known limitation, measured: with no pivots to check, a floating region with an inconsistent
    // right-hand side is "solved" by an x of norm about 1e17 whose componentwise backward error is
    // at rounding level, as Eigen SparseLU gives with min_pivot_ratio = 0 (|x| about 3e16). The
    // device layer rejects floating regions by topology before a solve (device.cpp), so this
    // backend leaves singularity to it.
    const SparseMatrix a = laplace_1d(random_conductances(99), false);
    const std::vector<double> b(100, 1.0);
    std::vector<double> x(100);
    LinearSolver solver = make();
    REQUIRE(solver.factorize(a).has_value());
    const auto r = solver.solve(b, x);
    REQUIRE(r.has_value());
    REQUIRE_FALSE(r->fell_back);
    REQUIRE(max_abs(x) > 1e10);

    LinearSolver direct = LinearSolver::create({.min_pivot_ratio = 0.0}).value();
    REQUIRE(direct.factorize(a).has_value());
    REQUIRE(direct.solve(b, x).has_value());
    REQUIRE(max_abs(x) > 1e10);

    // When GMRES stops short and falls back, the direct factorization's pivot check applies, and
    // the solver stays usable.
    LinearSolver one = make({.backend = SolverBackend::gmres_ilut, .max_iterations = 1,
                             .ilut_drop_tolerance = 0.5, .ilut_fill_factor = 1});
    REQUIRE(one.factorize(a).has_value());
    REQUIRE(one.solve(b, x).error().code == ErrorCode::singular_system);
    REQUIRE(one.factorize(nonsymmetric(100).a).has_value());
    REQUIRE(one.solve(std::vector<double>(100, 1.0), x).has_value());
}

TEST_CASE("gmres_ilut: an ILUT that cannot be built falls back, or fails in strict mode") {
    // ILUT does not pivot: the permutation has a zero diagonal, and with drop tolerance 0 Eigen's
    // zero-pivot shift (sqrt(drop) times the row norm) is 0, so the factor divides by zero.
    const SparseMatrix swap = from(2, {{0, 1, 1.0}, {1, 0, 1.0}});
    const std::array<double, 2> b{2.0, 3.0};
    std::array<double, 2> x{};

    LinearSolver none = make({.backend = SolverBackend::gmres_ilut, .ilut_drop_tolerance = 0.0,
                              .fallback = LinearFallback::none});
    const auto failed = none.factorize(swap);
    REQUIRE(failed.error().code == ErrorCode::inaccurate_solve);
    REQUIRE_FALSE(none.is_factorized());

    LinearSolver fallback = make({.backend = SolverBackend::gmres_ilut, .ilut_drop_tolerance = 0.0});
    const auto f = fallback.factorize(swap);
    REQUIRE(f.has_value());
    REQUIRE(f->pivot_ratio.has_value());  // the SparseLU factorization's
    const auto r = fallback.solve(b, x);
    REQUIRE(r.has_value());
    REQUIRE(r->fell_back);
    REQUIRE(*r->iterations == 0);
    REQUIRE(max_abs_diff(x, std::array<double, 2>{3.0, 2.0}) == 0.0);

    // With the default drop tolerance the shift makes a usable preconditioner, and GMRES solves it.
    LinearSolver shifted = make();
    REQUIRE(shifted.factorize(swap).has_value());
    const auto s = shifted.solve(b, x);
    REQUIRE(s.has_value());
    REQUIRE_FALSE(s->fell_back);
    REQUIRE(max_abs_diff(x, std::array<double, 2>{3.0, 2.0}) <= 1e-8);  // measured 4.3e-12
}

// --- reuse, reproducibility and edge cases -------------------------------------------------

TEST_CASE("gmres_ilut: one ordering per pattern, one ILUT per matrix") {
    LinearSolver solver = make();
    const SparseMatrix base = tridiagonal(500, -1.0, 4.0, -1.5);
    std::vector<double> b(500, 1.0), x(500);
    std::size_t iterations = 0;
    std::size_t solves = 0;
    for (int k = 0; k < 4; ++k) {
        SparseMatrix a = base;
        for (double& v : a.values()) v *= 1.0 + 0.1 * k;
        REQUIRE(solver.factorize(a).has_value());
        for (int r = 0; r < 3; ++r) {
            b[static_cast<std::size_t>(r)] += 1.0;
            const auto s = solver.solve(b, x);
            REQUIRE(s.has_value());
            iterations += *s->iterations;
            ++solves;
        }
    }
    REQUIRE(solver.analyses() == 1);
    REQUIRE(solver.backend_counts().analyses == 1);
    REQUIRE(solver.backend_counts().factorizations == 4);
    // One preconditioner application per iteration and one per cycle (here one cycle per solve).
    REQUIRE(solver.backend_counts().solves == iterations + solves);
    REQUIRE(solver.factorize(nonsymmetric(300).a).has_value());
    REQUIRE(solver.backend_counts().analyses == 2);
}

TEST_CASE("gmres_ilut: repeated solves are bit-identical") {
    const System s = badly_scaled_coupled(20);
    std::vector<double> first(s.b.size()), second(s.b.size());
    LinearSolver a = make();
    LinearSolver b = make();
    REQUIRE(a.factorize(s.a).has_value());
    REQUIRE(b.factorize(s.a).has_value());
    REQUIRE(a.solve(s.b, first).has_value());
    REQUIRE(b.solve(s.b, second).has_value());
    REQUIRE(a.solve(s.b, second).has_value());  // and the same solver again
    for (std::size_t i = 0; i < first.size(); ++i) {
        REQUIRE(std::bit_cast<std::uint64_t>(first[i]) == std::bit_cast<std::uint64_t>(second[i]));
    }
}

TEST_CASE("gmres_ilut: zero and huge right-hand sides, and no equilibration") {
    const System s = nonsymmetric(100);
    LinearSolver solver = make();
    REQUIRE(solver.factorize(s.a).has_value());
    std::vector<double> x(100, 1.0);
    const auto zero = solver.solve(std::vector<double>(100, 0.0), x);
    REQUIRE(zero.has_value());
    REQUIRE(zero->backward_error == 0.0);
    REQUIRE(*zero->iterations == 0);
    REQUIRE(std::ranges::all_of(x, [](double v) { return v == 0.0; }));

    std::vector<double> huge(s.b);
    for (double& v : huge) v *= 1e200;
    const auto big = solver.solve(huge, x);
    REQUIRE(big.has_value());
    for (std::size_t i = 0; i < x.size(); ++i) {
        REQUIRE(std::abs(x[i] / 1e200 - s.x_exact[i]) <= 1e-8);
    }

    LinearSolver raw = make({.backend = SolverBackend::gmres_ilut, .equilibrate = false,
                             .min_pivot_ratio = 0.0, .fallback = LinearFallback::none});
    REQUIRE(raw.factorize(s.a).has_value());
    const auto r = raw.solve(s.b, x);
    REQUIRE(r.has_value());
    REQUIRE(relative_error(x, s.x_exact) <= 1e-7);
}

TEST_CASE("gmres_ilut: a moved solver keeps its preconditioner") {
    const System s = nonsymmetric(200);
    LinearSolver source = make();
    REQUIRE(source.factorize(s.a).has_value());
    LinearSolver moved = std::move(source);
    std::vector<double> x(200);
    const auto r = moved.solve(s.b, x);
    REQUIRE(r.has_value());
    REQUIRE(relative_error(x, s.x_exact) <= 1e-7);
}
