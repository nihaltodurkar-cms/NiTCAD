// The MKL PARDISO backend (ARCHITECTURE.md 6.10, Unit 18), through the public interface only: the
// runtime loader, agreement with Eigen SparseLU on the shared systems, the singularity checks,
// complex matrices, the hard reuse gate (one analysis per pattern, one numeric factorization per
// matrix, one triangular solve per right-hand side) and thread counts. Results are compared with
// Eigen's within tolerances, never bit for bit.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <thread>
#include <vector>

#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "fixtures.hpp"
#include "mkl_test_runtime.hpp"

using namespace NiTCAD::linalg;
using NiTCAD::base::ErrorCode;
using namespace fixtures;

namespace {

constexpr SolverConfig pardiso{.backend = SolverBackend::mkl_pardiso};

template <class Solver = LinearSolver>
Solver make(SolverConfig config = pardiso) {
    return Solver::create(config).value();
}

// x from a fresh solver of the given configuration, with its reports.
struct Solved {
    std::vector<double> x;
    FactorizationReport factorization;
    SolveReport solve;
};

Solved solve_with(const SparseMatrix& a, std::span<const double> b, SolverConfig config) {
    LinearSolver solver = make(config);
    const auto f = solver.factorize(a);
    INFO("factorize: " << (f ? std::string("ok") : f.error().message) << " "
                       << (f || !f.error().context ? 0.0 : f.error().context->value.value_or(0.0)));
    REQUIRE(f.has_value());
    Solved s{std::vector<double>(b.size()), *f, {}};
    const auto r = solver.solve(b, s.x);
    INFO("solve: " << (r ? std::string("ok") : r.error().message) << " "
                   << (r || !r.error().context ? 0.0 : r.error().context->value.value_or(0.0)));
    REQUIRE(r.has_value());
    s.solve = *r;
    return s;
}

}  // namespace

TEST_CASE("pardiso: the runtime loader and the configuration checks", "[pardiso]") {
    // Without a runtime loaded: tests/linalg/mkl_absent_probe.cpp, in a process of its own.
    REQUIRE_MKL();
    const auto runtime = loaded_mkl_runtime();
    REQUIRE(runtime.has_value());
    REQUIRE(runtime->version.starts_with("Intel(R) oneAPI Math Kernel Library"));
    // A second load returns the runtime already loaded, whatever the path.
    REQUIRE(load_mkl_runtime("elsewhere/mkl_rt.dll")->path == runtime->path);
    for (const int threads : {0, -1, 1025}) {
        REQUIRE(LinearSolver::create({.backend = SolverBackend::mkl_pardiso, .threads = threads})
                    .error()
                    .code == ErrorCode::invalid_input);
    }
    REQUIRE(LinearSolver::create({.backend = SolverBackend::mkl_pardiso, .threads = 8}).has_value());
    REQUIRE(make().config().backend == SolverBackend::mkl_pardiso);
}

TEST_CASE("pardiso: agrees with Eigen SparseLU on the shared systems", "[pardiso]") {
    REQUIRE_MKL();
    const auto check = [](const char* name, const SparseMatrix& a, std::span<const double> b,
                          std::span<const double> exact, double tolerance) {
        const Solved e = solve_with(a, b, {});
        const Solved p = solve_with(a, b, pardiso);
        const double difference = max_abs_diff(p.x, e.x) / max_abs(e.x);
        std::printf("pardiso %-22s: |x_p - x_e| / |x_e| %.2e, backward error %.2e (Eigen %.2e), "
                    "perturbed %zu, error vs exact %.2e\n",
                    name, difference, p.solve.backward_error, e.solve.backward_error,
                    *p.factorization.perturbed_pivots, relative_error(p.x, exact));
        REQUIRE(!p.factorization.pivot_ratio.has_value());
        REQUIRE(*p.factorization.perturbed_pivots == 0);
        REQUIRE(difference <= tolerance);
        REQUIRE(relative_error(p.x, exact) <= tolerance);
    };
    {
        const System s = nonsymmetric(2000);
        check("nonsymmetric", s.a, s.b, s.x_exact, 1e-13);
    }
    {
        const SparseMatrix a = laplace_1d(random_conductances(999), true);
        std::vector<double> exact(1000), b(1000);
        for (std::size_t i = 0; i < exact.size(); ++i) exact[i] = 1.0 + hash01(70'000 + i);
        a.multiply(exact, b);
        check("contacted Laplacian", a, b, exact, 1e-9);
    }
    {
        const System s = badly_scaled_coupled(20);
        check("badly scaled coupled", s.a, s.b, s.x_exact, 1e-10);
    }
    // No equilibration on our side: PARDISO's own scaling and matching still solve it.
    {
        const System s = badly_scaled_coupled(20);
        const Solved p = solve_with(s.a, s.b,
                                    {.backend = SolverBackend::mkl_pardiso, .equilibrate = false,
                                     .min_pivot_ratio = 0.0});
        std::printf("pardiso unequilibrated  : backward error %.2e, error vs exact %.2e\n",
                    p.solve.backward_error, relative_error(p.x, s.x_exact));
        REQUIRE(relative_error(p.x, s.x_exact) <= 1e-10);
    }
}

TEST_CASE("pardiso: singular systems are singular_system, and the solver recovers", "[pardiso]") {
    REQUIRE_MKL();
    LinearSolver solver = make();
    const SparseMatrix singular = from(2, {{0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 2.0}, {1, 1, 4.0}});
    REQUIRE(solver.factorize(singular).error().code == ErrorCode::singular_system);
    REQUIRE_FALSE(solver.is_factorized());
    const SparseMatrix regular = from(2, {{0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 2.0}, {1, 1, 5.0}});
    REQUIRE(solver.factorize(regular).has_value());
    std::array<double, 2> x{};
    REQUIRE(solver.solve(std::array<double, 2>{3.0, 7.0}, x).has_value());
    REQUIRE(max_abs_diff(x, std::array<double, 2>{1.0, 1.0}) <= 1e-15);

    // Floating regions: exact (unit conductances) and at rounding level (random ones), alone and
    // among a contacted chain.
    for (const bool unit : {true, false}) {
        const SparseMatrix a =
            laplace_1d(unit ? std::vector<double>(999, 1.0) : random_conductances(999), false);
        const auto r = solver.factorize(a);
        REQUIRE(r.error().code == ErrorCode::singular_system);
    }
    {
        std::vector<Triplet> t;
        const std::vector<double> g = random_conductances(1000);
        for (Index i = 0; i + 1 < 1000; ++i) {
            if (i == 499) continue;
            const double w = g[static_cast<std::size_t>(i)];
            t.insert(t.end(), {{i, i, w}, {i + 1, i + 1, w}, {i, i + 1, -w}, {i + 1, i, -w}});
        }
        t.push_back({0, 0, 1.0});
        REQUIRE(solver.factorize(from(1000, t)).error().code == ErrorCode::singular_system);
    }
    // With the check disabled, PARDISO perturbs the pivot instead of failing, and says so.
    LinearSolver unchecked = make({.backend = SolverBackend::mkl_pardiso, .min_pivot_ratio = 0.0});
    const auto r = unchecked.factorize(laplace_1d(random_conductances(999), false));
    REQUIRE(r.has_value());
    REQUIRE(*r->perturbed_pivots >= 1);
    REQUIRE(!r->pivot_ratio.has_value());
}

TEST_CASE("pardiso: weakly anchored valid regions are accepted", "[pardiso]") {
    // The Eigen suite's known false positives (valid systems whose pivot ratio falls under
    // min_pivot_ratio: a zeroth-order anchor delta = 1e-14, a weak link w = 1e-12) are not
    // perturbed by PARDISO and solve accurately, with the default configuration.
    REQUIRE_MKL();
    for (const auto& [delta, w] : {std::pair{1e-12, 0.0}, std::pair{1e-14, 0.0},
                                   std::pair{0.0, 1e-9}, std::pair{0.0, 1e-12}}) {
        LinearSolver solver = make();
        const auto f = solver.factorize(anchored_chain(1000, delta, w));
        CAPTURE(delta, w);
        REQUIRE(f.has_value());
        REQUIRE(*f->perturbed_pivots == 0);
        std::vector<double> b(1000, 1.0), x(1000);
        const auto s = solver.solve(b, x);
        REQUIRE(s.has_value());
        REQUIRE(s->backward_error <= 1e-15);
    }
}

TEST_CASE("pardiso: the reuse gate, one analysis, one factorization per matrix, one solve per "
          "right-hand side",
          "[pardiso]") {
    REQUIRE_MKL();
    LinearSolver solver = make();
    const SparseMatrix base = tridiagonal(500, -1.0, 4.0, -1.5);
    std::vector<double> b(500, 1.0), x(500);
    for (int k = 0; k < 4; ++k) {
        SparseMatrix a = base;
        for (double& v : a.values()) v *= 1.0 + 0.1 * k;  // new values, same pattern
        REQUIRE(solver.factorize(a).has_value());
        for (int r = 0; r < 3; ++r) {
            b[static_cast<std::size_t>(r)] += 1.0;
            const auto s = solver.solve(b, x);
            REQUIRE(s.has_value());
            REQUIRE(s->refinement_steps == 0);
        }
    }
    REQUIRE(solver.analyses() == 1);
    REQUIRE(solver.factorizations() == 4);
    REQUIRE(solver.backend_counts().analyses == 1);
    REQUIRE(solver.backend_counts().factorizations == 4);
    REQUIRE(solver.backend_counts().solves == 12);

    // A changed pattern is analyzed again, once.
    REQUIRE(solver.factorize(tridiagonal(500, -1.0, 4.0, 0.0)).has_value());  // explicit zeros
    REQUIRE(solver.factorize(nonsymmetric(300).a).has_value());
    REQUIRE(solver.factorize(nonsymmetric(300).a).has_value());
    REQUIRE(solver.backend_counts().analyses == 2);  // the explicit zeros kept the pattern
    REQUIRE(solver.backend_counts().factorizations == 7);

    // analyze() alone does no backend work: PARDISO's analysis needs values (matching, scaling).
    LinearSolver fresh = make();
    REQUIRE(fresh.analyze(base).has_value());
    REQUIRE(fresh.backend_counts().analyses == 0);
    REQUIRE(fresh.factorize(base).has_value());
    REQUIRE(fresh.backend_counts().analyses == 1);
    REQUIRE(fresh.backend_counts().factorizations == 1);
}

TEST_CASE("pardiso: the same counts for Eigen SparseLU", "[pardiso]") {
    // The counters are backend-neutral: Eigen's analyzePattern, factorize and solve.
    LinearSolver solver = LinearSolver::create({}).value();
    const SparseMatrix a = tridiagonal(50, -1.0, 4.0, -1.5);
    std::vector<double> b(50, 1.0), x(50);
    for (int k = 0; k < 3; ++k) {
        REQUIRE(solver.factorize(a).has_value());
        REQUIRE(solver.solve(b, x).has_value());
    }
    REQUIRE(solver.backend_counts().analyses == 1);
    REQUIRE(solver.backend_counts().factorizations == 3);
    REQUIRE(solver.backend_counts().solves == 3);
}

TEST_CASE("pardiso: complex systems", "[pardiso]") {
    REQUIRE_MKL();
    // A known nonsymmetric system.
    {
        const Index n = 400;
        std::vector<ComplexTriplet> t;
        for (Index i = 0; i < n; ++i) {
            if (i > 0) t.push_back({i, i - 1, Complex{-1.3, 0.4}});
            t.push_back({i, i, Complex{2.5, 1.0 + 0.001 * i}});
            if (i + 1 < n) t.push_back({i, i + 1, Complex{-0.7, -0.2}});
        }
        const ComplexSparseMatrix a = complex_from(n, t);
        std::vector<Complex> exact, b(static_cast<std::size_t>(n)), x(b.size());
        for (Index i = 0; i < n; ++i) exact.push_back({1.0 + std::sin(0.01 * i), std::cos(0.07 * i)});
        a.multiply(exact, b);
        auto solver = make<ComplexLinearSolver>();
        const auto f = solver.factorize(a);
        REQUIRE(f.has_value());
        REQUIRE(!f->pivot_ratio.has_value());  // a complex factor's pivots are not read
        REQUIRE(*f->perturbed_pivots == 0);
        const auto s = solver.solve(b, x);
        REQUIRE(s.has_value());
        REQUIRE(s->backward_error < 4 * std::numeric_limits<double>::epsilon());
        REQUIRE(max_abs_diff(x, exact) / max_abs(exact) < 1e-14);
    }
    // G + i omega C against Eigen, over frequencies, reusing the analysis.
    {
        const SparseMatrix lap = laplace_1d(random_conductances(299), true);
        const auto n = static_cast<std::size_t>(lap.rows());
        std::vector<double> c(n);
        for (std::size_t i = 0; i < n; ++i) c[i] = 0.5 + hash01(5'000 + i);
        std::vector<Complex> b(n), xp(n), xe(n);
        for (std::size_t i = 0; i < n; ++i) b[i] = {hash01(i), hash01(10'000 + i) - 0.5};
        auto p = make<ComplexLinearSolver>();
        auto e = ComplexLinearSolver::create({}).value();
        double worst = 0.0;
        for (const double omega : {1e-6, 0.3, 30.0, 1e6}) {
            const ComplexSparseMatrix a = ac_matrix(lap, c, omega);
            REQUIRE(p.factorize(a).has_value());
            REQUIRE(e.factorize(a).has_value());
            REQUIRE(p.solve(b, xp).has_value());
            REQUIRE(e.solve(b, xe).has_value());
            worst = std::max(worst, max_abs_diff(xp, xe) / max_abs(xe));
        }
        std::printf("pardiso complex G + i omega C: worst |x_p - x_e| / |x_e| %.2e\n", worst);
        REQUIRE(worst < 1e-12);
        REQUIRE(p.backend_counts().analyses == 1);
        REQUIRE(p.backend_counts().factorizations == 4);
    }
    // A floating region (i omega times a pure-Neumann Laplacian): perturbed pivots, singular.
    {
        const SparseMatrix lap = laplace_1d(random_conductances(99), false);
        auto solver = make<ComplexLinearSolver>();
        REQUIRE(solver.factorize(to_complex(lap, 7.0)).error().code == ErrorCode::singular_system);
        auto unchecked = make<ComplexLinearSolver>(
            {.backend = SolverBackend::mkl_pardiso, .min_pivot_ratio = 0.0});
        const auto f = unchecked.factorize(to_complex(lap, 7.0));
        REQUIRE(f.has_value());
        REQUIRE(*f->perturbed_pivots >= 1);
    }
}

TEST_CASE("pardiso: thread counts agree within rounding, and solvers run concurrently",
          "[pardiso]") {
    REQUIRE_MKL();
    const System s = badly_scaled_coupled(40);
    const Solved one = solve_with(s.a, s.b, pardiso);
    const Solved four = solve_with(s.a, s.b, {.backend = SolverBackend::mkl_pardiso, .threads = 4});
    const double difference = max_abs_diff(four.x, one.x) / max_abs(one.x);
    std::printf("pardiso 4 threads vs 1: %.2e\n", difference);
    REQUIRE(difference <= 1e-10);
    // Repeated runs at one thread count agree within the same tolerance (not required bit for bit).
    const Solved again = solve_with(s.a, s.b, pardiso);
    REQUIRE(max_abs_diff(again.x, one.x) / max_abs(one.x) <= 1e-10);

    std::vector<double> first(s.b.size()), second(s.b.size());
    bool ok_first = false, ok_second = false;
    {
        const auto run = [&s](std::vector<double>& x, bool& ok, int threads) {
            LinearSolver solver =
                LinearSolver::create({.backend = SolverBackend::mkl_pardiso, .threads = threads})
                    .value();
            ok = solver.factorize(s.a).has_value() && solver.solve(s.b, x).has_value();
        };
        std::jthread t1([&] { run(first, ok_first, 1); });
        std::jthread t2([&] { run(second, ok_second, 2); });
    }
    REQUIRE(ok_first);
    REQUIRE(ok_second);
    REQUIRE(max_abs_diff(first, one.x) / max_abs(one.x) <= 1e-10);
    REQUIRE(max_abs_diff(second, one.x) / max_abs(one.x) <= 1e-10);
}

TEST_CASE("pardiso: invalid values and right-hand sides are errors before PARDISO runs",
          "[pardiso]") {
    REQUIRE_MKL();
    LinearSolver solver = make();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(solver.factorize(from(2, {{0, 0, nan}, {1, 1, 1.0}})).error().code ==
            ErrorCode::invalid_input);
    REQUIRE(solver.factorize(from(3, {{0, 0, 1.0}, {1, 0, 1.0}, {2, 2, 1.0}})).error().code ==
            ErrorCode::singular_system);
    REQUIRE(solver.backend_counts().analyses == 0);
    REQUIRE(solver.factorize(from(2, {{0, 0, 2.0}, {1, 1, 4.0}})).has_value());
    std::array<double, 2> x{};
    REQUIRE(solver.solve(std::array<double, 2>{1.0, nan}, x).error().code ==
            ErrorCode::invalid_input);
    REQUIRE(solver.solve(std::array<double, 2>{1.0, 2.0}, x).has_value());
    REQUIRE(x[0] == 0.5);
    REQUIRE(x[1] == 0.5);
}
