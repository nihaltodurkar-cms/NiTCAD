// The GMRES kernel of the gmres_ilut backend (Unit 44, LG2), on its own: without a preconditioner
// and with simple ones, so its iteration counts can be compared with what GMRES guarantees in
// exact arithmetic. This file includes the private src/linalg/gmres.hpp (it holds no Eigen and no
// solver state); everything else in the backend is tested through the public interface in
// gmres_ilut_test.cpp.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include "../../src/linalg/gmres.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "fixtures.hpp"

using namespace NiTCAD::linalg;
using namespace fixtures;
using detail::Gmres;
using detail::GmresResult;
using detail::GmresSettings;

namespace {

// Componentwise backward error max_i |c - A y|_i / (|A| |y| + |c|)_i (Oettli-Prager), written
// independently of the solver's.
template <MatrixScalar Scalar>
double backward_error(const BasicSparseMatrix<Scalar>& a, std::span<const Scalar> c,
                      std::span<const Scalar> y, std::span<Scalar> r) {
    double worst = 0.0;
    for (Index i = 0; i < a.rows(); ++i) {
        Scalar ay{};
        double abs_ay = 0.0;
        for (Index k = a.row_offsets()[i]; k < a.row_offsets()[i + 1]; ++k) {
            const Scalar term = a.values()[k] * y[static_cast<std::size_t>(a.col_indices()[k])];
            ay += term;
            abs_ay += std::abs(term);
        }
        const auto u = static_cast<std::size_t>(i);
        r[u] = c[u] - ay;
        const double denominator = abs_ay + std::abs(c[u]);
        const double row = denominator > 0.0 ? std::abs(r[u]) / denominator
                           : r[u] == Scalar{} ? 0.0
                                              : HUGE_VAL;
        worst = std::isnan(row) ? row : std::max(worst, row);
    }
    return worst;
}

// GMRES on A y = c with a diagonal right preconditioner M (identity when empty); counts the
// measures taken (one at the start and one per cycle).
template <MatrixScalar Scalar>
struct Run {
    GmresResult result;
    std::vector<Scalar> y;
    std::size_t measures = 0;
};

template <MatrixScalar Scalar>
Run<Scalar> run(const BasicSparseMatrix<Scalar>& a, std::span<const Scalar> c,
                GmresSettings settings, const std::vector<Scalar>& m_diag = {}) {
    Run<Scalar> out;
    out.y.resize(c.size());
    Gmres<Scalar> gmres;
    out.result = gmres.solve(
        settings, [&](std::span<const Scalar> v, std::span<Scalar> w) { a.multiply(v, w); },
        [&](std::span<const Scalar> v, std::span<Scalar> w) {
            for (std::size_t i = 0; i < v.size(); ++i) {
                w[i] = m_diag.empty() ? v[i] : v[i] / m_diag[i];
            }
        },
        [&](std::span<const Scalar> y, std::span<Scalar> r) {
            ++out.measures;
            return backward_error(a, c, y, r);
        },
        std::span<Scalar>(out.y));
    return out;
}

// Dense A = P D P^-1 with P unit lower bidiagonal (subdiagonal p), so A is diagonalizable with the
// eigenvalues of D and nonnormal; P^-1 has entries (-p)^(i-j) below the diagonal.
template <MatrixScalar Scalar>
BasicSparseMatrix<Scalar> similar_to_diagonal(const std::vector<Scalar>& d, double p) {
    const auto n = d.size();
    std::vector<std::vector<Scalar>> pd(n, std::vector<Scalar>(n));  // P D
    for (std::size_t i = 0; i < n; ++i) {
        pd[i][i] = d[i];
        if (i > 0) pd[i][i - 1] = p * d[i - 1];
    }
    std::vector<BasicTriplet<Scalar>> t;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            Scalar s{};
            for (std::size_t k = j; k <= i; ++k) {
                s += pd[i][k] * std::pow(-p, static_cast<double>(k - j));
            }
            t.push_back({static_cast<Index>(i), static_cast<Index>(j), s});
        }
    }
    return BasicSparseMatrix<Scalar>::from_triplets(static_cast<Index>(n), static_cast<Index>(n), t)
        .value();
}

// 2D convection-diffusion on an m x m grid, upwinded: nonsymmetric, needs many Krylov vectors.
SparseMatrix convection_diffusion(Index m, double peclet) {
    std::vector<Triplet> t;
    for (Index j = 0; j < m; ++j) {
        for (Index i = 0; i < m; ++i) {
            const Index k = j * m + i;
            t.push_back({k, k, 4.0 + peclet});
            if (i > 0) t.push_back({k, k - 1, -1.0 - peclet});
            if (i + 1 < m) t.push_back({k, k + 1, -1.0});
            if (j > 0) t.push_back({k, k - m, -1.0});
            if (j + 1 < m) t.push_back({k, k + m, -1.0});
        }
    }
    return from(m * m, t);
}

}  // namespace

TEST_CASE("gmres kernel: k distinct eigenvalues need at most k iterations") {
    // In exact arithmetic GMRES on a diagonalizable matrix terminates in at most as many steps as
    // the matrix has distinct eigenvalues. 60 unknowns, each of 6 eigenvalues 10 times, real and
    // complex, no preconditioner. The eigenvalues lie within one decade: over four (1 to 1e4) the
    // minimal polynomial's terms reach 1e19 before they cancel (measured), so termination at step 6
    // is below rounding and floating-point GMRES needs a seventh step.
    std::vector<double> d;
    std::vector<std::complex<double>> dc;
    const std::vector<double> values{1.0, 2.0, 3.0, 5.0, 8.0, 13.0};
    for (std::size_t i = 0; i < 60; ++i) {
        d.push_back(values[i % values.size()]);
        // i % 2 is fixed by i % 6, so the complex matrix also has 6 distinct eigenvalues.
        dc.push_back({values[i % values.size()], i % 2 == 0 ? 0.5 : -2.0});
    }
    const SparseMatrix a = similar_to_diagonal(d, 0.5);
    std::vector<double> c(60);
    for (std::size_t i = 0; i < c.size(); ++i) c[i] = 1.0 + hash01(300 + i);
    const auto r = run<double>(a, c, {.restart = 50, .max_iterations = 50, .tolerance = 1e-8});
    CAPTURE(r.result.iterations, r.result.measure);
    REQUIRE(r.result.converged);
    REQUIRE(r.result.iterations <= 6);
    REQUIRE(r.result.measure <= 1e-8);

    const ComplexSparseMatrix ac = similar_to_diagonal(dc, 0.5);
    std::vector<std::complex<double>> cc(60);
    for (std::size_t i = 0; i < cc.size(); ++i) cc[i] = {hash01(400 + i), hash01(500 + i) - 0.5};
    const auto rc =
        run<std::complex<double>>(ac, cc, {.restart = 50, .max_iterations = 50, .tolerance = 1e-8});
    CAPTURE(rc.result.iterations, rc.result.measure);
    REQUIRE(rc.result.converged);
    REQUIRE(rc.result.iterations <= 6);
}

TEST_CASE("gmres kernel: unpreconditioned nonsymmetric systems reach the target") {
    const System s = nonsymmetric(400);
    const auto r = run<double>(s.a, s.b, {.restart = 30, .max_iterations = 2000, .tolerance = 1e-12});
    CAPTURE(r.result.iterations, r.result.measure);
    REQUIRE(r.result.converged);
    REQUIRE(r.result.measure <= 1e-12);
    REQUIRE(relative_error(r.y, s.x_exact) <= 1e-9);

    // Restarted many times: the convection-diffusion operator with restart 10.
    const SparseMatrix a = convection_diffusion(30, 2.0);
    std::vector<double> c(900), y_exact(900);
    for (std::size_t i = 0; i < y_exact.size(); ++i) y_exact[i] = 1.0 + hash01(900 + i);
    a.multiply(y_exact, c);
    const auto restarted = run<double>(a, c, {.restart = 10, .max_iterations = 5000, .tolerance = 1e-10});
    CAPTURE(restarted.result.iterations, restarted.measures);
    REQUIRE(restarted.result.converged);
    REQUIRE(restarted.result.iterations > 10);  // several cycles
    REQUIRE(restarted.measures >= 3);
    REQUIRE(relative_error(restarted.y, y_exact) <= 1e-7);

    // Complex: (G + i omega C) of the AC system.
    const SparseMatrix lap = laplace_1d(random_conductances(199), true);
    std::vector<double> cap(200);
    for (std::size_t i = 0; i < cap.size(); ++i) cap[i] = 0.5 + hash01(5'000 + i);
    const ComplexSparseMatrix ac = ac_matrix(lap, cap, 3.0);
    std::vector<Complex> b(200);
    for (std::size_t i = 0; i < b.size(); ++i) b[i] = {hash01(i), hash01(10'000 + i) - 0.5};
    const auto rc = run<Complex>(ac, b, {.restart = 40, .max_iterations = 4000, .tolerance = 1e-10});
    CAPTURE(rc.result.iterations);
    REQUIRE(rc.result.converged);
    REQUIRE(rc.result.measure <= 1e-10);
}

TEST_CASE("gmres kernel: the restart length and the iteration limit are honoured") {
    // A hard system that no cycle of 5 can finish: cycles of 5, 5, 5, 5 and 3 products, then stop.
    const SparseMatrix a = convection_diffusion(30, 2.0);
    const std::vector<double> c(900, 1.0);
    const auto r = run<double>(a, c, {.restart = 5, .max_iterations = 23, .tolerance = 1e-10});
    REQUIRE_FALSE(r.result.converged);
    REQUIRE(r.result.iterations == 23);
    REQUIRE(r.measures == 1 + 5);
    REQUIRE(r.result.measure > 1e-10);
    // The returned iterate is the last one, and its measure is the one reported.
    std::vector<double> residual(900);
    REQUIRE(backward_error<double>(a, c, r.y, residual) == r.result.measure);

    // restart = 1 is steepest-descent-like GMRES(1): one product per cycle.
    const auto one = run<double>(a, c, {.restart = 1, .max_iterations = 7, .tolerance = 1e-10});
    REQUIRE(one.result.iterations == 7);
    REQUIRE(one.measures == 1 + 7);
}

TEST_CASE("gmres kernel: right preconditioning by the exact inverse takes one iteration") {
    // A diagonal with entries over 12 decades: unpreconditioned GMRES needs one step per distinct
    // entry; with M = A (right preconditioning) the first Krylov vector is the solution.
    const Index n = 200;
    std::vector<Triplet> t;
    std::vector<double> m_diag(n), c(n);
    for (Index i = 0; i < n; ++i) {
        const double v = std::pow(10.0, 12.0 * hash01(700 + static_cast<std::size_t>(i)) - 6.0);
        t.push_back({i, i, v});
        m_diag[static_cast<std::size_t>(i)] = v;
        c[static_cast<std::size_t>(i)] = 1.0 + hash01(800 + static_cast<std::size_t>(i));
    }
    const SparseMatrix a = from(n, t);
    const auto r = run<double>(a, c, {.restart = 30, .max_iterations = 30, .tolerance = 1e-12},
                               m_diag);
    REQUIRE(r.result.converged);
    REQUIRE(r.result.iterations == 1);
    for (std::size_t i = 0; i < c.size(); ++i) {
        REQUIRE(std::abs(r.y[i] * m_diag[i] - c[i]) <= 1e-14 * c[i]);
    }
    const auto plain = run<double>(a, c, {.restart = 30, .max_iterations = 30, .tolerance = 1e-12});
    REQUIRE_FALSE(plain.result.converged);

    // Columns scaled by 10^[-4, 4] (A = T D): right preconditioning by D undoes them exactly.
    const System s = nonsymmetric(300);
    std::vector<double> d(300);
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = std::pow(10.0, 8.0 * hash01(900 + i) - 4.0);
    std::vector<Triplet> col_scaled;
    for (Index i = 0; i < 300; ++i) {
        for (Index k = s.a.row_offsets()[i]; k < s.a.row_offsets()[i + 1]; ++k) {
            const Index j = s.a.col_indices()[k];
            col_scaled.push_back({i, j, s.a.values()[k] * d[static_cast<std::size_t>(j)]});
        }
    }
    const SparseMatrix td = from(300, col_scaled);
    const auto pre = run<double>(td, s.b, {.restart = 60, .max_iterations = 60, .tolerance = 1e-12}, d);
    const auto unpre = run<double>(s.a, s.b, {.restart = 60, .max_iterations = 60, .tolerance = 1e-12});
    CAPTURE(pre.result.iterations, unpre.result.iterations);
    REQUIRE(pre.result.converged);
    // T D M^-1 = T: the same Krylov space as the unscaled system, so the same iteration count.
    REQUIRE(pre.result.iterations == unpre.result.iterations);
    for (std::size_t i = 0; i < pre.y.size(); ++i) {
        REQUIRE(std::abs(pre.y[i] * d[i] - s.x_exact[i]) <= 1e-9 * std::abs(s.x_exact[i]));
    }
}

TEST_CASE("gmres kernel: a singular system is not detected by the backward error") {
    // Pure Neumann Laplacian with an inconsistent right-hand side: no y has a small residual, but
    // a y of huge norm has a small componentwise backward error (|A| |y| is huge), as a direct
    // solve without its pivot check finds. GMRES reaches one; the norm is the only sign.
    const SparseMatrix a = laplace_1d(random_conductances(99), false);
    const std::vector<double> c(100, 1.0);
    const auto r = run<double>(a, c, {.restart = 20, .max_iterations = 200, .tolerance = 1e-8});
    CAPTURE(r.result.iterations, max_abs(r.y));
    REQUIRE(r.result.converged);
    REQUIRE(max_abs(r.y) > 1e10);
}

TEST_CASE("gmres kernel: a zero right-hand side is converged with no iteration") {
    const System s = nonsymmetric(50);
    const std::vector<double> c(50, 0.0);
    const auto r = run<double>(s.a, c, {.restart = 10, .max_iterations = 10, .tolerance = 1e-8});
    REQUIRE(r.result.converged);
    REQUIRE(r.result.iterations == 0);
    REQUIRE(r.measures == 1);
    REQUIRE(std::ranges::all_of(r.y, [](double v) { return v == 0.0; }));
}
