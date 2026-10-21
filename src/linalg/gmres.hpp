// Restarted GMRES(m) with right preconditioning (Unit 44, LS2). Private to linalg: the kernel of the
// gmres_ilut backend, written against callables so it holds no matrix type and needs no Eigen.
//
// It solves A y = c from y = 0. Each cycle builds an orthonormal Krylov basis of A M^-1 by modified
// Gram-Schmidt with one full reorthogonalization pass ("twice is enough": without it the basis
// loses orthogonality as the residual falls by the condition number, and GMRES stalls), reduces the
// Hessenberg matrix with Givens rotations, and tracks the 2-norm of the residual as it goes. That
// estimate only decides when a cycle ends early; acceptance is always on the caller's measure of the
// true iterate (the componentwise backward error of x on the original system), taken at the end of
// every cycle. Right preconditioning makes the estimate the 2-norm of the true residual of A, not of
// M^-1 A, so the cycle ends near the point where the measure can pass: the estimate is aimed at the
// residual scaled by max_backward_error / measure, the factor the measure must still fall by if it
// is proportional to the residual, with a safety factor of 1/2. A cycle that ends on that estimate
// but fails the measure restarts from the true residual with a new target, so a wrong aim costs a
// restart, never a wrong answer.
//
// Single-threaded and with a fixed operation order, so a solve is bit-for-bit reproducible on one
// machine (6.8).
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace NiTCAD::linalg::detail {

struct GmresSettings {
    std::size_t restart;          // m: basis vectors per cycle, at least 1
    std::size_t max_iterations;   // products with A over all cycles, at least 1
    double tolerance;             // accepted value of the measure
};

struct GmresResult {
    bool converged;
    std::size_t iterations;  // products with A
    double measure;          // of the final iterate: at most tolerance when converged
};

namespace gmres_ops {

inline double conj(double v) noexcept { return v; }
inline std::complex<double> conj(std::complex<double> v) noexcept { return std::conj(v); }

// sum conj(u_i) v_i
template <MatrixScalar Scalar>
Scalar dot(std::span<const Scalar> u, std::span<const Scalar> v) noexcept {
    Scalar s{};
    for (std::size_t i = 0; i < u.size(); ++i) {
        s += conj(u[i]) * v[i];
    }
    return s;
}

// 2-norm scaled by the largest modulus, so huge or tiny vectors neither overflow nor underflow.
template <MatrixScalar Scalar>
double norm2(std::span<const Scalar> v) noexcept {
    double largest = 0.0;
    for (const Scalar x : v) {
        largest = std::max(largest, std::abs(x));
    }
    if (largest == 0.0 || !std::isfinite(largest)) {
        return largest;
    }
    double sum = 0.0;
    for (const Scalar x : v) {
        const double a = std::abs(x) / largest;
        sum += a * a;
    }
    return largest * std::sqrt(sum);
}

// The rotation G = [c, s; -conj(s), c] with real c that maps (a, b) to (rho, 0).
template <MatrixScalar Scalar>
struct Givens {
    double c;
    Scalar s;
};

template <MatrixScalar Scalar>
Givens<Scalar> givens(Scalar a, Scalar b) noexcept {
    const double abs_a = std::abs(a);
    const double abs_b = std::abs(b);
    if (abs_b == 0.0) {
        return {1.0, Scalar{}};
    }
    if (abs_a == 0.0) {
        return {0.0, conj(b) / abs_b};
    }
    const double r = std::hypot(abs_a, abs_b);
    return {abs_a / r, (a / abs_a) * conj(b) / r};
}

template <MatrixScalar Scalar>
void rotate(const Givens<Scalar>& g, Scalar& x, Scalar& y) noexcept {
    const Scalar top = g.c * x + g.s * y;
    y = -conj(g.s) * x + g.c * y;
    x = top;
}

}  // namespace gmres_ops

// Working storage, kept between solves so a run of solves allocates once.
template <MatrixScalar Scalar>
class Gmres {
public:
    // apply(v, w):         w = A v
    // precondition(v, w):  w = M^-1 v
    // measure(y, r):       returns the acceptance measure of y and writes r = c - A y
    // y is overwritten: it starts at zero and holds the last iterate on return.
    template <class Apply, class Precondition, class Measure>
    GmresResult solve(const GmresSettings& settings, Apply&& apply, Precondition&& precondition,
                      Measure&& measure, std::span<Scalar> y);

private:
    std::vector<Scalar> basis_;      // (m + 1) vectors of length n, one after the other
    std::vector<Scalar> hessenberg_; // (m + 1) x m, column-major
    std::vector<gmres_ops::Givens<Scalar>> rotations_;
    std::vector<Scalar> g_;          // rotated right-hand side of the least-squares problem
    std::vector<Scalar> t_;          // its solution
    std::vector<Scalar> residual_;
    std::vector<Scalar> z_;
    std::vector<Scalar> w_;
};

template <MatrixScalar Scalar>
template <class Apply, class Precondition, class Measure>
GmresResult Gmres<Scalar>::solve(const GmresSettings& settings, Apply&& apply,
                                 Precondition&& precondition, Measure&& measure,
                                 std::span<Scalar> y) {
    using gmres_ops::dot;
    using gmres_ops::norm2;
    const std::size_t n = y.size();
    const std::size_t m = settings.restart;
    basis_.resize((m + 1) * n);
    hessenberg_.resize((m + 1) * m);
    rotations_.resize(m);
    g_.resize(m + 1);
    t_.resize(m);
    residual_.resize(n);
    z_.resize(n);
    w_.resize(n);
    const auto vector = [&](std::size_t i) { return std::span<Scalar>(basis_).subspan(i * n, n); };
    const auto h = [&](std::size_t i, std::size_t j) -> Scalar& { return hessenberg_[j * (m + 1) + i]; };

    std::ranges::fill(y, Scalar{});
    double omega = measure(std::span<const Scalar>(y), std::span<Scalar>(residual_));
    std::size_t iterations = 0;
    while (!(omega <= settings.tolerance)) {
        const double beta = norm2<Scalar>(residual_);
        if (iterations >= settings.max_iterations || !std::isfinite(omega) || beta == 0.0 ||
            !std::isfinite(beta)) {
            return {false, iterations, omega};
        }
        // A new cycle from the true residual: nothing of the previous basis is kept.
        const std::span<Scalar> v0 = vector(0);
        for (std::size_t i = 0; i < n; ++i) {
            v0[i] = residual_[i] / beta;
        }
        std::ranges::fill(g_, Scalar{});
        g_[0] = beta;
        const double target = 0.5 * beta * settings.tolerance / omega;

        std::size_t columns = 0;
        for (std::size_t j = 0; j < m && iterations < settings.max_iterations; ++j) {
            precondition(std::span<const Scalar>(vector(j)), std::span<Scalar>(z_));
            apply(std::span<const Scalar>(z_), std::span<Scalar>(w_));
            ++iterations;
            for (std::size_t i = 0; i <= j + 1; ++i) {
                h(i, j) = Scalar{};
            }
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t i = 0; i <= j; ++i) {
                    const std::span<const Scalar> vi = vector(i);
                    const Scalar hij = dot<Scalar>(vi, w_);
                    h(i, j) += hij;
                    for (std::size_t k = 0; k < n; ++k) {
                        w_[k] -= hij * vi[k];
                    }
                }
            }
            const double next = norm2<Scalar>(w_);
            h(j + 1, j) = next;
            for (std::size_t i = 0; i < j; ++i) {
                gmres_ops::rotate(rotations_[i], h(i, j), h(i + 1, j));
            }
            rotations_[j] = gmres_ops::givens(h(j, j), h(j + 1, j));
            gmres_ops::rotate(rotations_[j], h(j, j), h(j + 1, j));
            if (h(j, j) == Scalar{}) {
                break;  // A M^-1 is singular on this Krylov space: the column adds nothing
            }
            gmres_ops::rotate(rotations_[j], g_[j], g_[j + 1]);
            columns = j + 1;
            if (!(next > 0.0) || !(std::abs(g_[j + 1]) > target)) {
                break;  // exact (breakdown) or estimated to pass: end the cycle and measure
            }
            const std::span<Scalar> vj1 = vector(j + 1);
            for (std::size_t k = 0; k < n; ++k) {
                vj1[k] = w_[k] / next;
            }
        }
        if (columns == 0) {
            return {false, iterations, omega};
        }
        // Back substitution on the triangular R t = g, then y += M^-1 V t.
        for (std::size_t i = columns; i-- > 0;) {
            Scalar s = g_[i];
            for (std::size_t k = i + 1; k < columns; ++k) {
                s -= h(i, k) * t_[k];
            }
            t_[i] = s / h(i, i);
        }
        std::ranges::fill(w_, Scalar{});
        for (std::size_t i = 0; i < columns; ++i) {
            const std::span<const Scalar> vi = vector(i);
            for (std::size_t k = 0; k < n; ++k) {
                w_[k] += t_[i] * vi[k];
            }
        }
        precondition(std::span<const Scalar>(w_), std::span<Scalar>(z_));
        for (std::size_t k = 0; k < n; ++k) {
            y[k] += z_[k];
        }
        omega = measure(std::span<const Scalar>(y), std::span<Scalar>(residual_));
    }
    return {true, iterations, omega};
}

}  // namespace NiTCAD::linalg::detail
