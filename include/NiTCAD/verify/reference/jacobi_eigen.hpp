// Small Hermitian eigenproblems (n <= 16) by the cyclic Jacobi method, in double (Unit 55, TVR-1, TC9).
//
// Source: G. E. Forsythe and P. Henrici, "The cyclic Jacobi method for computing the principal values of a complex
// matrix", Trans. Amer. Math. Soc. 94 (1960) 1-23, read in a scan (pages.mtu.edu, 2026-11-04): section 1, 2.1
// (Lemma 1, (10)-(14)), 2.5 (the necessity of (IV), example (42)) and 3.1 (Theorems 5 and 6); the proofs of 2.2-2.4
// and the principal-value case 3.2 were not read. The frozen citation, Golub and Van Loan section 8.5, was not
// available. Grade PA on Forsythe and Henrici for the rotation and its convergence.
//
// A rotation in the (p, q) plane is FH (6) with alpha = delta = 0, beta = theta, gamma = -theta, where
// a_pq = r e^{i theta}: U_pp = U_qq = c, U_pq = e^{i theta} s, U_qp = -e^{-i theta} s, which satisfies (5) and (43);
// A' = U A U^*. Then a'_pq = e^{i theta} [r cos 2phi - (a_pp - a_qq) sin 2phi / 2], zero by (44),
// tan 2phi = 2r / (a_pp - a_qq). phi is taken in [-pi/4, pi/4] (the smaller root t = tan phi of
// t^2 + 2 tau t - 1 = 0, tau = (a_pp - a_qq) / 2r), which is condition (IV) with J = [-pi/4, pi/4]; Theorem 3 then
// gives convergence of the cyclic method. Pairs are taken cyclic by rows (III_r). A zero a_pq is skipped (FH 2.5:
// rotating a zero element by pi/2 is what makes (42) cycle). After each rotation a'_pq is set to 0 and the
// diagonal to its real part (a Hermitian matrix's). Sweeps stop when off(A) = (sum_{p != q} |a_pq|^2)^{1/2} <=
// eps ||A||_F, eps = 2^-52.
//
// The error estimate is the final off(A): by Weyl's theorem each eigenvalue of the final matrix lies within
// off(A) of its diagonal. It does not include the rounding of the rotations (the frozen criterion bounds that by
// 8 n eps ||A||_2).
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numeric>
#include <vector>

namespace NiTCAD::verify::reference {

inline constexpr std::size_t jacobi_max_n = 16;
inline constexpr int jacobi_max_sweeps = 50;

enum class JacobiStatus {
    ok,
    invalid_input,   // n == 0, n > 16, a.size() != n * n, a non-finite entry, or A not exactly Hermitian
    not_converged,   // off(A) > eps ||A||_F after jacobi_max_sweeps sweeps
};

struct JacobiResult {
    std::vector<double> eigenvalues;                 // ascending
    std::vector<std::complex<double>> eigenvectors;  // row-major n x n: column j is the eigenvector of value j
    double off_norm = 0.0;                           // final off(A): the deflation error estimate
    int sweeps = 0;
    JacobiStatus status = JacobiStatus::ok;
};

// a: row-major n x n, exactly Hermitian (a_ij == conj(a_ji), Im a_ii == 0).
[[nodiscard]] inline JacobiResult jacobi_eigen(std::vector<std::complex<double>> a, std::size_t n) {
    using cd = std::complex<double>;
    JacobiResult out;
    if (n == 0 || n > jacobi_max_n || a.size() != n * n) {
        out.status = JacobiStatus::invalid_input;
        return out;
    }
    double fro2 = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const cd x = a[i * n + j];
            if (!std::isfinite(x.real()) || !std::isfinite(x.imag()) || x != std::conj(a[j * n + i])) {
                out.status = JacobiStatus::invalid_input;
                return out;
            }
            fro2 += std::norm(x);
        }
    const double eps = 0x1p-52;
    const double target = eps * std::sqrt(fro2);

    std::vector<cd> v(n * n, cd(0.0, 0.0));
    for (std::size_t i = 0; i < n; ++i) v[i * n + i] = cd(1.0, 0.0);

    const auto off = [&] {
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                if (i != j) s += std::norm(a[i * n + j]);
        return std::sqrt(s);
    };

    double o = off();
    while (o > target) {
        if (out.sweeps == jacobi_max_sweeps) {
            out.status = JacobiStatus::not_converged;
            break;
        }
        ++out.sweeps;
        for (std::size_t p = 0; p + 1 < n; ++p)
            for (std::size_t q = p + 1; q < n; ++q) {
                const cd apq = a[p * n + q];
                const double r = std::abs(apq);
                if (r == 0.0) continue;
                const cd w = apq / r;  // e^{i theta}
                const double tau = (a[p * n + p].real() - a[q * n + q].real()) / (2.0 * r);
                // If tau^2 overflows, t = 0: a rotation by 0, whose error r < |a_pp - a_qq| 1e-154 is below rounding.
                const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
                const double c = 1.0 / std::sqrt(1.0 + t * t), s = t * c;
                const cd ws = w * s, wbs = std::conj(w) * s;
                // Rows: A <- U A.
                for (std::size_t k = 0; k < n; ++k) {
                    const cd x = a[p * n + k], y = a[q * n + k];
                    a[p * n + k] = c * x + ws * y;
                    a[q * n + k] = -wbs * x + c * y;
                }
                // Columns: A <- A U^*, and V <- V U^*.
                for (std::size_t k = 0; k < n; ++k) {
                    const cd x = a[k * n + p], y = a[k * n + q];
                    a[k * n + p] = c * x + wbs * y;
                    a[k * n + q] = -ws * x + c * y;
                    const cd vx = v[k * n + p], vy = v[k * n + q];
                    v[k * n + p] = c * vx + wbs * vy;
                    v[k * n + q] = -ws * vx + c * vy;
                }
                a[p * n + q] = a[q * n + p] = cd(0.0, 0.0);
                a[p * n + p] = cd(a[p * n + p].real(), 0.0);
                a[q * n + q] = cd(a[q * n + q].real(), 0.0);
            }
        o = off();
    }
    out.off_norm = o;

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t i, std::size_t j) { return a[i * n + i].real() < a[j * n + j].real(); });
    out.eigenvalues.resize(n);
    out.eigenvectors.resize(n * n);
    for (std::size_t j = 0; j < n; ++j) {
        out.eigenvalues[j] = a[order[j] * n + order[j]].real();
        for (std::size_t i = 0; i < n; ++i) out.eigenvectors[i * n + j] = v[i * n + order[j]];
    }
    return out;
}

}  // namespace NiTCAD::verify::reference
