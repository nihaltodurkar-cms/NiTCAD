// Chebyshev collocation with domain splitting, in double-double (Unit 55, TVR-1, TC7).
//
// Sources, read before code (2026-11-04): L. N. Trefethen, Spectral Methods in MATLAB (SIAM, 2000), chapter 6
// (the author's sample-chapter PDF, people.maths.ox.ac.uk/trefethen/book.pdf, rendered page by page): Theorem 7,
// (6.3)-(6.5), the program cheb.m and the negative-sum diagonal (6.6), and the printed cheb(1) ... cheb(5); the
// online cheb.m and p13.m (the boundary-value problem u_xx = e^{4x}, u(+-1) = 0, with its exact solution), and the
// errata list (no entry touches chapter 6). J. P. Boyd, Chebyshev and Fourier Spectral Methods, 2nd ed. (Dover,
// 2000), sections 22.1-22.4: patching, (22.3)-(22.5), u and u_x continuous at the interface, the Dirichlet
// conditions at the outer ends. Provenance PA (both texts read).
//
// The differentiation matrix is cheb.m's: off-diagonal entries by (6.5), the diagonal by (6.6). Two deliberate
// changes, both exact identities, so that a double-double result is accurate to double-double: x_j = cos(j pi / N)
// is evaluated as sin(pi (N - 2j) / (2N)) (exactly antisymmetric), and x_i - x_j as -2 sin((i + j) pi / (2N))
// sin((i - j) pi / (2N)) (no cancellation between nearby nodes; cheb.m forms X - X'). sin is a Taylor series on
// [-pi/2, pi/2], internal to this library.
//
// The two-domain problem is -(eps u')' = f on [a, b] with eps = eps_left on [a, c] and eps_right on [c, b]
// (constant on each), u(a) = u_a, u(b) = u_b. Each subdomain carries N + 1 Chebyshev points; the collocation rows
// are -eps (2 / L)^2 (D^2 u)_j = f(x_j) at the interior points; the interface rows are Boyd's patching (22.5),
// with (22.5b) written for the flux, eps_left u_left'(c) = eps_right u_right'(c), which is what integrating the
// equation across x = c gives when eps jumps. The dense system is solved by Gaussian elimination with partial
// pivoting in double-double.
//
// Error estimate (U_ref): the magnitude of the last two Chebyshev coefficients of each subdomain's interpolant,
// the usual heuristic for a geometrically convergent expansion; reported, not a bound.
#pragma once

#include "NiTCAD/verify/reference/double_double.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace NiTCAD::verify::reference {

using DDVector = std::vector<DD>;

struct DDMatrix {
    std::size_t rows = 0, cols = 0;
    std::vector<DD> v;
    DDMatrix() = default;
    DDMatrix(std::size_t r, std::size_t c) : rows(r), cols(c), v(r * c, DD(0.0)) {}
    [[nodiscard]] DD& operator()(std::size_t i, std::size_t j) { return v[i * cols + j]; }
    [[nodiscard]] const DD& operator()(std::size_t i, std::size_t j) const { return v[i * cols + j]; }
};

[[nodiscard]] DDMatrix multiply(const DDMatrix& a, const DDMatrix& b);

// sin(m pi / d) for integers m, d > 0 (m reduced to [0, d / 2] by exact integer arithmetic first).
[[nodiscard]] DD sin_pi_ratio(long long m, long long d);

// x_j = cos(j pi / N), j = 0 ... N (x_0 = 1, as cheb.m).
[[nodiscard]] DDVector chebyshev_points(int n);

// cheb.m: the (N + 1) x (N + 1) differentiation matrix on chebyshev_points(N).
[[nodiscard]] DDMatrix chebyshev_differentiation(int n);

// The coefficients a_0 ... a_N of the interpolant sum a_k T_k(xi) through values at chebyshev_points(N).
[[nodiscard]] DDVector chebyshev_coefficients(const DDVector& values);

// Clenshaw's recurrence for sum a_k T_k(xi), xi in [-1, 1].
[[nodiscard]] DD chebyshev_evaluate(const DDVector& a, const DD& xi);

// Solves A x = b in place (b becomes x) by Gaussian elimination with partial pivoting; false if a pivot is zero.
[[nodiscard]] bool lu_solve(DDMatrix a, DDVector& b);

struct TwoDomainProblem {
    DD a = DD(-1.0), c = DD(0.0), b = DD(1.0);
    DD eps_left = DD(1.0), eps_right = DD(1.0);
    std::function<DD(const DD&)> f;
    DD u_a, u_b;
    int n = 16;  // Chebyshev degree on each subdomain
};

struct TwoDomainSolution {
    bool ok = false;
    DD a, c, b;
    DDVector x_left, x_right;  // physical nodes, x_left[0] = c, x_right[N] = c
    DDVector u_left, u_right;  // nodal values
    DDVector coef_left, coef_right;
    DD flux_left, flux_right;  // eps u'(c) from each side
    double tail_estimate = 0.0;
    // The interpolant at x in [a, b] (the left piece at x = c).
    [[nodiscard]] DD eval(const DD& x) const;
};

[[nodiscard]] TwoDomainSolution solve_two_domain(const TwoDomainProblem& p);

}  // namespace NiTCAD::verify::reference
