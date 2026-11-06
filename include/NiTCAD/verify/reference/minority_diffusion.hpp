// The linear minority-carrier diffusion problem of TVB-15 (DECISIONS Unit 55 16.2 (e), 16.3, 16.4), in
// double-double. Base-free: every input is a number the caller passes.
//
//     D u'' - u / tau = -G(x)   on [0, L],   u(0) = u(L) = 0,
//
// u = Delta n the excess minority density, D the minority diffusivity, tau the lifetime of the exact linearization of
// SRH about the equilibrium state, and G a Gaussian generation G0 exp(-(x - c)^2 / (2 s^2)) (the TVQ3 hook's G).
// Ohmic ends hold n = n0, so u = 0 there.
//
// Derived here (PD).
//   - The lifetime. R = (n p - n_ie^2) / (tau_p (n + n1) + tau_n (p + p1)) with n1 = p1 = n_ie vanishes at
//     equilibrium, so its first variation is the numerator's over the denominator D0. With quasi-neutral excess
//     (Delta p = Delta n), 1 / tau = (n0 + p0) / D0 and D0 = tau_p (n0 + n_ie) + tau_n (p0 + n_ie).
//   - The Green's function of -D u'' + u / tau with u(0) = u(L) = 0, k = 1 / sqrt(D tau):
//         g(x, s) = sinh(k x<) sinh(k (L - x>)) / (D k sinh(k L)),   x< = min(x, s), x> = max(x, s).
//     It is written without overflow as e^(-k |x - s|) (1 - e^(-2 k x<)) (1 - e^(-2 k (L - x>))) /
//     (2 D k (1 - e^(-2 k L))), with each 1 - e^(-y) as -expm1(-y).
//   - u(x) = integral_0^L g(x, s) G(s) ds. It is split at s = x (g has a kink there) and integrated by adaptive
//     Gauss-Kronrod (TC3) to a relative 0.5e-28.
// Second method (I3): Chebyshev collocation on P equal pieces of degree n. The collocation rows are
// D (2 / h_p)^2 (D^2 u)_j - u_j / tau = -G(x_j) at each piece's interior points. The patching rows are u and u'
// continuous (Boyd (22.5), as chebyshev.hpp). The ends are Dirichlet. The system is a dense LU in double-double,
// and the error estimate is the largest of each piece's last two coefficients.
#pragma once

#include <vector>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct MinorityProblem {
    double diffusivity_cm2_s = 0.0;  // D > 0
    double lifetime_s = 0.0;         // tau > 0
    double length_cm = 0.0;          // L > 0
    double peak_generation = 0.0;    // G0 [cm^-3 s^-1], finite
    double centre_cm = 0.0;          // c in [0, L]
    double width_cm = 0.0;           // s > 0
};

// 1 / tau of the linearized SRH rate (above): (n0 + p0) / (tau_p (n0 + n_ie) + tau_n (p0 + n_ie)).
// Throws std::invalid_argument unless every input is finite and positive.
[[nodiscard]] DD srh_linear_inverse_lifetime(double n0, double p0, double n_ie, double tau_n, double tau_p);

// The Green's function g(x, s) (above). Throws std::invalid_argument for a problem outside its ranges or x, s
// outside [0, L].
[[nodiscard]] DD minority_green(const MinorityProblem& p, double x, double s);

struct MinorityValue {
    DD value;               // u(x)
    double abs_error = 0;   // the quadrature's estimate (both halves)
    bool ok = false;        // both halves converged
};

// u(x) by the Green's function integral. Throws as minority_green.
[[nodiscard]] MinorityValue minority_excess(const MinorityProblem& p, double x);

struct MinorityCollocation {
    bool ok = false;
    std::vector<DD> values;     // u at the requested points
    double tail_estimate = 0;   // the largest of the pieces' last two Chebyshev coefficients
};

// The second method: u at each of xs (in [0, L]) by collocation on `pieces` pieces of degree `degree`.
// Throws std::invalid_argument for pieces < 1, degree < 4 or the ranges above.
[[nodiscard]] MinorityCollocation minority_collocation(const MinorityProblem& p, const std::vector<double>& xs,
                                                       int pieces, int degree);

}  // namespace NiTCAD::verify::reference
