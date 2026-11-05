// The equilibrium abrupt pn junction, Boltzmann statistics, constant permittivity: the exact solution of 1D
// Poisson-Boltzmann by its first integral (Unit 55, TVR-2, TVB-07; DECISIONS Unit 55 6, 14.4), with an independent
// cross-check by Chebyshev collocation (I3).
//
// The one-sided first integral and its profile are first_integral.hpp's (shared with TVB-09).
//
// Problem. p side x < 0 (acceptors N_A), n side x > 0 (donors N_D), the doping step at x = 0; reduced potential
// u = psi / V_T with n = n_i e^u and p = n_i e^-u, and
//     u'' = (q / (eps V_T)) (n - p - N_D)  (x > 0),     u'' = (q / (eps V_T)) (n - p + N_A)  (x < 0),
// u -> u_n and u -> u_p in the neutral bulks, u_n = asinh(N_D / (2 n_i)), u_p = -asinh(N_A / (2 n_i)) (the exact
// neutral levels; V_T ln(N_A N_D / n_i^2) is their majority-carrier limit, which differs by about (n_i / N)^2
// relative and so cannot serve as a reference to rounding).
//
// First integral (each side, b its bulk, t = u - u_b): (1/2) u'^2 = (q / (eps V_T)) Phi_b(u),
//     Phi_b(u) = n_b G(t) + p_b G(-t),   G(t) = e^t - 1 - t >= 0,
// with n_b, p_b the bulk's carrier densities (n_b - p_b its net doping). The junction potential u_j solves
// Phi_n(u_j) = Phi_p(u_j) on (u_p, u_n) (u and eps u' continuous at x = 0), by Brent's procedure; the field at the
// junction is E_max = V_T sqrt(2 q Phi(u_j) / (eps V_T)). The profile: on the n side, with t = u_n - u = e^s,
//     x(s) = int_s^{s_j} k(sigma) d sigma,   k(sigma) = e^sigma / u'(u_n - e^sigma),
// whose integrand tends to the bulk Debye length as sigma -> -infinity (no endpoint singularity), by tanh-sinh;
// psi(x) inverts it by safeguarded Newton on s (dx/ds = -k). The p side is the mirror, t = u - u_p.
// All of it in double-double; G by its Taylor series for |t| < 1/4 (no cancellation), else expm1(t) - t.
//
// The infinite-domain solution is the reference. A device of finite length with Dirichlet contacts at the bulk
// levels differs from it by about the tail V_T t at the contacts (the linearised problem's boundary mismatch,
// decaying into the device); pn_potential reports it per point (t), and the collocation solution, which solves the
// finite-domain problem, measures the difference directly.
//
// Cross-check (I3): pn_collocation solves the finite-domain problem u(x_left) = u_p, u(x_right) = u_n by Chebyshev
// collocation (chebyshev.hpp's cheb.m matrices, in double) on `subdomains` equal pieces per side, u and u'
// continuous at every piece boundary (Boyd's patching), by damped Newton with a dense LU in double. It shares no
// code with the first integral. Its error estimate is the last two Chebyshev coefficients per piece.
//
// Inputs are physical data, passed in by the caller (the benchmark takes them from NiTCAD's material set: I1).
#pragma once

#include <expected>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct PnJunctionInput {
    double q_C = 0.0;           // elementary charge
    double eps_F_per_cm = 0.0;  // permittivity
    double thermal_voltage_V = 0.0;
    double ni_cm3 = 0.0;         // intrinsic density
    double acceptors_cm3 = 0.0;  // x < 0
    double donors_cm3 = 0.0;     // x > 0
};

struct PnJunction {
    PnJunctionInput in;
    DD u_p, u_n;           // reduced bulk levels
    DD n_p, p_p, n_n, p_n; // bulk densities [cm^-3]
    DD psi_p, psi_n;       // V_T u_p, V_T u_n [V]
    DD built_in_V;         // psi_n - psi_p
    DD u_j, psi_j;         // the junction's
    DD field_max_V_per_cm; // |psi'(0)|
    DD depletion_p_cm, depletion_n_cm;  // eps E_max / (q N): the charge-balance depletion edges (reported)
    DD debye_p_cm, debye_n_cm;          // sqrt(eps V_T / (q N)) of each bulk
    long long brent_evaluations = 0;
};

// Errors: an input not finite and positive; Brent not converging.
[[nodiscard]] std::expected<PnJunction, std::string> pn_junction(const PnJunctionInput& in);

struct PnPoint {
    DD psi_V;
    DD tail;  // t = |u - u_bulk| at the point (reduced)
};

// psi at each x (any order, any sign; x = 0 is the junction). Errors: a Newton or quadrature failure.
[[nodiscard]] std::expected<std::vector<PnPoint>, std::string> pn_potential(const PnJunction& j,
                                                                            std::span<const double> x_cm);

struct PnCollocation {
    std::vector<double> psi_V;   // at the query points
    double tail_estimate_V = 0;  // max over pieces of |a_{N-1}| + |a_N|, times V_T
    int newton_iterations = 0;
    double last_update = 0.0;    // max |du| of the last Newton step
};

// The finite-domain problem on [x_left, x_right] (x_left < 0 < x_right), `subdomains` pieces per side of
// Chebyshev degree `degree`; psi at the query points. Errors: Newton not converging in 50 iterations, a singular
// matrix, bad arguments.
[[nodiscard]] std::expected<PnCollocation, std::string> pn_collocation(const PnJunction& j, double x_left_cm,
                                                                       double x_right_cm, int subdomains, int degree,
                                                                       std::span<const double> query_cm);

}  // namespace NiTCAD::verify::reference
