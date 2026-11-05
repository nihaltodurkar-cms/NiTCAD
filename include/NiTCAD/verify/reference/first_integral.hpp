// The first integral of 1D Poisson-Boltzmann toward a neutral bulk, in double-double (Unit 55, TVR-2; shared by
// TVB-07's pn junction and TVB-09's MOS capacitor, DECISIONS Unit 55 14.4).
//
// A neutral bulk of net doping C (n_b - p_b = C, n_b p_b = n_i^2) has the reduced level u_b = asinh(C / (2 n_i))
// (u = psi / V_T, n = n_i e^u, p = n_i e^-u). Toward it, (1/2) u'^2 = (q / (eps V_T)) Phi_b(u) with
//     Phi_b(u) = n_b G(t) + p_b G(-t),   t = u - u_b,   G(t) = e^t - 1 - t
// (G by its Taylor series for |t| < 1/4, else expm1(t) - t). From a surface (or junction) level u_s the profile
// at distance d into the bulk follows by inverting d(t) = int k: with t = |u - u_b| = e^s and x_hat = x / L_D,
// L_D = sqrt(eps V_T / (q |C|)),
//     d_hat(s) = int_s^{s_s} k(sigma) d sigma,   k(sigma) = e^sigma / sqrt(2 Phi_b / |C|) at t = e^sigma,
// by tanh-sinh (the integrand tends to 1 in the bulk, so no endpoint singularity) and safeguarded Newton on s
// (d d_hat / ds = -k). Phi is evaluated from t itself, never from u - u_b, which cancels when t is small.
#pragma once

#include <expected>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

// e^t - 1 - t without cancellation.
[[nodiscard]] DD poisson_boltzmann_G(const DD& t);

struct NeutralBulk {
    DD u_b;     // asinh(C / (2 n_i))
    DD n_b, p_b;
    DD doping;  // |C| [cm^-3]
};

// Preconditions: n_i > 0 and C != 0, both finite (the caller checks).
[[nodiscard]] NeutralBulk neutral_bulk(double ni_cm3, double net_doping_cm3);

// Phi_b(u) [cm^-3].
[[nodiscard]] DD first_integral_phi(const NeutralBulk& b, const DD& u);

struct ProfilePoint {
    DD u;     // reduced potential
    DD tail;  // t = |u - u_b|
};

// u at each distance d >= 0 [cm] from the surface level u_s into the bulk (any order); L_D is the bulk's Debye
// length [cm]. u_s = u_b gives u_b everywhere. Errors: a Newton or quadrature failure, or a negative distance.
[[nodiscard]] std::expected<std::vector<ProfilePoint>, std::string> bulk_profile(const NeutralBulk& b,
                                                                                 const DD& u_surface,
                                                                                 const DD& debye_cm,
                                                                                 std::span<const double> distance_cm);

}  // namespace NiTCAD::verify::reference
