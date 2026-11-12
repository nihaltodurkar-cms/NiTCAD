// Avalanche multiplication in an abrupt pn junction: the ionization integrals of the van Overstraeten-de Man
// coefficients over the quasi-equilibrium (QE) field, in double-double (Unit 55, TVR-5a, TVB-23; DECISIONS Unit 55
// 17.3, 17.12, 17.13 (b); IG-Q3 A). Base-free.
//
// The field. The QE junction of qe_junction.hpp (p side x < 0, n side x > 0; anode bias V on the p side; m = n_i
// e^(V / 2 V_T)), written here so that it holds at any reverse bias (IG-Q3 A): qe_junction forms m^2, which leaves
// the normal range near -15 V. Each side is described from its own bulk, with the minority density carried as a
// logarithm:
//     p side, s = u - u_p >= 0:  Phi_p = p_p G(-s) + n_p G(s),   ln n_p = 2 ln m - ln p_p,
//     n side, r = u_n - u >= 0:  Phi_n = n_n G(-r) + p_n G(r),   ln p_n = 2 ln m - ln n_n,
// G(t) = e^t - 1 - t, and the minority term taken as e^(ln n_p + s) - n_p (1 + s) once n_p itself is below e^-300
// (it is then exact to rounding or zero; nothing is dropped by hand). The barrier is B = u_n - u_p = ln n_n +
// ln p_p - 2 ln m, the junction level s_j solves Phi_p(s_j) = Phi_n(B - s_j) by Brent, and the field at a level is
// E = V_T sqrt(2 q Phi / (eps V_T)), rising monotonically from each bulk to the junction.
//
// The integrals (17.13 (b); the form as J. M. Park, TU Wien 2004, eqs (3.26)-(3.28), PB; derived here for both
// injections, PD). Electrons move toward +x and enter the depletion layer from the p side. With phi(x) =
// int (alpha_n - alpha_p) dx from the p side,
//     I_n = int alpha_n e^-phi dx,   I_p = int alpha_p e^-phi dx,   1 - I_n = e^-phi(W) - I_p  (an identity),
//     M_n = 1 / (1 - I_n)  (electrons injected from the p side),   M_p = e^-phi(W) / (1 - I_n)  (holes from the n side).
// In the level coordinate dx = V_T ds / E on the p side and -V_T dr / E on the n side, so no profile is inverted.
// alpha = gamma A exp(-gamma B / E) with the branch (low below the switch field, high from it on) of each carrier
// (TVR-3's model_generation::vodm). The integration covers the levels where E >= gamma B_min / 200, B_min the
// smallest B of the four branches: outside, every alpha is below A e^-200, and the omitted part of each integral is
// below A e^-200 times the depletion width (about 1e-80; stated, not estimated). The range is split at each
// doubling of that cut field and at each carrier's switch field, where alpha keeps its value but not its slope.
//
// Method: Dormand-Prince 5(4) (dopri5.hpp) on d(phi, I_n, I_p) / d level, piece by piece, at relative tolerance
// 1e-13. `integration_error` is the difference of I_n from the run at 1e-11 (a measured tolerance error, not a
// bound), and `identity_residual` the identity above, which the three integrated quantities do not share.
// Second method (multiplication_quadrature): on equal panels of each piece, an outer Kronrod-15 rule of
// alpha e^-phi with phi at each node by an inner Kronrod-15 rule from the panel's start. It shares the field and
// the coefficients, not the integration. It converges slowly (about a factor 200 per doubling of the panels,
// measured), which is why it is the cross-check and not the method: 16 panels agree to about 5e-10.
#pragma once

#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct IonizationBranches {
    double A_low_per_cm = 0.0;
    double B_low_V_per_cm = 0.0;
    double A_high_per_cm = 0.0;
    double B_high_V_per_cm = 0.0;
    double switch_V_per_cm = 0.0;
};

struct MultiplicationInput {
    double q_C = 0.0;
    double eps_F_per_cm = 0.0;
    double thermal_voltage_V = 0.0;
    double ni_cm3 = 0.0;
    double acceptors_cm3 = 0.0;  // x < 0
    double donors_cm3 = 0.0;     // x > 0
    IonizationBranches electron, hole;
    double gamma = 1.0;          // the temperature factor (1 at 300 K)
};

struct Multiplication {
    DD barrier_V;              // V_bi - V = V_T B
    DD bending_p_V, bending_n_V;  // V_T s_j and V_T (B - s_j)
    DD field_max_V_per_cm;
    DD depletion_p_cm, depletion_n_cm;  // eps E_max / (q N), reported
    DD ionization_integral;    // I_n
    DD hole_integral;          // I_p
    DD phi_total;              // phi(W)
    bool below_breakdown = false;  // I_n < 1
    DD multiplication_n, multiplication_p;  // 0 unless below_breakdown
    DD inverse_n, inverse_p;   // 1 / M_n = 1 - I_n and 1 / M_p = (1 - I_n) e^phi(W) (negative past breakdown)
    double integration_error = 0.0;  // |I_n(rtol 1e-13) - I_n(rtol 1e-11)|
    double identity_residual = 0.0;  // |(1 - I_n) - (e^-phi(W) - I_p)|
    long long steps = 0;             // accepted steps of the 1e-13 run
};

// Errors: an input not finite and positive (a coefficient set with A = 0 everywhere is allowed: no ionization), a
// bias not finite or forward (V > 0), Brent not converging, or the integrator not finishing.
[[nodiscard]] std::expected<Multiplication, std::string> multiplication(const MultiplicationInput& in, double bias_V);

struct MultiplicationQuadrature {
    DD ionization_integral, hole_integral, phi_total;
};

// The second method, on `panels` >= 1 equal panels per piece. Errors as multiplication.
[[nodiscard]] std::expected<MultiplicationQuadrature, std::string> multiplication_quadrature(
    const MultiplicationInput& in, double bias_V, int panels = 16);

// Where the field falls to `fraction` of its peak on each side (TVB-23's Q2 window, DECISIONS Unit 55 17.13 (b)):
// x = -int V_T ds / E from that level to the junction's on the p side, and the mirror on the n side, the junction at
// x = 0. The level by Brent, the integral by adaptive Kronrod-15 (relative 1e-25); E >= fraction E_max on the range,
// so the integrand is bounded. `abs_error` is the quadrature's summed estimate.
struct FieldWindow {
    DD field_max_V_per_cm;
    DD x_p_cm, x_n_cm;  // x_p < 0 < x_n
    double abs_error = 0.0;
};
// Errors: as multiplication; fraction not in (0, 1); Brent or the quadrature not converging.
[[nodiscard]] std::expected<FieldWindow, std::string> multiplication_field_window(const MultiplicationInput& in,
                                                                                 double bias_V, double fraction);

// The bias in [reverse_V, forward_V] (reverse_V < forward_V <= 0) at which 1 / M of the chosen carrier equals
// `inverse` (0: breakdown, where I_n = 1), by Brent to 1e-10 V. 1 / M falls as the reverse bias grows. Errors: as
// multiplication; the bracket does not hold the value.
[[nodiscard]] std::expected<double, std::string> bias_at_inverse_multiplication(const MultiplicationInput& in,
                                                                               bool electron, double inverse,
                                                                               double reverse_V, double forward_V);

}  // namespace NiTCAD::verify::reference
