// The abrupt pn junction under bias in the quasi-equilibrium (QE) model: depletion charge and capacitance (Unit 55,
// TVR-4, TVB-17; DECISIONS Unit 55 16.2 (c), 16.3). Base-free, in double-double.
//
// The QE model (16.2 (c)). Electrons are referred to the n side's quasi-Fermi level and holes to the p side's, both
// flat across the junction. With V = V_p - V_n the anode voltage and v = V / V_T, the n side's level is the
// reference: n = n_i e^u and p = n_i e^(v - u). Writing u = u' + v/2 gives n = m e^u' and p = m e^-u', with
//     m = n_i e^(v/2).
// So the QE junction at bias V is TVR-2's equilibrium junction (pn_junction.hpp, first_integral.hpp) with n_i
// replaced by m and the potential shifted by v/2 (derived here, PD). The bias enters only through m.
//
// Bulks (each exact): n_n = N_D/2 + sqrt(N_D^2/4 + m^2), p_n = m^2 / n_n, u_n = ln(n_n / m); the mirror on the p
// side. The first integral Phi_b(u) = n_b G(u - u_b) + p_b G(u_b - u), G(t) = e^t - 1 - t. The junction level u_j
// solves Phi_n(u_j) = Phi_p(u_j) (Brent, in double-double).
// The depletion charge per area on the n side is Gauss's Q = eps E_max = sqrt(2 q eps V_T Phi_n(u_j)).
// The capacitance is C = -dQ/dV (Q falls as V rises), taken exactly:
//   - one Newton step on Phi_n - Phi_p = 0 in dual numbers over double-double, from the converged u_j, with
//     dm/dv = m/2 seeded;
//   - this gives du_j/dv = -F_v / F_u to rounding, and Q's dual part is dQ/dv.
// Second method: complex step (Squire and Trapp; complex_step.hpp's rule) on the same Q(v) in complex<double>, with
// u_j(v + ih) by complex Newton from the real root. It shares the formulas but not the arithmetic or the
// differentiation.
// Reported, not judged: the depletion approximation's C_dep = eps / W, W = sqrt(2 eps (V_bi - V - 2 V_T) / q
// (1/N_A + 1/N_D)) (the -2 V_T majority-tail correction), with V_bi - V = V_T (u_n - u_p) at m. Derived here: dropping
// the junction's minority terms p_n e^(u_n - u_j) and n_p e^(u_j - u_p) (the carrier densities at x = 0) and taking
// n_n = N_D, p_p = N_A turns the first integral into exactly this form. So C_dep differs from C only by terms
// exponentially small in the barrier and by m^2 / N^2: 9e-8 relative at -0.5 V, rounding by -5 V.
#pragma once

#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/pn_junction.hpp"

namespace NiTCAD::verify::reference {

struct QeJunction {
    DD m_cm3;                    // n_i e^(V / (2 V_T))
    DD u_n, u_p, u_j;            // reduced levels in the shifted frame (u')
    DD barrier_V;                // V_bi - V = V_T (u_n - u_p)
    DD charge_C_per_cm2;         // Q = eps E_max
    DD capacitance_F_per_cm2;    // -dQ/dV
    DD depletion_form_F_per_cm2; // eps / W with the -2 V_T correction (reported)
    long long brent_evaluations = 0;
};

// Errors: an input of `in` not finite and positive, a bias not finite or so far reverse that m < 1e-140 cm^-3 (m^2
// would leave the normal range; about -15 V at 300 K for n_i ~ 1e10), or Brent not converging.
[[nodiscard]] std::expected<QeJunction, std::string> qe_junction(const PnJunctionInput& in, double bias_V);

// -dQ/dV by complex step in complex<double> (the second method). Errors as qe_junction, and Newton not converging.
[[nodiscard]] std::expected<double, std::string> qe_capacitance_complex_step(const PnJunctionInput& in,
                                                                            double bias_V);

}  // namespace NiTCAD::verify::reference
