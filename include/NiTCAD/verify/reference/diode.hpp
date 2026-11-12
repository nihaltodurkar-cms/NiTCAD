// The ideal-diode current of an abrupt pn junction: Shockley's law with the quasi-equilibrium (QE) depletion edges,
// and the Sah-Noyce-Shockley (SNS) depletion-region recombination current (Unit 55, TVR-4, TVB-16; DECISIONS Unit 55
// 16.2 (d), 16.3, 16.18 (f)). Base-free, in double-double.
//
// The junction at bias V (anode p side, V = V_p - V_n, v = V / V_T) is TVR-2's equilibrium junction at
// m = n_i e^(v/2) (qe_junction.hpp): n = m e^u, p = m e^-u in the shifted frame, np = m^2 = n_i^2 e^v everywhere in
// the depletion region. Its charge-balance edges x_p = Q / (q N_A), x_n = Q / (q N_D) (Q = eps E_max) bound the
// depletion region, and the neutral bases are W'_p = L_p - x_p, W'_n = L_n - x_n (contact to edge).
//
// Diffusion current (Shockley, BSTJ 28 (1949) eqs 4.23-4.24, PA; the finite base with an ohmic contact derived here,
// PD). In each neutral base the minority excess obeys D u'' - u / tau = 0 with u = n_0 (e^v - 1) at the edge (the QE
// law of the junction) and u = 0 at the ohmic contact. So
//     J_diff = q (e^v - 1) [ D_n n_p0 / L_n coth(W'_p / L_n) + D_p p_n0 / L_p coth(W'_n / L_p) ],
// with L = sqrt(D tau) and 1 / tau SRH's exact linearization at the bulk's equilibrium (minority_diffusion.hpp). The
// long base is coth -> 1 (W' >> L), the short base L coth(W'/L) -> W' (W' << L). Without recombination (tau -> inf)
// the minority current is constant across the base, J_diff = q (e^v - 1) [D_n n_p0 / W'_p + D_p p_n0 / W'_n],
// exactly the short-base law. n_p0, p_n0 are the exact equilibrium bulk minorities (pn_junction at n_i).
//
// SNS current (Sah, Noyce and Shockley, Proc. IRE 45 (1957); used as the definition, PD, DT-Q18): with SRH on,
// J_SNS = q int_{-x_p}^{x_n} U dx, U = (np - n_i^2) / (tau_p (n + n_i) + tau_n (p + n_i)) (NiTCAD's mid-gap SRH
// form), each side with its own lifetimes, over the reference's QE potential (never the solver's: C-2). In u-space,
// dx = du / |u'(u)| with (1/2) u'^2 = (q / (eps V_T)) Phi_b(u), so on each side
//     int U dx = int_{u_edge}^{u_j} U(u) / sqrt(2 q Phi_b(u) / (eps V_T)) du,
// whose integrand is finite on the closed interval (Phi_b > 0 off the bulk). Gauss-Kronrod (qag15) in u; the second
// method is tanh-sinh on the same integrals (16.3), and their difference is reported. The change of variable itself is
// checked at L0 against U(x) integrated in x with u(x) from pn_potential.
//
// J = J_diff + J_SNS (positive forward). Low injection, Boltzmann statistics, constant D and tau per side, abrupt
// doping, and the edge definition are the approximations; their first-order sizes are measured on the solved state by
// the benchmark (16.18 (f)), not here.
#pragma once

#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/pn_junction.hpp"

namespace NiTCAD::verify::reference {

struct DiodeSide {
    double length_cm = 0.0;               // junction to ohmic contact
    double minority_diffusivity_cm2_s = 0.0;  // electrons on the p side, holes on the n side
    double tau_n_s = 0.0, tau_p_s = 0.0;  // SRH lifetimes in this side (read only with srh)
};

struct DiodeInput {
    PnJunctionInput junction;  // q, eps, V_T, n_i, N_A, N_D
    DiodeSide p_side, n_side;
    bool srh = false;
};

struct DiodeCurrent {
    DD m_cm3;                          // n_i e^(v/2)
    DD edge_p_cm, edge_n_cm;           // x_p, x_n (charge balance at V)
    DD base_p_cm, base_n_cm;           // W'_p, W'_n
    DD debye_p_cm, debye_n_cm;         // sqrt(eps V_T / (q N)) of each bulk
    DD n_p0_cm3, p_n0_cm3;             // equilibrium bulk minorities
    DD diffusion_length_p_cm, diffusion_length_n_cm;  // sqrt(D tau) (0 without srh)
    DD diffusion_A_per_cm2;            // J_diff
    DD sns_A_per_cm2;                  // J_SNS (0 without srh)
    DD total_A_per_cm2;                // J_diff + J_SNS
    double sns_quadrature_error = 0.0;     // Gauss-Kronrod's own estimate, absolute [A/cm^2]
    double sns_second_method_difference = 0.0;  // |GK in u - tanh-sinh in x| [A/cm^2]
};

// Errors: an input not finite and positive (the lifetimes only with srh); a base W' <= 0 (the depletion edge at or
// beyond a contact); a bias outside qe_junction's range; a Brent, Newton or quadrature failure.
[[nodiscard]] std::expected<DiodeCurrent, std::string> diode_current(const DiodeInput& in, double bias_V);

}  // namespace NiTCAD::verify::reference
