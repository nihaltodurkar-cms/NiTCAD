// The small-signal admittance of the ideal diode at forward bias: Shockley's eq 4.23 with the finite-base (coth) form,
// the depletion capacitance, and the depletion edges' motion (Unit 55, TVR-4, TVB-18; DECISIONS Unit 55 16.2 (l),
// 16.3, 16.18 (g)). Base-free, in double-double.
//
// The operating point is diode.hpp's: the QE junction at m = n_i e^(v/2), the edges x_p, x_n, the bases
// W'_p = L_p - x_p, W'_n = L_n - x_n, the minority excess U0 = n_0 (e^v - 1) at each edge and 0 at the contact.
// A perturbation V + dV e^(i w t) obeys, in each base (y measured from the contact, the edge at y = W'),
//     D du'' = (1/tau + i w) du,  du(0) = 0,
// with 1/tau SRH's exact linearization at the bulk (0 without srh). Write kappa = sqrt((1/tau + i w) / D) (principal
// root, Re kappa > 0) and L_eff = 1 / kappa.
//
// Fixed edges (Shockley, BSTJ 28 (1949) eq 4.23, PA; the finite base derived here, PD): du(W') = n_0 e^v dV / V_T,
// so the minority current entering the base is
//     Y_fixed = q D n_0 e^v / V_T * kappa coth(kappa W'),
// which is G_0 sqrt(1 + i w tau) for W' >> |L_eff| and the short-base form for W' << |L_eff|.
//
// Edge motion (derived here, PD; not in eq 4.23): the edge moves with V, dW'/dV = C_j / (q N) on each side (C_j =
// -dQ/dV from qe_junction, the edge x = Q / (q N)). Linearizing the boundary condition at the moving edge,
// du(W') + u0'(W') dW' = n_0 e^v dV / V_T, and the current there, q D (du'(W') + u0''(W') dW'), gives
//     Y_edge = q D U0 dW'/dV [ u0''/U0 - kappa coth(kappa W') u0'/U0 ],
// with u0'/U0 = coth(W'/L)/L and u0''/U0 = 1/L^2 (srh), or 1/W' and 0 (without). At w = 0, Y_fixed + Y_edge is
// d J_diff / dV of diode.hpp exactly (its W(V) edges); the L0 test checks this against a double-double difference.
//
// Capacitive: i w C_j, the QE junction's capacitance at V.
//
// Total: Y = sum over the sides (Y_fixed + Y_edge) + i w C_j. Shockley's form is reported without Y_edge. Not here:
// the depletion region's recombination (SNS), the neutral bases' ohmic impedance, high injection; their sizes are
// the benchmark's (16.18 (g)).
//
// Second method for Y_fixed: the Riccati form of the same ODE, w(y) = du / du', w' = 1 - kappa^2 w^2, w(0) = 0,
// integrated from the contact to the edge by Dormand-Prince in double-double (stable in that direction: the root
// w = 1/kappa attracts), so kappa coth(kappa W') = 1 / w(W'). Integrated in y / |L_eff|, starting at most 40 |L_eff|
// before the edge from w = 1/kappa when the base is longer (the contact's influence there is e^-80). The
// difference is reported.
#pragma once

#include <complex>
#include <expected>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/diode.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct AdmittanceSide {
    std::complex<double> fixed_edge;          // Y_fixed [S/cm^2]
    std::complex<double> fixed_edge_riccati;  // Y_fixed by the second method
    std::complex<double> edge_motion;         // Y_edge
    std::complex<double> effective_length_cm; // L_eff = 1 / kappa
};

struct AdmittancePoint {
    double omega_rad_s = 0.0;
    AdmittanceSide p, n;
    std::complex<double> capacitive;  // i w C_j
    std::complex<double> shockley;    // sum Y_fixed + i w C_j (eq 4.23 at the edges of V)
    std::complex<double> total;       // shockley + sum Y_edge
    double second_method_difference = 0.0;  // |sum Y_fixed - sum Y_fixed,riccati| [S/cm^2]
};

struct DiodeAdmittance {
    DiodeCurrent dc;                     // diode_current at the bias
    DD capacitance_F_per_cm2;            // C_j
    double base_rate_p_cm_per_V = 0.0;   // dW'_p / dV
    double base_rate_n_cm_per_V = 0.0;   // dW'_n / dV
    double lifetime_p_s = 0.0;           // the minority lifetime of each base (0 without srh)
    double lifetime_n_s = 0.0;
    std::vector<AdmittancePoint> points;  // one per frequency, in order
};

// Errors: those of diode_current and qe_junction; a frequency not finite or negative; w = 0 without srh is allowed
// (the short-base limit); the second method not converging.
[[nodiscard]] std::expected<DiodeAdmittance, std::string> diode_admittance(const DiodeInput& in, double bias_V,
                                                                          std::span<const double> omegas_rad_s);

}  // namespace NiTCAD::verify::reference
