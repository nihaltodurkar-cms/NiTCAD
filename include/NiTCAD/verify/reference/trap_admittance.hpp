// A MOS capacitor with one interface trap level: the operating point and the small-signal admittance of the
// single-level equivalent circuit (Unit 55, TVR-5a, TVB-22; DECISIONS Unit 55 17.3, 17.12, 17.13 (a), 17.14).
// Base-free: the operating point in double-double, the admittance in complex double.
//
// Operating point (thermal equilibrium; u = psi / V_T with n = n_i e^u, p = n_i e^-u; u_b the bulk level, t = u_s -
// u_b). The semiconductor charge is the Kingston-Neustadter charge of first_integral.hpp, Q_s = -sgn(t)
// sqrt(2 q eps V_T Phi_b(u_s)). The level, at tau = (E_t - E_i) / kT, has the Fermi occupancy f_0 = 1 / (1 +
// e^(tau - u_s)) and the sheet charge Q_it = q N_t (1 - f_0) (donor-like) or -q N_t f_0 (acceptor-like). Gauss's law
// across the oxide (derived here, PD):
//     V_G = V_FB + V_T (u_s - u_b) - (Q_s + Q_it) / C_ox,
// V_FB the flat-band voltage of the same structure without the level. The right side rises with u_s, so u_s is the
// one root, by Brent in double-double.
//
// The level's admittance (17.13 (a), derived from the rate equation df/dt = c_n n (1 - f) - c_n n_1 f - c_p p f +
// c_p p_1 (1 - f), n_1 = n_i e^tau, p_1 = n_i e^-tau; PD). With the majority carriers in equilibrium at the surface
// and the minority carriers not following,
//     Y_it = i w C_T / (1 + i w tau),  C_T = q N_t f_0 (1 - f_0) / V_T,  tau = 1 / (c_maj (m_0 + m_1)) = f_0 / (c_n n_0)
// (n-type; (1 - f_0) / (c_p p_0) on p-type), m the majority density at the surface. The same series C-R branch is
// Goetzberger's Fig. 1B (1974; PA for the circuit), whose eq (2) writes tau = 1 / (sigma v n_s0), without the
// occupancy factor: the two agree only as f_0 -> 1, and the model's own linearization is the reference.
//
// The circuit (per unit area):  1 / Y = 1 / (i w C_ox) + 1 / (i w C_s + Y_it) + Z_b,
//   - C_s = -dQ_s / dpsi_s, the Kingston-Neustadter capacitance (analytic);
//   - Z_b = R_b / (1 + i w tau_d), the neutral substrate: R_b = L_n / (q mu N), tau_d = eps / (q mu N), L_n the
//     neutral length the caller gives. The depletion edge is not sharp; the caller takes the difference between
//     L_n = L and L_n = L - W (W = |Q_s| / (q N), reported) as the reference's uncertainty.
//
// Predicates, computed here from the inputs (17.13 (a); their limits are frozen in 17.14):
//   P1 = c_min (m'_0 + m'_1) / (c_maj (m_0 + m_1)): the minority carriers' share of the level's exchange;
//   P2 = N_t c_maj / (mu E_s): the majority quasi-Fermi drop across the depletion layer that the level's current
//        causes, relative to the signal, at w tau = 1 (E_s = |Q_s| / eps);
//   P3 = the minority carriers' share of Phi_b(u_s): below it the low- and high-frequency C_s agree.
//
// Second method (trap_admittance_direct): the linearized equations for (d psi_s, d f, d psi_ox) at unit gate signal,
// solved as a 3 x 3 complex system, with the level's rate written from n_0 and n_1 (not from f_0 (1 - f_0)), so it
// shares neither the circuit's algebra nor the detailed-balance identity n_0 / (n_0 + n_1) = f_0.
#pragma once

#include <complex>
#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct TrapAdmittanceInput {
    double q_C = 0.0;
    double thermal_voltage_V = 0.0;
    double eps_F_per_cm = 0.0;      // the semiconductor's permittivity
    double ni_cm3 = 0.0;
    double net_doping_cm3 = 0.0;    // N_D - N_A of the substrate, not zero
    double oxide_capacitance_F_per_cm2 = 0.0;
    double flat_band_V = 0.0;       // of the structure without the level
    double gate_V = 0.0;
    bool donor = false;             // donor-like (else acceptor-like)
    double density_cm2 = 0.0;       // N_t >= 0
    double energy_eV = 0.0;         // E_t - E_i
    double capture_n_cm3_s = 0.0;   // c_n = sigma_n v_th,n
    double capture_p_cm3_s = 0.0;   // c_p = sigma_p v_th,p
    double mobility_cm2_Vs = 0.0;   // of the substrate's majority carriers
};

struct TrapOperatingPoint {
    DD u_s, u_b;                    // reduced surface and bulk levels
    double surface_potential_V = 0.0;  // V_T (u_s - u_b)
    double occupancy = 0.0;         // f_0
    double semiconductor_charge_C_per_cm2 = 0.0;  // Q_s
    double trap_charge_C_per_cm2 = 0.0;           // Q_it
    double semiconductor_capacitance_F_per_cm2 = 0.0;  // C_s
    double trap_capacitance_F_per_cm2 = 0.0;           // C_T
    double time_constant_s = 0.0;   // tau
    double surface_field_V_per_cm = 0.0;
    double depletion_width_cm = 0.0;  // |Q_s| / (q N)
    double majority_surface_cm3 = 0.0, majority_level_cm3 = 0.0;  // m_0 and m_1
    double minority_surface_cm3 = 0.0, minority_level_cm3 = 0.0;
    double p1_minority_exchange = 0.0, p2_fermi_drop = 0.0, p3_minority_charge = 0.0;
    double sheet_resistance_Ohm_cm = 0.0;  // 1 / (q mu N): R_b per unit neutral length
    double relaxation_time_s = 0.0;        // tau_d
    long long brent_evaluations = 0;
};

// Errors: an input not finite, or not positive where it must be (N_t may be 0); a level whose n_1 or p_1 leaves the
// normal range; Brent not converging.
[[nodiscard]] std::expected<TrapOperatingPoint, std::string> trap_operating_point(const TrapAdmittanceInput& in);

// The circuit's admittance [S/cm^2] at angular frequency omega > 0, with neutral length L_n >= 0 [cm].
[[nodiscard]] std::complex<double> trap_admittance(const TrapAdmittanceInput& in, const TrapOperatingPoint& op,
                                                   double omega, double neutral_length_cm);

// The second method. An empty result if the 3 x 3 system is singular.
[[nodiscard]] std::expected<std::complex<double>, std::string> trap_admittance_direct(
    const TrapAdmittanceInput& in, const TrapOperatingPoint& op, double omega, double neutral_length_cm);

}  // namespace NiTCAD::verify::reference
