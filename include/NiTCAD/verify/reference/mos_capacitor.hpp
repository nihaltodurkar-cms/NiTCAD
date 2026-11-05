// The ideal MOS capacitor on a uniform substrate, Boltzmann statistics, at thermal equilibrium (the quasi-static
// state): the exact 1D charge-voltage relation (Unit 55, TVR-2, TVB-09; DECISIONS Unit 55 6, 14.4).
//
// Semiconductor charge (Kingston and Neustadter, J. Appl. Phys. 26, 718 (1955); the F-function of Sze and Ng,
// Physics of Semiconductor Devices, 3rd ed., 2007, 4.2): with u = psi / V_T, the bulk level u_b and the first
// integral of first_integral.hpp,
//     Q_s(u_s) = -sgn(u_s - u_b) sqrt(2 q eps V_T Phi_b(u_s))      [C/cm^2]
// (Phi_b in cm^-3: q Phi_b V_T has the units of an energy density). The gate bias holding u_s (Gauss across the
// lumped oxide; Q_G + q Q_f + Q_s = 0):
//     V_G = V_FB + V_T (u_s - u_b) - Q_s / C_ox,
// monotone in u_s, so u_s(V_G) by Brent in double-double on [u_b - a, u_b + a], a = |V_G - V_FB| / V_T + 1 (since
// |V_T (u_s - u_b)| <= |V_G - V_FB|). The low-frequency (quasi-static) capacitance is dQ_G / dV_G = C_s C_ox /
// (C_s + C_ox), C_s = -dQ_s / dpsi_s, taken analytically,
//     dPhi_b / du = n_b expm1(t) - p_b expm1(-t),   dQ_s / dpsi_s = -sgn(t) sqrt(2 q eps V_T) Phi_b' / (2 V_T sqrt Phi_b),
// with C_s = sqrt(q eps (n_b + p_b) / V_T) at t = 0 (its limit), and by complex step on a complex-double Q_s(u) as the
// second method (TVB-09's "C_LF by complex step"), reported for |t| >= 1/2, where the complex Phi does not cancel.
//
// Flat band (the gate's work function phi_m against the semiconductor's: Sze and Ng 4.2):
//     V_FB = phi_m - chi - (E_g / 2 + (V_T / 2) ln(N_c / N_v)) + V_T u_b - q Q_f / C_ox,
// E_c - E_i in the brackets (Boltzmann), V_T u_b = E_i - E_F at the bulk, Q_f the fixed oxide charge [cm^-2].
#pragma once

#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"

namespace NiTCAD::verify::reference {

struct FlatBandInput {
    double work_function_eV = 0.0;
    double electron_affinity_eV = 0.0;
    double band_gap_eV = 0.0;
    double conduction_dos_cm3 = 0.0;
    double valence_dos_cm3 = 0.0;
    double thermal_voltage_V = 0.0;
    double ni_cm3 = 0.0;
    double net_doping_cm3 = 0.0;  // N_D - N_A of the substrate, nonzero
    double fixed_charge_cm2 = 0.0;
    double q_C = 0.0;
    double oxide_capacitance_F_per_cm2 = 0.0;
};

// Errors: an input not finite, or not positive where it must be.
[[nodiscard]] std::expected<DD, std::string> flat_band_voltage(const FlatBandInput& in);

struct MosInput {
    double q_C = 0.0;
    double eps_F_per_cm = 0.0;  // the semiconductor's
    double thermal_voltage_V = 0.0;
    double ni_cm3 = 0.0;
    double net_doping_cm3 = 0.0;  // nonzero
    double oxide_capacitance_F_per_cm2 = 0.0;
    double flat_band_V = 0.0;
};

struct MosPoint {
    DD gate_V;
    DD u_s;                 // reduced surface potential
    DD band_bending_V;      // V_T (u_s - u_b)
    DD charge_C_per_cm2;    // Q_s
    DD surface_capacitance; // C_s [F/cm^2]
    DD capacitance;         // C_LF [F/cm^2]
    double capacitance_complex_step = 0.0;  // C_LF with C_s by complex step; 0 when |t| < 1/2
    bool complex_step_evaluated = false;
    long long brent_evaluations = 0;
};

// Errors: bad inputs; Brent not converging.
[[nodiscard]] std::expected<MosPoint, std::string> mos_point(const MosInput& in, double gate_V);

[[nodiscard]] NeutralBulk mos_bulk(const MosInput& in);
[[nodiscard]] DD mos_debye_cm(const MosInput& in);

}  // namespace NiTCAD::verify::reference
