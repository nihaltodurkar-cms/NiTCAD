// The ideal MOS capacitor on a uniform p-type substrate away from equilibrium in two limits (Unit 55, TVR-4, TVB-10;
// DECISIONS Unit 55 16.2 (g), 16.18 (b), DT-Q17): the high-frequency (HF) capacitance in the redistribution limit,
// and the deep-depletion charge with no minority carriers. Boltzmann statistics, lumped oxide, as mos_capacitor.hpp.
//
// The redistribution limit (Brews, J. Appl. Phys. 45, 1276 (1974), premise from the abstract, PC; the model derived
// here, PD). The minority electrons share one quasi-Fermi level, shifted by V = V_T v from the holes' (which stay at
// the bulk's), uniformly across the layer, n = n_b e^{t - v}; t = u - u_b as first_integral.hpp. The first integral
// of Poisson toward the bulk is then
//     Phi(t, v) = p_b G(-t) + n_b e^{-v} G(t),
// (DECISIONS 16.17's reproduction (2) with chi = 0: the electrons' bulk density n_b e^{-v} is neutralized with them,
// so Phi has a double zero at t = 0 for every v; reproduction (1) does not, and is not used). With k = sqrt(2 q eps V_T)
// and g(t, v) = |dt/dx| = sqrt(2 q Phi / (eps V_T)):
//     Q_s(t, v) = -sgn(t) k sqrt(Phi),
//     Q_n(t, v) = -q sgn(t) int_0^t n_b e^{-v} expm1(s) / g(s, v) ds          (the electrons' excess charge)
// and the partials
//     Phi_t = -p_b expm1(-t) + n_b e^{-v} expm1(t),   Phi_v = -n_b e^{-v} G(t),
//     dQ_s/dt = -sgn k Phi_t / (2 sqrt Phi),   dQ_s/dv = -sgn k Phi_v / (2 sqrt Phi),
//     dQ_n/dt = -q sgn n_b e^{-v} expm1(t) / g(t, v),
//     dQ_n/dv = q sgn int_0^t n_b e^{-v} expm1(s) / g (1 + Phi_v / (2 Phi)) ds
// (the integrands are finite at s = 0: Phi ~ (p_b + n_b e^{-v}) s^2 / 2 there). HF keeps Q_n at its DC value (v = 0
// at DC: the operating point is mos_point's), so dv/dt = -(dQ_n/dt) / (dQ_n/dv) and
//     C_s,HF = -(1/V_T) (dQ_s/dt + dQ_s/dv dv/dt),   C_HF = C_s,HF C_ox / (C_s,HF + C_ox).
// With v held at 0 instead it is mos_point's C_LF (the electrons following the substrate). In depletion Q_n -> 0 and
// the v-term cancels the electrons' part of dQ_s/dt: C_s,HF -> the holes' depletion capacitance.
// Second method (16.18 (b)): Q_s along the constrained path, Q_n(t, v(t)) = Q_n(t_s, 0) held by Brent on v at
// t_s +- delta, delta / 2, the central differences Richardson-extrapolated (O(delta^4)); delta = 1e-3.
//
// Deep depletion (16.2 (g)): the state of a gate ramp faster than generation, with no minority carriers. The holes
// are at the bulk's level, p_b = N_A, and Phi_dd(t) = N_A G(-t), Q_dd = -sgn(t) k sqrt(Phi_dd), V_G = V_FB + V_T t -
// Q_dd / C_ox, t by Brent. Exact for that state. Second method: Q_dd as the quadrature of the charge density
// q N_A expm1(-t) over the profile, sgn(t) int_0^t q N_A expm1(-s) / g_dd(s) ds (Gauss-Kronrod in double-double).
#pragma once

#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/mos_capacitor.hpp"

namespace NiTCAD::verify::reference {

// Q_s(t, v) and Q_n(t, v) of the redistribution model, with the quadrature's estimated absolute error.
struct MosCharges {
    DD surface_C_per_cm2;
    DD minority_C_per_cm2;
    double minority_abs_error = 0.0;
};

// Preconditions as mos_point, and net_doping < 0 (p-type; else an error).
[[nodiscard]] std::expected<MosCharges, std::string> mos_redistribution_charges(const MosInput& in, const DD& t,
                                                                              const DD& v);

struct MosHfPoint {
    MosPoint dc;                      // the operating point (v = 0) and C_LF
    DD minority_charge_C_per_cm2;     // Q_n at DC
    DD dv_dt;                         // the quasi-Fermi shift per unit surface potential, both reduced
    DD surface_capacitance_hf;        // C_s,HF [F/cm^2]
    DD capacitance_hf;                // C_HF [F/cm^2]
    double quadrature_error = 0.0;    // relative, of dQ_n/dv
    double capacitance_hf_second_method = 0.0;
    double second_method_difference = 0.0;  // |C_HF - second method| [F/cm^2]
    bool second_method_evaluated = false;
};

// Errors: bad inputs, n-type substrate, t_s = 0 exactly (flat band: the sign is undefined), Brent or quadrature
// failure. The second method runs only when asked (it costs about 100 quadratures).
[[nodiscard]] std::expected<MosHfPoint, std::string> mos_hf_point(const MosInput& in, double gate_V,
                                                                  bool second_method = false);

struct DeepDepletionPoint {
    DD gate_V;
    DD t_s;                   // reduced band bending
    DD charge_C_per_cm2;      // Q_dd
    double charge_second_method = 0.0;
    double second_method_difference = 0.0;  // |Q_dd - quadrature| [C/cm^2]
    long long brent_evaluations = 0;
};

// Errors: as mos_hf_point.
[[nodiscard]] std::expected<DeepDepletionPoint, std::string> mos_deep_depletion_point(const MosInput& in,
                                                                                      double gate_V);

}  // namespace NiTCAD::verify::reference
