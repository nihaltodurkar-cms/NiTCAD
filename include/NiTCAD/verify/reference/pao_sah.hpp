// The long-channel MOSFET under the gradual-channel approximation: Pao and Sah's double integral, Brews' charge sheet
// (a diagnostic) and the subthreshold swing, for TVB-20 and TVB-21 (Unit 55 TVR-4; DECISIONS 16.3, 16.17, 16.18 (d),
// (e), 16.46). Base-free: every constant and material value is an input (R-7).
//
// Pao and Sah (Solid-State Electronics 9, 927 (1966); PB through two reproductions, 16.17): on a uniform p-substrate,
// Boltzmann statistics, the electrons at the channel's quasi-Fermi level V across each slice, with u = psi / V_T and
// v = V / V_T,
//     E^2(psi, V) = (2 q V_T / eps) [p_b (e^-u + u - 1) + n_b B(u, v)],
//     form 1 (Wong's reproduction):            B = e^-v (e^u - 1) - u,
//     form 2 (Galup-Montoro's eq 5, chi = 0):  B = e^-v (e^u - u - 1);
// the surface potential psi_s(V) solves V_G - V_FB = psi_s + eps E(psi_s, V) / C_ox (Brent in double-double; the
// field's root taken as 0 where form 1's E^2 is negative, psi below about 1e-15 V);
//     Q_I(V) = -q int_{c psi_B}^{psi_s(V)} n_b e^{u - v} / E(psi, V) dpsi,   psi_B = V_T ln(p_b / n_i),
// c the cutoff fraction (1: Wong's lower limit psi_B; 16.18 (d) also takes c = 1/2), and
//     Q = I_D L / W = mu int_{V_S}^{V_D} (-Q_I(V)) dV        [A, per unit W / L].
// At c psi_B > 0 the inner integrand is finite (E > 0 there): no singularity at the inner lower limit (16.18 (n)).
// Both integrals by adaptive Gauss-Kronrod in double-double.
// Second method (16.18 (d)): g_d = dQ/dV_D = mu (-Q_I(V_D)) (Galup-Montoro's eq 4), against a fourth-order
// difference of Q in V_D in double-double (a dyadic step, so that V_D +- h are exact).
// Brews' charge sheet (Solid-State Electronics 21, 345 (1978); PC, used as PD, 16.18 (e)): with psi_s at the source
// and drain from the same gate equation (form 2), Q_I(psi) = -C_ox (V_G - V_FB - psi) + sqrt(2 q eps N_A (psi - V_T)),
//     Q_cs = mu [int_{psi_0}^{psi_L} (-Q_I) dpsi + V_T (Q_I(psi_L) - Q_I(psi_0))],
// in closed form, with Gauss-Kronrod on the first integral as its second method. A diagnostic, never gated.
#pragma once

#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct PaoSahInput {
    double q_C = 0.0;
    double eps_F_per_cm = 0.0;  // the semiconductor's
    double thermal_voltage_V = 0.0;
    double ni_cm3 = 0.0;
    double acceptors_cm3 = 0.0;  // the uniform substrate
    double oxide_capacitance_F_per_cm2 = 0.0;
    double flat_band_V = 0.0;
    double mobility_cm2_Vs = 0.0;  // constant
    int form = 1;                  // 1 or 2 (above)
    double cutoff_fraction = 1.0;  // c in (0, 1]
};

struct PaoSahPoint {
    DD current;              // Q = I_D L / W [A] at (V_G, V_D), V_S = 0
    DD surface_source_V;     // psi_s(0)
    DD surface_drain_V;      // psi_s(V_D)
    DD bulk_potential_V;     // psi_B
    double quadrature_error = 0.0;  // the outer integral's estimate plus mu times the worst inner estimate, relative
};

// Errors: an input not finite or out of range, V_G <= V_FB, V_D < 0; Brent or a quadrature failing.
[[nodiscard]] std::expected<PaoSahPoint, std::string> pao_sah(const PaoSahInput& in, double gate_V, double drain_V);

// Q_I(V) [C/cm^2] and psi_s(V) at one channel potential.
struct InversionCharge {
    DD charge;
    DD surface_V;
    double quadrature_error = 0.0;  // relative
};
[[nodiscard]] std::expected<InversionCharge, std::string> pao_sah_inversion_charge(const PaoSahInput& in,
                                                                                   double gate_V, double channel_V);

// The second method: |dQ/dV_D (fourth-order central difference, step h) - mu (-Q_I(V_D))| / |mu Q_I(V_D)|.
[[nodiscard]] std::expected<double, std::string> pao_sah_conductance_identity(const PaoSahInput& in, double gate_V,
                                                                              double drain_V, double h_V);

// Brews' charge sheet (a diagnostic): Q_cs [A] and its second method's relative difference.
struct ChargeSheet {
    DD current;
    double second_method = 0.0;
};
[[nodiscard]] std::expected<ChargeSheet, std::string> charge_sheet(const PaoSahInput& in, double gate_V,
                                                                   double drain_V);

}  // namespace NiTCAD::verify::reference
