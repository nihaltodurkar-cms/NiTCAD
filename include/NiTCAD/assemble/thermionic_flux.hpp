// Thermionic-emission fluxes across an abrupt heterointerface, in scaled variables with exact
// partials (Unit 15; legacy device1d.cpp, the M33-S2 te_edge branch of residual_jacobian). With
// PhysicsModels::thermionic_emission they replace the Scharfetter-Gummel fluxes of sg_flux.hpp on
// every edge between two materials, in the same sign convention and EdgeFlux layout, for an edge
// from node 1 to node 2 with driving term delta (the edge's delta_n or delta_p):
//     electrons  Jn =  K (n2 g2 - n1 g1),  g1 = min(1, e^u),  g2 = r min(1, e^-u),
//                      u = delta - ln(N2 / N1),  r = N1 / N2   (N = Nc);
//     holes      Jp = -K (p2 h2 - p1 h1),  h1 = min(1, e^w),  h2 = r min(1, e^-w),
//                      w = -delta - ln(N2 / N1),  r = N1 / N2   (N = Nv).
// u is the conduction-band edge's fall from node 1 to node 2 in units of kT (w the valence band's
// rise): a carrier crossing up a step is cut by e^-step and one crossing down is not, and each
// side emits at the emission velocity (physics/thermionic_emission.hpp). K is the harmonic mean of
// the two sides' velocities times the scaled interface area (the legacy K = hmean(v) L_D / D0 in
// 1D). At equilibrium n2 / n1 = e^delta, so n2 g2 = n1 g1 and the flux vanishes exactly, for either
// sign of the step. The partials in delta take the flat side's slope, 0, at u = 0 or w = 0, where
// min(1, e^x) has its kink (legacy).
#pragma once

#include <cmath>

#include "NiTCAD/assemble/sg_flux.hpp"

namespace NiTCAD::assemble {

// log_ratio = ln(N2 / N1), ratio = N1 / N2.
[[nodiscard]] inline EdgeFlux thermionic_electron_flux(double K, double delta, double log_ratio,
                                                       double ratio, double n1,
                                                       double n2) noexcept {
    const double u = delta - log_ratio;
    const double g1 = u < 0.0 ? std::exp(u) : 1.0;
    const double g2 = ratio * (u > 0.0 ? std::exp(-u) : 1.0);
    const double dg1 = u < 0.0 ? g1 : 0.0;   // d g1 / d delta
    const double dg2 = u > 0.0 ? -g2 : 0.0;  // d g2 / d delta
    const double d_delta = K * (n2 * dg2 - n1 * dg1);
    return {K * (n2 * g2 - n1 * g1), -d_delta, d_delta, -K * g1, K * g2};
}

[[nodiscard]] inline EdgeFlux thermionic_hole_flux(double K, double delta, double log_ratio,
                                                   double ratio, double p1, double p2) noexcept {
    const double w = -delta - log_ratio;
    const double h1 = w < 0.0 ? std::exp(w) : 1.0;
    const double h2 = ratio * (w > 0.0 ? std::exp(-w) : 1.0);
    const double dh1 = w < 0.0 ? -h1 : 0.0;  // d h1 / d delta (dw / d delta = -1)
    const double dh2 = w > 0.0 ? h2 : 0.0;   // d h2 / d delta
    const double d_delta = -K * (p2 * dh2 - p1 * dh1);
    return {-K * (p2 * h2 - p1 * h1), -d_delta, d_delta, K * h1, -K * h2};
}

}  // namespace NiTCAD::assemble
