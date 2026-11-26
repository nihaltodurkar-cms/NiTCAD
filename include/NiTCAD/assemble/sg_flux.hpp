// Scharfetter-Gummel edge fluxes in scaled variables, with exact partials (legacy device1d.cpp
// residual, "Jn" / "Jp" and their derivatives).
//
// For an edge from node 1 (mesh Edge::first) to node 2 (Edge::second), with delta = psi2 - psi1:
//     electrons  Jn =  a (n2 B(delta) - n1 B(-delta))
//     holes      Jp = -a (p2 B(-delta) - p1 B(delta))
// where a is the edge's scaled diffusivity times its scaled coupling area over its scaled length
// (the assembler's business; the legacy 1D form is dn_edge / h). Both are the particle-flux
// expressions the legacy continuity rows use, in its sign convention. At equilibrium
// (n = n_ie e^psi, p = n_ie e^-psi) both vanish identically, since B(-delta) = e^delta B(delta).
// Swapping the two ends negates the flux exactly.
#pragma once

#include "NiTCAD/assemble/bernoulli.hpp"

namespace NiTCAD::assemble {

struct EdgeFlux {
    double flux;
    double d_psi1;  // d flux / d psi at node 1
    double d_psi2;  // = -d_psi1
    double d_c1;    // d flux / d density at node 1
    double d_c2;
    // d flux / d delta alone, the driving term's own partial (Unit 62: density gradient's Lambda
    // enters delta only); set by SystemData::edge_fluxes with the model on, else 0.
    double d_delta = 0.0;
    // d flux / d the other carrier's density at node 1 and node 2 (Unit 66: the carrier-dependent
    // Klaassen mobility); set by SystemData::edge_fluxes with that model on, else 0.
    double d_other1 = 0.0;
    double d_other2 = 0.0;
};

[[nodiscard]] inline EdgeFlux sg_electron_flux(double a, double psi1, double psi2, double n1,
                                               double n2) noexcept {
    const double delta = psi2 - psi1;
    const double bp = bernoulli(delta);
    const double bm = bernoulli(-delta);
    const double d_psi2 =
        a * (n2 * bernoulli_derivative(delta) + n1 * bernoulli_derivative(-delta));
    return {a * (n2 * bp - n1 * bm), -d_psi2, d_psi2, -a * bm, a * bp};
}

[[nodiscard]] inline EdgeFlux sg_hole_flux(double a, double psi1, double psi2, double p1,
                                           double p2) noexcept {
    const double delta = psi2 - psi1;
    const double bp = bernoulli(delta);
    const double bm = bernoulli(-delta);
    const double d_psi2 =
        a * (p2 * bernoulli_derivative(-delta) + p1 * bernoulli_derivative(delta));
    return {-a * (p2 * bm - p1 * bp), -d_psi2, d_psi2, a * bp, -a * bm};
}

}  // namespace NiTCAD::assemble
