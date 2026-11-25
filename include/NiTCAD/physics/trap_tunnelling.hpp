// Trap-assisted tunnelling (Unit 67): Hurkx's field-enhancement factor of the SRH lifetimes.
//
// G. A. M. Hurkx, D. B. M. Klaassen and M. P. G. Knuvers, IEEE Trans. Electron Devices 39, 331
// (1992), as given by the ATLAS manual (eqs 3-71 to 3-75): a lifetime tau becomes tau / (1 + Gamma),
//
//     Gamma = a int_0^1 exp(a u - K u^(3/2)) du,   a = dE / kT,
//     K = (4/3) sqrt(2 m m0 dE^3) / (q hbar F),
//
// dE the tunnelling energy range of the carrier (Unit 67: band edge to the intrinsic level), m
// its tunnelling mass in m0, F the field. Evaluated with u = v^2,
//
//     Gamma = 2a int_0^1 v exp(a v^2 - K v^3) dv,   dGamma/dF = (K/F) 2a int_0^1 v^4 exp(...) dv,
//
// whose integrand is analytic, by 16-panel 8-point Gauss-Legendre on [0, v_max], v_max = 1 or the
// point past the maximum where the exponent has fallen 50 below it (the rest is below e^-50 of
// the peak). Gamma -> a Gamma(5/3) K^(-2/3) as F -> 0 (so dGamma/dF is unbounded there) and
// -> 2 sqrt(3 pi) (F/F_G) exp((F/F_G)^2), F_G = sqrt(24 m m0 (kT)^3) / (q hbar), where the peak
// lies well inside (Laplace's method).
#pragma once

namespace NiTCAD::physics {

struct TrapTunnellingFactor {
    double gamma;  // Gamma, dimensionless
    double d_dF;   // dGamma/dF [per V/cm]
};

// Gamma and dGamma/dF at field F [V/cm]. Both 0 at F = 0 and where K overflows (F so small that
// Gamma would be below the smallest double's scale of interest).
// Preconditions (NITCAD_EXPECTS): a finite in (0, 700]; mass_ratio and T_K finite and > 0; F
// finite and >= 0.
[[nodiscard]] TrapTunnellingFactor hurkx_enhancement(double a, double mass_ratio, double T_K,
                                                     double field_V_per_cm);

}  // namespace NiTCAD::physics
