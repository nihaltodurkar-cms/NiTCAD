// The energy-balance boundary layer of a uniform resistor (Unit 55 TVR-4; DECISIONS 16.18 (a), DT-Q15 A, 16.42;
// TVB-13). Derived here (PD); base-free: every constant and parameter is an input.
//
// The uniform state: electrons at density N, field E0 > 0, temperature T_h (with E_eff(T_h) = E0), mobility mu at
// T_h and its slope mu' = d mu / dT, current J = q N mu E0. Its perturbations ~ e^{s x} at fixed current
// (delta J = 0), unknowns (delta psi, delta n, delta T), with kappa = c (k^2 / q) N mu T_h (c = 5/2, ED3) and
// G = (3/2) k N / tau_w, satisfy (ED2, Poisson, ED3)
//     current: -mu N q s dpsi + (mu q E0 + mu k T_h s) dn + (mu' q N E0 + mu k N s) dT = 0,
//     Poisson: eps s^2 dpsi - q dn = 0,
//     energy:  J s dpsi + (3/2) k (T_h - T_L) / tau_w dn + (-kappa s^2 - c (k / q) J s + G) dT = 0.
// Their determinant is s P(s) up to a constant, with the quartic (eliminating dpsi)
//     P(s) = [mu k T_h s^2 + mu q E0 s - mu N q^2 / eps] (-kappa s^2 - c (k / q) J s + G)
//            - (mu' q N E0 + mu k N s) (J q / eps + (3/2) k (T_h - T_L) s / tau_w).
// lambda = 2 / min |Re s| over P's roots: the slowest mode's length, doubled for safety.
//
// The roots: Durand-Kerner (Weierstrass) on the monic quartic in complex double-double, after scaling s by the
// geometric mean of the roots' magnitudes. Second method: the determinant of the unreduced 3 x 3 system at each
// root, relative to the product of its rows' norms (independent of the elimination and of the root finder).
#pragma once

#include <array>
#include <complex>
#include <expected>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct ComplexDD {
    DD re, im;
};

struct PolynomialRoots {
    std::vector<ComplexDD> roots;
    int iterations = 0;
    double residual = 0.0;  // max |p(z)| / sum |c_k| |z|^k over the roots
};

// The roots of s^n + c[n-1] s^{n-1} + ... + c[0] (c.size() = n >= 1, real coefficients, c[0] != 0), by
// Durand-Kerner in double-double. Fails if it does not converge in 500 iterations.
[[nodiscard]] std::expected<PolynomialRoots, std::string> durand_kerner(const std::vector<DD>& c);

struct LayerInput {
    double q_C = 0.0, k_J_per_K = 0.0, eps_F_per_cm = 0.0;
    double density_cm3 = 0.0;          // N
    double field_V_per_cm = 0.0;       // E0 > 0
    double mobility_cm2_Vs = 0.0;      // mu at T_h
    double mobility_slope = 0.0;       // d mu / dT at T_h [cm^2 / (V s K)]
    double temperature_K = 0.0;        // T_h
    double lattice_K = 0.0;            // T_L
    double tau_w_s = 0.0;
    double heat_coefficient = 2.5;     // c (ED3)
};

struct LayerResult {
    std::array<std::complex<double>, 4> roots;  // [1/cm]
    double slowest_rate_per_cm = 0.0;           // min |Re s|
    double length_cm = 0.0;                     // lambda = 2 / min |Re s|
    double root_residual = 0.0;                 // Durand-Kerner's
    double determinant_residual = 0.0;          // the second method: max over the roots
    int iterations = 0;
};

[[nodiscard]] std::expected<LayerResult, std::string> energy_layer(const LayerInput& in);

}  // namespace NiTCAD::verify::reference
