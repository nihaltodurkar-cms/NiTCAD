// Shockley-Read-Hall recombination and the Scharfetter doping-dependent lifetimes
// (legacy materials.lifetime_scharfetter and materials.recombination / kernels.hpp
// recombination_boltzmann and recombination_fd, SRH part).
//
// SRH with a mid-gap trap (n1 = p1 = n_ie):
//
//     R = (n p - E) / (tau_p (n + n_ie) + tau_n (p + n_ie))
//
// where E = n_eq p_eq is the equilibrium product from statistics.hpp (n_ie^2 under Boltzmann
// statistics), so R vanishes at equilibrium. R > 0 is net recombination, R < 0 net generation.
// Returned partials are exact and include the chain through E:
//
//     dR/dn = ((p - dE/dn) D - (n p - E) tau_p) / D^2,   D the denominator,
//
// and likewise for p. With the Boltzmann E (partials zero) this is the legacy Boltzmann form
// operation for operation; with a Fermi-Dirac E it is the legacy recombination_fd form.
// Auger recombination (Unit 11) and radiative recombination (Unit 15) are below; the assembler
// adds them to SRH.
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::physics {

// tau = tau0 / (1 + N / N_ref) [s], with N the total ionised impurity N_A + N_D.
// Preconditions (NITCAD_EXPECTS): total_impurity_cm3 is finite and >= 0. The legacy clamp of N
// to at least 1 cm^-3 is not carried over (see mobility.hpp).
[[nodiscard]] double scharfetter_lifetime(const Semiconductor& m, Carrier carrier,
                                          double total_impurity_cm3);

struct RecombinationRate {
    double rate;  // R [cm^-3 s^-1] for inputs in cm^-3 and s
    double d_dn;  // dR/dn [s^-1]
    double d_dp;  // dR/dp [s^-1]
};

// A per-node, per-Newton-iteration kernel, so it has no checks. The result is the formula above;
// it is finite when the denominator is positive, which holds for n, p >= 0 and positive n_ie,
// tau_n, tau_p. Keeping n and p in that range during Newton overshoot is the caller's job.
// Homogeneous of degree one in the concentrations: scaling n, p, n_ie by s and E by s^2 scales R
// by s and leaves the partials unchanged.
[[nodiscard]] constexpr RecombinationRate srh_recombination(double n, double p,
                                                            EquilibriumProduct np_eq,
                                                            double n_ie, double tau_n,
                                                            double tau_p) noexcept {
    const double excess = n * p - np_eq.value;
    const double den = tau_p * (n + n_ie) + tau_n * (p + n_ie);
    return {
        excess / den,
        ((p - np_eq.d_dn) * den - excess * tau_p) / (den * den),
        ((n - np_eq.d_dp) * den - excess * tau_n) / (den * den),
    };
}

// Auger recombination (legacy materials.recombination, auger branch):
//
//     R = (Cn n + Cp p) (n p - E),   dR/dn = Cn (n p - E) + (Cn n + Cp p) (p - dE/dn),
//
// and likewise for p, with E the equilibrium product (so R vanishes at equilibrium). Cubic in the
// concentrations, not homogeneous of degree one like SRH: call it with physical densities
// (cm^-3, and E in cm^-6) and the coefficients in cm^6/s; R is then in cm^-3 s^-1.
[[nodiscard]] constexpr RecombinationRate auger_recombination(double n, double p,
                                                              EquilibriumProduct np_eq,
                                                              double Cn, double Cp) noexcept {
    const double excess = n * p - np_eq.value;
    const double C = Cn * n + Cp * p;
    return {C * excess, Cn * excess + C * (p - np_eq.d_dn), Cp * excess + C * (n - np_eq.d_dp)};
}

// Radiative (band-to-band) recombination (Unit 15; not in the legacy):
//
//     R = B (n p - E),   dR/dn = B (p - dE/dn),   dR/dp = B (n - dE/dp),
//
// with B the material's radiative coefficient and E the equilibrium product. Quadratic in the
// concentrations: call it with physical densities (cm^-3, E in cm^-6) and B in cm^3/s.
[[nodiscard]] constexpr RecombinationRate radiative_recombination(double n, double p,
                                                                  EquilibriumProduct np_eq,
                                                                  double B) noexcept {
    return {B * (n * p - np_eq.value), B * (p - np_eq.d_dn), B * (n - np_eq.d_dp)};
}

}  // namespace NiTCAD::physics
