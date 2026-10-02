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
// Auger recombination is deferred (ARCHITECTURE.md section 11, Units 11+).
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

}  // namespace NiTCAD::physics
