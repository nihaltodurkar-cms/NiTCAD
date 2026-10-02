// Carrier statistics (ARCHITECTURE.md section 5, item 7; R4). Boltzmann statistics only.
//
// Statistics enter the equations in three places, and each has a function here, so the
// assembler never writes a statistics formula itself and Fermi-Dirac can later be added beside
// these functions without changing any caller's signature:
// 1. carrier density from the reduced potential (boltzmann_density);
// 2. the equilibrium product n_eq p_eq that recombination drives towards (EquilibriumProduct,
//    consumed by recombination.hpp; Fermi-Dirac makes it depend on n and p, hence the partials);
// 3. the charge-neutral equilibrium of an ohmic contact or a neutral bulk node
//    (boltzmann_neutral_equilibrium; Fermi-Dirac needs a root solve there).
//
// Reduced potentials are in units of V_T: eta_n = (psi - phi_n) / V_T for electrons and
// eta_p = (phi_p - psi) / V_T for holes, with psi the electrostatic potential and phi_n, phi_p the
// quasi-Fermi potentials, all referenced to the intrinsic level so that n = p = n_ie at
// eta = 0. Legacy: n = nie exp(psi + s), p = nie exp(-psi - s) in scaled units
// (device1d.cpp); its band-offset gauge shift s belongs to heterojunctions and is not carried.
//
// Every function here is homogeneous in the concentrations: multiplying n_ie (and the doping)
// by a factor multiplies every returned concentration by the same factor and leaves eta
// unchanged. So the assembler may pass concentrations already divided by its scale Ns.
#pragma once

#include <cmath>

namespace NiTCAD::physics {

struct DensityResult {
    double density;  // same unit as n_ie
    double d_eta;    // d density / d eta
};

// n = n_ie exp(eta), for either carrier with its own reduced potential. Not clipped: exp
// overflows to +inf above eta ~ 709.78, and bounding a Newton overshoot is the caller's job
// (the legacy clips eta to +-700 inside its residual, ARCHITECTURE.md 6.7).
[[nodiscard]] inline DensityResult boltzmann_density(double n_ie, double eta) noexcept {
    const double density = n_ie * std::exp(eta);
    return {density, density};
}

// The equilibrium carrier product n_eq p_eq and its partials with respect to n and p.
struct EquilibriumProduct {
    double value;  // concentration^2
    double d_dn;
    double d_dp;
};

// Boltzmann: n_eq p_eq = n_ie^2, independent of n and p. n_ie^2 underflows (to subnormal, then 0)
// for n_ie below about 1.5e-154; then SRH no longer vanishes at equilibrium to full precision.
[[nodiscard]] constexpr EquilibriumProduct boltzmann_equilibrium_product(double n_ie) noexcept {
    return {n_ie * n_ie, 0.0, 0.0};
}

struct NeutralEquilibrium {
    double n;    // same unit as n_ie
    double p;
    double eta;  // reduced equilibrium potential psi / V_T, with n = n_ie exp(eta)
};

// Equilibrium densities of a charge-neutral node: n - p = net_doping (N_D - N_A) and
// n p = n_ie^2. The majority carrier comes from 0.5 (|C| + sqrt(C^2 + 4 n_ie^2)) and the minority
// from mass action as n_ie (n_ie / majority), so neither cancels (legacy Device1D::contact_value,
// which forms n_ie^2 / majority); eta = asinh(C / (2 n_ie)) (the legacy initial guess). No
// intermediate overflows or underflows where the result is representable: n, p and eta are finite
// whenever the majority density is, and p loses precision only near the subnormal range.
// Preconditions (NITCAD_EXPECTS): net_doping is finite; n_ie is finite and positive.
[[nodiscard]] NeutralEquilibrium boltzmann_neutral_equilibrium(double net_doping,
                                                               double n_ie);

}  // namespace NiTCAD::physics
