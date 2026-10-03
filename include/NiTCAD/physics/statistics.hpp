// Carrier statistics (ARCHITECTURE.md section 5, item 9; R4): Boltzmann (Unit 5) and Fermi-Dirac
// for parabolic bands (Unit 14).
//
// Statistics enter the equations in three places, and each has a function here for each
// statistics, so the assembler never writes a statistics formula itself:
// 1. carrier density from the reduced potential (boltzmann_density, fermi_dirac_density);
// 2. the equilibrium product n_eq p_eq that recombination drives towards (EquilibriumProduct,
//    consumed by recombination.hpp; Fermi-Dirac makes it depend on n and p, hence the partials);
// 3. the charge-neutral equilibrium of an ohmic contact or a neutral bulk node
//    (boltzmann_neutral_equilibrium; fermi_dirac_neutral_equilibrium, a root solve).
// Fermi-Dirac adds a fourth, the degeneracy factor of a density (fermi_dirac_degeneracy), which
// the transport fluxes need (the assembler's business, assemble/drift_diffusion.hpp).
//
// Reduced potentials are in units of V_T: eta_n = (psi - phi_n) / V_T for electrons and
// eta_p = (phi_p - psi) / V_T for holes, with psi the electrostatic potential and phi_n, phi_p the
// quasi-Fermi potentials, all referenced to the intrinsic level so that n = p = n_ie at
// eta = 0. Legacy: n = nie exp(psi + s), p = nie exp(-psi - s) in scaled units
// (device1d.cpp); its band-offset gauge shift s belongs to heterojunctions and is not carried.
//
// Fermi-Dirac keeps that gauge. With N the band's effective density of states (Nc for electrons,
// Nv for holes) and g = ln(N / n_ie) (the band edge sits g V_T from the reference level),
//     n = N F_{1/2}(eta - g) = n_ie e^eta gamma(eta - g),   gamma = F_{1/2}(x) / e^x
// (fermi_dirac.hpp), which tends to the Boltzmann n_ie e^eta for eta - g << 0. With n_ie the
// effective one (band-gap narrowing), g_n + g_p = ln(Nc Nv / n_ie^2) is the narrowed gap over kT,
// so the narrowing is shared between the bands as in the Boltzmann case.
//
// Every function here is homogeneous in the concentrations: multiplying n_ie (and the doping)
// by a factor multiplies every returned concentration by the same factor and leaves eta
// unchanged (g is a ratio, so it does not change). So the assembler may pass concentrations
// already divided by its scale Ns.
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

// Fermi-Dirac (Unit 14). log_dos is g = ln(N / n_ie) of the carrier's band: ln(Nc / n_ie) for
// electrons, ln(Nv / n_ie) for holes.

// n = n_ie e^(eta + ln gamma(eta - g)) and dn/deta = n (1 + d ln gamma). A per-node kernel: no
// checks (NaN in, NaN out).
[[nodiscard]] DensityResult fermi_dirac_density(double n_ie, double log_dos, double eta) noexcept;

struct Degeneracy {
    double log_gamma;  // ln gamma at the density's own reduced energy: n = n_ie e^(eta + log_gamma)
    double d_density;  // d log_gamma / d density (<= 0), in 1 / (the unit of n_ie)
};

// The degeneracy factor of a carrier density (legacy fd_node_factors: L = ln nu, w = dL/dn): the
// x with N F_{1/2}(x) = density by inverse_fermi_half, then ln gamma(x), with
// d ln gamma / d density = d ln gamma / dx / (density (1 + d ln gamma / dx)).
// Below density / N = 1e-6 the series ln gamma = -log1p(a v + b v^2), v = density / N,
// a = 2^(-3/2), b = 1/4 - 3^(-3/2) is used instead (error below 1e-18; no inversion). It is
// analytic through v = 0, so it also accepts a density at or slightly below zero, as a finite-
// difference probe of a minority density may produce.
// A per-node, per-iteration kernel: no checks. A density that is not finite gives a non-finite
// result (so a diverging Newton iterate surfaces as a non-finite residual); n_ie must be positive
// and log_dos finite, as the assembler's are.
[[nodiscard]] Degeneracy fermi_dirac_degeneracy(double n_ie, double log_dos,
                                                double density) noexcept;

// n_eq p_eq = n_ie^2 gamma_n gamma_p at the carriers' own degeneracy factors (legacy npq_args):
// with a common Fermi level n p = Nc Nv F_{1/2}(eta_n) F_{1/2}(eta_p), which this is; partials
// value d log_gamma_n / dn and value d log_gamma_p / dp. Boltzmann limit: n_ie^2.
[[nodiscard]] inline EquilibriumProduct fermi_dirac_equilibrium_product(
    double n_ie, const Degeneracy& n, const Degeneracy& p) noexcept {
    const double value = n_ie * n_ie * std::exp(n.log_gamma + p.log_gamma);
    return {value, value * n.d_density, value * p.d_density};
}

// Equilibrium of a charge-neutral node under Fermi-Dirac statistics: n - p = net_doping with
// n = n_ie e^(eta + ln gamma(eta - g_n)) and p = n_ie e^(-eta + ln gamma(-eta - g_p)); a
// safeguarded Newton solve for the majority carrier, ln n = ln(C + p) (C >= 0, and the mirror for
// C < 0), from the Boltzmann solution. n, p and eta as in boltzmann_neutral_equilibrium.
// Preconditions (NITCAD_EXPECTS): net_doping finite; n_ie finite and positive; log_dos_n and
// log_dos_p finite.
[[nodiscard]] NeutralEquilibrium fermi_dirac_neutral_equilibrium(double net_doping, double n_ie,
                                                                 double log_dos_n,
                                                                 double log_dos_p);

}  // namespace NiTCAD::physics
