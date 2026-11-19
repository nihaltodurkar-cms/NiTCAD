// Quantum confinement by the modified local density approximation (Unit 47; DECISIONS.md,
// Unit 47): G. Paasch, H. Ubensee, phys. stat. sol. (b) 113, 165 (1982), in its nondegenerate
// form (W. Hansch et al., Solid-State Electron. 32, 839 (1989)), as reproduced by M. Ieong et al.
// (SISPAD 1998) and S. Entner (TU Wien thesis, eqs. 2.106-2.108); the papers themselves were not
// read. Near an insulator interface the density of states of a band is reduced by
//
//     f(z) = 1 - exp(-z^2 / lambda^2),    lambda = hbar / sqrt(2 m k T),
//
// z the distance to the interface and m the carrier's quantization mass. Under Boltzmann
// statistics this is MLDA's j0-kernel integral exactly; under Fermi-Dirac it is Hansch's
// approximation of it.
#pragma once

namespace NiTCAD::physics {

// Silicon's quantization masses [m0] (Ioffe: the longitudinal electron mass and the heavy-hole
// mass), the defaults of PhysicsModels::mlda_mass_n and mlda_mass_p.
inline constexpr double mlda_silicon_electron_mass = 0.98;
inline constexpr double mlda_silicon_hole_mass = 0.49;

// lambda = hbar / sqrt(2 m m0 k T) [cm]. Preconditions (NITCAD_EXPECTS): mass_m0 finite and
// positive; temperature_K finite and positive.
[[nodiscard]] double mlda_wavelength_cm(double mass_m0, double temperature_K);

// f(z) at distance_cm >= 0 (finite); 0 at the interface. Precondition: wavelength_cm finite and
// positive.
[[nodiscard]] double mlda_factor(double distance_cm, double wavelength_cm);

// The mean of f over z in [0, a]: 1 - (sqrt(pi) lambda / (2 a)) erf(a / lambda), by its series
// where that form cancels. Positive for a > 0; the factor of a node on the interface (a lumped
// gate), a its depth. Preconditions: depth_cm and wavelength_cm finite and positive.
[[nodiscard]] double mlda_mean_factor(double depth_cm, double wavelength_cm);

// Density-gradient confinement (Unit 62; DECISIONS.md QP1-QP6). The quantum potential of carrier c is
// Lambda_c / kT = -beta_c (laplacian sqrt c) / sqrt c, with
//
//     beta = gamma hbar^2 / (6 m m0 k T)  [cm^2],
//
// and an insulator barrier dE in front of a semiconductor lets sqrt c decay as exp(-x / x_p) into it,
//
//     x_p = hbar / sqrt(2 m_ox m0 q dE)  [cm].
//
// The defaults are placeholders until the source check QP0 records the theoretical coefficient and
// the masses it pairs with: gamma = 1 and MLDA's masses.
inline constexpr double dg_default_gamma = 1.0;
inline constexpr double dg_default_electron_mass = mlda_silicon_electron_mass;
inline constexpr double dg_default_hole_mass = mlda_silicon_hole_mass;

// beta [cm^2]. Preconditions (NITCAD_EXPECTS): gamma and mass_m0 finite and positive; temperature_K
// finite and positive.
[[nodiscard]] double dg_coefficient_cm2(double gamma, double mass_m0, double temperature_K);

// x_p [cm]. Preconditions: mass_m0 and barrier_eV finite and positive.
[[nodiscard]] double barrier_decay_length_cm(double mass_m0, double barrier_eV);

}  // namespace NiTCAD::physics
