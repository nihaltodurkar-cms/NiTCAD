// Zincblende III-V parameter sets (Unit 68, slice 3, MAD10; DECISIONS.md amendment 7).
//
// Sources (the owner's values and rulings of 2026-11-27, amendment 7, from):
//   M1, I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, J. Appl. Phys. 89, 5815 (2001): gaps and Varshni
//     (the fundamental valley: X for GaP and AlAs), electron masses, Luttinger gamma_1 and gamma_2,
//     valence-band offsets, alloy bowing;
//   M. Sotoodeh, A. H. Khalid, A. A. Rezazadeh, J. Appl. Phys. 87, 2890 (2000), Table I: the
//     Caughey-Thomas mobility of InP, InAs, GaP, AlAs and GaAs (GaSb and InSb from a temperature-
//     extended fit, the InSb exponents provisional);
//   V. Palankovski (thesis, TU Wien 2000), Tables 3.30-3.32, and Quay et al. (2000), Table 1: saturation
//     velocities and their temperature coefficient A; the alloy C_mu;
//   handbook and owner values: permittivity, lifetimes, Auger, radiative, thermal data.
// beta 2 (velocity saturation) is the engine's model default, not a sourced value. Isotropic.
//
// Band alignment (owner ruling 1): each set's affinity comes from M1's valence-band offsets with
// GaAs's 4.07 eV as the anchor, at each set's 300 K fundamental gap:
//     chi = 4.07 + Eg_GaAs(300) + VBO_GaAs - VBO - Eg(300),
// a derived, internally consistent lineup (not measured affinities), held at every temperature.
// Densities of states: the conduction band's from the fundamental valley's density-of-states mass
// (GaP and AlAs: their X ellipsoids, three valleys, (m_l m_t^2)^(1/3) 3^(2/3)); the valence band's
// from the [001] heavy- and light-hole masses 1 / (gamma_1 -+ 2 gamma_2), combined as
// (m_hh^(3/2) + m_lh^(3/2))^(2/3). Dopants fully ionized; band-gap narrowing and impact ionization
// off. A lifetime of 0 is not given: SRH refuses the material (InSb, the alloys).
#pragma once

#include <expected>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

inline constexpr double iii_v_affinity_anchor_eV = 4.07;  // GaAs

// GaAs on M1's band parameters with Sotoodeh's mobility (`gaas_2001`), the alloys' GaAs end; the
// legacy `gaas` set stays as it was (owner ruling 2).
[[nodiscard]] SemiconductorParameters gaas_2001_parameters();
[[nodiscard]] SemiconductorParameters inp_parameters();
[[nodiscard]] SemiconductorParameters inas_parameters();
[[nodiscard]] SemiconductorParameters gap_parameters();
[[nodiscard]] SemiconductorParameters alas_parameters();
[[nodiscard]] SemiconductorParameters gasb_parameters();
[[nodiscard]] SemiconductorParameters insb_parameters();

// `ingaas <x>`, In(x)Ga(1-x)As, and `inalas <x>`, In(x)Al(1-x)As, x the In fraction in [0, 1]; x = 0
// is gaas_2001_parameters() or alas_parameters() and x = 1 inas_parameters(), exactly. Between:
//     P(x) = (1 - x) P_A + x P_B - b x (1 - x)
// for the Gamma gap (b 0.477 and 0.70 eV), the X gap (InAlAs, b 0), the Gamma electron mass (0.0091
// and 0.049), the [001] hole masses (InGaAs -0.145 and 0.0202; InAlAs 0) and the valence-band
// offset (-0.38 and -0.64 eV); the conduction band the lower of the Gamma and X valleys at 300 K
// (InAlAs: X below x of about 0.37); the gap's temperature change linear in x (AlloyGapParameters),
// a satellite valley shifting with the fundamental one; mobility by AlloyMobilityParameters (C_mu
// 1e6 cm^2/(V s), both carriers); thermal conductivity (1 - x) k_A + x k_B - C_k x (1 - x), C_k 1.41
// and 2.27 W/(cm K), (T / 300)^-1.17, none (electrothermal refused) where that is not positive
// (InAlAs, x about 0.54 to 0.74); lifetimes not given; permittivity, saturation, Auger, radiative
// and heat capacity linear. InGaAs carries the direct-gap tunnelling masses (Gamma electron, [001]
// light hole). Errors: invalid_input for x outside [0, 1].
[[nodiscard]] std::expected<SemiconductorParameters, base::Error> ingaas_parameters(double x);
[[nodiscard]] std::expected<SemiconductorParameters, base::Error> inalas_parameters(double x);

}  // namespace NiTCAD::physics
