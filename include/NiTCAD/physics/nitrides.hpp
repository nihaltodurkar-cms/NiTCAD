// Wurtzite nitride parameter sets (Unit 68, slice 2, MAD8; DECISIONS.md amendment 6).
//
// Sources: M2', Levinshtein, Rumyantsev and Shur, eds., Properties of Advanced Semiconductor
// Materials (Wiley 2001), ch. 1-3 (permittivity, GaN gap and Varshni, masses, radiative
// coefficient, 300 K thermal conductivity, specific heat, density); and the owner's values of
// 2026-11-27 (amendment 6: InN gap, mobility, saturation, lifetimes, Auger, dopant levels,
// GaN impact ionization, thermal exponents, conduction-band offset ratios); the AlN gap of J. Piprek,
// ed., Nitride Semiconductor Devices (Wiley-VCH 2007), Section 2. Isotropic (MAQ6).
//
// Band alignment (amendment 6, option A): GaN's affinity is 4.10 eV and every other nitride's
// conduction edge sits Q_c of the gap difference from GaN's at 300 K:
//     chi = chi_GaN - Q_c [Eg(300) - Eg_GaN(300)],  Q_c 0.70 toward AlN, 0.58 toward InN,
// an effective TCAD alignment, held at its 300 K value at every temperature.
#pragma once

#include <expected>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

inline constexpr double gan_affinity_eV = 4.10;
inline constexpr double algan_conduction_share = 0.70;  // Q_c, AlGaN/GaN
inline constexpr double ingan_conduction_share = 0.58;  // Q_c, InGaN/GaN

[[nodiscard]] SemiconductorParameters gan_parameters();
[[nodiscard]] SemiconductorParameters aln_parameters();
[[nodiscard]] SemiconductorParameters inn_parameters();

}  // namespace NiTCAD::physics

namespace NiTCAD::physics {

// `algan <x>` (Al fraction) and `ingan <x>` (In fraction), 0 <= x <= 1; x = 0 is gan_parameters()
// and x = 1 aln_parameters() or inn_parameters(), exactly. Between (amendment 6 and the owner's
// rulings of 2026-11-27):
//     Eg(x, 300) = (1 - x) Eg_GaN + x Eg_B - b x (1 - x),  b 0.8 (AlGaN), 1.4 (InGaN) eV;
//     the gap's temperature change linear in x (AlloyGapParameters); chi by Q_c (above);
//     m_e(x) = (1 - x) m_GaN + x m_B - b_m x (1 - x), b_m 4.8e-3 (AlGaN), 0 (InGaN); Nv linear;
//     mobility by AlloyMobilityParameters (C_mu 40 and 1e6 cm^2/(V s)); permittivity, lifetimes,
//     Auger, radiative, saturation, dopant levels and heat capacity linear; thermal conductivity
//     off (electrothermal refused: alloy scattering makes a linear rule wrong); impact ionization
//     off. Errors: invalid_input for x outside [0, 1].
[[nodiscard]] std::expected<SemiconductorParameters, base::Error> algan_parameters(double x);
[[nodiscard]] std::expected<SemiconductorParameters, base::Error> ingan_parameters(double x);

}  // namespace NiTCAD::physics
