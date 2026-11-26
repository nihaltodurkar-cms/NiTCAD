// Total-dose oxide charge (Unit 69; DECISIONS.md TDD1, TDQ1-TDQ4).
//
// Ionizing dose D [rad(SiO2)] leaves holes trapped in a gate oxide of thickness t_ox, taken as a
// sheet at the semiconductor interface (TDQ1 a):
//
//     N_ot = f_ot f_y(E) g0 D t_ox  [cm^-2],
//
// f_ot the hole trapping fraction (0..1), g0 the pairs per cm^3 per rad, and f_y the Co-60 charge
// yield at the oxide field E [V/cm] during irradiation (Oldham and McLean, IEEE TNS 50(3), 2003;
// quoted, the paper not in hand: TDQ3):
//
//     f_y(E) = ((E + E0) / (E + E1))^m.
//
// The deck adds N_ot to the interface's or gate's fixed charge, and builds the interface traps
// k_it D as two bands of half-width dose_trap_half_width_eV about E_i (TDQ4).
#pragma once

namespace NiTCAD::physics {

// SiO2: 2.2 g/cm^3 * 6.2415e13 eV/(g rad) / 17 eV per pair = 8.08e12, the literature's 8.1e12.
inline constexpr double dose_pair_density_cm3_rad = 8.1e12;
inline constexpr double dose_yield_e0_V_cm = 0.1e6;
inline constexpr double dose_yield_e1_V_cm = 1.35e6;
inline constexpr double dose_yield_exponent = 0.9;
inline constexpr double dose_trap_half_width_eV = 0.5;

// f_y(E); E finite and >= 0.
[[nodiscard]] double charge_yield(double field_V_cm);

// N_ot [cm^-2]; dose > 0, field >= 0, thickness > 0, all finite; trapping fraction in [0, 1].
[[nodiscard]] double oxide_trapped_charge_cm2(double dose_rad, double field_V_cm, double thickness_cm,
                                              double trapping_fraction);

}  // namespace NiTCAD::physics
