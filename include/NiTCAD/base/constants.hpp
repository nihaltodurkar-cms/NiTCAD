// Physical constants and the unit convention used throughout NiTCAD.
//
// Unit convention (the standard TCAD convention, unchanged from the legacy project):
//     length           cm
//     concentration    cm^-3
//     potential        V
//     current density  A/cm^2
//     time             s
// Quantities that do not follow it carry their unit in the name (eps0_F_per_cm).
//
// Values are CODATA 2018, pinned deliberately. CODATA 2022 changed eps0 and the electron
// mass at the 1e-9 level; adopting it would silently move validated results.
#pragma once

namespace NiTCAD::base {

inline constexpr double q_C = 1.602176634e-19;             // elementary charge [C] (exact)
inline constexpr double k_B_J_per_K = 1.380649e-23;        // Boltzmann constant [J/K] (exact)
inline constexpr double k_B_eV_per_K = k_B_J_per_K / q_C;  // Boltzmann constant [eV/K]
inline constexpr double eps0_F_per_cm = 8.8541878128e-14;  // vacuum permittivity [F/cm]
// The reduced Planck constant h / 2 pi [J s], h = 6.62607015e-34 exact: the double nearest to
// 1.054 571 817 646 156 5...e-34. CODATA prints it to 10 figures, 1.054571817e-34, which is 6.1e-10 low
// (DECISIONS Unit 55 12, FQ1; finding (i) of TVB-01).
inline constexpr double hbar_J_s = 1.0545718176461565e-34;
inline constexpr double m0_kg = 9.1093837015e-31;          // electron rest mass [kg]

// Thermal voltage kT/q [V]. 25.852 mV at T = 300 K.
[[nodiscard]] constexpr double thermal_voltage(double temperature_K) noexcept {
    return k_B_J_per_K * temperature_K / q_C;
}

}  // namespace NiTCAD::base
