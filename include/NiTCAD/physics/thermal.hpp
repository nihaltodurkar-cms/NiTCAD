// Lattice heat and the temperature derivatives of the band, mobility and generation models (Unit 23,
// electrothermal coupling; legacy pytcad/thermal.py, thermal_grid.py, core/src/thermal/grid.cpp).
//
// Heat conduction: kappa(T) = kappa_300 (T / 300)^-a (legacy materials.kappa_th, silicon 1.48 W/(cm K)
// and a = 1.33). Heat storage rho c (thermal transients), temperature independent.
//
// Thermopower (Wachutka, IEEE Trans. CAD 9, 1141 (1990); Lundstrom, Fundamentals of Carrier
// Transport, sec. 3.4): with a momentum relaxation time tau ~ E^r, a current carries per electron
// the energy E_c + (r + 5/2) kT F_{r+3/2}(x) / F_{r+1/2}(x) (x = (E_Fn - E_c) / kT; Boltzmann:
// E_c + (r + 5/2) kT), holes likewise below E_v. In the drift-diffusion current the thermopower and
// the temperature dependence of Nc ~ T^(3/2) combine (assemble/electrothermal.hpp) into a thermal
// diffusion term k h n grad T with
//     h = (r + 5/2) F_{r+3/2}(x) / F_{r+1/2}(x) - 3/2,
// so 1 + r for Boltzmann statistics and, for Fermi-Dirac with r = -1/2 (acoustic phonons, the one
// exponent Fermi-Dirac is implemented for), 2 F_1(x) / F_0(x) - 3/2. F_j are the normalized
// complete Fermi-Dirac integrals, F_j -> e^x as x -> -inf; F_0 = ln(1 + e^x), F_1 = -Li_2(-e^x),
// dF_1/dx = F_0, dF_0/dx = F_-1 = 1 / (1 + e^-x).
//
// Temperature derivatives (d / dT at fixed doping and field) of the legacy formulas the other
// headers evaluate; each is the exact derivative of that formula:
//     Eg(T) = Eg0 - alpha T^2 / (T + beta)        dEg/dT = -alpha T (T + 2 beta) / (T + beta)^2
//     Nc, Nv ~ T^(3/2)                            d ln N / dT = 3 / (2 T)
//     n_i = sqrt(Nc Nv) exp(-Eg / 2kT)            d ln n_i / dT = 3 / (2 T) + (Eg - T dEg/dT) / (2 k T^2)
//     n_ie = n_i exp(dEg_BGN / 2kT)               d ln n_ie / dT = d ln n_i / dT - dEg_BGN / (2 k T^2)
//     mu_CT = mu_min + (mu_max (T/300)^b - mu_min) / (1 + (N / N_ref)^alpha)
//                                                 dmu/dT = b mu_max (T/300)^b / T / (1 + (N/N_ref)^alpha)
//     gamma_ii = tanh(h / 300) / tanh(h / T), h = hbar w / 2k
//                                                 dgamma/dT = tanh(h / 300) h / (T^2 sinh^2(h / T))
//     v_emission ~ T^(1/2) (both forms)           dv/dT = v / (2 T)
// The electron affinity is temperature independent (as in semiconductor.hpp), so in the vacuum-level
// gauge the conduction-band edge does not move with T and the valence-band edge carries all of the
// gap's change.
#pragma once

#include <expected>
#include <string_view>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

// Whether a parameter set carries thermal data (a conductivity above 0).
[[nodiscard]] constexpr bool has_thermal_data(const ThermalParameters& t) noexcept {
    return t.conductivity_W_cmK > 0.0;
}

// Errors (invalid_input; the message names `owner` and the field, the context value is it): every
// value not finite; a conductivity below 0; with a conductivity above 0, a heat capacity not
// positive.
[[nodiscard]] std::expected<void, base::Error> check_thermal_parameters(const ThermalParameters& t,
                                                                        std::string_view owner);

struct ValueSlope {
    double value;
    double d_dT;  // d value / d T [per K]
};

// kappa(T) [W/(cm K)] and its derivative. Precondition (NITCAD_EXPECTS): T finite and positive.
[[nodiscard]] ValueSlope thermal_conductivity(const ThermalParameters& t, double temperature_K);

// dEg/dT [eV/K] of band_gap_eV. Precondition (NITCAD_EXPECTS): T finite and positive.
[[nodiscard]] double band_gap_slope_eV_per_K(const Semiconductor& m, double temperature_K);

// d ln n_i / dT [1/K] of intrinsic_density, and of effective_intrinsic_density (band-gap
// narrowing dEg_BGN [eV] given). Precondition (NITCAD_EXPECTS): T finite and positive.
[[nodiscard]] double intrinsic_density_log_slope(const Semiconductor& m, double temperature_K,
                                                 double narrowing_eV = 0.0);

// dmu/dT [cm^2/(V s K)] of caughey_thomas_mobility. Preconditions as that function.
[[nodiscard]] double caughey_thomas_mobility_slope(const Semiconductor& m, Carrier carrier,
                                                   double total_impurity_cm3,
                                                   double temperature_K);

// d gamma / dT [1/K] of impact_ionization_temperature_factor (0 for a phonon energy of 0).
// Preconditions as that function.
[[nodiscard]] double impact_ionization_temperature_factor_slope(double phonon_energy_eV,
                                                                double temperature_K);

// d v / dT [cm/(s K)] of emission_velocity_cm_s(m, carrier, T). Preconditions as that function.
[[nodiscard]] double emission_velocity_slope(const Semiconductor& m, Carrier carrier,
                                             double temperature_K);

// F_-1, F_0 and F_1 at x. Below x = -1 from the alternating series in e^x (F_j = e^x S_j, the
// series S_j summed so that e^x itself never multiplies them: the ratios below are finite where
// e^x underflows); on [-1, 1] from F_0 = ln 2 + x / 2 + ln cosh(x / 2) and its integral, the Taylor
// series of ln cosh y = sum_n (-1)^(n+1) (4^n - 1) zeta(2n) y^(2n) / (n pi^(2n)); above 1 by the
// reflections F_0(x) = x + F_0(-x), F_1(x) = x^2 / 2 + pi^2 / 6 - F_1(-x). F_0 and F_1 underflow
// with e^x below x of about -745; F_-1 = 1 / (1 + e^-x).
struct FermiOrdersZeroOne {
    double minus_one, zero, one;
};
[[nodiscard]] FermiOrdersZeroOne fermi_orders_zero_one(double x) noexcept;

// The thermal diffusion factor of Fermi-Dirac statistics with r = -1/2 (see the header comment):
// h(x) = 2 F_1(x) / F_0(x) - 3/2 and dh/dx = 2 - 2 F_1 F_-1 / F_0^2. Tends to 1/2 (the Boltzmann
// 1 + r) as x -> -inf, with dh/dx -> 0, and to x - 3/2 for large x. Finite for every finite x.
struct ThermalDiffusionFactor {
    double value;
    double d_x;
};
[[nodiscard]] ThermalDiffusionFactor fermi_dirac_thermal_diffusion(double x) noexcept;

}  // namespace NiTCAD::physics
