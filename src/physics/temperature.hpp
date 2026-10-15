// Private to the physics layer: the one temperature precondition, the one temperature scaling of
// the Caughey-Thomas mu_max, and the Unit 35 band forms (Passler gap, density-of-states mass
// polynomials) with their temperature slopes, shared by the model functions, check_temperature
// and thermal.cpp.
#pragma once

#include <cmath>
#include <numbers>

#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics::detail {

inline void expect_temperature(double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(temperature_K) && temperature_K > 0.0);
}

// mu_max (T / 300)^T_exponent [cm^2/(V s)]. No check; callers check the temperature.
[[nodiscard]] inline double mu_max_at(const CaugheyThomasParameters& ct, double temperature_K) {
    return ct.mu_max * std::pow(temperature_K / 300.0, ct.T_exponent);
}

struct BandSlope {
    double value;
    double d_dT;  // d value / d T [per K]
};

// Passler's Eg(T) [eV] and dEg/dT [eV/K] (semiconductor.hpp, PasslerGapParameters). No check.
[[nodiscard]] inline BandSlope passler_gap(double Eg0_eV, const PasslerGapParameters& g,
                                           double T) {
    const double D2 = g.Delta * g.Delta;
    const double u = g.Theta_K / T;
    const double em = std::expm1(u);  // exp(u) - 1
    const double c = (1.0 - 3.0 * D2) / em;
    // dc/dT = (1 - 3 Delta^2) (u / T) exp(u) / (exp(u) - 1)^2, written so as not to overflow.
    const double dc = (1.0 - 3.0 * D2) * (u / T) / (em * -std::expm1(-u));
    const double x = 2.0 * T / g.Theta_K;
    const double a2 = std::numbers::pi * std::numbers::pi / (3.0 * (1.0 + D2));
    const double a3 = (3.0 * D2 - 1.0) / 4.0;
    const double x2 = x * x, x3 = x2 * x, x4 = x3 * x, x5 = x4 * x, x6 = x5 * x;
    const double S = 1.0 + a2 * x2 + a3 * x3 + 8.0 / 3.0 * x4 + x6;
    const double dS = 2.0 * a2 * x + 3.0 * a3 * x2 + 32.0 / 3.0 * x3 + 6.0 * x5;  // dS / dx
    const double R = std::pow(S, 1.0 / 6.0);
    const double dR = R / (6.0 * S) * dS * (2.0 / g.Theta_K);  // dR / dT
    const double k = g.alpha_eV_per_K * g.Theta_K;
    return {Eg0_eV - k * (c + 1.5 * D2 * (R - 1.0)), -k * (dc + 1.5 * D2 * dR)};
}

// (m_d / m0)^(3/2) = c[0] T^3 + c[1] T^2 + c[2] T + c[3] and its T derivative. No check.
[[nodiscard]] inline BandSlope dos_mass(const double (&c)[4], double T) {
    return {((c[0] * T + c[1]) * T + c[2]) * T + c[3], (3.0 * c[0] * T + 2.0 * c[1]) * T + c[2]};
}

// 2 (2 pi m0 k / h^2)^(3/2) [cm^-3 K^-3/2] as Couderc et al. (2014) print it, eqs (4) and (5).
inline constexpr double dos_constant = 4.83e15;

}  // namespace NiTCAD::physics::detail
