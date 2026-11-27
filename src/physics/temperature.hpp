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
// Unit 68 slice 2: mu_min, N_ref and alpha at T (unchanged with their exponents 0).
[[nodiscard]] inline double mu_min_at(const CaugheyThomasParameters& ct, double temperature_K) {
    return ct.mu_min_T_exponent == 0.0 ? ct.mu_min : ct.mu_min * std::pow(temperature_K / 300.0, ct.mu_min_T_exponent);
}

// One Caughey-Thomas set's mobility and its temperature slope at N [cm^-3] and T [K]. With the
// slice-2 exponents 0 both are the legacy expressions, operation for operation.
[[nodiscard]] inline double ct_mobility(const CaugheyThomasParameters& ct, double N, double T) {
    const double mu_max = mu_max_at(ct, T);
    const double mu_min = mu_min_at(ct, T);
    const double N_ref = ct.N_ref_T_exponent == 0.0 ? ct.N_ref : ct.N_ref * std::pow(T / 300.0, ct.N_ref_T_exponent);
    const double alpha = ct.alpha_T_exponent == 0.0 ? ct.alpha : ct.alpha * std::pow(T / 300.0, ct.alpha_T_exponent);
    // Legacy operation order (materials.mobility_caughey_thomas).
    return mu_min + (mu_max - mu_min) / (1.0 + std::pow(N / N_ref, alpha));
}
[[nodiscard]] inline double ct_mobility_slope(const CaugheyThomasParameters& ct, double N, double T) {
    const double mu_max = mu_max_at(ct, T);
    if (ct.mu_min_T_exponent == 0.0 && ct.N_ref_T_exponent == 0.0 && ct.alpha_T_exponent == 0.0) {
        return ct.T_exponent * mu_max / T / (1.0 + std::pow(N / ct.N_ref, ct.alpha));
    }
    // mu = mu_min + (mu_max - mu_min) s, s = 1 / (1 + u), u = (N / N_ref)^alpha, each of mu_min,
    // mu_max, N_ref, alpha a power of T: d ln u / dT = (gamma3 ln(N / N_ref) - gamma2) alpha / T.
    const double mu_min = mu_min_at(ct, T);
    const double N_ref = ct.N_ref * std::pow(T / 300.0, ct.N_ref_T_exponent);
    const double alpha = ct.alpha * std::pow(T / 300.0, ct.alpha_T_exponent);
    const double u = N > 0.0 ? std::pow(N / N_ref, alpha) : 0.0;
    const double s = 1.0 / (1.0 + u);
    const double dln_u = N > 0.0 ? (ct.alpha_T_exponent * std::log(N / N_ref) - ct.N_ref_T_exponent) * alpha / T : 0.0;
    const double d_min = ct.mu_min_T_exponent * mu_min / T;
    const double d_max = ct.T_exponent * mu_max / T;
    return d_min + (d_max - d_min) * s - (mu_max - mu_min) * s * s * u * dln_u;
}

// A carrier's low-field mobility with the alloy rule (AlloyMobilityParameters), and its slope.
[[nodiscard]] inline double alloy_mobility(const CaugheyThomasParameters& a, const AlloyMobilityParameters& m,
                                           bool electron, double N, double T) {
    const double mu_a = ct_mobility(a, N, T);
    if (m.fraction == 0.0) return mu_a;
    const double x = m.fraction;
    const double mu_b = ct_mobility(electron ? m.electron_b : m.hole_b, N, T);
    const double C_mu = electron ? m.C_mu_electron : m.C_mu_hole;
    const double scatter = C_mu > 0.0 ? x * (1.0 - x) / C_mu : 0.0;
    return 1.0 / ((1.0 - x) / mu_a + x / mu_b + scatter);
}
[[nodiscard]] inline double alloy_mobility_slope(const CaugheyThomasParameters& a, const AlloyMobilityParameters& m,
                                                 bool electron, double N, double T) {
    if (m.fraction == 0.0) return ct_mobility_slope(a, N, T);
    const double x = m.fraction;
    const CaugheyThomasParameters& b = electron ? m.electron_b : m.hole_b;
    const double mu_a = ct_mobility(a, N, T), mu_b = ct_mobility(b, N, T);
    const double mu = alloy_mobility(a, m, electron, N, T);
    // d(1/mu) = -(1 - x) mu_a' / mu_a^2 - x mu_b' / mu_b^2; the scattering term is constant.
    return mu * mu * ((1.0 - x) * ct_mobility_slope(a, N, T) / (mu_a * mu_a) + x * ct_mobility_slope(b, N, T) / (mu_b * mu_b));
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

// The alloy gap correction (semiconductor.hpp, AlloyGapParameters) and its T derivative. No check.
[[nodiscard]] inline BandSlope alloy_gap(const AlloyGapParameters& g, double T) {
    if (g.fraction == 0.0) return {0.0, 0.0};
    const auto d = [T](double a, double b) -> BandSlope {
        const double s = T + b;
        return {a * 300.0 * 300.0 / (300.0 + b) - a * T * T / s, -a * T * (T + 2.0 * b) / (s * s)};
    };
    const BandSlope da = d(g.alpha_a_eV_per_K, g.beta_a_K), db = d(g.alpha_b_eV_per_K, g.beta_b_K);
    return {g.fraction * (db.value - da.value), g.fraction * (db.d_dT - da.d_dT)};
}

// 2 (2 pi m0 k / h^2)^(3/2) [cm^-3 K^-3/2] as Couderc et al. (2014) print it, eqs (4) and (5).
inline constexpr double dos_constant = 4.83e15;

}  // namespace NiTCAD::physics::detail
