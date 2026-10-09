#include "NiTCAD/physics/thermal.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "temperature.hpp"

namespace NiTCAD::physics {

namespace {

base::Error thermal_error(std::string_view owner, std::string_view field,
                          std::string_view requirement, double value) {
    std::string message{owner};
    message += " thermal parameter ";
    message += field;
    message += requirement;
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = std::nullopt, .value = value}};
}

// zeta(s), s >= 4: the sum to N - 1 = 99 and the Euler-Maclaurin tail N^(1-s) / (s - 1) + N^-s / 2
// + s N^(-s-1) / 12 - s (s+1) (s+2) N^(-s-3) / 720 + s (s+1) (s+2) (s+3) (s+4) N^(-s-5) / 30240,
// whose remainder is below 1e-20 at s = 4.
double zeta_value(double s) {
    constexpr double N = 100.0;
    double sum = 0.0;
    for (int k = 99; k >= 1; --k) sum += std::pow(static_cast<double>(k), -s);
    const double t = std::pow(N, -s);
    return sum + N * t / (s - 1.0) + 0.5 * t + s * t / (12.0 * N) -
           s * (s + 1.0) * (s + 2.0) * t / (720.0 * N * N * N) +
           s * (s + 1.0) * (s + 2.0) * (s + 3.0) * (s + 4.0) * t / (30240.0 * std::pow(N, 5.0));
}

// The Taylor coefficients of ln cosh y = sum_n c_n y^(2n), c_n = (-1)^(n+1) (4^n - 1) zeta(2n) /
// (n pi^(2n)), n = 1..30 (zeta(2) = pi^2 / 6, so c_1 = 1/2; c_2 = -1/12, c_3 = 1/45).
const std::array<double, 30>& log_cosh_coefficients() {
    static const std::array<double, 30> c = [] {
        std::array<double, 30> a{};
        constexpr double pi2 = std::numbers::pi * std::numbers::pi;
        for (int n = 1; n <= 30; ++n) {
            const double zeta = n == 1 ? pi2 / 6.0 : zeta_value(2.0 * n);
            const double sign = n % 2 == 1 ? 1.0 : -1.0;
            a[static_cast<std::size_t>(n - 1)] =
                sign * (std::pow(4.0, n) - 1.0) * zeta / (n * std::pow(pi2, n));
        }
        return a;
    }();
    return c;
}

// x <= -1, z = e^x: S_j = sum_k (-1)^(k+1) z^(k-1) / k^(j+1) for j = 0, 1 (F_j = z S_j) and
// S_-1 = 1 / (1 + z) (F_-1 = z S_-1).
struct Series {
    double s_minus_one, s_zero, s_one;
};

Series series(double x) {
    const double z = std::exp(x);
    double s0 = 1.0, s1 = 1.0, power = 1.0, sign = 1.0;
    for (int k = 2; k <= 80; ++k) {
        power *= z;
        sign = -sign;
        const double kd = static_cast<double>(k);
        s0 += sign * power / kd;
        s1 += sign * power / (kd * kd);
        if (!(power > 1e-18)) break;
    }
    return {1.0 / (1.0 + z), s0, s1};
}

// |x| <= 1: F_0 = ln 2 + x / 2 + ln cosh(x / 2), F_1 = pi^2 / 12 + x ln 2 + x^2 / 4
// + 2 integral_0^(x/2) ln cosh, from the Taylor series of ln cosh (|x / 2| <= 1/2, a quarter of its
// radius pi / 2: the terms fall by about 0.1 each).
FermiOrdersZeroOne taylor(double x) {
    const auto& c = log_cosh_coefficients();
    const double y = 0.5 * x, y2 = y * y;
    double log_cosh = 0.0, integral = 0.0, power = y2;  // power = y^(2n)
    for (std::size_t k = 0; k < c.size(); ++k) {
        const double n = static_cast<double>(k + 1);
        log_cosh += c[k] * power;
        integral += c[k] * power * y / (2.0 * n + 1.0);
        power *= y2;
    }
    constexpr double ln2 = std::numbers::ln2;
    constexpr double pi2 = std::numbers::pi * std::numbers::pi;
    return {1.0 / (1.0 + std::exp(-x)), ln2 + 0.5 * x + log_cosh,
            pi2 / 12.0 + x * ln2 + 0.25 * x * x + 2.0 * integral};
}

}  // namespace

std::expected<void, base::Error> check_thermal_parameters(const ThermalParameters& t,
                                                          std::string_view owner) {
    const struct {
        std::string_view field;
        double value;
    } all[] = {{"conductivity_W_cmK", t.conductivity_W_cmK},
               {"conductivity_exponent", t.conductivity_exponent},
               {"heat_capacity_J_cm3K", t.heat_capacity_J_cm3K},
               {"thermopower_exponent_n", t.thermopower_exponent_n},
               {"thermopower_exponent_p", t.thermopower_exponent_p}};
    for (const auto& f : all) {
        if (!std::isfinite(f.value)) {
            return std::unexpected(thermal_error(owner, f.field, " must be finite", f.value));
        }
    }
    if (t.conductivity_W_cmK < 0.0) {
        return std::unexpected(thermal_error(owner, "conductivity_W_cmK", " must not be negative",
                                             t.conductivity_W_cmK));
    }
    if (has_thermal_data(t) && !(t.heat_capacity_J_cm3K > 0.0)) {
        return std::unexpected(thermal_error(owner, "heat_capacity_J_cm3K",
                                             " must be positive with a conductivity",
                                             t.heat_capacity_J_cm3K));
    }
    return {};
}

ValueSlope thermal_conductivity(const ThermalParameters& t, double temperature_K) {
    detail::expect_temperature(temperature_K);
    const double kappa =
        t.conductivity_W_cmK * std::pow(temperature_K / 300.0, -t.conductivity_exponent);
    return {kappa, -t.conductivity_exponent * kappa / temperature_K};
}

double band_gap_slope_eV_per_K(const Semiconductor& m, double temperature_K) {
    detail::expect_temperature(temperature_K);
    const SemiconductorParameters& p = m.parameters();
    const double T = temperature_K, s = T + p.varshni_beta_K;
    return -p.varshni_alpha_eV_per_K * T * (T + 2.0 * p.varshni_beta_K) / (s * s);
}

double intrinsic_density_log_slope(const Semiconductor& m, double temperature_K,
                                   double narrowing_eV) {
    const double T = temperature_K;
    const double kT = base::thermal_voltage(T);  // [eV]
    const double Eg = band_gap_eV(m, T);
    return 1.5 / T + (Eg - T * band_gap_slope_eV_per_K(m, T) - narrowing_eV) / (2.0 * kT * T);
}

double caughey_thomas_mobility_slope(const Semiconductor& m, Carrier carrier,
                                     double total_impurity_cm3, double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(total_impurity_cm3) && total_impurity_cm3 >= 0.0);
    detail::expect_temperature(temperature_K);
    const CaugheyThomasParameters& ct = carrier == Carrier::electron
                                            ? m.parameters().electron_mobility
                                            : m.parameters().hole_mobility;
    const double mu_max = detail::mu_max_at(ct, temperature_K);
    return ct.T_exponent * mu_max / temperature_K /
           (1.0 + std::pow(total_impurity_cm3 / ct.N_ref, ct.alpha));
}

double impact_ionization_temperature_factor_slope(double phonon_energy_eV,
                                                  double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(temperature_K) && temperature_K > 0.0);
    NITCAD_EXPECTS(std::isfinite(phonon_energy_eV) && phonon_energy_eV >= 0.0);
    if (phonon_energy_eV == 0.0) return 0.0;
    const double half = 0.5 * phonon_energy_eV / base::k_B_eV_per_K;
    const double s = std::sinh(half / temperature_K);
    return std::tanh(half / 300.0) * half / (temperature_K * temperature_K * s * s);
}

double emission_velocity_slope(const Semiconductor& m, Carrier carrier, double temperature_K) {
    return 0.5 * emission_velocity_cm_s(m, carrier, temperature_K) / temperature_K;
}

FermiOrdersZeroOne fermi_orders_zero_one(double x) noexcept {
    if (std::isnan(x)) return {x, x, x};
    if (x <= -1.0) {
        const Series s = series(x);
        const double z = std::exp(x);
        return {z * s.s_minus_one, z * s.s_zero, z * s.s_one};
    }
    if (x <= 1.0) return taylor(x);
    constexpr double pi2 = std::numbers::pi * std::numbers::pi;
    const FermiOrdersZeroOne r = fermi_orders_zero_one(-x);  // the series, -x < -1
    return {1.0 / (1.0 + std::exp(-x)), x + r.zero, 0.5 * x * x + pi2 / 6.0 - r.one};
}

ThermalDiffusionFactor fermi_dirac_thermal_diffusion(double x) noexcept {
    // Below -1 from the series ratios (e^x cancels); elsewhere from the integrals themselves.
    double ratio, product;  // F_1 / F_0 and F_1 F_-1 / F_0^2
    if (x <= -1.0) {
        const Series s = series(x);
        ratio = s.s_one / s.s_zero;
        product = ratio * s.s_minus_one / s.s_zero;
    } else {
        const FermiOrdersZeroOne f = fermi_orders_zero_one(x);
        ratio = f.one / f.zero;
        product = ratio * f.minus_one / f.zero;
    }
    return {2.0 * ratio - 1.5, 2.0 - 2.0 * product};
}

}  // namespace NiTCAD::physics
