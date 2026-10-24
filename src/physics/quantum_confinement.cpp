#include "NiTCAD/physics/quantum_confinement.hpp"

#include <cmath>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "temperature.hpp"

namespace NiTCAD::physics {

double mlda_wavelength_cm(double mass_m0, double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(mass_m0) && mass_m0 > 0.0);
    detail::expect_temperature(temperature_K);
    const double m = mass_m0 * base::m0_kg;
    return base::hbar_J_s / std::sqrt(2.0 * m * base::k_B_J_per_K * temperature_K) * 100.0;
}

double mlda_factor(double distance_cm, double wavelength_cm) {
    NITCAD_EXPECTS(std::isfinite(distance_cm) && distance_cm >= 0.0);
    NITCAD_EXPECTS(std::isfinite(wavelength_cm) && wavelength_cm > 0.0);
    const double x = distance_cm / wavelength_cm;
    return -std::expm1(-x * x);
}

double mlda_mean_factor(double depth_cm, double wavelength_cm) {
    NITCAD_EXPECTS(std::isfinite(depth_cm) && depth_cm > 0.0);
    NITCAD_EXPECTS(std::isfinite(wavelength_cm) && wavelength_cm > 0.0);
    const double x = depth_cm / wavelength_cm;
    if (x >= 0.5) return 1.0 - std::sqrt(std::numbers::pi) / (2.0 * x) * std::erf(x);
    // (1/x) int_0^x (1 - e^(-u^2)) du = sum_k>=1 (-1)^(k+1) x^(2k) / (k! (2k + 1)); at x < 0.5 the
    // terms fall by at least x^2 / 2 per step, so 20 of them reach rounding.
    const double x2 = x * x;
    double term = 1.0;  // (-1)^(k+1) x^(2k) / k!
    double sum = 0.0;
    for (int k = 1; k <= 20; ++k) {
        term *= (k == 1 ? x2 : -x2 / k);
        sum += term / (2 * k + 1);
    }
    return sum;
}

}  // namespace NiTCAD::physics
