#include "NiTCAD/physics/impact_ionization.hpp"

#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

double impact_ionization_temperature_factor(double phonon_energy_eV, double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(temperature_K) && temperature_K > 0.0);
    NITCAD_EXPECTS(std::isfinite(phonon_energy_eV) && phonon_energy_eV >= 0.0);
    if (phonon_energy_eV == 0.0) return 1.0;
    const double half = 0.5 * phonon_energy_eV / base::k_B_eV_per_K;
    return std::tanh(half / 300.0) / std::tanh(half / temperature_K);
}

ImpactIonizationRate impact_ionization_coefficient(const ImpactIonizationCoefficients& c,
                                                   double gamma, double field_V_per_cm) {
    NITCAD_EXPECTS(std::isfinite(field_V_per_cm) && field_V_per_cm >= 0.0);
    NITCAD_EXPECTS(std::isfinite(gamma) && gamma > 0.0);
    const bool high = field_V_per_cm >= c.switch_V_per_cm;
    const double A = high ? c.A_high_per_cm : c.A_low_per_cm;
    const double B = high ? c.B_high_V_per_cm : c.B_low_V_per_cm;
    if (A == 0.0 || field_V_per_cm == 0.0) return {0.0, 0.0};
    const double x = gamma * B / field_V_per_cm;
    if (x > 700.0) return {0.0, 0.0};  // exp(-x) underflows; so would the derivative
    const double alpha = gamma * A * std::exp(-x);
    return {alpha, alpha * x / field_V_per_cm};
}

}  // namespace NiTCAD::physics
