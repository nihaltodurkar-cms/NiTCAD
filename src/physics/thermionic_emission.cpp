#include "NiTCAD/physics/thermionic_emission.hpp"

#include <cmath>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

double emission_velocity_cm_s(double dos_cm3, double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(dos_cm3) && dos_cm3 > 0.0);
    NITCAD_EXPECTS(std::isfinite(temperature_K) && temperature_K > 0.0);
    const double h = 2.0 * std::numbers::pi * base::hbar_J_s;  // [J s]
    const double kT = base::k_B_J_per_K * temperature_K;        // [J]
    return kT / h * std::cbrt(2.0 / dos_cm3);                  // [1/s] [cm]
}

double emission_velocity_cm_s(const Semiconductor& m, Carrier carrier, double temperature_K) {
    const double N = carrier == Carrier::electron ? conduction_band_dos(m, temperature_K)
                                                  : valence_band_dos(m, temperature_K);
    const double A = carrier == Carrier::electron ? m.parameters().richardson.electron
                                                  : m.parameters().richardson.hole;
    if (A > 0.0) return A * temperature_K * temperature_K / (base::q_C * N);  // [cm/s]
    return emission_velocity_cm_s(N, temperature_K);
}

}  // namespace NiTCAD::physics
