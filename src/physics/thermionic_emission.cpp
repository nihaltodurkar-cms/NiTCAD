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

}  // namespace NiTCAD::physics
