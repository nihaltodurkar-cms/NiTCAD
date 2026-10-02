// Private to the physics layer: the one temperature precondition and the one temperature
// scaling of the Caughey-Thomas mu_max, shared by the model functions and check_temperature.
#pragma once

#include <cmath>

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

}  // namespace NiTCAD::physics::detail
