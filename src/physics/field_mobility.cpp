#include "NiTCAD/physics/field_mobility.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

FieldMobility canali_mobility(double mu0, double field_V_per_cm, const CanaliParameters& c) {
    NITCAD_EXPECTS(std::isfinite(mu0) && mu0 > 0.0);
    NITCAD_EXPECTS(std::isfinite(field_V_per_cm) && field_V_per_cm >= 0.0);
    // Legacy operation order (materials.mobility_field).
    const double x = mu0 * field_V_per_cm / c.v_sat_cm_s;
    const double s = 1.0 + std::pow(x, c.beta);
    const double mu = mu0 / std::pow(s, 1.0 / c.beta);
    // dmu/dE = -mu x^(beta - 1) (mu0 / v_sat) / s; x^(beta - 1) is 1 at x = 0 when beta = 1.
    const double x_power = c.beta == 1.0 ? 1.0 : std::pow(x, c.beta - 1.0);
    return {mu, -mu * x_power * (mu0 / c.v_sat_cm_s) / s};
}

}  // namespace NiTCAD::physics
