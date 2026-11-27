#include "NiTCAD/physics/drift_velocity.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

DriftVelocity caughey_thomas_velocity(double mu0, double field_V_per_cm, double v_sat_cm_s, double beta) {
    NITCAD_EXPECTS(std::isfinite(mu0) && mu0 > 0.0);
    NITCAD_EXPECTS(std::isfinite(v_sat_cm_s) && v_sat_cm_s > 0.0);
    NITCAD_EXPECTS(std::isfinite(beta) && beta > 0.0);
    NITCAD_EXPECTS(std::isfinite(field_V_per_cm) && field_V_per_cm >= 0.0);
    const double s = 1.0 + std::pow(mu0 * field_V_per_cm / v_sat_cm_s, beta);
    const double factor = std::pow(s, -1.0 / beta);
    return {mu0 * field_V_per_cm * factor, mu0 * factor / s};
}

double piecewise_switch_field(double mu01, const HighFieldLaw& law) {
    const double mu02 = law.v_sat_cm_s / law.E_c_V_per_cm;
    NITCAD_EXPECTS(std::isfinite(mu01) && mu01 > 0.0 && mu01 < mu02);
    return law.E_c_V_per_cm * std::pow(std::pow(mu02 / mu01, law.beta) - 1.0, 1.0 / law.beta);
}

DriftVelocity piecewise_velocity(double mu01, double field_V_per_cm, const HighFieldLaw& law) {
    NITCAD_EXPECTS(std::isfinite(field_V_per_cm) && field_V_per_cm >= 0.0);
    if (field_V_per_cm < piecewise_switch_field(mu01, law)) return {mu01 * field_V_per_cm, mu01};
    return caughey_thomas_velocity(law.v_sat_cm_s / law.E_c_V_per_cm, field_V_per_cm, law.v_sat_cm_s, law.beta);
}

}  // namespace NiTCAD::physics
