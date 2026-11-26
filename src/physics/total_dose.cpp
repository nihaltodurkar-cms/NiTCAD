#include "NiTCAD/physics/total_dose.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

double charge_yield(double field_V_cm) {
    NITCAD_EXPECTS(std::isfinite(field_V_cm) && field_V_cm >= 0.0);
    return std::pow((field_V_cm + dose_yield_e0_V_cm) / (field_V_cm + dose_yield_e1_V_cm), dose_yield_exponent);
}

double oxide_trapped_charge_cm2(double dose_rad, double field_V_cm, double thickness_cm, double trapping_fraction) {
    NITCAD_EXPECTS(std::isfinite(dose_rad) && dose_rad > 0.0);
    NITCAD_EXPECTS(std::isfinite(thickness_cm) && thickness_cm > 0.0);
    NITCAD_EXPECTS(trapping_fraction >= 0.0 && trapping_fraction <= 1.0);
    return trapping_fraction * charge_yield(field_V_cm) * dose_pair_density_cm3_rad * dose_rad * thickness_cm;
}

}  // namespace NiTCAD::physics
