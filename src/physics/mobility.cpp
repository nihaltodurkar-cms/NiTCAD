#include "NiTCAD/physics/mobility.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

double caughey_thomas_mobility(const Semiconductor& m, Carrier carrier,
                               double total_impurity_cm3, double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(total_impurity_cm3) && total_impurity_cm3 >= 0.0);
    NITCAD_EXPECTS(std::isfinite(temperature_K) && temperature_K > 0.0);
    const CaugheyThomasParameters& ct = carrier == Carrier::electron
                                            ? m.parameters().electron_mobility
                                            : m.parameters().hole_mobility;
    // Legacy operation order (materials.mobility_caughey_thomas).
    const double mu_max = ct.mu_max * std::pow(temperature_K / 300.0, ct.T_exponent);
    return ct.mu_min +
           (mu_max - ct.mu_min) / (1.0 + std::pow(total_impurity_cm3 / ct.N_ref, ct.alpha));
}

}  // namespace NiTCAD::physics
