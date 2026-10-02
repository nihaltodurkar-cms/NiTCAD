#include "NiTCAD/physics/mobility.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"
#include "temperature.hpp"

namespace NiTCAD::physics {

double caughey_thomas_mobility(const Semiconductor& m, Carrier carrier,
                               double total_impurity_cm3, double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(total_impurity_cm3) && total_impurity_cm3 >= 0.0);
    detail::expect_temperature(temperature_K);
    const CaugheyThomasParameters& ct = carrier == Carrier::electron
                                            ? m.parameters().electron_mobility
                                            : m.parameters().hole_mobility;
    // Legacy operation order (materials.mobility_caughey_thomas).
    const double mu_max = detail::mu_max_at(ct, temperature_K);
    return ct.mu_min +
           (mu_max - ct.mu_min) / (1.0 + std::pow(total_impurity_cm3 / ct.N_ref, ct.alpha));
}

}  // namespace NiTCAD::physics
