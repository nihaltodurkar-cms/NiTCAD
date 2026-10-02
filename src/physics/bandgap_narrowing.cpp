#include "NiTCAD/physics/bandgap_narrowing.hpp"

#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "temperature.hpp"

namespace NiTCAD::physics {

double bandgap_narrowing_eV(const Semiconductor& m, double total_impurity_cm3) {
    NITCAD_EXPECTS(std::isfinite(total_impurity_cm3) && total_impurity_cm3 >= 0.0);
    const SlotboomParameters& s = m.parameters().bandgap_narrowing;
    if (!(total_impurity_cm3 > s.N0)) return 0.0;
    const double x = std::log(total_impurity_cm3 / s.N0);
    return s.E0_eV * (x + std::sqrt(x * x + 0.5) - std::sqrt(0.5));
}

double effective_intrinsic_density(const Semiconductor& m, double total_impurity_cm3,
                                   double temperature_K) {
    detail::expect_temperature(temperature_K);
    const double narrowing = bandgap_narrowing_eV(m, total_impurity_cm3);
    return intrinsic_density(m, temperature_K) *
           std::exp(narrowing / (2.0 * base::thermal_voltage(temperature_K)));
}

}  // namespace NiTCAD::physics
