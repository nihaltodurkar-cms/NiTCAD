#include "NiTCAD/physics/recombination.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

double scharfetter_lifetime(const Semiconductor& m, Carrier carrier, double total_impurity_cm3) {
    NITCAD_EXPECTS(std::isfinite(total_impurity_cm3) && total_impurity_cm3 >= 0.0);
    const ScharfetterLifetimeParameters& lt = m.parameters().lifetime;
    const double tau0 = carrier == Carrier::electron ? lt.tau_n0 : lt.tau_p0;
    return tau0 / (1.0 + total_impurity_cm3 / lt.N_ref);
}

}  // namespace NiTCAD::physics
