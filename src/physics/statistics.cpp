#include "NiTCAD/physics/statistics.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

NeutralEquilibrium boltzmann_neutral_equilibrium(double net_doping, double n_ie) {
    NITCAD_EXPECTS(std::isfinite(net_doping));
    NITCAD_EXPECTS(std::isfinite(n_ie) && n_ie > 0.0);
    const double C = net_doping;
    // hypot avoids overflow of C^2 for any finite C.
    const double root = std::hypot(C, 2.0 * n_ie);
    const double n2 = n_ie * n_ie;
    double n = 0.0;
    double p = 0.0;
    if (C >= 0.0) {
        n = 0.5 * (C + root);
        p = n2 / n;
    } else {
        p = 0.5 * (-C + root);
        n = n2 / p;
    }
    return {n, p, std::asinh(C / (2.0 * n_ie))};
}

}  // namespace NiTCAD::physics
