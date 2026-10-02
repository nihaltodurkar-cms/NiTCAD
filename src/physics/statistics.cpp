#include "NiTCAD/physics/statistics.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

NeutralEquilibrium boltzmann_neutral_equilibrium(double net_doping, double n_ie) {
    NITCAD_EXPECTS(std::isfinite(net_doping));
    NITCAD_EXPECTS(std::isfinite(n_ie) && n_ie > 0.0);
    const double C = net_doping;
    // No intermediate overflows or underflows where the result is representable:
    // - the majority is 0.5 |C| + hypot(0.5 C, n_ie), a sum of two terms no larger than itself;
    //   neither C^2, n_ie^2 nor 2 n_ie is formed (the legacy 0.5 (|C| + sqrt(C^2 + 4 n_ie^2))
    //   overflows for |C| near DBL_MAX);
    // - the minority is n_ie (n_ie / majority), not n_ie^2 / majority, since n_ie^2 underflows
    //   for n_ie below about 1.5e-154 (cryogenic temperatures, or concentrations divided by Ns);
    // - eta = asinh(C / (2 n_ie)), with the quotient as 0.5 (C / n_ie): equal to C / (2 n_ie)
    //   except at overflow or underflow. On overflow |eta| = ln|C| - ln n_ie, which is asinh to
    //   double precision for any quotient above 1e8.
    const double majority = 0.5 * std::abs(C) + std::hypot(0.5 * C, n_ie);
    const double minority = n_ie * (n_ie / majority);
    const double ratio = 0.5 * (C / n_ie);
    const double eta = std::isfinite(ratio)
                           ? std::asinh(ratio)
                           : std::copysign(std::log(std::abs(C)) - std::log(n_ie), C);
    if (C >= 0.0) return {majority, minority, eta};
    return {minority, majority, eta};
}

}  // namespace NiTCAD::physics
