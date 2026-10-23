#include "NiTCAD/physics/surface_mobility.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"
#include "temperature.hpp"

namespace NiTCAD::physics {

SurfaceMobility lombardi_mobility(const LombardiParameters& p, double bulk_mobility,
                                  double total_impurity_cm3, double temperature_K,
                                  double field_V_per_cm) {
    NITCAD_EXPECTS(std::isfinite(bulk_mobility) && bulk_mobility > 0.0);
    NITCAD_EXPECTS(std::isfinite(total_impurity_cm3) && total_impurity_cm3 >= 0.0);
    detail::expect_temperature(temperature_K);
    NITCAD_EXPECTS(std::isfinite(field_V_per_cm) && field_V_per_cm > 0.0);
    const double E = field_V_per_cm, T = temperature_K;
    const double cube = std::cbrt(E);
    // mu_ac = B / E + k / (T E^(1/3)), k = C N^lambda (0 at N = 0, lambda > 0).
    const double k = total_impurity_cm3 == 0.0 ? 0.0 : p.C * std::pow(total_impurity_cm3, p.lambda);
    const double phonon = k / (T * cube);
    const double mu_ac = p.B / E + phonon;
    const double d_ac_dE = -p.B / (E * E) - phonon / (3.0 * E);
    const double d_ac_dT = -phonon / T;
    // mu_sr = delta / E^2.
    const double mu_sr = p.delta / (E * E);
    const double d_sr_dE = -2.0 * mu_sr / E;
    const double mu = 1.0 / (1.0 / bulk_mobility + 1.0 / mu_ac + 1.0 / mu_sr);
    // d mu / d mu_k = (mu / mu_k)^2 for each term of the sum.
    const double w_ac = (mu / mu_ac) * (mu / mu_ac);
    const double w_sr = (mu / mu_sr) * (mu / mu_sr);
    return {mu, w_ac * d_ac_dE + w_sr * d_sr_dE, w_ac * d_ac_dT,
            (mu / bulk_mobility) * (mu / bulk_mobility)};
}

}  // namespace NiTCAD::physics
