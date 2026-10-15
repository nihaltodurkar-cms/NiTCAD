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

namespace {

// 1 / (1 + (N / C_r)^alpha), the factor of the mu_max term.
double masetti_lattice_share(const MasettiParameters& p, double N) {
    return 1.0 / (1.0 + std::pow(N / p.C_r, p.alpha));
}

}  // namespace

double masetti_mobility(const MasettiParameters& p, double total_impurity_cm3,
                        double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(total_impurity_cm3) && total_impurity_cm3 >= 0.0);
    detail::expect_temperature(temperature_K);
    const double N = total_impurity_cm3;
    const double mu_max = p.mu_max * std::pow(temperature_K / 300.0, p.T_exponent);
    if (N == 0.0) return mu_max - p.mu_min2;  // exp(-P_c / 0) and 1 / (1 + (C_s / 0)^beta) vanish
    return p.mu_min1 * std::exp(-p.P_c / N) + (mu_max - p.mu_min2) * masetti_lattice_share(p, N) -
           p.mu_1 / (1.0 + std::pow(p.C_s / N, p.beta));
}

double masetti_mobility_slope(const MasettiParameters& p, double total_impurity_cm3,
                              double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(total_impurity_cm3) && total_impurity_cm3 >= 0.0);
    detail::expect_temperature(temperature_K);
    const double mu_max = p.mu_max * std::pow(temperature_K / 300.0, p.T_exponent);
    return p.T_exponent * mu_max / temperature_K * masetti_lattice_share(p, total_impurity_cm3);
}

}  // namespace NiTCAD::physics
