#include "NiTCAD/physics/ionization.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"
#include "increasing_root.hpp"

namespace NiTCAD::physics {

IonizedDensity ionized_density(double N, double eta, double depth_over_kT,
                               double degeneracy) noexcept {
    if (!(depth_over_kT > 0.0)) return {N, 0.0};
    // s = 1 / (1 + e^y), y = eta + depth + ln g, with e = e^-|y| so nothing overflows:
    // s = e / (1 + e) for y > 0, 1 / (1 + e) otherwise, and s (1 - s) = e / (1 + e)^2 either way.
    const double y = eta + depth_over_kT + std::log(degeneracy);
    const double e = std::exp(-std::abs(y));
    const double d = 1.0 / (1.0 + e);
    const double s = y > 0.0 ? e * d : d;
    return {N * s, -N * e * d * d};
}

NeutralEquilibrium ionized_neutral_equilibrium(double donors, double acceptors, double n_ie,
                                               double log_dos_n, double log_dos_p,
                                               const DopantLevels& levels, bool fermi_dirac) {
    NITCAD_EXPECTS(std::isfinite(donors) && donors >= 0.0);
    NITCAD_EXPECTS(std::isfinite(acceptors) && acceptors >= 0.0);
    NITCAD_EXPECTS(std::isfinite(n_ie) && n_ie > 0.0);
    NITCAD_EXPECTS(std::isfinite(log_dos_n) && std::isfinite(log_dos_p));
    const auto density = [&](double g, double eta) {
        return fermi_dirac ? fermi_dirac_density(n_ie, g, eta) : boltzmann_density(n_ie, eta);
    };
    // f(eta) = n - p - N_D+ + N_A-, increasing in eta; scaled to order one.
    const double scale = donors + acceptors + n_ie;
    const auto charge = [&](double eta) {
        const DensityResult n = density(log_dos_n, eta);
        const DensityResult p = density(log_dos_p, -eta);
        const IonizedDensity nd = ionized_density(donors, eta - log_dos_n, levels.donor_kT,
                                                  levels.donor_degeneracy);
        const IonizedDensity na = ionized_density(acceptors, -eta - log_dos_p, levels.acceptor_kT,
                                                  levels.acceptor_degeneracy);
        return detail::ValueSlope{(n.density - p.density - nd.value + na.value) / scale,
                                  (n.d_eta + p.d_eta - nd.d_eta - na.d_eta) / scale};
    };
    const double C = donors - acceptors;
    const double guess = fermi_dirac
                             ? fermi_dirac_neutral_equilibrium(C, n_ie, log_dos_n, log_dos_p).eta
                             : boltzmann_neutral_equilibrium(C, n_ie).eta;
    const double eta = detail::increasing_root(charge, guess);
    return {density(log_dos_n, eta).density, density(log_dos_p, -eta).density, eta};
}

}  // namespace NiTCAD::physics
