#include "NiTCAD/physics/statistics.hpp"

#include <cmath>
#include <limits>

#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/fermi_dirac.hpp"
#include "increasing_root.hpp"

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

DensityResult fermi_dirac_density(double n_ie, double log_dos, double eta) noexcept {
    const double x = eta - log_dos;
    if (eta > 700.0) {  // e^eta would overflow; N F_{1/2}(x) need not
        const double N = n_ie * std::exp(log_dos);
        const FermiIntegral f = fermi_half(x);
        return {N * f.value, N * f.derivative};
    }
    // n_ie e^eta gamma keeps the correctly rounded e^eta, as boltzmann_density does.
    const LogDegeneracy L = log_degeneracy(x);
    const double density = n_ie * std::exp(eta) * std::exp(L.value);
    return {density, density * (1.0 + L.d_eta)};
}

Degeneracy fermi_dirac_degeneracy(double n_ie, double log_dos, double density) noexcept {
    const double N = n_ie * std::exp(log_dos);
    const double v = density / N;
    if (v < 1e-6) {
        // e^x = v (1 + a v + b v^2 + O(v^3)) inverts F_{1/2}(x) = e^x - e^2x / 2^(3/2) + ...
        constexpr double a = 0.35355339059327376220;  // 2^(-3/2)
        constexpr double b = 0.05754991027012474516;  // 1/4 - 3^(-3/2)
        const double q = v * (a + b * v);
        return {-std::log1p(q), -(a + 2.0 * b * v) / ((1.0 + q) * N)};
    }
    if (!std::isfinite(v)) {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan, nan};
    }
    const LogDegeneracy L = log_degeneracy(inverse_fermi_half(v));
    return {L.value, L.d_eta / (density * (1.0 + L.d_eta))};
}

NeutralEquilibrium fermi_dirac_neutral_equilibrium(double net_doping, double n_ie,
                                                   double log_dos_n, double log_dos_p) {
    NITCAD_EXPECTS(std::isfinite(net_doping));
    NITCAD_EXPECTS(std::isfinite(n_ie) && n_ie > 0.0);
    NITCAD_EXPECTS(std::isfinite(log_dos_n) && std::isfinite(log_dos_p));
    // In the majority carrier's orientation y (eta for C >= 0, -eta for C < 0), with densities in
    // units of n_ie: ln major = y + ln gamma(y - g_major), minor = e^(-y + ln gamma(-y - g_minor)),
    // and the root of ln major - ln(c + minor), c = |C| / n_ie, increasing in y. ln(c + minor) is
    // formed as ln c + log1p(minor / c), so c itself never overflows. The majority is then
    // |C| + minor, so neutrality holds to rounding (the root is good to a few 1e-14 in y).
    const bool donors = net_doping >= 0.0;
    const double g_major = donors ? log_dos_n : log_dos_p;
    const double g_minor = donors ? log_dos_p : log_dos_n;
    const double log_c = net_doping != 0.0 ? std::log(std::abs(net_doping)) - std::log(n_ie)
                                           : 0.0;
    const auto residual = [&](double y) {
        const LogDegeneracy major = log_degeneracy(y - g_major);
        const LogDegeneracy minor = log_degeneracy(-y - g_minor);
        const double log_minor = -y + minor.value;
        const double d_log_minor = -(1.0 + minor.d_eta);
        double log_rhs = log_minor, d_log_rhs = d_log_minor;
        if (net_doping != 0.0) {
            const double r = std::exp(log_minor - log_c);  // minor / c
            log_rhs = log_c + std::log1p(r);
            d_log_rhs = r / (1.0 + r) * d_log_minor;
        }
        return detail::ValueSlope{y + major.value - log_rhs, 1.0 + major.d_eta - d_log_rhs};
    };
    const double guess = boltzmann_neutral_equilibrium(std::abs(net_doping), n_ie).eta;
    const double y = detail::increasing_root(residual, guess);
    const double minor = fermi_dirac_density(n_ie, g_minor, -y).density;
    const double major = std::abs(net_doping) + minor;
    if (donors) return {major, minor, y};
    return {minor, major, -y};
}

}  // namespace NiTCAD::physics
