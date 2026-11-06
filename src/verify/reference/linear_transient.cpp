#include "NiTCAD/verify/reference/linear_transient.hpp"

#include <algorithm>
#include <cmath>

#include "NiTCAD/verify/reference/dopri5.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

}  // namespace

std::expected<DD, std::string> dielectric_relaxation_time(double eps_F_per_cm, double q_C, double mobility_cm2_Vs,
                                                          double density_cm3) {
    if (!positive_finite(eps_F_per_cm) || !positive_finite(q_C) || !positive_finite(mobility_cm2_Vs) ||
        !positive_finite(density_cm3)) {
        return std::unexpected("dielectric_relaxation_time: an input is not finite and positive");
    }
    return DD(eps_F_per_cm) / (DD(q_C) * DD(mobility_cm2_Vs) * DD(density_cm3));
}

std::expected<double, std::string> relaxation_second_method(const DD& tau_d, double t_end_s) {
    if (!(tau_d.hi > 0.0) || !positive_finite(t_end_s)) {
        return std::unexpected("relaxation_second_method: tau_d and t_end must be positive");
    }
    const DD rate = DD(1.0) / tau_d;
    double worst = 0.0;
    const auto f = [&](const DD&, const DDVector& y, DDVector& dy) { dy[0] = -rate * y[0]; };
    const long long n = 4000;
    const DD h = DD(t_end_s) / DD(static_cast<double>(n));
    DDVector y{DD(1.0)};
    long long k = 0;
    const auto out = [&](const Dopri5Dense&) {};
    for (k = 0; k < n; ++k) {
        y = dopri5_fixed(f, DD(static_cast<double>(k)) * h, y, DD(static_cast<double>(k + 1)) * h, 1, out);
        const DD exact = exp(-(DD(static_cast<double>(k + 1)) * h) * rate);
        worst = std::max(worst, std::abs((y[0] - exact).to_double()));
    }
    return worst;
}

std::expected<PulseMoments, std::string> haynes_shockley_moments(const PulseInput& in, double t_s) {
    if (!positive_finite(in.diffusivity_cm2_s) || !std::isfinite(in.velocity_cm_s) || !(in.lifetime_s > 0.0) ||
        !std::isfinite(in.centre_cm) || !positive_finite(in.sigma_cm) || !positive_finite(in.rate_cm3_s) ||
        !std::isfinite(in.t_on_s) || !std::isfinite(in.t_off_s) || !(in.t_on_s < in.t_off_s) || !std::isfinite(t_s)) {
        return std::unexpected("haynes_shockley_moments: an input is not finite, or out of range");
    }
    if (t_s < in.t_off_s) return std::unexpected("haynes_shockley_moments: t before the source ends");
    const DD a1 = DD(t_s) - DD(in.t_off_s), a2 = DD(t_s) - DD(in.t_on_s), span = a2 - a1;
    const DD lambda = std::isfinite(in.lifetime_s) ? DD(1.0) / DD(in.lifetime_s) : DD(0.0);
    const DD x = lambda * span;
    // The truncated exponential's moments of s = a - a1 on [0, span].
    DD m0, m1, var;  // int e^{-lambda s} ds, <s>, Var(s)
    if (x.hi == 0.0) {
        m0 = span;
        m1 = ldexp(span, -1);
        var = span * span / DD(12.0);
    } else {
        const DD em = expm1(x);  // e^x - 1
        m0 = span * (-expm1(-x)) / x;
        m1 = span * (DD(1.0) / x - DD(1.0) / em);
        var = span * span * (DD(1.0) / (x * x) - exp(x) / (em * em));
    }
    const DD w1 = exp(-lambda * a1);
    const DD amp = DD(in.rate_cm3_s) * DD(in.sigma_cm) * sqrt(ldexp(pi(), 1));
    const DD mean_a = a1 + m1;
    const DD v(in.velocity_cm_s), s(in.sigma_cm), D(in.diffusivity_cm2_s);
    PulseMoments p;
    p.area_cm2 = amp * w1 * m0;
    p.centre_cm = DD(in.centre_cm) + v * mean_a;
    p.width_cm = sqrt(s * s + ldexp(D, 1) * mean_a + v * v * var);
    // The second method: the age integrals by Gauss-Kronrod.
    const auto w = [&](const DD& a) { return exp(-lambda * a); };
    const AdaptiveResult I0 = qag15(w, a1, a2, 0.0, 1e-26);
    const AdaptiveResult I1 = qag15([&](const DD& a) { return a * w(a); }, a1, a2, 0.0, 1e-26);
    const AdaptiveResult I2 = qag15([&](const DD& a) { return a * a * w(a); }, a1, a2, 0.0, 1e-26);
    for (const AdaptiveResult* r : {&I0, &I1, &I2}) {
        if (r->status != QuadratureStatus::ok && r->status != QuadratureStatus::roundoff) {
            return std::unexpected("haynes_shockley_moments: a quadrature failed");
        }
    }
    const DD ma = I1.value / I0.value, va = I2.value / I0.value - ma * ma;
    const DD area2 = amp * I0.value, centre2 = DD(in.centre_cm) + v * ma;
    const DD width2 = sqrt(s * s + ldexp(D, 1) * ma + v * v * va);
    p.area_difference = std::abs(((area2 - p.area_cm2) / p.area_cm2).to_double());
    const double drift = std::abs((v * a2).to_double());
    p.centre_difference = drift > 0.0 ? std::abs((centre2 - p.centre_cm).to_double()) / drift : 0.0;
    p.width_difference = std::abs(((width2 - p.width_cm) / p.width_cm).to_double());
    return p;
}

}  // namespace NiTCAD::verify::reference
