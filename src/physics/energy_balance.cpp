#include "NiTCAD/physics/energy_balance.hpp"

#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "increasing_root.hpp"

namespace NiTCAD::physics {

namespace {

// beta = 2: y = x^2 = a s / 2 and dy/da = s^2 / (2 r), r = sqrt(a^2 + 4), s = a + r (formed as
// 4 / (r - a) for a < 0, where a + r cancels).
struct Squared {
    double y, d_a;
};
Squared squared_root(double a) {
    const double r = std::sqrt(a * a + 4.0);
    const double s = a >= 0.0 ? a + r : 4.0 / (r - a);
    return {0.5 * a * s, s * s / (2.0 * r)};
}

// The exact root for a > 0: closed forms for beta = 1 and 2, Newton otherwise.
ReducedField exact_root(double a, double beta) {
    if (beta == 1.0) {
        // x^2 = a (1 + x): x = (a + sqrt(a^2 + 4 a)) / 2, dx/da = (1 + x) / sqrt(a^2 + 4 a).
        const double root = std::sqrt(a * a + 4.0 * a);
        const double x = 0.5 * (a + root);
        return {x, (1.0 + x) / root};
    }
    if (beta == 2.0) {
        const Squared q = squared_root(a);
        const double x = std::sqrt(q.y);
        return {x, q.d_a / (2.0 * x)};
    }
    return reduced_field_newton(a, beta);
}

}  // namespace

ReducedField reduced_field_newton(double a, double beta) {
    NITCAD_EXPECTS(std::isfinite(a) && a > 0.0);
    NITCAD_EXPECTS(std::isfinite(beta) && beta >= 1.0);
    // f(u) = 2 u - ln(1 + e^(beta u)) / beta - ln a with u = ln x; f' = 2 - s in [1, 2),
    // s = x^beta / (1 + x^beta). x ~ sqrt(a) for small a and ~ a for large a.
    const double target = std::log(a);
    const auto share = [beta](double u) {
        const double z = beta * u;
        return z > 0.0 ? 1.0 / (1.0 + std::exp(-z)) : std::exp(z) / (1.0 + std::exp(z));
    };
    const double u = detail::increasing_root(
        [&](double v) {
            const double z = beta * v;
            const double log1pe = z > 0.0 ? z + std::log1p(std::exp(-z)) : std::log1p(std::exp(z));
            return detail::ValueSlope{2.0 * v - log1pe / beta - target, 2.0 - share(v)};
        },
        a < 1.0 ? 0.5 * target : target);
    const double x = std::exp(u);
    return {x, x / (a * (2.0 - share(u)))};
}

ReducedField reduced_field(double a, double beta) {
    NITCAD_EXPECTS(std::isfinite(a));
    NITCAD_EXPECTS(std::isfinite(beta) && beta >= 1.0);
    if (a <= 0.0) return {0.0, 0.0};
    constexpr double a0 = energy_balance_regularization;
    if (a >= a0) return exact_root(a, beta);
    // The cubic c2 a^2 + c3 a^3 with the exact value and slope at a0: x and dx/da vanish at 0.
    const ReducedField e = exact_root(a0, beta);
    const double c2 = (3.0 * e.x - a0 * e.d_a) / (a0 * a0);
    const double c3 = (a0 * e.d_a - 2.0 * e.x) / (a0 * a0 * a0);
    return {a * a * (c2 + c3 * a), a * (2.0 * c2 + 3.0 * c3 * a)};
}

CarrierTemperatureMobility carrier_temperature_mobility(double mu0, double rise_K, double tau_w_s,
                                                        const CanaliParameters& c) {
    NITCAD_EXPECTS(std::isfinite(mu0) && mu0 > 0.0);
    NITCAD_EXPECTS(std::isfinite(rise_K));
    NITCAD_EXPECTS(std::isfinite(tau_w_s) && tau_w_s > 0.0);
    const double v2 = c.v_sat_cm_s * c.v_sat_cm_s;
    // a = mu0 P / v_sat^2 with P = (3/2) (k/q) rise / tau_w; a is linear in mu0 and in the rise.
    const double a_d_rise = 1.5 * (base::k_B_J_per_K / base::q_C) * mu0 / (tau_w_s * v2);
    const double a = a_d_rise * rise_K;
    const ReducedField r = reduced_field(a, c.beta);
    // The mobility factor m and dm/da.
    double m = 1.0, m_a = 0.0;
    if (c.beta == 2.0) {
        const Squared q = squared_root(a);
        m = 1.0 / std::sqrt(1.0 + q.y);
        m_a = -m * q.d_a / (2.0 * (1.0 + q.y));
    } else if (r.x > 0.0) {
        const double xb = std::pow(r.x, c.beta);
        const double x_power = c.beta == 1.0 ? 1.0 : std::pow(r.x, c.beta - 1.0);
        m = std::pow(1.0 + xb, -1.0 / c.beta);
        m_a = -m * x_power / (1.0 + xb) * r.d_a;
    }
    // mu = mu0 m(a), E = x v_sat / mu0; a / mu0 = da/dmu0.
    const double field = r.x * c.v_sat_cm_s / mu0;
    return {.mobility = mu0 * m,
            .d_mu0 = m + m_a * a,
            .d_rise = mu0 * m_a * a_d_rise,
            .field = field,
            .field_d_mu0 = c.v_sat_cm_s * (r.d_a * a - r.x) / (mu0 * mu0),
            .field_d_rise = c.v_sat_cm_s * r.d_a * a_d_rise / mu0};
}

}  // namespace NiTCAD::physics
