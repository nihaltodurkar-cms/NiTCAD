#include "NiTCAD/physics/klaassen_mobility.hpp"

#include <cmath>

#include "NiTCAD/base/contract.hpp"
#include "temperature.hpp"

namespace NiTCAD::physics {

namespace {

// A value with its partials in n, p and T (forward mode): every operation below is exact
// differentiation of the recorded formulas, so the partials are analytic to rounding.
struct D {
    double v, n = 0.0, p = 0.0, T = 0.0;
};
D operator+(D a, D b) { return {a.v + b.v, a.n + b.n, a.p + b.p, a.T + b.T}; }
D operator-(D a, D b) { return {a.v - b.v, a.n - b.n, a.p - b.p, a.T - b.T}; }
D operator*(D a, D b) {
    return {a.v * b.v, a.n * b.v + a.v * b.n, a.p * b.v + a.v * b.p, a.T * b.v + a.v * b.T};
}
D operator*(double k, D a) { return {k * a.v, k * a.n, k * a.p, k * a.T}; }
D operator/(D a, D b) {
    const double q = a.v / b.v;
    return {q, (a.n - q * b.n) / b.v, (a.p - q * b.p) / b.v, (a.T - q * b.T) / b.v};
}
D operator/(double k, D b) { return D{k} / b; }
// a^e for a >= 0. At a = 0 with e > 0 the value is 0 and the partials are taken as 0: the
// header's non-differentiable point, and elsewhere a factor that multiplies a vanishing one.
D pow(D a, double e) {
    if (a.v == 0.0) return {e == 0.0 ? 1.0 : 0.0};
    const double y = std::pow(a.v, e);
    const double k = e * y / a.v;
    return {y, k * a.n, k * a.p, k * a.T};
}

const KlaassenCarrierParameters& own(const KlaassenParameters& K, Carrier c) {
    return c == Carrier::electron ? K.electron : K.hole;
}
const KlaassenCarrierParameters& other(const KlaassenParameters& K, Carrier c) {
    return c == Carrier::electron ? K.hole : K.electron;
}

// G's two scaled arguments: a = ((m0 / m)(T / 300))^s4, b = ((m / m0)(300 / T))^s7.
struct GScales {
    D a, b;
};
GScales g_scales(const KlaassenParameters& K, double m, D t300) {
    return {pow((1.0 / m) * t300, K.s[3]), pow(m * (1.0 / t300), K.s[6])};
}
D g_raw(const KlaassenParameters& K, const GScales& g, D P) {
    return D{1.0} - K.s[0] / pow(D{K.s[1]} + g.a * P, K.s[2]) + K.s[4] / pow(g.b * P, K.s[5]);
}

// dG/dP at P (values only): s1 s3 a (s2 + a P)^(-s3-1) - s5 s6 b^-s6 P^(-s6-1).
double g_slope(const KlaassenParameters& K, double a, double b, double P) {
    return K.s[0] * K.s[2] * a * std::pow(K.s[1] + a * P, -K.s[2] - 1.0) -
           K.s[4] * K.s[5] * std::pow(b, -K.s[5]) * std::pow(P, -K.s[5] - 1.0);
}

}  // namespace

KlaassenAtTemperature klaassen_at(const KlaassenParameters& K, Carrier carrier,
                                  double temperature_K) {
    detail::expect_temperature(temperature_K);
    const double m = own(K, carrier).mass_ratio;
    const GScales g = g_scales(K, m, D{temperature_K / 300.0});
    // G' < 0 at small P (the s5 term), > 0 at large P: bisect on ln P. 60 halvings of
    // [1e-8, 1e8] leave 3e-17 of ln P; G is flat at its minimum, so the held value is exact.
    double lo = std::log(1e-8), hi = std::log(1e8);
    NITCAD_EXPECTS(g_slope(K, g.a.v, g.b.v, std::exp(lo)) < 0.0 &&
                  g_slope(K, g.a.v, g.b.v, std::exp(hi)) > 0.0);
    for (int k = 0; k < 60; ++k) {
        const double mid = 0.5 * (lo + hi);
        (g_slope(K, g.a.v, g.b.v, std::exp(mid)) < 0.0 ? lo : hi) = mid;
    }
    return {&K, carrier, temperature_K, std::exp(0.5 * (lo + hi))};
}

KlaassenMobility klaassen_mobility(const KlaassenAtTemperature& at, double N_D, double N_A,
                                   double n, double p) {
    NITCAD_EXPECTS(std::isfinite(N_D) && N_D >= 0.0 && std::isfinite(N_A) && N_A >= 0.0);
    NITCAD_EXPECTS(std::isfinite(n) && n >= 0.0 && std::isfinite(p) && p >= 0.0);
    const KlaassenParameters& K = *at.parameters;
    const bool electron = at.carrier == Carrier::electron;
    const KlaassenCarrierParameters& c = own(K, at.carrier);
    const D t300{at.temperature_K / 300.0, 0.0, 0.0, 1.0 / 300.0};
    const D dn{n, 1.0, 0.0, 0.0}, dp{p, 0.0, 1.0, 0.0};

    const auto cluster = [](double N, double c_I, double N_ref) {
        if (N == 0.0) return 0.0;
        const double x = N_ref / N;
        return N * (1.0 + 1.0 / (c_I + x * x));
    };
    const double Ds = cluster(N_D, K.c_D, K.N_ref_D), As = cluster(N_A, K.c_A, K.N_ref_A);
    const double same = electron ? Ds : As, opposite = electron ? As : Ds;
    const D screen = electron ? dp : dn;
    const D N_sc = D{Ds + As} + screen;
    const D carriers = dn + dp;

    const D mu_L = c.mu_max * pow(t300, -c.theta);
    if (N_sc.v == 0.0 && carriers.v == 0.0) return {mu_L.v, 0.0, 0.0, mu_L.T};

    // P = [f_CW / P_CW + f_BH / P_BH]^-1; 1 / P is finite (0 only at the excluded point above).
    const D t2 = t300 * t300;
    const D inv_P = (K.f_CW / 3.97e13) * (pow(N_sc, 2.0 / 3.0) / t2) +
                    (K.f_BH / (1.36e20 * c.mass_ratio)) * (carriers / t2);
    const D P = 1.0 / inv_P;
    const double ratio = c.mass_ratio / other(K, at.carrier).mass_ratio;
    const D Pr = pow(P, K.r[5]);
    const D F = (K.r[0] * Pr + D{K.r[1] + K.r[2] * ratio}) / (Pr + D{K.r[3] + K.r[4] * ratio});
    // G held at its minimum below P_min: there dG/dP = 0, so G's T-partial is the explicit one.
    const GScales gs = g_scales(K, c.mass_ratio, t300);
    const D G = g_raw(K, gs, P.v < at.P_min ? D{at.P_min} : P);

    const D N_eff = D{same} + opposite * G + screen / F;
    const D mu_N = (c.mu_max * c.mu_max / (c.mu_max - c.mu_min)) *
                   pow(t300, 3.0 * c.alpha - 1.5);
    const D mu_c = (c.mu_max * c.mu_min / (c.mu_max - c.mu_min)) * pow(t300, -0.5);
    // 1 / mu_DAP = N_eff / [mu_N N_ref^alpha N_sc^(1 - alpha) + mu_c (n + p)].
    const D inv_dap = N_eff / (std::pow(c.N_ref, c.alpha) * (mu_N * pow(N_sc, 1.0 - c.alpha)) +
                               mu_c * carriers);
    const D mu = 1.0 / (1.0 / mu_L + inv_dap);
    return {mu.v, mu.n, mu.p, mu.T};
}

KlaassenMobility klaassen_mobility(const KlaassenParameters& K, Carrier carrier, double N_D,
                                   double N_A, double n, double p, double temperature_K) {
    return klaassen_mobility(klaassen_at(K, carrier, temperature_K), N_D, N_A, n, p);
}

EquilibriumDensities neutral_equilibrium(double N_D, double N_A, double n_i) {
    NITCAD_EXPECTS(std::isfinite(N_D) && N_D >= 0.0 && std::isfinite(N_A) && N_A >= 0.0);
    NITCAD_EXPECTS(std::isfinite(n_i) && n_i > 0.0);
    const double h = 0.5 * std::abs(N_D - N_A);
    const double root = std::hypot(h, n_i);
    const double major = h + root;  // the majority density (n_i when N_D = N_A)
    const double minor = n_i * n_i / major;
    const double d_major = n_i / root;
    const double d_minor = 2.0 * n_i / major - minor / major * d_major;
    return N_D >= N_A ? EquilibriumDensities{major, minor, d_major, d_minor}
                      : EquilibriumDensities{minor, major, d_minor, d_major};
}

KlaassenMobility klaassen_equilibrium_mobility(const KlaassenAtTemperature& at, double N_D,
                                               double N_A, double n_i, double dlnni_dT) {
    const EquilibriumDensities eq = neutral_equilibrium(N_D, N_A, n_i);
    KlaassenMobility m = klaassen_mobility(at, N_D, N_A, eq.n, eq.p);
    const double dni_dT = n_i * dlnni_dT;
    m.d_dT += (m.d_dn * eq.dn_dni + m.d_dp * eq.dp_dni) * dni_dT;
    return m;
}

}  // namespace NiTCAD::physics
