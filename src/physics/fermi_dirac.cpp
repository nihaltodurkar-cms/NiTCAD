#include "NiTCAD/physics/fermi_dirac.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

#include "NiTCAD/base/contract.hpp"
#include "increasing_root.hpp"

namespace NiTCAD::physics {

namespace {

constexpr double series_limit = -2.0;  // the series below, the table above
constexpr double table_limit = 40.0;   // the table below, Sommerfeld above
constexpr double step = 0.025;
constexpr auto table_size = static_cast<std::size_t>((table_limit - series_limit) / step + 0.5) + 1;

// ln gamma, its slope, and ln F_{1/2} = eta + ln gamma (each formed where it is exact).
struct Core {
    double log_f;
    double log_gamma;
    double d_log_gamma;
};

// eta <= -2, x = e^eta: F_{1/2} = x S1 and F_{-1/2} = x S0, with
// S1 = sum_k (-1)^(k+1) x^(k-1) / k^(3/2) and S0 the same over k^(1/2). Both start at 1, so the
// tails are summed on their own: ln gamma = log1p(S1 - 1), d ln gamma = (S0 - S1) / S1.
Core series(double eta) {
    const double x = std::exp(eta);
    double tail1 = 0.0, tail01 = 0.0;  // S1 - 1 and S0 - S1
    double power = 1.0, sign = 1.0;
    for (int k = 2; k <= 60; ++k) {
        power *= x;
        sign = -sign;
        const double root = std::sqrt(static_cast<double>(k));
        const double k32 = 1.0 / (k * root);
        tail1 += sign * power * k32;
        tail01 += sign * power * (1.0 / root - k32);
        if (!(power > 1e-17 * x)) break;  // also ends at once for x = 0 and on NaN
    }
    const double log_gamma = std::log1p(tail1);
    return {eta + log_gamma, log_gamma, tail01 / (1.0 + tail1)};
}

// eta > 40: F_{1/2} = A eta^(3/2) S, S = 1 + sum_k c_k u^k, u = eta^-2, with
// c_k = 2 (1 - 2^(1-2k)) zeta(2k) prod_{m=0}^{2k-1} (3/2 - m) (Sommerfeld; c_1 = pi^2 / 8,
// c_2 = 7 pi^4 / 640), the values of that closed form (zeta(2k) through the Bernoulli numbers);
// F_{1/2} above eta = 40 is checked against the double-double quadrature of the tests.
// F_{-1/2} = A eta^(1/2) (3/2 S - sum 2k c_k u^k).
Core sommerfeld(double eta) {
    static constexpr std::array<double, 10> c{
        1.233700550136169827354,   1.065411933184401657274,  9.701518554959126805969,
        242.715048146678318261,    11865.69174544562355668,  958843.3951441106180105,
        115801357.4749065758053,   19542370099.72680073539,  4392197855056.560499082,
        1268250756252705.624397};
    constexpr double A = 0.7522527780636750492641;  // 4 / (3 sqrt(pi))
    const double u = 1.0 / (eta * eta);
    double S = 0.0, D = 0.0;  // sum c_k u^k and sum 2k c_k u^k, by Horner
    for (std::size_t k = c.size(); k-- > 0;) {
        S = (S + c[k]) * u;
        D = (D + 2.0 * static_cast<double>(k + 1) * c[k]) * u;
    }
    S += 1.0;
    const double log_f = std::log(A) + 1.5 * std::log(eta) + std::log(S);
    return {log_f, log_f - eta, (1.5 * S - D) / (eta * S) - 1.0};
}

// The table: ln gamma and its first two derivatives at eta_i = -2 + i step, from Gauss-Legendre
// quadratures of F_{1/2}, F_{-1/2} and F_{-3/2} (d^2/deta^2 F_{1/2}):
//     F_j = (1 / Gamma(j + 1)) Integral_0^inf t^j s(t - eta) dt for j = 1/2, -1/2, and
//     F_{-3/2} = (1 / sqrt(pi)) Integral_0^inf t^(-1/2) s (1 - s) dt, with s(x) = 1 / (1 + e^x).
// t in [0, 1] is integrated in s = sqrt(t) (one panel; the substitution removes the t^(-1/2)
// endpoint singularity), and t in [1, max(eta, 0) + 40] in panels of width at most 2, so the
// Fermi edge (width ~1, poles at distance pi) is resolved for every eta. The neglected tail is
// below e^-40 relative. 32 nodes per panel: the quadrature error is far below rounding.
struct Table {
    std::vector<double> value, slope, curvature;
};

struct Rule {
    std::array<double, 32> x, w;  // on [-1, 1]
};

Rule gauss_legendre() {
    Rule r{};
    constexpr int n = 32;
    for (int i = 0; i < n / 2; ++i) {
        double z = std::cos(std::numbers::pi * (i + 0.75) / (n + 0.5));
        double derivative = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = z;  // Legendre recurrence up to P_n
            for (int k = 2; k <= n; ++k) {
                const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
                p0 = p1;
                p1 = p2;
            }
            derivative = n * (z * p1 - p0) / (z * z - 1.0);
            const double dz = p1 / derivative;
            z -= dz;
            if (std::abs(dz) <= 1e-16) break;
        }
        const double w = 2.0 / ((1.0 - z * z) * derivative * derivative);
        r.x[static_cast<std::size_t>(i)] = -z;
        r.x[static_cast<std::size_t>(n - 1 - i)] = z;
        r.w[static_cast<std::size_t>(i)] = w;
        r.w[static_cast<std::size_t>(n - 1 - i)] = w;
    }
    return r;
}

// Neumaier's compensated sum: a node's ~700 positive terms would otherwise leave a few 1e-15.
struct Sum {
    double sum = 0.0, carry = 0.0;
    void add(double x) {
        const double t = sum + x;
        carry += std::abs(sum) >= std::abs(x) ? (sum - t) + x : (x - t) + sum;
        sum = t;
    }
    [[nodiscard]] double value() const { return sum + carry; }
};

Table build_table() {
    const Rule rule = gauss_legendre();
    constexpr double rpi = std::numbers::inv_sqrtpi;
    // s(x) and s (1 - s) without overflow.
    const auto fermi = [](double x, double& s, double& s1s) {
        const double e = std::exp(-std::abs(x));
        const double d = 1.0 / (1.0 + e);
        s = x > 0.0 ? e * d : d;
        s1s = e * d * d;
    };
    Table t;
    t.value.resize(table_size);
    t.slope.resize(table_size);
    t.curvature.resize(table_size);
    for (std::size_t i = 0; i < table_size; ++i) {
        const double eta = series_limit + static_cast<double>(i) * step;
        Sum f12, fm12, fm32;
        double s = 0.0, s1s = 0.0;
        for (std::size_t q = 0; q < rule.x.size(); ++q) {  // t = sigma^2, sigma in [0, 1]
            const double sigma = 0.5 * (1.0 + rule.x[q]);
            const double w = 0.5 * rule.w[q];
            fermi(sigma * sigma - eta, s, s1s);
            f12.add(w * 4.0 * rpi * sigma * sigma * s);
            fm12.add(w * 2.0 * rpi * s);
            fm32.add(w * 2.0 * rpi * s1s);
        }
        const double t_max = std::max(eta, 0.0) + 40.0;
        const int panels = static_cast<int>(std::ceil((t_max - 1.0) / 2.0));
        const double width = (t_max - 1.0) / panels;
        for (int k = 0; k < panels; ++k) {
            const double mid = 1.0 + (k + 0.5) * width;
            for (std::size_t q = 0; q < rule.x.size(); ++q) {
                const double tt = mid + 0.5 * width * rule.x[q];
                const double w = 0.5 * width * rule.w[q];
                const double root = std::sqrt(tt);
                fermi(tt - eta, s, s1s);
                f12.add(w * 2.0 * rpi * root * s);
                fm12.add(w * rpi * s / root);
                fm32.add(w * rpi * s1s / root);
            }
        }
        const double f = f12.value();
        const double r1 = fm12.value() / f;
        t.value[i] = std::log(f) - eta;
        t.slope[i] = r1 - 1.0;
        t.curvature[i] = fm32.value() / f - r1 * r1;
    }
    return t;
}

const Table& table() {
    static const Table t = build_table();  // immutable once built; C++ makes the build thread-safe
    return t;
}

// -2 < eta <= 40: quintic Hermite interpolation on the cell holding eta.
Core interpolate(double eta) {
    const Table& t = table();
    const double position = (eta - series_limit) / step;
    const auto i = std::min(static_cast<std::size_t>(position), table_size - 2);
    const double u = (eta - (series_limit + static_cast<double>(i) * step)) / step;
    const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, u5 = u4 * u;
    const double h0 = 1.0 - 10.0 * u3 + 15.0 * u4 - 6.0 * u5;
    const double h1 = u - 6.0 * u3 + 8.0 * u4 - 3.0 * u5;
    const double h2 = 0.5 * (u2 - 3.0 * u3 + 3.0 * u4 - u5);
    const double h3 = 1.0 - h0;
    const double h4 = -4.0 * u3 + 7.0 * u4 - 3.0 * u5;
    const double h5 = 0.5 * (u3 - 2.0 * u4 + u5);
    const double d0 = -30.0 * u2 + 60.0 * u3 - 30.0 * u4;
    const double d1 = 1.0 - 18.0 * u2 + 32.0 * u3 - 15.0 * u4;
    const double d2 = 0.5 * (2.0 * u - 9.0 * u2 + 12.0 * u3 - 5.0 * u4);
    const double d4 = -12.0 * u2 + 28.0 * u3 - 15.0 * u4;
    const double d5 = 0.5 * (3.0 * u2 - 8.0 * u3 + 5.0 * u4);
    const double y0 = t.value[i], y1 = t.value[i + 1];
    const double s0 = step * t.slope[i], s1 = step * t.slope[i + 1];
    const double c0 = step * step * t.curvature[i], c1 = step * step * t.curvature[i + 1];
    const double log_gamma = h0 * y0 + h1 * s0 + h2 * c0 + h3 * y1 + h4 * s1 + h5 * c1;
    const double slope = (d0 * (y0 - y1) + d1 * s0 + d2 * c0 + d4 * s1 + d5 * c1) / step;
    return {eta + log_gamma, log_gamma, slope};
}

Core core(double eta) {
    if (eta <= series_limit) return series(eta);
    if (eta <= table_limit) return interpolate(eta);
    if (eta > table_limit) return sommerfeld(eta);
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    return {nan, nan, nan};
}

}  // namespace

FermiIntegral fermi_half(double eta) noexcept {
    const Core c = core(eta);
    // e^eta gamma keeps the correctly rounded e^eta (exp(eta + ln gamma) would lose |eta| ulp in
    // the sum); above the table e^eta may overflow where F_{1/2} does not.
    const double value =
        eta <= table_limit ? std::exp(eta) * std::exp(c.log_gamma) : std::exp(c.log_f);
    return {value, value * (1.0 + c.d_log_gamma)};
}

LogDegeneracy log_degeneracy(double eta) noexcept {
    const Core c = core(eta);
    return {c.log_gamma, c.d_log_gamma};
}

double inverse_fermi_half(double nu) {
    NITCAD_EXPECTS(std::isfinite(nu) && nu > 0.0);
    const double target = std::log(nu);
    // Below nu = 1 ln nu is within nu / 2^(3/2) of the root; above, the leading Sommerfeld term
    // (3 sqrt(pi) nu / 4)^(2/3) is within a few units of it.
    const double guess =
        nu <= 1.0 ? target : std::pow(0.75 * std::sqrt(std::numbers::pi) * nu, 2.0 / 3.0);
    return detail::increasing_root(
        [&](double eta) {
            const Core c = core(eta);
            return detail::ValueSlope{c.log_f - target, 1.0 + c.d_log_gamma};
        },
        guess);
}

}  // namespace NiTCAD::physics
