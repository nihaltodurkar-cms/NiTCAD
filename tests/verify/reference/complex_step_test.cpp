// TC8 (DECISIONS Unit 55 section 11.3; readings ruled in 11.13): L0 tests of complex-step differentiation in double.
//   - derivatives of e^x, sin x, x^{3/2} (as z sqrt z) and the Lyness-Moler composite e^x / sqrt(sin^3 x + cos^3 x)
//     at h = 1e-200 to 4 eps relative (eps = 2^-52) against the analytic derivative evaluated in double-double;
//     the composite at its published point x = 1.5 (Martins, Kroo and Alonso 2000, (15)) and 10^3 points of each
//     function on the ruled ranges; e^x and x^{3/2} (no zero of f') judged at every point, sin and the composite
//     where |f'| >= 0.01 max |f'| on the grid and printed elsewhere (rulings after measurement, 11.13);
//   - the Richardson truncation estimate (with F_s = 1.25) bounds the true error at h = 1e-3 without exceeding
//     1.25 (1 + 1e-3) times it, and is at rounding level at 1e-200;
//   - Table 5's new definitions (abs, max, min, sign, dim) give the branch derivatives, and std::abs does not;
//   - refusals, and the papers' pitfalls (arcsin by its log formula, sin by exponentials) shown to lose h.
// sin and cos in double-double are test-local: reduction by k pi/2 with the library's pi, then Taylor series on
// |r| <= pi/4, checked below by sin^2 + cos^2 = 1 and sin 2x = 2 sin x cos x.
#include "NiTCAD/verify/reference/complex_step.hpp"

#include "NiTCAD/verify/reference/double_double.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::complex_d;
using ref::DD;
using ref::u_dd;

namespace {

constexpr double eps = 0x1p-52;
constexpr double pi_d = std::numbers::pi;

DD sin_series(const DD& x) {
    const DD x2 = x * x;
    DD term = x, sum = x;
    for (int k = 1; k < 40; ++k) {
        term = -(term * x2) / DD((2.0 * k) * (2.0 * k + 1.0));
        sum += term;
        if (std::abs(term.hi) <= 0x1p-112 * std::abs(sum.hi)) break;
    }
    return sum;
}
DD cos_series(const DD& x) {
    const DD x2 = x * x;
    DD term(1.0), sum(1.0);
    for (int k = 1; k < 40; ++k) {
        term = -(term * x2) / DD((2.0 * k - 1.0) * (2.0 * k));
        sum += term;
        if (std::abs(term.hi) < 0x1p-112) break;
    }
    return sum;
}

// sin and cos of x by r = x - k pi/2, |r| <= pi/4 (|x| a few multiples of 2 pi here).
struct SinCos {
    DD s, c;
};
SinCos sincos_dd(const DD& x) {
    const DD half_pi = ldexp(ref::pi(), -1);
    const double k = std::nearbyint((x / half_pi).to_double());
    const DD r = x - DD(k) * half_pi;
    const DD s = sin_series(r), c = cos_series(r);
    const int q = static_cast<int>(((static_cast<long long>(k) % 4) + 4) % 4);
    switch (q) {
        case 0: return {s, c};
        case 1: return {c, -s};
        case 2: return {-s, -c};
        default: return {-c, s};
    }
}

struct Grid {
    std::vector<double> x;
};
Grid linear(double a, double b, int n) {
    Grid g;
    for (int i = 0; i < n; ++i) g.x.push_back(a + (b - a) * i / (n - 1));
    return g;
}
Grid logspaced(double e0, double e1, int n) {
    Grid g;
    for (int i = 0; i < n; ++i) g.x.push_back(std::pow(10.0, e0 + (e1 - e0) * i / (n - 1)));
    return g;
}

// The test functions in complex form, and their analytic derivatives in double-double.
complex_d f_exp(const complex_d& z) { return std::exp(z); }
complex_d f_sin(const complex_d& z) { return std::sin(z); }
complex_d f_x32(const complex_d& z) { return z * std::sqrt(z); }
complex_d f_x32_pow(const complex_d& z) { return std::pow(z, 1.5); }
complex_d f_lm(const complex_d& z) {
    const complex_d s = std::sin(z), c = std::cos(z);
    return std::exp(z) / std::sqrt(s * s * s + c * c * c);
}
DD d_exp(double x) { return ref::exp(DD(x)); }
DD d_sin(double x) { return sincos_dd(DD(x)).c; }
DD d_x32(double x) { return DD(1.5) * ref::sqrt(DD(x)); }
// ln f = x - (1/2) ln(s^3 + c^3), so f' = f (1 - (3/2) s c (s - c) / (s^3 + c^3)).
DD d_lm(double x) {
    const auto [s, c] = sincos_dd(DD(x));
    const DD q = s * s * s + c * c * c;
    const DD f = ref::exp(DD(x)) / ref::sqrt(q);
    return f * (DD(1.0) - DD(1.5) * s * c * (s - c) / q);
}

double rel_err_eps(double d, const DD& r) { return std::abs((DD(d) - r).to_double()) / std::abs(r.to_double()) / eps; }

struct Sweep {
    double worst_judged = 0.0;    // in eps, over points with |f'| >= 0.01 max |f'|
    double worst_printed = 0.0;   // in eps, over the rest
    double x_worst = 0.0;
    int judged = 0, printed = 0, not_ok = 0;
};

// f_prime_has_zero = false: f' has no zero on the range, and every point is judged.
template <class F, class D>
Sweep sweep(F f, D dref, const Grid& g, bool f_prime_has_zero) {
    std::vector<DD> r;
    double maxd = 0.0;
    for (double x : g.x) {
        r.push_back(dref(x));
        maxd = std::max(maxd, std::abs(r.back().to_double()));
    }
    Sweep s;
    for (std::size_t i = 0; i < g.x.size(); ++i) {
        const auto res = ref::complex_step(f, g.x[i]);
        if (res.status != ref::ComplexStepStatus::ok) {
            ++s.not_ok;
            continue;
        }
        const double e = rel_err_eps(res.derivative, r[i]);
        if (!f_prime_has_zero || std::abs(r[i].to_double()) >= 0.01 * maxd) {
            ++s.judged;
            if (e > s.worst_judged) {
                s.worst_judged = e;
                s.x_worst = g.x[i];
            }
        } else {
            ++s.printed;
            s.worst_printed = std::max(s.worst_printed, e);
        }
    }
    return s;
}

}  // namespace

TEST_CASE("TC8: the test-local double-double sin and cos (identities on [-2 pi, 2 pi])") {
    double worst_pyth = 0.0, worst_double = 0.0;
    for (double x : linear(-2.0 * pi_d, 2.0 * pi_d, 1001).x) {
        const auto [s, c] = sincos_dd(DD(x));
        worst_pyth = std::max(worst_pyth, std::abs((s * s + c * c - DD(1.0)).to_double()) / u_dd);
        const auto [s2, c2] = sincos_dd(ldexp(DD(x), 1));
        worst_double = std::max(worst_double, std::abs((s2 - ldexp(s * c, 1)).to_double()) / u_dd);
    }
    INFO("sin^2 + cos^2 - 1: " << worst_pyth << " u_dd; sin 2x - 2 sin x cos x: " << worst_double << " u_dd");
    CHECK(worst_pyth <= 4.0);
    CHECK(worst_double <= 4.0);
}

TEST_CASE("TC8: the Lyness-Moler composite at x = 1.5 to 4 eps at h = 1e-200") {
    const auto res = ref::complex_step(f_lm, 1.5);
    REQUIRE(res.status == ref::ComplexStepStatus::ok);
    const DD r = d_lm(1.5);
    const double e = rel_err_eps(res.derivative, r);
    const auto [s, c] = sincos_dd(DD(1.5));
    const DD fv = ref::exp(DD(1.5)) / ref::sqrt(s * s * s + c * c * c);
    const double ef = rel_err_eps(res.value, fv);
    INFO("f'(1.5): complex step " << res.derivative << ", reference " << r.to_double() << ", error " << e
                                  << " eps; value Re f(x + ih) error " << ef << " eps");
    CHECK(e <= 4.0);
    CHECK(ef <= 4.0);  // addition: the reported value is f(x) to O(h^2), judged at TC8's unit
}

TEST_CASE("TC8: e^x, sin x, x^{3/2} and the composite on their grids to 4 eps at h = 1e-200") {
    struct Case {
        const char* name;
        Sweep s;
    };
    const Case cases[] = {
        {"e^x on [-200, 200]", sweep(f_exp, d_exp, linear(-200.0, 200.0, 1000), false)},
        {"sin x on [-2 pi, 2 pi]", sweep(f_sin, d_sin, linear(-2.0 * pi_d, 2.0 * pi_d, 1000), true)},
        {"x^{3/2} (z sqrt z) on [1e-50, 1e50]", sweep(f_x32, d_x32, logspaced(-50.0, 50.0, 1000), false)},
        {"Lyness-Moler on [0, 1.5]", sweep(f_lm, d_lm, linear(0.0, 1.5, 1000), true)},
    };
    for (const auto& c : cases) {
        INFO(c.name << ": worst judged " << c.s.worst_judged << " eps at x = " << c.s.x_worst << " (" << c.s.judged
                    << " points); printed only: worst " << c.s.worst_printed << " eps (" << c.s.printed
                    << " points); status not ok: " << c.s.not_ok);
        CHECK(c.s.not_ok == 0);
        CHECK(c.s.judged >= 500);
        CHECK(c.s.worst_judged <= 4.0);
    }
    // pow(z, 1.5) is printed only (11.13: exp(1.5 log z) carries the condition number |1.5 ln x|).
    const Sweep p = sweep(f_x32_pow, d_x32, logspaced(-50.0, 50.0, 1000), false);
    UNSCOPED_INFO("printed only: x^{3/2} as pow(z, 1.5): worst " << p.worst_judged << " eps at x = " << p.x_worst);
    CHECK(p.not_ok == 0);
}

TEST_CASE("TC8: the Richardson truncation estimate bounds the true error at h = 1e-3, rounding level at 1e-200") {
    struct Case {
        const char* name;
        complex_d (*f)(const complex_d&);
        DD (*d)(double);
        double x;
    };
    const Case cases[] = {
        {"e^x at 1", f_exp, d_exp, 1.0},
        {"sin x at 1", f_sin, d_sin, 1.0},
        {"x^{3/2} at 2", f_x32, d_x32, 2.0},
        {"Lyness-Moler at 1.5", f_lm, d_lm, 1.5},
    };
    for (const auto& c : cases) {
        const auto big = ref::complex_step(c.f, c.x, 1e-3);
        const double true_err = std::abs((DD(big.derivative) - c.d(c.x)).to_double());
        const auto tiny = ref::complex_step(c.f, c.x);
        INFO(c.name << ": h = 1e-3 true " << true_err << ", estimate " << big.truncation_estimate << " (ratio "
                    << big.truncation_estimate / true_err << ")"
                    << "; h = 1e-200 estimate " << tiny.truncation_estimate / (eps * std::abs(tiny.derivative))
                    << " eps |d|");
        CHECK(big.status == ref::ComplexStepStatus::ok);
        CHECK(true_err <= big.truncation_estimate);
        CHECK(big.truncation_estimate <= 1.25 * (1.0 + 1e-3) * true_err);  // the ruled F_s, not the library's
        CHECK(tiny.truncation_estimate <= 4.0 * eps * std::abs(tiny.derivative));
    }
}

TEST_CASE("TC8: the truncation error is second order in h (14)") {
    const DD r = d_lm(1.5);
    double e[3];
    const double hs[3] = {0x1p-4, 0x1p-5, 0x1p-6};
    for (int i = 0; i < 3; ++i) e[i] = std::abs((DD(ref::complex_step(f_lm, 1.5, hs[i]).derivative) - r).to_double());
    const double p1 = std::log2(e[0] / e[1]), p2 = std::log2(e[1] / e[2]);
    INFO("observed orders " << p1 << ", " << p2);
    CHECK(std::abs(p1 - 2.0) <= 0.1);
    CHECK(std::abs(p2 - 2.0) <= 0.1);
}

TEST_CASE("TC8: invalid steps and points are refused; a non-finite f is reported") {
    const double bad_h[] = {0.0, -1e-200, NAN, INFINITY, 1e-310, DBL_MAX};
    for (double h : bad_h) {
        INFO("h = " << h);
        CHECK(ref::complex_step(f_exp, 1.0, h).status == ref::ComplexStepStatus::invalid_input);
    }
    CHECK(ref::complex_step(f_exp, NAN).status == ref::ComplexStepStatus::invalid_input);
    CHECK(ref::complex_step(f_exp, INFINITY).status == ref::ComplexStepStatus::invalid_input);
    CHECK(ref::complex_step(f_exp, 1.0, DBL_MIN).status == ref::ComplexStepStatus::ok);
    CHECK(ref::complex_step([](const complex_d&) { return complex_d(NAN, 0.0); }, 1.0).status ==
          ref::ComplexStepStatus::non_finite);
    // The papers' lower limit: when f' h underflows inside f the derivative is silently lost (here e^-700 * 1e-200).
    const auto lost = ref::complex_step(f_exp, -700.0);
    CHECK(lost.status == ref::ComplexStepStatus::ok);
    CHECK(lost.derivative == 0.0);
}

TEST_CASE("TC8: Table 5's new definitions give the branch derivatives; std::abs does not") {
    const auto d = [](auto f, double x) { return ref::complex_step(f, x).derivative; };
    const auto close = [](double a, double b) { return std::abs(a - b) <= 4.0 * eps * std::abs(b); };
    namespace cs = ref::cs;
    CHECK(close(d([](const complex_d& z) { return cs::abs(z); }, -2.0), -1.0));
    CHECK(close(d([](const complex_d& z) { return cs::abs(z); }, 2.0), 1.0));
    CHECK(close(d([](const complex_d& z) { return cs::abs(z); }, 0.0), 1.0));  // x >= 0 takes +z: right derivative
    CHECK(d([](const complex_d& z) { return complex_d(std::abs(z), 0.0); }, 2.0) == 0.0);  // the modulus loses it

    const auto mx = [](const complex_d& z) { return cs::max(z * z, 3.0 * z); };
    CHECK(close(d(mx, 1.0), 3.0));
    CHECK(close(d(mx, 4.0), 8.0));
    CHECK(close(d(mx, 3.0), 6.0));  // tie: z1
    const auto mn = [](const complex_d& z) { return cs::min(z * z, 3.0 * z); };
    CHECK(close(d(mn, 1.0), 2.0));
    CHECK(close(d(mn, 4.0), 3.0));
    CHECK(close(d(mn, 3.0), 6.0));  // tie: z1
    const auto sg = [](const complex_d& z) { return cs::sign(z * z, z - 1.0); };
    CHECK(close(d(sg, 2.0), 4.0));
    CHECK(close(d(sg, 0.5), -1.0));
    CHECK(close(d([](const complex_d& z) { return cs::sign(z, complex_d(1.0, 0.0)); }, -2.0), -1.0));
    const auto dm = [](const complex_d& z) { return cs::dim(z * z, 3.0 * z); };
    CHECK(close(d(dm, 4.0), 5.0));
    CHECK(d(dm, 1.0) == 0.0);
}

TEST_CASE("TC8: the papers' pitfalls lose the h term (Martins, Sturdza and Alonso 2001, (2)-(5))") {
    const complex_d i(0.0, 1.0);
    const double exact_asin = 1.0 / std::sqrt(0.75);
    const auto bad_asin = ref::complex_step([&](const complex_d& z) { return -i * std::log(i * z + std::sqrt(1.0 - z * z)); }, 0.5);
    const auto bad_sin = ref::complex_step([&](const complex_d& z) { return (std::exp(i * z) - std::exp(-i * z)) / (2.0 * i); }, 0.5);
    const auto std_asin = ref::complex_step([](const complex_d& z) { return std::asin(z); }, 0.5);
    UNSCOPED_INFO("printed only: std::asin at 0.5: " << std_asin.derivative << " against " << exact_asin);
    INFO("log-formula arcsin " << bad_asin.derivative << ", exponential sin " << bad_sin.derivative);
    CHECK(std::abs(bad_asin.derivative - exact_asin) > 0.5 * exact_asin);
    CHECK(std::abs(bad_sin.derivative - std::cos(0.5)) > 0.5 * std::cos(0.5));
}
