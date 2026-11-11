// TVR-3 step 2 (DECISIONS Unit 55 section 15.11, frozen before code): L0 tests of the double-double special
// functions. Three kinds of evidence:
//   - printed values: erf(1) (OEIS A099286) and tanh(1) (OEIS A073744), 40 digits typed from the DATA fields
//     (read 2026-11-06), compared exactly by TC1's big-integer method;
//   - exact special values from pi and ln 2 (DLMF 25.12.4, 25.12.6 at -1 and 1/2; Li2(1) = zeta(2));
//   - a second method: adaptive Gauss-Kronrod (TC3) of each function's defining integral (DLMF 7.2.1 for erf
//     and erfc, 4.24.1 for atan, 25.12.2 for Li2), with the quadrature's own error estimate as U_ref;
// plus the identities of 15.11. Bounds are 15.11's; none is widened after a run (TVM2).
#include "NiTCAD/verify/reference/special_functions.hpp"

#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

#include "exact_integer.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using ref::u_dd;

namespace {

// |a - b| / |b| in u_dd (b != 0), evaluated in double-double.
double rel_udd(const DD& a, const DD& b) { return std::abs(((a - b) / b).hi) / u_dd; }
double abs_udd(const DD& a, const DD& b) { return std::abs((a - b).hi) / u_dd; }

// n points from lo to hi, both included.
std::vector<double> inclusive(double lo, double hi, int n) {
    std::vector<double> v;
    for (int i = 0; i < n; ++i) v.push_back(lo + (hi - lo) * i / (n - 1));
    return v;
}

std::vector<double> linear(double lo, double hi, int n) {
    std::vector<double> v;
    for (int i = 0; i < n; ++i) v.push_back(lo + (hi - lo) * (i + 1) / n);  // excludes lo, includes hi
    return v;
}

std::vector<double> log_spaced(double lo_exp, double hi_exp, int n) {
    std::vector<double> v;
    for (int i = 0; i < n; ++i) v.push_back(std::pow(10.0, lo_exp + (hi_exp - lo_exp) * i / (n - 1)));
    return v;
}

constexpr double quad_epsrel = 0.5e-28;  // the smallest relative tolerance qag15 accepts

struct QuadCheck {
    double worst_excess_udd = -1e300;  // max over points of |f - I| - (16 u_dd |I| + U_ref), in u_dd
    double worst_rel_udd = 0.0;         // the measured |f - I| / |I|, for the record
    int not_ok = 0;                     // quadrature status other than ok
};

template <class Integrand, class Function, class Limits>
QuadCheck quadrature_check(const std::vector<double>& xs, Integrand&& integrand, Function&& function,
                           Limits&& limits, double sign = 1.0) {
    QuadCheck q;
    for (const double x : xs) {
        const auto [a, b] = limits(x);
        const auto r = ref::qag15(integrand, a, b, 0.0, quad_epsrel);
        if (r.status != ref::QuadratureStatus::ok) ++q.not_ok;
        const DD I = r.value * DD(sign);
        const DD f = function(DD(x));
        const double diff = std::abs((f - I).hi);
        const double bound = 16.0 * u_dd * std::abs(I.hi) + r.abs_error;
        q.worst_excess_udd = std::max(q.worst_excess_udd, (diff - bound) / u_dd);
        q.worst_rel_udd = std::max(q.worst_rel_udd, diff / std::abs(I.hi) / u_dd);
    }
    return q;
}

}  // namespace

TEST_CASE("TVR-3 step 2: erf(1) and tanh(1) within 4 u_dd of OEIS's digits") {
    const auto e = exact::compare_with_literal(ref::erf(DD(1.0)), "0.8427007929497148693412206350826092592960", 4);
    INFO("erf(1): " << e.measured_udd << " u_dd");
    CHECK(e.within);
    const auto t = exact::compare_with_literal(ref::tanh(DD(1.0)), "0.7615941559557648881194582826047935904127", 4);
    INFO("tanh(1): " << t.measured_udd << " u_dd");
    CHECK(t.within);
}

TEST_CASE("TVR-3 step 2: exact special values to 8 u_dd") {
    const DD pi = ref::pi();
    const DD pi2 = pi * pi;
    const DD l2 = ref::ln2();
    struct Case {
        const char* name;
        DD value;
        DD exact;
    };
    const Case cases[] = {
        {"li2(-1) = -pi^2/12", ref::li2(DD(-1.0)), -(pi2 / DD(12.0))},
        {"li2(1/2) = pi^2/12 - ln^2 2 / 2", ref::li2(DD(0.5)), pi2 / DD(12.0) - ref::ldexp(l2 * l2, -1)},
        {"li2(1) = pi^2/6", ref::li2(DD(1.0)), pi2 / DD(6.0)},
        {"atan(1) = pi/4", ref::atan(DD(1.0)), ref::ldexp(pi, -2)},
        {"asin(1/2) = pi/6", ref::asin(DD(0.5)), pi / DD(6.0)},
        {"sin(pi/6) = 1/2", ref::sin(pi / DD(6.0)), DD(0.5)},
        {"cos(pi/4) = sqrt2/2", ref::cos(ref::ldexp(pi, -2)), ref::ldexp(ref::sqrt2(), -1)},
    };
    for (const auto& c : cases) {
        const double r = rel_udd(c.value, c.exact);
        INFO(c.name << ": " << r << " u_dd");
        CHECK(r <= 8.0);
    }
}

TEST_CASE("TVR-3 step 2: erf, erfc, atan and li2 against Gauss-Kronrod of their integrals") {
    const DD two_over_sqrt_pi = DD(2.0) / ref::sqrt(ref::pi());
    const auto gauss = [&](const DD& t) { return two_over_sqrt_pi * ref::exp(-(t * t)); };
    struct Named {
        const char* name;
        QuadCheck q;
    };
    const Named checks[] = {
        {"erf on (0, 6]",
         quadrature_check(linear(0.0, 6.0, 64), gauss, [](const DD& x) { return ref::erf(x); },
                          [](double x) { return std::pair{DD(0.0), DD(x)}; })},
        {"erfc on [3, 25]",
         quadrature_check(inclusive(3.0, 25.0, 64), gauss, [](const DD& x) { return ref::erfc(x); },
                          [](double x) { return std::pair{DD(x), DD(x + 15.0)}; })},
        {"atan on (0, 20]",
         quadrature_check(linear(0.0, 20.0, 64), [](const DD& t) { return DD(1.0) / (DD(1.0) + t * t); },
                          [](const DD& x) { return ref::atan(x); },
                          [](double x) { return std::pair{DD(0.0), DD(x)}; })},
        {"li2 on [-20, 0.99]",
         quadrature_check(inclusive(-20.0, 0.99, 64),
                          [](const DD& t) { return t.hi == 0.0 ? DD(-1.0) : ref::log1p(-t) / t; },
                          [](const DD& x) { return ref::li2(x); },
                          [](double x) { return std::pair{DD(0.0), DD(x)}; }, -1.0)},
    };
    for (const auto& c : checks) {
        INFO(c.name << ": worst |f - I| / |I| = " << c.q.worst_rel_udd << " u_dd, worst excess over the bound "
                    << c.q.worst_excess_udd << " u_dd, quadrature not ok at " << c.q.not_ok << " points");
        CHECK(c.q.not_ok == 0);
        CHECK(c.q.worst_excess_udd <= 0.0);
    }
}

TEST_CASE("TVR-3 step 2: log1p against its Taylor polynomial for |x| <= 1e-9 (4 u_dd)") {
    double worst = 0.0;
    for (const double m : log_spaced(-30.0, -9.0, 200)) {
        for (const double x : {m, -m}) {
            const DD X(x);
            const DD x2 = X * X;
            const DD taylor = X - ref::ldexp(x2, -1) + x2 * X / DD(3.0) - ref::ldexp(x2 * x2, -2);
            worst = std::max(worst, rel_udd(ref::log1p(X), taylor));
        }
    }
    INFO("worst " << worst << " u_dd");
    CHECK(worst <= 4.0);
}

TEST_CASE("TVR-3 step 2: identities (15.11)") {
    double expm1_log1p = 0.0;
    for (const double m : log_spaced(-25.0, 1.0, 300)) {
        for (const double x : {m, -m}) {
            if (!(x > -0.9)) continue;
            expm1_log1p = std::max(expm1_log1p, rel_udd(ref::expm1(ref::log1p(DD(x))), DD(x)));
        }
    }
    INFO("expm1(log1p x) = x: " << expm1_log1p << " u_dd");
    CHECK(expm1_log1p <= 16.0);

    double softplus = 0.0;
    for (const double x : inclusive(-600.0, 600.0, 1001)) {
        const DD d = ref::log1p_exp(DD(x)) - ref::log1p_exp(DD(-x));
        softplus = std::max(softplus, abs_udd(d, DD(x)) / std::max(1.0, std::abs(x)));
    }
    INFO("log1p_exp(x) - log1p_exp(-x) = x: " << softplus << " u_dd max(1, |x|)");
    CHECK(softplus <= 8.0);

    double tanh_double = 0.0;
    for (const double m : log_spaced(-20.0, std::log10(40.0), 300)) {
        for (const double x : {m, -m}) {
            const DD t = ref::tanh(DD(x));
            const DD rhs = ref::ldexp(t, 1) / (DD(1.0) + t * t);
            tanh_double = std::max(tanh_double, rel_udd(ref::tanh(ref::ldexp(DD(x), 1)), rhs));
        }
    }
    INFO("tanh(2x) = 2t / (1 + t^2): " << tanh_double << " u_dd");
    CHECK(tanh_double <= 16.0);

    double pythagoras = 0.0, sin_asin = 0.0, atan_pair = 0.0;
    for (const double y : inclusive(-1.0, 1.0, 401)) {
        const DD s = ref::sin(DD(y)), c = ref::cos(DD(y));
        pythagoras = std::max(pythagoras, abs_udd(s * s + c * c, DD(1.0)));
        if (y != 0.0) sin_asin = std::max(sin_asin, rel_udd(ref::sin(ref::asin(DD(y))), DD(y)));
    }
    for (const double x : log_spaced(-12.0, 12.0, 400)) {
        atan_pair = std::max(atan_pair, abs_udd(ref::atan(DD(x)) + ref::atan(DD(1.0) / DD(x)),
                                                ref::ldexp(ref::pi(), -1)));
    }
    INFO("sin^2 + cos^2 = 1: " << pythagoras << "; sin(asin x) = x: " << sin_asin
                               << "; atan x + atan 1/x = pi/2: " << atan_pair << " u_dd");
    CHECK(pythagoras <= 8.0);
    CHECK(sin_asin <= 16.0);
    CHECK(atan_pair <= 8.0);
}

// Added after the mutation run (DECISIONS 15.12; checks added, none loosened): atan2's x < 0 branch (mutant M16)
// and erfc below 3, where erfc = 1 - erf is absolutely, not relatively, accurate (mutant M20).
TEST_CASE("TVR-3 step 2: atan2 in the left half-plane (8 u_dd)") {
    const DD pi = ref::pi();
    const DD three_quarters = ref::ldexp(pi * DD(3.0), -2);
    const double a = rel_udd(ref::atan2(DD(1.0), DD(-1.0)), three_quarters);
    const double b = rel_udd(ref::atan2(DD(-1.0), DD(-1.0)), -three_quarters);
    const double c = rel_udd(ref::atan2(DD(0.5), -ref::sqrt(DD(3.0)) * DD(0.5)), pi * DD(5.0) / DD(6.0));
    INFO("atan2(1, -1): " << a << ", atan2(-1, -1): " << b << ", atan2(1/2, -sqrt3/2): " << c << " u_dd");
    CHECK(a <= 8.0);
    CHECK(b <= 8.0);
    CHECK(c <= 8.0);
}

TEST_CASE("TVR-3 step 2: erfc on [0, 3) against Gauss-Kronrod, absolutely (16 u_dd + U_ref)") {
    const DD two_over_sqrt_pi = DD(2.0) / ref::sqrt(ref::pi());
    const auto gauss = [&](const DD& t) { return two_over_sqrt_pi * ref::exp(-(t * t)); };
    double worst_excess = -1e300, worst_abs = 0.0;
    int not_ok = 0;
    for (const double x : inclusive(0.0, 2.99, 64)) {
        const auto r = ref::qag15(gauss, DD(x), DD(x + 15.0), 0.0, quad_epsrel);
        if (r.status != ref::QuadratureStatus::ok) ++not_ok;
        const double diff = std::abs((ref::erfc(DD(x)) - r.value).hi);
        worst_abs = std::max(worst_abs, diff / u_dd);
        worst_excess = std::max(worst_excess, (diff - 16.0 * u_dd - r.abs_error) / u_dd);
    }
    INFO("worst |erfc - I| = " << worst_abs << " u_dd absolute; excess over the bound " << worst_excess);
    CHECK(not_ok == 0);
    CHECK(worst_excess <= 0.0);
}
