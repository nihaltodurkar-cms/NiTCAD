// TC2 (DECISIONS Unit 55 section 11.3; TVR1-Q1): the new double-double library and the existing test
// helper tests/physics/reference_arithmetic.hpp agree to 2 u_dd relative on the TC1 grids (log absolutely
// and pow scaled by its condition number, owner ruling 2026-11-04). The helper is
// included here only (never by the library, TVD-14) and is not changed. The two share the published
// add/multiply/divide algorithms but not sqrt (Newton on x^-2 here, one step on x there), exp (k = 256 here,
// 1024 there), log (an atanh branch near 1 here) or the constants' formulas (Gauss and acoth 7, 17 here;
// Machin and atanh(1/3) there), so agreement of the functions is evidence about both.
#include "NiTCAD/verify/reference/double_double.hpp"

#include "../../physics/reference_arithmetic.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace ref = NiTCAD::verify::reference;

namespace {

constexpr int n_points = 10000;
constexpr double tolerance_udd = 2.0;

reference::DD to_helper(const ref::DD& x) { return {x.hi, x.lo}; }

// |a - b| / |b| in u_dd, computed in the new library (b nonzero).
double diff_udd(const ref::DD& a, const reference::DD& b) {
    const ref::DD d = a - ref::DD(b.hi, b.lo);
    return std::abs(d.to_double()) / std::abs(b.hi) / ref::u_dd;
}

ref::DD make_dd(double x, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> t(-0.5, 0.5);
    return ref::quick_two_sum(x, x * 0x1p-53 * t(rng));
}

std::vector<double> log_spaced(double lo, double hi, int n) {
    std::vector<double> v;
    for (int i = 0; i < n; ++i) v.push_back(std::pow(10.0, lo + (hi - lo) * i / (n - 1)));
    return v;
}

struct Worst {
    double value = 0.0;
    double at = 0.0;
    void update(double e, double x) {
        if (e > value) {
            value = e;
            at = x;
        }
    }
};

}  // namespace

TEST_CASE("TC2: +, -, *, / agree with the test helper to 2 u_dd") {
    std::mt19937_64 rng(11);
    const auto xs = log_spaced(-140.0, 140.0, n_points);
    auto ys = xs;
    std::shuffle(ys.begin(), ys.end(), rng);
    Worst add, sub, mul, div;
    for (int i = 0; i < n_points; ++i) {
        const ref::DD x = make_dd(xs[i], rng);
        const ref::DD y = make_dd(i % 2 == 0 ? ys[i] : -ys[i], rng);
        const reference::DD hx = to_helper(x), hy = to_helper(y);
        add.update(diff_udd(x + y, hx + hy), xs[i]);
        sub.update(diff_udd(x - y, hx - hy), xs[i]);
        mul.update(diff_udd(x * y, hx * hy), xs[i]);
        div.update(diff_udd(x / y, hx / hy), xs[i]);
    }
    INFO("worst +: " << add.value << ", -: " << sub.value << ", *: " << mul.value << ", /: " << div.value
                     << " u_dd");
    CHECK(add.value <= tolerance_udd);
    CHECK(sub.value <= tolerance_udd);
    CHECK(mul.value <= tolerance_udd);
    CHECK(div.value <= tolerance_udd);
}

TEST_CASE("TC2: sqrt agrees with the test helper to 2 u_dd") {
    std::mt19937_64 rng(12);
    Worst w;
    for (const double v : log_spaced(-280.0, 280.0, n_points)) {
        const ref::DD x = make_dd(v, rng);
        w.update(diff_udd(ref::sqrt(x), reference::sqrt(to_helper(x))), v);
    }
    INFO("worst " << w.value << " u_dd at x = " << w.at);
    CHECK(w.value <= tolerance_udd);
}

TEST_CASE("TC2: exp and expm1 agree with the test helper to 2 u_dd") {
    std::mt19937_64 rng(13);
    std::bernoulli_distribution sign(0.5);
    Worst we, wm;
    for (const double v : log_spaced(-10.0, std::log10(670.0), n_points)) {
        const ref::DD a = make_dd(sign(rng) ? v : -v, rng);
        we.update(diff_udd(ref::exp(a), reference::exp(to_helper(a))), a.hi);
    }
    for (const double v : log_spaced(-30.0, std::log10(2.0), n_points)) {
        const ref::DD a = make_dd(sign(rng) ? v : -v, rng);
        wm.update(diff_udd(ref::expm1(a), reference::expm1(to_helper(a))), a.hi);
    }
    INFO("worst exp " << we.value << " u_dd at a = " << we.at << "; expm1 " << wm.value << " u_dd at a = " << wm.at);
    CHECK(we.value <= tolerance_udd);
    CHECK(wm.value <= tolerance_udd);
}

// Owner ruling 2026-11-04 (DECISIONS 11.13): log is compared absolutely, |difference| <= 2 u_dd max(1, |ln x|),
// which is the helper's guarantee (its Newton log is absolute-accurate and loses relative precision near 1,
// TVR-0 finding I-11); pow relative to 2 u_dd max(1, |b|, |b ln a|): the helper's absolute log error times |b|
// (second ruling, after a measured 8.3 miss at a = 1.024 under max(1, |b ln a|), also I-11).
TEST_CASE("TC2: log (absolute) and pow (condition-scaled) agree with the test helper to 2 u_dd") {
    std::mt19937_64 rng(14);
    Worst wl, wp;
    for (const double v : log_spaced(-280.0, 280.0, n_points)) {
        const ref::DD x = make_dd(v, rng);
        const ref::DD mine = ref::log(x);
        const reference::DD theirs = reference::log(to_helper(x));
        const ref::DD d = mine - ref::DD(theirs.hi, theirs.lo);
        wl.update(std::abs(d.to_double()) / std::max(1.0, std::abs(mine.hi)) / ref::u_dd, v);
    }
    const auto bases = log_spaced(-3.0, 3.0, n_points);
    std::uniform_real_distribution<double> expo(-40.0, 40.0);
    for (const double v : bases) {
        const ref::DD x = make_dd(v, rng);
        const ref::DD b = make_dd(expo(rng), rng);
        const double scale = std::max({1.0, std::abs(b.hi), std::abs(b.hi * std::log(v))});
        wp.update(diff_udd(ref::pow(x, b), reference::pow(to_helper(x), to_helper(b))) / scale, v);
    }
    INFO("worst log " << wl.value << " u_dd at x = " << wl.at << "; pow " << wp.value << " u_dd at x = " << wp.at);
    CHECK(wl.value <= tolerance_udd);
    CHECK(wp.value <= tolerance_udd);
}

TEST_CASE("TC2: pi and ln 2 agree with the test helper to 2 u_dd") {
    const double dpi = diff_udd(ref::pi(), reference::pi());
    const double dln2 = diff_udd(ref::ln2(), reference::ln2());
    INFO("pi " << dpi << " u_dd, ln 2 " << dln2 << " u_dd");
    CHECK(dpi <= tolerance_udd);
    CHECK(dln2 <= tolerance_udd);
}
