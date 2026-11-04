// TC5 (DECISIONS Unit 55 section 11.3): L0 tests of Brent's procedure zero.
//   - the root of x^3 - 2x - 5 (Wallis's number; 64 digits, OEIS A007493, PA) and of cos x - x (the Dottie
//     number; 68 digits, OEIS A003957, PA), each to 4 u_dd relative, compared exactly with the printed digits;
//   - the root 0 of e^x - 1 to 4 u_dd absolute (a relative error has no meaning at 0; f'(0) = 1 sets the scale);
//     with t = u_dd, Brent's (2.18) guarantees 2t = 2 u_dd there;
//   - each count of function evaluations within (3.4), (k + 1)^2 - 2;
//   - an interval without a sign change is refused.
// cos is not in the reference library (TC1); a test-local Taylor series serves, |x| <= 1, its terms decreasing
// with no cancellation beyond the alternation (cos x >= 0.54 there), and is checked below by two identities,
// cos^2 x + sin^2 x = 1 (sin by its own series) and cos x = 2 cos^2(x/2) - 1, on a grid of [-1, 1].
#include "NiTCAD/verify/reference/brent.hpp"

#include "exact_integer.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using ref::u_dd;

namespace {

// cos and sin by their Taylor series, |x| <= 1.
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
DD sin_series(const DD& x) {
    const DD x2 = x * x;
    DD term = x, sum = x;
    for (int k = 1; k < 40; ++k) {
        term = -(term * x2) / DD((2.0 * k) * (2.0 * k + 1.0));
        sum += term;
        if (std::abs(term.hi) < 0x1p-112 * std::abs(sum.hi)) break;
    }
    return sum;
}

}  // namespace

TEST_CASE("TC5: the test-local cos is accurate to 2 u_dd on [-1, 1] (identities)") {
    double worst_pyth = 0.0, worst_double = 0.0;
    for (int i = -100; i <= 100; ++i) {
        const DD x = DD(i / 100.0);
        const DD c = cos_series(x), s = sin_series(x);
        worst_pyth = std::max(worst_pyth, std::abs((c * c + s * s - DD(1.0)).to_double()) / u_dd);
        const DD h = ldexp(x, -1);
        const DD ch = cos_series(h);
        worst_double = std::max(worst_double, std::abs((ldexp(ch * ch, 1) - DD(1.0) - c).to_double()) / u_dd);
    }
    INFO("cos^2 + sin^2 - 1: " << worst_pyth << " u_dd; 2 cos^2(x/2) - 1 - cos x: " << worst_double << " u_dd");
    CHECK(worst_pyth <= 2.0);
    CHECK(worst_double <= 2.0);
}

TEST_CASE("TC5: roots of x^3 - 2x - 5 and cos x - x to 4 u_dd, within Brent's evaluation bound") {
    const double t = u_dd;
    struct Case {
        const char* name;
        DD (*f)(const DD&);
        DD a, b;
        const char* digits;
    };
    const Case cases[] = {
        {"x^3 - 2x - 5 on [2, 3]", [](const DD& x) { return x * x * x - ldexp(x, 1) - DD(5.0); }, DD(2.0), DD(3.0),
         "2.094551481542326591482386540579302963857306105628239180304128529"},
        {"cos x - x on [0, 1]", [](const DD& x) { return cos_series(x) - x; }, DD(0.0), DD(1.0),
         "0.73908513321516064165531208767387340401341175890075746496568063577328"},
    };
    for (const auto& c : cases) {
        const auto res = ref::brent_zero(c.f, c.a, c.b, t);
        const long long bound = ref::brent_evaluation_bound(c.a, c.b, t);
        const auto cmp = exact::compare_with_literal(res.root, c.digits, 4);
        INFO(c.name << ": " << cmp.measured_udd << " u_dd, " << res.evaluations << " evaluations, bound " << bound);
        CHECK(res.status == ref::BrentStatus::ok);
        CHECK(cmp.within);
        CHECK(res.evaluations <= bound);
    }
}

TEST_CASE("TC5: the root 0 of e^x - 1 to 4 u_dd absolute, within Brent's evaluation bound") {
    const double t = u_dd;
    const DD a(-1.0), b(2.0);
    const auto res = ref::brent_zero([](const DD& x) { return ref::expm1(x); }, a, b, t);
    const long long bound = ref::brent_evaluation_bound(a, b, t);
    const double err = std::abs(res.root.to_double()) / u_dd;
    INFO("|root| " << err << " u_dd, " << res.evaluations << " evaluations, bound " << bound);
    CHECK(res.status == ref::BrentStatus::ok);
    CHECK(err <= 4.0);
    CHECK(res.evaluations <= bound);
}

TEST_CASE("TC5: superlinear on the smooth roots: at most a quarter of bisection's k + 1 evaluations") {
    // Section 3: bisection to the same tolerance takes exactly k + 1 evaluations.
    const double t = u_dd;
    auto k_plus_1 = [t](const DD& a, const DD& b) {
        const double dm = 2.0 * u_dd * std::min(std::abs(a.hi), std::abs(b.hi)) + t;
        return std::ceil(std::log2((b - a).hi / dm)) + 1.0;
    };
    const auto w = ref::brent_zero([](const DD& x) { return x * x * x - ldexp(x, 1) - DD(5.0); }, DD(2.0), DD(3.0), t);
    const auto d = ref::brent_zero([](const DD& x) { return cos_series(x) - x; }, DD(0.0), DD(1.0), t);
    INFO("Wallis " << w.evaluations << " of k + 1 = " << k_plus_1(DD(2.0), DD(3.0)) << "; Dottie " << d.evaluations
                   << " of " << k_plus_1(DD(0.0), DD(1.0)));
    CHECK(static_cast<double>(w.evaluations) <= k_plus_1(DD(2.0), DD(3.0)) / 4.0);
    CHECK(static_cast<double>(d.evaluations) <= k_plus_1(DD(0.0), DD(1.0)) / 4.0);
}

TEST_CASE("TC5: Brent's example (2.1), on which Dekker's algorithm takes ~(b - a)/delta steps, within (3.4)") {
    // f = 2^(x/delta) on [a + delta, b]; f(a) = -((b - a - delta)/delta) 2^(b/delta); "arbitrary" on (a, a + delta),
    // here linear between those two values. a = 0, b = 1, delta = 2^-8: Dekker's steps of delta would need ~256
    // evaluations; (3.4) allows (8 + 1)^2 - 2 = 79 with t = delta (2 u_dd |b| adds less than one part in 2^95).
    const double delta = 0x1p-8;
    const DD a(0.0), b(1.0);
    const DD fa = -DD((1.0 - delta) / delta) * ldexp(DD(1.0), 256);
    const DD f_left = DD(2.0);  // 2^((a + delta)/delta)
    auto f = [&](const DD& x) {
        if (!(x < DD(delta))) return ref::exp(x / DD(delta) * ref::ln2());
        return fa + (f_left - fa) * (x / DD(delta));
    };
    const auto res = ref::brent_zero(f, a, b, delta);
    const long long bound = ref::brent_evaluation_bound(a, b, delta);
    INFO(res.evaluations << " evaluations, bound " << bound << ", root " << res.root.hi);
    CHECK(res.status == ref::BrentStatus::ok);
    CHECK(res.evaluations <= bound);
    // The root of the linear piece, x = -fa delta / (f_left - fa), within Brent's (2.18) 6 eps |x| + 2t.
    const DD root = -fa * DD(delta) / (f_left - fa);
    CHECK(std::abs((res.root - root).to_double()) <= 6.0 * u_dd * std::abs(root.hi) + 2.0 * delta);
}

TEST_CASE("TC5: Brent's Table 4.1 (x^9, a multiple root): evaluation counts within 5% of the published ones") {
    // Table 4.1 (p. 55; procedure zero2, macheps = 16^-13 on an IBM 360): x^9 on [-1, 1.1], t = 1e-9: 81
    // evaluations; on [-1, 4], t = 1e-20: 189. Here macheps = u_dd, so the counts differ by rounding only; the 5%
    // band was set after measuring 78 and 191 (DECISIONS 11.14) and separates the four interpolation and
    // bookkeeping mutations the frozen criteria do not. The x^19 row is not used: x^19 underflows in double
    // below |x| ~ 1e-16, so f vanishes on an interval there (Brent's own remark on Dekker's procedure), and Brent
    // ran it with zero2's extended exponent.
    auto x9 = [](const DD& x) {
        DD p = x;
        for (int i = 1; i < 9; ++i) p = p * x;
        return p;
    };
    const auto r1 = ref::brent_zero(x9, DD(-1.0), DD(1.1), 1e-9);
    const auto r2 = ref::brent_zero(x9, DD(-1.0), DD(4.0), 1e-20);
    INFO(r1.evaluations << " (published 81), " << r2.evaluations << " (published 189); roots " << r1.root.hi << ", "
                        << r2.root.hi);
    CHECK(r1.status == ref::BrentStatus::ok);
    CHECK(r2.status == ref::BrentStatus::ok);
    CHECK(std::abs(static_cast<double>(r1.evaluations) - 81.0) <= 0.05 * 81.0);
    CHECK(std::abs(static_cast<double>(r2.evaluations) - 189.0) <= 0.05 * 189.0);
    // (2.18): |root| <= 6 eps |root| + 2t, the true root being 0.
    CHECK(std::abs(r1.root.hi) <= 2e-9 * (1.0 + 1e-20));
    CHECK(std::abs(r2.root.hi) <= 2e-20 * (1.0 + 1e-20));
}

TEST_CASE("TC5: an interval without a sign change is refused; the bound's formula") {
    const auto res = ref::brent_zero([](const DD& x) { return x * x + DD(1.0); }, DD(-1.0), DD(1.0), u_dd);
    CHECK(res.status == ref::BrentStatus::not_bracketed);
    CHECK(res.evaluations == 2);
    const auto bad = ref::brent_zero([](const DD& x) { return x; }, DD(-1.0), DD(1.0), 0.0);
    CHECK(bad.status == ref::BrentStatus::invalid_input);
    // (3.2)-(3.4) by hand: [2, 3], t = 2^-10, macheps = 2^-20: delta_m = 2^-19 * 2 + 2^-10 ~ 2^-10, k = 10, 119.
    CHECK(ref::brent_evaluation_bound(DD(2.0), DD(3.0), 0x1p-10, 0x1p-20) == 119);
    CHECK(ref::brent_evaluation_bound(DD(-1.0), DD(1.0), 4.0) == 2);  // k <= 0
}
