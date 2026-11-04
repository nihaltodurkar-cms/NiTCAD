// TC1 (DECISIONS Unit 55 section 11.3, as amended 2026-11-04): L0 tests of the double-double library.
//   - error-free transformations exact, checked by integer arithmetic on random integer-scaled inputs;
//   - pi, e, ln 2, sqrt 2 within 1 u_dd relative of OEIS's 40 digits (A000796, A001113, A002162,
//     A002193; read 2026-11-04; A002162 re-read against its b-file after a first transcription
//     dropped a digit), compared exactly with big integers (no decimal parsing in floating point);
//     the 1 u_dd tolerance is the owner's ruling of 2026-11-04 on "to 32 digits";
//   - identities on 10^4 log-spaced double-double points: (sqrt x)^2 = x and (x / y) y = x to 8 u_dd;
//     exp(log x) = x and exp(a + b) = exp a exp b to 8 u_dd times the identity's condition number
//     (max(1, |ln x|), max(1, |a| + |b|); owner ruling 2026-11-04, since no double-double log can carry
//     more than about 2^-106 |ln x| absolute precision).
// Additional identities (not in the frozen list; nothing is loosened by adding them): expm1(2x) =
// expm1(x) (expm1(x) + 2) and expm1(log(1 + d)) = d, both conditioned by 1, to 8 u_dd relative; and
// double-double + (under cancellation) and * against the exact integer result, to 8 u_dd.
// Ranges are where double-double precision exists (lo normal): |x| in [1e-280, 1e280], quotients too.
#include "NiTCAD/verify/reference/double_double.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using ref::u_dd;

namespace {

// ---- exact integers --------------------------------------------------------------------------

// Non-negative big integer in base 2^32 (little-endian limbs).
struct Big {
    std::vector<std::uint32_t> w;

    static Big from_u64(std::uint64_t v) {
        Big b;
        while (v != 0) {
            b.w.push_back(static_cast<std::uint32_t>(v));
            v >>= 32;
        }
        return b;
    }
    void trim() {
        while (!w.empty() && w.back() == 0) w.pop_back();
    }
    void mul_small(std::uint32_t m) {
        std::uint64_t carry = 0;
        for (auto& limb : w) {
            const std::uint64_t t = static_cast<std::uint64_t>(limb) * m + carry;
            limb = static_cast<std::uint32_t>(t);
            carry = t >> 32;
        }
        if (carry != 0) w.push_back(static_cast<std::uint32_t>(carry));
    }
    void add_small(std::uint32_t a) {
        std::uint64_t carry = a;
        for (std::size_t i = 0; carry != 0; ++i) {
            if (i == w.size()) w.push_back(0);
            const std::uint64_t t = static_cast<std::uint64_t>(w[i]) + carry;
            w[i] = static_cast<std::uint32_t>(t);
            carry = t >> 32;
        }
    }
    void shl(int bits) {
        if (w.empty()) return;
        const int limbs = bits / 32;
        const int rest = bits % 32;
        w.insert(w.begin(), static_cast<std::size_t>(limbs), 0u);
        if (rest != 0) {
            std::uint32_t carry = 0;
            for (std::size_t i = static_cast<std::size_t>(limbs); i < w.size(); ++i) {
                const std::uint32_t next = w[i] >> (32 - rest);
                w[i] = (w[i] << rest) | carry;
                carry = next;
            }
            if (carry != 0) w.push_back(carry);
        }
    }
    [[nodiscard]] double to_double() const {
        double v = 0.0;
        for (std::size_t i = w.size(); i-- > 0;) v = v * 4294967296.0 + w[i];
        return v;
    }
};

int cmp(const Big& a, const Big& b) {
    if (a.w.size() != b.w.size()) return a.w.size() < b.w.size() ? -1 : 1;
    for (std::size_t i = a.w.size(); i-- > 0;)
        if (a.w[i] != b.w[i]) return a.w[i] < b.w[i] ? -1 : 1;
    return 0;
}
Big add(const Big& a, const Big& b) {
    Big r;
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < std::max(a.w.size(), b.w.size()) || carry != 0; ++i) {
        const std::uint64_t t = carry + (i < a.w.size() ? a.w[i] : 0u) + (i < b.w.size() ? b.w[i] : 0u);
        r.w.push_back(static_cast<std::uint32_t>(t));
        carry = t >> 32;
    }
    r.trim();
    return r;
}
// a - b for a >= b.
Big sub(const Big& a, const Big& b) {
    Big r = a;
    std::int64_t borrow = 0;
    for (std::size_t i = 0; i < r.w.size(); ++i) {
        std::int64_t t = static_cast<std::int64_t>(r.w[i]) - borrow - (i < b.w.size() ? b.w[i] : 0u);
        borrow = t < 0 ? 1 : 0;
        if (t < 0) t += 4294967296LL;
        r.w[i] = static_cast<std::uint32_t>(t);
    }
    r.trim();
    return r;
}

// Signed big integer.
struct SBig {
    Big mag;
    bool neg = false;
};
SBig sadd(const SBig& a, const SBig& b) {
    if (a.neg == b.neg) return {add(a.mag, b.mag), a.neg};
    const int c = cmp(a.mag, b.mag);
    if (c == 0) return {};
    return c > 0 ? SBig{sub(a.mag, b.mag), a.neg} : SBig{sub(b.mag, a.mag), b.neg};
}
bool operator==(const SBig& a, const SBig& b) {
    if (a.mag.w.empty() && b.mag.w.empty()) return true;
    return a.neg == b.neg && cmp(a.mag, b.mag) == 0;
}

// d = M 2^k with M an integer of at most 53 bits (M = 0 for d = 0).
struct Scaled {
    std::int64_t m = 0;
    int k = 0;
};
Scaled decompose(double d) {
    if (d == 0.0) return {};
    int k = 0;
    const double f = std::frexp(d, &k);
    Scaled s{static_cast<std::int64_t>(std::ldexp(f, 53)), k - 53};
    while (s.m % 2 == 0) {  // odd M: the largest k with d on the grid 2^k
        s.m /= 2;
        ++s.k;
    }
    return s;
}
// M 2^(k - base) as a signed big integer, k >= base.
SBig to_big(std::int64_t m, int k, int base) {
    SBig r{Big::from_u64(static_cast<std::uint64_t>(m < 0 ? -m : m)), m < 0};
    r.mag.shl(k - base);
    return r;
}
// |m1 m2| exactly.
Big product(std::int64_t m1, std::int64_t m2) {
    Big r = Big::from_u64(static_cast<std::uint64_t>(m1 < 0 ? -m1 : m1));
    const auto b = static_cast<std::uint64_t>(m2 < 0 ? -m2 : m2);
    // r * b by two 32-bit halves.
    Big lo = r;
    lo.mul_small(static_cast<std::uint32_t>(b));
    Big hi = r;
    hi.mul_small(static_cast<std::uint32_t>(b >> 32));
    hi.shl(32);
    return add(lo, hi);
}

// ---- helpers ---------------------------------------------------------------------------------

// A double-double with a nonzero lo part near x: hi = x, lo = x 2^-53 t, t in (-1/2, 1/2).
DD make_dd(double x, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> t(-0.5, 0.5);
    return ref::quick_two_sum(x, x * 0x1p-53 * t(rng));
}

// n log-spaced doubles 10^lo .. 10^hi.
std::vector<double> log_spaced(double lo, double hi, int n) {
    std::vector<double> v;
    for (int i = 0; i < n; ++i) v.push_back(std::pow(10.0, lo + (hi - lo) * i / (n - 1)));
    return v;
}

// |r - x| / |x| in units of u_dd.
double rel_err_udd(const DD& r, const DD& x) {
    const DD d = r - x;
    return std::abs(d.to_double()) / std::abs(x.hi) / u_dd;
}

constexpr int n_points = 10000;

}  // namespace

TEST_CASE("TC1: two_sum is exact (integer check)") {
    std::mt19937_64 rng(20261104);
    std::uniform_int_distribution<std::int64_t> mant(-(std::int64_t{1} << 53) + 1, (std::int64_t{1} << 53) - 1);
    std::uniform_int_distribution<int> exp_a(-60, 60);
    std::uniform_int_distribution<int> offset(-40, 40);
    int failures = 0;
    for (int i = 0; i < n_points; ++i) {
        const std::int64_t ma = mant(rng);
        const std::int64_t mb = mant(rng);
        const int ea = exp_a(rng);
        const int eb = ea + offset(rng);
        const double a = std::ldexp(static_cast<double>(ma), ea);
        const double b = std::ldexp(static_cast<double>(mb), eb);
        const DD s = ref::two_sum(a, b);
        const Scaled sh = decompose(s.hi);
        const Scaled sl = decompose(s.lo);
        const int base = std::min(ea, eb);
        // s.hi and s.lo lie on the grid 2^base (DECISIONS 11.3 TC1 note: exact sum on that grid).
        const bool on_grid = (s.hi == 0.0 || sh.k >= base) && (s.lo == 0.0 || sl.k >= base);
        if (!on_grid) {
            ++failures;
            continue;
        }
        const SBig exact = sadd(to_big(ma, ea, base), to_big(mb, eb, base));
        const SBig got = sadd(to_big(sh.m, sh.k, base), to_big(sl.m, sl.k, base));
        if (!(exact == got)) ++failures;
    }
    CHECK(failures == 0);
}

TEST_CASE("TC1: quick_two_sum is exact for |a| >= |b| (integer check)") {
    std::mt19937_64 rng(20261105);
    std::uniform_int_distribution<std::int64_t> mant(-(std::int64_t{1} << 53) + 1, (std::int64_t{1} << 53) - 1);
    std::uniform_int_distribution<int> exp_a(-60, 60);
    std::uniform_int_distribution<int> offset(0, 40);
    int failures = 0;
    for (int i = 0; i < n_points; ++i) {
        std::int64_t ma = mant(rng);
        std::int64_t mb = mant(rng);
        const int ea = exp_a(rng);
        const int eb = ea - offset(rng);
        if (eb == ea && std::abs(mb) > std::abs(ma)) std::swap(ma, mb);
        const double a = std::ldexp(static_cast<double>(ma), ea);
        const double b = std::ldexp(static_cast<double>(mb), eb);
        if (std::abs(a) < std::abs(b)) continue;  // only a few draws, mantissa of a short
        const DD s = ref::quick_two_sum(a, b);
        const Scaled sh = decompose(s.hi);
        const Scaled sl = decompose(s.lo);
        const int base = eb;
        if ((s.hi != 0.0 && sh.k < base) || (s.lo != 0.0 && sl.k < base)) {
            ++failures;
            continue;
        }
        const SBig exact = sadd(to_big(ma, ea, base), to_big(mb, eb, base));
        const SBig got = sadd(to_big(sh.m, sh.k, base), to_big(sl.m, sl.k, base));
        if (!(exact == got)) ++failures;
    }
    CHECK(failures == 0);
}

TEST_CASE("TC1: two_prod is exact (integer check)") {
    std::mt19937_64 rng(20261106);
    std::uniform_int_distribution<std::int64_t> mant(-(std::int64_t{1} << 53) + 1, (std::int64_t{1} << 53) - 1);
    std::uniform_int_distribution<int> expo(-200, 200);
    int failures = 0;
    for (int i = 0; i < n_points; ++i) {
        const std::int64_t ma = mant(rng);
        const std::int64_t mb = mant(rng);
        const int ea = expo(rng);
        const int eb = expo(rng);
        const double a = std::ldexp(static_cast<double>(ma), ea);
        const double b = std::ldexp(static_cast<double>(mb), eb);
        const DD p = ref::two_prod(a, b);
        const Scaled ph = decompose(p.hi);
        const Scaled pl = decompose(p.lo);
        const int base = ea + eb;
        if ((p.hi != 0.0 && ph.k < base) || (p.lo != 0.0 && pl.k < base)) {
            ++failures;
            continue;
        }
        const SBig exact{product(ma, mb), (ma < 0) != (mb < 0)};
        const SBig got = sadd(to_big(ph.m, ph.k, base), to_big(pl.m, pl.k, base));
        if (!(exact == got)) ++failures;
    }
    CHECK(failures == 0);
}

namespace {

// A random double-double on an integer grid: hi = A 2^e with |A| in [2^52, 2^53) (so ulp(hi) = 2^e),
// lo = B 2^(e - 53 - k), k in [1, 10], |B| < 2^53 (so |lo| < ulp(hi) / 2).
struct ExactDD {
    std::int64_t a = 0;
    int ea = 0;
    std::int64_t b = 0;
    int eb = 0;
    [[nodiscard]] DD value() const {
        return {std::ldexp(static_cast<double>(a), ea), std::ldexp(static_cast<double>(b), eb)};
    }
};
ExactDD random_exact(std::mt19937_64& rng, int e) {
    std::uniform_int_distribution<std::int64_t> top(std::int64_t{1} << 52, (std::int64_t{1} << 53) - 1);
    std::uniform_int_distribution<std::int64_t> any(-(std::int64_t{1} << 53) + 1, (std::int64_t{1} << 53) - 1);
    std::uniform_int_distribution<int> k(1, 10);
    std::bernoulli_distribution sign(0.5);
    ExactDD x;
    x.a = sign(rng) ? top(rng) : -top(rng);
    x.ea = e;
    x.b = any(rng);
    x.eb = e - 53 - k(rng);
    return x;
}
// |r - exact| <= 8 u_dd |exact|, i.e. 2^104 |r - exact| <= 8 |exact|, everything on the grid 2^base.
bool within_8_udd(const DD& r, const SBig& exact, int base) {
    const Scaled h = decompose(r.hi);
    const Scaled l = decompose(r.lo);
    if ((r.hi != 0.0 && h.k < base) || (r.lo != 0.0 && l.k < base)) {
        // r carries bits below the exact result's grid: compare on the finer grid instead.
        const int fine = std::min(r.hi != 0.0 ? h.k : base, r.lo != 0.0 ? l.k : base);
        SBig ex = exact;
        ex.mag.shl(base - fine);
        return within_8_udd(r, ex, fine);
    }
    const SBig got = sadd(to_big(h.m, h.k, base), to_big(l.m, l.k, base));
    SBig neg_exact = exact;
    neg_exact.neg = !neg_exact.neg;
    Big diff = sadd(got, neg_exact).mag;
    if (exact.mag.w.empty()) return diff.w.empty();
    diff.shl(104);
    Big bound = exact.mag;
    bound.shl(3);
    return cmp(diff, bound) <= 0;
}

}  // namespace

// The IEEE-style addition must keep the low-order error terms under cancellation; the identity tests
// cannot see them (they are below 1 u_dd of operands that do not cancel). The exact sum is integer.
TEST_CASE("TC1 (added): double-double + under cancellation within 8 u_dd of the exact sum (integer check)") {
    std::mt19937_64 rng(7);
    std::uniform_int_distribution<int> expo(-100, 100);
    std::uniform_int_distribution<int> mode(0, 2);
    std::uniform_int_distribution<std::int64_t> nudge(-1000, 1000);
    int failures = 0;
    for (int i = 0; i < n_points; ++i) {
        const ExactDD x = random_exact(rng, expo(rng));
        ExactDD y = random_exact(rng, x.ea);
        switch (mode(rng)) {
            case 0: y.a = -x.a; break;  // hi parts cancel exactly
            case 1:  // hi parts nearly cancel
                y.a = -x.a + nudge(rng);
                if (std::abs(y.a) < (std::int64_t{1} << 52)) y.a = -x.a;
                break;
            default: y = random_exact(rng, x.ea - 20 + expo(rng) % 20); break;
        }
        const int base = std::min({x.ea, x.eb, y.ea, y.eb});
        SBig exact = sadd(sadd(to_big(x.a, x.ea, base), to_big(x.b, x.eb, base)),
                          sadd(to_big(y.a, y.ea, base), to_big(y.b, y.eb, base)));
        if (!within_8_udd(x.value() + y.value(), exact, base)) ++failures;
    }
    CHECK(failures == 0);
}

TEST_CASE("TC1 (added): double-double * within 8 u_dd of the exact product (integer check)") {
    std::mt19937_64 rng(8);
    std::uniform_int_distribution<int> expo(-100, 100);
    int failures = 0;
    for (int i = 0; i < n_points; ++i) {
        const ExactDD x = random_exact(rng, expo(rng));
        const ExactDD y = random_exact(rng, expo(rng));
        const int base = x.eb + y.eb;
        auto term = [&](std::int64_t m1, int e1, std::int64_t m2, int e2) {
            SBig t{product(m1, m2), (m1 < 0) != (m2 < 0)};
            t.mag.shl(e1 + e2 - base);
            if (t.mag.w.empty()) t.neg = false;
            return t;
        };
        const SBig exact = sadd(sadd(term(x.a, x.ea, y.a, y.ea), term(x.a, x.ea, y.b, y.eb)),
                                sadd(term(x.b, x.eb, y.a, y.ea), term(x.b, x.eb, y.b, y.eb)));
        if (!within_8_udd(x.value() * y.value(), exact, base)) ++failures;
    }
    CHECK(failures == 0);
}

namespace {

// |c - L| <= 2^-104 |L| exactly, where L is the decimal literal (digits with one '.').
// With c = C 2^E (E < 0) and L = N / 10^d: 2^104 |C 10^d - N 2^-E| <= N 2^-E. Returns the measured
// |c - L| / |L| in u_dd (approximate, for the record) through `measured`.
bool within_one_udd(const DD& c, std::string_view literal, double& measured) {
    Big n;
    int d = 0;
    bool after_point = false;
    for (const char ch : literal) {
        if (ch == '.') {
            after_point = true;
            continue;
        }
        n.mul_small(10);
        n.add_small(static_cast<std::uint32_t>(ch - '0'));
        if (after_point) ++d;
    }
    const Scaled h = decompose(c.hi);
    const Scaled l = decompose(c.lo);
    const int base = std::min(h.k, c.lo == 0.0 ? h.k : l.k);
    REQUIRE(base < 0);
    SBig cbig = sadd(to_big(h.m, h.k, base), to_big(l.m, l.k, base));
    REQUIRE_FALSE(cbig.neg);
    for (int i = 0; i < d; ++i) cbig.mag.mul_small(10);  // C 10^d
    Big n_scaled = n;
    n_scaled.shl(-base);  // N 2^-E
    const int c_vs_n = cmp(cbig.mag, n_scaled);
    Big diff = c_vs_n >= 0 ? sub(cbig.mag, n_scaled) : sub(n_scaled, cbig.mag);
    measured = diff.to_double() / n_scaled.to_double() / u_dd;
    diff.shl(104);
    return cmp(diff, n_scaled) <= 0;
}

}  // namespace

TEST_CASE("TC1: constants within 1 u_dd of OEIS's digits (exact comparison)") {
    struct Case {
        const char* name;
        const DD& value;
        const char* oeis;
    };
    const Case cases[] = {
        {"pi (A000796)", ref::pi(), "3.1415926535897932384626433832795028841971"},
        {"e (A001113)", ref::e(), "2.7182818284590452353602874713526624977572"},
        {"ln 2 (A002162)", ref::ln2(), "0.6931471805599453094172321214581765680755"},
        {"sqrt 2 (A002193)", ref::sqrt2(), "1.4142135623730950488016887242096980785696"},
    };
    for (const auto& c : cases) {
        double measured = 0.0;
        const bool ok = within_one_udd(c.value, c.oeis, measured);
        INFO(c.name << ": |c - L| / L = " << measured << " u_dd");
        CHECK(ok);
    }
}

TEST_CASE("TC1: (sqrt x)^2 = x to 8 u_dd") {
    std::mt19937_64 rng(1);
    double worst = 0.0;
    for (const double v : log_spaced(-280.0, 280.0, n_points)) {
        const DD x = make_dd(v, rng);
        const DD s = ref::sqrt(x);
        worst = std::max(worst, rel_err_udd(s * s, x));
    }
    INFO("worst " << worst << " u_dd");
    CHECK(worst <= 8.0);
}

TEST_CASE("TC1: (x / y) y = x to 8 u_dd") {
    std::mt19937_64 rng(2);
    const auto xs = log_spaced(-140.0, 140.0, n_points);
    auto ys = xs;
    std::shuffle(ys.begin(), ys.end(), rng);
    double worst = 0.0;
    for (int i = 0; i < n_points; ++i) {
        const DD x = make_dd(xs[i], rng);
        const DD y = make_dd(i % 2 == 0 ? ys[i] : -ys[i], rng);
        worst = std::max(worst, rel_err_udd((x / y) * y, x));
    }
    INFO("worst " << worst << " u_dd");
    CHECK(worst <= 8.0);
}

TEST_CASE("TC1: exp(log x) = x to 8 max(1, |ln x|) u_dd") {
    std::mt19937_64 rng(3);
    double worst = 0.0;  // in units of the scaled tolerance
    for (const double v : log_spaced(-280.0, 280.0, n_points)) {
        const DD x = make_dd(v, rng);
        const double scale = std::max(1.0, std::abs(std::log(v)));
        worst = std::max(worst, rel_err_udd(ref::exp(ref::log(x)), x) / scale);
    }
    INFO("worst " << worst << " u_dd per unit condition number");
    CHECK(worst <= 8.0);
}

TEST_CASE("TC1: exp(a + b) = exp a exp b to 8 max(1, |a| + |b|) u_dd") {
    std::mt19937_64 rng(4);
    const auto mags = log_spaced(-10.0, std::log10(335.0), n_points);
    auto mags_b = mags;
    std::shuffle(mags_b.begin(), mags_b.end(), rng);
    std::bernoulli_distribution sign(0.5);
    double worst = 0.0;
    for (int i = 0; i < n_points; ++i) {
        const DD a = make_dd(sign(rng) ? mags[i] : -mags[i], rng);
        const DD b = make_dd(sign(rng) ? mags_b[i] : -mags_b[i], rng);
        const double scale = std::max(1.0, std::abs(a.hi) + std::abs(b.hi));
        const DD lhs = ref::exp(a + b);
        worst = std::max(worst, rel_err_udd(ref::exp(a) * ref::exp(b), lhs) / scale);
    }
    INFO("worst " << worst << " u_dd per unit condition number");
    CHECK(worst <= 8.0);
}

TEST_CASE("TC1 (added): expm1(2x) = expm1(x) (expm1(x) + 2) to 8 u_dd") {
    std::mt19937_64 rng(5);
    std::bernoulli_distribution sign(0.5);
    double worst = 0.0;
    for (const double v : log_spaced(-30.0, std::log10(2.0), n_points)) {
        const DD x = make_dd(sign(rng) ? v : -v, rng);
        const DD m = ref::expm1(x);
        worst = std::max(worst, rel_err_udd(m * (m + DD(2.0)), ref::expm1(ldexp(x, 1))));
    }
    INFO("worst " << worst << " u_dd");
    CHECK(worst <= 8.0);
}

// x = 1 + d is exact (two_sum of two doubles), so x - 1 = d exactly and the identity is conditioned
// by 1: it measures log's accuracy relative to its own value near its zero (the atanh branch).
TEST_CASE("TC1 (added): expm1(log(1 + d)) = d for |d| <= 0.1 to 8 u_dd relative") {
    std::mt19937_64 rng(6);
    std::bernoulli_distribution sign(0.5);
    double worst = 0.0;
    for (const double v : log_spaced(-30.0, -1.0, n_points)) {
        const double d = sign(rng) ? v : -v;
        const DD x = ref::two_sum(1.0, d);
        worst = std::max(worst, rel_err_udd(ref::expm1(ref::log(x)), DD(d)));
    }
    INFO("worst " << worst << " u_dd");
    CHECK(worst <= 8.0);
}

TEST_CASE("TC1: special values") {
    CHECK(ref::exp(DD(0.0)) == DD(1.0));
    CHECK(ref::exp(DD(800.0)).hi == std::numeric_limits<double>::infinity());
    CHECK(ref::exp(DD(-800.0)).hi == 0.0);
    CHECK(std::isnan(ref::log(DD(0.0)).hi));
    CHECK(std::isnan(ref::log(DD(-1.0)).hi));
    CHECK(ref::log(DD(1.0)).hi == 0.0);
    CHECK(ref::sqrt(DD(0.0)).hi == 0.0);
    CHECK(std::isnan(ref::sqrt(DD(-1.0)).hi));
    CHECK(ref::expm1(DD(0.0)).hi == 0.0);
}
