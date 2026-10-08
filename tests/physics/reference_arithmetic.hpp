// Test helper: double-double arithmetic (a value as an unevaluated sum hi + lo of two doubles,
// about 32 significant digits; Dekker 1971, Bailey's QD library), so the tests compute their own
// high-precision references in C++ instead of taking literals from an outside tool. Elementary
// functions to about 1e-30 relative: exp (range reduction by ln 2 and 2^-10, Taylor series,
// squaring), log (Newton on exp), sqrt (one Newton step), pow through them. ln 2 and pi are
// computed here (series for atanh(1/3) and Machin's formula), not typed in. Every operation is
// exact rounding-error arithmetic (two_sum, and two_prod by std::fma), so the result does not
// depend on the compiler contracting or reassociating (MSVC /fp:precise does neither).
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace reference {

struct DD {
    double hi = 0.0, lo = 0.0;
    constexpr DD() = default;
    constexpr DD(double h) : hi(h) {}  // NOLINT: a double is exactly a DD
    constexpr DD(double h, double l) : hi(h), lo(l) {}
    [[nodiscard]] double value() const { return hi + lo; }
};

namespace detail {

inline DD two_sum(double a, double b) {
    const double s = a + b;
    const double bb = s - a;
    return {s, (a - (s - bb)) + (b - bb)};
}
inline DD quick_two_sum(double a, double b) {
    const double s = a + b;
    return {s, b - (s - a)};
}
inline DD two_prod(double a, double b) {
    const double p = a * b;
    return {p, std::fma(a, b, -p)};
}

}  // namespace detail

inline DD operator+(const DD& a, const DD& b) {
    DD s = detail::two_sum(a.hi, b.hi);
    const DD t = detail::two_sum(a.lo, b.lo);
    s.lo += t.hi;
    s = detail::quick_two_sum(s.hi, s.lo);
    s.lo += t.lo;
    return detail::quick_two_sum(s.hi, s.lo);
}
inline DD operator-(const DD& a) { return {-a.hi, -a.lo}; }
inline DD operator-(const DD& a, const DD& b) { return a + (-b); }
inline DD operator*(const DD& a, const DD& b) {
    DD p = detail::two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return detail::quick_two_sum(p.hi, p.lo);
}
inline DD operator/(const DD& a, const DD& b) {
    const double q1 = a.hi / b.hi;
    DD r = a - b * DD(q1);
    const double q2 = r.hi / b.hi;
    r = r - b * DD(q2);
    const double q3 = r.hi / b.hi;
    return detail::quick_two_sum(q1, q2) + DD(q3);
}
inline DD& operator+=(DD& a, const DD& b) { return a = a + b; }
inline DD& operator-=(DD& a, const DD& b) { return a = a - b; }
inline DD& operator*=(DD& a, const DD& b) { return a = a * b; }
inline bool operator<(const DD& a, const DD& b) {
    return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}
inline bool operator>(const DD& a, const DD& b) { return b < a; }
inline DD abs(const DD& a) { return a.hi < 0.0 ? -a : a; }
inline DD ldexp(const DD& a, int e) { return {std::ldexp(a.hi, e), std::ldexp(a.lo, e)}; }

inline DD sqrt(const DD& a) {
    if (a.hi <= 0.0) return {0.0};
    const DD y = std::sqrt(a.hi);
    return y + (a - y * y) / (DD(2.0) * y);
}

// atanh(1/m) = sum 1 / ((2k + 1) m^(2k+1)), for the constants.
inline DD atanh_inverse(double m) {
    DD sum = 0.0, power = DD(1.0) / DD(m);
    const DD m2 = DD(m) * DD(m);
    for (int k = 0; k < 200; ++k) {
        const DD term = power / DD(2.0 * k + 1.0);
        sum += term;
        if (std::abs(term.hi) < 1e-40 * std::abs(sum.hi)) break;
        power = power / m2;
    }
    return sum;
}
// atan(1/m) = sum (-1)^k / ((2k + 1) m^(2k+1)).
inline DD atan_inverse(double m) {
    DD sum = 0.0, power = DD(1.0) / DD(m);
    const DD m2 = DD(m) * DD(m);
    for (int k = 0; k < 200; ++k) {
        const DD term = power / DD(2.0 * k + 1.0);
        sum = (k % 2 == 0) ? sum + term : sum - term;
        if (std::abs(term.hi) < 1e-40 * std::abs(sum.hi)) break;
        power = power / m2;
    }
    return sum;
}
inline const DD& ln2() {
    static const DD v = DD(2.0) * atanh_inverse(3.0);  // ln 2 = 2 atanh(1/3)
    return v;
}
inline const DD& pi() {
    static const DD v = DD(16.0) * atan_inverse(5.0) - DD(4.0) * atan_inverse(239.0);  // Machin
    return v;
}

inline DD exp(const DD& a) {
    if (a.hi > 709.7) return {std::numeric_limits<double>::infinity()};
    if (a.hi < -745.2) return {0.0};
    const double k = std::nearbyint(a.hi / ln2().hi);
    const DD r = ldexp(a - ln2() * DD(k), -10);  // |r| <= ln2 / 2048
    // e^r - 1 by Taylor; |r|^n / n! < 1e-36 by n = 9.
    DD term = r, sum = r;
    for (int n = 2; n <= 12; ++n) {
        term = term * r / DD(static_cast<double>(n));
        sum += term;
    }
    // (1 + s)^(2^10) by s -> 2 s + s^2, ten times (keeps e^r - 1 without cancellation).
    for (int i = 0; i < 10; ++i) sum = DD(2.0) * sum + sum * sum;
    return ldexp(sum + DD(1.0), static_cast<int>(k));
}
// e^x - 1, accurate near 0.
inline DD expm1(const DD& a) {
    if (std::abs(a.hi) > 0.5) return exp(a) - DD(1.0);
    DD term = a, sum = a;
    for (int n = 2; n <= 40; ++n) {
        term = term * a / DD(static_cast<double>(n));
        sum += term;
        if (std::abs(term.hi) < 1e-36 * std::abs(sum.hi)) break;
    }
    return sum;
}
inline DD log(const DD& a) {
    DD y = std::log(a.hi);
    for (int i = 0; i < 3; ++i) y = y + a * exp(-y) - DD(1.0);
    return y;
}
inline DD pow(const DD& a, const DD& b) { return exp(b * log(a)); }

// Gauss-Legendre nodes and weights of order n on [-1, 1], computed by Newton on P_n in DD
// arithmetic (the n / 2 positive nodes, ascending, and their weights).
struct GaussLegendre {
    std::vector<DD> x, w;
};
inline GaussLegendre gauss_legendre(int n) {
    GaussLegendre g;
    for (int i = n / 2; i >= 1; --i) {
        DD x = std::cos(pi().hi * (i - 0.25) / (n + 0.5));
        DD dp = 0.0;
        for (int it = 0; it < 100; ++it) {
            DD p0 = 1.0, p1 = x;  // P_0, P_1
            for (int k = 2; k <= n; ++k) {
                const DD p2 = (DD(2.0 * k - 1.0) * x * p1 - DD(k - 1.0) * p0) / DD(k);
                p0 = p1;
                p1 = p2;
            }
            dp = DD(n) * (x * p1 - p0) / (x * x - DD(1.0));  // P_n'
            const DD step = p1 / dp;
            x -= step;
            if (std::abs(step.hi) < 1e-34) break;
        }
        DD p0 = 1.0, p1 = x;
        for (int k = 2; k <= n; ++k) {
            const DD p2 = (DD(2.0 * k - 1.0) * x * p1 - DD(k - 1.0) * p0) / DD(k);
            p0 = p1;
            p1 = p2;
        }
        dp = DD(n) * (x * p1 - p0) / (x * x - DD(1.0));
        g.x.push_back(x);
        g.w.push_back(DD(2.0) / ((DD(1.0) - x * x) * dp * dp));
    }
    return g;
}

// The Bernoulli numbers B_2 .. B_20 (exact rationals).
inline DD bernoulli_number(int k2) {
    static constexpr double num[] = {1, -1, 1, -1, 5, -691, 7, -3617, 43867, -174611};
    static constexpr double den[] = {6, 30, 42, 30, 66, 2730, 6, 510, 798, 330};
    return DD(num[k2 / 2 - 1]) / DD(den[k2 / 2 - 1]);
}

// zeta(s) for real s > 1 by Euler-Maclaurin: sum_{k<N} k^-s + N^(1-s)/(s-1) + N^-s / 2
// + sum_j B_2j / (2j)! s (s+1) ... (s+2j-2) N^(-s-2j+1), N = 20, j up to 10.
inline DD zeta(const DD& s) {
    const int N = 20;
    DD sum = 0.0;
    for (int k = 1; k < N; ++k) sum += exp(-s * log(DD(k)));
    const DD lnN = log(DD(N));
    sum += exp((DD(1.0) - s) * lnN) / (s - DD(1.0));
    sum += DD(0.5) * exp(-s * lnN);
    DD rising = s;  // s (s+1) ... (s+2j-2)
    DD factorial = 2.0;  // (2j)!
    for (int j = 1; j <= 10; ++j) {
        sum += bernoulli_number(2 * j) / factorial * rising *
               exp((-s - DD(2.0 * j - 1.0)) * lnN);
        rising = rising * (s + DD(2.0 * j - 1.0)) * (s + DD(2.0 * j));
        factorial = factorial * DD(2.0 * j + 1.0) * DD(2.0 * j + 2.0);
    }
    return sum;
}

}  // namespace reference
