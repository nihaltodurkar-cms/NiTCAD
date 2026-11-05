// Double-double arithmetic for the independent reference layer (Unit 55, TVR-1, TC1).
//
// A value is the unevaluated sum hi + lo of two doubles with |lo| <= ulp(hi) / 2, about 106
// significant bits. The error-free transformations are Algorithms 3, 4 and 7 of Hida, Li and
// Bailey, "Library for Double-Double and Quad-Double Arithmetic" (8 May 2008); the double-double
// addition is the QD library's ieee_add (due to Briggs and Kahan), the product its FMA product and
// the quotient its accurate_div (QD source, include/qd/dd_inline.h); sqrt is the Newton iteration on
// x^-2 of the report's section 4.2, exp the reduction e^(kr + m ln2) = 2^m (e^r)^k with k = 256 of
// section 5.1, log the Newton iteration on e^x - a of section 5.2. Dekker (Numer. Math. 18, 1971)
// is the origin of the method; it was not read for this unit (paywalled), so the citations above
// are the sources (provenance PB, DECISIONS Unit 55 section 11).
//
// This library depends on the standard library only (TVD-19). It is not the test helper
// tests/physics/reference_arithmetic.hpp, which it neither includes nor changes (TVR1-Q1); the
// verify tests cross-check the two (TC2).
//
// Domain: double-double precision exists only where the lo part is a normal double, so every bound
// below is for values (operands, results and intermediates) of magnitude in [1e-280, 1e280]; below
// that the lo part goes subnormal and the precision falls towards a double's.
//
// Accuracy, in u_dd = 2^-104 (the bound each operation's L0 test asserts is in DECISIONS 11.3 as
// amended): +, -, *, / and sqrt a few u_dd relative; exp a few u_dd plus |x| / 8 u_dd (the
// double-double ln 2 used in the reduction has relative error about 2^-107, multiplied by m);
// log a few u_dd absolute for |x - 1| >= 1/8 and a few u_dd relative below (atanh series). Every
// operation uses two_sum and std::fma only, so the result does not depend on the compiler
// contracting or reassociating (MSVC /fp:precise does neither).
#pragma once

#include <cmath>

namespace NiTCAD::verify::reference {

// The double-double unit roundoff the tolerances are expressed in (DECISIONS 11.3).
inline constexpr double u_dd = 0x1p-104;

struct DD {
    double hi = 0.0;
    double lo = 0.0;
    constexpr DD() = default;
    constexpr DD(double h) : hi(h) {}  // NOLINT: a double is exactly a double-double
    constexpr DD(double h, double l) : hi(h), lo(l) {}
    [[nodiscard]] constexpr double to_double() const { return hi + lo; }
};

// Error-free transformations: the returned pair's hi + lo is exactly a + b (two_sum), a + b for
// |a| >= |b| (quick_two_sum), and a * b barring underflow (two_prod).
[[nodiscard]] inline DD two_sum(double a, double b) {
    const double s = a + b;
    const double v = s - a;
    return {s, (a - (s - v)) + (b - v)};
}
[[nodiscard]] inline DD quick_two_sum(double a, double b) {
    const double s = a + b;
    return {s, b - (s - a)};
}
[[nodiscard]] inline DD two_prod(double a, double b) {
    const double p = a * b;
    return {p, std::fma(a, b, -p)};
}

[[nodiscard]] inline DD operator-(const DD& a) { return {-a.hi, -a.lo}; }

// QD ieee_add.
[[nodiscard]] inline DD operator+(const DD& a, const DD& b) {
    DD s = two_sum(a.hi, b.hi);
    const DD t = two_sum(a.lo, b.lo);
    s = quick_two_sum(s.hi, s.lo + t.hi);
    return quick_two_sum(s.hi, s.lo + t.lo);
}
[[nodiscard]] inline DD operator-(const DD& a, const DD& b) { return a + (-b); }

// QD operator* with the FMA two_prod.
[[nodiscard]] inline DD operator*(const DD& a, const DD& b) {
    const DD p = two_prod(a.hi, b.hi);
    return quick_two_sum(p.hi, p.lo + (a.hi * b.lo + a.lo * b.hi));
}

// QD accurate_div: three quotient digits by long division.
[[nodiscard]] inline DD operator/(const DD& a, const DD& b) {
    const double q1 = a.hi / b.hi;
    DD r = a - DD(q1) * b;
    const double q2 = r.hi / b.hi;
    r = r - DD(q2) * b;
    const double q3 = r.hi / b.hi;
    return quick_two_sum(q1, q2) + DD(q3);
}

inline DD& operator+=(DD& a, const DD& b) { return a = a + b; }
inline DD& operator-=(DD& a, const DD& b) { return a = a - b; }
inline DD& operator*=(DD& a, const DD& b) { return a = a * b; }
inline DD& operator/=(DD& a, const DD& b) { return a = a / b; }

[[nodiscard]] inline bool operator==(const DD& a, const DD& b) { return a.hi == b.hi && a.lo == b.lo; }
[[nodiscard]] inline bool operator<(const DD& a, const DD& b) {
    return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}
[[nodiscard]] inline bool operator>(const DD& a, const DD& b) { return b < a; }

[[nodiscard]] inline DD abs(const DD& a) { return (a.hi < 0.0 || (a.hi == 0.0 && a.lo < 0.0)) ? -a : a; }
// Exact multiplication by 2^e (barring overflow and underflow of either part).
[[nodiscard]] inline DD ldexp(const DD& a, int e) { return {std::ldexp(a.hi, e), std::ldexp(a.lo, e)}; }

// Square root of a >= 0 (report section 4.2): x -> x + x (1 - a x^2) / 2 three times from the double
// 1 / sqrt(a.hi), then sqrt(a) = a x. Domain: a = 0, or a.hi in [1e-280, 1e280].
[[nodiscard]] DD sqrt(const DD& a);

// e^a (report section 5.1). Returns +inf above 709.78 and 0 below -745.13; within u_dd-level accuracy
// for |a| <= 670 (where the result's lo part stays a normal double).
[[nodiscard]] DD exp(const DD& a);

// e^a - 1 without cancellation near 0: for |a| < 1/2 the reduced series of exp with the squaring
// s -> 2 s + s^2 (which keeps e^r - 1 instead of e^r), otherwise exp(a) - 1.
[[nodiscard]] DD expm1(const DD& a);

// Natural logarithm of a > 0 (report section 5.2: three Newton steps on e^x - a from std::log(a.hi));
// for |a - 1| < 1/8, 2 atanh((a - 1) / (a + 1)) by its series instead, so the result is accurate
// relative to itself near its zero.
[[nodiscard]] DD log(const DD& a);

// a^b = e^(b log a) for a > 0; its condition number is |b log a|.
[[nodiscard]] DD pow(const DD& a, const DD& b);

// Constants computed here by series (no typed digits): pi by Gauss's formula
// pi = 48 atan(1/18) + 32 atan(1/57) - 20 atan(1/239); ln 2 = 4 acoth 7 + 2 acoth 17;
// e = exp(1); sqrt 2 = sqrt(2). The TC1 test checks each against OEIS's digits.
[[nodiscard]] const DD& pi();
[[nodiscard]] const DD& ln2();
[[nodiscard]] const DD& e();
[[nodiscard]] const DD& sqrt2();

}  // namespace NiTCAD::verify::reference
