// Double-double elementary functions and constants (Unit 55, TVR-1, TC1). Sources and accuracy are in
// the header.
#include "NiTCAD/verify/reference/double_double.hpp"

#include <cmath>
#include <limits>

namespace NiTCAD::verify::reference {
namespace {

// The series below stop when a term is below 2^-110 of the sum (a sixteenth of u_dd).
constexpr double series_stop = 0x1p-110;

// sum_{k >= 0} sign^k / ((2k + 1) m^(2k + 1)): atan(1/m) for sign = -1, atanh(1/m) = acoth(m) for
// sign = +1. m^2 is exact in a double for every m used here (m < 2^26).
DD inverse_series(double m, bool alternating) {
    const DD m2 = two_prod(m, m);
    DD power = DD(1.0) / DD(m);
    DD sum;
    for (int k = 0; k < 200; ++k) {
        const DD term = power / DD(2.0 * k + 1.0);
        sum = (alternating && k % 2 == 1) ? sum - term : sum + term;
        if (std::abs(term.hi) < series_stop * std::abs(sum.hi)) break;
        power = power / m2;
    }
    return sum;
}

// e^r - 1 for |r| <= 2^-8 / 2 (or so) by Taylor: terms r^n / n! until below series_stop of the sum.
DD expm1_small(const DD& r) {
    DD term = r;
    DD sum = r;
    for (int n = 2; n <= 40; ++n) {
        term = term * r / DD(static_cast<double>(n));
        sum += term;
        if (std::abs(term.hi) <= series_stop * std::abs(sum.hi)) break;
    }
    return sum;
}

// (1 + s)^(2^8) - 1 from s = e^r - 1 by s -> 2 s + s^2, eight times.
DD square_up(DD s) {
    for (int i = 0; i < 8; ++i) s = ldexp(s, 1) + s * s;
    return s;
}

constexpr int reduction_log2 = 8;  // k = 2^8 = 256 (report section 5.1)

}  // namespace

DD sqrt(const DD& a) {
    if (a.hi == 0.0) return DD(0.0);
    if (a.hi < 0.0) return DD(std::numeric_limits<double>::quiet_NaN());
    DD x = 1.0 / std::sqrt(a.hi);
    for (int i = 0; i < 3; ++i) x = x + ldexp(x * (DD(1.0) - a * (x * x)), -1);
    return a * x;
}

DD exp(const DD& a) {
    if (a.hi > 709.78) return DD(std::numeric_limits<double>::infinity());
    if (a.hi < -745.13) return DD(0.0);
    const double m = std::nearbyint(a.hi / ln2().hi);
    const DD r = ldexp(a - ln2() * DD(m), -reduction_log2);  // |r| <= ln 2 / 512
    const DD s = square_up(expm1_small(r));
    return ldexp(s + DD(1.0), static_cast<int>(m));
}

DD expm1(const DD& a) {
    if (std::abs(a.hi) >= 0.5) return exp(a) - DD(1.0);
    return square_up(expm1_small(ldexp(a, -reduction_log2)));
}

DD log(const DD& a) {
    if (!(a.hi > 0.0)) return DD(std::numeric_limits<double>::quiet_NaN());
    const DD d = a - DD(1.0);
    if (std::abs(d.hi) < 0.125) {
        // log a = 2 atanh(z), z = (a - 1) / (a + 1), |z| < 1/15.
        const DD z = d / (a + DD(1.0));
        const DD z2 = z * z;
        DD power = z;
        DD sum = z;
        for (int k = 1; k < 100; ++k) {
            power = power * z2;
            const DD term = power / DD(2.0 * k + 1.0);
            sum += term;
            if (std::abs(term.hi) <= series_stop * std::abs(sum.hi)) break;
        }
        return ldexp(sum, 1);
    }
    DD x = std::log(a.hi);
    for (int i = 0; i < 3; ++i) x = x + a * exp(-x) - DD(1.0);
    return x;
}

DD pow(const DD& a, const DD& b) { return exp(b * log(a)); }

const DD& pi() {
    static const DD v = DD(48.0) * inverse_series(18.0, true) + DD(32.0) * inverse_series(57.0, true) -
                        DD(20.0) * inverse_series(239.0, true);
    return v;
}

const DD& ln2() {
    // acoth m = ln((m + 1) / (m - 1)) / 2, so 4 acoth 7 + 2 acoth 17 = 2 ln(4/3) + ln(9/8) = ln 2.
    static const DD v = DD(4.0) * inverse_series(7.0, false) + DD(2.0) * inverse_series(17.0, false);
    return v;
}

const DD& e() {
    static const DD v = exp(DD(1.0));
    return v;
}

const DD& sqrt2() {
    static const DD v = sqrt(DD(2.0));
    return v;
}

}  // namespace NiTCAD::verify::reference
