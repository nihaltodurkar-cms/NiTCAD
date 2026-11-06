// Double-double special functions (Unit 55, TVR-3 step 2; DECISIONS 15.11). Sources in the header.
#include "NiTCAD/verify/reference/special_functions.hpp"

#include <cmath>
#include <stdexcept>

namespace NiTCAD::verify::reference {

namespace {

// Series stop when a term falls below 2^-110 of the sum (a sixteenth of u_dd), as in double_double.cpp.
constexpr double series_stop = 0x1p-110;
constexpr int max_terms = 10000;

bool small_term(const DD& term, const DD& sum) {
    return std::abs(term.hi) <= series_stop * std::abs(sum.hi);
}

// atanh(s) = sum s^(2k+1) / (2k+1), |s| small.
DD atanh_series(const DD& s) {
    const DD s2 = s * s;
    DD power = s;
    DD sum = s;
    for (int k = 1; k < max_terms; ++k) {
        power = power * s2;
        const DD term = power / DD(2.0 * k + 1.0);
        sum = sum + term;
        if (small_term(term, sum)) return sum;
    }
    throw std::logic_error("atanh_series did not converge");
}

DD half_pi() { return ldexp(pi(), -1); }
DD pi_squared_over_6() { return pi() * pi() / DD(6.0); }

// erf for 0 <= x < 3 by DLMF 7.6.2 (positive terms).
DD erf_series(const DD& x) {
    const DD x2 = x * x;
    const DD two_x2 = ldexp(x2, 1);
    DD term = x;
    DD sum = x;
    for (int n = 1; n < max_terms; ++n) {
        term = term * two_x2 / DD(2.0 * n + 1.0);
        sum = sum + term;
        if (small_term(term, sum)) {
            return DD(2.0) / sqrt(pi()) * exp(-x2) * sum;
        }
    }
    throw std::logic_error("erf_series did not converge");
}

// The DLMF 7.9.2 continued fraction 2x / (2x^2 + 1 - 1*2 / (2x^2 + 5 - 3*4 / (2x^2 + 9 - ...))), evaluated
// from depth `depth` upwards.
DD erfc_fraction(const DD& x, int depth) {
    const DD b0 = ldexp(x * x, 1);  // 2x^2
    DD t = b0 + DD(4.0 * depth + 1.0);
    for (int k = depth; k >= 1; --k) {
        const double a = (2.0 * k - 1.0) * (2.0 * k);
        t = b0 + DD(4.0 * (k - 1) + 1.0) - DD(a) / t;
    }
    return ldexp(x, 1) / t;
}

// erfc for x >= 3: the fraction at doubling depths until two agree to series_stop.
DD erfc_large(const DD& x) {
    if (x.hi > 26.5) return DD(0.0);
    DD previous = erfc_fraction(x, 16);
    for (int depth = 32; depth <= 1 << 16; depth *= 2) {
        const DD current = erfc_fraction(x, depth);
        if (std::abs((current - previous).hi) <= series_stop * std::abs(current.hi)) {
            return current * exp(-(x * x)) / sqrt(pi());
        }
        previous = current;
    }
    throw std::logic_error("erfc continued fraction did not converge");
}

// Li2 by its series (DLMF 25.12.1) for |x| <= 1/2.
DD li2_series(const DD& x) {
    DD power = x;
    DD sum = x;
    for (int k = 2; k < max_terms; ++k) {
        power = power * x;
        const DD term = power / DD(static_cast<double>(k) * k);
        sum = sum + term;
        if (small_term(term, sum)) return sum;
    }
    throw std::logic_error("li2_series did not converge");
}

}  // namespace

DD log1p(const DD& x) {
    if (!(x.hi > -1.0)) throw std::domain_error("log1p: x <= -1");
    if (std::abs(x.hi) < 0.125) {
        if (x.hi == 0.0 && x.lo == 0.0) return DD(0.0);
        return ldexp(atanh_series(x / (DD(2.0) + x)), 1);
    }
    return log(DD(1.0) + x);
}

DD log1p_exp(const DD& x) {
    if (x.hi > 0.0) return x + log1p(exp(-x));
    if (x.hi < -745.0) return DD(0.0);
    return log1p(exp(x));
}

DD tanh(const DD& x) {
    const bool negative = x.hi < 0.0;
    const DD ax = negative ? -x : x;
    if (ax.hi > 40.0) return DD(negative ? -1.0 : 1.0);
    if (ax.hi == 0.0) return DD(0.0);
    const DD t = expm1(ldexp(ax, 1));
    const DD r = t / (t + DD(2.0));
    return negative ? -r : r;
}

DD erf(const DD& x) {
    const bool negative = x.hi < 0.0;
    const DD ax = negative ? -x : x;
    if (ax.hi == 0.0) return DD(0.0);
    const DD r = ax.hi < 3.0 ? erf_series(ax) : DD(1.0) - erfc_large(ax);
    return negative ? -r : r;
}

DD erfc(const DD& x) {
    if (x.hi >= 3.0) return erfc_large(x);
    if (x.hi >= 0.0) return DD(1.0) - erf_series(x);
    return DD(1.0) + erf(-x);
}

DD sin(const DD& x) {
    if (!(std::abs(x.hi) <= 2.0)) throw std::domain_error("sin: |x| > 2");
    const DD x2 = x * x;
    DD term = x;
    DD sum = x;
    if (x.hi == 0.0) return x;
    for (int k = 1; k < max_terms; ++k) {
        term = -(term * x2) / DD((2.0 * k) * (2.0 * k + 1.0));
        sum = sum + term;
        if (small_term(term, sum)) return sum;
    }
    throw std::logic_error("sin series did not converge");
}

DD cos(const DD& x) {
    if (!(std::abs(x.hi) <= 2.0)) throw std::domain_error("cos: |x| > 2");
    const DD x2 = x * x;
    DD term(1.0);
    DD sum(1.0);
    for (int k = 1; k < max_terms; ++k) {
        term = -(term * x2) / DD((2.0 * k - 1.0) * (2.0 * k));
        sum = sum + term;
        if (small_term(term, sum)) return sum;
    }
    throw std::logic_error("cos series did not converge");
}

DD atan(const DD& x) {
    const bool negative = x.hi < 0.0;
    const DD ax = negative ? -x : x;
    if (ax.hi == 0.0) return DD(0.0);
    DD r;
    if (ax.hi > 1.0) {
        r = half_pi() - atan(DD(1.0) / ax);
    } else {
        // Newton on f(y) = sin y - x cos y (zero at tan y = x), f'(y) = cos y + x sin y; |y| <= pi/4.
        DD y(std::atan(ax.hi));
        for (int i = 0; i < 3; ++i) {
            const DD s = sin(y);
            const DD c = cos(y);
            y = y - (s - ax * c) / (c + ax * s);
        }
        r = y;
    }
    return negative ? -r : r;
}

DD atan2(const DD& y, const DD& x) {
    if (x.hi == 0.0 && y.hi == 0.0) throw std::domain_error("atan2(0, 0)");
    if (x.hi == 0.0) return y.hi > 0.0 ? half_pi() : -half_pi();
    const DD base = atan(y / x);
    if (x.hi > 0.0) return base;
    return y.hi >= 0.0 ? base + pi() : base - pi();
}

DD asin(const DD& x) {
    if (!(std::abs(x.hi) <= 1.0)) throw std::domain_error("asin: |x| > 1");
    const DD c = sqrt((DD(1.0) - x) * (DD(1.0) + x));
    if (c.hi == 0.0) return x.hi > 0.0 ? half_pi() : -half_pi();
    return atan2(x, c);
}

DD li2(const DD& x) {
    if (!(x.hi <= 1.0) || (x.hi == 1.0 && x.lo > 0.0)) throw std::domain_error("li2: x > 1");
    if (x.hi == 1.0 && x.lo == 0.0) return pi_squared_over_6();
    if (x.hi == 0.0) return DD(0.0);
    if (x.hi < -1.0) {
        const DD l = log(-x);
        return -pi_squared_over_6() - ldexp(l * l, -1) - li2(DD(1.0) / x);
    }
    if (x.hi < -0.5) {
        const DD l = log1p(-x);  // ln(1 - x)
        return -li2_series(x / (x - DD(1.0))) - ldexp(l * l, -1);
    }
    if (x.hi <= 0.5) return li2_series(x);
    const DD one_minus = DD(1.0) - x;
    return pi_squared_over_6() - log(x) * log(one_minus) - li2_series(one_minus);
}

}  // namespace NiTCAD::verify::reference
