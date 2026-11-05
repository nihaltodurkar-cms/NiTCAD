// tanh-sinh (double-exponential) quadrature in double-double (Unit 55, TVR-1, TC4).
//
// Source: Takahasi and Mori, Publ. RIMS 9 (1974) 721-741 (paywalled; graded PB under TVR1-Q2 from two open
// reproductions that agree on everything used here): Bailey, Jeyabalan and Li, "A comparison of three
// high-precision quadrature schemes", Experimental Mathematics 14 (2005) 317-329, sections 3-5 (QUADTS), and
// Bailey, "Tanh-sinh high-precision quadrature" (19 Jan 2006), sections 2-4 (both read 2026-11-04).
//
// The method as both describe it:
//   int_{-1}^{1} f(x) dx ~ h sum_j w_j f(x_j),  x_j = tanh(u2), w_j = u1 / cosh^2(u2),
//   u1 = pi/2 cosh(hj), u2 = pi/2 sinh(hj);
// level k uses h = 2^-k, starting at level 1; the pairs of a level are those of the finest level taken with a
// stride, and each level after the first evaluates f only at its odd-indexed abscissae (the even ones give the
// previous level's sum). 1 - x_j is stored, as y_j = 1 / (e^u2 cosh u2), so that abscissae near an endpoint keep
// their relative accuracy; [a, b] is mapped linearly, the points near a as a + r y_j and near b as b - r y_j
// (r = (b - a) / 2), in double-double (the sources' "secondary precision" scaling). Pairs are generated for
// j >= 0 until w_j < eps2 = eps1^3, eps1 = u_dd here (the sources' 10^-p). The secondary epsilon is the
// note's eps1^n with n = 3 (it allows any n >= 2; BJL use n = 2): with n = 2 the endpoint term d4 below stays
// at ~7e-29 for x^-1/2 ln x, above a 16 u_dd target (owner ruling 2026-11-04, after measurement).
//
// The error estimate is the sources' heuristic: with S_k the level-k result, E_n = 1 for n <= 2,
// 0 if S_n = S_{n-1}, else 10^d, d = min(0, max(d1^2/d2, 2 d1, d3, d4)), where d1 = log10|S_n - S_{n-1}|,
// d2 = log10|S_n - S_{n-2}|, d3 = log10(eps1 max_j |w_j f(x_j)|), d4 = log10 max(|w_l f(x_l)|, |w_r f(x_r)|)
// over the abscissae nearest the two endpoints. The w_j f(x_j) are the summands without the factor h, as in
// the sources' formula (1), with the interval's r; d is not rounded to an integer (the sources allow, not
// require, it). E is absolute, and capped at 1: the scheme is meant for integrals of order one, as BJL's are.
// One deviation (owner ruling 2026-11-04, after measurement): for n >= 3, E_n is raised to at least
// 8 u_dd h sum_j |w_j f(x_j)|, the 8 u_dd TC1 guarantees for the library's functions times the computed
// int |f|. The sources' floor d3 assumes f exact to eps1; with f's own rounding (a few u_dd for log and
// sqrt) their E fell below the true error for ln x and x^-1/2 ln x. So no epsrel below ~8 u_dd is reachable.
// The iteration stops at the first level n >= 3 with E_n <= max(epsabs, epsrel |S_n|).
//
// Infinite intervals are not built in: as BJL do (problems 11-13), the caller maps them to a finite one.
// Domain: finite a, b; f finite at every abscissa (f is never evaluated at a or b); the endpoint distances
// r y_j reach ~1e-63, which double-double holds with full relative precision.
#pragma once

#include "NiTCAD/verify/reference/double_double.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace NiTCAD::verify::reference {

// The finest level's pairs, j = 0 .. size - 1 at h = 2^-max_level; y[j] = 1 - x_j, w[j] = g'(hj).
struct TanhSinhTable {
    int max_level = 0;
    std::vector<DD> y;
    std::vector<DD> w;
};

// Computed once (thread-safe static initialisation), max_level = 10.
[[nodiscard]] const TanhSinhTable& tanh_sinh_table();

enum class TanhSinhStatus {
    ok,             // E_n <= max(epsabs, epsrel |S_n|) at some level n >= 3
    limit_reached,  // the last level was completed without that
    bad_integrand,  // f returned a non-finite value
    invalid_input,  // non-finite a or b, or a negative or all-zero tolerance
};

struct TanhSinhResult {
    DD value;
    double abs_error = 1.0;  // the sources' E_n
    int levels = 0;
    int evaluations = 0;
    TanhSinhStatus status = TanhSinhStatus::ok;
};

template <class F>
[[nodiscard]] TanhSinhResult tanh_sinh(F&& f, const DD& a, const DD& b, double epsabs, double epsrel,
                                       int max_level = 10) {
    TanhSinhResult out;
    const TanhSinhTable& tab = tanh_sinh_table();
    if (!std::isfinite(a.hi) || !std::isfinite(b.hi) || !(epsabs >= 0.0) || !(epsrel >= 0.0) ||
        (epsabs == 0.0 && epsrel == 0.0) || max_level < 1 || max_level > tab.max_level) {
        out.status = TanhSinhStatus::invalid_input;
        return out;
    }
    const DD c = ldexp(a + b, -1);
    const DD r = ldexp(b - a, -1);
    const std::size_t last = tab.y.size() - 1;
    DD sum;                  // sum of r w_j f(x_j) over the abscissae evaluated so far
    double max_term = 0.0;   // max |r w_j f(x_j)|
    double abs_sum = 0.0;    // sum of |r w_j f(x_j)|
    std::size_t far_j = 0;   // the largest j evaluated, and its two terms
    double far_left = 0.0, far_right = 0.0;
    std::vector<DD> s;       // S_1, S_2, ...
    auto add = [&](const DD& term) {
        sum += term;
        max_term = std::max(max_term, std::abs(term.hi));
        abs_sum += std::abs(term.hi);
        return std::isfinite(term.hi);
    };
    for (int k = 1; k <= max_level; ++k) {
        const std::size_t stride = std::size_t{1} << (tab.max_level - k);
        const std::size_t first = k == 1 ? 0 : stride;
        const std::size_t step = k == 1 ? stride : 2 * stride;
        for (std::size_t j = first; j <= last; j += step) {
            if (j == 0) {
                ++out.evaluations;
                if (!add(r * tab.w[0] * f(c))) {
                    out.status = TanhSinhStatus::bad_integrand;
                    return out;
                }
                continue;
            }
            const DD dx = r * tab.y[j];
            const DD left = r * tab.w[j] * f(a + dx);
            const DD right = r * tab.w[j] * f(b - dx);
            out.evaluations += 2;
            if (!add(left) || !add(right)) {
                out.status = TanhSinhStatus::bad_integrand;
                return out;
            }
            if (j > far_j) {
                far_j = j;
                far_left = left.hi;
                far_right = right.hi;
            }
        }
        s.push_back(ldexp(sum, -k));
        out.levels = k;
        out.value = s.back();
        const std::size_t n = s.size();
        if (n <= 2) {
            out.abs_error = 1.0;
        } else if (s[n - 1] == s[n - 2]) {
            out.abs_error = 0.0;
        } else {
            const double inf = std::numeric_limits<double>::infinity();
            const double d1 = std::log10(std::abs((s[n - 1] - s[n - 2]).to_double()));
            const double d2 = std::log10(std::abs((s[n - 1] - s[n - 3]).to_double()));
            const double d3 = std::log10(u_dd * max_term);
            const double d4 = std::log10(std::max(std::abs(far_left), std::abs(far_right)));
            const double t1 = d2 == -inf ? -inf : (d2 == 0.0 ? inf : d1 * d1 / d2);
            const double d = std::min(0.0, std::max({t1, 2.0 * d1, d3, d4}));
            out.abs_error = std::pow(10.0, d);
        }
        if (n >= 3) out.abs_error = std::max(out.abs_error, 8.0 * u_dd * std::ldexp(abs_sum, -k));  // the deviation
        if (n >= 3 && out.abs_error <= std::max(epsabs, epsrel * std::abs(out.value.to_double()))) return out;
    }
    out.status = TanhSinhStatus::limit_reached;
    return out;
}

}  // namespace NiTCAD::verify::reference
