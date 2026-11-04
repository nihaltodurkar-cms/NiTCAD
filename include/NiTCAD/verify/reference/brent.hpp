// Bracketed root finding by Brent's procedure zero, in double-double (Unit 55, TVR-1, TC5).
//
// Source: R. P. Brent, Algorithms for Minimization without Derivatives (Prentice-Hall, 1973), chapter 4, read in
// the author's scan (maths-people.anu.edu.au/~brent/pd/rpb011i.pdf, 2026-11-04): the ALGOL 60 procedure zero
// (section 6, pp. 58-59), the tolerance (2.9) delta = 2 eps |b| + t, the error bound (2.18) |zeta^ - zeta| <=
// 6 eps |zeta| + 2t, and the bound (3.2)-(3.4) on the number of function evaluations, (k + 1)^2 - 2 with
// k = ceil(log2((b - a) / delta_m)), delta_m the minimum over [a, b] of 2 macheps |x| + t (k <= 0: two
// evaluations). netlib's go/zeroin.f, "a slightly modified translation" of zero, was compared and is not
// followed where it differs (it uses 2 eps |b| + tol / 2).
//
// The translation is line by line, with the ALGOL's labels int and ext as the two loop entries, in double-double
// with macheps = u_dd (the book's "relative machine precision" for rounded arithmetic; it says macheps "should not
// be decreased"). The only additions: refusal of an interval whose end values have the same sign (the procedure
// "assumes" they differ; zeroin.f checks and aborts), and a stop after twice the book's evaluation bound, so that a
// count above the bound is reported rather than looping.
#pragma once

#include "NiTCAD/verify/reference/double_double.hpp"

#include <algorithm>
#include <cmath>

namespace NiTCAD::verify::reference {

// (3.4): the most function evaluations zero can take on [a, b] with absolute tolerance t > 0.
[[nodiscard]] inline long long brent_evaluation_bound(const DD& a, const DD& b, double t, double macheps = u_dd) {
    const double lo = std::min(a.to_double(), b.to_double());
    const double hi = std::max(a.to_double(), b.to_double());
    const double delta_m = (lo <= 0.0 && hi >= 0.0) ? t : 2.0 * macheps * std::min(std::abs(lo), std::abs(hi)) + t;
    const double k = std::ceil(std::log2(std::abs((b - a).to_double()) / delta_m));
    if (k <= 0.0) return 2;
    const long long kk = static_cast<long long>(k);
    return (kk + 1) * (kk + 1) - 2;
}

enum class BrentStatus {
    ok,             // terminated by zero's own test: |m| <= tol or f(b) = 0
    not_bracketed,  // f(a) and f(b) are nonzero with the same sign
    bound_exceeded, // more than twice (3.4)'s evaluations (never expected; reported, not looped)
    invalid_input,  // t <= 0 or a non-finite end point
};

struct BrentResult {
    DD root;
    long long evaluations = 0;
    BrentStatus status = BrentStatus::ok;
};

template <class F>
[[nodiscard]] BrentResult brent_zero(F&& f, DD a, DD b, double t, double macheps = u_dd) {
    BrentResult out;
    if (!(t > 0.0) || !std::isfinite(a.hi) || !std::isfinite(b.hi)) {
        out.status = BrentStatus::invalid_input;
        return out;
    }
    const long long cap = 2 * brent_evaluation_bound(a, b, t, macheps);
    DD fa = f(a), fb = f(b);
    out.evaluations = 2;
    if (fa.hi != 0.0 && fb.hi != 0.0 && (fa.hi > 0.0) == (fb.hi > 0.0)) {
        out.status = BrentStatus::not_bracketed;
        return out;
    }
    DD c, fc, d, e;
    bool to_int = true;  // which label the loop enters: int (c := a ...) or ext
    for (;;) {
        if (to_int) {  // int:
            c = a;
            fc = fa;
            d = e = b - a;
        }
        // ext:
        if (abs(fc) < abs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        const DD tol = DD(2.0 * macheps) * abs(b) + DD(t);
        const DD m = ldexp(c - b, -1);
        if (!(tol < abs(m) && fb.hi != 0.0)) break;
        // See if a bisection is forced.
        if (abs(e) < tol || !(abs(fb) < abs(fa))) {
            d = e = m;
        } else {
            DD s = fb / fa, p, q;
            if (a == c) {  // linear interpolation
                p = ldexp(m * s, 1);
                q = DD(1.0) - s;
            } else {  // inverse quadratic interpolation
                q = fa / fc;
                const DD r = fb / fc;
                p = s * (ldexp(m * q * (q - r), 1) - (b - a) * (r - DD(1.0)));
                q = (q - DD(1.0)) * (r - DD(1.0)) * (s - DD(1.0));
            }
            if (p.hi > 0.0)
                q = -q;
            else
                p = -p;
            s = e;
            e = d;
            if (ldexp(p, 1) < DD(3.0) * m * q - abs(tol * q) && p < abs(ldexp(s * q, -1)))
                d = p / q;
            else
                d = e = m;
        }
        a = b;
        fa = fb;
        b = b + (tol < abs(d) ? d : (m.hi > 0.0 ? tol : -tol));
        fb = f(b);
        ++out.evaluations;
        if (out.evaluations > cap) {
            out.status = BrentStatus::bound_exceeded;
            break;
        }
        to_int = (fb.hi > 0.0) == (fc.hi > 0.0);
    }
    out.root = b;
    return out;
}

}  // namespace NiTCAD::verify::reference
