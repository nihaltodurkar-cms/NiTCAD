// Gauss-Kronrod 7/15 quadrature in double-double, with QUADPACK's error estimate and adaptive bisection
// (Unit 55, TVR-1, TC3).
//
// Source: QUADPACK (Piessens, de Doncker-Kapenga, Ueberhuber and Kahaner, 1983), netlib's public dqk15.f and
// dqage.f (read 2026-11-04). The rule's nodes and weights are *computed* here, not typed: the 7 Gauss nodes
// by Newton on P_7, their weights 2 / ((1 - x^2) P_7'(x)^2); the 8 Kronrod nodes as the roots of the
// Stieltjes polynomial E_8 = P_8 + c_6 P_6 + c_4 P_4 + c_2 P_2 + c_0 (even, since P_7 is odd), fixed by
// int P_7 E_8 P_k = 0 for k = 1, 3, 5, 7 (the orthogonality that defines it), each root bracketed between
// consecutive Gauss nodes (they interlace) and polished by Newton; the 15 Kronrod weights by making the rule
// exact for P_0 .. P_14 (the even ones; the odd ones hold by symmetry). The test compares them with the
// 33 digits dqk15.f prints (Fullerton's 80-digit values).
//
// The error estimate is dqk15's: |(K15 - G7) h| scaled by resasc min(1, (200 |K15 - G7| h / resasc)^1.5),
// floored at 50 eps resabs, with eps = u_dd = 2^-104 here (QUADPACK's eps is the working precision's). The
// adaptive driver is dqage's (bisection of the interval of largest error; roundoff counters iroff1 >= 6 or
// iroff2 >= 20; stop when the summed error <= max(epsabs, epsrel |result|)), except that the interval of
// largest error is found by a linear search instead of dqpsrt's partially sorted list (the same choice
// whenever the list is not truncated near the limit), and that the running error sum is kept in double-double.
//
// The integrand is any callable DD -> DD. Domain: finite a < b or a > b, integrand finite at the 15 interior
// nodes of every subinterval (the endpoints are never evaluated).
#pragma once

#include "NiTCAD/verify/reference/double_double.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace NiTCAD::verify::reference {

// Positive abscissae on [-1, 1] in dqk15's order: xgk[0] > xgk[1] > ... > xgk[7] = 0; xgk[1], xgk[3], xgk[5]
// (and 0) are the Gauss nodes. wgk are the Kronrod weights, wg the Gauss weights of xgk[1], [3], [5], [7].
struct KronrodRule15 {
    std::array<DD, 8> xgk;
    std::array<DD, 8> wgk;
    std::array<DD, 4> wg;
};

// Computed once (thread-safe static initialisation).
[[nodiscard]] const KronrodRule15& kronrod15();

struct RuleResult {
    DD value;               // the Kronrod estimate
    double abs_error = 0.0;  // dqk15's estimate of |I - value|
    double resabs = 0.0;     // approximation of the integral of |f|
    double resasc = 0.0;     // approximation of the integral of |f - I / (b - a)|
};

// dqk15 over [a, b].
template <class F>
[[nodiscard]] RuleResult qk15(F&& f, const DD& a, const DD& b) {
    const KronrodRule15& r = kronrod15();
    const DD centr = ldexp(a + b, -1);
    const DD hlgth = ldexp(b - a, -1);
    const double dhlgth = std::abs(hlgth.to_double());
    std::array<DD, 7> fv1{}, fv2{};
    const DD fc = f(centr);
    DD resg = fc * r.wg[3];
    DD resk = fc * r.wgk[7];
    DD resabs = abs(resk);
    for (int j = 0; j < 3; ++j) {
        const int jtw = 2 * j + 1;  // Gauss nodes
        const DD absc = hlgth * r.xgk[jtw];
        const DD f1 = f(centr - absc);
        const DD f2 = f(centr + absc);
        fv1[jtw] = f1;
        fv2[jtw] = f2;
        resg += r.wg[j] * (f1 + f2);
        resk += r.wgk[jtw] * (f1 + f2);
        resabs += r.wgk[jtw] * (abs(f1) + abs(f2));
    }
    for (int j = 0; j < 4; ++j) {
        const int jtwm1 = 2 * j;  // Kronrod nodes
        const DD absc = hlgth * r.xgk[jtwm1];
        const DD f1 = f(centr - absc);
        const DD f2 = f(centr + absc);
        fv1[jtwm1] = f1;
        fv2[jtwm1] = f2;
        resk += r.wgk[jtwm1] * (f1 + f2);
        resabs += r.wgk[jtwm1] * (abs(f1) + abs(f2));
    }
    const DD reskh = ldexp(resk, -1);
    DD resasc = r.wgk[7] * abs(fc - reskh);
    for (int j = 0; j < 7; ++j) resasc += r.wgk[j] * (abs(fv1[j] - reskh) + abs(fv2[j] - reskh));
    RuleResult out;
    out.value = resk * hlgth;
    out.resabs = resabs.to_double() * dhlgth;
    out.resasc = resasc.to_double() * dhlgth;
    double abserr = std::abs(((resk - resg) * hlgth).to_double());
    if (out.resasc != 0.0 && abserr != 0.0) abserr = out.resasc * std::min(1.0, std::pow(200.0 * abserr / out.resasc, 1.5));
    const double uflow = std::numeric_limits<double>::min();
    if (out.resabs > uflow / (50.0 * u_dd)) abserr = std::max(50.0 * u_dd * out.resabs, abserr);
    out.abs_error = abserr;
    return out;
}

enum class QuadratureStatus {
    ok,             // the summed error estimate is within max(epsabs, epsrel |result|)
    limit_reached,  // dqage ier = 1
    roundoff,       // ier = 2
    bad_integrand,  // ier = 3 (an interval shrank to the precision of its endpoints)
    invalid_input,  // ier = 6
};

struct AdaptiveResult {
    DD value;
    double abs_error = 0.0;
    int intervals = 0;
    int evaluations = 0;
    QuadratureStatus status = QuadratureStatus::ok;
};

// dqage with the 7/15 rule (key = 1).
template <class F>
[[nodiscard]] AdaptiveResult qag15(F&& f, const DD& a, const DD& b, double epsabs, double epsrel, int limit = 2000) {
    AdaptiveResult out;
    if (limit < 1 || (epsabs <= 0.0 && epsrel < std::max(50.0 * u_dd, 0.5e-28))) {
        out.status = QuadratureStatus::invalid_input;
        return out;
    }
    const double uflow = std::numeric_limits<double>::min();
    struct Interval {
        DD a, b, area;
        double error;
    };
    std::vector<Interval> list;
    list.reserve(static_cast<std::size_t>(limit));
    const RuleResult first = qk15(f, a, b);
    out.evaluations = 15;
    list.push_back({a, b, first.value, first.abs_error});
    double errbnd = std::max(epsabs, epsrel * std::abs(first.value.to_double()));
    if (first.abs_error <= 50.0 * u_dd * first.resabs && first.abs_error > errbnd) out.status = QuadratureStatus::roundoff;
    if (limit == 1) out.status = QuadratureStatus::limit_reached;
    if (out.status != QuadratureStatus::ok || (first.abs_error <= errbnd && first.abs_error != first.resasc) ||
        first.abs_error == 0.0) {
        out.value = first.value;
        out.abs_error = first.abs_error;
        out.intervals = 1;
        return out;
    }
    DD area = first.value;
    DD errsum(first.abs_error);  // in double-double: in double its drift (~2^-53 of the first, largest
                                 // errors) would swamp a 1e-24 target and can go negative
    int iroff1 = 0, iroff2 = 0;
    for (int last = 2; last <= limit; ++last) {
        std::size_t maxerr = 0;
        for (std::size_t i = 1; i < list.size(); ++i)
            if (list[i].error > list[maxerr].error) maxerr = i;
        const Interval parent = list[maxerr];
        const DD b1 = ldexp(parent.a + parent.b, -1);
        const RuleResult r1 = qk15(f, parent.a, b1);
        const RuleResult r2 = qk15(f, b1, parent.b);
        out.evaluations += 30;
        const DD area12 = r1.value + r2.value;
        const double erro12 = r1.abs_error + r2.abs_error;
        errsum = errsum + DD(r1.abs_error) + DD(r2.abs_error) - DD(parent.error);
        area = area + area12 - parent.area;
        if (r1.resasc != r1.abs_error && r2.resasc != r2.abs_error) {
            if (std::abs((parent.area - area12).to_double()) <= 1e-5 * std::abs(area12.to_double()) &&
                erro12 >= 0.99 * parent.error)
                ++iroff1;
            if (last > 10 && erro12 > parent.error) ++iroff2;
        }
        list[maxerr] = {parent.a, b1, r1.value, r1.abs_error};
        list.push_back({b1, parent.b, r2.value, r2.abs_error});
        errbnd = std::max(epsabs, epsrel * std::abs(area.to_double()));
        if (errsum.to_double() > errbnd) {
            if (iroff1 >= 6 || iroff2 >= 20) out.status = QuadratureStatus::roundoff;
            if (last == limit) out.status = QuadratureStatus::limit_reached;
            // In double-double: the endpoints of a DD interval can differ below a double's resolution.
            const DD ends = std::max(abs(parent.a), abs(parent.b));
            if (!(DD(1.0, 100.0 * u_dd) * (abs(b1) + DD(1000.0 * uflow)) < ends))
                out.status = QuadratureStatus::bad_integrand;
        }
        if (out.status != QuadratureStatus::ok || errsum.to_double() <= errbnd) break;
    }
    DD result;
    for (const Interval& iv : list) result += iv.area;
    out.value = result;
    out.abs_error = errsum.to_double();
    out.intervals = static_cast<int>(list.size());
    return out;
}

}  // namespace NiTCAD::verify::reference
