// The tanh-sinh abscissa-weight table, computed in double-double (Unit 55, TVR-1, TC4). The method is in the
// header.
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

namespace NiTCAD::verify::reference {
namespace {

TanhSinhTable compute_table(int max_level) {
    TanhSinhTable tab;
    tab.max_level = max_level;
    const double eps2 = u_dd * u_dd * u_dd;
    const DD half_pi = ldexp(pi(), -1);
    const DD one(1.0);
    for (int j = 0;; ++j) {
        const DD t = ldexp(DD(static_cast<double>(j)), -max_level);  // exact
        const DD et = exp(t);
        const DD cosh_t = ldexp(et + one / et, -1);
        const DD sinh_t = ldexp(expm1(t) - expm1(-t), -1);  // no cancellation: the two have opposite signs
        const DD u1 = half_pi * cosh_t;
        const DD u2 = half_pi * sinh_t;
        const DD eu = exp(u2);
        const DD cosh_u = ldexp(eu + one / eu, -1);
        const DD w = u1 / (cosh_u * cosh_u);
        if (w.hi < eps2) break;
        tab.y.push_back(one / (eu * cosh_u));
        tab.w.push_back(w);
    }
    return tab;
}

}  // namespace

const TanhSinhTable& tanh_sinh_table() {
    static const TanhSinhTable tab = compute_table(10);
    return tab;
}

}  // namespace NiTCAD::verify::reference
