// The Gauss-Kronrod 7/15 nodes and weights, computed in double-double (Unit 55, TVR-1, TC3). The method is in
// the header.
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace NiTCAD::verify::reference {
namespace {

// P_0(x) .. P_n(x) by (k + 1) P_{k+1} = (2k + 1) x P_k - k P_{k-1}.
std::vector<DD> legendre_all(int n, const DD& x) {
    std::vector<DD> p(static_cast<std::size_t>(n) + 1);
    p[0] = DD(1.0);
    if (n >= 1) p[1] = x;
    for (int k = 1; k < n; ++k)
        p[static_cast<std::size_t>(k) + 1] =
            (DD(2.0 * k + 1.0) * x * p[static_cast<std::size_t>(k)] - DD(static_cast<double>(k)) * p[static_cast<std::size_t>(k) - 1]) /
            DD(k + 1.0);
    return p;
}

// P_n'(x) = n (x P_n - P_{n-1}) / (x^2 - 1), |x| < 1; x^2 - 1 as -(1 - x)(1 + x), which does not cancel near |x| = 1.
DD legendre_derivative(int n, const DD& x, const std::vector<DD>& p) {
    return DD(static_cast<double>(n)) * (x * p[static_cast<std::size_t>(n)] - p[static_cast<std::size_t>(n) - 1]) /
           -((DD(1.0) - x) * (DD(1.0) + x));
}

// Gauss-Legendre: the nonnegative nodes of P_n (descending) and their weights, Newton from the classical
// estimate cos(pi (i - 1/4) / (n + 1/2)).
struct Gauss {
    std::vector<DD> x, w;
};
Gauss gauss_legendre(int n) {
    Gauss g;
    for (int i = 1; i <= (n + 1) / 2; ++i) {
        DD x = std::cos(pi().hi * (i - 0.25) / (n + 0.5));
        if (2 * i - 1 == n) x = DD(0.0);  // the middle node of an odd rule
        for (int it = 0; it < 100 && !(2 * i - 1 == n); ++it) {
            const auto p = legendre_all(n, x);
            const DD step = p[static_cast<std::size_t>(n)] / legendre_derivative(n, x, p);
            x -= step;
            if (std::abs(step.hi) < 1e-34) break;
        }
        const auto p = legendre_all(n, x);
        const DD dp = legendre_derivative(n, x, p);
        g.x.push_back(x);
        g.w.push_back(DD(2.0) / ((DD(1.0) - x) * (DD(1.0) + x) * dp * dp));
    }
    return g;
}

// Dense solve A c = r by Gaussian elimination with partial pivoting (small systems only).
std::vector<DD> solve(std::vector<std::vector<DD>> a, std::vector<DD> r) {
    const std::size_t n = r.size();
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t piv = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (abs(a[piv][k]) < abs(a[i][k])) piv = i;
        std::swap(a[k], a[piv]);
        std::swap(r[k], r[piv]);
        for (std::size_t i = k + 1; i < n; ++i) {
            const DD m = a[i][k] / a[k][k];
            for (std::size_t j = k; j < n; ++j) a[i][j] -= m * a[k][j];
            r[i] -= m * r[k];
        }
    }
    std::vector<DD> c(n);
    for (std::size_t k = n; k-- > 0;) {
        DD s = r[k];
        for (std::size_t j = k + 1; j < n; ++j) s -= a[k][j] * c[j];
        c[k] = s / a[k][k];
    }
    return c;
}

KronrodRule15 compute_rule() {
    constexpr int n = 7;
    const Gauss g7 = gauss_legendre(n);  // x: 3 positive descending, then 0
    // E_8 = P_8 + sum_{j = 0, 2, 4, 6} c_j P_j with int P_7 E_8 P_k = 0, k = 1, 3, 5, 7. The integrands have
    // degree <= 22, so a 16-point Gauss-Legendre rule (exact to degree 31) evaluates them exactly up to rounding.
    const Gauss g16 = gauss_legendre(16);
    auto triple = [&](int j, int k) {  // int_{-1}^{1} P_7 P_j P_k, by symmetry 2 * the positive half
        DD s;
        for (std::size_t i = 0; i < g16.x.size(); ++i) {
            const auto p = legendre_all(16, g16.x[i]);
            s += g16.w[i] * p[7] * p[static_cast<std::size_t>(j)] * p[static_cast<std::size_t>(k)];
        }
        return ldexp(s, 1);
    };
    std::vector<std::vector<DD>> a(4, std::vector<DD>(4));
    std::vector<DD> rhs(4);
    for (int row = 0; row < 4; ++row) {
        const int k = 2 * row + 1;
        for (int col = 0; col < 4; ++col) a[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] = triple(2 * col, k);
        rhs[static_cast<std::size_t>(row)] = -triple(8, k);
    }
    const std::vector<DD> c = solve(a, rhs);
    auto e8 = [&](const DD& x) {
        const auto p = legendre_all(8, x);
        DD s = p[8];
        for (int col = 0; col < 4; ++col) s += c[static_cast<std::size_t>(col)] * p[2 * static_cast<std::size_t>(col)];
        return s;
    };
    auto e8_derivative = [&](const DD& x) {
        const auto p = legendre_all(8, x);
        DD s = legendre_derivative(8, x, p);
        for (int col = 1; col < 4; ++col) {
            const int j = 2 * col;
            s += c[static_cast<std::size_t>(col)] * legendre_derivative(j, x, p);
        }
        return s;
    };
    // The 4 positive roots interlace with the Gauss nodes: in (g1, 1), (g2, g1), (g3, g2), (0, g3).
    std::vector<DD> bounds{DD(1.0), g7.x[0], g7.x[1], g7.x[2], DD(0.0)};
    std::vector<DD> kronrod;
    for (std::size_t i = 0; i + 1 < bounds.size(); ++i) {
        DD hi = bounds[i], lo = bounds[i + 1];
        const bool lo_negative = e8(lo).hi < 0.0;
        for (int it = 0; it < 60; ++it) {  // bisection to a double's resolution
            const DD mid = ldexp(lo + hi, -1);
            ((e8(mid).hi < 0.0) == lo_negative ? lo : hi) = mid;
        }
        DD x = ldexp(lo + hi, -1);
        for (int it = 0; it < 6; ++it) x -= e8(x) / e8_derivative(x);
        kronrod.push_back(x);
    }
    KronrodRule15 rule;
    for (int i = 0; i < 4; ++i) rule.xgk[2 * static_cast<std::size_t>(i)] = kronrod[static_cast<std::size_t>(i)];
    for (int i = 0; i < 3; ++i) rule.xgk[2 * static_cast<std::size_t>(i) + 1] = g7.x[static_cast<std::size_t>(i)];
    rule.xgk[7] = DD(0.0);
    for (int i = 0; i < 4; ++i) rule.wg[static_cast<std::size_t>(i)] = g7.w[static_cast<std::size_t>(i)];
    // Kronrod weights: the rule integrates P_j exactly for even j = 0 .. 14: w_7 P_j(0) + sum_{i < 7} 2 w_i P_j(x_i)
    // = 2 delta_{j0}.
    std::vector<std::vector<DD>> m(8, std::vector<DD>(8));
    std::vector<DD> moments(8);
    for (std::size_t i = 0; i < 8; ++i) {
        const auto p = legendre_all(14, rule.xgk[i]);
        for (std::size_t row = 0; row < 8; ++row) m[row][i] = (i == 7 ? DD(1.0) : DD(2.0)) * p[2 * row];
    }
    moments[0] = DD(2.0);
    const std::vector<DD> w = solve(m, moments);
    for (std::size_t i = 0; i < 8; ++i) rule.wgk[i] = w[i];
    return rule;
}

}  // namespace

const KronrodRule15& kronrod15() {
    static const KronrodRule15 rule = compute_rule();
    return rule;
}

}  // namespace NiTCAD::verify::reference
