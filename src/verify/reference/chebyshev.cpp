// Chebyshev collocation with domain splitting (Unit 55, TVR-1, TC7); see the header for sources and choices.
#include "NiTCAD/verify/reference/chebyshev.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace NiTCAD::verify::reference {

namespace {

// Taylor series of sin on [-pi/2, pi/2]: the terms decrease from the first, |t|^(2k+1)/(2k+1)! <= 1.57.
DD sin_series(const DD& t) {
    if (t.hi == 0.0) return DD(0.0);
    const DD t2 = t * t;
    DD term = t, sum = t;
    for (int k = 1; k < 60; ++k) {
        term = -(term * t2) / DD((2.0 * k) * (2.0 * k + 1.0));
        sum = sum + term;
        if (std::abs(term.hi) < 0x1p-112 * std::abs(sum.hi)) break;
    }
    return sum;
}

}  // namespace

DD sin_pi_ratio(long long m, long long d) {
    long long r = m % (2 * d);
    if (r < 0) r += 2 * d;
    double sign = 1.0;
    if (r >= d) {  // sin(theta + pi) = -sin(theta)
        r -= d;
        sign = -1.0;
    }
    if (2 * r > d) r = d - r;  // sin(pi - theta) = sin(theta)
    const DD s = sin_series(pi() * DD(static_cast<double>(r)) / DD(static_cast<double>(d)));
    return sign < 0.0 ? -s : s;
}

DDMatrix multiply(const DDMatrix& a, const DDMatrix& b) {
    DDMatrix c(a.rows, b.cols);
    for (std::size_t i = 0; i < a.rows; ++i)
        for (std::size_t j = 0; j < b.cols; ++j) {
            DD s(0.0);
            for (std::size_t k = 0; k < a.cols; ++k) s = s + a(i, k) * b(k, j);
            c(i, j) = s;
        }
    return c;
}

DDVector chebyshev_points(int n) {
    DDVector x(static_cast<std::size_t>(n) + 1);
    if (n == 0) {
        x[0] = DD(1.0);
        return x;
    }
    for (int j = 0; j <= n; ++j) x[static_cast<std::size_t>(j)] = sin_pi_ratio(n - 2LL * j, 2LL * n);
    return x;
}

DDMatrix chebyshev_differentiation(int n) {
    const std::size_t m = static_cast<std::size_t>(n) + 1;
    DDMatrix d(m, m);
    if (n == 0) return d;
    auto c = [n](int i) { return (i == 0 || i == n) ? 2.0 : 1.0; };
    for (int i = 0; i <= n; ++i) {
        DD row(0.0);
        for (int j = 0; j <= n; ++j) {
            if (i == j) continue;
            // x_i - x_j = cos(i pi / N) - cos(j pi / N) = -2 sin((i + j) pi / 2N) sin((i - j) pi / 2N)
            const DD dx = -ldexp(sin_pi_ratio(i + j, 2LL * n) * sin_pi_ratio(i - j, 2LL * n), 1);
            const double sgn = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            const DD e = DD(sgn * c(i) / c(j)) / dx;  // c(i) / c(j) is 1, 2 or 1/2: exact
            d(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = e;
            row = row + e;
        }
        d(static_cast<std::size_t>(i), static_cast<std::size_t>(i)) = -row;  // (6.6)
    }
    return d;
}

DDVector chebyshev_coefficients(const DDVector& values) {
    const int n = static_cast<int>(values.size()) - 1;
    DDVector a(values.size(), DD(0.0));
    if (n == 0) {
        a[0] = values[0];
        return a;
    }
    const DDVector x = chebyshev_points(n);
    auto cos_pi = [&](long long m) {  // cos(m pi / N) = x_m, reduced
        long long r = m % (2LL * n);
        if (r > n) r = 2LL * n - r;
        return x[static_cast<std::size_t>(r)];
    };
    for (int k = 0; k <= n; ++k) {
        DD s(0.0);
        for (int j = 0; j <= n; ++j) {
            DD t = values[static_cast<std::size_t>(j)] * cos_pi(static_cast<long long>(k) * j);
            if (j == 0 || j == n) t = ldexp(t, -1);
            s = s + t;
        }
        s = s * DD(2.0) / DD(static_cast<double>(n));
        if (k == 0 || k == n) s = ldexp(s, -1);
        a[static_cast<std::size_t>(k)] = s;
    }
    return a;
}

DD chebyshev_evaluate(const DDVector& a, const DD& xi) {
    DD b1(0.0), b2(0.0);
    const DD two_xi = ldexp(xi, 1);
    for (std::size_t k = a.size(); k-- > 1;) {
        const DD b0 = a[k] + two_xi * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return a.empty() ? DD(0.0) : a[0] + xi * b1 - b2;
}

bool lu_solve(DDMatrix a, DDVector& b) {
    const std::size_t n = a.rows;
    // Row equilibration by powers of two (exact): each row's largest entry in [1, 2), so partial pivoting compares
    // rows on one scale (collocation rows carry entries of order N^4, the interface and boundary rows of order 1).
    for (std::size_t i = 0; i < n; ++i) {
        double big = 0.0;
        for (std::size_t j = 0; j < n; ++j) big = std::max(big, std::abs(a(i, j).hi));
        if (big == 0.0) return false;
        int e = 0;
        (void)std::frexp(big, &e);
        for (std::size_t j = 0; j < n; ++j) a(i, j) = ldexp(a(i, j), 1 - e);
        b[i] = ldexp(b[i], 1 - e);
    }
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (abs(a(p, k)) < abs(a(i, k))) p = i;
        if (a(p, k).hi == 0.0) return false;
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(a(k, j), a(p, j));
            std::swap(b[k], b[p]);
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            const DD l = a(i, k) / a(k, k);
            if (l.hi == 0.0) continue;
            for (std::size_t j = k + 1; j < n; ++j) a(i, j) = a(i, j) - l * a(k, j);
            b[i] = b[i] - l * b[k];
        }
    }
    for (std::size_t k = n; k-- > 0;) {
        DD s = b[k];
        for (std::size_t j = k + 1; j < n; ++j) s = s - a(k, j) * b[j];
        b[k] = s / a(k, k);
    }
    return true;
}

TwoDomainSolution solve_two_domain(const TwoDomainProblem& p) {
    TwoDomainSolution s;
    s.a = p.a;
    s.c = p.c;
    s.b = p.b;
    const int n = p.n;
    if (n < 2 || !p.f || !(p.a < p.c) || !(p.c < p.b)) return s;
    const std::size_t m = static_cast<std::size_t>(n) + 1;
    const DDVector xi = chebyshev_points(n);
    const DDMatrix d = chebyshev_differentiation(n);
    const DDMatrix d2 = multiply(d, d);
    const DD sl = DD(2.0) / (p.c - p.a), sr = DD(2.0) / (p.b - p.c);
    s.x_left.resize(m);
    s.x_right.resize(m);
    for (std::size_t j = 0; j < m; ++j) {
        s.x_left[j] = ldexp(p.a + p.c, -1) + ldexp(p.c - p.a, -1) * xi[j];
        s.x_right[j] = ldexp(p.c + p.b, -1) + ldexp(p.b - p.c, -1) * xi[j];
    }
    s.x_left[0] = p.c;
    s.x_left[m - 1] = p.a;
    s.x_right[0] = p.b;
    s.x_right[m - 1] = p.c;

    DDMatrix a(2 * m, 2 * m);
    DDVector rhs(2 * m, DD(0.0));
    const DD kl = -(p.eps_left * sl * sl), kr = -(p.eps_right * sr * sr);
    for (std::size_t j = 1; j + 1 < m; ++j) {
        for (std::size_t k = 0; k < m; ++k) {
            a(j, k) = kl * d2(j, k);
            a(m + j, m + k) = kr * d2(j, k);
        }
        rhs[j] = p.f(s.x_left[j]);
        rhs[m + j] = p.f(s.x_right[j]);
    }
    // Outer Dirichlet conditions (22.3), (22.4).
    a(m - 1, m - 1) = DD(1.0);
    rhs[m - 1] = p.u_a;
    a(m, m) = DD(1.0);
    rhs[m] = p.u_b;
    // Patching (22.5): u continuous, and the flux eps u' continuous.
    a(0, 0) = DD(1.0);
    a(0, 2 * m - 1) = DD(-1.0);
    const DD fl = p.eps_left * sl, fr = p.eps_right * sr;
    for (std::size_t k = 0; k < m; ++k) {
        a(2 * m - 1, k) = fl * d(0, k);
        a(2 * m - 1, m + k) = -(fr * d(m - 1, k));
    }
    if (!lu_solve(std::move(a), rhs)) return s;
    s.u_left.assign(rhs.begin(), rhs.begin() + static_cast<std::ptrdiff_t>(m));
    s.u_right.assign(rhs.begin() + static_cast<std::ptrdiff_t>(m), rhs.end());
    DD gl(0.0), gr(0.0);
    for (std::size_t k = 0; k < m; ++k) {
        gl = gl + d(0, k) * s.u_left[k];
        gr = gr + d(m - 1, k) * s.u_right[k];
    }
    s.flux_left = fl * gl;
    s.flux_right = fr * gr;
    s.coef_left = chebyshev_coefficients(s.u_left);
    s.coef_right = chebyshev_coefficients(s.u_right);
    auto tail = [](const DDVector& c) {
        return std::abs(c[c.size() - 1].to_double()) + std::abs(c[c.size() - 2].to_double());
    };
    s.tail_estimate = std::max(tail(s.coef_left), tail(s.coef_right));
    s.ok = true;
    return s;
}

DD TwoDomainSolution::eval(const DD& x) const {
    if (!(c < x)) return chebyshev_evaluate(coef_left, (ldexp(x, 1) - a - c) / (c - a));
    return chebyshev_evaluate(coef_right, (ldexp(x, 1) - c - b) / (b - c));
}

}  // namespace NiTCAD::verify::reference
