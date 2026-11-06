#include "NiTCAD/verify/reference/minority_diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "NiTCAD/verify/reference/chebyshev.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

void check(const MinorityProblem& p) {
    if (!positive(p.diffusivity_cm2_s) || !positive(p.lifetime_s) || !positive(p.length_cm) ||
        !positive(p.width_cm) || !std::isfinite(p.peak_generation) || !std::isfinite(p.centre_cm) ||
        p.centre_cm < 0.0 || p.centre_cm > p.length_cm) {
        throw std::invalid_argument("minority problem: D, tau, L, s finite and positive; G0 finite; c in [0, L]");
    }
}

void check_point(const MinorityProblem& p, double x) {
    if (!std::isfinite(x) || x < 0.0 || x > p.length_cm) throw std::invalid_argument("minority: x outside [0, L]");
}

DD k_of(const MinorityProblem& p) { return DD(1.0) / sqrt(DD(p.diffusivity_cm2_s) * DD(p.lifetime_s)); }

// 1 - e^(-y), y >= 0.
DD one_minus_exp(const DD& y) { return -expm1(-y); }

DD generation(const MinorityProblem& p, const DD& s) {
    const DD d = (s - DD(p.centre_cm)) / DD(p.width_cm);
    return DD(p.peak_generation) * exp(-ldexp(d * d, -1));
}

DD green(const MinorityProblem& p, const DD& k, const DD& x, const DD& s) {
    const DD lo = x < s ? x : s;
    const DD hi = x < s ? s : x;
    const DD L(p.length_cm);
    const DD num = exp(-k * (hi - lo)) * one_minus_exp(ldexp(k * lo, 1)) * one_minus_exp(ldexp(k * (L - hi), 1));
    return num / (ldexp(DD(p.diffusivity_cm2_s) * k, 1) * one_minus_exp(ldexp(k * L, 1)));
}

constexpr double quad_epsrel = 0.5e-28;  // the smallest relative tolerance qag15 accepts

}  // namespace

DD srh_linear_inverse_lifetime(double n0, double p0, double n_ie, double tau_n, double tau_p) {
    if (!positive(n0) || !positive(p0) || !positive(n_ie) || !positive(tau_n) || !positive(tau_p)) {
        throw std::invalid_argument("srh_linear_inverse_lifetime: inputs finite and positive");
    }
    const DD D0 = DD(tau_p) * (DD(n0) + DD(n_ie)) + DD(tau_n) * (DD(p0) + DD(n_ie));
    return (DD(n0) + DD(p0)) / D0;
}

DD minority_green(const MinorityProblem& p, double x, double s) {
    check(p);
    check_point(p, x);
    check_point(p, s);
    return green(p, k_of(p), DD(x), DD(s));
}

MinorityValue minority_excess(const MinorityProblem& p, double x) {
    check(p);
    check_point(p, x);
    const DD k = k_of(p), X(x);
    const auto f = [&](const DD& s) { return green(p, k, X, s) * generation(p, s); };
    MinorityValue out;
    out.ok = true;
    for (const auto& [a, b] : {std::pair{DD(0.0), X}, std::pair{X, DD(p.length_cm)}}) {
        if (!(a < b)) continue;
        const AdaptiveResult r = qag15(f, a, b, 0.0, quad_epsrel);
        out.value += r.value;
        out.abs_error += r.abs_error;
        out.ok = out.ok && r.status == QuadratureStatus::ok;
    }
    return out;
}

MinorityCollocation minority_collocation(const MinorityProblem& p, const std::vector<double>& xs, int pieces,
                                         int degree) {
    check(p);
    if (pieces < 1 || degree < 4) throw std::invalid_argument("minority_collocation: pieces >= 1, degree >= 4");
    for (const double x : xs) check_point(p, x);
    const auto P = static_cast<std::size_t>(pieces);
    const auto n = static_cast<std::size_t>(degree);
    const std::size_t w = n + 1, size = P * w;
    const DD L(p.length_cm), h = L / DD(static_cast<double>(pieces));
    const DD scale = DD(2.0) / h;  // d/dx = (2 / h) d/dxi
    const DDVector xi = chebyshev_points(degree);
    const DDMatrix D1 = chebyshev_differentiation(degree);
    const DDMatrix D2 = multiply(D1, D1);
    const auto x_of = [&](std::size_t q, std::size_t j) {
        return h * DD(static_cast<double>(q)) + ldexp((xi[j] + DD(1.0)) * h, -1);
    };
    DDMatrix A(size, size);
    DDVector b(size, DD(0.0));
    std::size_t row = 0;
    const DD Dc(p.diffusivity_cm2_s), inv_tau = DD(1.0) / DD(p.lifetime_s);
    for (std::size_t q = 0; q < P; ++q) {
        for (std::size_t j = 1; j < n; ++j, ++row) {  // interior collocation
            for (std::size_t m = 0; m < w; ++m) A(row, q * w + m) = Dc * scale * scale * D2(j, m);
            A(row, q * w + j) -= inv_tau;
            b[row] = -generation(p, x_of(q, j));
        }
    }
    A(row, 0 * w + n) = DD(1.0);  // u(0) = 0: piece 0's xi = -1 point
    ++row;
    A(row, (P - 1) * w + 0) = DD(1.0);  // u(L) = 0: the last piece's xi = 1 point
    ++row;
    for (std::size_t q = 0; q + 1 < P; ++q) {  // piece q's xi = 1 meets piece q + 1's xi = -1
        A(row, q * w + 0) = DD(1.0);
        A(row, (q + 1) * w + n) = DD(-1.0);
        ++row;
        for (std::size_t m = 0; m < w; ++m) {
            A(row, q * w + m) = D1(0, m);
            A(row, (q + 1) * w + m) = -D1(n, m);
        }
        ++row;
    }
    MinorityCollocation out;
    if (!lu_solve(std::move(A), b)) return out;
    std::vector<DDVector> coef(P);
    for (std::size_t q = 0; q < P; ++q) {
        coef[q] = chebyshev_coefficients(DDVector(b.begin() + static_cast<std::ptrdiff_t>(q * w),
                                                  b.begin() + static_cast<std::ptrdiff_t>((q + 1) * w)));
        out.tail_estimate = std::max({out.tail_estimate, std::abs(coef[q][n].hi), std::abs(coef[q][n - 1].hi)});
    }
    for (const double x : xs) {
        const auto q = std::min(static_cast<std::size_t>(x / h.hi), P - 1);
        const DD t = ldexp((DD(x) - h * DD(static_cast<double>(q))) / h, 1) - DD(1.0);
        out.values.push_back(chebyshev_evaluate(coef[q], t));
    }
    out.ok = true;
    return out;
}

}  // namespace NiTCAD::verify::reference
