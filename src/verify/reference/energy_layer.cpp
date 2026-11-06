#include "NiTCAD/verify/reference/energy_layer.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace NiTCAD::verify::reference {

namespace {

using CD = ComplexDD;

CD operator+(const CD& a, const CD& b) { return {a.re + b.re, a.im + b.im}; }
CD operator-(const CD& a, const CD& b) { return {a.re - b.re, a.im - b.im}; }
CD operator*(const CD& a, const CD& b) { return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}; }
CD operator*(const DD& s, const CD& a) { return {s * a.re, s * a.im}; }
CD operator/(const CD& a, const CD& b) {
    const DD d = b.re * b.re + b.im * b.im;
    return {(a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d};
}
DD magnitude(const CD& a) { return sqrt(a.re * a.re + a.im * a.im); }
CD real(const DD& v) { return {v, DD(0.0)}; }

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

// Polynomial products on coefficient vectors, lowest degree first.
std::vector<DD> times(const std::vector<DD>& a, const std::vector<DD>& b) {
    std::vector<DD> r(a.size() + b.size() - 1, DD(0.0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) r[i + j] = r[i + j] + a[i] * b[j];
    }
    return r;
}

}  // namespace

std::expected<PolynomialRoots, std::string> durand_kerner(const std::vector<DD>& c) {
    const std::size_t n = c.size();
    if (n == 0) return std::unexpected("durand_kerner: no coefficients");
    for (const DD& v : c) {
        if (!std::isfinite(v.hi)) return std::unexpected("durand_kerner: a coefficient is not finite");
    }
    if (c[0].hi == 0.0) return std::unexpected("durand_kerner: zero constant term (a root at 0)");
    // Scale z = s / sigma, sigma = |c0|^{1/n}: the roots' geometric mean is then 1.
    const DD sigma = exp(log(abs(c[0])) / DD(static_cast<double>(n)));
    std::vector<DD> d(n);
    for (std::size_t k = 0; k < n; ++k) d[k] = c[k] / pow(sigma, DD(static_cast<double>(n - k)));
    const auto p = [&](const CD& z) {
        CD v = real(DD(1.0));
        for (std::size_t k = n; k-- > 0;) v = v * z + real(d[k]);
        return v;
    };
    const auto scale = [&](const CD& z) {
        const DD m = magnitude(z);
        DD s(0.0), zk(1.0);
        for (std::size_t k = 0; k < n; ++k) {
            s = s + abs(d[k]) * zk;
            zk = zk * m;
        }
        return s + zk;
    };
    std::vector<CD> z(n);
    const CD seed{DD(0.4), DD(0.9)};
    CD w = real(DD(1.0));
    for (std::size_t i = 0; i < n; ++i) {
        z[i] = w;
        w = w * seed;
    }
    PolynomialRoots out;
    int settled = 0;
    for (int it = 1; it <= 500; ++it) {
        double worst = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            CD den = real(DD(1.0));
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i) den = den * (z[i] - z[j]);
            }
            const CD step = p(z[i]) / den;
            z[i] = z[i] - step;
            worst = std::max(worst, (magnitude(step) / std::max(DD(1.0), magnitude(z[i]))).hi);
        }
        out.iterations = it;
        if (!std::isfinite(worst)) return std::unexpected("durand_kerner: the iteration diverged");
        if (worst <= 1e-30) {
            if (++settled == 2) break;
        }
        if (it == 500) return std::unexpected("durand_kerner: no convergence in 500 iterations");
    }
    for (std::size_t i = 0; i < n; ++i) {
        out.residual = std::max(out.residual, (magnitude(p(z[i])) / scale(z[i])).hi);
        out.roots.push_back(sigma * z[i]);
    }
    return out;
}

std::expected<LayerResult, std::string> energy_layer(const LayerInput& in) {
    if (!positive_finite(in.q_C) || !positive_finite(in.k_J_per_K) || !positive_finite(in.eps_F_per_cm) ||
        !positive_finite(in.density_cm3) || !positive_finite(in.field_V_per_cm) ||
        !positive_finite(in.mobility_cm2_Vs) || !std::isfinite(in.mobility_slope) ||
        !positive_finite(in.temperature_K) || !positive_finite(in.lattice_K) || !positive_finite(in.tau_w_s) ||
        !positive_finite(in.heat_coefficient) || !(in.temperature_K >= in.lattice_K)) {
        return std::unexpected("energy_layer: an input is not finite, or out of range");
    }
    const DD q(in.q_C), k(in.k_J_per_K), eps(in.eps_F_per_cm), N(in.density_cm3), E0(in.field_V_per_cm);
    const DD mu(in.mobility_cm2_Vs), mup(in.mobility_slope), Th(in.temperature_K), TL(in.lattice_K);
    const DD tau(in.tau_w_s), c(in.heat_coefficient);
    const DD J = q * N * mu * E0;
    const DD kappa = c * k * k / q * N * mu * Th;
    const DD G = DD(1.5) * k * N / tau;
    const DD H = DD(1.5) * k * (Th - TL) / tau;  // the relaxation term's density coefficient
    // P = A D - B C, lowest degree first.
    const std::vector<DD> A{-mu * N * q * q / eps, mu * q * E0, mu * k * Th};
    const std::vector<DD> D{G, -c * k / q * J, -kappa};
    const std::vector<DD> B{mup * q * N * E0, mu * k * N};
    const std::vector<DD> C{J * q / eps, H};
    std::vector<DD> P = times(A, D);
    const std::vector<DD> BC = times(B, C);
    for (std::size_t i = 0; i < BC.size(); ++i) P[i] = P[i] - BC[i];
    const DD lead = P[4];
    if (lead.hi == 0.0) return std::unexpected("energy_layer: the quartic's leading coefficient is 0");
    std::vector<DD> monic(4);
    for (std::size_t i = 0; i < 4; ++i) monic[i] = P[i] / lead;
    const auto r = durand_kerner(monic);
    if (!r) return std::unexpected("energy_layer: " + r.error());
    LayerResult out;
    out.iterations = r->iterations;
    out.root_residual = r->residual;
    double slowest = INFINITY;
    for (std::size_t i = 0; i < 4; ++i) {
        const CD& s = r->roots[i];
        out.roots[i] = {s.re.to_double(), s.im.to_double()};
        slowest = std::min(slowest, std::abs(s.re.to_double()));
        // The second method: det of the unreduced system at s, relative to the product of its rows' norms.
        const CD s2 = s * s;
        const CD m[3][3] = {
            {DD(-1.0) * mu * N * q * s, real(mu * q * E0) + (mu * k * Th) * s, real(mup * q * N * E0) + (mu * k * N) * s},
            {eps * s2, real(-q), real(DD(0.0))},
            {J * s, real(H), real(G) - (c * k / q * J) * s - kappa * s2}};
        const CD det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                       m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                       m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        // Relative to the sum of the six expansion terms' magnitudes: the cancellation the zero takes.
        DD terms(0.0);
        for (const auto& [a, b, cc] : {std::array{0, 1, 2}, std::array{0, 2, 1}, std::array{1, 0, 2},
                                       std::array{1, 2, 0}, std::array{2, 0, 1}, std::array{2, 1, 0}}) {
            terms = terms + magnitude(m[0][a]) * magnitude(m[1][b]) * magnitude(m[2][cc]);
        }
        out.determinant_residual = std::max(out.determinant_residual, (magnitude(det) / terms).hi);
    }
    if (!(slowest > 0.0)) return std::unexpected("energy_layer: a root on the imaginary axis (no decay)");
    out.slowest_rate_per_cm = slowest;
    out.length_cm = 2.0 / slowest;
    return out;
}

}  // namespace NiTCAD::verify::reference
