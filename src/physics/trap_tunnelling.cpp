#include "NiTCAD/physics/trap_tunnelling.hpp"

#include <algorithm>
#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

namespace {

// 8-point Gauss-Legendre on [-1, 1]: abscissae (positive half) and weights.
constexpr double gl_x[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                            0.9602898564975363};
constexpr double gl_w[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745,
                            0.1012285362903763};
constexpr int min_panels = 16;
constexpr double tail = 50.0;  // the exponent's fall past its maximum where the rule stops

}  // namespace

TrapTunnellingFactor hurkx_enhancement(double a, double mass_ratio, double T_K,
                                       double field_V_per_cm) {
    NITCAD_EXPECTS(std::isfinite(a) && a > 0.0 && a <= 700.0);
    NITCAD_EXPECTS(std::isfinite(mass_ratio) && mass_ratio > 0.0);
    NITCAD_EXPECTS(std::isfinite(T_K) && T_K > 0.0);
    NITCAD_EXPECTS(std::isfinite(field_V_per_cm) && field_V_per_cm >= 0.0);
    if (field_V_per_cm == 0.0) return {0.0, 0.0};
    const double dE = a * base::k_B_J_per_K * T_K;
    const double K = 4.0 / 3.0 * std::sqrt(2.0 * mass_ratio * base::m0_kg * dE * dE * dE) /
                     (base::q_C * base::hbar_J_s * field_V_per_cm * 100.0);
    if (!std::isfinite(K)) return {0.0, 0.0};
    // phi(v) = a v^2 - K v^3 on [0, 1]: maximum at v* = 2a / (3K) (or at 1).
    const double vs = 2.0 * a / (3.0 * K);
    const auto phi = [&](double v) { return v * v * (a - K * v); };
    const double top = vs < 1.0 ? phi(vs) : phi(1.0);
    double v_max = 1.0;
    if (vs < 1.0 && phi(1.0) < top - tail) {
        // phi is decreasing on [v*, 1]: bisect phi(v) = top - tail.
        double lo = vs, hi = 1.0;
        for (int it = 0; it < 200 && hi - lo > 1e-16 * hi; ++it) {
            const double mid = 0.5 * (lo + hi);
            (phi(mid) > top - tail ? lo : hi) = mid;
        }
        v_max = hi;
    }
    // Panels: at least 16, and enough that phi changes by at most 4 across one (|phi'| <= 2a v +
    // 3K v^2 <= 2a v_max + 3K v_max^2).
    const double slope = 2.0 * a * v_max + 3.0 * K * v_max * v_max;
    const int n = std::max(min_panels, static_cast<int>(std::ceil(slope * v_max / 4.0)));
    const double h = v_max / n;
    double I1 = 0.0, I4 = 0.0;  // int v e^phi, int v^4 e^phi
    for (int p = 0; p < n; ++p) {
        const double c = (p + 0.5) * h;
        for (int k = 0; k < 4; ++k) {
            for (const double s : {-1.0, 1.0}) {
                const double v = c + s * 0.5 * h * gl_x[k];
                const double e = gl_w[k] * std::exp(phi(v));
                const double v2 = v * v;
                I1 += v * e;
                I4 += v2 * v2 * e;
            }
        }
    }
    I1 *= 0.5 * h;
    I4 *= 0.5 * h;
    return {2.0 * a * I1, 2.0 * a * (K * I4) / field_V_per_cm};
}

}  // namespace NiTCAD::physics
