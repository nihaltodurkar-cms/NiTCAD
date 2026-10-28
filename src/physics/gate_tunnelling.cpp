#include "NiTCAD/physics/gate_tunnelling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

namespace {

// ln(1 + e^x) without overflow, and its derivative 1 / (1 + e^-x), from one exponential e^-|x|.
struct Softplus {
    double value;
    double derivative;
};
Softplus softplus(double x) noexcept {
    const double e = std::exp(-std::abs(x));
    return {std::max(x, 0.0) + std::log1p(e), x >= 0.0 ? 1.0 / (1.0 + e) : e / (1.0 + e)};
}

// 4-point Gauss-Legendre on [-1, 1].
constexpr std::array<double, 4> gl_x{-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                                     0.8611363115940526};
constexpr std::array<double, 4> gl_w{0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                                     0.3478548451374538};

// The supply decays as e^(-(E - F)/kT): beyond 40 kT above the highest level it is below 4e-18 of
// its peak.
constexpr double supply_cut_kT = 40.0;

}  // namespace

double wkb_action_scale() noexcept {
    return std::sqrt(2.0 * base::m0_kg * base::q_C) / base::hbar_J_s * 1e-2;
}

namespace {

// segment_action with the ends' roots sqrt(a+), sqrt(b+) given (adjacent segments share an end, so a
// barrier's roots are taken once per point), one division per segment.
SegmentAction action_of_roots(double a, double b, double sa, double sb, double L) noexcept {
    if (a <= 0.0 && b <= 0.0) return {0.0, 0.0, 0.0};
    if (a > 0.0 && b > 0.0) {
        // Both ends inside: dividing a^(3/2) - b^(3/2) by a - b = (sqrt a - sqrt b)(sqrt a + sqrt b)
        // leaves forms without cancellation, exact at a == b:
        //     g = (2/3) (a + sqrt(a b) + b) / (sqrt a + sqrt b)
        //     dg/da = (sqrt a + 2 sqrt b) / (3 (sqrt a + sqrt b)^2), dg/db symmetrically.
        const double s = sa + sb;
        return {L * (2.0 / 3.0) * (a + sa * sb + b) / s, L * (sa + 2.0 * sb) / (3.0 * s * s),
                L * (2.0 * sa + sb) / (3.0 * s * s)};
    }
    // One end outside: the barrier ends at the turning point, g = (2/3) c^(3/2) / (c - o) with c > 0
    // the inside end's value and o <= 0 the outside one; c - o >= c, so nothing cancels.
    //     dg/dc = sqrt(c) (c - 3 o) / (3 (c - o)^2),   dg/do = g / (c - o)
    const bool a_inside = a > 0.0;
    const double c = a_inside ? a : b;
    const double o = a_inside ? b : a;
    const double sc = a_inside ? sa : sb;
    const double w = c - o;
    const double g = (2.0 / 3.0) * c * sc / w;
    const double value = L * g;
    const double d_inside = L * sc * (c - 3.0 * o) / (3.0 * w * w);
    const double d_outside = L * g / w;
    return a_inside ? SegmentAction{value, d_inside, d_outside}
                    : SegmentAction{value, d_outside, d_inside};
}

double root(double a) noexcept { return a > 0.0 ? std::sqrt(a) : 0.0; }

}  // namespace

SegmentAction segment_action(double a, double b, double L) noexcept {
    return action_of_roots(a, b, root(a), root(b), L);
}

double wkb_exponent(std::span<const BarrierSegment> barrier, double E, double* d_E,
                    std::span<double> d_W, double weight) {
    NITCAD_EXPECTS(d_W.empty() || d_W.size() == 2 * barrier.size());
    const double scale = wkb_action_scale();
    double S = 0.0;
    double dS_dE = 0.0;
    for (std::size_t k = 0; k < barrier.size(); ++k) {
        const BarrierSegment& s = barrier[k];
        const SegmentAction g = segment_action(s.W0_eV - E, s.W1_eV - E, s.length_cm);
        const double c = scale * std::sqrt(s.mass_m0);
        S += c * g.value;
        dS_dE -= c * (g.d_a + g.d_b);
        if (!d_W.empty()) {
            d_W[2 * k] += weight * c * g.d_a;
            d_W[2 * k + 1] += weight * c * g.d_b;
        }
    }
    if (d_E != nullptr) *d_E = dS_dE;
    return S;
}

double wkb_transmission(std::span<const BarrierSegment> barrier, double E) {
    return std::exp(-2.0 * wkb_exponent(barrier, E));
}

Supply tsu_esaki_supply(double E, double F1, double F2, double kT) noexcept {
    const double x1 = (F1 - E) / kT;
    const double x2 = (F2 - E) / kT;
    const Softplus s1 = softplus(x1);
    const Softplus s2 = softplus(x2);
    return {kT * (s1.value - s2.value), s1.derivative, -s2.derivative, s2.derivative - s1.derivative};
}

TunnelCurrent tsu_esaki_current(std::span<const BarrierSegment> barrier, TunnelRange range,
                                double F1, double F2, double temperature_K, double richardson,
                                std::span<double> d_W) {
    NITCAD_EXPECTS(!barrier.empty());
    NITCAD_EXPECTS(d_W.empty() || d_W.size() == 2 * barrier.size());
    NITCAD_EXPECTS(std::isfinite(range.E_low_eV) && temperature_K > 0.0);
    const double kT = base::k_B_eV_per_K * temperature_K;

    // The upper end: the given one, the barrier top or the supply cut, whichever is lowest.
    std::size_t top_at = 0;
    for (std::size_t k = 1; k < 2 * barrier.size(); ++k) {
        const double w = (k % 2 == 0) ? barrier[k / 2].W0_eV : barrier[k / 2].W1_eV;
        const double w_top = (top_at % 2 == 0) ? barrier[top_at / 2].W0_eV : barrier[top_at / 2].W1_eV;
        if (w > w_top) top_at = k;
    }
    const double top = (top_at % 2 == 0) ? barrier[top_at / 2].W0_eV : barrier[top_at / 2].W1_eV;
    enum class Cut { given, top, F1, F2, low };
    const double level = std::max({F1, F2, range.E_low_eV});
    const Cut level_is = level == F1 ? Cut::F1 : level == F2 ? Cut::F2 : Cut::low;
    double Eh = range.E_high_eV;
    Cut high_is = Cut::given;
    if (top < Eh) { Eh = top; high_is = Cut::top; }
    if (level + supply_cut_kT * kT < Eh) { Eh = level + supply_cut_kT * kT; high_is = level_is; }

    TunnelCurrent out{0.0, 0.0, 0.0, 0.0, 0.0};
    const double El = range.E_low_eV;
    if (!(Eh > El)) return out;

    const double span = Eh - El;
    const auto n = static_cast<std::size_t>(std::max(1.0, std::ceil(span / (0.5 * kT))));
    const double h = span / static_cast<double>(n);
    const double prefactor = richardson / (base::k_B_eV_per_K * base::k_B_eV_per_K);

    // wkb_exponent's sum inlined, its per-segment factor scale sqrt(m) taken once.
    const std::size_t m = barrier.size();
    const double scale = wkb_action_scale();
    // Scratch on the stack for the usual few segments (the flux runs on several threads, GR11).
    constexpr std::size_t inline_segments = 32;
    std::array<double, 3 * inline_segments> stack;
    std::vector<double> heap;
    if (m > inline_segments) heap.resize(3 * m);
    const std::span<double> scratch = m > inline_segments ? std::span<double>(heap)
                                                          : std::span<double>(stack).first(3 * m);
    const std::span<double> c = scratch.first(m);
    const std::span<double> dS_dW = scratch.subspan(m, d_W.empty() ? 0 : 2 * m);
    for (std::size_t k = 0; k < m; ++k) c[k] = scale * std::sqrt(barrier[k].mass_m0);
    double J = 0.0, dF1 = 0.0, dF2 = 0.0, dEl = 0.0, dEh = 0.0;
    for (std::size_t p = 0; p < n; ++p) {
        for (std::size_t j = 0; j < gl_x.size(); ++j) {
            const double t = (static_cast<double>(p) + 0.5 * (1.0 + gl_x[j])) / static_cast<double>(n);
            const double E = El + span * t;
            const double w = prefactor * 0.5 * h * gl_w[j];
            double S = 0.0, dS_dE = 0.0;
            double r1 = 0.0;
            for (std::size_t k = 0; k < m; ++k) {
                const BarrierSegment& s = barrier[k];
                const double a = s.W0_eV - E, b = s.W1_eV - E;
                const double r0 = k > 0 && s.W0_eV == barrier[k - 1].W1_eV ? r1 : root(a);
                r1 = root(b);
                const SegmentAction g = action_of_roots(a, b, r0, r1, s.length_cm);
                S += c[k] * g.value;
                dS_dE -= c[k] * (g.d_a + g.d_b);
                if (!dS_dW.empty()) {
                    dS_dW[2 * k] = c[k] * g.d_a;
                    dS_dW[2 * k + 1] = c[k] * g.d_b;
                }
            }
            const double TC = std::exp(-2.0 * S);
            const Supply N = tsu_esaki_supply(E, F1, F2, kT);
            const double f = TC * N.value;
            const double df_dE = TC * (N.d_E - 2.0 * dS_dE * N.value);
            J += w * f;
            dF1 += w * TC * N.d_F1;
            dF2 += w * TC * N.d_F2;
            dEl += w * df_dE * (1.0 - t);
            dEh += w * df_dE * t;
            for (std::size_t k = 0; k < dS_dW.size(); ++k) d_W[k] += -2.0 * w * f * dS_dW[k];
        }
    }
    // h = span / n: the sum scales with the span at fixed nodes' fractions.
    dEh += J / span;
    dEl -= J / span;

    switch (high_is) {
        case Cut::given: out.d_E_high = dEh; break;
        case Cut::top:
            if (!d_W.empty()) d_W[top_at] += dEh;
            break;
        case Cut::F1: dF1 += dEh; break;
        case Cut::F2: dF2 += dEh; break;
        case Cut::low: dEl += dEh; break;
    }
    out.J = J;
    out.d_F1 = dF1;
    out.d_F2 = dF2;
    out.d_E_low = dEl;
    return out;
}

}  // namespace NiTCAD::physics
