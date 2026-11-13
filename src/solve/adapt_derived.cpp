// Derived adaptation goals (Unit 17c, DG1-DG3, DG12). The layer table forbids solve -> analysis,
// so the three figures are restated here from analysis/dc.cpp (threshold_max_gm,
// max_transconductance, subthreshold_swing at its default ceiling 1e-2) and curve.cpp (derivative);
// tests/solve/adapt_derived_test.cpp holds them equal to dc.hpp's, bit for bit.
#include "NiTCAD/solve/adapt_derived.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <utility>

namespace NiTCAD::solve {

namespace {

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

constexpr double swing_ceiling = 1e-2;  // subthreshold_swing's default ceiling_fraction

// curve.cpp derivative at point k, and its partials with respect to y[k-1], y[k], y[k+1] (zero
// where a neighbour is absent).
struct Slope {
    double value = 0.0;
    double d_prev = 0.0;
    double d_self = 0.0;
    double d_next = 0.0;
};

Slope slope_at(const std::vector<double>& x, const std::vector<double>& y, std::size_t k) {
    const std::size_t n = x.size();
    Slope s;
    if (k == 0) {
        const double h = x[1] - x[0];
        s.value = (y[1] - y[0]) / h;
        s.d_self = -1.0 / h;
        s.d_next = 1.0 / h;
    } else if (k == n - 1) {
        const double h = x[n - 1] - x[n - 2];
        s.value = (y[n - 1] - y[n - 2]) / h;
        s.d_prev = -1.0 / h;
        s.d_self = 1.0 / h;
    } else {
        const double hl = x[k] - x[k - 1], hr = x[k + 1] - x[k];
        const double den = hl * hr * (hl + hr);
        s.value = (hl * hl * y[k + 1] - hr * hr * y[k - 1] + (hr * hr - hl * hl) * y[k]) / den;
        s.d_prev = -hr * hr / den;
        s.d_self = (hr * hr - hl * hl) / den;
        s.d_next = hl * hl / den;
    }
    return s;
}

bool strictly_monotone(const std::vector<double>& x) {
    if (x.size() < 2) return false;
    const bool up = x[1] > x[0];
    for (std::size_t k = 1; k < x.size(); ++k) {
        if (up ? !(x[k] > x[k - 1]) : !(x[k] < x[k - 1])) return false;
    }
    return true;
}

}  // namespace

std::expected<DerivedGoal, base::Error> derived_goal(GoalKind kind, const TransferCurve& curve) {
    if (kind == GoalKind::current) return std::unexpected(invalid("a current is not a derived goal"));
    const std::size_t n = curve.gate_V.size();
    if (curve.current.size() != n || curve.resolution.size() != n) {
        return std::unexpected(invalid("a transfer curve needs one current and resolution per gate bias"));
    }
    if (n < 2) return std::unexpected(invalid("a transfer curve needs at least 2 points"));
    for (std::size_t k = 0; k < n; ++k) {
        if (!std::isfinite(curve.gate_V[k]) || !std::isfinite(curve.current[k])) {
            return std::unexpected(invalid("a transfer curve's value is not finite"));
        }
    }
    // DG1: a p-channel curve (its largest |I| negative) read with currents and biases times -1.
    std::size_t largest = 0;
    for (std::size_t k = 1; k < n; ++k) {
        if (std::abs(curve.current[k]) > std::abs(curve.current[largest])) largest = k;
    }
    const double s = curve.current[largest] < 0.0 ? -1.0 : 1.0;
    std::vector<double> x(n), y(n);
    for (std::size_t k = 0; k < n; ++k) {
        x[k] = s * curve.gate_V[k];
        y[k] = s * curve.current[k];
    }
    if (!strictly_monotone(x)) return std::unexpected(invalid("x is not strictly monotone"));
    DerivedGoal g;
    if (kind == GoalKind::vt || kind == GoalKind::gm) {
        // peak_transconductance: the largest derivative, the first on a tie.
        std::size_t k = 0;
        Slope best = slope_at(x, y, 0);
        for (std::size_t j = 1; j < n; ++j) {
            const Slope sj = slope_at(x, y, j);
            if (sj.value > best.value) {
                best = sj;
                k = j;
            }
        }
        if (!(best.value > 0.0)) return std::unexpected(invalid("the transconductance is nowhere positive"));
        g.first = k == 0 ? 0 : k - 1;
        g.last = std::min(k + 1, n - 1);
        g.point = k;
        // dG/dy at first..last: the slope's partials laid out on the window.
        std::vector<double> dg(g.last - g.first + 1, 0.0);
        const auto put = [&](std::size_t at, double v) {
            if (at >= g.first && at <= g.last) dg[at - g.first] += v;
        };
        if (k > 0) put(k - 1, best.d_prev);
        put(k, best.d_self);
        if (k + 1 < n) put(k + 1, best.d_next);
        if (kind == GoalKind::gm) {
            g.value = best.value;
            // dG/dI = s dG/dy (y = s I, x = s V: g_m = dy/dx is dI/dV).
            for (double& d : dg) d *= s;
        } else {
            // threshold_max_gm: x_k - y_k / g_m - drain / 2, in the read sign; V_T = s times it.
            const double gm = best.value;
            const double drain = s * curve.drain_V;
            g.value = s * (x[k] - y[k] / gm - 0.5 * drain);
            // d(x_k - y_k/g)/dy_j = -delta_jk / g + y_k / g^2 dg/dy_j; dV_T/dI_j = s * s * that.
            for (double& d : dg) d *= y[k] / (gm * gm);
            dg[k - g.first] -= 1.0 / gm;
        }
        g.sensitivity = std::move(dg);
    } else {
        // subthreshold_swing: the steepest segment between usable points, ceiling 1e-2 of the largest.
        const double ceiling = swing_ceiling * *std::max_element(y.begin(), y.end());
        const auto usable = [&](std::size_t k) {
            const double r = curve.resolution[k];
            return y[k] > 0.0 && y[k] > r && y[k] <= ceiling;
        };
        std::optional<std::size_t> best;
        double swing = 0.0;
        for (std::size_t k = 0; k + 1 < n; ++k) {
            if (!usable(k) || !usable(k + 1)) continue;
            const double sw = (x[k + 1] - x[k]) / (std::log10(y[k + 1]) - std::log10(y[k]));
            if (!(sw > 0.0) || !std::isfinite(sw)) continue;
            if (!best || sw < swing) {
                best = k;
                swing = sw;
            }
        }
        if (!best) return std::unexpected(invalid("no subthreshold segment with a rising current"));
        const std::size_t k = *best;
        g.value = 1000.0 * swing;
        g.first = k;
        g.last = k + 1;
        g.point = k;
        // S = 1000 dx / L, L = log10 y1 - log10 y0: dS/dy1 = -S / (L y1 ln 10), dS/dy0 = S / (L y0 ln 10).
        const double L = std::log10(y[k + 1]) - std::log10(y[k]);
        const double ln10 = std::log(10.0);
        g.sensitivity = {s * g.value / (L * y[k] * ln10), -s * g.value / (L * y[k + 1] * ln10)};
    }
    // DG5: the resolution through the sensitivities.
    for (std::size_t j = g.first; j <= g.last; ++j) {
        g.resolution += std::abs(g.sensitivity[j - g.first]) * curve.resolution[j];
    }
    return g;
}

}  // namespace NiTCAD::solve
