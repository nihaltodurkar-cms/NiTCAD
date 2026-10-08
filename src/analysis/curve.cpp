#include "NiTCAD/analysis/curve.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <utility>

#include "common.hpp"

namespace NiTCAD::analysis {

namespace {

using detail::below;
using detail::invalid;

// Reads one value per bias point of a sweep: the swept bias, value(point) and, when every point
// has one, resolution(point).
template <class Value, class Resolution>
std::expected<Curve, base::Error> read_sweep(const results::Sweep& sweep, std::size_t swept,
                                             Value value, Resolution resolution) {
    std::vector<double> x, y, r;
    bool resolved = true;
    for (std::size_t k = 0; k < sweep.points.size(); ++k) {
        const results::BiasPoint& p = sweep.points[k];
        const auto v = value(p);
        if (swept >= p.bias_V.size() || !v) return std::unexpected(invalid("contact out of range", k));
        x.push_back(p.bias_V[swept]);
        y.push_back(*v);
        const auto b = resolution(p);
        resolved = resolved && b.has_value();
        if (b) r.push_back(*b);
    }
    if (!resolved) r.clear();
    return make_curve(std::move(x), std::move(y), std::move(r));
}

}  // namespace

std::expected<Curve, base::Error> make_curve(std::vector<double> x, std::vector<double> y,
                                             std::vector<double> resolution) {
    if (x.size() != y.size()) return std::unexpected(invalid("x and y differ in size"));
    if (x.size() < 2) return std::unexpected(invalid("a curve needs at least 2 points"));
    if (!resolution.empty() && resolution.size() != x.size()) {
        return std::unexpected(invalid("resolution is not one value per point"));
    }
    for (std::size_t k = 0; k < x.size(); ++k) {
        if (!std::isfinite(x[k]) || !std::isfinite(y[k])) {
            return std::unexpected(invalid("a point is not finite", k));
        }
        if (!resolution.empty() && !(std::isfinite(resolution[k]) && resolution[k] >= 0.0)) {
            return std::unexpected(invalid("a resolution is not finite and non-negative", k));
        }
    }
    return Curve{std::move(x), std::move(y), std::move(resolution)};
}

std::expected<Curve, base::Error> current_curve(const results::Sweep& sweep, std::size_t swept,
                                                std::size_t sense) {
    return read_sweep(
        sweep, swept,
        [&](const results::BiasPoint& p) -> std::optional<double> {
            if (sense >= p.terminal_current.size()) return std::nullopt;
            return p.terminal_current[sense];
        },
        [&](const results::BiasPoint& p) -> std::optional<double> {
            if (sense >= p.terminal_current_resolution.size()) return std::nullopt;
            return p.terminal_current_resolution[sense];
        });
}

std::expected<Curve, base::Error> gate_charge_curve(const results::Sweep& sweep,
                                                    std::size_t swept, std::size_t gate) {
    return read_sweep(
        sweep, swept,
        [&](const results::BiasPoint& p) -> std::optional<double> {
            if (gate >= p.gate_charge.size()) return std::nullopt;
            return p.gate_charge[gate];
        },
        [](const results::BiasPoint&) -> std::optional<double> { return std::nullopt; });
}

std::expected<Curve, base::Error> capacitance_curve(const results::SmallSignal& small_signal,
                                                    std::size_t frequency, std::size_t row,
                                                    std::size_t column, std::size_t swept) {
    if (frequency >= small_signal.frequency_Hz.size()) {
        return std::unexpected(invalid("frequency out of range"));
    }
    if (!(small_signal.frequency_Hz[frequency] > 0.0)) {
        return std::unexpected(invalid("a capacitance needs a positive frequency"));
    }
    if (row >= small_signal.contacts || column >= small_signal.contacts) {
        return std::unexpected(invalid("contact out of range"));
    }
    std::vector<double> x, y;
    for (std::size_t k = 0; k < small_signal.points.size(); ++k) {
        const results::SmallSignalPoint& p = small_signal.points[k];
        if (swept >= p.dc.bias_V.size() || frequency >= p.admittance.size()) {
            return std::unexpected(invalid("contact out of range", k));
        }
        x.push_back(p.dc.bias_V[swept]);
        y.push_back(small_signal.capacitance(k, frequency, row, column));
    }
    return make_curve(std::move(x), std::move(y));
}

Curve scaled(Curve curve, double factor) {
    for (double& v : curve.y) v *= factor;
    for (double& v : curve.resolution) v *= std::abs(factor);
    return curve;
}

std::expected<Curve, base::Error> slice(const Curve& curve, double from_x, double to_x) {
    Curve out;
    for (std::size_t k = 0; k < curve.x.size(); ++k) {
        if (!detail::between(curve.x[k], from_x, to_x)) continue;
        out.x.push_back(curve.x[k]);
        out.y.push_back(curve.y[k]);
        if (!curve.resolution.empty()) out.resolution.push_back(curve.resolution[k]);
    }
    if (out.x.size() < 2) return std::unexpected(invalid("fewer than 2 points in the range"));
    return out;
}

bool strictly_monotone(std::span<const double> x) noexcept {
    if (x.size() < 2) return false;
    const bool rising = x[1] > x[0];
    for (std::size_t k = 0; k + 1 < x.size(); ++k) {
        if (rising ? !(x[k + 1] > x[k]) : !(x[k + 1] < x[k])) return false;
    }
    return true;
}

std::expected<std::vector<double>, base::Error> derivative(const Curve& curve) {
    const std::vector<double>& x = curve.x;
    const std::vector<double>& y = curve.y;
    if (!strictly_monotone(x)) return std::unexpected(invalid("x is not strictly monotone"));
    const std::size_t n = x.size();
    std::vector<double> d(n);
    d[0] = (y[1] - y[0]) / (x[1] - x[0]);
    d[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
    for (std::size_t k = 1; k + 1 < n; ++k) {
        const double hl = x[k] - x[k - 1], hr = x[k + 1] - x[k];
        d[k] = (hl * hl * y[k + 1] - hr * hr * y[k - 1] + (hr * hr - hl * hl) * y[k]) /
               (hl * hr * (hl + hr));
    }
    return d;
}

std::expected<Extraction, base::Error> value_at(const Curve& curve, double x, Scale scale) {
    if (!strictly_monotone(curve.x)) return std::unexpected(invalid("x is not strictly monotone"));
    const std::size_t n = curve.x.size();
    for (std::size_t k = 0; k + 1 < n; ++k) {
        const double x0 = curve.x[k], x1 = curve.x[k + 1];
        if (!(std::min(x0, x1) <= x && x <= std::max(x0, x1))) continue;
        // On a point: that point alone.
        if (x == x0 || x == x1) {
            const std::size_t j = x == x0 ? k : k + 1;
            return Extraction{.value = curve.y[j], .unit = {}, .window = {j, j},
                              .extrapolated = false, .below_resolution = below(curve, j, j)};
        }
        const double t = (x - x0) / (x1 - x0);
        const double y0 = curve.y[k], y1 = curve.y[k + 1];
        double v = 0.0;
        if (scale == Scale::linear) {
            v = y0 + t * (y1 - y0);
        } else {
            if (y0 == 0.0 || y1 == 0.0 || (y0 > 0.0) != (y1 > 0.0)) {
                return std::unexpected(
                    invalid("logarithmic interpolation across zero or a sign change", k));
            }
            const double sign = y0 > 0.0 ? 1.0 : -1.0;
            v = sign * std::exp(std::log(std::abs(y0)) +
                                t * (std::log(std::abs(y1)) - std::log(std::abs(y0))));
        }
        return Extraction{.value = v,
                          .unit = {},
                          .window = {k, k + 1},
                          .extrapolated = false,
                          .below_resolution = below(curve, k, k + 1)};
    }
    return std::unexpected(invalid("x outside the curve"));
}

std::expected<Extraction, base::Error> crossing(const Curve& curve, double target, Scale scale) {
    const std::vector<double>& x = curve.x;
    const std::vector<double>& y = curve.y;
    if (y[0] == target) {
        return Extraction{.value = x[0], .unit = "V", .window = {0, 0},
                          .extrapolated = false, .below_resolution = below(curve, 0, 0)};
    }
    for (std::size_t k = 0; k + 1 < x.size(); ++k) {
        const double a = y[k] - target, b = y[k + 1] - target;
        if (!(b == 0.0 || (a < 0.0) != (b < 0.0))) continue;
        double t = 0.0;
        if (scale == Scale::linear) {
            t = (target - y[k]) / (y[k + 1] - y[k]);
        } else {
            if (target == 0.0 || y[k] == 0.0 || (y[k] > 0.0) != (target > 0.0) ||
                (y[k + 1] > 0.0) != (target > 0.0)) {
                return std::unexpected(
                    invalid("logarithmic interpolation across zero or a sign change", k));
            }
            const double l0 = std::log(std::abs(y[k]));
            t = (std::log(std::abs(target)) - l0) / (std::log(std::abs(y[k + 1])) - l0);
        }
        return Extraction{.value = x[k] + t * (x[k + 1] - x[k]),
                          .unit = "V",
                          .window = {k, k + 1},
                          .extrapolated = false,
                          .below_resolution = below(curve, k, k + 1)};
    }
    return std::unexpected(invalid("the curve never reaches the target"));
}

}  // namespace NiTCAD::analysis
