#include "NiTCAD/analysis/dc.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "common.hpp"

namespace NiTCAD::analysis {

namespace {

using detail::below;
using detail::invalid;

struct Line {
    double slope;
    double intercept;
};

// Least squares through (u, v), about their means.
Line fit(const std::vector<double>& u, const std::vector<double>& v) {
    const double n = static_cast<double>(u.size());
    double mu = 0.0, mv = 0.0;
    for (std::size_t k = 0; k < u.size(); ++k) {
        mu += u[k];
        mv += v[k];
    }
    mu /= n;
    mv /= n;
    double suu = 0.0, suv = 0.0;
    for (std::size_t k = 0; k < u.size(); ++k) {
        suu += (u[k] - mu) * (u[k] - mu);
        suv += (u[k] - mu) * (v[k] - mv);
    }
    const double slope = suv / suu;
    return {slope, mv - slope * mu};
}

// The indices of the points with x between the two bounds, which must be consecutive.
std::expected<Window, base::Error> window_of(const Curve& c, double from_x, double to_x,
                                             std::size_t minimum) {
    std::optional<std::size_t> first;
    std::size_t last = 0, count = 0;
    for (std::size_t k = 0; k < c.x.size(); ++k) {
        if (!detail::between(c.x[k], from_x, to_x)) continue;
        if (first && k != last + 1) return std::unexpected(invalid("the window is not contiguous", k));
        if (!first) first = k;
        last = k;
        ++count;
    }
    if (count < minimum) return std::unexpected(invalid("too few points in the window"));
    return Window{*first, last};
}

struct Peak {
    std::size_t point;
    Extraction gm;
};

// The largest transconductance (the first, on a tie, as numpy.argmax) and the point it is at.
std::expected<Peak, base::Error> peak_transconductance(const Curve& transfer) {
    const auto gm = derivative(transfer);
    if (!gm) return std::unexpected(gm.error());
    const auto k = static_cast<std::size_t>(std::max_element(gm->begin(), gm->end()) - gm->begin());
    if (!((*gm)[k] > 0.0)) return std::unexpected(invalid("the transconductance is nowhere positive"));
    const Window w{k == 0 ? 0 : k - 1, std::min(k + 1, transfer.x.size() - 1)};
    return Peak{k, {.value = (*gm)[k],
                    .unit = "A cm^(D-3)/V",
                    .window = w,
                    .extrapolated = false,
                    .below_resolution = below(transfer, w.first, w.last)}};
}

}  // namespace

std::expected<Extraction, base::Error> max_transconductance(const Curve& transfer) {
    const auto p = peak_transconductance(transfer);
    if (!p) return std::unexpected(p.error());
    return p->gm;
}

std::expected<Extraction, base::Error> threshold_max_gm(const Curve& transfer, double drain_V) {
    const auto p = peak_transconductance(transfer);
    if (!p) return std::unexpected(p.error());
    Extraction e = p->gm;
    const std::size_t k = p->point;
    e.value = transfer.x[k] - transfer.y[k] / p->gm.value - 0.5 * drain_V;
    e.unit = "V";
    e.extrapolated = true;  // the tangent's zero lies outside the points it came from
    return e;
}

std::expected<Extraction, base::Error> threshold_constant_current(const Curve& transfer,
                                                                  double criterion) {
    if (!(std::isfinite(criterion) && criterion > 0.0)) {
        return std::unexpected(invalid("the criterion is not finite and positive"));
    }
    // Only a rising crossing from positive current counts: find the first point at or above.
    for (std::size_t k = 0; k < transfer.y.size(); ++k) {
        if (transfer.y[k] < criterion) continue;
        if (k == 0) {
            return Extraction{.value = transfer.x[0], .unit = "V", .window = {0, 0},
                              .extrapolated = false, .below_resolution = below(transfer, 0, 0)};
        }
        const double y0 = transfer.y[k - 1], y1 = transfer.y[k];
        double t = 0.0;
        if (y0 > 0.0) {
            t = std::log(criterion / y0) / std::log(y1 / y0);
        } else {
            t = (criterion - y0) / (y1 - y0);  // from zero or negative current: linear
        }
        return Extraction{.value = transfer.x[k - 1] + t * (transfer.x[k] - transfer.x[k - 1]),
                          .unit = "V",
                          .window = {k - 1, k},
                          .extrapolated = false,
                          .below_resolution = below(transfer, k - 1, k)};
    }
    return std::unexpected(invalid("the current never reaches the criterion"));
}

std::expected<Extraction, base::Error> subthreshold_swing(const Curve& transfer,
                                                          double ceiling_fraction) {
    if (!(ceiling_fraction > 0.0 && ceiling_fraction <= 1.0)) {
        return std::unexpected(invalid("ceiling_fraction is not in (0, 1]"));
    }
    const double ceiling = ceiling_fraction * *std::max_element(transfer.y.begin(), transfer.y.end());
    const auto usable = [&](std::size_t k) {
        const double r = transfer.resolution.empty() ? 0.0 : transfer.resolution[k];
        return transfer.y[k] > 0.0 && transfer.y[k] > r && transfer.y[k] <= ceiling;
    };
    std::optional<std::size_t> best;
    double swing = 0.0;
    for (std::size_t k = 0; k + 1 < transfer.y.size(); ++k) {
        if (!usable(k) || !usable(k + 1)) continue;
        const double s = (transfer.x[k + 1] - transfer.x[k]) /
                         (std::log10(transfer.y[k + 1]) - std::log10(transfer.y[k]));
        if (!(s > 0.0) || !std::isfinite(s)) continue;  // the current must rise with the gate
        if (!best || s < swing) {
            best = k;
            swing = s;
        }
    }
    if (!best) return std::unexpected(invalid("no subthreshold segment with a rising current"));
    return Extraction{.value = 1000.0 * swing,
                      .unit = "mV/decade",
                      .window = {*best, *best + 1},
                      .extrapolated = false,
                      .below_resolution = false};
}

std::expected<Extraction, base::Error> dibl(const Curve& low_drain, const Curve& high_drain,
                                            double low_drain_V, double high_drain_V,
                                            double criterion) {
    if (!(high_drain_V != low_drain_V)) return std::unexpected(invalid("the drain biases are equal"));
    const auto low = threshold_constant_current(low_drain, criterion);
    if (!low) return low;
    auto high = threshold_constant_current(high_drain, criterion);
    if (!high) return high;
    high->below_resolution = high->below_resolution || low->below_resolution;
    high->value = (low->value - high->value) / (high_drain_V - low_drain_V);
    high->unit = "V/V";
    return high;
}

std::expected<Extraction, base::Error> on_off_ratio(const Curve& transfer, double on_V,
                                                    double off_V) {
    const auto on = value_at(transfer, on_V, Scale::logarithmic);
    if (!on) return on;
    const auto off = value_at(transfer, off_V, Scale::logarithmic);
    if (!off) return off;
    if ((on->value > 0.0) != (off->value > 0.0)) {
        return std::unexpected(invalid("the on and off currents have opposite signs"));
    }
    return Extraction{.value = on->value / off->value,
                      .unit = "",
                      .window = {std::min(on->window.first, off->window.first),
                                 std::max(on->window.last, off->window.last)},
                      .extrapolated = false,
                      .below_resolution = on->below_resolution || off->below_resolution};
}

std::expected<Extraction, base::Error> output_conductance(const Curve& output,
                                                          double fit_fraction) {
    if (!(fit_fraction > 0.0 && fit_fraction <= 1.0)) {
        return std::unexpected(invalid("fit_fraction is not in (0, 1]"));
    }
    if (!strictly_monotone(output.x)) return std::unexpected(invalid("x is not strictly monotone"));
    const double end = output.x.back();
    const double span = std::abs(end - output.x.front());
    const auto w = window_of(output, end, end - std::copysign(fit_fraction * span,
                                                              end - output.x.front()), 3);
    if (!w) return std::unexpected(w.error());
    const std::vector<double> u(output.x.begin() + static_cast<std::ptrdiff_t>(w->first),
                                output.x.begin() + static_cast<std::ptrdiff_t>(w->last) + 1);
    const std::vector<double> v(output.y.begin() + static_cast<std::ptrdiff_t>(w->first),
                                output.y.begin() + static_cast<std::ptrdiff_t>(w->last) + 1);
    return Extraction{.value = fit(u, v).slope,
                      .unit = "A cm^(D-3)/V",
                      .window = *w,
                      .extrapolated = false,
                      .below_resolution = below(output, w->first, w->last)};
}

std::expected<DiodeFit, base::Error> diode_fit(const Curve& forward, double from_V, double to_V,
                                               double temperature_K) {
    if (!(temperature_K > 0.0)) return std::unexpected(invalid("the temperature is not positive"));
    const auto w = window_of(forward, from_V, to_V, 2);
    if (!w) return std::unexpected(w.error());
    std::vector<double> u, v;
    for (std::size_t k = w->first; k <= w->last; ++k) {
        if (!(forward.y[k] > 0.0)) return std::unexpected(invalid("a current is not positive", k));
        u.push_back(forward.x[k]);
        v.push_back(std::log(forward.y[k]));
    }
    const Line line = fit(u, v);
    if (!(line.slope > 0.0)) return std::unexpected(invalid("ln I does not rise with the bias"));
    const bool noise = below(forward, w->first, w->last);
    const double V_T = base::thermal_voltage(temperature_K);
    return DiodeFit{
        .ideality = {.value = 1.0 / (line.slope * V_T), .unit = "", .window = *w,
                     .extrapolated = false, .below_resolution = noise},
        .saturation_current = {.value = std::exp(line.intercept), .unit = "A cm^(D-3)",
                               .window = *w, .extrapolated = true, .below_resolution = noise}};
}

std::expected<std::vector<double>, base::Error> local_ideality(const Curve& forward,
                                                               double temperature_K) {
    if (!(temperature_K > 0.0)) return std::unexpected(invalid("the temperature is not positive"));
    Curve log_curve{forward.x, forward.y, {}};
    for (std::size_t k = 0; k < log_curve.y.size(); ++k) {
        if (!(forward.y[k] > 0.0)) return std::unexpected(invalid("a current is not positive", k));
        log_curve.y[k] = std::log(forward.y[k]);
    }
    auto d = derivative(log_curve);
    if (!d) return d;
    const double V_T = base::thermal_voltage(temperature_K);
    for (double& v : *d) v = 1.0 / (V_T * v);
    return d;
}

std::expected<SeriesResistanceFit, base::Error> series_resistance(const Curve& forward,
                                                                  double from_V, double to_V,
                                                                  double temperature_K) {
    if (!(temperature_K > 0.0)) return std::unexpected(invalid("the temperature is not positive"));
    const auto w = window_of(forward, from_V, to_V, 3);
    if (!w) return std::unexpected(w.error());
    // The curve's end points have only a one-sided (first-order) derivative: left out. The others
    // need their two neighbours.
    const std::size_t first = std::max<std::size_t>(w->first, 1);
    const std::size_t last = std::min(w->last, forward.y.size() - 2);
    if (last < first + 2) return std::unexpected(invalid("too few interior points in the window"));
    // V against ln I, differentiated directly: dV / d ln I is the quantity fitted.
    Curve by_log;
    for (std::size_t k = first - 1; k <= last + 1; ++k) {
        if (!(forward.y[k] > 0.0)) return std::unexpected(invalid("a current is not positive", k));
        by_log.x.push_back(std::log(forward.y[k]));
        by_log.y.push_back(forward.x[k]);
    }
    const auto slope = derivative(by_log);
    if (!slope) return std::unexpected(invalid("ln I is not strictly monotone in the bias"));
    const double V_T = base::thermal_voltage(temperature_K);
    std::vector<double> u, v;
    for (std::size_t k = first; k <= last; ++k) {
        u.push_back(forward.y[k]);
        v.push_back((*slope)[k - first + 1]);
    }
    const Line line = fit(u, v);
    const Window used{first, last};
    const bool noise = below(forward, first, last);
    return SeriesResistanceFit{
        .ideality = {.value = line.intercept / V_T, .unit = "", .window = used,
                     .extrapolated = true, .below_resolution = noise},
        .resistance = {.value = line.slope, .unit = "Ohm cm^(3-D)", .window = used,
                       .extrapolated = false, .below_resolution = noise}};
}

std::expected<Extraction, base::Error> breakdown_at_current(const Curve& curve, double limit) {
    if (!(std::isfinite(limit) && limit > 0.0)) {
        return std::unexpected(invalid("the limit is not finite and positive"));
    }
    for (std::size_t k = 0; k < curve.y.size(); ++k) {
        const double a1 = std::abs(curve.y[k]);
        if (a1 < limit) continue;
        if (k == 0) {
            return Extraction{.value = curve.x[0], .unit = "V", .window = {0, 0},
                              .extrapolated = false, .below_resolution = below(curve, 0, 0)};
        }
        const double a0 = std::abs(curve.y[k - 1]);
        const double t = a0 > 0.0 ? std::log(limit / a0) / std::log(a1 / a0)
                                  : (limit - a0) / (a1 - a0);
        return Extraction{.value = curve.x[k - 1] + t * (curve.x[k] - curve.x[k - 1]),
                          .unit = "V",
                          .window = {k - 1, k},
                          .extrapolated = false,
                          .below_resolution = below(curve, k - 1, k)};
    }
    return std::unexpected(invalid("|I| never reaches the limit"));
}

std::expected<Extraction, base::Error> turning_point(const Curve& curve) {
    const std::vector<double>& x = curve.x;
    std::size_t start = 1;
    while (start < x.size() && x[start] == x[0]) ++start;
    if (start == x.size()) return std::unexpected(invalid("the bias never moves"));
    const double direction = x[start] > x[0] ? 1.0 : -1.0;
    for (std::size_t k = start; k + 1 < x.size(); ++k) {
        if (!((x[k + 1] - x[k]) * direction < 0.0)) continue;
        Extraction e{.value = x[k], .unit = "V", .window = {k, k},
                     .extrapolated = false, .below_resolution = below(curve, k, k)};
        // V(t), t = ln |I|, through the three points: its vertex.
        const double y0 = std::abs(curve.y[k - 1]), y1 = std::abs(curve.y[k]),
                     y2 = std::abs(curve.y[k + 1]);
        if (y0 > 0.0 && y1 > 0.0 && y2 > 0.0) {
            const double t0 = std::log(y0), t1 = std::log(y1), t2 = std::log(y2);
            if ((t1 - t0) * (t2 - t1) > 0.0) {
                const double d01 = (x[k] - x[k - 1]) / (t1 - t0);
                const double d12 = (x[k + 1] - x[k]) / (t2 - t1);
                const double a = (d12 - d01) / (t2 - t0);
                const double b = d01 - a * (t0 + t1);
                if (a != 0.0) {
                    const double tv = -b / (2.0 * a);
                    if (std::min(t0, t2) <= tv && tv <= std::max(t0, t2)) {
                        e.value = x[k - 1] + (tv - t0) * (d01 + a * (tv - t1));
                        e.window = {k - 1, k + 1};
                        e.below_resolution = below(curve, k - 1, k + 1);
                    }
                }
            }
        }
        return e;
    }
    return std::unexpected(invalid("the bias never turns back"));
}

}  // namespace NiTCAD::analysis
