#include "NiTCAD/solve/waveform.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <string>

namespace NiTCAD::solve {

namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt,
                    std::optional<double> value = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = value}};
}

}  // namespace

Waveform Waveform::constant(double value_V) {
    Waveform w;
    w.corners_ = {{0.0, value_V}};
    return w;
}

std::expected<Waveform, base::Error> Waveform::piecewise_linear(
    std::vector<std::pair<double, double>> corners) {
    if (corners.empty()) return std::unexpected(invalid("a waveform needs at least one corner"));
    for (std::size_t k = 0; k < corners.size(); ++k) {
        const auto [t, v] = corners[k];
        if (!std::isfinite(t)) return std::unexpected(invalid("waveform time not finite", k, t));
        if (!std::isfinite(v)) return std::unexpected(invalid("waveform value not finite", k, v));
        if (k == 0) continue;
        if (t < corners[k - 1].first) {
            return std::unexpected(invalid("waveform times must not decrease", k, t));
        }
        if (k >= 2 && t == corners[k - 2].first) {
            return std::unexpected(invalid("at most two waveform corners at one time", k, t));
        }
    }
    Waveform w;
    w.corners_ = std::move(corners);
    return w;
}

std::expected<Waveform, base::Error> Waveform::step(double v0, double v1, double t) {
    return piecewise_linear({{t, v0}, {t, v1}});
}

std::expected<Waveform, base::Error> Waveform::ramp(double v0, double v1, double t0, double t1) {
    if (!(t1 > t0)) return std::unexpected(invalid("a ramp needs t1 > t0", std::nullopt, t1));
    return piecewise_linear({{t0, v0}, {t1, v1}});
}

std::expected<Waveform, base::Error> Waveform::pulse(double v_base, double v_pulse,
                                                     double t_start, double width, double rise,
                                                     double fall) {
    if (!(width > 0.0)) {
        return std::unexpected(invalid("a pulse needs a positive width", std::nullopt, width));
    }
    if (!(rise >= 0.0) || !(fall >= 0.0)) {
        return std::unexpected(invalid("pulse rise and fall times must not be negative"));
    }
    const double t1 = t_start + rise, t2 = t1 + width, t3 = t2 + fall;
    return piecewise_linear({{t_start, v_base}, {t1, v_pulse}, {t2, v_pulse}, {t3, v_base}});
}

std::expected<Waveform, base::Error> Waveform::sine(double offset_V, double amplitude_V,
                                                    double frequency_Hz, double delay_s,
                                                    double phase_rad) {
    for (const double v : {offset_V, amplitude_V, frequency_Hz, delay_s, phase_rad}) {
        if (!std::isfinite(v)) {
            return std::unexpected(invalid("sine parameter not finite", std::nullopt, v));
        }
    }
    if (!(frequency_Hz > 0.0)) {
        return std::unexpected(
            invalid("sine frequency must be positive", std::nullopt, frequency_Hz));
    }
    Waveform w;
    w.kind_ = Kind::sine;
    w.offset_ = offset_V;
    w.amplitude_ = amplitude_V;
    w.frequency_ = frequency_Hz;
    w.delay_ = delay_s;
    w.phase_ = phase_rad;
    return w;
}

double Waveform::value(double t) const noexcept {
    if (kind_ == Kind::sine) {
        const double arg =
            t < delay_ ? phase_ : 2.0 * std::numbers::pi * frequency_ * (t - delay_) + phase_;
        return offset_ + amplitude_ * std::sin(arg);
    }
    // The last corner at or before t (the second of a jump), else the first value.
    const auto after = std::upper_bound(corners_.begin(), corners_.end(), t,
                                        [](double t, const auto& c) { return t < c.first; });
    if (after == corners_.begin()) return corners_.front().second;
    if (after == corners_.end()) return corners_.back().second;
    const auto& [t0, v0] = *(after - 1);
    const auto& [t1, v1] = *after;
    return v0 + (v1 - v0) * ((t - t0) / (t1 - t0));
}

double Waveform::left_value(double t) const noexcept {
    if (kind_ == Kind::sine) return value(t);  // continuous
    // The first corner at or after t (the first of a jump), else the last value.
    const auto at = std::lower_bound(corners_.begin(), corners_.end(), t,
                                     [](const auto& c, double t) { return c.first < t; });
    if (at == corners_.end()) return corners_.back().second;
    if (at == corners_.begin()) return corners_.front().second;
    const auto& [t0, v0] = *(at - 1);
    const auto& [t1, v1] = *at;
    if (t1 == t) return v1;
    return v0 + (v1 - v0) * ((t - t0) / (t1 - t0));
}

std::vector<double> Waveform::breakpoints(double t0, double t1) const {
    std::vector<double> b;
    if (kind_ == Kind::sine) {
        if (delay_ > t0 && delay_ <= t1) b.push_back(delay_);
        return b;
    }
    for (const auto& [t, v] : corners_) {
        if (t > t0 && t <= t1 && (b.empty() || b.back() != t)) b.push_back(t);
    }
    return b;
}

std::vector<double> Waveform::parameters() const {
    if (kind_ == Kind::sine) return {offset_, amplitude_, frequency_, delay_, phase_};
    std::vector<double> p;
    for (const auto& [t, v] : corners_) {
        p.push_back(t);
        p.push_back(v);
    }
    return p;
}

}  // namespace NiTCAD::solve
