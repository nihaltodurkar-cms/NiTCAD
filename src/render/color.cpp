#include "NiTCAD/render/color.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <string>
#include <utility>

namespace NiTCAD::render {

namespace {

struct Entry {
    double r, g, b;
};

#include "colormap_data.inc"

// The workspace reference's softer map (VQ14): deep blue, teal, green, gold, muted red. The page
// names the colors without values; these are the values chosen for them (DECISIONS.md, Unit 25a).
constexpr std::array<Entry, 5> soft_table{{
    {0.145, 0.247, 0.561},  // deep blue  #253f8f
    {0.165, 0.616, 0.612},  // teal       #2a9d9c
    {0.427, 0.733, 0.353},  // green      #6dbb5a
    {0.890, 0.698, 0.235},  // gold       #e3b23c
    {0.753, 0.314, 0.302},  // muted red  #c0504d
}};

template <std::size_t N>
Rgb lookup(const std::array<Entry, N>& table, double t) noexcept {
    const double x = std::clamp(t, 0.0, 1.0) * static_cast<double>(N - 1);
    const auto i = std::min(static_cast<std::size_t>(x), N - 2);
    const double f = x - static_cast<double>(i);
    const Entry& a = table[i];
    const Entry& b = table[i + 1];
    return {static_cast<float>(a.r + f * (b.r - a.r)), static_cast<float>(a.g + f * (b.g - a.g)),
            static_cast<float>(a.b + f * (b.b - a.b))};
}

base::Error invalid(std::string message, double value) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = std::nullopt, .value = value}};
}

// UTF-8 of the minus sign U+2212 and of the superscripts.
constexpr std::string_view minus = "\xE2\x88\x92";
constexpr std::array<std::string_view, 10> superscript{
    "\xE2\x81\xB0", "\xC2\xB9",     "\xC2\xB2",     "\xC2\xB3",     "\xE2\x81\xB4",
    "\xE2\x81\xB5", "\xE2\x81\xB6", "\xE2\x81\xB7", "\xE2\x81\xB8", "\xE2\x81\xB9"};
constexpr std::string_view superscript_minus = "\xE2\x81\xBB";

std::string with_minus(std::string text) {
    if (!text.empty() && text.front() == '-') text.replace(0, 1, minus);
    return text;
}

// Integer k such that 10^k <= x < 10^(k+1), robust at exact powers.
int decade_of(double x) noexcept {
    int k = static_cast<int>(std::floor(std::log10(x)));
    if (std::pow(10.0, k + 1) <= x) ++k;
    if (std::pow(10.0, k) > x) --k;
    return k;
}

}  // namespace

Rgb sample(Colormap map, double t) noexcept {
    if (std::isnan(t)) return no_data;
    switch (map) {
        case Colormap::viridis: return lookup(viridis_table, t);
        case Colormap::cividis: return lookup(cividis_table, t);
        case Colormap::coolwarm: return lookup(coolwarm_table, t);
        case Colormap::soft: return lookup(soft_table, t);
        case Colormap::turbo:
        default: return lookup(turbo_table, t);
    }
}

std::string_view colormap_name(Colormap map) noexcept {
    switch (map) {
        case Colormap::viridis: return "viridis";
        case Colormap::cividis: return "cividis";
        case Colormap::coolwarm: return "coolwarm";
        case Colormap::soft: return "soft";
        case Colormap::turbo:
        default: return "turbo";
    }
}

std::optional<Colormap> find_colormap(std::string_view name) noexcept {
    for (const Colormap m : {Colormap::turbo, Colormap::viridis, Colormap::cividis,
                             Colormap::coolwarm, Colormap::soft}) {
        if (colormap_name(m) == name) return m;
    }
    return std::nullopt;
}

// --- ColorScale -------------------------------------------------------------------------------

ColorScale::ColorScale(ScaleKind kind, double lo, double hi, double floor) noexcept
    : kind_(kind), lo_(lo), hi_(hi), floor_(floor), s_lo_(0.0), s_hi_(0.0) {
    s_lo_ = transform(lo);
    s_hi_ = transform(hi);
}

std::expected<ColorScale, base::Error> ColorScale::create(ScaleKind kind, double lo, double hi,
                                                          double floor) {
    if (!std::isfinite(lo) || !std::isfinite(hi) || !(lo < hi)) {
        return std::unexpected(invalid("a scale needs finite lo < hi", lo));
    }
    if (kind == ScaleKind::log10 && !(lo > 0.0)) {
        return std::unexpected(invalid("a log scale needs lo > 0", lo));
    }
    if (kind == ScaleKind::signed_log10 && !(std::isfinite(floor) && floor > 0.0)) {
        return std::unexpected(invalid("a signed log scale needs a positive floor", floor));
    }
    if (kind != ScaleKind::linear && kind != ScaleKind::log10 && kind != ScaleKind::signed_log10) {
        return std::unexpected(invalid("unknown scale", static_cast<double>(kind)));
    }
    return ColorScale(kind, lo, hi, floor);
}

std::expected<ColorScale, base::Error> ColorScale::fit(ScaleKind kind,
                                                       std::span<const double> values,
                                                       double floor) {
    double lo = std::numeric_limits<double>::infinity();
    double hi = -lo;
    for (const double v : values) {
        if (!std::isfinite(v) || (kind == ScaleKind::log10 && !(v > 0.0))) continue;
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    if (!(lo <= hi)) return std::unexpected(invalid("no value to fit a scale to", 0.0));
    switch (kind) {
        case ScaleKind::log10: {
            const int a = decade_of(lo);
            int b = decade_of(hi);
            if (std::pow(10.0, b) < hi) ++b;
            if (a == b) return create(kind, std::pow(10.0, a - 1), std::pow(10.0, b + 1), floor);
            return create(kind, std::pow(10.0, a), std::pow(10.0, b), floor);
        }
        case ScaleKind::signed_log10:
            if (lo == hi) {
                const double m = std::max(std::abs(lo) * 10.0, floor * 10.0);
                return create(kind, lo < 0.0 ? -m : std::min(lo / 10.0, -floor), std::max(m, lo),
                              floor);
            }
            return create(kind, lo, hi, floor);
        case ScaleKind::linear:
        default:
            if (lo == hi) {
                const double w = std::max(std::abs(lo) / 2.0, 1.0);
                return create(kind, lo - w, hi + w, floor);
            }
            return create(kind, lo, hi, floor);
    }
}

double ColorScale::transform(double v) const noexcept {
    switch (kind_) {
        case ScaleKind::log10: return v > 0.0 ? std::log10(v) : std::numeric_limits<double>::quiet_NaN();
        case ScaleKind::signed_log10: return std::copysign(std::log10(1.0 + std::abs(v) / floor_), v);
        case ScaleKind::linear:
        default: return v;
    }
}

double ColorScale::inverse(double s) const noexcept {
    switch (kind_) {
        case ScaleKind::log10: return std::pow(10.0, s);
        case ScaleKind::signed_log10: return std::copysign(floor_ * (std::pow(10.0, std::abs(s)) - 1.0), s);
        case ScaleKind::linear:
        default: return s;
    }
}

double ColorScale::position(double v) const noexcept {
    const double s = transform(v);
    if (std::isnan(s)) return s;
    return std::clamp((s - s_lo_) / (s_hi_ - s_lo_), 0.0, 1.0);
}

double ColorScale::value_at(double t) const noexcept {
    return inverse(s_lo_ + std::clamp(t, 0.0, 1.0) * (s_hi_ - s_lo_));
}

std::vector<Tick> ColorScale::ticks(int target) const {
    switch (kind_) {
        case ScaleKind::log10: return log_ticks(lo_, hi_, target);
        case ScaleKind::linear: return linear_ticks(lo_, hi_, target);
        case ScaleKind::signed_log10:
        default: break;
    }
    // 0 and the decades of |v| at least 10 times the floor, every k-th when there are too many.
    const int first = decade_of(floor_ * 10.0);
    const int last = decade_of(std::max(std::abs(lo_), std::abs(hi_)));
    const int count = std::max(0, last - first + 1);
    const int sides = (lo_ < 0.0 ? 1 : 0) + (hi_ > 0.0 ? 1 : 0);
    const int stride = std::max(1, (count * sides + 1 + target - 1) / std::max(target, 1));
    std::vector<Tick> out;
    const auto add = [&](double v, std::string label) {
        if (v < lo_ || v > hi_) return;
        out.push_back({v, position(v), std::move(label), true});
    };
    for (int k = last; k >= first; --k) {
        if ((k - first) % stride == 0) add(-std::pow(10.0, k), std::string(minus) + power_of_ten(k));
    }
    add(0.0, "0");
    for (int k = first; k <= last; ++k) {
        if ((k - first) % stride == 0) add(std::pow(10.0, k), power_of_ten(k));
    }
    return out;
}

// --- ticks and labels -------------------------------------------------------------------------

std::vector<Tick> linear_ticks(double lo, double hi, int target) {
    std::vector<Tick> out;
    if (!(std::isfinite(lo) && std::isfinite(hi) && lo < hi)) return out;
    target = std::max(target, 2);
    const double span = hi - lo;
    const double raw = span / static_cast<double>(target);
    const double magnitude = std::pow(10.0, decade_of(raw));
    double step = magnitude;
    double best = std::numeric_limits<double>::infinity();
    for (const double m : {1.0, 2.0, 5.0, 10.0}) {
        const double s = m * magnitude;
        const double miss = std::abs(span / s - static_cast<double>(target));
        if (miss < best) {
            best = miss;
            step = s;
        }
    }
    const double scale = std::max(std::abs(lo), std::abs(hi));
    const auto first = static_cast<long long>(std::ceil(lo / step - 1e-9));
    const auto last = static_cast<long long>(std::floor(hi / step + 1e-9));
    for (long long k = first; k <= last; ++k) {
        double v = static_cast<double>(k) * step;
        if (std::abs(v) < step * 1e-9) v = 0.0;
        out.push_back({v, std::clamp((v - lo) / span, 0.0, 1.0), format_tick(v, step, scale), true});
    }
    return out;
}

std::vector<Tick> log_ticks(double lo, double hi, int target) {
    std::vector<Tick> out;
    if (!(std::isfinite(lo) && std::isfinite(hi) && lo > 0.0 && lo < hi)) return out;
    target = std::max(target, 2);
    const double l0 = std::log10(lo);
    const double l1 = std::log10(hi);
    const auto position = [&](double logv) { return std::clamp((logv - l0) / (l1 - l0), 0.0, 1.0); };
    int first = decade_of(lo);
    if (std::pow(10.0, first) < lo) ++first;
    const int last = decade_of(hi);
    const int count = std::max(0, last - first + 1);
    const int stride = std::max(1, (count + target - 1) / target);
    for (int k = first; k <= last; ++k) {
        if (((k % stride) + stride) % stride != 0) continue;
        out.push_back({std::pow(10.0, k), position(static_cast<double>(k)), power_of_ten(k), true});
    }
    if (l1 - l0 <= 3.0 + 1e-12) {
        for (int k = decade_of(lo); k <= last; ++k) {
            for (int m = 2; m <= 9; ++m) {
                const double v = m * std::pow(10.0, k);
                if (v < lo || v > hi) continue;
                out.push_back({v, position(std::log10(v)), std::string(), false});
            }
        }
        std::ranges::sort(out, {}, &Tick::value);
    }
    return out;
}

std::string format_tick(double value, double step, double scale) {
    if (value == 0.0) return "0";
    const bool scientific = scale >= 1e4 || scale < 1e-3;
    if (!scientific) {
        const int decimals = std::max(0, -decade_of(step * (1.0 + 1e-9)));
        return with_minus(std::format("{:.{}f}", value, decimals));
    }
    const int e = decade_of(std::abs(value) * (1.0 + 1e-12));
    const double mantissa = value / std::pow(10.0, e);
    const int decimals = std::max(0, e - decade_of(step * (1.0 + 1e-9)));
    std::string m = std::format("{:.{}f}", mantissa, std::min(decimals, 6));
    if (m == "1" || m == "-1") return (m == "-1" ? std::string(minus) : std::string()) + power_of_ten(e);
    return with_minus(std::move(m)) + "\xC3\x97" + power_of_ten(e);  // U+00D7 multiplication sign
}

std::string power_of_ten(int k) {
    std::string out = "10";
    if (k < 0) out += superscript_minus;
    const std::string digits = std::to_string(std::abs(k));
    for (const char c : digits) out += superscript[static_cast<std::size_t>(c - '0')];
    return out;
}

}  // namespace NiTCAD::render
