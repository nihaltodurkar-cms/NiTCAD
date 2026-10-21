// Colors, colormaps, value scales and their ticks for the views (Unit 25a; DECISIONS.md VR3, VQ5,
// VQ14).
//
// A colormap turns a position in [0, 1] into a color; a ColorScale turns a field value into that
// position, on a linear, a log10 or a signed log10 scale. Positions are what the scene stores per
// vertex, so drawing interpolates the position and looks the color up afterwards: a log scale is
// interpolated in log space and a color never comes from blending two colors of the table.
//
// The tables are published data: Turbo (Google, Apache 2.0), viridis and cividis (CC0) as
// matplotlib 3.9.2 ships them in 256 entries, and Moreland's cool-warm diverging map in its 33
// control points. "soft" is the five-color map of the owner's workspace reference (VQ14): deep blue,
// teal, green, gold and muted red, equally spaced. Every table is interpolated linearly in sRGB.
#pragma once

#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::render {

// sRGB components in [0, 1].
struct Rgb {
    float r = 0.0f;
    float g = 0.0f;
    float b = 0.0f;
};

// sRGB with straight (not premultiplied) alpha.
struct Rgba {
    float r = 0.0f;
    float g = 0.0f;
    float b = 0.0f;
    float a = 1.0f;
};

enum class Colormap : std::uint8_t {
    turbo,     // the default (VQ5)
    viridis,
    cividis,
    coolwarm,  // diverging, for signed quantities
    soft,      // the workspace reference's five colors (VQ14)
};

// The color at position t, clamped into [0, 1]; no_data for a NaN position.
[[nodiscard]] Rgb sample(Colormap map, double t) noexcept;
[[nodiscard]] std::string_view colormap_name(Colormap map) noexcept;
// The colormap of this name ("turbo", "viridis", "cividis", "coolwarm", "soft").
[[nodiscard]] std::optional<Colormap> find_colormap(std::string_view name) noexcept;

// The color of a value with no position: NaN, or not positive on a log scale.
inline constexpr Rgb no_data{0.62f, 0.62f, 0.62f};

enum class ScaleKind : std::uint8_t {
    linear,
    log10,         // positive values; a decade is an equal step
    signed_log10,  // sign(v) log10(1 + |v| / floor): signed quantities over many decades (doping)
};

// One tick of an axis or a colorbar: the value, its position along the axis in [0, 1], and its
// label in UTF-8 (powers of ten with superscript digits, the minus sign U+2212).
struct Tick {
    double value = 0.0;
    double position = 0.0;
    std::string label;
    bool major = true;  // minor ticks (2..9 within a decade of a short log axis) have no label
};

class ColorScale {
public:
    // Errors (invalid_input): lo or hi not finite, or lo >= hi; log10 with lo <= 0; signed_log10
    // with a floor that is not positive and finite.
    [[nodiscard]] static std::expected<ColorScale, base::Error> create(ScaleKind kind, double lo,
                                                                      double hi,
                                                                      double floor = 1.0);
    // The range of the finite values (for log10 the positive ones, widened to the whole decades
    // holding them; a single power of ten to one decade each side). A single value v is widened to
    // v +- max(|v| / 2, 1) on a linear scale and to ten times its size on the signed log scale.
    // Errors (invalid_input): no usable value; a floor as in create.
    [[nodiscard]] static std::expected<ColorScale, base::Error> fit(ScaleKind kind,
                                                                   std::span<const double> values,
                                                                   double floor = 1.0);

    // The position of v in [0, 1] (clamped), or NaN for NaN and, on log10, for v <= 0.
    [[nodiscard]] double position(double v) const noexcept;
    // The value at position t in [0, 1] (the inverse of position on the range).
    [[nodiscard]] double value_at(double t) const noexcept;
    // About `target` ticks: linear steps of 1, 2 or 5 times a power of ten; decades on log10 (every
    // k-th decade when there are too many, with the minor 2..9 ticks when there are at most
    // three); 0 and decades of |v| on signed_log10.
    [[nodiscard]] std::vector<Tick> ticks(int target = 6) const;

    [[nodiscard]] ScaleKind kind() const noexcept { return kind_; }
    [[nodiscard]] double lo() const noexcept { return lo_; }
    [[nodiscard]] double hi() const noexcept { return hi_; }
    [[nodiscard]] double floor() const noexcept { return floor_; }

private:
    ColorScale(ScaleKind kind, double lo, double hi, double floor) noexcept;
    [[nodiscard]] double transform(double v) const noexcept;
    [[nodiscard]] double inverse(double s) const noexcept;

    ScaleKind kind_;
    double lo_;
    double hi_;
    double floor_;
    double s_lo_;  // transform(lo_) and transform(hi_)
    double s_hi_;
};

// Ticks of a linear axis over [lo, hi] (lo < hi, finite), steps of 1, 2 or 5 times a power of ten
// chosen so that about `target` ticks fall inside; positions relative to [lo, hi].
[[nodiscard]] std::vector<Tick> linear_ticks(double lo, double hi, int target = 6);
// Ticks of a log10 axis over [lo, hi] (0 < lo < hi, finite); see ColorScale::ticks.
[[nodiscard]] std::vector<Tick> log_ticks(double lo, double hi, int target = 6);
// The label of a linear tick value on a grid of this step, on an axis whose largest |value| is
// `scale`: fixed point with as many decimals as the step needs, or a mantissa times a power of ten
// when the axis reaches 1e4 or stays below 1e-3.
[[nodiscard]] std::string format_tick(double value, double step, double scale);
// 10 to the power k as a label: "10" with superscript digits ("10⁻³").
[[nodiscard]] std::string power_of_ten(int k);

}  // namespace NiTCAD::render
