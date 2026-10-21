// Colormaps, value scales and ticks (Unit 25a, VG3): table endpoints and interpolation, log scales
// exact at decades, ticks of 1, 2 or 5 times a power of ten, and their labels.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "NiTCAD/render/color.hpp"

using namespace NiTCAD::render;
using NiTCAD::base::ErrorCode;

namespace {

bool near(Rgb c, double r, double g, double b, double tol = 1e-6) {
    return std::abs(c.r - r) <= tol && std::abs(c.g - g) <= tol && std::abs(c.b - b) <= tol;
}

}  // namespace

TEST_CASE("color: colormap tables end where the published maps end") {
    // Turbo #30123b .. #7a0402, viridis #440154 .. #fde725, cividis, cool-warm (Moreland).
    CHECK(near(sample(Colormap::turbo, 0.0), 0.18995, 0.07176, 0.23217));
    CHECK(near(sample(Colormap::turbo, 1.0), 0.47960, 0.01583, 0.01055));
    CHECK(near(sample(Colormap::viridis, 0.0), 0.267004, 0.004874, 0.329415));
    CHECK(near(sample(Colormap::viridis, 1.0), 0.993248, 0.906157, 0.143936));
    CHECK(near(sample(Colormap::cividis, 0.0), 0.0, 0.135112, 0.304751));
    CHECK(near(sample(Colormap::cividis, 1.0), 0.995737, 0.909344, 0.217772));
    CHECK(near(sample(Colormap::coolwarm, 0.0), 0.2298057, 0.298717966, 0.753683153));
    CHECK(near(sample(Colormap::coolwarm, 1.0), 0.705673158, 0.01555616, 0.150232812));
    CHECK(near(sample(Colormap::coolwarm, 0.5), 0.865395197, 0.86541021, 0.865395561));
    CHECK(near(sample(Colormap::soft, 0.0), 0.145, 0.247, 0.561));
    CHECK(near(sample(Colormap::soft, 0.5), 0.427, 0.733, 0.353));
    // Clamped outside [0, 1]; NaN is no data.
    CHECK(near(sample(Colormap::turbo, -3.0), 0.18995, 0.07176, 0.23217));
    CHECK(near(sample(Colormap::turbo, 7.0), 0.47960, 0.01583, 0.01055));
    const Rgb nd = sample(Colormap::viridis, std::numeric_limits<double>::quiet_NaN());
    CHECK(near(nd, no_data.r, no_data.g, no_data.b));
}

TEST_CASE("color: tables are interpolated linearly between entries") {
    // Between viridis entries 0 and 1 (t = 1/255): halfway at t = 0.5 / 255.
    const Rgb half = sample(Colormap::viridis, 0.5 / 255.0);
    CHECK(near(half, 0.5 * (0.267004 + 0.268510), 0.5 * (0.004874 + 0.009605),
               0.5 * (0.329415 + 0.335427), 1e-6));
    // Turbo entry 128 exactly at t = 128 / 255.
    const Rgb e = sample(Colormap::turbo, 128.0 / 255.0);
    CHECK(near(e, 0.64362, 0.98999, 0.23356, 1e-6));
}

TEST_CASE("color: names") {
    for (const Colormap m : {Colormap::turbo, Colormap::viridis, Colormap::cividis,
                             Colormap::coolwarm, Colormap::soft}) {
        REQUIRE(find_colormap(colormap_name(m)) == m);
    }
    CHECK_FALSE(find_colormap("jet").has_value());
}

TEST_CASE("color: linear and log scales; log exact at decades") {
    const ColorScale lin = *ColorScale::create(ScaleKind::linear, -2.0, 6.0);
    CHECK(lin.position(-2.0) == 0.0);
    CHECK(lin.position(6.0) == 1.0);
    CHECK(lin.position(2.0) == 0.5);
    CHECK(lin.position(100.0) == 1.0);  // clamped
    CHECK(lin.value_at(0.25) == 0.0);

    const ColorScale log = *ColorScale::create(ScaleKind::log10, 1e2, 1e20);
    for (int k = 2; k <= 20; ++k) {
        CAPTURE(k);
        CHECK(log.position(std::pow(10.0, k)) == static_cast<double>(k - 2) / 18.0);
    }
    CHECK(std::isnan(log.position(0.0)));
    CHECK(std::isnan(log.position(-5.0)));
    CHECK(std::abs(log.value_at(0.5) / 1e11 - 1.0) < 1e-14);

    // Signed log: symmetric about 0, odd.
    const ColorScale s = *ColorScale::create(ScaleKind::signed_log10, -1e20, 1e20, 1e10);
    CHECK(s.position(0.0) == 0.5);
    CHECK(std::abs(s.position(1e15) + s.position(-1e15) - 1.0) < 1e-15);
    CHECK(s.position(1e20) == 1.0);
}

TEST_CASE("color: scale errors and fitting") {
    CHECK(ColorScale::create(ScaleKind::linear, 1.0, 1.0).error().code == ErrorCode::invalid_input);
    CHECK(ColorScale::create(ScaleKind::linear, 0.0, std::nan("")).error().code ==
          ErrorCode::invalid_input);
    CHECK(ColorScale::create(ScaleKind::log10, 0.0, 1.0).error().code == ErrorCode::invalid_input);
    CHECK(ColorScale::create(ScaleKind::signed_log10, -1.0, 1.0, 0.0).error().code ==
          ErrorCode::invalid_input);

    const std::vector<double> n{3e5, 2e17, std::nan(""), 0.0, -4.0, 7e10};
    const ColorScale f = *ColorScale::fit(ScaleKind::log10, n);
    CHECK(f.lo() == 1e5);   // widened to whole decades
    CHECK(f.hi() == 1e18);
    const ColorScale l = *ColorScale::fit(ScaleKind::linear, n);
    CHECK(l.lo() == -4.0);
    CHECK(l.hi() == 2e17);
    const std::vector<double> flat{5.0, 5.0};
    const ColorScale c = *ColorScale::fit(ScaleKind::linear, flat);
    CHECK(c.lo() == 2.5);
    CHECK(c.hi() == 7.5);
    const ColorScale d = *ColorScale::fit(ScaleKind::log10, flat);
    CHECK(d.lo() == 1.0);  // the decade holding 5
    CHECK(d.hi() == 10.0);
    const std::vector<double> decade{100.0};
    const ColorScale p = *ColorScale::fit(ScaleKind::log10, decade);  // a power of ten: widened
    CHECK(p.lo() == 10.0);
    CHECK(p.hi() == 1000.0);
    const std::vector<double> none{std::nan(""), -1.0};
    CHECK(ColorScale::fit(ScaleKind::log10, none).error().code == ErrorCode::invalid_input);
}

TEST_CASE("color: linear ticks are 1, 2 or 5 times a power of ten") {
    const auto check = [](double lo, double hi, int target) {
        const auto t = linear_ticks(lo, hi, target);
        REQUIRE(t.size() >= 2);
        const double step = t[1].value - t[0].value;
        const double m = step / std::pow(10.0, std::floor(std::log10(step) + 1e-9));
        CAPTURE(lo, hi, step, m);
        CHECK((std::abs(m - 1.0) < 1e-9 || std::abs(m - 2.0) < 1e-9 || std::abs(m - 5.0) < 1e-9));
        for (std::size_t i = 0; i < t.size(); ++i) {
            CHECK(t[i].value >= lo - 1e-12 * std::abs(step));
            CHECK(t[i].value <= hi + 1e-12 * std::abs(step));
            CHECK(std::abs(t[i].position - (t[i].value - lo) / (hi - lo)) < 1e-12);
            if (i > 0) CHECK(std::abs(t[i].value - t[i - 1].value - step) <= 1e-9 * step);
        }
        CHECK(static_cast<int>(t.size()) >= target / 2);
        CHECK(static_cast<int>(t.size()) <= 2 * target + 1);
    };
    check(0.0, 1.0, 6);
    check(-0.7, 2.3, 6);
    check(1e-9, 7.3e-9, 5);
    check(-3e18, 9e18, 6);
    check(0.31, 0.32, 4);
    const auto t = linear_ticks(0.0, 1.0, 6);
    REQUIRE(t.size() == 6);  // 0, 0.2, ..., 1
    CHECK(t[1].label == "0.2");
    CHECK(t[0].label == "0");
    const auto neg = linear_ticks(-1.0, 1.0, 4);
    CHECK(neg.front().label == "\xE2\x88\x92" "1.0");  // U+2212, the decimals of the 0.5 step
}

TEST_CASE("color: log ticks at decades, with minor ticks on short axes") {
    const auto t = log_ticks(1e2, 1e20, 6);
    for (const Tick& k : t) {
        const double e = std::log10(k.value);
        CHECK(e == std::round(e));
        CHECK(k.major);
    }
    CHECK(t.size() <= 7);
    CHECK(t.front().value >= 1e2);
    const auto short_axis = log_ticks(1.0, 1000.0, 6);
    std::size_t majors = 0, minors = 0;
    for (const Tick& k : short_axis) (k.major ? majors : minors) += 1;
    CHECK(majors == 4);
    CHECK(minors == 24);
    CHECK(power_of_ten(18) == "10\xC2\xB9\xE2\x81\xB8");
    CHECK(power_of_ten(-3) == "10\xE2\x81\xBB\xC2\xB3");
    CHECK(power_of_ten(0) == "10\xE2\x81\xB0");
}

TEST_CASE("color: labels of large and small values use powers of ten") {
    CHECK(format_tick(2e18, 1e18, 9e18) == "2\xC3\x97" "10\xC2\xB9\xE2\x81\xB8");
    CHECK(format_tick(1e18, 1e18, 9e18) == "10\xC2\xB9\xE2\x81\xB8");
    CHECK(format_tick(2.5e-6, 5e-7, 3e-6) == "2.5\xC3\x97" "10\xE2\x81\xBB\xE2\x81\xB6");
    CHECK(format_tick(0.25, 0.05, 1.0) == "0.25");
    CHECK(format_tick(1500.0, 500.0, 2000.0) == "1500");
    const ColorScale s = *ColorScale::create(ScaleKind::signed_log10, -1e19, 1e19, 1e10);
    const auto t = s.ticks(8);
    bool zero = false;
    for (const Tick& k : t) zero = zero || (k.value == 0.0 && k.label == "0");
    CHECK(zero);
    CHECK(t.front().value < 0.0);
    CHECK(t.back().value > 0.0);
}
