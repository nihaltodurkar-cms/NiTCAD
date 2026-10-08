// DC extraction on synthetic curves whose figures are known exactly (ARCHITECTURE.md section 11,
// Unit 24, gate 1): a piecewise-linear transfer curve (max-g_m threshold), exponential
// subthreshold currents (constant-current threshold, swing, DIBL, on/off), a linear saturation
// tail (g_ds), the diode law with and without series resistance, and a parabolic fold.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/base/constants.hpp"

using namespace NiTCAD;
using analysis::Curve;
using base::ErrorCode;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> v(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) v[static_cast<std::size_t>(k)] = a + (b - a) * k / (n - 1);
    return v;
}

template <class F>
Curve sample(const std::vector<double>& x, F f) {
    std::vector<double> y;
    for (const double v : x) y.push_back(f(v));
    return *analysis::make_curve(x, std::move(y));
}

// I = I0 10^((V - V0) / S) with S in V/decade.
Curve subthreshold(double S, double V0, double shift = 0.0) {
    return sample(linspace(-0.5, 0.5, 21),
                  [&](double v) { return 1e-7 * std::pow(10.0, (v + shift - V0) / S); });
}

}  // namespace

TEST_CASE("dc: max-g_m threshold of a piecewise-linear curve is exact") {
    // I = 0 below V_t, g (V - V_t) above, V_t on a grid point: the largest numpy.gradient is g, first
    // reached one point above V_t, and its tangent meets zero at V_t.
    const double Vt = 0.4, g = 2e-4, Vds = 0.05;
    const Curve c = sample(linspace(0.0, 1.0, 26), [&](double v) { return v > Vt ? g * (v - Vt) : 0.0; });
    const auto gm = analysis::max_transconductance(c);
    REQUIRE(close(gm->value, g, 1e-12));
    const auto vt = analysis::threshold_max_gm(c, Vds);
    REQUIRE(vt.has_value());
    REQUIRE(std::abs(vt->value - (Vt - Vds / 2.0)) < 1e-12);
    REQUIRE(vt->extrapolated);
    const Curve flat = sample(linspace(0.0, 1.0, 5), [](double v) { return 1.0 - v; });
    REQUIRE(analysis::threshold_max_gm(flat, 0.0).error().code == ErrorCode::invalid_input);
}

TEST_CASE("dc: constant-current threshold, swing, DIBL and on/off on exponential currents") {
    const double S = 0.075, V0 = 0.1;
    const Curve c = subthreshold(S, V0);
    // 1e-7 is reached exactly at V0, between samples.
    const auto vt = analysis::threshold_constant_current(c, 1e-7);
    REQUIRE(std::abs(vt->value - V0) < 1e-13);
    REQUIRE(vt->window.last == vt->window.first + 1);
    REQUIRE(!analysis::threshold_constant_current(c, 1.0).has_value());
    REQUIRE(!analysis::threshold_constant_current(c, 0.0).has_value());
    const auto ss = analysis::subthreshold_swing(c, 1.0);
    REQUIRE(close(ss->value, 1000.0 * S, 1e-12));
    REQUIRE(ss->unit == "mV/decade");
    // DIBL: a drain bias of 1 V lowers the threshold by eta (1 V - 0.05 V).
    const double eta = 0.08;
    const Curve high = subthreshold(S, V0, eta * 0.95);
    const auto d = analysis::dibl(c, high, 0.05, 1.0, 1e-7);
    REQUIRE(close(d->value, eta, 1e-10));
    REQUIRE(!analysis::dibl(c, high, 1.0, 1.0, 1e-7).has_value());
    const auto r = analysis::on_off_ratio(c, 0.43, -0.37);
    REQUIRE(close(r->value, std::pow(10.0, 0.8 / S), 1e-11));
    REQUIRE(!analysis::on_off_ratio(c, 0.6, 0.0).has_value());
}

TEST_CASE("dc: the swing is the steepest segment below the ceiling and above the resolution") {
    // 90 mV/decade at low current, a steeper 60 mV/decade segment through 1e-9 to 1e-8 (the
    // device's swing), 80 mV/decade to 1e-7, and above 1e-2 I_max a far steeper (unphysical) rise
    // the ceiling excludes; the lowest points carry a resolution above their current, and their
    // (steep) segment is noise.
    const std::vector<double> x{0.0, 0.01, 0.1, 0.19, 0.25, 0.33, 0.34, 0.35};
    const std::vector<double> y{1e-13, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-4, 1e-3};
    std::vector<double> res(x.size(), 0.0);
    res[0] = res[1] = 2e-11;
    const Curve c = *analysis::make_curve(x, y, res);
    const auto ss = analysis::subthreshold_swing(c);
    REQUIRE(close(ss->value, 60.0, 1e-12));
    REQUIRE(ss->window.first == 3);
    const auto all = analysis::subthreshold_swing(c, 1.0);
    REQUIRE(close(all->value, 10.0 / 3.0, 1e-12));
    REQUIRE(!analysis::subthreshold_swing(c, 0.0).has_value());
    const Curve falling = sample(linspace(0.0, 1.0, 5), [](double v) { return std::exp(-v); });
    REQUIRE(!analysis::subthreshold_swing(falling, 1.0).has_value());
}

TEST_CASE("dc: output conductance of a linear saturation tail") {
    const double g = 3e-6;
    const Curve c = sample(linspace(0.0, 2.0, 41), [&](double v) {
        return v < 1.0 ? 1e-3 * (2.0 * v - v * v) : 1e-3 + g * (v - 1.0);
    });
    const auto gds = analysis::output_conductance(c, 0.25);
    REQUIRE(close(gds->value, g, 1e-9));
    REQUIRE(gds->window.first == 30);
    REQUIRE(gds->window.last == 40);
    REQUIRE(!analysis::output_conductance(c, 0.01).has_value());  // fewer than 3 points
    REQUIRE(!analysis::output_conductance(c, 1.5).has_value());
}

TEST_CASE("dc: the diode law, its ideality, saturation current and series resistance") {
    const double T = 300.0, V_T = base::thermal_voltage(T), n = 1.37, Is = 3e-14;
    const Curve c = sample(linspace(0.0, 0.7, 36), [&](double v) { return Is * std::exp(v / (n * V_T)); });
    const auto fit = analysis::diode_fit(c, 0.2, 0.6, T);
    REQUIRE(close(fit->ideality.value, n, 1e-12));
    REQUIRE(close(fit->saturation_current.value, Is, 1e-10));
    REQUIRE(fit->saturation_current.extrapolated);
    const auto local = analysis::local_ideality(c, T);
    for (const double v : *local) REQUIRE(close(v, n, 1e-12));
    REQUIRE(!analysis::diode_fit(c, 0.2, 0.6, 0.0).has_value());
    REQUIRE(!analysis::diode_fit(c, 0.71, 0.9, T).has_value());
    const Curve reverse = sample(linspace(-1.0, 0.0, 5), [](double v) { return -1e-12 * (1.0 - v); });
    REQUIRE(!analysis::diode_fit(reverse, -1.0, 0.0, T).has_value());

    // V = n V_T ln(I / I_s) + I R_s, sampled on log-spaced currents.
    const double Rs = 2.5;
    std::vector<double> V, I;
    for (int k = 0; k <= 160; ++k) {
        I.push_back(1e-6 * std::pow(10.0, k / 40.0));
        V.push_back(n * V_T * std::log(I.back() / Is) + I.back() * Rs);
    }
    const Curve r = *analysis::make_curve(V, I);
    const auto sr = analysis::series_resistance(r, V.front(), V.back(), T);
    REQUIRE(sr.has_value());
    CAPTURE(sr->ideality.value, sr->resistance.value);
    REQUIRE(close(sr->ideality.value, n, 1e-3));
    REQUIRE(close(sr->resistance.value, Rs, 1e-3));
}

TEST_CASE("dc: breakdown at a current limit and at a turning point") {
    // A reverse curve: the bias falls, the current grows negative exponentially.
    const Curve c = sample(linspace(0.0, -40.0, 41), [](double v) { return -1e-12 * std::exp(-v / 2.0); });
    const auto bd = analysis::breakdown_at_current(c, 1e-6);
    REQUIRE(std::abs(bd->value - (-2.0 * std::log(1e6))) < 1e-12);
    REQUIRE(!analysis::breakdown_at_current(c, 1.0).has_value());
    REQUIRE(!analysis::breakdown_at_current(c, -1.0).has_value());

    // A snapback: along the trace ln|I| rises steadily while the bias reaches -V_b at ln|I| = t_b
    // and turns back: V = -V_b + c (t - t_b)^2, sampled unevenly. The parabola is exact.
    const double Vb = 52.3, tb = -12.0;
    std::vector<double> x, y;
    for (const double t : {-30.0, -22.0, -16.0, -13.1, -11.6, -9.0, -6.5}) {
        x.push_back(-Vb + 0.05 * (t - tb) * (t - tb));
        y.push_back(-std::exp(t));
    }
    const Curve fold = *analysis::make_curve(x, y);
    const auto tp = analysis::turning_point(fold);
    REQUIRE(tp.has_value());
    REQUIRE(std::abs(tp->value + Vb) < 1e-11);
    REQUIRE(tp->window.first == 3);
    REQUIRE(tp->window.last == 5);
    REQUIRE(analysis::turning_point(c).error().code == ErrorCode::invalid_input);
}
