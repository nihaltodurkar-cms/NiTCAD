// Curves read from results and the shared tools (ARCHITECTURE.md section 11, Unit 24): reading
// sweeps and small-signal runs, input checks, the numpy.gradient derivative, interpolation and
// crossings, exact on the functions they are built for.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/small_signal.hpp"

using namespace NiTCAD;
using analysis::Curve;
using analysis::Scale;
using base::ErrorCode;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

results::BiasPoint point(std::vector<double> bias, std::vector<double> current,
                         std::vector<double> resolution, std::vector<double> gate_charge = {}) {
    results::BiasPoint p;
    p.bias_V = std::move(bias);
    p.terminal_current = std::move(current);
    p.terminal_current_resolution = std::move(resolution);
    p.gate_charge = std::move(gate_charge);
    return p;
}

}  // namespace

TEST_CASE("curve: make_curve checks its input") {
    REQUIRE(analysis::make_curve({0.0, 1.0}, {1.0, 2.0}).has_value());
    REQUIRE(analysis::make_curve({0.0}, {1.0}).error().code == ErrorCode::invalid_input);
    REQUIRE(analysis::make_curve({0.0, 1.0}, {1.0}).error().code == ErrorCode::invalid_input);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto bad = analysis::make_curve({0.0, 1.0, 2.0}, {1.0, nan, 2.0});
    REQUIRE(bad.error().code == ErrorCode::invalid_input);
    REQUIRE(bad.error().context->index == 1u);
    REQUIRE(!analysis::make_curve({0.0, 1.0}, {1.0, 2.0}, {1.0}).has_value());
    REQUIRE(!analysis::make_curve({0.0, 1.0}, {1.0, 2.0}, {1.0, -1.0}).has_value());
}

TEST_CASE("curve: reading a sweep's current, resolution and gate charge") {
    results::Sweep s;
    s.points.push_back(point({0.0, 0.1}, {1.0, -1.0}, {1e-3, 2e-3}, {0.0, 5.0}));
    s.points.push_back(point({0.0, 0.2}, {2.0, -2.0}, {1e-3, 3e-3}, {0.0, 6.0}));
    const auto c = analysis::current_curve(s, 1, 0);
    REQUIRE(c.has_value());
    REQUIRE(c->x == std::vector<double>{0.1, 0.2});
    REQUIRE(c->y == std::vector<double>{1.0, 2.0});
    REQUIRE(c->resolution == std::vector<double>{1e-3, 1e-3});
    const auto q = analysis::gate_charge_curve(s, 1, 1);
    REQUIRE(q->y == std::vector<double>{5.0, 6.0});
    REQUIRE(q->resolution.empty());
    // A point without a resolution (a quasi-static sweep) leaves the whole curve without one.
    s.points[1].terminal_current_resolution.clear();
    REQUIRE(analysis::current_curve(s, 1, 0)->resolution.empty());
    REQUIRE(!analysis::current_curve(s, 2, 0).has_value());
    REQUIRE(!analysis::current_curve(s, 1, 2).has_value());
    REQUIRE(!analysis::gate_charge_curve(s, 1, 2).has_value());
}

TEST_CASE("curve: capacitance from a small-signal run") {
    results::SmallSignal r;
    r.contacts = 2;
    r.frequency_Hz = {0.0, 1e6};
    const double w = 2.0 * std::numbers::pi * 1e6;
    for (const double v : {-1.0, 0.0, 1.0}) {
        results::SmallSignalPoint p;
        p.dc.bias_V = {v, 0.0};
        const std::complex<double> y{1e-3, w * (2e-7 + 1e-8 * v)};
        p.admittance = {std::vector<std::complex<double>>(4), {y, -y, -y, y}};
        r.points.push_back(p);
    }
    const auto c = analysis::capacitance_curve(r, 1, 0, 0, 0);
    REQUIRE(c.has_value());
    REQUIRE(c->x == std::vector<double>{-1.0, 0.0, 1.0});
    for (std::size_t k = 0; k < 3; ++k) REQUIRE(close(c->y[k], 2e-7 + 1e-8 * c->x[k], 1e-14));
    REQUIRE(!analysis::capacitance_curve(r, 0, 0, 0, 0).has_value());  // zero frequency
    REQUIRE(!analysis::capacitance_curve(r, 2, 0, 0, 0).has_value());
    REQUIRE(!analysis::capacitance_curve(r, 1, 2, 0, 0).has_value());
}

TEST_CASE("curve: scaled and slice") {
    const Curve c = *analysis::make_curve({0.0, 1.0, 2.0, 3.0}, {1.0, 2.0, 3.0, 4.0},
                                          {0.1, 0.1, 0.1, 0.1});
    const Curve s = analysis::scaled(c, -2.0);
    REQUIRE(s.y == std::vector<double>{-2.0, -4.0, -6.0, -8.0});
    REQUIRE(s.resolution == std::vector<double>{0.2, 0.2, 0.2, 0.2});
    const auto part = analysis::slice(c, 2.5, 0.5);
    REQUIRE(part->x == std::vector<double>{1.0, 2.0});
    REQUIRE(part->resolution.size() == 2);
    REQUIRE(!analysis::slice(c, 0.5, 1.5).has_value());
}

TEST_CASE("curve: the derivative is numpy.gradient, exact for a quadratic inside") {
    // Uneven, decreasing x: the three-point formula is exact for a quadratic at interior points;
    // the ends are one-sided differences.
    const std::vector<double> x{3.0, 2.5, 1.75, 1.0, 0.2, -0.5};
    std::vector<double> y;
    for (const double v : x) y.push_back(2.0 * v * v - 3.0 * v + 1.0);
    const auto d = analysis::derivative(*analysis::make_curve(x, y));
    REQUIRE(d.has_value());
    for (std::size_t k = 1; k + 1 < x.size(); ++k) REQUIRE(close((*d)[k], 4.0 * x[k] - 3.0, 1e-13));
    REQUIRE(close((*d)[0], (y[1] - y[0]) / (x[1] - x[0]), 1e-15));
    REQUIRE(close((*d)[5], (y[5] - y[4]) / (x[5] - x[4]), 1e-15));
    const auto bad = analysis::derivative(*analysis::make_curve({0.0, 1.0, 0.5}, {0.0, 1.0, 2.0}));
    REQUIRE(bad.error().code == ErrorCode::invalid_input);
    REQUIRE(!analysis::strictly_monotone(std::vector<double>{0.0, 1.0, 1.0}));
    REQUIRE(analysis::strictly_monotone(std::vector<double>{1.0, 0.0, -1.0}));
}

TEST_CASE("curve: value_at, linear and logarithmic") {
    std::vector<double> x, y;
    for (int k = 0; k <= 10; ++k) {
        x.push_back(0.1 * k);
        y.push_back(1e-9 * std::exp(x.back() / 0.04));
    }
    const Curve c = *analysis::make_curve(x, y);
    const auto e = analysis::value_at(c, 0.537, Scale::logarithmic);
    REQUIRE(close(e->value, 1e-9 * std::exp(0.537 / 0.04), 1e-12));
    REQUIRE(e->window.first == 5);
    REQUIRE(e->window.last == 6);
    const Curve line = *analysis::make_curve({0.0, 1.0, 2.0}, {1.0, 3.0, 4.0});
    REQUIRE(analysis::value_at(line, 0.25)->value == 1.5);
    REQUIRE(analysis::value_at(line, 2.0)->value == 4.0);
    REQUIRE(analysis::value_at(line, 1.0)->window.first == 1);  // on a point: that point alone
    REQUIRE(analysis::value_at(line, 1.0)->window.last == 1);
    REQUIRE(!analysis::value_at(line, 2.1).has_value());
    const Curve across = *analysis::make_curve({0.0, 1.0}, {-1.0, 1.0});
    REQUIRE(!analysis::value_at(across, 0.5, Scale::logarithmic).has_value());
}

TEST_CASE("curve: crossings from either side, the resolution flag") {
    const Curve c =
        *analysis::make_curve({0.0, -1.0, -2.0, -3.0}, {1.0, 3.0, 7.0, 2.0}, {0.0, 0.0, 0.0, 9.0});
    const auto up = analysis::crossing(c, 5.0);
    REQUIRE(up->value == -1.5);
    REQUIRE(!up->below_resolution);
    // The first crossing in the curve's order: the falling one after the rise is not reached.
    const auto fall = analysis::crossing(c, 2.0);
    REQUIRE(fall->value == -0.5);
    const auto late = analysis::crossing(*analysis::make_curve({0.0, 1.0, 2.0}, {7.0, 7.0, 2.0},
                                                               {0.0, 0.0, 9.0}),
                                         4.0);
    REQUIRE(close(late->value, 1.6, 1e-15));
    REQUIRE(late->below_resolution);
    REQUIRE(analysis::crossing(c, 1.0)->value == 0.0);
    REQUIRE(analysis::crossing(c, 1.0)->window.last == 0);
    REQUIRE(analysis::crossing(c, 8.0).error().code == ErrorCode::invalid_input);
    // Logarithmic: exact for an exponential.
    const Curve e = *analysis::make_curve({0.0, 1.0}, {1e-6, 1e-2});
    REQUIRE(close(analysis::crossing(e, 1e-4, Scale::logarithmic)->value, 0.5, 1e-14));
    REQUIRE(!analysis::crossing(e, -1e-4, Scale::logarithmic).has_value());
}
