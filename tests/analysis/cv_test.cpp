// C-V and interface extraction on synthetic data with exact answers (ARCHITECTURE.md section 11,
// Unit 24, gate 1): a one-sided junction's depletion capacitance, alone and behind an oxide; the
// flat-band capacitance; the conductance method on a single trap level behind an oxide.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/cv.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/results/small_signal.hpp"

using namespace NiTCAD;
using analysis::Curve;
using Complex = std::complex<double>;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

constexpr double eps_si = 11.7 * base::eps0_F_per_cm;

// A one-sided junction: C = sqrt(q eps N / (2 (V_bi - V))) per unit area, in series with C_s
// (none when zero).
Curve junction(double N, double V_bi, double C_s = 0.0) {
    std::vector<double> V, C;
    for (int k = 0; k <= 20; ++k) {
        V.push_back(-5.0 + 0.25 * k);
        const double c = std::sqrt(base::q_C * eps_si * N / (2.0 * (V_bi - V.back())));
        C.push_back(C_s > 0.0 ? 1.0 / (1.0 / c + 1.0 / C_s) : c);
    }
    return *analysis::make_curve(V, C);
}

}  // namespace

TEST_CASE("cv: quasi-static capacitance and the accumulation capacitance") {
    std::vector<double> V, Q;
    for (int k = 0; k <= 10; ++k) {
        V.push_back(-1.0 + 0.2 * k);
        Q.push_back(3e-7 * V.back() - 4e-8 * V.back() * V.back());
    }
    const auto C = analysis::quasi_static_capacitance(*analysis::make_curve(V, Q));
    REQUIRE(C.has_value());
    for (std::size_t k = 1; k + 1 < V.size(); ++k) REQUIRE(close(C->y[k], 3e-7 - 8e-8 * V[k], 1e-12));
    const auto acc = analysis::accumulation_capacitance(*C);
    REQUIRE(acc->value == C->y[0]);
    REQUIRE(acc->window.first == 0);
    const Curve negative = analysis::scaled(*C, -1.0);
    REQUIRE(!analysis::accumulation_capacitance(negative).has_value());
}

TEST_CASE("cv: the doping profile of a uniform junction, bare and behind an oxide") {
    // 1/C^2 is linear in V, so the derivative and the profile are exact.
    const double N = 3e16, V_bi = 0.8, C_ox = 6.9e-7;
    for (const double C_s : {0.0, C_ox}) {
        const auto p = analysis::doping_profile(junction(N, V_bi, C_s), eps_si, C_s);
        REQUIRE(p.has_value());
        for (std::size_t k = 0; k < p->doping_cm3.size(); ++k) {
            const double V = -5.0 + 0.25 * static_cast<double>(k);
            REQUIRE(close(p->doping_cm3[k], N, 1e-9));
            REQUIRE(close(p->depth_cm[k], std::sqrt(2.0 * eps_si * (V_bi - V) / (base::q_C * N)), 1e-12));
        }
    }
    REQUIRE(!analysis::doping_profile(junction(N, V_bi), 0.0).has_value());
    REQUIRE(!analysis::doping_profile(junction(N, V_bi), eps_si, -1.0).has_value());
    // A capacitance at or above the series capacitance cannot be in series with it.
    REQUIRE(!analysis::doping_profile(junction(N, V_bi), eps_si, 1e-9).has_value());
    const Curve flat = *analysis::make_curve({0.0, 1.0, 2.0}, {1e-7, 1e-7, 1e-7});
    REQUIRE(!analysis::doping_profile(flat, eps_si).has_value());
}

TEST_CASE("cv: the flat-band capacitance and voltage") {
    const double C_ox = 6.9e-7, N = 1e17, T = 300.0;
    const double L_D = std::sqrt(eps_si * base::thermal_voltage(T) / (base::q_C * N));
    const double C_fb = analysis::flat_band_capacitance(C_ox, eps_si, N, T);
    REQUIRE(close(C_fb, C_ox * (eps_si / L_D) / (C_ox + eps_si / L_D), 1e-14));
    // A falling high-frequency-like curve through C_fb at -0.8 V.
    std::vector<double> V, C;
    for (int k = 0; k <= 20; ++k) {
        V.push_back(-2.0 + 0.1 * k);
        C.push_back(C_fb - 1e-7 * (V.back() + 0.8));
    }
    const auto vfb = analysis::flat_band_voltage(*analysis::make_curve(V, C), C_fb);
    REQUIRE(std::abs(vfb->value + 0.8) < 1e-12);
}

TEST_CASE("cv: the conductance method finds a trap level's peak") {
    // Behind C_ox: C_d in parallel with a trap level, Y_it = C_it i omega / (1 + i omega tau), whose
    // G_p / omega = C_it omega tau / (1 + omega^2 tau^2) peaks at C_it / 2 at omega tau = 1.
    const double C_ox = 3.45e-7, C_d = 5e-8, C_it = 2e-8, tau = 1e-5;
    const auto run = [&](const std::vector<double>& wt) {
        results::SmallSignal r;
        r.contacts = 1;
        results::SmallSignalPoint p;
        for (const double x : wt) {
            const double w = x / tau;
            r.frequency_Hz.push_back(w / (2.0 * std::numbers::pi));
            const Complex ys = Complex{0.0, w * C_d} + C_it * Complex{0.0, w} / Complex{1.0, w * tau};
            p.admittance.push_back({1.0 / (1.0 / ys + 1.0 / Complex{0.0, w * C_ox})});
        }
        r.points.push_back(p);
        return r;
    };
    // Symmetric samples about the peak in ln f: exact.
    std::vector<double> wt;
    for (int k = -8; k <= 8; ++k) wt.push_back(std::pow(10.0, 0.25 * k));
    const auto peak = analysis::conductance_peak(run(wt), 0, 0, C_ox);
    REQUIRE(peak.has_value());
    REQUIRE(close(peak->conductance_over_omega.value, C_it / 2.0, 1e-10));
    REQUIRE(close(peak->frequency.value, 1.0 / (2.0 * std::numbers::pi * tau), 1e-10));
    REQUIRE(peak->frequency.window.first == 7);
    // Off-centre samples (a quarter decade apart, 0.1 decade off): the parabola recovers it within
    // its fourth-order error (measured 0.31% low in the peak, 0.41% in the frequency).
    std::vector<double> skew;
    for (int k = -8; k <= 8; ++k) skew.push_back(std::pow(10.0, 0.25 * k + 0.1));
    const auto off = analysis::conductance_peak(run(skew), 0, 0, C_ox);
    CAPTURE(off->conductance_over_omega.value / (C_it / 2.0) - 1.0,
            off->frequency.value * 2.0 * std::numbers::pi * tau - 1.0);
    REQUIRE(close(off->conductance_over_omega.value, C_it / 2.0, 5e-3));
    REQUIRE(close(off->frequency.value, 1.0 / (2.0 * std::numbers::pi * tau), 1e-2));
    // The peak at the edge of the range, bad input.
    REQUIRE(!analysis::conductance_peak(run({1.0, 3.0, 10.0}), 0, 0, C_ox).has_value());
    REQUIRE(!analysis::conductance_peak(run({0.3, 1.0}), 0, 0, C_ox).has_value());
    REQUIRE(!analysis::conductance_peak(run(wt), 1, 0, C_ox).has_value());
    REQUIRE(!analysis::conductance_peak(run(wt), 0, 1, C_ox).has_value());
    REQUIRE(!analysis::conductance_peak(run(wt), 0, 0, 0.0).has_value());
    REQUIRE(close(analysis::interface_trap_density(1e-8), 2.5e-8 / base::q_C, 1e-15));
}
