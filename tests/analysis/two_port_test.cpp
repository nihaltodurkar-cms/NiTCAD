// Two-port parameters and RF figures on networks with known answers (ARCHITECTURE.md section 11,
// Unit 24, gate 1): the conversions round-trip, S of matched and open ports, Mason's U on a
// unilateral amplifier and its invariance under lossless reciprocal embedding, and f_T and f_max
// of the intrinsic transistor's pi model (g_m, C_gs, C_gd, g_ds) with a gate resistance.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <vector>

#include "NiTCAD/analysis/two_port.hpp"
#include "NiTCAD/results/small_signal.hpp"

using namespace NiTCAD;
using analysis::Complex;
using analysis::TwoPort;
using base::ErrorCode;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

double distance(const TwoPort& a, const TwoPort& b) {
    return std::abs(a.p11 - b.p11) + std::abs(a.p12 - b.p12) + std::abs(a.p21 - b.p21) +
           std::abs(a.p22 - b.p22);
}

double size(const TwoPort& a) {
    return std::abs(a.p11) + std::abs(a.p12) + std::abs(a.p21) + std::abs(a.p22);
}

constexpr double two_pi = 2.0 * std::numbers::pi;

struct PiModel {
    double gm = 5e-3, gds = 2e-4, Cgs = 4e-14, Cgd = 1e-14, Cds = 2e-14, rg = 0.0;

    // Intrinsic Y, then the gate resistance in series with port 1 (z11 += r_g).
    [[nodiscard]] TwoPort y(double f) const {
        const double w = two_pi * f;
        const Complex j{0.0, 1.0};
        const TwoPort yi{j * w * (Cgs + Cgd), -j * w * Cgd, gm - j * w * Cgd, gds + j * w * (Cgd + Cds)};
        if (rg == 0.0) return yi;
        TwoPort z = *analysis::inverse(yi);
        z.p11 += rg;
        return *analysis::inverse(z);
    }
};

// A small-signal run of a three-contact device (gate 0, drain 1, source 2) at one operating point:
// the source row and column close the admittance so every row and column sums to zero.
results::SmallSignal run(const PiModel& m, const std::vector<double>& f) {
    results::SmallSignal r;
    r.contacts = 3;
    r.frequency_Hz = f;
    results::SmallSignalPoint p;
    for (const double v : f) {
        const TwoPort y = m.y(v);
        const Complex a[2][2] = {{y.p11, y.p12}, {y.p21, y.p22}};
        std::vector<Complex> Y(9);
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                Y[i * 3 + j] = a[i][j];
                Y[i * 3 + 2] -= a[i][j];
                Y[2 * 3 + j] -= a[i][j];
                Y[2 * 3 + 2] += a[i][j];
            }
        }
        p.admittance.push_back(Y);
    }
    r.points.push_back(p);
    return r;
}

std::vector<double> log_grid(double from, double to, int per_decade) {
    std::vector<double> f;
    const int n = static_cast<int>(std::round(std::log10(to / from) * per_decade));
    for (int k = 0; k <= n; ++k) f.push_back(from * std::pow(10.0, static_cast<double>(k) / per_decade));
    return f;
}

}  // namespace

TEST_CASE("two-port: reading the ports from the admittance matrix") {
    const PiModel m;
    const auto r = run(m, {1e9});
    const auto y = analysis::admittance_two_port(r, 0, 0, 0, 1);
    REQUIRE(distance(*y, m.y(1e9)) == 0.0);
    const auto swapped = analysis::admittance_two_port(r, 0, 0, 1, 0);
    REQUIRE(swapped->p11 == m.y(1e9).p22);
    REQUIRE(swapped->p21 == m.y(1e9).p12);
    REQUIRE(!analysis::admittance_two_port(r, 0, 0, 1, 1).has_value());
    REQUIRE(!analysis::admittance_two_port(r, 1, 0, 0, 1).has_value());
    REQUIRE(!analysis::admittance_two_port(r, 0, 1, 0, 1).has_value());
    REQUIRE(!analysis::admittance_two_port(r, 0, 0, 0, 3).has_value());
}

TEST_CASE("two-port: Y, Z, h and S round-trip") {
    const TwoPort y{Complex{1e-3, 2e-3}, Complex{-1e-5, -3e-4}, Complex{4e-2, -3e-4}, Complex{2e-4, 6e-4}};
    const TwoPort z = *analysis::inverse(y);
    REQUIRE(distance(*analysis::inverse(z), y) < 1e-15 * size(y));
    const TwoPort h = *analysis::y_to_h(y);
    REQUIRE(h.p21 == y.p21 / y.p11);
    REQUIRE(distance(*analysis::h_to_y(h), y) < 1e-15 * size(y));
    const TwoPort s = *analysis::y_to_s(y, 50.0);
    REQUIRE(distance(*analysis::s_to_y(s, 50.0), y) < 1e-14 * size(y));
    // Matched ports reflect nothing; open ports everything.
    const TwoPort matched{Complex{0.02, 0.0}, 0.0, 0.0, Complex{0.02, 0.0}};
    REQUIRE(size(*analysis::y_to_s(matched, 50.0)) < 1e-16);
    const TwoPort open{0.0, 0.0, 0.0, 0.0};
    REQUIRE(distance(*analysis::y_to_s(open, 50.0), TwoPort{1.0, 0.0, 0.0, 1.0}) == 0.0);
    // S of a series resistor R between the ports: S11 = R / (R + 2 z0), S21 = 2 z0 / (R + 2 z0).
    const double R = 30.0, z0 = 50.0;
    const TwoPort series{1.0 / R, -1.0 / R, -1.0 / R, 1.0 / R};
    const TwoPort ss = *analysis::y_to_s(series, z0);
    REQUIRE(close(ss.p11.real(), R / (R + 2 * z0), 1e-14));
    REQUIRE(close(ss.p21.real(), 2 * z0 / (R + 2 * z0), 1e-14));
    REQUIRE(analysis::inverse(open).error().code == ErrorCode::singular_system);
    REQUIRE(analysis::y_to_h(open).error().code == ErrorCode::singular_system);
    REQUIRE(analysis::y_to_s(y, 0.0).error().code == ErrorCode::invalid_input);
    REQUIRE(distance(analysis::scaled(y, 2.0), TwoPort{2.0 * y.p11, 2.0 * y.p12, 2.0 * y.p21, 2.0 * y.p22}) == 0.0);
}

TEST_CASE("two-port: Mason's U, the stability factor and the maximum gain") {
    // A unilateral amplifier with real input and output conductances: U = g_m^2 / (4 g_1 g_2).
    const TwoPort u{1e-3, 0.0, 4e-2, 2e-4};
    REQUIRE(close(*analysis::unilateral_gain(u), 4e-2 * 4e-2 / (4.0 * 1e-3 * 2e-4), 1e-14));
    REQUIRE(std::isinf(analysis::stability_factor(u)));
    REQUIRE(!analysis::maximum_gain(u).has_value());
    // U is invariant under lossless reciprocal embedding: a shunt susceptance at the input and a
    // series reactance at the output.
    PiModel m;
    m.rg = 15.0;
    const TwoPort y = m.y(2e9);
    TwoPort embedded = y;
    embedded.p11 += Complex{0.0, 3e-3};
    TwoPort z = *analysis::inverse(embedded);
    z.p22 += Complex{0.0, 40.0};
    embedded = *analysis::inverse(z);
    REQUIRE(close(*analysis::unilateral_gain(embedded), *analysis::unilateral_gain(y), 1e-10));
    // k and MAG by hand: a passive attenuator (pi of conductances) is unconditionally stable and has
    // gain below 1.
    const double a = 0.01, b = 0.02;
    const TwoPort att{a + b, -b, -b, a + b};
    const double k = analysis::stability_factor(att);
    REQUIRE(close(k, (2.0 * (a + b) * (a + b) - b * b) / (b * b), 1e-14));
    const double mag = *analysis::maximum_gain(att);
    REQUIRE(close(mag, k - std::sqrt(k * k - 1.0), 1e-14));
    REQUIRE(mag < 1.0);
    // Below k = 1 the maximum stable gain |y21 / y12|: the unloaded pi model at low frequency.
    const TwoPort low = PiModel{}.y(1e6);
    REQUIRE(analysis::stability_factor(low) < 1.0);
    REQUIRE(close(*analysis::maximum_gain(low), std::abs(low.p21 / low.p12), 1e-14));
    REQUIRE(!analysis::unilateral_gain(TwoPort{0.0, 0.0, 1.0, 0.0}).has_value());
}

TEST_CASE("two-port: the unity-gain frequency, interpolated and extrapolated") {
    // gain = f_T / f exactly: the log-log interpolation is exact.
    const double fT = 7.3e9;
    std::vector<double> f, g;
    for (const double v : log_grid(1e8, 1e11, 5)) {
        f.push_back(v);
        g.push_back(fT / v);
    }
    const auto e = analysis::unity_gain_frequency(f, g);
    REQUIRE(close(e->value, fT, 1e-13));
    REQUIRE(!e->extrapolated);
    // Data ending below f_T: extrapolated at -20 dB/decade, exact here too.
    const std::vector<double> f_low(f.begin(), f.begin() + 8), g_low(g.begin(), g.begin() + 8);
    const auto x = analysis::unity_gain_frequency(f_low, g_low);
    REQUIRE(close(x->value, fT, 1e-13));
    REQUIRE(x->extrapolated);
    // A zero frequency or an undefined gain is skipped.
    std::vector<double> f0{0.0}, g0{std::nan("")};
    f0.insert(f0.end(), f.begin(), f.end());
    g0.insert(g0.end(), g.begin(), g.end());
    REQUIRE(close(analysis::unity_gain_frequency(f0, g0)->value, fT, 1e-13));
    REQUIRE(analysis::unity_gain_frequency(f0, g0)->window.first >= 1);
    // Errors.
    const std::vector<double> high(f.end() - 3, f.end()), gh(g.end() - 3, g.end());
    REQUIRE(!analysis::unity_gain_frequency(high, gh).has_value());
    REQUIRE(!analysis::unity_gain_frequency(f, g_low).has_value());
    const std::vector<double> back{2e9, 1e9}, gb{3.0, 2.0};
    REQUIRE(!analysis::unity_gain_frequency(back, gb).has_value());
}

TEST_CASE("two-port: f_T and f_max of the pi model") {
    // |h21|^2 = (g_m^2 + w^2 C_gd^2) / (w^2 (C_gs + C_gd)^2) = 1 at
    // w_T = g_m / sqrt((C_gs + C_gd)^2 - C_gd^2), with or without the gate resistance: h21 = -z21 / z22,
    // and a series input resistance changes only z11.
    PiModel m;
    m.rg = 20.0;
    const double wT = m.gm / std::sqrt((m.Cgs + m.Cgd) * (m.Cgs + m.Cgd) - m.Cgd * m.Cgd);
    const auto r = run(m, log_grid(1e8, 1e12, 20));
    const auto fT = analysis::transition_frequency(r, 0, 0, 1);
    REQUIRE(fT.has_value());
    CAPTURE(fT->value, wT / two_pi);
    REQUIRE(close(fT->value, wT / two_pi, 1e-4));
    REQUIRE(!fT->extrapolated);
    // f_max where U = 1, found by bisection on U itself (U verified above) against the extractor's
    // interpolation of sqrt(U) on the grid.
    const auto U = [&](double f) { return *analysis::unilateral_gain(m.y(f)); };
    double lo = 1e9, hi = 1e12;
    for (int k = 0; k < 200; ++k) {
        const double mid = std::sqrt(lo * hi);
        (U(mid) > 1.0 ? lo : hi) = mid;
    }
    const auto fmax = analysis::maximum_oscillation_frequency(r, 0, 0, 1);
    REQUIRE(fmax.has_value());
    CAPTURE(fmax->value, lo);
    REQUIRE(close(fmax->value, lo, 1e-3));
    // The textbook estimate f_max ~ f_T / (2 sqrt(r_g (g_ds + w_T C_gd))) is within 30%.
    const double estimate = (wT / two_pi) / (2.0 * std::sqrt(m.rg * (m.gds + wT * m.Cgd)));
    CAPTURE(estimate);
    REQUIRE(close(fmax->value, estimate, 0.3));
}
