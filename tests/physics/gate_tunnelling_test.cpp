// Gate tunnelling functions (Unit 51, DECISIONS.md GC2-GC5): the WKB segment action against
// quadrature and its closed forms, the Fowler-Nordheim and direct-tunnelling limits, the supply and
// equilibrium, the energy quadrature against a fine reference, the partials against central
// differences, and the insulator band data's validation.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/gate_tunnelling.hpp"
#include "NiTCAD/physics/insulator.hpp"

using namespace NiTCAD::physics;
namespace base = NiTCAD::base;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

constexpr double nm = 1e-7;  // [cm]
constexpr double inf = std::numeric_limits<double>::infinity();

// Composite Gauss-Legendre (5 points) on [a, b] with m panels.
template <class F>
double gauss(F f, double a, double b, int m) {
    static constexpr double x[5] = {-0.9061798459386640, -0.5384693101056831, 0.0,
                                    0.5384693101056831, 0.9061798459386640};
    static constexpr double w[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                    0.4786286704993665, 0.2369268850561891};
    const double h = (b - a) / m;
    double s = 0.0;
    for (int p = 0; p < m; ++p)
        for (int j = 0; j < 5; ++j) s += w[j] * f(a + h * (p + 0.5 * (1.0 + x[j])));
    return 0.5 * h * s;
}

// S by quadrature of the defining integrand, the turning points located exactly per segment. At a
// turning point x_t the integrand goes as sqrt(|x - x_t|); the map x = x_t + (x_e - x_t) s^2 makes
// it smooth in s (dx = 2 (x_e - x_t) s ds).
double reference_exponent(const std::vector<BarrierSegment>& barrier, double E) {
    double S = 0.0;
    for (const BarrierSegment& s : barrier) {
        const auto f = [&](double x) {
            const double W = s.W0_eV + (s.W1_eV - s.W0_eV) * x / s.length_cm;
            return std::sqrt(std::max(W - E, 0.0));
        };
        const double a0 = s.W0_eV - E, b0 = s.W1_eV - E;
        if (a0 <= 0.0 && b0 <= 0.0) continue;
        double I = 0.0;
        if (a0 >= 0.0 && b0 >= 0.0) {
            I = gauss(f, 0.0, s.length_cm, 2000);
        } else {
            const double xt = s.length_cm * a0 / (a0 - b0);
            const double xe = a0 > 0.0 ? 0.0 : s.length_cm;  // the end inside the barrier
            I = gauss([&](double u) { return f(xt + (xe - xt) * u * u) * 2.0 * std::abs(xe - xt) * u; },
                      0.0, 1.0, 2000);
        }
        S += wkb_action_scale() * std::sqrt(s.mass_m0) * I;
    }
    return S;
}

}  // namespace

TEST_CASE("gate tunnelling: GC2 the WKB action against quadrature and closed forms") {
    const std::vector<std::vector<BarrierSegment>> barriers{
        {{3.0, 2.0, 1.5 * nm, 0.5}},                                       // trapezoid
        {{3.0, -1.0, 2.0 * nm, 0.5}},                                      // turning point inside
        {{-0.5, 2.5, 2.0 * nm, 0.5}},                                      // turning point at the start side
        {{3.0, 2.9, 0.5 * nm, 0.5}, {1.5, 1.0, 1.5 * nm, 0.2}, {0.4, -0.3, 1.0 * nm, 0.4}},
    };
    for (const auto& b : barriers) {
        for (double E : {-0.2, 0.0, 0.3, 0.8, 1.2}) {
            const double S = wkb_exponent(b, E);
            const double R = reference_exponent(b, E);
            if (R == 0.0) {
                REQUIRE(S == 0.0);
            } else {
                REQUIRE(close(S, R, 1e-10));
            }
        }
    }
    // A rectangular barrier: exp(-2 sqrt(2 m (W - E)) d / hbar).
    const std::vector<BarrierSegment> rect{{3.1, 3.1, 1.2 * nm, 0.5}};
    const double kappa = std::sqrt(2.0 * 0.5 * base::m0_kg * (3.1 - 0.4) * base::q_C) / base::hbar_J_s;
    REQUIRE(close(wkb_transmission(rect, 0.4), std::exp(-2.0 * kappa * 1.2e-9), 1e-13));
    // Two identical layers equal one of the summed thickness.
    const std::vector<BarrierSegment> one{{3.0, 2.0, 2.0 * nm, 0.5}};
    const std::vector<BarrierSegment> two{{3.0, 2.5, 1.0 * nm, 0.5}, {2.5, 2.0, 1.0 * nm, 0.5}};
    for (double E : {0.0, 1.0, 2.2}) REQUIRE(close(wkb_exponent(two, E), wkb_exponent(one, E), 1e-13));
}

TEST_CASE("gate tunnelling: the segment action's partials and equal ends") {
    const double L = 1.0 * nm;
    for (auto [a, b] : {std::pair{3.0, 2.0}, {3.0, -1.0}, {-0.5, 2.5}, {2.0, 2.0 + 1e-7},
                        {2.0, 2.0 + 3e-5}, {1.0, 1.0}, {0.3, -0.2}}) {
        const SegmentAction g = segment_action(a, b, L);
        const double h = 1e-5;
        const double da = (segment_action(a + h, b, L).value - segment_action(a - h, b, L).value) / (2 * h);
        const double db = (segment_action(a, b + h, L).value - segment_action(a, b - h, L).value) / (2 * h);
        REQUIRE(close(g.d_a, da, 1e-6));
        REQUIRE(close(g.d_b, db, 1e-6));
    }
    // Equal and nearly equal ends: L sqrt(a), and the partials 1 / (4 sqrt a) each.
    REQUIRE(close(segment_action(2.0, 2.0, L).value, L * std::sqrt(2.0), 1e-15));
    const SegmentAction near = segment_action(2.0, 2.0 * (1.0 - 1e-12), L);
    REQUIRE(close(near.value, L * std::sqrt(2.0) * (1.0 - 0.25e-12), 1e-15));
    REQUIRE(close(near.d_a, L / (4.0 * std::sqrt(2.0)), 1e-11));
}

TEST_CASE("gate tunnelling: GC3 the Fowler-Nordheim slope") {
    // A metal-like source (F1 = 0) under a triangular SiO2 barrier of height 3.0 eV, m = 0.5.
    const double phi = 3.0, m = 0.5, T = 300.0;
    std::vector<double> x, y;
    for (double F_MV = 10.0; F_MV <= 14.0001; F_MV += 0.5) {
        const double F = F_MV * 1e6;  // [V/cm]
        const double t = 10.0 * nm;   // long enough that the barrier is triangular at every energy
        const std::vector<BarrierSegment> b{{phi, phi - F * t, t, m}};
        const double J = tsu_esaki_current(b, {-2.0, inf}, 0.0, -5.0, T, 120.0).J;
        x.push_back(1.0 / F);
        y.push_back(std::log(J / (F * F)));
    }
    const std::size_t n = x.size();
    double mx = 0, my = 0;
    for (std::size_t i = 0; i < n; ++i) { mx += x[i] / n; my += y[i] / n; }
    double sxy = 0, sxx = 0;
    for (std::size_t i = 0; i < n; ++i) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) * (x[i] - mx); }
    const double B = (4.0 / 3.0) * wkb_action_scale() * std::sqrt(m) * std::pow(phi, 1.5);
    REQUIRE(close(-sxy / sxx, B, 0.02));
}

TEST_CASE("gate tunnelling: GC3 the direct-tunnelling thickness slope") {
    // Source band edge at 0, Boltzmann supply (F1 = -0.2 eV), 0.5 V across the oxide.
    const double phi = 3.0, m = 0.5, T = 300.0, V = 0.5;
    const auto J = [&](double t) {
        const std::vector<BarrierSegment> b{{phi, phi - V, t, m}};
        return tsu_esaki_current(b, {0.0, inf}, -0.2, -0.2 - V, T, 120.0).J;
    };
    const double t1 = 1.2 * nm, t2 = 2.0 * nm;
    const double slope = (std::log(J(t2)) - std::log(J(t1))) / (t2 - t1);
    const SegmentAction g = segment_action(phi, phi - V, 1.0);
    const double kappa = wkb_action_scale() * std::sqrt(m) * g.value;
    REQUIRE(close(slope, -2.0 * kappa, 0.03));
}

TEST_CASE("gate tunnelling: GC4 the supply and equilibrium") {
    const double kT = base::k_B_eV_per_K * 300.0;
    // Both levels 10 kT below E: the Boltzmann closed form kT (e^((F1-E)/kT) - e^((F2-E)/kT)).
    const double E = 1.0, F1 = E - 10.0 * kT, F2 = E - 12.0 * kT;
    const double boltz = kT * (std::exp((F1 - E) / kT) - std::exp((F2 - E) / kT));
    REQUIRE(close(tsu_esaki_supply(E, F1, F2, kT).value, boltz, 5e-5));  // ln(1+x) = x - x^2/2
    REQUIRE(tsu_esaki_supply(E, 0.3, 0.3, kT).value == 0.0);
    const std::vector<BarrierSegment> b{{3.0, 2.0, 1.5 * nm, 0.5}};
    const double J1 = tsu_esaki_current(b, {0.0, inf}, 0.1, 0.1 - 1.0, 300.0, 120.0).J;
    const double J0 = tsu_esaki_current(b, {0.0, inf}, 0.1, 0.1, 300.0, 120.0).J;
    REQUIRE(J1 > 0.0);
    REQUIRE(std::abs(J0) <= 1e-14 * J1);
    // No over-barrier current (GD3): a flat barrier at E_c leaves an empty range.
    const std::vector<BarrierSegment> none{{0.0, 0.0, 1.0 * nm, 0.5}};
    const double T = 300.0, A = 120.0, Fs = -0.3;
    const double Jt = tsu_esaki_current(none, {0.0, 40.0}, Fs, -10.0, T, A).J;
    REQUIRE(Jt == 0.0);  // the barrier top (0) caps the range at E_low: no over-barrier current
}

TEST_CASE("gate tunnelling: GC5 the energy quadrature against a fine reference") {
    const double T = 300.0, A = 120.0, kT = base::k_B_eV_per_K * T;
    const std::vector<std::vector<BarrierSegment>> barriers{
        {{3.0, 2.0, 1.5 * nm, 0.5}},
        {{3.0, -6.0, 10.0 * nm, 0.5}},
        {{2.9, 2.6, 0.5 * nm, 0.5}, {1.4, 0.6, 1.5 * nm, 0.2}},
    };
    for (const auto& b : barriers) {
        for (auto [F1, F2] : {std::pair{-0.2, -1.2}, {0.1, -0.9}, {-0.3, 0.2}}) {
            const double J = tsu_esaki_current(b, {0.0, inf}, F1, F2, T, A).J;
            double top = -inf;
            for (const auto& s : b) top = std::max({top, s.W0_eV, s.W1_eV});
            const double hi = std::min(top, std::max({F1, F2, 0.0}) + 40.0 * kT);
            const double R = A / (base::k_B_eV_per_K * base::k_B_eV_per_K) *
                             gauss([&](double E) {
                                 return wkb_transmission(b, E) * tsu_esaki_supply(E, F1, F2, kT).value;
                             }, 0.0, hi, 4000);
            REQUIRE(close(J, R, 1e-6));
        }
    }
}

TEST_CASE("gate tunnelling: the current's partials against central differences") {
    const double T = 300.0, A = 120.0;
    std::vector<BarrierSegment> b{{2.9, 2.6, 0.5 * nm, 0.5}, {1.4, 0.6, 1.5 * nm, 0.2}};
    struct Case { double El, Eh, F1, F2; };
    for (const Case c : {Case{0.0, inf, -0.2, -1.2}, Case{0.0, 0.5, 0.1, -0.9},
                         Case{-0.4, inf, -0.3, 0.2}, Case{0.0, inf, 1.0, 0.9}}) {
        std::vector<double> dW(4, 0.0);
        const TunnelCurrent J = tsu_esaki_current(b, {c.El, c.Eh}, c.F1, c.F2, T, A, dW);
        const auto at = [&](double El, double Eh, double F1, double F2, const std::vector<BarrierSegment>& bb) {
            return tsu_esaki_current(bb, {El, Eh}, F1, F2, T, A).J;
        };
        const double h = 1e-6;
        const double scale = std::abs(J.J) / 0.0259;
        const auto check = [&](double analytic, double fd) {
            REQUIRE(std::abs(analytic - fd) <= 1e-6 * std::max(std::abs(fd), scale));
        };
        check(J.d_F1, (at(c.El, c.Eh, c.F1 + h, c.F2, b) - at(c.El, c.Eh, c.F1 - h, c.F2, b)) / (2 * h));
        check(J.d_F2, (at(c.El, c.Eh, c.F1, c.F2 + h, b) - at(c.El, c.Eh, c.F1, c.F2 - h, b)) / (2 * h));
        check(J.d_E_low, (at(c.El + h, c.Eh, c.F1, c.F2, b) - at(c.El - h, c.Eh, c.F1, c.F2, b)) / (2 * h));
        if (std::isfinite(c.Eh))
            check(J.d_E_high, (at(c.El, c.Eh + h, c.F1, c.F2, b) - at(c.El, c.Eh - h, c.F1, c.F2, b)) / (2 * h));
        for (std::size_t k = 0; k < 4; ++k) {
            auto bp = b, bm = b;
            double& wp = (k % 2 == 0) ? bp[k / 2].W0_eV : bp[k / 2].W1_eV;
            double& wm = (k % 2 == 0) ? bm[k / 2].W0_eV : bm[k / 2].W1_eV;
            wp += h;
            wm -= h;
            check(dW[k], (at(c.El, c.Eh, c.F1, c.F2, bp) - at(c.El, c.Eh, c.F1, c.F2, bm)) / (2 * h));
        }
    }
}

TEST_CASE("gate tunnelling: insulator band data") {
    REQUIRE_FALSE(has_band_data(silicon_dioxide_parameters));
    REQUIRE(has_band_data(silicon_dioxide_tunnelling_parameters));
    REQUIRE(evb_mass(silicon_dioxide_tunnelling_parameters) == 0.50);
    REQUIRE(Insulator::create(silicon_dioxide_tunnelling_parameters).has_value());
    InsulatorParameters p = silicon_dioxide_parameters;
    p.band_gap_eV = 5.7;  // gap without masses
    auto e = Insulator::create(p);
    REQUIRE_FALSE(e.has_value());
    REQUIRE(e.error().code == base::ErrorCode::invalid_input);
    p = silicon_dioxide_tunnelling_parameters;
    p.tunnel_mass_p = -0.1;
    REQUIRE_FALSE(Insulator::create(p).has_value());
    p = silicon_dioxide_tunnelling_parameters;
    p.electron_affinity_eV = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_FALSE(Insulator::create(p).has_value());
}
