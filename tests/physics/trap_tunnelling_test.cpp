// Hurkx trap-assisted tunnelling field-enhancement factor (Unit 67, DECISIONS.md TAC2-TAC5): the
// zero-field and monotone limits, the quadrature against an independent reference, the
// low-field and Laplace asymptotes and the field derivative.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/trap_tunnelling.hpp"

using namespace NiTCAD::physics;
namespace base = NiTCAD::base;

namespace {

constexpr double mass = 0.25;  // S4 Table 3-7 default
const double as[] = {10.0, 21.7, 40.0, 62.0};
const double Ts[] = {250.0, 300.0, 400.0};

// K of the exponent at field F [V/cm]: (4/3) sqrt(2 m m0 dE^3) / (q hbar F), dE = a kT.
double kappa(double a, double m, double T, double F_V_per_cm) {
    const double dE = a * base::k_B_J_per_K * T;
    return 4.0 / 3.0 * std::sqrt(2.0 * m * base::m0_kg * dE * dE * dE) /
           (base::q_C * base::hbar_J_s * F_V_per_cm * 100.0);
}

// The field [V/cm] at which K takes the given value.
double field_for(double a, double m, double T, double K) {
    return kappa(a, m, T, 1.0) / K;
}

// F_Gamma = sqrt(24 m m0 (kT)^3) / (q hbar) [V/cm].
double f_gamma(double m, double T) {
    const double kT = base::k_B_J_per_K * T;
    return std::sqrt(24.0 * m * base::m0_kg * kT * kT * kT) / (base::q_C * base::hbar_J_s) / 100.0;
}

// Independent reference: Gamma = 2a int_0^1 v exp(a v^2 - K v^3) dv (u = v^2) and
// dGamma/dF = (K/F) 2a int_0^1 v^4 exp(...) dv, composite Simpson with 2^16 intervals.
struct Ref {
    double gamma, d_dF;
};
Ref reference(double a, double m, double T, double F) {
    const double K = kappa(a, m, T, F);
    const int n = 1 << 16;
    const double h = 1.0 / n;
    double s0 = 0.0, s4 = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double v = i * h;
        const double w = (i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
        const double e = std::exp(a * v * v - K * v * v * v);
        s0 += w * v * e;
        s4 += w * v * v * v * v * e;
    }
    return {2.0 * a * s0 * h / 3.0, K / F * 2.0 * a * s4 * h / 3.0};
}

std::vector<double> field_grid() {  // 40 log points 1e3 .. 1e7 V/cm
    std::vector<double> F;
    for (int i = 0; i < 40; ++i) F.push_back(1e3 * std::pow(1e4, i / 39.0));
    return F;
}

}  // namespace

TEST_CASE("hurkx zero field", "[trap_tunnelling]") {
    for (const double a : as) {
        for (const double T : Ts) {
            const TrapTunnellingFactor g = hurkx_enhancement(a, mass, T, 0.0);
            CHECK(g.gamma == 0.0);
            CHECK(g.d_dF == 0.0);
        }
    }
}

TEST_CASE("hurkx monotone and bounded", "[trap_tunnelling]") {
    for (const double a : as) {
        for (const double T : Ts) {
            double last = 0.0;
            for (const double F : field_grid()) {
                const TrapTunnellingFactor g = hurkx_enhancement(a, mass, T, F);
                INFO("a " << a << " T " << T << " F " << F);
                CHECK(g.gamma > last);
                CHECK(g.d_dF > 0.0);
                CHECK(g.gamma <= std::expm1(a));
                last = g.gamma;
            }
        }
    }
}

TEST_CASE("hurkx quadrature reference", "[trap_tunnelling]") {
    double worst = 0.0;
    for (const double a : as) {
        for (const double T : Ts) {
            for (const double F : field_grid()) {
                const TrapTunnellingFactor g = hurkx_enhancement(a, mass, T, F);
                const Ref r = reference(a, mass, T, F);
                INFO("a " << a << " T " << T << " F " << F);
                const double eg = std::abs(g.gamma - r.gamma) / r.gamma;
                const double ed = std::abs(g.d_dF - r.d_dF) / r.d_dF;
                CHECK(eg <= 1e-12);
                CHECK(ed <= 1e-12);
                worst = std::max({worst, eg, ed});
            }
        }
    }
    WARN("TAC3 largest relative deviation " << worst);
}

TEST_CASE("hurkx low-field limit", "[trap_tunnelling]") {
    // K -> infinity: Gamma -> a Gamma(5/3) K^(-2/3); taken where a K^(-2/3) = 1e-4.
    for (const double a : as) {
        for (const double T : Ts) {
            const double K = std::pow(1e4 * a, 1.5);
            const double F = field_for(a, mass, T, K);
            const double limit = a * std::tgamma(5.0 / 3.0) * std::pow(K, -2.0 / 3.0);
            const TrapTunnellingFactor g = hurkx_enhancement(a, mass, T, F);
            INFO("a " << a << " T " << T << " F " << F);
            CHECK(std::abs(g.gamma - limit) <= 1e-3 * limit);
        }
    }
}

TEST_CASE("hurkx laplace asymptote", "[trap_tunnelling]") {
    // Peak at v* = 2a/(3K) = 0.7: Gamma -> 2 sqrt(3 pi) (F/F_G) exp((F/F_G)^2) as a grows.
    CHECK(std::abs(f_gamma(mass, 300.0) - 3.69e5) <= 0.005 * 3.69e5);
    const double T = 300.0;
    double last = 1e300;
    for (const double a : {40.0, 62.0, 100.0, 200.0}) {
        const double F = field_for(a, mass, T, 2.0 * a / (3.0 * 0.7));
        const double x = F / f_gamma(mass, T);
        const double laplace = 2.0 * std::sqrt(3.0 * std::numbers::pi) * x * std::exp(x * x);
        const double dev =
            std::abs(hurkx_enhancement(a, mass, T, F).gamma - laplace) / laplace;
        INFO("a " << a);
        CHECK(dev < last);
        WARN("TAC4 a " << a << " Laplace deviation " << dev);
        last = dev;
    }
    CHECK(last <= 0.05);
}

TEST_CASE("hurkx derivative", "[trap_tunnelling]") {
    for (const double a : as) {
        for (const double T : Ts) {
            for (const double F : field_grid()) {
                const TrapTunnellingFactor g = hurkx_enhancement(a, mass, T, F);
                if (g.gamma <= 1e-12) continue;
                // Fourth-order central difference: Gamma ~ exp((F/F_G)^2) is steep enough that
                // the second-order one's truncation reaches 1e-5 on this grid.
                const double h = 1e-4 * F;
                const auto G = [&](double f) { return hurkx_enhancement(a, mass, T, f).gamma; };
                const double fd =
                    (-G(F + 2.0 * h) + 8.0 * G(F + h) - 8.0 * G(F - h) + G(F - 2.0 * h)) /
                    (12.0 * h);
                INFO("a " << a << " T " << T << " F " << F);
                CHECK(std::abs(g.d_dF - fd) <= 1e-6 * std::abs(fd));
            }
        }
    }
}

TEST_CASE("tat extreme field finite", "[trap_tunnelling]") {
    for (const double a : as) {
        for (const double F : {1e7, 1e8, 1e10}) {
            const TrapTunnellingFactor g = hurkx_enhancement(a, mass, 300.0, F);
            CHECK(std::isfinite(g.gamma));
            CHECK(std::isfinite(g.d_dF));
            CHECK(g.gamma <= std::expm1(a));
        }
        // A field so small that K overflows: Gamma 0, no NaN.
        const TrapTunnellingFactor tiny = hurkx_enhancement(a, mass, 300.0, 1e-300);
        CHECK(tiny.gamma == 0.0);
        CHECK(tiny.d_dF == 0.0);
    }
}
