// Bernoulli function and Scharfetter-Gummel edge fluxes (ARCHITECTURE.md section 11, Unit 7 gate:
// Bernoulli limits and symmetry). References computed here in double-double arithmetic
// (../physics/references.hpp) from x / expm1(x) and its derivative, independently of this code.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

#include "NiTCAD/assemble/bernoulli.hpp"
#include "NiTCAD/assemble/sg_flux.hpp"
#include "../physics/references.hpp"

using namespace NiTCAD::assemble;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

constexpr double eps = std::numeric_limits<double>::epsilon();

constexpr double points[] = {-1000, -800,  -40,    -1, -0.01, -0.0099, -1e-4, 1e-6, 1e-4,
                             0.0099, 0.01, 0.5,    1,  3,     40,      700,   709};

}  // namespace

TEST_CASE("bernoulli: values and derivatives against double-double references") {
    for (const double x : points) {
        CAPTURE(x);
        const reference::Bernoulli r = reference::bernoulli(x);
        REQUIRE(close(bernoulli(x), r.value.value(), 4.0 * eps));
        REQUIRE(close(bernoulli_derivative(x), r.derivative.value(), 1e-13));
    }
    REQUIRE(bernoulli(0.0) == 1.0);
    REQUIRE(bernoulli_derivative(0.0) == -0.5);
}

TEST_CASE("bernoulli: symmetry B(-x) = B(x) + x and B'(x) + B'(-x) = -1") {
    for (double x = 1e-8; x < 700.0; x *= 1.37) {
        CAPTURE(x);
        REQUIRE(close(bernoulli(-x), bernoulli(x) + x, 4.0 * eps));
        REQUIRE(std::abs(bernoulli_derivative(x) + bernoulli_derivative(-x) + 1.0) <= 1e-13);
    }
}

TEST_CASE("bernoulli: limits, with no clipping") {
    // The legacy clipped x to [-700, 700]: B(-1000) was 700.
    REQUIRE(bernoulli(-1000.0) == 1000.0);
    REQUIRE(bernoulli(-1e6) == 1e6);
    REQUIRE(bernoulli_derivative(-1e6) == -1.0);
    REQUIRE(bernoulli(1000.0) == 0.0);
    REQUIRE(bernoulli_derivative(1000.0) == 0.0);
    // Positive and decreasing everywhere.
    double previous = bernoulli(-50.0);
    for (double x = -50.0; x <= 50.0; x += 0.0625) {
        const double b = bernoulli(x);
        REQUIRE(b > 0.0);
        REQUIRE(b <= previous);
        previous = b;
    }
}

TEST_CASE("bernoulli: continuous across the series switch at |x| = 1e-2") {
    for (const double edge : {1e-2, -1e-2}) {
        const double below = std::nextafter(edge, 0.0);
        REQUIRE(close(bernoulli(below), bernoulli(edge), 4.0 * eps));
        REQUIRE(close(bernoulli_derivative(below), bernoulli_derivative(edge), 1e-13));
    }
}

TEST_CASE("sg flux: zero at equilibrium for both carriers") {
    const double ni = 7.3e-7;  // a scaled n_ie
    for (const double psi1 : {-15.0, -2.0, 0.0, 4.0}) {
        for (const double d : {-30.0, -1.0, -1e-3, 0.0, 2e-5, 0.7, 25.0}) {
            const double psi2 = psi1 + d;
            const double n1 = ni * std::exp(psi1), n2 = ni * std::exp(psi2);
            const double p1 = ni * std::exp(-psi1), p2 = ni * std::exp(-psi2);
            const double a = 3.0;
            // Relative to the larger of the two one-sided terms.
            const double scale_n = a * std::max(n2 * bernoulli(d), n1 * bernoulli(-d));
            const double scale_p = a * std::max(p2 * bernoulli(-d), p1 * bernoulli(d));
            CAPTURE(psi1, d);
            REQUIRE(std::abs(sg_electron_flux(a, psi1, psi2, n1, n2).flux) <= 1e-13 * scale_n);
            REQUIRE(std::abs(sg_hole_flux(a, psi1, psi2, p1, p2).flux) <= 1e-13 * scale_p);
        }
    }
}

TEST_CASE("sg flux: diffusion and drift limits, and reversal") {
    // No field: pure diffusion, Jn = a (n2 - n1), Jp = -a (p2 - p1).
    REQUIRE(sg_electron_flux(2.0, 1.0, 1.0, 3.0, 5.0).flux == 4.0);
    REQUIRE(sg_hole_flux(2.0, 1.0, 1.0, 3.0, 5.0).flux == -4.0);
    // Strong field (delta = 60): upwinded drift, electrons carried from node 1 (Jn -> -a n1 delta)
    // and holes from node 2 (Jp -> -a p2 delta); the downwind density drops out.
    REQUIRE(close(sg_electron_flux(1.0, 0.0, 60.0, 2.0, 1e3).flux, -2.0 * 60.0, 1e-12));
    REQUIRE(close(sg_hole_flux(1.0, 0.0, 60.0, 1e3, 2.0).flux, -2.0 * 60.0, 1e-12));
    // Swapping the ends negates the flux exactly.
    const EdgeFlux f = sg_electron_flux(1.5, 0.3, -2.0, 0.7, 4.0);
    const EdgeFlux r = sg_electron_flux(1.5, -2.0, 0.3, 4.0, 0.7);
    REQUIRE(r.flux == -f.flux);
    const EdgeFlux g = sg_hole_flux(1.5, 0.3, -2.0, 0.7, 4.0);
    REQUIRE(sg_hole_flux(1.5, -2.0, 0.3, 4.0, 0.7).flux == -g.flux);
}

TEST_CASE("sg flux: partials match finite differences") {
    struct State {
        double psi1, psi2, c1, c2;
    };
    constexpr State states[] = {
        {0.0, 0.0, 1.0, 2.0},  {0.1, -0.3, 0.5, 4.0}, {-5.0, 7.0, 1e-3, 3.0},
        {2.0, 2.005, 1.0, 1.0}, {10.0, -12.0, 2.0, 1e-6},
    };
    double worst = 0.0;
    for (const auto& flux : {sg_electron_flux, sg_hole_flux}) {
        for (const State& s : states) {
            const EdgeFlux f = flux(1.7, s.psi1, s.psi2, s.c1, s.c2);
            const double analytic[] = {f.d_psi1, f.d_psi2, f.d_c1, f.d_c2};
            const double scale = std::max({std::abs(f.d_psi1), std::abs(f.d_c1), std::abs(f.d_c2)});
            for (int k = 0; k < 4; ++k) {
                double u[] = {s.psi1, s.psi2, s.c1, s.c2};
                const double h = 1e-6 * std::max(std::abs(u[k]), 1.0);
                const double base = u[k];
                u[k] = base + h;
                const double up = u[k];
                const double fp = flux(1.7, u[0], u[1], u[2], u[3]).flux;
                u[k] = base - h;
                const double down = u[k];
                const double fm = flux(1.7, u[0], u[1], u[2], u[3]).flux;
                worst = std::max(worst, std::abs((fp - fm) / (up - down) - analytic[k]) / scale);
            }
            REQUIRE(f.d_psi2 == -f.d_psi1);
        }
    }
    REQUIRE(worst <= 1e-8);
}
