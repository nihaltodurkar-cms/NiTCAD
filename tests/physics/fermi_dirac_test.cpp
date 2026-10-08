// The Fermi-Dirac integral F_{1/2} and its inverse (ARCHITECTURE.md section 11, Unit 14; legacy
// tests/test_m13_fermi.py, gates G1-G3).
//
// References computed here in double-double arithmetic (../physics/references.hpp): for
// eta < -1 the polylogarithm series F_j = -Li_{j+1}(-e^eta), otherwise composite Gauss-Legendre
// quadrature of the defining integrals after t = u^2. The two methods are checked against each
// other where both converge, and the quadrature against the closed form
// F_{1/2}(0) = (1 - 2^(-1/2)) zeta(3/2), zeta by Euler-Maclaurin. Independent of this code
// (different arithmetic, different discretization). The points avoid the table nodes, so they
// measure the interpolation and not only the quadrature it is built from.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include "NiTCAD/physics/fermi_dirac.hpp"
#include "references.hpp"

using namespace NiTCAD::physics;

namespace {

struct Reference {
    double eta, f12, fm12;  // F_{1/2}(eta), F_{-1/2}(eta)
};

constexpr double points[] = {-700,     -40,     -12.5, -5,      -2.000001, -1.987654,
                             -1.2345,  -0.31,   0,     0.0123,  0.777,     1.5,
                             2.71828,  4.444,   7.3,   12.34567, 19.99,    27.1828,
                             39.98765, 40.5,    60,    100,     1000};

const std::vector<Reference>& references() {
    static const std::vector<Reference> r = [] {
        std::vector<Reference> v;
        for (const double eta : points) {
            v.push_back({eta, reference::fermi_half(eta).value(),
                         reference::fermi_minus_half(eta).value()});
        }
        return v;
    }();
    return r;
}

double relative(double a, double b) { return std::abs(a - b) / std::abs(b); }

constexpr double eps = std::numeric_limits<double>::epsilon();

}  // namespace

TEST_CASE("Fermi-Dirac: the reference methods agree with each other and the closed form") {
    // Series and quadrature where both converge; the quadrature at 0 against
    // (1 - 2^(-1/2)) zeta(3/2).
    for (const double eta : {-12.5, -5.0, -2.000001, -1.2345, -1.0001}) {
        CAPTURE(eta);
        const reference::DD s = reference::fermi_half_series(eta);
        const reference::DD q = reference::fermi_half_quadrature(eta);
        REQUIRE(std::abs((s - q).value()) <= 1e-26 * std::abs(s.value()));
    }
    const reference::DD zeta = reference::zeta(1.5);
    const reference::DD f0 =
        (reference::DD(1.0) - reference::DD(1.0) / reference::sqrt(2.0)) * zeta;
    CAPTURE(zeta.value(), f0.value());
    REQUIRE(std::abs((reference::fermi_half(0.0) - f0).value()) <= 1e-26);
}

TEST_CASE("Fermi-Dirac: F_1/2 and F_-1/2 against double-double references") {
    double worst_value = 0.0, worst_derivative = 0.0;
    for (const Reference& r : references()) {
        CAPTURE(r.eta);
        const FermiIntegral f = fermi_half(r.eta);
        worst_value = std::max(worst_value, relative(f.value, r.f12));
        worst_derivative = std::max(worst_derivative, relative(f.derivative, r.fm12));
        REQUIRE(relative(f.value, r.f12) <= 2e-15);
        // F_{-1/2} is the table interpolant's own derivative: interpolation and the rounding of the
        // stored ln gamma, over a slope F_{-1/2} / F_{1/2} that falls to 0.04 at eta = 40.
        REQUIRE(relative(f.derivative, r.fm12) <= (r.eta <= 15.0 ? 2e-13 : 5e-12));
    }
    UNSCOPED_INFO("worst relative error: F_1/2 " << worst_value << ", F_-1/2 " << worst_derivative);
    // The exact anchor at eta = 0: (1 - 2^(-1/2)) zeta(3/2), zeta by Euler-Maclaurin.
    const double anchor =
        ((reference::DD(1.0) - reference::DD(1.0) / reference::sqrt(2.0)) * reference::zeta(1.5))
            .value();
    REQUIRE(relative(fermi_half(0.0).value, anchor) <= 2e-15);
}

TEST_CASE("Fermi-Dirac: the returned derivative is that of the value, across both joins") {
    // Central differences of F_{1/2} against F_{-1/2}, on a grid that crosses the series/table join
    // at -2, table cells, and the table/Sommerfeld join at 40. Step 1e-4: truncation below
    // h^2 / 6 = 2e-9; rounding up to about 2e-9 near eta = 40 (ln gamma = -35 there, so its
    // rounding over a slope of 0.04).
    double worst = 0.0;
    for (double eta = -30.0; eta <= 120.0; eta += 0.0731) {
        for (const double at : {eta, -2.0 + 1e-7 * eta, 40.0 + 1e-7 * eta}) {
            const double h = 1e-4;
            const double fd = (fermi_half(at + h).value - fermi_half(at - h).value) / (2.0 * h);
            const double error = relative(fd, fermi_half(at).derivative);
            worst = std::max(worst, error);
            CAPTURE(at);
            REQUIRE(error <= 1e-8);
        }
    }
    UNSCOPED_INFO("worst derivative mismatch " << worst);
}

TEST_CASE("Fermi-Dirac: ln gamma is F_1/2 / e^eta, continuous through the joins") {
    for (const Reference& r : references()) {
        if (r.eta < -600.0) continue;  // F underflows towards the subnormals there
        CAPTURE(r.eta);
        const LogDegeneracy L = log_degeneracy(r.eta);
        const FermiIntegral f = fermi_half(r.eta);
        REQUIRE(std::abs(L.value - (std::log(f.value) - r.eta)) <= 1e-14 * std::max(1.0, r.eta));
        REQUIRE(std::abs((1.0 + L.d_eta) * f.value - f.derivative) <= 1e-15 * f.derivative);
        REQUIRE(L.value <= 0.0);
        REQUIRE(L.d_eta <= 0.0);
        REQUIRE(L.d_eta > -1.0);
    }
    // ln gamma deep in the Boltzmann regime is e^eta / 2^(3/2) below zero, not rounded to zero.
    REQUIRE(relative(-log_degeneracy(-100.0).value, std::exp(-100.0) / std::pow(2.0, 1.5)) <=
            1e-14);
    REQUIRE(log_degeneracy(-std::numeric_limits<double>::infinity()).value == 0.0);
    // The joins: value and slope just below and just above.
    for (const double join : {-2.0, 40.0}) {
        CAPTURE(join);
        const double below = std::nextafter(join, -1e9), above = std::nextafter(join, 1e9);
        const LogDegeneracy a = log_degeneracy(below), b = log_degeneracy(above);
        // ln gamma = ln F - eta carries the rounding of eta (Sommerfeld forms it by that sum).
        REQUIRE(std::abs(a.value - b.value) <= 4.0 * eps * std::max(1.0, std::abs(join)));
        REQUIRE(std::abs(a.d_eta - b.d_eta) <= 4e-16);
    }
}

TEST_CASE("Fermi-Dirac: increasing, log-concave and smooth on a dense grid") {
    // Legacy G1 smoothness guard: second differences of ln F_{1/2} on a 1e-3 grid. Here they must
    // also be negative to rounding (ln F_{1/2} is concave), which a wiggle in the interpolant
    // would break.
    double previous = fermi_half(-50.0 + 1e-3).value;
    double lf0 = std::log(fermi_half(-50.0).value), lf1 = std::log(previous);
    for (int k = 2; k <= 150000; ++k) {
        const double eta = -50.0 + 1e-3 * k;
        const double f = fermi_half(eta).value;
        const double lf2 = std::log(f);
        const double rounding = 8.0 * std::numeric_limits<double>::epsilon() *
                                std::max(1.0, std::abs(lf2));
        CAPTURE(eta);
        REQUIRE(f > previous);
        REQUIRE(lf2 - 2.0 * lf1 + lf0 <= rounding);
        REQUIRE(lf2 - 2.0 * lf1 + lf0 > -1e-6);
        previous = f;
        lf0 = lf1;
        lf1 = lf2;
    }
}

TEST_CASE("Fermi-Dirac: the Boltzmann and degenerate limits") {
    // G2: F_{1/2} = e^eta (1 - e^eta / 2^(3/2) + e^2eta / 3^(3/2) - ...), so
    // ln gamma = -e^eta / 2^(3/2) (1 - 0.367 e^eta + ...), and F_{1/2} = e^eta to rounding below
    // eta = -36.
    for (const double eta : {-30.0, -20.0, -15.0, -10.0}) {
        CAPTURE(eta);
        const double expected = -std::exp(eta) / std::pow(2.0, 1.5);
        const double error = relative(log_degeneracy(eta).value, expected);
        REQUIRE(error <= 0.37 * std::exp(eta) + 1e-15);
        REQUIRE(error >= 0.36 * std::exp(eta) - 1e-15);
    }
    REQUIRE(relative(fermi_half(-40.0).value, std::exp(-40.0)) <= 2e-16);
    // G3: the two-term Sommerfeld form leaves the 7 pi^4 / (640 eta^4) term, so its deviation
    // falls by 16 per doubling of eta.
    const auto sommerfeld2 = [](double eta) {
        return 4.0 / (3.0 * std::sqrt(std::numbers::pi)) * std::pow(eta, 1.5) *
               (1.0 + std::numbers::pi * std::numbers::pi / (8.0 * eta * eta));
    };
    const double d20 = relative(sommerfeld2(20.0), fermi_half(20.0).value);
    const double d40 = relative(sommerfeld2(40.0), fermi_half(40.0).value);
    const double third = 7.0 * std::pow(std::numbers::pi, 4) / 640.0;
    REQUIRE(relative(d40, third / std::pow(40.0, 4)) <= 0.02);
    REQUIRE(d20 / d40 > 15.0);
    REQUIRE(d20 / d40 < 17.0);
}

TEST_CASE("Fermi-Dirac: the inverse recovers eta") {
    double worst = 0.0;
    for (double eta = -700.0; eta <= 1000.0; eta += eta < -50.0 || eta > 50.0 ? 13.7 : 0.0917) {
        const double nu = fermi_half(eta).value;
        const double back = inverse_fermi_half(nu);
        const double error = std::abs(back - eta) / std::max(1.0, std::abs(eta));
        worst = std::max(worst, error);
        CAPTURE(eta);
        REQUIRE(error <= 1e-14);
    }
    UNSCOPED_INFO("worst inverse error " << worst);
    // Exactly at a reference and at the extremes of double range.
    REQUIRE(std::abs(inverse_fermi_half(0.7651470246254079453673)) <= 2e-15);
    const double tiny = inverse_fermi_half(1e-300);
    REQUIRE(std::abs(tiny - std::log(1e-300)) <= 1e-13);
    const double huge = inverse_fermi_half(1e30);
    REQUIRE(relative(fermi_half(huge).value, 1e30) <= 1e-14);
}

TEST_CASE("Fermi-Dirac: a NaN argument gives NaN") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(std::isnan(fermi_half(nan).value));
    REQUIRE(std::isnan(fermi_half(nan).derivative));
    REQUIRE(std::isnan(log_degeneracy(nan).value));
}


