// The Fermi-Dirac integral F_{1/2} and its inverse (ARCHITECTURE.md section 11, Unit 14; legacy
// tests/test_m13_fermi.py, gates G1-G3).
//
// Reference values, 40-digit mpmath: for eta < 0 the polylogarithm F_j = -Li_{j+1}(-e^eta); for
// eta >= 0 quadratures of the defining integrals subdivided at the Fermi edge t ~ eta. The two
// agree to 1e-31 where both were evaluated (-40 <= eta < 0; at -700 the quadrature loses 7
// digits, hence the polylogarithm), and the quadrature gives F_{1/2}(0) = (1 - 2^(-1/2)) zeta(3/2)
// and F_{-1/2}(0) = (1 - 2^(1/2)) zeta(1/2) to all digits.
// Independent of this code (different arithmetic, different discretization). The points avoid the
// table nodes, so they measure the interpolation and not only the quadrature it is built from.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include "NiTCAD/physics/fermi_dirac.hpp"

using namespace NiTCAD::physics;

namespace {

struct Reference {
    double eta, f12, fm12;  // F_{1/2}(eta), F_{-1/2}(eta)
};

constexpr Reference references[] = {
    {-700, 9.859676543759770856705e-305, 9.859676543759770856705e-305},
    {-40, 4.248354255291588988948e-18, 4.248354255291588982567e-18},
    {-12.5, 0.000003726648261958989501984, 0.000003726643351849268323077},
    {-5, 0.006721954314505912707786, 0.006706019989268209127845},
    {-2.000001, 0.1292983896545103617568, 0.1236655086574104040811},
    {-1.987654, 0.1308339415807036505649, 0.1250696910250291204103},
    {-1.2345, 0.265043325548207017367, 0.2424872355199710840049},
    {-0.31, 0.59527185592054674347, 0.4931650088039558883891},
    {0, 0.7651470246254079453673, 0.6048986434216303702473},
    {0.0123, 0.7726160676922649119224, 0.6095828828964767049883},
    {0.777, 1.35765955632693248621, 0.9281156611207317590773},
    {1.5, 2.144860877583114035964, 1.24932334785271220083},
    {2.71828, 3.979869169740514112218, 1.749864442627717048176},
    {4.444, 7.511764021998324437939, 2.318135866755477192745},
    {7.3, 15.18751192060517394408, 3.022500415460546063469},
    {12.34567, 32.89714466814245236305, 3.953672833234494470665},
    {19.99, 67.44110836493078802666, 5.039752784653832820401},
    {27.1828, 106.7899664623740525882, 5.879751022508753663206},
    {39.98765, 190.3652841674622841409, 7.133554596611505259675},
    {40.5, 194.0318552090013066966, 7.179155892384824068736},
    {60, 349.7353379459762089049, 8.73938781383138152102},
    {100, 752.3455915521961188446, 11.28332744292768060324},
    {1000, 23788.35089639434090471, 35.68246764915936962222},
};

double relative(double a, double b) { return std::abs(a - b) / std::abs(b); }

constexpr double eps = std::numeric_limits<double>::epsilon();

}  // namespace

TEST_CASE("Fermi-Dirac: F_1/2 and F_-1/2 against 40-digit references") {
    double worst_value = 0.0, worst_derivative = 0.0;
    for (const Reference& r : references) {
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
    // The exact anchor at eta = 0, with zeta(3/2) = 2.612375348685488343348567...
    REQUIRE(relative(fermi_half(0.0).value, (1.0 - 1.0 / std::numbers::sqrt2) *
                                                 2.6123753486854883433) <= 2e-15);
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
    for (const Reference& r : references) {
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


