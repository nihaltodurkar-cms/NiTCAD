// Boltzmann carrier statistics (ARCHITECTURE.md section 11, Unit 5; R4).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

constexpr double eps = std::numeric_limits<double>::epsilon();

}  // namespace

TEST_CASE("statistics: Boltzmann density and its derivative") {
    const double ni = 1.0674e10;
    REQUIRE(boltzmann_density(ni, 0.0).density == ni);
    for (const double eta : {-40.0, -3.0, -1e-3, 0.5, 7.0, 40.0}) {
        CAPTURE(eta);
        const DensityResult r = boltzmann_density(ni, eta);
        REQUIRE(r.density == ni * std::exp(eta));
        REQUIRE(r.d_eta == r.density);
        // Central finite difference, divided by the step actually taken (eta +- h rounds).
        const double up = eta + 1e-6;
        const double down = eta - 1e-6;
        const double fd =
            (boltzmann_density(ni, up).density - boltzmann_density(ni, down).density) /
            (up - down);
        REQUIRE(close(r.d_eta, fd, 1e-9));
    }
}

TEST_CASE("statistics: at equilibrium n p = n_ie^2 for any potential") {
    const double ni = 1.0674e10;
    for (const double eta : {-35.0, -5.0, 0.0, 2.5, 35.0}) {
        CAPTURE(eta);
        const double n = boltzmann_density(ni, eta).density;   // eta_n = psi / V_T
        const double p = boltzmann_density(ni, -eta).density;  // eta_p = -psi / V_T
        REQUIRE(close(n * p, ni * ni, 4.0 * eps));
    }
}

TEST_CASE("statistics: the Boltzmann equilibrium product has no carrier dependence") {
    constexpr EquilibriumProduct e = boltzmann_equilibrium_product(3.0);
    static_assert(e.value == 9.0 && e.d_dn == 0.0 && e.d_dp == 0.0);
    REQUIRE(boltzmann_equilibrium_product(1.0674e10).value == 1.0674e10 * 1.0674e10);
}

TEST_CASE("statistics: charge-neutral equilibrium satisfies neutrality and mass action") {
    const double ni = 1.0674e10;
    for (const double C : {0.0, 1.0, 1e10, -1e10, 1e15, -1e15, 1e17, -1e17, 1e20, -1e20, 1e22}) {
        CAPTURE(C);
        const NeutralEquilibrium e = boltzmann_neutral_equilibrium(C, ni);
        REQUIRE(e.n > 0.0);
        REQUIRE(e.p > 0.0);
        // Neutrality n - p = C, measured against the majority density (the subtraction cancels).
        REQUIRE(std::abs((e.n - e.p) - C) <= 4.0 * eps * std::max(e.n, e.p));
        REQUIRE(close(e.n * e.p, ni * ni, 4.0 * eps));
        REQUIRE(std::abs(e.eta - std::log(e.n / ni)) <= 4.0 * eps * std::max(1.0, std::abs(e.eta)));
    }
    const NeutralEquilibrium intrinsic = boltzmann_neutral_equilibrium(0.0, ni);
    REQUIRE(intrinsic.n == ni);
    REQUIRE(intrinsic.p == ni);
    REQUIRE(intrinsic.eta == 0.0);
}

TEST_CASE("statistics: n-type and p-type are mirror images, and the result is scale-free") {
    const double ni = 1.0674e10;
    for (const double C : {3e14, 1e17, 5e19}) {
        const NeutralEquilibrium nt = boltzmann_neutral_equilibrium(C, ni);
        const NeutralEquilibrium pt = boltzmann_neutral_equilibrium(-C, ni);
        REQUIRE(nt.n == pt.p);
        REQUIRE(nt.p == pt.n);
        REQUIRE(nt.eta == -pt.eta);
        // The assembler passes concentrations divided by its scale Ns (6.1); a power-of-two
        // scale is exact, so the result must scale exactly.
        const double s = 0x1p-60;
        const NeutralEquilibrium scaled = boltzmann_neutral_equilibrium(C * s, ni * s);
        REQUIRE(scaled.n == nt.n * s);
        REQUIRE(scaled.p == nt.p * s);
        REQUIRE(scaled.eta == nt.eta);
    }
}

TEST_CASE("statistics: built-in potential of a 1e17/1e17 junction") {
    // V_bi = V_T (eta(N_D) - eta(-N_A)) equals V_T ln(N_D N_A / n_i^2) up to n_i^2 / N^2.
    const double ni = intrinsic_density(silicon(), 300.0);
    const double VT = NiTCAD::base::thermal_voltage(300.0);
    const double N = 1e17;
    const double vbi = VT * (boltzmann_neutral_equilibrium(N, ni).eta -
                             boltzmann_neutral_equilibrium(-N, ni).eta);
    REQUIRE(close(vbi, VT * std::log(N * N / (ni * ni)), 1e-14));
    REQUIRE(close(vbi, 0.82999866587, 1e-10));
}

TEST_CASE("statistics: extreme doping does not overflow") {
    const NeutralEquilibrium e = boltzmann_neutral_equilibrium(1e300, 1.0);
    REQUIRE(e.n == 1e300);
    REQUIRE(close(e.p, 1e-300, eps));
    REQUIRE(std::isfinite(e.eta));
}

TEST_CASE("statistics: neutral equilibrium with |C| near DBL_MAX does not overflow") {
    // The legacy 0.5 (C + sqrt(C^2 + 4 n_ie^2)) overflows here.
    constexpr double big = std::numeric_limits<double>::max();
    for (const double C : {big, -big, 1e308, -1e308}) {
        CAPTURE(C);
        const NeutralEquilibrium e = boltzmann_neutral_equilibrium(C, 1.0);
        const double majority = C > 0.0 ? e.n : e.p;
        const double minority = C > 0.0 ? e.p : e.n;
        REQUIRE(majority == std::abs(C));
        REQUIRE(minority > 0.0);
        REQUIRE(close(minority, 1.0 / std::abs(C), 1e-14));  // subnormal, so 1e-14
        REQUIRE(close(e.eta, std::copysign(std::log(std::abs(C)), C), 4.0 * eps));
    }
}

TEST_CASE("statistics: eta when C / (2 n_ie) overflows") {
    const NeutralEquilibrium e = boltzmann_neutral_equilibrium(1e300, 1e-10);
    REQUIRE(e.n == 1e300);
    REQUIRE(std::isfinite(e.eta));
    REQUIRE(close(e.eta, std::log(1e300) - std::log(1e-10), 4.0 * eps));
    REQUIRE(boltzmann_neutral_equilibrium(-1e300, 1e-10).eta == -e.eta);
    // The fallback continues the asinh branch: where the quotient (5e306) is still finite, asinh
    // gives the same value as ln|C| - ln n_ie.
    const double C = 1e300;
    const double ni = 1e-7;
    REQUIRE(close(boltzmann_neutral_equilibrium(C, ni).eta, std::log(C) - std::log(ni), 4.0 * eps));
}

TEST_CASE("statistics: neutral equilibrium when n_ie^2 underflows") {
    // n_ie = 1e-170: n_ie^2 is 0 in double, yet p = n_ie^2 / C = 1e-300 is representable.
    // The legacy n_ie^2 / n returned p = 0.
    const double ni = 1e-170;
    REQUIRE(boltzmann_equilibrium_product(ni).value == 0.0);  // the documented limit
    const NeutralEquilibrium e = boltzmann_neutral_equilibrium(1e-40, ni);
    REQUIRE(e.n == 1e-40);
    REQUIRE(close(e.p, 1e-300, 4.0 * eps));
    REQUIRE(close((e.n / ni) * (e.p / ni), 1.0, 8.0 * eps));  // mass action in units of n_ie
    REQUIRE(close(e.eta, std::log(1e-40 / ni), 4.0 * eps));
    const NeutralEquilibrium h = boltzmann_neutral_equilibrium(-1e-40, ni);
    REQUIRE(h.n == e.p);
    REQUIRE(h.p == e.n);
    REQUIRE(h.eta == -e.eta);
}
