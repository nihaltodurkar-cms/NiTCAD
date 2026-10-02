// SRH recombination and Scharfetter lifetimes (ARCHITECTURE.md section 11, Unit 5): SRH vanishes
// at equilibrium, analytic limits, and the exact partials against finite differences.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

constexpr double eps = std::numeric_limits<double>::epsilon();
constexpr double ni = 1.0674e10;
constexpr double tau_n = 1e-5 / 3.0;  // silicon at 1e17 total impurity
constexpr double tau_p = 3e-6 / 3.0;

RecombinationRate srh(double n, double p) {
    return srh_recombination(n, p, boltzmann_equilibrium_product(ni), ni, tau_n, tau_p);
}

// Largest finite-difference error over the two partials, relative to the larger partial.
template <class F>
double fd_error(F rate, double n, double p, double d_dn, double d_dp) {
    const double hn = 1e-5 * std::max(n, ni);
    const double hp = 1e-5 * std::max(p, ni);
    // Divided by the steps actually taken, since n +- hn rounds.
    const double fd_n = (rate(n + hn, p) - rate(n - hn, p)) / ((n + hn) - (n - hn));
    const double fd_p = (rate(n, p + hp) - rate(n, p - hp)) / ((p + hp) - (p - hp));
    const double scale = std::max(std::abs(d_dn), std::abs(d_dp));
    return std::max(std::abs(fd_n - d_dn), std::abs(fd_p - d_dp)) / scale;
}

// Operating points: equilibrium, low and high injection on both sides, depletion (generation).
struct Point {
    double n;
    double p;
};
constexpr Point points[] = {
    {1e17, ni * ni / 1e17},  // n-type at equilibrium
    {1e17, 1e12},            // n-type, low injection
    {1e5, 1e16},             // p-type, low injection
    {3e18, 3e18},            // high injection
    {1e14, 2e13},            // moderate injection
    {1e3, 1e4},              // depletion, generation
    {0.0, 0.0},              // no carriers
    {ni, ni},                // intrinsic
};

}  // namespace

TEST_CASE("recombination: Scharfetter lifetime") {
    const Semiconductor si = silicon();
    REQUIRE(scharfetter_lifetime(si, Carrier::electron, 0.0) == 1e-5);
    REQUIRE(scharfetter_lifetime(si, Carrier::hole, 0.0) == 3e-6);
    REQUIRE(scharfetter_lifetime(si, Carrier::electron, 5e16) == 0.5e-5);
    REQUIRE(close(scharfetter_lifetime(si, Carrier::hole, 1e17), tau_p, 1e-15));
    REQUIRE(close(scharfetter_lifetime(si, Carrier::electron, 1e17), tau_n, 1e-15));
}

TEST_CASE("recombination: SRH vanishes at equilibrium (section 11 gate)") {
    // np = n_ie^2 exactly: a factor of two each way is exact in floating point.
    static_assert(srh_recombination(2.0 * ni, 0.5 * ni, boltzmann_equilibrium_product(ni), ni,
                                    tau_n, tau_p)
                      .rate == 0.0);
    // Equilibrium carriers from the Boltzmann statistics, for any doping or potential:
    // |R| is at the rounding level of n p.
    for (const double C : {0.0, 1e14, -1e14, 1e17, -1e17, 1e20, -1e20}) {
        CAPTURE(C);
        const NeutralEquilibrium e = boltzmann_neutral_equilibrium(C, ni);
        const double den = tau_p * (e.n + ni) + tau_n * (e.p + ni);
        REQUIRE(std::abs(srh(e.n, e.p).rate) <= 4.0 * eps * e.n * e.p / den);
    }
    for (const double eta : {-30.0, -2.0, 0.0, 3.0, 30.0}) {
        CAPTURE(eta);
        const double n = boltzmann_density(ni, eta).density;
        const double p = boltzmann_density(ni, -eta).density;
        const double den = tau_p * (n + ni) + tau_n * (p + ni);
        REQUIRE(std::abs(srh(n, p).rate) <= 4.0 * eps * n * p / den);
    }
}

TEST_CASE("recombination: low-injection limits give the minority-carrier lifetime") {
    // n-type: R -> dp / tau_p; p-type: R -> dn / tau_n (to n_ie / N and tau ratio * p / n).
    const double N = 1e17;
    const double dp = 1e10;
    const double p0 = ni * ni / N;
    REQUIRE(close(srh(N, p0 + dp).rate, dp / tau_p, 1e-6));
    const double dn = 1e10;
    REQUIRE(close(srh(ni * ni / N + dn, N).rate, dn / tau_n, 1e-6));
}

TEST_CASE("recombination: generation with no carriers, and the sign of R") {
    REQUIRE(close(srh(0.0, 0.0).rate, -ni / (tau_n + tau_p), 1e-15));
    REQUIRE(srh(1e3, 1e4).rate < 0.0);     // np < n_ie^2: generation
    REQUIRE(srh(1e14, 2e13).rate > 0.0);   // np > n_ie^2: recombination
}

TEST_CASE("recombination: Boltzmann form matches the legacy operation order") {
    // kernels.hpp recombination_boltzmann, SRH part, with Auger off.
    for (const Point& pt : points) {
        const double ni2 = ni * ni;
        const double excess = pt.n * pt.p - ni2;
        const double den = tau_p * (pt.n + ni) + tau_n * (pt.p + ni);
        const RecombinationRate r = srh(pt.n, pt.p);
        REQUIRE(r.rate == excess / den);
        REQUIRE(r.d_dn == (pt.p * den - excess * tau_p) / (den * den));
        REQUIRE(r.d_dp == (pt.n * den - excess * tau_n) / (den * den));
    }
}

TEST_CASE("recombination: SRH partials match finite differences") {
    // Measured worst error 6.4e-11 (central differences, step 1e-5 of the density).
    // The section 10 gate for an assembled Jacobian is 5e-5.
    for (const Point& pt : points) {
        CAPTURE(pt.n, pt.p);
        const RecombinationRate r = srh(pt.n, pt.p);
        const double err = fd_error([](double n, double p) { return srh(n, p).rate; }, pt.n, pt.p,
                                    r.d_dn, r.d_dp);
        REQUIRE(err <= 1e-8);
    }
}

TEST_CASE("recombination: partials include the chain through a carrier-dependent product") {
    // A stand-in for a Fermi-Dirac equilibrium product (R4): E(n, p) = n_ie^2 (1 + n/K)(1 + p/K).
    // The returned partials must be the total derivatives of R(n, p, E(n, p)).
    const double K = 1e18;
    const auto product = [&](double n, double p) {
        return EquilibriumProduct{ni * ni * (1.0 + n / K) * (1.0 + p / K),
                                  ni * ni * (1.0 + p / K) / K, ni * ni * (1.0 + n / K) / K};
    };
    const auto rate = [&](double n, double p) {
        return srh_recombination(n, p, product(n, p), ni, tau_n, tau_p).rate;
    };
    for (const Point& pt : points) {
        CAPTURE(pt.n, pt.p);
        const RecombinationRate r =
            srh_recombination(pt.n, pt.p, product(pt.n, pt.p), ni, tau_n, tau_p);
        const double err = fd_error(rate, pt.n, pt.p, r.d_dn, r.d_dp);
        REQUIRE(err <= 1e-8);
    }
}

TEST_CASE("recombination: SRH is homogeneous in the concentrations") {
    // The assembler may pass densities divided by its scale Ns (6.1); a power-of-two scale is
    // exact, so R scales exactly and the partials do not change.
    const double s = 0x1p-50;
    for (const Point& pt : points) {
        const RecombinationRate r = srh(pt.n, pt.p);
        const RecombinationRate scaled = srh_recombination(
            pt.n * s, pt.p * s, boltzmann_equilibrium_product(ni * s), ni * s, tau_n, tau_p);
        REQUIRE(scaled.rate == r.rate * s);
        REQUIRE(scaled.d_dn == r.d_dn);
        REQUIRE(scaled.d_dp == r.d_dp);
    }
}
