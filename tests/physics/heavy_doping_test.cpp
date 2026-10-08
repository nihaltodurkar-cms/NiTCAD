// Slotboom band-gap narrowing and Auger recombination (ARCHITECTURE.md section 11, Unit 11).
// References computed here in double-double arithmetic (references.hpp) from the legacy formulas,
// independently of this code.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "references.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"

using namespace NiTCAD::physics;
using NiTCAD::base::ErrorCode;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

constexpr double Cn = 2.8e-31, Cp = 9.9e-32;

}  // namespace

TEST_CASE("bgn: Slotboom narrowing of silicon") {
    const Semiconductor si = silicon();
    REQUIRE(bandgap_narrowing_eV(si, 0.0) == 0.0);
    REQUIRE(bandgap_narrowing_eV(si, 1e17) == 0.0);
    REQUIRE(bandgap_narrowing_eV(si, 1.3e17) == 0.0);  // zero at N0, so no step there
    REQUIRE(bandgap_narrowing_eV(si, 1.3e17 * (1 + 1e-9)) < 1e-11);
    for (const double N : {1e18, 1e19, 1e20}) {
        CAPTURE(N);
        REQUIRE(close(bandgap_narrowing_eV(si, N),
                      reference::slotboom_narrowing(silicon_parameters, N).value(), 1e-13));
    }
    // Legacy test_slotboom_bgn_positive_and_monotonic.
    const double lo = bandgap_narrowing_eV(si, 1e18), hi = bandgap_narrowing_eV(si, 1e20);
    REQUIRE((hi > lo && lo >= 0.0));
    double previous = 0.0;
    for (double N = 1.3e17; N < 1e21; N *= 1.3) {
        const double d = bandgap_narrowing_eV(si, N);
        REQUIRE(d >= previous);
        previous = d;
    }
}

TEST_CASE("bgn: effective intrinsic density") {
    const Semiconductor si = silicon();
    REQUIRE(effective_intrinsic_density(si, 1e16, 300.0) == intrinsic_density(si, 300.0));
    REQUIRE(close(effective_intrinsic_density(si, 1e18, 300.0), 17033941184.25466429, 1e-12));
    REQUIRE(close(effective_intrinsic_density(si, 1e19, 300.0), 31289337616.097783699, 1e-12));
    REQUIRE(close(effective_intrinsic_density(si, 1e20, 300.0), 57800322297.905986582, 1e-12));
}

TEST_CASE("auger: vanishes at equilibrium and grows as the cube of the injection") {
    const double ni = effective_intrinsic_density(silicon(), 1e15, 300.0);
    // n p = n_ie^2 exactly (factors of two).
    static_assert(auger_recombination(4.0, 0.25, boltzmann_equilibrium_product(1.0), Cn, Cp).rate ==
                  0.0);
    REQUIRE(auger_recombination(2.0 * ni, 0.5 * ni, boltzmann_equilibrium_product(ni), Cn, Cp)
                .rate == 0.0);
    // Far above equilibrium R = (Cn n + Cp p) n p, so doubling both densities multiplies it by 8.
    const auto R = [&](double n, double p) {
        return auger_recombination(n, p, boltzmann_equilibrium_product(ni), Cn, Cp).rate;
    };
    REQUIRE(close(R(2e18, 2e18), 8.0 * R(1e18, 1e18), 1e-12));
    REQUIRE(close(R(1e18, 1e18), (Cn + Cp) * 1e36 * 1e18, 1e-12));
    // Legacy test_auger_grows_quadratically_with_carrier_density: SRH + Auger, doubled densities
    // give more than twice the rate.
    const auto total = [&](double n, double p) {
        const EquilibriumProduct e = boltzmann_equilibrium_product(ni);
        return srh_recombination(n, p, e, ni, 1e-5, 3e-6).rate + R(n, p);
    };
    REQUIRE(total(2e18, 2e18) > 2.0 * total(1e18, 1e18));
}

TEST_CASE("auger: partials match finite differences, including the equilibrium-product chain") {
    const double ni = 1.07e10, K = 1e19;
    const auto product = [&](double n, double p) {
        return EquilibriumProduct{ni * ni * (1.0 + n / K) * (1.0 + p / K),
                                  ni * ni * (1.0 + p / K) / K, ni * ni * (1.0 + n / K) / K};
    };
    struct Point {
        double n, p;
    };
    const Point points[] = {{1e18, 1e18}, {1e19, 1e5}, {1e3, 1e17}, {1e10, 2e9}};
    for (const Point pt : points) {
        for (const bool chain : {false, true}) {
            const auto rate = [&](double n, double p) {
                const EquilibriumProduct e =
                    chain ? product(n, p) : boltzmann_equilibrium_product(ni);
                return auger_recombination(n, p, e, Cn, Cp).rate;
            };
            const EquilibriumProduct e =
                chain ? product(pt.n, pt.p) : boltzmann_equilibrium_product(ni);
            const RecombinationRate r = auger_recombination(pt.n, pt.p, e, Cn, Cp);
            const double hn = 1e-5 * pt.n, hp = 1e-5 * pt.p;
            const double fd_n = (rate(pt.n + hn, pt.p) - rate(pt.n - hn, pt.p)) /
                                ((pt.n + hn) - (pt.n - hn));
            const double fd_p = (rate(pt.n, pt.p + hp) - rate(pt.n, pt.p - hp)) /
                                ((pt.p + hp) - (pt.p - hp));
            const double scale = std::max(std::abs(r.d_dn), std::abs(r.d_dp));
            CAPTURE(pt.n, pt.p, chain);
            REQUIRE(std::abs(fd_n - r.d_dn) <= 1e-8 * scale);
            REQUIRE(std::abs(fd_p - r.d_dp) <= 1e-8 * scale);
        }
    }
}

TEST_CASE("bgn and auger: parameters are validated") {
    using P = SemiconductorParameters;
    struct Case {
        std::string name;
        std::function<void(P&)> spoil;
    };
    const std::vector<Case> cases{
        {"auger.Cn", [](P& p) { p.auger.Cn = -1e-31; }},
        {"auger.Cp", [](P& p) { p.auger.Cp = std::numeric_limits<double>::infinity(); }},
        {"bandgap_narrowing.E0_eV", [](P& p) { p.bandgap_narrowing.E0_eV = -1e-3; }},
        {"bandgap_narrowing.N0", [](P& p) { p.bandgap_narrowing.N0 = 0.0; }},
    };
    for (const Case& c : cases) {
        CAPTURE(c.name);
        P p = silicon_parameters;
        c.spoil(p);
        const auto m = Semiconductor::create(p);
        REQUIRE_FALSE(m.has_value());
        REQUIRE(m.error().code == ErrorCode::invalid_input);
        REQUIRE(m.error().message.find(c.name) != std::string::npos);
    }
    P off = silicon_parameters;  // switching the models off through the parameters is allowed
    off.auger = {.Cn = 0.0, .Cp = 0.0};
    off.bandgap_narrowing.E0_eV = 0.0;
    REQUIRE(Semiconductor::create(off).has_value());
}
