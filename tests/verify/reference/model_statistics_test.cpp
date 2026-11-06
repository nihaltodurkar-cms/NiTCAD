// TVR-3 step 6 (DECISIONS Unit 55 15.9, frozen before code): L0 tests of the ionization and equilibrium references
// of TVB-31.17-31.19 (verify/reference/model_statistics.hpp).
//   - d N_ion / d eta against complex step (TC8) under 15.9's partial rule (MI-Q13); N_ion at b = 1 is N s exactly;
//     the documented complete ionization at depth <= 0.
//   - Altermatt: E_dop(N_ref) = E_dop0 / 2 and b(N_b) = 1 / 2, exactly where the power is exact (c = 2).
//   - The equilibrium root: the Boltzmann root against its closed form eta = ln((C + sqrt(C^2 + 4 n_ie^2)) / (2 n_ie))
//     (the second method, 15.3); the Fermi-Dirac root's n - p = C to the charge's rounding; F_1/2's series and
//     quadrature agree at the switch x = -2; the implicit condition sums against condition_sum of re-solved roots.
#include "NiTCAD/verify/reference/model_statistics.hpp"

#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/reference/condition.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>

namespace ref = NiTCAD::verify::reference;
namespace ms = NiTCAD::verify::reference::model_statistics;
using ref::complex_d;
using ref::DD;

namespace {
constexpr double eps = std::numeric_limits<double>::epsilon();

ms::EquilibriumInputs inputs(double ND, double NA, double kD = 0.0, double kA = 0.0, double bD = 1.0,
                             double bA = 1.0) {
    const double ni = 1.0e10, gn = std::log(2.86e19 / ni), gp = std::log(3.10e19 / ni);
    return {DD(ND), DD(NA), DD(ni), DD(gn), DD(gp), DD(kD), DD(kA), DD(2.0), DD(4.0), DD(bD), DD(bA)};
}
}  // namespace

TEST_CASE("model_statistics: ionized density", "[verify][model_statistics]") {
    for (double eta = -50.0; eta <= 50.0; eta += 1.25) {
        for (const double depth : {1e-3, 1.74, 11.6}) {
            for (const double b : {1.0, 0.5, 1e-3}) {
                const std::array<double, 5> x{1e17, eta, depth, 2.0, b};
                auto g = [&](complex_d t) {
                    return ms::ionized(std::array<complex_d, 5>{x[0], t, x[2], x[3], x[4]})[0];
                };
                const auto cs = ref::complex_step(g, eta, 1e-30 * std::fmax(1.0, std::abs(eta)));
                REQUIRE(cs.status == ref::ComplexStepStatus::ok);
                auto f = [](const std::array<DD, 5>& y) { return ms::ionized(y); };
                const std::array<DD, 5> xd{DD(x[0]), DD(x[1]), DD(x[2]), DD(x[3]), DD(x[4])};
                const DD d = f(xd)[1];
                const auto S = ref::condition_sums<2>(f, xd);
                INFO("eta " << eta << " depth " << depth << " b " << b);
                CHECK(std::abs(d.to_double() - cs.derivative) <=
                      std::fmax(1e-13 * std::abs(d.to_double()), 64.0 * eps * S[1]));
                if (b == 1.0) {
                    const DD s = ref::model_recombination::fermi_occupancy(DD(eta) + DD(depth) + log(DD(2.0)))[0];
                    CHECK(f(xd)[0] == DD(1e17) * s);
                }
            }
        }
        // Complete ionization at depth 0 and below.
        const auto c = ms::ionized(std::array<DD, 5>{DD(1e17), DD(eta), DD(0.0), DD(2.0), DD(1.0)});
        CHECK(c[0] == DD(1e17));
        CHECK(c[1] == DD(0.0));
    }
}

TEST_CASE("model_statistics: Altermatt level and bound fraction", "[verify][model_statistics]") {
    const DD E0(0.0455), Nref(2.2e18), Nb(6e18);
    CHECK(ms::altermatt_level(std::array<DD, 4>{E0, Nref, DD(2.0), Nref}) == E0 / DD(2.0));
    CHECK(ms::altermatt_level(std::array<DD, 4>{E0, Nref, DD(2.0), DD(0.0)}) == E0);
    CHECK(ms::altermatt_bound(std::array<DD, 3>{Nb, DD(2.0), Nb}) == DD(0.5));
    CHECK(ms::altermatt_bound(std::array<DD, 3>{Nb, DD(2.3), DD(0.0)}) == DD(1.0));
    // c = 2 at N = 3 N_ref: E0 / 10.
    CHECK(std::abs((ms::altermatt_level(std::array<DD, 4>{E0, Nref, DD(2.0), DD(3.0) * Nref}) * DD(10.0) - E0).hi) <=
          8.0 * ref::u_dd * 0.0455);
}

TEST_CASE("model_statistics: F_1/2 series and quadrature agree at the switch", "[verify][model_statistics]") {
    for (const double x : {-2.0, -2.5, -3.0}) {
        const DD s = ref::fermi_series(ref::FermiOrder::half, DD(x)).value;
        const auto q = ref::fermi_gauss_kronrod(ref::FermiOrder::half, DD(x));
        REQUIRE(q.ok);
        CHECK(std::abs(((s - q.value) / s).hi) <= 1e-26);
    }
    // fermi_half above the switch (tanh-sinh) against Gauss-Kronrod.
    for (const double x : {-1.999, -1.0, 0.0, 3.0, 20.0, 60.0}) {
        const DD t = ms::fermi_half(DD(x));
        const auto q = ref::fermi_gauss_kronrod(ref::FermiOrder::half, DD(x));
        REQUIRE(q.ok);
        INFO("x " << x);
        CHECK(std::abs(((t - q.value) / q.value).hi) <= 1e-26);
    }
}

TEST_CASE("model_statistics: equilibrium roots", "[verify][model_statistics]") {
    const double ni = 1.0e10;
    for (const double C : {0.0, 1.0, 1e5, 1e10, 1e15, 1e19, 1e21, -1e5, -1e15, -1e21}) {
        const ms::EquilibriumInputs x = C >= 0.0 ? inputs(C, 0.0) : inputs(0.0, -C);
        INFO("C " << C);
        // Boltzmann: the closed form.
        const ms::Equilibrium b = ms::equilibrium(x, false, std::asinh(C / (2.0 * ni)));
        REQUIRE(b.ok);
        const DD Cd(std::abs(C)), nid(ni);
        const DD major = (Cd + sqrt(Cd * Cd + DD(4.0) * nid * nid)) / DD(2.0);
        const DD eta = C >= 0.0 ? log(major / nid) : -log(major / nid);
        CHECK(std::abs((b.eta - eta).hi) <= 1e-29 * std::fmax(1.0, std::abs(eta.hi)));
        CHECK(std::abs(((C >= 0.0 ? b.n : b.p) - major).hi) <= 1e-28 * major.hi);
        // Fermi-Dirac: the balance.
        const ms::Equilibrium f = ms::equilibrium(x, true, std::asinh(C / (2.0 * ni)));
        REQUIRE(f.ok);
        CHECK(std::abs((f.n - f.p - DD(C)).hi) <= 1e-27 * (std::abs(C) + ni));
        if (std::abs(C) >= 1e19) CHECK(std::abs(f.eta.hi) > std::abs(b.eta.hi));  // degeneracy lifts eta
    }
    // The implicit condition sums against re-solved roots (Boltzmann; n, p, eta).
    for (const double C : {1e5, 1e15, -1e17}) {
        const ms::EquilibriumInputs x = C >= 0.0 ? inputs(C, 0.0, 1.74, 1.74) : inputs(0.0, -C, 1.74, 1.74);
        const ms::Equilibrium e = ms::equilibrium(x, false, std::asinh(C / (2.0 * ni)));
        REQUIRE(e.ok);
        auto solve = [&](const ms::EquilibriumInputs& y) {
            const ms::Equilibrium r = ms::equilibrium(y, false, std::asinh(C / (2.0 * ni)));
            return std::array<DD, 3>{r.n, r.p, r.eta};
        };
        const auto S = ref::condition_sums<3>(solve, x);
        for (int k = 0; k < 3; ++k) {
            INFO("C " << C << " k " << k << ": implicit " << e.S[k] << " re-solved " << S[k]);
            CHECK(std::abs(e.S[k] - S[k]) <= 1e-6 * S[k]);
        }
    }
}
