// The density-dependent dopant levels of Altermatt et al. (Unit 37, physics/ionization.hpp): the
// ionized density against eqs (1), (3)-(5) and Table I as recorded in DECISIONS.md, evaluated here
// independently (II2), and the non-ionized fraction of neutral silicon at 300 K (II3).
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/ionization.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;

namespace {

// Written out from DECISIONS.md (Unit 37), not from the implementation's constants.
struct Set {
    double E0_eV, N_ref, c, N_b, d, g;  // g as the paper's Table I (the reciprocal degeneracy)
};
constexpr Set phosphorus{0.0455, 2.2e18, 2.0, 6e18, 2.3, 0.5};
constexpr Set boron{0.04439, 1.3e18, 1.4, 4.5e18, 2.4, 0.25};

// N+ / N = 1 - b(N) x / (x + g), x = n / n1 (eq (3)), at a given x.
double ionized_fraction(const Set& s, double N, double x) {
    const double b = 1.0 / (1.0 + std::pow(N / s.N_b, s.d));
    return 1.0 - b * x / (x + s.g);
}

double level_eV(const Set& s, double N) { return s.E0_eV / (1.0 + std::pow(N / s.N_ref, s.c)); }

}  // namespace

TEST_CASE("ionization density: Altermatt's levels as published (II2)", "[physics][ionization]") {
    for (const auto& [set, model] :
         {std::pair{phosphorus, altermatt_phosphorus}, std::pair{boron, altermatt_boron}}) {
        const double degeneracy = 1.0 / set.g;  // the legacy g_D = 2, g_A = 4
        for (int e = 14; e <= 21; ++e) {
            const double N = std::pow(10.0, e);
            CAPTURE(set.E0_eV, N);
            CHECK(std::abs(dopant_level_eV(model, N) / level_eV(set, N) - 1.0) < 1e-12);
            const double depth = level_eV(set, N) / NiTCAD::base::thermal_voltage(300.0);
            const double b = bound_fraction(model, N);
            for (const double x : {1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e3}) {
                // x = n / n1 = e^(eta + depth) (Boltzmann), so eta = ln x - depth.
                const double eta = std::log(x) - depth;
                const IonizedDensity r = ionized_density(N, eta, depth, degeneracy, b);
                CAPTURE(x, r.value / N, ionized_fraction(set, N, x));
                CHECK(std::abs(r.value / (N * ionized_fraction(set, N, x)) - 1.0) < 1e-12);
                // The slope: the paper's formula differentiated, d/deta = -N b g x / (x + g)^2
                // (dx/deta = x), everywhere; a central difference of the implementation wherever
                // it resolves the slope (its rounding, about eps N / h, below 1e-8 of it: where
                // nearly every atom is ionized the value is N (1 - 1e-6) and a difference of it is
                // rounding).
                const double bp = 1.0 / (1.0 + std::pow(N / set.N_b, set.d));
                const double analytic = -N * bp * set.g * x / ((x + set.g) * (x + set.g));
                CHECK(std::abs(r.d_eta / analytic - 1.0) < 1e-10);
                const double h = 1e-5;
                if (2.2e-16 * N / h < 1e-8 * std::abs(r.d_eta)) {
                    const double fd = (ionized_density(N, eta + h, depth, degeneracy, b).value -
                                       ionized_density(N, eta - h, depth, degeneracy, b).value) /
                                      (2.0 * h);
                    CHECK(std::abs(r.d_eta / fd - 1.0) < 1e-6);
                }
            }
        }
    }
    // b = 1 is the fixed-level kernel, bit for bit.
    const IonizedDensity a = ionized_density(1e17, -2.0, 1.7, 2.0);
    const IonizedDensity b = ionized_density(1e17, -2.0, 1.7, 2.0, 1.0);
    CHECK(a.value == b.value);
    CHECK(a.d_eta == b.d_eta);
}

TEST_CASE("ionization density: non-ionized fraction of neutral silicon at 300 K (II3)",
          "[physics][ionization]") {
    // The paper: up to 25% of P and B atoms non-ionized near the Mott density (3.74e18 cm^-3) at
    // room temperature, none to speak of far above it. Neutral uniform silicon (the published band
    // set, Boltzmann statistics, no band-gap narrowing), solved for the Fermi level.
    const Semiconductor si = silicon();
    const double T = 300.0, VT = NiTCAD::base::thermal_voltage(T);
    const double ni = intrinsic_density(si, T);
    const double log_nc = std::log(conduction_band_dos(si, T) / ni);
    const double log_nv = std::log(valence_band_dos(si, T) / ni);
    const auto non_ionized = [&](bool donor, double N) {
        DopantLevels l{0.0, 0.0, 2.0, 4.0};
        const DensityDependentLevel& set = donor ? altermatt_phosphorus : altermatt_boron;
        (donor ? l.donor_kT : l.acceptor_kT) = dopant_level_eV(set, N) / VT;
        (donor ? l.donor_bound : l.acceptor_bound) = bound_fraction(set, N);
        const NeutralEquilibrium q = ionized_neutral_equilibrium(
            donor ? N / ni : 0.0, donor ? 0.0 : N / ni, 1.0, log_nc, log_nv, l, false);
        // Neutral: the ionized dopant equals the majority density less the minority one.
        return 1.0 - (donor ? q.n - q.p : q.p - q.n) * ni / N;
    };
    for (const bool donor : {true, false}) {
        double worst = 0.0, at = 0.0;
        for (double lg = 17.0; lg <= 20.0001; lg += 0.02) {
            const double N = std::pow(10.0, lg);
            const double f = non_ionized(donor, N);
            if (f > worst) {
                worst = f;
                at = N;
            }
        }
        const double high = non_ionized(donor, 1e20), low = non_ionized(donor, 1e15);
        std::printf("%s: largest non-ionized fraction %.3f at %.3g cm^-3; %.4f at 1e20, %.5f at "
                    "1e15\n",
                    donor ? "phosphorus" : "boron", worst, at, high, low);
        CAPTURE(donor, worst, at, high, low);
        CHECK(worst >= 0.20);
        CHECK(worst <= 0.30);
        CHECK(at >= 3.74e18 / 3.0);
        CHECK(at <= 3.74e18 * 3.0);
        CHECK(high < 0.03);
        CHECK(low < 0.01);
    }
}
