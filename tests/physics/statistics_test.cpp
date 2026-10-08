// Boltzmann (ARCHITECTURE.md section 11, Unit 5; R4) and Fermi-Dirac (Unit 14) carrier statistics.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/fermi_dirac.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "references.hpp"

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

// Unit 14: Fermi-Dirac statistics in the n_ie gauge.

namespace {

// Silicon at 300 K: n_i and g = ln(N / n_i) for both bands.
struct Bands {
    double ni, gn, gp;
};

Bands silicon_bands() {
    const Semiconductor si = silicon();
    const double ni = intrinsic_density(si, 300.0);
    return {ni, std::log(conduction_band_dos(si, 300.0) / ni),
            std::log(valence_band_dos(si, 300.0) / ni)};
}

}  // namespace

TEST_CASE("statistics: the Fermi-Dirac density is N F_1/2(eta - g), Boltzmann when nondegenerate") {
    const Bands b = silicon_bands();
    const double Nc = b.ni * std::exp(b.gn);
    for (const double eta : {-30.0, -5.0, 0.0, 10.0, 20.0, 21.0, 24.0, 30.0, 60.0}) {
        CAPTURE(eta);
        const DensityResult r = fermi_dirac_density(b.ni, b.gn, eta);
        REQUIRE(close(r.density, Nc * fermi_half(eta - b.gn).value, 1e-13));
        // Below the band edge the deviation from n_i e^eta is the series term e^(eta-g) / 2^(3/2),
        // and nothing more than rounding.
        const double boltzmann = boltzmann_density(b.ni, eta).density;
        const double x = std::exp(eta - b.gn);
        if (eta - b.gn < -10.0) {
            REQUIRE(std::abs(r.density / boltzmann - 1.0) <= x / std::pow(2.0, 1.5) + 4.0 * eps);
        }
        REQUIRE(r.density <= boltzmann);
        const double up = eta + 1e-4, down = eta - 1e-4;
        const double fd = (fermi_dirac_density(b.ni, b.gn, up).density -
                           fermi_dirac_density(b.ni, b.gn, down).density) /
                          (up - down);
        REQUIRE(close(r.d_eta, fd, 1e-8));
    }
}

TEST_CASE("statistics: the degeneracy factor of a density, and its derivative") {
    const Bands b = silicon_bands();
    const double Nc = b.ni * std::exp(b.gn);
    // From the density back to ln gamma at its reduced energy, through the series branch
    // (density / Nc below 1e-6) and the inversion.
    for (const double eta : {-10.0, 5.0, 7.8, 7.9, 15.0, 21.7, 25.0, 40.0, 70.0}) {
        CAPTURE(eta);
        const double n = fermi_dirac_density(b.ni, b.gn, eta).density;
        const Degeneracy d = fermi_dirac_degeneracy(b.ni, b.gn, n);
        REQUIRE(std::abs(d.log_gamma - log_degeneracy(eta - b.gn).value) <=
                1e-14 * std::max(1.0, std::abs(d.log_gamma)) + 1e-18);
        const double h = 1e-6 * n;
        const double fd = (fermi_dirac_degeneracy(b.ni, b.gn, n + h).log_gamma -
                           fermi_dirac_degeneracy(b.ni, b.gn, n - h).log_gamma) /
                          (2.0 * h);
        REQUIRE(close(d.d_density, fd, 1e-7));
        REQUIRE(d.d_density <= 0.0);
    }
    // The series and the inversion meet at density / Nc = 1e-6.
    const Degeneracy below = fermi_dirac_degeneracy(b.ni, b.gn, Nc * std::nextafter(1e-6, 0.0));
    const Degeneracy above = fermi_dirac_degeneracy(b.ni, b.gn, Nc * std::nextafter(1e-6, 1.0));
    REQUIRE(std::abs(below.log_gamma - above.log_gamma) <= 1e-18);
    REQUIRE(close(below.d_density, above.d_density, 1e-6));
    // A finite-difference probe may take a minority density through zero: still finite, and
    // continuous there.
    const Degeneracy zero = fermi_dirac_degeneracy(b.ni, b.gn, 0.0);
    REQUIRE(zero.log_gamma == 0.0);
    REQUIRE(close(zero.d_density, -1.0 / (std::pow(2.0, 1.5) * Nc), 1e-15));
    const Degeneracy negative = fermi_dirac_degeneracy(b.ni, b.gn, -1e3);
    REQUIRE(negative.log_gamma > 0.0);
    REQUIRE(std::isfinite(negative.d_density));
}

TEST_CASE("statistics: generalized mass action under Fermi-Dirac statistics") {
    // With a common Fermi level, n p = n_ie^2 gamma_n gamma_p (legacy G4(c)); the product's
    // partials match finite differences of the product.
    const Bands b = silicon_bands();
    for (const double eta : {-24.0, -10.0, 0.0, 10.0, 20.0, 23.0, 26.0}) {
        CAPTURE(eta);
        const double n = fermi_dirac_density(b.ni, b.gn, eta).density;
        const double p = fermi_dirac_density(b.ni, b.gp, -eta).density;
        const Degeneracy dn = fermi_dirac_degeneracy(b.ni, b.gn, n);
        const Degeneracy dp = fermi_dirac_degeneracy(b.ni, b.gp, p);
        const EquilibriumProduct e = fermi_dirac_equilibrium_product(b.ni, dn, dp);
        REQUIRE(close(e.value, n * p, 1e-13));
        const double hn = 1e-6 * n, hp = 1e-6 * p;
        const auto product = [&](double nn, double pp) {
            return fermi_dirac_equilibrium_product(b.ni, fermi_dirac_degeneracy(b.ni, b.gn, nn),
                                                   fermi_dirac_degeneracy(b.ni, b.gp, pp))
                .value;
        };
        // Only where a 1e-6 step moves the product by more than rounding (judged on the
        // difference, not on the partial under test): deep in the Boltzmann regime the partials
        // are below what a difference resolves. Degenerate electrons (eta >= 20 here) always are.
        const double fd_n = (product(n + hn, p) - product(n - hn, p)) / (2.0 * hn);
        const double fd_p = (product(n, p + hp) - product(n, p - hp)) / (2.0 * hp);
        if (std::abs(fd_n * hn) > 1e-9 * e.value) REQUIRE(close(e.d_dn, fd_n, 1e-6));
        if (std::abs(fd_p * hp) > 1e-9 * e.value) REQUIRE(close(e.d_dp, fd_p, 1e-6));
        if (eta >= 20.0) REQUIRE(std::abs(fd_n * hn) > 1e-9 * e.value);
    }
    // Nondegenerate: the Boltzmann n_ie^2 to the series term.
    const Degeneracy light = fermi_dirac_degeneracy(b.ni, b.gn, 1e10);
    REQUIRE(close(fermi_dirac_equilibrium_product(b.ni, light, light).value, b.ni * b.ni, 1e-9));
}

TEST_CASE("statistics: Fermi-Dirac neutral equilibrium against double-double roots") {
    // Silicon, 300 K, no band-gap narrowing: the root of Nc F(eta - g_n) - Nv F(-eta - g_p) = C by
    // bisection in double-double arithmetic with the reference F_{1/2} (references.hpp), with this
    // code's n_i, Nc and Nv.
    const Bands b = silicon_bands();
    for (const double C : {1e20, -1e20, 1e17, 1e19, -3e19, 0.0}) {
        CAPTURE(C);
        reference::Dopants d;
        (C > 0.0 ? d.donors : d.acceptors) = std::abs(C);
        const reference::Neutral r = reference::neutral_equilibrium(d, b.ni, b.gn, b.gp, true);
        const NeutralEquilibrium e = fermi_dirac_neutral_equilibrium(C, b.ni, b.gn, b.gp);
        CAPTURE(e.eta, r.eta);
        REQUIRE(std::abs(e.eta - r.eta) <= 1e-13 * std::max(1.0, std::abs(r.eta)));
        REQUIRE(close(e.n, r.n.value(), 1e-12));
        REQUIRE(close(e.p, r.p.value(), 1e-12));
        REQUIRE(std::abs((e.n - e.p) - C) <= 4.0 * eps * std::max(e.n, e.p));
        // Consistent with the density functions at its own eta.
        REQUIRE(close(e.n, fermi_dirac_density(b.ni, b.gn, e.eta).density, 1e-13));
        REQUIRE(close(e.p, fermi_dirac_density(b.ni, b.gp, -e.eta).density, 1e-13));
    }
    // 1e20 cm^-3 is degenerate: the Fermi level 2.4 kT above the band edge (legacy G7(a): > 2).
    REQUIRE(fermi_dirac_neutral_equilibrium(1e20, b.ni, b.gn, b.gp).eta - b.gn > 2.0);
}

TEST_CASE("statistics: Fermi-Dirac neutral equilibrium is scale-free and Boltzmann when light") {
    const Bands b = silicon_bands();
    for (const double C : {1e12, -1e15, 3e17, -5e19, 2e21}) {
        CAPTURE(C);
        const NeutralEquilibrium e = fermi_dirac_neutral_equilibrium(C, b.ni, b.gn, b.gp);
        const double s = 1e-18;
        const NeutralEquilibrium scaled =
            fermi_dirac_neutral_equilibrium(C * s, b.ni * s, b.gn, b.gp);
        REQUIRE(std::abs(scaled.eta - e.eta) <= 1e-14 * std::abs(e.eta));
        REQUIRE(close(scaled.n, e.n * s, 1e-13));
        REQUIRE(close(scaled.p, e.p * s, 1e-13));
        // FD needs a higher Fermi level than Boltzmann for the same majority density.
        REQUIRE(std::abs(e.eta) >= std::abs(boltzmann_neutral_equilibrium(C, b.ni).eta));
    }
    // At 1e12 the series deviation of the majority is e^(eta-g) / 2^(3/2) ~ 1e-8 relative.
    const NeutralEquilibrium light = fermi_dirac_neutral_equilibrium(1e12, b.ni, b.gn, b.gp);
    const NeutralEquilibrium classic = boltzmann_neutral_equilibrium(1e12, b.ni);
    REQUIRE(std::abs(light.eta - classic.eta) <= 1e-7);
}
