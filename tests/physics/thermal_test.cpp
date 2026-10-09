// Thermal data, the Fermi-Dirac integrals F_0 and F_1, the thermal diffusion factor and the
// temperature derivatives of the material models (ARCHITECTURE.md section 11, Unit 23).
//
// References in double-double arithmetic (references.hpp): F_0 = ln(1 + e^eta) directly; F_1 by
// its alternating series below eta = -1 and by Gauss-Legendre quadrature of its defining integral
// above, the two checked against each other where both converge. The temperature derivatives are
// checked against fourth-order central differences of the material formulas evaluated in DD (step
// 1e-3 K: truncation about (1e-3 / T)^4, rounding about 1e-32 / 1e-3, both far below double
// precision), so they test the closed forms of thermal.hpp against the formulas they differentiate.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "references.hpp"

using namespace NiTCAD::physics;
using reference::DD;

namespace {

constexpr double eps = std::numeric_limits<double>::epsilon();

double relative(double a, double b) { return std::abs(a - b) / std::abs(b); }

// d f / dT at T by the fourth-order central difference in DD.
double slope(const std::function<DD(const DD&)>& f, double T) {
    const DD h = 1e-3;
    const DD t(T);
    const DD d = (f(t - DD(2.0) * h) - DD(8.0) * f(t - h) + DD(8.0) * f(t + h) -
                  f(t + DD(2.0) * h)) /
                 (DD(12.0) * h);
    return d.hi + d.lo;
}

DD kT_eV(const DD& T) {
    return DD(NiTCAD::base::k_B_J_per_K) * T / DD(NiTCAD::base::q_C);
}
DD band_gap(const SemiconductorParameters& p, const DD& T) {
    return DD(p.Eg0_eV) - DD(p.varshni_alpha_eV_per_K) * T * T / (T + DD(p.varshni_beta_K));
}
DD dos_factor(const DD& T) {
    const DD r = T / DD(300.0);
    return r * sqrt(r);
}
DD log_intrinsic(const SemiconductorParameters& p, const DD& T, double narrowing_eV) {
    return log(sqrt(DD(p.Nc300) * DD(p.Nv300)) * dos_factor(T)) -
           (band_gap(p, T) - DD(narrowing_eV)) / (DD(2.0) * kT_eV(T));
}
DD tanh_dd(const DD& a) {
    const DD e = exp(DD(2.0) * a);
    return (e - DD(1.0)) / (e + DD(1.0));
}

}  // namespace

TEST_CASE("F_0 and F_1: the reference methods agree with each other") {
    for (const double eta : {-12.5, -5.0, -2.000001, -1.2345, -1.0001}) {
        CAPTURE(eta);
        const DD s = reference::fermi_one_series(eta);
        const DD q = reference::fermi_one_quadrature(eta);
        REQUIRE(std::abs((s - q).value()) <= 1e-26 * std::abs(s.value()));
    }
    // F_1(0) = pi^2 / 12.
    const DD q = reference::fermi_one_quadrature(0.0);
    const DD exact = reference::pi() * reference::pi() / DD(12.0);
    REQUIRE(std::abs((q - exact).value()) <= 1e-28);
}

TEST_CASE("F_-1, F_0 and F_1 match double-double references") {
    for (const double x : {-700.0, -40.0, -12.5, -3.3, -1.0000001, -1.0, -0.999, -0.5, -1e-9,
                           0.0, 1e-9, 0.25, 0.7777, 1.0, 1.0000001, 1.5, 2.0, 3.14159, 7.3,
                           19.99, 40.5, 100.0, 1000.0}) {
        CAPTURE(x);
        const FermiOrdersZeroOne f = fermi_orders_zero_one(x);
        const double f0 = reference::fermi_zero(x).value();
        const double f1 = reference::fermi_one(x).value();
        const double fm1 = 1.0 / (1.0 + std::exp(-x));
        CHECK(relative(f.zero, f0) <= 4.0 * eps);
        CHECK(relative(f.one, f1) <= 4.0 * eps);
        CHECK(relative(f.minus_one, fm1) <= 2.0 * eps);
    }
    // Underflow with e^x, and NaN propagation.
    REQUIRE(fermi_orders_zero_one(-800.0).one == 0.0);
    REQUIRE(std::isnan(fermi_orders_zero_one(std::numeric_limits<double>::quiet_NaN()).one));
}

TEST_CASE("Fermi-Dirac thermal diffusion factor and its slope") {
    // h = 2 F_1 / F_0 - 3/2 and dh/dx = 2 - 2 F_1 F_-1 / F_0^2, against DD values of the same.
    for (const double x : {-800.0, -40.0, -12.5, -3.3, -1.0, -0.5, 0.0, 0.25, 1.0, 1.5, 3.14159,
                           7.3, 19.99, 40.5, 100.0}) {
        CAPTURE(x);
        const ThermalDiffusionFactor h = fermi_dirac_thermal_diffusion(x);
        // Below -40 F_0 and F_1 in DD would need their own e^x scaling; the Boltzmann limit.
        if (x < -40.0) {
            REQUIRE(h.value == 0.5);
            REQUIRE(std::abs(h.d_x) <= 1e-300);
            continue;
        }
        const DD f0 = reference::fermi_zero(x), f1 = reference::fermi_one(x);
        const DD fm1 = DD(1.0) / (DD(1.0) + exp(DD(-x)));
        const DD value = DD(2.0) * f1 / f0 - DD(1.5);
        const DD d = DD(2.0) - DD(2.0) * f1 * fm1 / (f0 * f0);
        CHECK(std::abs(h.value - value.value()) <= 8.0 * eps * std::max(1.0, std::abs(value.value())));
        // The slope cancels to its own size near the Boltzmann limit: absolute.
        CHECK(std::abs(h.d_x - d.value()) <= 8.0 * eps);
    }
    // Limits: 1/2 (Boltzmann, r = -1/2) and, for large x, F_1 -> x^2 / 2 + pi^2 / 6 and F_0 -> x:
    // h -> x - 3/2 + pi^2 / (3 x).
    REQUIRE(std::abs(fermi_dirac_thermal_diffusion(-40.0).value - 0.5) <= 1e-15);
    const double pi2 = std::numbers::pi * std::numbers::pi;
    REQUIRE(std::abs(fermi_dirac_thermal_diffusion(1000.0).value - (1000.0 - 1.5 + pi2 / 3000.0)) <=
            1e-12);
}

TEST_CASE("Thermal data of the material sets and its validation") {
    const ThermalParameters& si = silicon_parameters.thermal;
    REQUIRE(has_thermal_data(si));
    REQUIRE(si.conductivity_W_cmK == 1.48);  // legacy materials.kappa_th
    REQUIRE(si.conductivity_exponent == 1.33);
    REQUIRE(std::abs(si.heat_capacity_J_cm3K - 2.329 * 0.7) <= 1e-12);
    REQUIRE(has_thermal_data(silicon_carbide_4h_parameters.thermal));
    REQUIRE(silicon_carbide_4h_parameters.thermal.conductivity_W_cmK == 3.7);
    // No invented values: the other sets carry none.
    for (const SemiconductorParameters* p :
         {&germanium_parameters, &gallium_arsenide_parameters,
          &indium_gallium_arsenide_parameters}) {
        REQUIRE_FALSE(has_thermal_data(p->thermal));
    }
    REQUIRE_FALSE(has_thermal_data(algaas_parameters(0.3)->thermal));
    REQUIRE(silicon_dioxide_parameters.thermal.conductivity_W_cmK == 0.014);
    REQUIRE(silicon_dioxide_parameters.thermal.conductivity_exponent == 0.0);
    REQUIRE(silicon_dioxide().parameters().thermal.heat_capacity_J_cm3K == 1.628);

    // kappa(300 K) is the 300 K value; its slope against DD differences.
    REQUIRE(thermal_conductivity(si, 300.0).value == 1.48);
    for (const double T : {250.0, 300.0, 412.5, 600.0}) {
        CAPTURE(T);
        const ValueSlope k = thermal_conductivity(si, T);
        const auto f = [&](const DD& t) {
            return DD(si.conductivity_W_cmK) * pow(t / DD(300.0), DD(-si.conductivity_exponent));
        };
        CHECK(relative(k.value, f(DD(T)).value()) <= 4.0 * eps);
        CHECK(relative(k.d_dT, slope(f, T)) <= 1e-13);
    }

    // Validation: finite values, conductivity >= 0, heat capacity > 0 with a conductivity.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto rejects = [](ThermalParameters t) { return !check_thermal_parameters(t, "x"); };
    REQUIRE(check_thermal_parameters({}, "x"));  // no data
    REQUIRE(check_thermal_parameters(si, "x"));
    ThermalParameters bad = si;
    bad.conductivity_W_cmK = -1.0;
    REQUIRE(rejects(bad));
    bad = si;
    bad.heat_capacity_J_cm3K = 0.0;
    REQUIRE(rejects(bad));
    bad = si;
    bad.conductivity_exponent = nan;
    REQUIRE(rejects(bad));
    bad = si;
    bad.thermopower_exponent_p = std::numeric_limits<double>::infinity();
    REQUIRE(rejects(bad));
    SemiconductorParameters p = silicon_parameters;
    p.thermal.heat_capacity_J_cm3K = -1.0;
    const auto created = Semiconductor::create(p);
    REQUIRE_FALSE(created);
    REQUIRE(created.error().message.find("heat_capacity_J_cm3K") != std::string::npos);
    InsulatorParameters ox = silicon_dioxide_parameters;
    ox.thermal.conductivity_W_cmK = nan;
    REQUIRE_FALSE(Insulator::create(ox));
}

TEST_CASE("Temperature derivatives of the material models") {
    const Semiconductor si = silicon();
    const Semiconductor sic = *Semiconductor::create(silicon_carbide_4h_parameters);
    for (const Semiconductor* m : {&si, &sic}) {
        const SemiconductorParameters& p = m->parameters();
        for (const double T : {200.0, 300.0, 377.7, 520.0}) {
            CAPTURE(T);
            CHECK(relative(band_gap_slope_eV_per_K(*m, T),
                           slope([&](const DD& t) { return band_gap(p, t); }, T)) <= 1e-12);
            CHECK(relative(intrinsic_density_log_slope(*m, T),
                           slope([&](const DD& t) { return log_intrinsic(p, t, 0.0); }, T)) <=
                  1e-12);
            const double narrowing = bandgap_narrowing_eV(*m, 3e19);
            CHECK(relative(intrinsic_density_log_slope(*m, T, narrowing),
                           slope([&](const DD& t) { return log_intrinsic(p, t, narrowing); }, T)) <=
                  1e-12);
            for (const Carrier c : {Carrier::electron, Carrier::hole}) {
                const CaugheyThomasParameters& ct =
                    c == Carrier::electron ? p.electron_mobility : p.hole_mobility;
                for (const double N : {0.0, 1e15, 3e17, 1e20}) {
                    const auto mu = [&](const DD& t) {
                        const DD mu_max = DD(ct.mu_max) * pow(t / DD(300.0), DD(ct.T_exponent));
                        const DD x = N == 0.0 ? DD(0.0) : pow(DD(N) / DD(ct.N_ref), DD(ct.alpha));
                        return DD(ct.mu_min) + (mu_max - DD(ct.mu_min)) / (DD(1.0) + x);
                    };
                    CHECK(relative(caughey_thomas_mobility_slope(*m, c, N, T), slope(mu, T)) <=
                          1e-12);
                }
                // v = k T / h (2 / N)^(1/3) with N = N300 (T / 300)^(3/2) (no Richardson constant).
                const double N300 = c == Carrier::electron ? p.Nc300 : p.Nv300;
                const auto v = [&](const DD& t) {
                    const DD h = DD(2.0) * reference::pi() * DD(NiTCAD::base::hbar_J_s);
                    return DD(NiTCAD::base::k_B_J_per_K) * t / h *
                           exp(log(DD(2.0) / (DD(N300) * dos_factor(t))) / DD(3.0));
                };
                CHECK(relative(emission_velocity_cm_s(*m, c, T), v(DD(T)).value()) <= 1e-14);
                CHECK(relative(emission_velocity_slope(*m, c, T), slope(v, T)) <= 1e-13);
            }
        }
    }
    // The impact-ionization temperature factor gamma = tanh(h / 300) / tanh(h / T).
    const double phonon = silicon_parameters.impact_ionization.phonon_energy_eV;
    const double half = 0.5 * phonon / NiTCAD::base::k_B_eV_per_K;
    for (const double T : {150.0, 300.0, 450.0, 700.0}) {
        CAPTURE(T);
        const auto gamma = [&](const DD& t) {
            return tanh_dd(DD(half) / DD(300.0)) / tanh_dd(DD(half) / t);
        };
        CHECK(relative(impact_ionization_temperature_factor_slope(phonon, T), slope(gamma, T)) <=
              1e-12);
    }
    REQUIRE(impact_ionization_temperature_factor_slope(0.0, 400.0) == 0.0);
}
