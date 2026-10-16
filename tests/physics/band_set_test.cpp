// The default silicon band set (Unit 35, DECISIONS.md S1-S2): Couderc, Amara and Lemiti, J. Appl.
// Phys. 115, 093705 (2014), against the formulas and coefficients as recorded in DECISIONS.md,
// evaluated here independently (DC2), n_i against the paper's eq (15) and Misiakos's measurement
// (DC3), the temperature slopes the electrothermal model uses (DC6), and the validation of the
// band forms.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"

using namespace NiTCAD::physics;

namespace {

double relative(double a, double b) { return std::abs(a / b - 1.0); }

// Written out from DECISIONS.md (Unit 35), not from the parameter struct.
double gap(double T) {
    const double E0 = 1.17, alpha = 3.23e-4, Theta = 446.0, D = 0.51;
    const double chi = 2.0 * T / Theta;
    const double pi2 = std::numbers::pi * std::numbers::pi;
    const double root = std::pow(1.0 + pi2 * chi * chi / (3.0 * (1.0 + D * D)) +
                                     (3.0 * D * D - 1.0) * std::pow(chi, 3) / 4.0 +
                                     8.0 * std::pow(chi, 4) / 3.0 + std::pow(chi, 6),
                                 1.0 / 6.0);
    return E0 - alpha * Theta *
                    ((1.0 - 3.0 * D * D) / (std::exp(Theta / T) - 1.0) +
                     1.5 * D * D * (root - 1.0));
}
double nc(double T) {
    const double m = -4.609e-10 * T * T * T + 6.753e-7 * T * T - 1.312e-5 * T + 1.094;
    return 4.83e15 * m * std::pow(T, 1.5);
}
double nv(double T) {
    const double m = 2.525e-9 * T * T * T - 4.689e-6 * T * T + 3.376e-3 * T + 3.426e-1;
    return 4.83e15 * m * std::pow(T, 1.5);
}
double kT(double T) { return NiTCAD::base::thermal_voltage(T); }
// Couderc et al., eq (15).
double ni_eq15(double T) { return 1.541e15 * std::pow(T, 1.712) * std::exp(-gap(T) / (2.0 * kT(T))); }

// The central difference of f at T with step 1e-4 T (truncation about 1e-8 relative).
template <class F>
double central(F f, double T) {
    const double h = 1e-4 * T;
    return (f(T + h) - f(T - h)) / (2.0 * h);
}

}  // namespace

TEST_CASE("band set: Couderc's Eg, Nc and Nv as published (DC2)", "[physics][band_set]") {
    const Semiconductor si = silicon();
    for (const double T : {50.0, 77.0, 100.0, 200.0, 300.0, 400.0}) {
        CAPTURE(T);
        CHECK(relative(band_gap_eV(si, T), gap(T)) < 1e-12);
        CHECK(relative(conduction_band_dos(si, T), nc(T)) < 1e-12);
        CHECK(relative(valence_band_dos(si, T), nv(T)) < 1e-12);
        CHECK(relative(intrinsic_density(si, T),
                       std::sqrt(nc(T) * nv(T)) * std::exp(-gap(T) / (2.0 * kT(T)))) < 1e-12);
    }
    // The values DECISIONS.md states at 300 K.
    CHECK(std::abs(band_gap_eV(si, 300.0) - 1.12416) < 5e-6);
    CHECK(relative(conduction_band_dos(si, 300.0), 2.857e19) < 5e-4);
    CHECK(relative(valence_band_dos(si, 300.0), 2.514e19) < 5e-4);
}

TEST_CASE("band set: n_i against eq (15) and the measurement (DC3)", "[physics][band_set]") {
    const Semiconductor si = silicon();
    const double ni = intrinsic_density(si, 300.0);
    CAPTURE(ni);
    CHECK(ni >= 9.6e9);  // Misiakos and Tsamakis: (9.7 +/- 0.1)e9 cm^-3
    CHECK(ni <= 9.8e9);
    CHECK(relative(ni, 9.673e9) < 5e-4);
    CHECK(relative(ni_eq15(300.0), 9.68e9) < 1e-3);  // the paper's stated 300 K value
    for (double T = 50.0; T <= 400.0; T += 10.0) {
        CAPTURE(T, intrinsic_density(si, T), ni_eq15(T));
        CHECK(relative(intrinsic_density(si, T), ni_eq15(T)) < 0.015);  // the paper's precision
    }
}

TEST_CASE("band set: the temperature slopes (DC6)", "[physics][band_set]") {
    const Semiconductor si = silicon();
    const Semiconductor legacy = legacy_silicon();
    for (const Semiconductor* m : {&si, &legacy}) {
        for (const double T : {100.0, 300.0, 400.0}) {
            CAPTURE(T, m == &si);
            CHECK(relative(band_gap_slope_eV_per_K(*m, T),
                           central([&](double t) { return band_gap_eV(*m, t); }, T)) < 1e-6);
            CHECK(relative(dos_log_slope(*m, Carrier::electron, T),
                           central([&](double t) { return std::log(conduction_band_dos(*m, t)); },
                                   T)) < 1e-6);
            CHECK(relative(dos_log_slope(*m, Carrier::hole, T),
                           central([&](double t) { return std::log(valence_band_dos(*m, t)); },
                                   T)) < 1e-6);
            CHECK(relative(intrinsic_density_log_slope(*m, T),
                           central([&](double t) { return std::log(intrinsic_density(*m, t)); },
                                   T)) < 1e-6);
            for (const Carrier c : {Carrier::electron, Carrier::hole}) {
                CHECK(relative(emission_velocity_slope(*m, c, T),
                               central([&](double t) { return emission_velocity_cm_s(*m, c, t); },
                                       T)) < 1e-6);
            }
        }
    }
    // With a Richardson constant, v = A T^2 / (q N).
    SemiconductorParameters p = silicon_parameters;
    p.richardson = {.electron = 252.0, .hole = 32.0};
    const Semiconductor r = *Semiconductor::create(p);
    for (const Carrier c : {Carrier::electron, Carrier::hole}) {
        CHECK(relative(emission_velocity_slope(r, c, 300.0),
                       central([&](double t) { return emission_velocity_cm_s(r, c, t); }, 300.0)) <
              1e-6);
    }
    // Without a mass fit the mass term is exactly 0.
    CHECK(dos_mass_log(legacy, Carrier::electron, 350.0).value == 0.0);
    CHECK(dos_mass_log(legacy, Carrier::hole, 350.0).d_dT == 0.0);
}

TEST_CASE("band set: a mixed or incomplete band form is refused", "[physics][band_set]") {
    const auto refused = [](const SemiconductorParameters& p, const char* name) {
        const auto m = Semiconductor::create(p);
        REQUIRE_FALSE(m.has_value());
        CHECK(m.error().code == NiTCAD::base::ErrorCode::invalid_input);
        CHECK(m.error().message.find(name) != std::string::npos);
    };
    SemiconductorParameters p = silicon_parameters;
    p.varshni_beta_K = 636.0;
    refused(p, "varshni");
    p = silicon_parameters;
    p.passler_gap.Theta_K = 0.0;
    refused(p, "Theta_K");
    p = silicon_parameters;
    p.passler_gap.Delta = 1.0;
    refused(p, "Delta");
    p = legacy_silicon_parameters;
    p.passler_gap.Delta = 0.5;  // without alpha
    refused(p, "alpha_eV_per_K");
    p = silicon_parameters;
    p.Nc300 = 2.86e19;
    refused(p, "Nc300");
    p = silicon_parameters;
    p.dos_mass.valence[0] = p.dos_mass.valence[1] = p.dos_mass.valence[2] =
        p.dos_mass.valence[3] = 0.0;
    refused(p, "dos_mass.valence");
    p = legacy_silicon_parameters;
    p.Nv300 = 0.0;
    refused(p, "Nv300");
    // A mass fit that turns negative is refused at that temperature: silicon's conduction cubic
    // near 2030 K (its mobility is refused from 857 K, so it is tested with a temperature-flat one).
    p = silicon_parameters;
    p.electron_mobility.T_exponent = p.hole_mobility.T_exponent = 0.0;
    const Semiconductor flat = *Semiconductor::create(p);
    CHECK(check_temperature(flat, 1000.0).has_value());
    const auto hot = check_temperature(flat, 2500.0);
    REQUIRE_FALSE(hot.has_value());
    CHECK(hot.error().message.find("mass fit") != std::string::npos);
}
