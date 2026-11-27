// Unit 68 slice 4b (MAD11; DECISIONS.md amendment 9): `ga2o3_beta`, the Barnes velocity law and the
// "not given" saturation law (owner ruling 1).
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <variant>

#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/gallium_oxide.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/materials.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;

TEST_CASE("ga2o3_beta: the owner's [010] values", "[ga2o3][unit68]") {
    const MaterialEntry* e = find_material("ga2o3_beta");
    REQUIRE(e != nullptr);
    CHECK(e->recorded);
    CHECK(e->n_type_only);
    const SemiconductorParameters p = std::get<SemiconductorParameters>(*e->make({}));
    CHECK(p == ga2o3_beta_parameters());
    const auto m = Semiconductor::create(p);
    REQUIRE(m.has_value());
    CHECK(p.eps_r == 10.87);
    CHECK(band_gap_eV(*m, 300.0) == 4.43);  // the 300 K gap at every temperature
    CHECK(band_gap_eV(*m, 450.0) == 4.43);
    CHECK(p.electron_affinity_eV == 4.00);
    CHECK(p.ionization.donor_eV == 0.030);
    CHECK(p.ionization.acceptor_eV == 1.10);
    CHECK(p.hole_saturation == CanaliParameters{.v_sat_cm_s = 0.0, .beta = 0.0});  // not given
    CHECK(p.impact_ionization.hole == ImpactIonizationCoefficients{});
    const double alpha = impact_ionization_coefficient(impact_ionization(*m, Carrier::electron), 1.0, 3e6).alpha;
    CHECK(alpha == 2.16e6 * std::exp(-1.77e7 / 3e6));
    CHECK(std::abs(p.thermal.heat_capacity_J_cm3K - 2.9155) < 1e-12);
}

TEST_CASE("Barnes law: Ghosh's form, its derivative and its velocity peak", "[ga2o3][unit68]") {
    const CanaliParameters c = ga2o3_beta_parameters().electron_saturation;
    const double mu0 = 140.0;  // Ghosh's own mu0, so the curve is his Table I y curve
    CHECK(canali_mobility(mu0, 0.0, c).mobility == mu0);
    CHECK(canali_mobility(mu0, 0.0, c).d_dE == 0.0);  // gamma > 2
    double peak_v = 0.0, peak_F = 0.0;
    for (double F = 1e4; F <= 5e5; F += 1e4) {
        CAPTURE(F);
        const FieldMobility f = canali_mobility(mu0, F, c);
        const double r = std::pow(F / 1.54e5, 2.47);
        CHECK(std::abs(f.mobility / ((mu0 + 1.5e7 / F * r) / (1.0 + r)) - 1.0) < 1e-14);
        const double h = F * 1e-6;
        const double fd = (canali_mobility(mu0, F + h, c).mobility - canali_mobility(mu0, F - h, c).mobility) / (2 * h);
        CHECK(std::abs(f.d_dE - fd) <= 1e-6 * std::abs(fd) + 1e-12);
        if (f.mobility * F > peak_v) peak_v = f.mobility * F, peak_F = F;
    }
    // Ghosh Fig. 6, y: a peak near 2e7 cm/s between 200 and 300 kV/cm, then negative differential mobility.
    CHECK(peak_F > 1.5e5);
    CHECK(peak_F < 3.0e5);
    CHECK(peak_v > 1.8e7);
    CHECK(peak_v < 2.2e7);
    CHECK(canali_mobility(mu0, 5e5, c).mobility * 5e5 < peak_v);
}

TEST_CASE("Saturation law: all zero is not given; the Barnes law's checks", "[ga2o3][unit68]") {
    SemiconductorParameters p = silicon_parameters;
    p.hole_saturation = {.v_sat_cm_s = 0.0, .beta = 0.0};
    CHECK(Semiconductor::create(p).has_value());
    const auto refused = [&p](CanaliParameters c, const std::string& field) {
        p.electron_saturation = c;
        const auto made = Semiconductor::create(p);
        REQUIRE_FALSE(made.has_value());
        CHECK(made.error().message.find(field) != std::string::npos);
    };
    refused({.v_sat_cm_s = 0.0, .beta = 2.0}, "v_sat_cm_s");  // half given
    refused({.v_sat_cm_s = 1.5e7, .beta = 0.0, .F_c_V_per_cm = 1.54e5, .gamma = 1.5}, "gamma");
    refused({.v_sat_cm_s = 1.5e7, .beta = 0.0, .F_c_V_per_cm = 0.0, .gamma = 2.47}, "F_c_V_per_cm");
    refused({.v_sat_cm_s = 1.5e7, .beta = 1.0, .F_c_V_per_cm = 1.54e5, .gamma = 2.47}, "beta");
}
