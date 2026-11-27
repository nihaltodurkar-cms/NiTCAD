// Unit 68 slice 4a (MAD11; DECISIONS.md amendment 8): `sic_4h_2004`, 4H-SiC with the <0001>
// impact-ionization law of Hatakeyama et al. (2004); `sic_4h` unchanged.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <variant>

#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/materials.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;

TEST_CASE("sic_4h_2004: sic_4h plus the Hatakeyama <0001> law, nothing else", "[sic][unit68]") {
    const MaterialEntry* e = find_material("sic_4h_2004");
    REQUIRE(e != nullptr);
    CHECK(e->recorded);
    const auto made = e->make({});
    REQUIRE(made.has_value());
    SemiconductorParameters p = std::get<SemiconductorParameters>(*made);
    CHECK(p.ionization.donor_eV == 0.070);     // incomplete ionization kept (owner ruling 9)
    CHECK(p.ionization.acceptor_eV == 0.220);
    p.impact_ionization = silicon_carbide_4h_parameters.impact_ionization;
    CHECK(p == silicon_carbide_4h_parameters);
    // `sic_4h` itself still does not ionize, and its set is the one it was.
    CHECK(std::get<SemiconductorParameters>(*find_material("sic_4h")->make({})) == silicon_carbide_4h_parameters);
    CHECK(silicon_carbide_4h_parameters.impact_ionization == ImpactIonizationParameters{});
}

TEST_CASE("sic_4h_2004: the coefficients at 1, 2 and 3 MV/cm (owner check 7)", "[sic][unit68]") {
    const auto made = find_material("sic_4h_2004")->make({});
    const Semiconductor m = *Semiconductor::create(std::get<SemiconductorParameters>(*made));
    const ImpactIonizationCoefficients& n = impact_ionization(m, Carrier::electron);
    const ImpactIonizationCoefficients& h = impact_ionization(m, Carrier::hole);
    struct Row {
        double F, alpha_n, alpha_p;  // owner's table, from the unrounded Table I values
    };
    for (const Row& r : {Row{1e6, 8.19964e-7, 4.73579e-3}, Row{2e6, 1.20131e1, 1.27079e3}, Row{3e6, 2.93950e3, 8.19660e4}}) {
        CAPTURE(r.F);
        const double an = impact_ionization_coefficient(n, 1.0, r.F).alpha;
        const double ap = impact_ionization_coefficient(h, 1.0, r.F).alpha;
        CHECK(std::abs(an / r.alpha_n - 1.0) < 1e-5);  // the table's six digits
        CHECK(std::abs(ap / r.alpha_p - 1.0) < 1e-5);
        CHECK(an == 1.76e8 * std::exp(-3.30e7 / r.F));  // the source's expression, bit for bit
        CHECK(ap == 3.41e8 * std::exp(-2.50e7 / r.F));
    }
    // No temperature factor: the 300 K law at every temperature.
    CHECK(impact_ionization_temperature_factor(0.0, 450.0) == 1.0);
}
