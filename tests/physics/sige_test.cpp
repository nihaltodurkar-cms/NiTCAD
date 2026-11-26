// Unit 68 slice 1 (MAD7, MAC9; DECISIONS.md amendment 5): SiGe sets.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>

#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/sige.hpp"
#include "NiTCAD/physics/thermal.hpp"

using namespace NiTCAD::physics;

namespace {
Semiconductor make(const SemiconductorParameters& p) { return *Semiconductor::create(p); }
double varshni_change(double a, double b, double T) { return a * 300.0 * 300.0 / (300.0 + b) - a * T * T / (T + b); }
}  // namespace

TEST_CASE("SiGe: x = 0 is NiTCAD silicon's band and transport set (S2, F3a)", "[sige][unit68]") {
    const SemiconductorParameters p = *sige_parameters(0.0, false);
    const SemiconductorParameters& s = silicon_parameters;
    REQUIRE(p.eps_r == s.eps_r);
    REQUIRE(p.Eg0_eV == s.Eg0_eV);
    REQUIRE(p.electron_affinity_eV == s.electron_affinity_eV);
    REQUIRE(p.electron_mobility == s.electron_mobility);
    REQUIRE(p.hole_mobility == s.hole_mobility);
    REQUIRE(p.dos_mass == s.dos_mass);
    REQUIRE(p.lifetime == s.lifetime);
    REQUIRE(p.auger == s.auger);
    REQUIRE(p.bandgap_narrowing == s.bandgap_narrowing);
    REQUIRE(p.electron_saturation == s.electron_saturation);
    REQUIRE(p.ionization == s.ionization);
    REQUIRE(p.alloy_gap.fraction == 0.0);
    for (const double T : {250.0, 300.0, 400.0}) REQUIRE(band_gap_eV(make(p), T) == band_gap_eV(make(s), T));
    REQUIRE(*sige_parameters(0.0, true) == p);
    REQUIRE(strained_silicon_parameters(0.0)->Eg0_eV == s.Eg0_eV);
}

TEST_CASE("SiGe: relaxed band edges, gap fit and mobility lines (M7, S2)", "[sige][unit68]") {
    const SiGeBands b = *sige_bands(0.2, 0.2);
    REQUIRE(b.in_plane_strain == 0.0);
    REQUIRE(std::abs(b.gap_change - (-0.41 * 0.2 + 0.008 * 0.04)) < 1e-15);
    REQUIRE(std::abs(b.valence_edge_offset - (0.458 * 0.2 + 0.246 * 0.2 / 3.0)) < 1e-15);
    const SemiconductorParameters p = *sige_parameters(0.2, false);
    REQUIRE(std::abs(band_gap_eV(make(p), 300.0) - band_gap_eV(make(silicon_parameters), 300.0) - b.gap_change) < 1e-12);
    REQUIRE(std::abs(p.electron_mobility.mu_max - 1360.0 * (1396.0 - 863.0) / 1396.0) < 1e-9);
    REQUIRE(std::abs(p.hole_mobility.mu_max - 495.0 * (450.0 - 173.0) / 450.0) < 1e-9);
    REQUIRE(p.eps_r == 11.7 + 0.9);
    REQUIRE(std::abs(p.ionization.acceptor_eV - (0.045 - 0.00692)) < 1e-15);
    REQUIRE(std::abs(p.bandgap_narrowing.E0_eV - 0.9 * silicon_parameters.bandgap_narrowing.E0_eV) < 1e-18);
}

TEST_CASE("SiGe: the strained offset MAC9 checks, and the strained band edges", "[sige][unit68]") {
    // MAC9 (S10): sige_strained 0.2 on Si, the Rieger-Vogl weighted-average valence offset
    // (0.47 - 0.06 x) x = 0.0916 eV as printed.
    const SiGeBands b = *sige_bands(0.2, 0.0);
    REQUIRE(std::abs(b.average_valence_offset - 0.0916) < 0.00005);
    std::printf("MAC9 sige_strained 0.2: average valence offset %.4f eV, valence edge offset %.4f eV, gap change "
                "%.4f eV, in-plane strain %.5f\n",
                b.average_valence_offset, b.valence_edge_offset, b.gap_change, b.in_plane_strain);
    const SiGeBands r = *sige_bands(0.2, 0.2);
    REQUIRE(b.in_plane_strain < 0.0);                        // compressive on silicon
    REQUIRE(b.valence_edge_offset > r.valence_edge_offset);  // the heavy hole rises
    REQUIRE(b.gap_change < r.gap_change);                    // strain narrows the gap further
    // Strained silicon on relaxed SiGe: tensile, the Delta2 valleys drop below relaxed silicon's.
    const SiGeBands t = *sige_bands(0.0, 0.3);
    REQUIRE(t.in_plane_strain > 0.0);
    REQUIRE(t.conduction_edge_offset < 0.0);
    REQUIRE(t.gap_change < 0.0);
    const SemiconductorParameters s = *strained_silicon_parameters(0.3);
    REQUIRE(s.electron_mobility == silicon_parameters.electron_mobility);  // S11
}

TEST_CASE("SiGe: the gap's temperature correction (F3a) and its slope", "[sige][unit68]") {
    const double x = 0.25;
    const Semiconductor m = make(*sige_parameters(x, false));
    const Semiconductor si = make(silicon_parameters);
    const double gap300 = band_gap_eV(m, 300.0) - band_gap_eV(si, 300.0);
    for (const double T : {200.0, 350.0, 450.0}) {
        const double expected = gap300 + x * (varshni_change(4.774e-4, 235.0, T) - varshni_change(4.73e-4, 636.0, T));
        REQUIRE(std::abs(band_gap_eV(m, T) - band_gap_eV(si, T) - expected) < 1e-12);
        const double h = 1e-3;
        const double fd = (band_gap_eV(m, T + h) - band_gap_eV(m, T - h)) / (2.0 * h);
        REQUIRE(std::abs(band_gap_slope_eV_per_K(m, T) - fd) < 1e-9);
    }
}

TEST_CASE("SiGe: refusals (S3)", "[sige][unit68]") {
    REQUIRE_FALSE(sige_parameters(0.31, false).has_value());
    REQUIRE_FALSE(sige_parameters(-0.01, true).has_value());
    REQUIRE_FALSE(strained_silicon_parameters(0.5).has_value());
    REQUIRE(sige_parameters(0.30, true).has_value());
    REQUIRE(Semiconductor::create(*sige_parameters(0.30, true)).has_value());
}
