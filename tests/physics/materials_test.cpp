// Material sets, the intrinsic-level depth and the thermionic emission velocity (ARCHITECTURE.md
// section 11, Unit 15; legacy materials.py M11-S1 and the 4H-SiC set, device.py emission_velocity,
// tests/test_m33_interface.py S2).
//
// References computed at 40 digits (mpmath) from the legacy formulas and parameters, independently
// of this code.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"

using namespace NiTCAD;
using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

struct Reference {
    const char* name;
    SemiconductorParameters parameters;
    double Eg, ni, depth;  // at 300 K: eV, cm^-3, eV
};

}  // namespace

TEST_CASE("materials: every legacy set validates and gives its band quantities at 300 K") {
    const Reference references[] = {
        {"Si", silicon_parameters, 1.1245192307692307692, 10673775147.815624025,
         4.6112180320226613287},
        {"Ge", germanium_parameters, 0.6636897196261682243, 21030229406096.730196,
         4.4689547587053674452},
        {"GaAs", gallium_arsenide_parameters, 1.4224821428571428571, 2043025.6913016548437,
         4.7463288152059162278},
        {"In0.53Ga0.47As", indium_gallium_arsenide_parameters, 0.72971812080536912752,
         944044715427.80413982, 4.8683013140223912968},
        // The legacy comment says Eg(300 K) is about 3.23 eV; its own Varshni numbers give 3.201.
        {"4H-SiC", silicon_carbide_4h_parameters, 3.201, 2.6728700597565474886e-8,
         4.7655149268142061874},
        {"Al0.3Ga0.7As", *algaas_parameters(0.3), 1.7965821428571428571, 1506.058167826050866,
         4.6814044378909079971},
    };
    for (const Reference& r : references) {
        CAPTURE(r.name);
        const auto m = Semiconductor::create(r.parameters);
        REQUIRE(m.has_value());
        REQUIRE(check_temperature(*m, 300.0).has_value());
        REQUIRE(close(band_gap_eV(*m, 300.0), r.Eg, 1e-14));
        REQUIRE(close(intrinsic_density(*m, 300.0), r.ni, 1e-12));
        REQUIRE(close(intrinsic_level_depth_eV(*m, 300.0), r.depth, 1e-14));
    }
}

TEST_CASE("materials: AlGaAs follows the legacy interpolation inside the direct-gap regime") {
    // x = 0 has GaAs's band parameters (the other fields are the silicon defaults, as in the
    // legacy, so it is not GaAs's full set).
    const SemiconductorParameters g = *algaas_parameters(0.0);
    const SemiconductorParameters& gaas = gallium_arsenide_parameters;
    REQUIRE(g.eps_r == gaas.eps_r);
    REQUIRE(g.Eg0_eV == gaas.Eg0_eV);
    REQUIRE(g.electron_affinity_eV == gaas.electron_affinity_eV);
    REQUIRE(g.Nc300 == gaas.Nc300);
    REQUIRE(g.Nv300 == gaas.Nv300);
    REQUIRE(g.electron_mobility.mu_min == silicon_parameters.electron_mobility.mu_min);
    // The conduction band takes 0.85 / 1.247 of the gap step.
    const SemiconductorParameters a = *algaas_parameters(0.45);
    REQUIRE(close((g.electron_affinity_eV - a.electron_affinity_eV) / (a.Eg0_eV - g.Eg0_eV),
                  0.85 / 1.247, 1e-14));
    REQUIRE(Semiconductor::create(a).has_value());
    for (const double x : {-0.01, 0.46, 1.0, std::nan("")}) {
        CAPTURE(x);
        const auto bad = algaas_parameters(x);
        REQUIRE_FALSE(bad.has_value());
        REQUIRE(bad.error().code == base::ErrorCode::invalid_input);
    }
}

TEST_CASE("materials: the intrinsic-level depth is chi plus E_c - E_i") {
    // The band offsets in the equations' reference: across a junction of two materials the step of
    // the depth is the step of the potential at which each side holds n = n_i. For two materials
    // that differ only in chi the step is the affinity step.
    SemiconductorParameters shifted = silicon_parameters;
    shifted.electron_affinity_eV += 0.3;
    const Semiconductor si = silicon();
    const Semiconductor s2 = *Semiconductor::create(shifted);
    for (const double T : {200.0, 300.0, 400.0}) {
        CAPTURE(T);
        REQUIRE(close(intrinsic_level_depth_eV(s2, T) - intrinsic_level_depth_eV(si, T), 0.3,
                      1e-13));
        const double ec_ei = 0.5 * band_gap_eV(si, T) +
                             0.5 * base::thermal_voltage(T) *
                                 std::log(conduction_band_dos(si, T) / valence_band_dos(si, T));
        REQUIRE(close(intrinsic_level_depth_eV(si, T), 4.05 + ec_ei, 1e-15));
    }
}

TEST_CASE("thermionic: the emission velocity is sqrt(kT / 2 pi m) with m from N") {
    // Legacy test_s2_emission_velocity_matches_the_closed_form: m_DOS from
    // N = 2 (2 pi m kT / h^2)^1.5, at 40 digits. Silicon Nc gives the legacy 2.575e6 cm/s.
    struct Case {
        double N, v;
    };
    for (const Case& c : {Case{2.86e19, 2575351.3998235678652}, Case{4.7e17, 10129606.016810534997},
                          Case{7.0e18, 4117110.2172364250936}}) {
        CAPTURE(c.N);
        REQUIRE(close(emission_velocity_cm_s(c.N, 300.0), c.v, 1e-14));
    }
    // Scales as T (N fixed) and N^(-1/3).
    REQUIRE(close(emission_velocity_cm_s(1e19, 600.0), 2.0 * emission_velocity_cm_s(1e19, 300.0),
                  1e-15));
    REQUIRE(close(emission_velocity_cm_s(8e19, 300.0), 0.5 * emission_velocity_cm_s(1e19, 300.0),
                  1e-15));
}
