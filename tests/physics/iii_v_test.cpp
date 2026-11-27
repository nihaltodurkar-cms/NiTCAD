// Unit 68 slice 3 (MAD10; DECISIONS.md amendment 7): the III-V binaries and alloys: MAC3 (alloy
// ends, the owner's tabulated compositions, refusals), MAC4 (Sotoodeh mobility, binaries gated, the
// lattice-matched alloy rows reported), gaps and alignment, the saturation temperature law, alloy
// hole scattering and thermal conductivity, the direct-gap tunnelling masses.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <string>

#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/iii_v.hpp"
#include "NiTCAD/physics/materials.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;

namespace {
Semiconductor make(const SemiconductorParameters& p) {
    auto m = Semiconductor::create(p);
    if (!m) FAIL(m.error().message);
    return *m;
}
SemiconductorParameters ingaas(double x) { return *ingaas_parameters(x); }
SemiconductorParameters inalas(double x) { return *inalas_parameters(x); }
double close(double a, double b) { return std::abs(a - b) / std::abs(b); }
// The density-of-states mass back from Nc or Nv at 300 K (Nc = 2.509e19 (m / m0)^(3/2)).
double mass(const SemiconductorParameters& reference, double N300) {
    return std::pow(N300 / reference.Nc300, 2.0 / 3.0) * 0.067;  // against gaas_2001's 0.067
}
// Sotoodeh's form as printed (owner amendment 7), independent of the engine's expression.
double sotoodeh(double mu_max, double mu_min, double N_ref, double alpha, double theta1, double theta2, double N,
                double T) {
    return mu_min + (mu_max * std::pow(300.0 / T, theta1) - mu_min) / (1.0 + std::pow(N / (N_ref * std::pow(T / 300.0, theta2)), alpha));
}
}  // namespace

TEST_CASE("III-V: every set builds, gaps at 300 K and the VBO alignment", "[iii_v][unit68]") {
    struct Row {
        const char* name;
        SemiconductorParameters p;
        double gap300;  // owner amendment 7 (rounded)
        double chi;     // the owner's derived affinity (from rounded gaps)
    };
    const Row rows[] = {{"InP", inp_parameters(), 1.353, 4.281},  {"InAs", inas_parameters(), 0.354, 4.930},
                        {"GaP", gap_parameters(), 2.273, 3.69},   {"AlAs", alas_parameters(), 2.164, 3.86},
                        {"GaSb", gasb_parameters(), 0.727, 4.00}, {"InSb", insb_parameters(), 0.173, 4.52}};
    const Semiconductor gaas = make(gaas_2001_parameters());
    CHECK(gaas_2001_parameters().electron_affinity_eV == 4.07);
    for (const Row& r : rows) {
        CAPTURE(r.name);
        const Semiconductor m = make(r.p);
        const double Eg = band_gap_eV(m, 300.0);
        CHECK(std::abs(Eg - r.gap300) < 1.5e-3);
        // chi = 4.07 + Eg_GaAs + VBO_GaAs - VBO - Eg: the offsets reproduce, the gaps exactly.
        CHECK(std::abs(r.p.electron_affinity_eV - r.chi) < 0.01);
        std::printf("III-V %-5s Eg(300) %.4f eV, chi %.4f eV, Nc %.4g Nv %.4g cm^-3\n", r.name, Eg,
                    r.p.electron_affinity_eV, r.p.Nc300, r.p.Nv300);
    }
    // Indirect GaP and AlAs: the X ellipsoids, three valleys.
    CHECK(close(mass(gaas_2001_parameters(), gap_parameters().Nc300), std::cbrt(2.00 * 0.253 * 0.253) * std::pow(3.0, 2.0 / 3.0)) < 1e-12);
    CHECK(close(mass(gaas_2001_parameters(), alas_parameters().Nc300), std::cbrt(0.97 * 0.22 * 0.22) * std::pow(3.0, 2.0 / 3.0)) < 1e-12);
    CHECK_FALSE(gap_parameters().band_to_band.direct_gap);
    CHECK_FALSE(alas_parameters().band_to_band.direct_gap);
    // Holes: the [001] heavy and light holes combined; InP's (0.532, 0.121) -> 0.570 (owner).
    CHECK(std::abs(mass(gaas_2001_parameters(), inp_parameters().Nv300) - 0.570) < 1e-3);
    CHECK(std::abs(mass(gaas_2001_parameters(), inas_parameters().Nv300) - 0.338) < 1e-3);
    // SRH: InSb and the alloys give no lifetimes.
    CHECK(insb_parameters().lifetime.tau_n0 == 0.0);
    CHECK(ingaas(0.53).lifetime.tau_p0 == 0.0);
    CHECK(gaas_2001_parameters().lifetime.tau_n0 == 1.8e-7);
    CHECK(gaas_2001_parameters().lifetime.tau_p0 == 3.0e-8);
}

TEST_CASE("III-V MAC3: alloy ends, the owner's compositions, refusals", "[iii_v][unit68]") {
    CHECK(ingaas(0.0) == gaas_2001_parameters());
    CHECK(ingaas(1.0) == inas_parameters());
    CHECK(inalas(0.0) == alas_parameters());
    CHECK(inalas(1.0) == inas_parameters());
    for (const double x : {-0.01, 1.01, std::nan("")}) {
        for (const auto& made : {ingaas_parameters(x), inalas_parameters(x)}) {
            REQUIRE_FALSE(made.has_value());
            CHECK(made.error().message.find("[0, 1]") != std::string::npos);
        }
    }
    // Continuity at the ends (the gap, the affinity), within the bowing's first-order term.
    for (const double x : {1e-9, 1.0 - 1e-9}) {
        const SemiconductorParameters& end = x < 0.5 ? gaas_2001_parameters() : inas_parameters();
        const Semiconductor a = make(ingaas(x)), b = make(end);
        CHECK(std::abs(band_gap_eV(a, 300.0) - band_gap_eV(b, 300.0)) < 1e-8);
        CHECK(std::abs(ingaas(x).electron_affinity_eV - end.electron_affinity_eV) < 1e-8);
    }
    // In0.53Ga0.47As (owner D23): m_e 0.043, m_lh 0.052; their reduced mass 0.0235 (the owner's
    // table prints 0.0237, which those two masses do not give).
    const SemiconductorParameters lm = ingaas(0.53);
    CHECK(std::abs(lm.band_to_band.electron_mass - 0.043) < 5e-4);
    CHECK(std::abs(lm.band_to_band.hole_mass - 0.052) < 5e-4);
    const double mr = lm.band_to_band.electron_mass * lm.band_to_band.hole_mass /
                      (lm.band_to_band.electron_mass + lm.band_to_band.hole_mass);
    CHECK(std::abs(mr - 0.043 * 0.052 / 0.095) < 2e-4);
    // InAlAs: X lowest on the Al-rich side, Gamma beyond the crossover (owner: x_In about 0.36).
    double crossover = 0.0;
    for (double x = 0.30; x <= 0.45; x += 1e-4) {
        if (inalas(x).band_to_band.direct_gap) {
            crossover = x;
            break;
        }
    }
    std::printf("III-V InAlAs Gamma-X crossover x_In = %.4f; In0.53Ga0.47As Eg(300) %.4f eV, "
                "In0.52Al0.48As %.4f eV\n",
                crossover, band_gap_eV(make(lm), 300.0), band_gap_eV(make(inalas(0.52)), 300.0));
    CHECK(std::abs(crossover - 0.36) < 0.02);
    CHECK_FALSE(inalas(0.2).band_to_band.direct_gap);
    // The gap is continuous through the crossover (the lower valley), and so is the affinity.
    CHECK(std::abs(band_gap_eV(make(inalas(crossover)), 300.0) - band_gap_eV(make(inalas(crossover - 1e-4)), 300.0)) <
          1e-3);
}

TEST_CASE("III-V MAC4: Sotoodeh mobility, binaries within 2 %, the alloy rows reported", "[iii_v][unit68]") {
    struct Fit {
        const char* name;
        SemiconductorParameters p;
        double e[6], h[6];  // mu_max mu_min N_ref alpha theta1 theta2 (owner amendment 7)
    };
    const Fit fits[] = {
        {"GaAs", gaas_2001_parameters(), {9400, 500, 6.0e16, 0.394, 2.10, 3.0}, {491.5, 20, 1.48e17, 0.38, 2.20, 3.0}},
        {"InP", inp_parameters(), {5200, 400, 3.0e17, 0.47, 2.00, 3.25}, {170, 10, 4.87e17, 0.62, 2.00, 3.0}},
        {"InAs", inas_parameters(), {34000, 1000, 1.1e18, 0.32, 1.57, 3.0}, {530, 20, 1.1e17, 0.46, 2.30, 3.0}},
        {"GaP", gap_parameters(), {152, 10, 4.4e18, 0.80, 1.60, 0.71}, {147, 10, 1.0e18, 0.85, 1.98, 0.0}},
        {"AlAs", alas_parameters(), {400, 10, 5.46e17, 1.00, 2.10, 3.0}, {200, 10, 3.84e17, 0.488, 2.24, 3.0}},
        {"GaSb", gasb_parameters(), {5650, 1050, 2.8e17, 1.05, 2.00, 2.8}, {875, 190, 9.0e17, 0.65, 1.70, 2.7}},
        {"InSb", insb_parameters(), {78000, 5000, 7.0e16, 0.70, 1.66, 3.0}, {750, 100, 6.0e17, 0.60, 1.80, 3.0}}};
    for (const Fit& f : fits) {
        CAPTURE(f.name);
        const Semiconductor m = make(f.p);
        for (const double T : {300.0, 400.0}) {
            for (const double N : {1e15, 1e16, 1e17, 1e18, 1e19}) {
                CAPTURE(T, N);
                CHECK(close(caughey_thomas_mobility(m, Carrier::electron, N, T),
                            sotoodeh(f.e[0], f.e[1], f.e[2], f.e[3], f.e[4], f.e[5], N, T)) < 0.02);
                CHECK(close(caughey_thomas_mobility(m, Carrier::hole, N, T),
                            sotoodeh(f.h[0], f.h[1], f.h[2], f.h[3], f.h[4], f.h[5], N, T)) < 0.02);
            }
        }
    }
    // The lattice-matched rows (owner): the reciprocal rule on the binaries, reported, not gated
    // (owner ruling 2026-11-27, MAC4 amended).
    struct Row {
        const char* name;
        SemiconductorParameters p;
        Carrier c;
        double fit[6];
    };
    const Row rows[] = {{"In0.53Ga0.47As e", ingaas(0.53), Carrier::electron, {14000, 300, 1.3e17, 0.48, 1.59, 3.68}},
                        {"In0.53Ga0.47As h", ingaas(0.53), Carrier::hole, {320, 10, 4.9e17, 0.403, 1.59, 3.0}},
                        {"In0.52Al0.48As e", inalas(0.52), Carrier::electron, {4800, 800, 3.0e16, 1.10, 1.85, 3.0}}};
    for (const Row& r : rows) {
        const Semiconductor m = make(r.p);
        std::printf("III-V MAC4 reference %s:", r.name);
        for (const double N : {1e15, 1e16, 1e17, 1e18, 1e19}) {
            const double ours = caughey_thomas_mobility(m, r.c, N, 300.0);
            const double ref = sotoodeh(r.fit[0], r.fit[1], r.fit[2], r.fit[3], r.fit[4], r.fit[5], N, 300.0);
            std::printf(" %.0e: %.0f vs %.0f (%+.0f %%)", N, ours, ref, 100.0 * (ours / ref - 1.0));
            CHECK(ours > 0.0);
        }
        std::printf("\n");
    }
}

TEST_CASE("III-V: saturation velocity's temperature law (Quay)", "[iii_v][unit68]") {
    const Semiconductor m = make(inp_parameters());
    CHECK(saturation(m, Carrier::electron, 300.0).v_sat_cm_s == 6.7e6);
    CHECK(close(saturation(m, Carrier::electron, 400.0).v_sat_cm_s, 6.7e6 / (1.0 - 0.68 + 0.68 * 400.0 / 300.0)) < 1e-15);
    CHECK(close(saturation(m, Carrier::hole, 500.0).v_sat_cm_s, 5.0e6 / (1.0 - 0.30 + 0.30 * 500.0 / 300.0)) < 1e-15);
    // A = 0 (GaSb, every earlier set): temperature independent.
    CHECK(saturation(make(gasb_parameters()), Carrier::electron, 450.0).v_sat_cm_s == 6.0e6);
    CHECK(saturation(make(silicon_parameters), Carrier::electron, 450.0) == silicon_parameters.electron_saturation);
    SemiconductorParameters bad = inp_parameters();
    bad.hole_saturation.temperature_A = 1.5;
    const auto refused = Semiconductor::create(bad);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().message.find("temperature_A") != std::string::npos);
}

TEST_CASE("III-V: alloy hole scattering, thermal conductivity, tunnelling masses", "[iii_v][unit68]") {
    // C_mu 1e6 on both carriers: the hole mobility sits just below the reciprocal rule's.
    SemiconductorParameters p = ingaas(0.5);
    REQUIRE(p.alloy_mobility.C_mu_hole == 1.0e6);
    const double with = caughey_thomas_mobility(make(p), Carrier::hole, 1e17, 300.0);
    p.alloy_mobility.C_mu_hole = 0.0;
    const double without = caughey_thomas_mobility(make(p), Carrier::hole, 1e17, 300.0);
    CHECK(close(1.0 / with - 1.0 / without, 0.25 / 1.0e6) < 1e-9);
    // Thermal: the fitted bowing; In0.53Ga0.47As about 5 and In0.52Al0.48As about 1 W/(m K);
    // InAlAs between about x 0.54 and 0.74 not positive: no thermal data.
    CHECK(std::abs(ingaas(0.53).thermal.conductivity_W_cmK - 0.050) < 0.002);
    CHECK(std::abs(inalas(0.52).thermal.conductivity_W_cmK - 0.0106) < 0.001);
    CHECK(inalas(0.60).thermal.conductivity_W_cmK == 0.0);
    CHECK(inalas(0.80).thermal.conductivity_W_cmK > 0.0);
    CHECK(gaas_2001_parameters().thermal.conductivity_W_cmK == 0.55);
    // Direct-gap tunnelling masses (MAD10): Gamma electron, [001] light hole 1 / (gamma1 + 2 gamma2).
    CHECK(gaas_2001_parameters().band_to_band.electron_mass == 0.067);
    CHECK(close(gaas_2001_parameters().band_to_band.hole_mass, 1.0 / (6.98 + 2.0 * 2.06)) < 1e-15);
    CHECK(close(inp_parameters().band_to_band.hole_mass, 1.0 / (5.08 + 3.2)) < 1e-15);
    CHECK(close(inas_parameters().band_to_band.hole_mass, 1.0 / 37.0) < 1e-15);
    CHECK(gasb_parameters().band_to_band.electron_mass == 0.0);
    CHECK(inalas(0.52).band_to_band.electron_mass == 0.0);
}

TEST_CASE("III-V: the registry resolves ingaas by its argument count", "[iii_v][unit68]") {
    const MaterialEntry* legacy = find_material("ingaas", 0);
    const MaterialEntry* alloy = find_material("ingaas", 1);
    REQUIRE(legacy != nullptr);
    REQUIRE(alloy != nullptr);
    CHECK(legacy != alloy);
    CHECK_FALSE(legacy->recorded);
    CHECK(alloy->recorded);
    CHECK(std::get<SemiconductorParameters>(*legacy->make({})) == indium_gallium_arsenide_parameters);
    const double x[] = {0.53};
    CHECK(std::get<SemiconductorParameters>(*alloy->make(x)) == ingaas(0.53));
    CHECK(find_material("ingaas", 2) == alloy);  // its usage refuses the count
    CHECK(find_material("nope", 0) == nullptr);
}
