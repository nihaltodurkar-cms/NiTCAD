// Material sets, the intrinsic-level depth and the thermionic emission velocity (ARCHITECTURE.md
// section 11, Unit 15; legacy materials.py M11-S1 and the 4H-SiC set, device.py emission_velocity,
// tests/test_m33_interface.py S2).
//
// References computed here in double-double arithmetic (references.hpp) from the legacy formulas
// and parameters, independently of this code.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/ionization.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "references.hpp"

using namespace NiTCAD;
using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

}  // namespace

TEST_CASE("materials: every legacy set validates and gives its band quantities at 300 K") {
    // References from the legacy formulas in double-double arithmetic (references.hpp): the
    // Varshni gap, n_i = sqrt(Nc Nv) e^(-Eg / 2kT) and the depth chi + Eg/2 + (kT/2) ln(Nc/Nv).
    struct Set {
        const char* name;
        SemiconductorParameters parameters;
    };
    const Set sets[] = {
        {"Si", silicon_parameters},
        {"Ge", germanium_parameters},
        {"GaAs", gallium_arsenide_parameters},
        {"In0.53Ga0.47As", indium_gallium_arsenide_parameters},
        {"4H-SiC", silicon_carbide_4h_parameters},
        {"Al0.3Ga0.7As", *algaas_parameters(0.3)},
    };
    for (const Set& r : sets) {
        CAPTURE(r.name);
        const auto m = Semiconductor::create(r.parameters);
        REQUIRE(m.has_value());
        REQUIRE(check_temperature(*m, 300.0).has_value());
        REQUIRE(close(band_gap_eV(*m, 300.0), reference::band_gap(r.parameters, 300.0).value(),
                      1e-14));
        REQUIRE(close(intrinsic_density(*m, 300.0),
                      reference::intrinsic_density(r.parameters, 300.0).value(), 1e-12));
        REQUIRE(close(intrinsic_level_depth_eV(*m, 300.0),
                      reference::intrinsic_level_depth(r.parameters, 300.0).value(), 1e-14));
    }
    // The legacy comment says Eg(300 K) of 4H-SiC is about 3.23 eV; its own Varshni numbers give
    // 3.201.
    REQUIRE(close(band_gap_eV(*Semiconductor::create(silicon_carbide_4h_parameters), 300.0),
                  3.201, 1e-14));
}

TEST_CASE("materials: AlGaAs follows the legacy interpolation inside the direct-gap regime") {
    // x = 0 is GaAs except the radiative coefficient (AlGaAs's) and the dopant levels (complete
    // ionization): the legacy took silicon's defaults for the fields it did not interpolate.
    const SemiconductorParameters g = *algaas_parameters(0.0);
    SemiconductorParameters gaas = gallium_arsenide_parameters;
    gaas.radiative_cm3_s = 1.8e-10;
    gaas.ionization.donor_eV = 0.0;
    gaas.ionization.acceptor_eV = 0.0;
    REQUIRE(g == gaas);
    // The conduction-band share sets chi: 0.62 of the 1.247 x gap step.
    REQUIRE(close(algaas_parameters(0.3, 0.62)->electron_affinity_eV, 4.07 - 0.62 * 1.247 * 0.3,
                  1e-15));
    REQUIRE_FALSE(algaas_parameters(0.3, 1.2).has_value());
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
    // N = 2 (2 pi m kT / h^2)^1.5, then v = sqrt(kT / (2 pi m)), in double-double arithmetic
    // (references.hpp) by that route, not the code's simplified kT / h (2 / N)^(1/3). Silicon Nc
    // gives the legacy 2.575e6 cm/s.
    using reference::DD;
    const auto velocity = [](double N_cm3, double T) {
        const DD kT = DD(base::k_B_J_per_K) * DD(T);
        const DD h = DD(2.0) * reference::pi() * DD(base::hbar_J_s);
        const DD N = DD(N_cm3) * DD(1e6);  // m^-3
        const DD m = reference::pow(N / DD(2.0), DD(2.0) / DD(3.0)) * h * h /
                     (DD(2.0) * reference::pi() * kT);
        return (reference::sqrt(kT / (DD(2.0) * reference::pi() * m)) * DD(100.0)).value();
    };
    for (const double N : {2.86e19, 4.7e17, 7.0e18}) {
        CAPTURE(N);
        REQUIRE(close(emission_velocity_cm_s(N, 300.0), velocity(N, 300.0), 1e-14));
    }
    REQUIRE(close(velocity(2.86e19, 300.0), 2.575e6, 2e-4));
    // Scales as T (N fixed) and N^(-1/3).
    REQUIRE(close(emission_velocity_cm_s(1e19, 600.0), 2.0 * emission_velocity_cm_s(1e19, 300.0),
                  1e-15));
    REQUIRE(close(emission_velocity_cm_s(8e19, 300.0), 0.5 * emission_velocity_cm_s(1e19, 300.0),
                  1e-15));
}

// Unit 15 fixes: incomplete ionization, radiative recombination, Richardson constants.

TEST_CASE("ionization: the ionized fraction, its derivative and its limits") {
    // N / (1 + g e^(eta + E/kT)) with an exact eta derivative; complete for a zero level depth;
    // no overflow at either end.
    const double N = 3e17, depth = 1.7, g = 4.0;
    for (const double eta : {-40.0, -5.0, -1.7 - std::log(4.0), 0.0, 3.0, 40.0}) {
        CAPTURE(eta);
        const IonizedDensity d = ionized_density(N, eta, depth, g);
        REQUIRE(close(d.value, N / (1.0 + g * std::exp(eta + depth)), 1e-14));
        const double h = 1e-6;
        const double fd = (ionized_density(N, eta + h, depth, g).value -
                           ionized_density(N, eta - h, depth, g).value) /
                          (2.0 * h);
        REQUIRE(std::abs(d.d_eta - fd) <= 1e-8 * N);
        REQUIRE(d.d_eta <= 0.0);
    }
    REQUIRE(ionized_density(N, -1000.0, depth, g).value == N);
    REQUIRE(ionized_density(N, 1000.0, depth, g).value == 0.0);
    REQUIRE(ionized_density(N, 2.0, 0.0, g).value == N);  // no level modelled
    REQUIRE(ionized_density(N, 2.0, 0.0, g).d_eta == 0.0);
}

TEST_CASE("ionization: freeze-out of boron in silicon against double-double roots (legacy G7)") {
    // N_A = 1e16, Fermi-Dirac statistics, no band-gap narrowing: the ionized fraction from the
    // neutral root, against bisection in double-double arithmetic with the reference F_{1/2}
    // (references.hpp), and the legacy's literature bands (Sze and Ng freeze-out curves;
    // Altermatt et al. 2002): 77 K 15-45%, 150 K 70-98%, 250 K >= 85%, 300 K >= 95%.
    struct Case {
        double T, lo, hi;
    };
    const Semiconductor si = silicon();
    for (const Case& c : {Case{77.0, 0.15, 0.45}, Case{150.0, 0.70, 0.98},
                          Case{250.0, 0.85, 1.01}, Case{300.0, 0.95, 1.01}}) {
        CAPTURE(c.T);
        const double ni = intrinsic_density(si, c.T);
        const double VT = base::thermal_voltage(c.T);
        const DopantLevels levels{0.045 / VT, 0.045 / VT, 2.0, 4.0};
        const double gn = std::log(conduction_band_dos(si, c.T) / ni);
        const double gp = std::log(valence_band_dos(si, c.T) / ni);
        const NeutralEquilibrium e =
            ionized_neutral_equilibrium(0.0, 1e16, ni, gn, gp, levels, true);
        const reference::Neutral r = reference::neutral_equilibrium(
            {.acceptors = 1e16, .donor_kT = levels.donor_kT, .acceptor_kT = levels.acceptor_kT},
            ni, gn, gp, true);
        const double fraction = (e.p - e.n) / 1e16;
        CAPTURE(fraction, ((r.p - r.n) / reference::DD(1e16)).value());
        REQUIRE(close(fraction, ((r.p - r.n) / reference::DD(1e16)).value(), 1e-10));
        REQUIRE(fraction >= c.lo);
        REQUIRE(fraction <= c.hi);
    }
}

TEST_CASE("ionization: 4H-SiC dopants at 300 K against double-double roots") {
    // 1e17 cm^-3: nitrogen (70 meV) is 87% ionized, aluminium (220 meV) only 10.6%, the reason
    // p-type SiC needs incomplete ionization at room temperature.
    const Semiconductor sic = *Semiconductor::create(silicon_carbide_4h_parameters);
    const double T = 300.0, VT = base::thermal_voltage(T);
    const double ni = intrinsic_density(sic, T);
    const double gn = std::log(conduction_band_dos(sic, T) / ni);
    const double gp = std::log(valence_band_dos(sic, T) / ni);
    const DopantLevels levels{0.070 / VT, 0.220 / VT, 2.0, 4.0};
    const auto n_type = ionized_neutral_equilibrium(1e17, 0.0, ni, gn, gp, levels, false);
    const auto p_type = ionized_neutral_equilibrium(0.0, 1e17, ni, gn, gp, levels, false);
    const auto p_fd = ionized_neutral_equilibrium(0.0, 1e17, ni, gn, gp, levels, true);
    const reference::Dopants nd{.donors = 1e17, .donor_kT = levels.donor_kT,
                                .acceptor_kT = levels.acceptor_kT};
    const reference::Dopants na{.acceptors = 1e17, .donor_kT = levels.donor_kT,
                                .acceptor_kT = levels.acceptor_kT};
    const reference::Neutral rn = reference::neutral_equilibrium(nd, ni, gn, gp, false);
    const reference::Neutral rp = reference::neutral_equilibrium(na, ni, gn, gp, false);
    const reference::Neutral rf = reference::neutral_equilibrium(na, ni, gn, gp, true);
    const reference::DD N = 1e17;
    REQUIRE(close((n_type.n - n_type.p) / 1e17, ((rn.n - rn.p) / N).value(), 1e-10));
    REQUIRE(close((p_type.p - p_type.n) / 1e17, ((rp.p - rp.n) / N).value(), 1e-10));
    REQUIRE(close((p_fd.p - p_fd.n) / 1e17, ((rf.p - rf.n) / N).value(), 1e-10));
    REQUIRE(close((n_type.n - n_type.p) / 1e17, 0.867, 1e-3));
    REQUIRE(close((p_type.p - p_type.n) / 1e17, 0.106, 1e-2));
    // Complete-ionization limit: a zero level depth gives the ordinary neutral root.
    const DopantLevels none{0.0, 0.0, 2.0, 4.0};
    const auto complete = ionized_neutral_equilibrium(1e17, 0.0, ni, gn, gp, none, false);
    REQUIRE(std::abs(complete.eta - boltzmann_neutral_equilibrium(1e17, ni).eta) <=
            1e-12 * complete.eta);
}

TEST_CASE("recombination: radiative is B (n p - E) with the partials of E") {
    const EquilibriumProduct E{4.0e20, 3.0e3, 2.0e3};
    const double B = 7.2e-10, n = 1e17, p = 5e10;
    const RecombinationRate r = radiative_recombination(n, p, E, B);
    REQUIRE(r.rate == B * (n * p - E.value));
    REQUIRE(r.d_dn == B * (p - E.d_dn));
    REQUIRE(r.d_dp == B * (n - E.d_dp));
    REQUIRE(radiative_recombination(2e10, 2e10, {4e20, 0.0, 0.0}, B).rate == 0.0);
}

TEST_CASE("thermionic: a set Richardson constant gives A* T^2 / (q N)") {
    // The legacy Schottky table's silicon electron value, 252 A/(cm^2 K^2), gives 4.950e6 cm/s at
    // 300 K (legacy emission_velocity docstring), 1.92 times the density-of-states velocity.
    SemiconductorParameters p = silicon_parameters;
    p.richardson = {.electron = 252.0, .hole = 0.0};
    const Semiconductor m = *Semiconductor::create(p);
    const double Nc = conduction_band_dos(m, 300.0);
    const double v = emission_velocity_cm_s(m, Carrier::electron, 300.0);
    REQUIRE(close(v, 252.0 * 300.0 * 300.0 / (base::q_C * Nc), 1e-15));
    REQUIRE(close(v, 4.950e6, 1e-3));
    REQUIRE(close(v / emission_velocity_cm_s(Nc, 300.0), 1.92, 0.01));
    // An unset constant falls back to the density of states.
    REQUIRE(emission_velocity_cm_s(m, Carrier::hole, 300.0) ==
            emission_velocity_cm_s(valence_band_dos(m, 300.0), 300.0));
    SemiconductorParameters bad = silicon_parameters;
    bad.richardson.hole = -1.0;
    REQUIRE_FALSE(Semiconductor::create(bad).has_value());
    bad = silicon_parameters;
    bad.ionization.donor_degeneracy = 0.0;
    REQUIRE_FALSE(Semiconductor::create(bad).has_value());
    bad = silicon_parameters;
    bad.radiative_cm3_s = -1e-10;
    REQUIRE_FALSE(Semiconductor::create(bad).has_value());
}
