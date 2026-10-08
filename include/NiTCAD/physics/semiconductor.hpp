// Bulk semiconductor parameter set and its temperature-dependent band quantities.
//
// Formulas and silicon values are the legacy ones (materials.py, Semiconductor dataclass;
// core/include/tcad/physics/materials.hpp):
//     Eg(T)  = Eg0 - alpha T^2 / (T + beta)        Varshni [eV]
//     Nc(T)  = Nc300 (T / 300)^1.5                 [cm^-3], Nv likewise
//     n_i(T) = sqrt(Nc Nv) exp(-Eg / (2 k T))      non-degenerate (Boltzmann) limit [cm^-3]
// Carried: permittivity, band structure, Caughey-Thomas mobility and Scharfetter SRH lifetimes
// (Unit 5), Auger coefficients and Slotboom band-gap narrowing (Unit 11), electron affinity (Unit
// 12, for gate work functions; Unit 15, band offsets), Canali velocity saturation (Unit 13).
// Effective masses come with the unit that uses them. Material sets: silicon (Unit 5); germanium,
// GaAs, In0.53Ga0.47As, 4H-SiC and Al_x Ga_1-x As (Unit 15, legacy materials.py M11-S1).
//
// A Semiconductor is validated once, at construction, so the model functions can rely on it.
#pragma once

#include <expected>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::physics {

enum class Carrier { electron, hole };

// Caughey-Thomas doping-dependent mobility for one carrier (mobility.hpp).
struct CaugheyThomasParameters {
    double mu_min;      // [cm^2/(V s)], temperature independent
    double mu_max;      // [cm^2/(V s)] at 300 K
    double N_ref;       // [cm^-3]
    double alpha;       // exponent on N / N_ref
    double T_exponent;  // mu_max(T) = mu_max (T / 300)^T_exponent

    bool operator==(const CaugheyThomasParameters&) const = default;
};

// Scharfetter doping-dependent SRH lifetimes: tau = tau0 / (1 + N / N_ref) (recombination.hpp).
struct ScharfetterLifetimeParameters {
    double tau_n0;  // [s]
    double tau_p0;  // [s]
    double N_ref;   // [cm^-3], shared by both carriers

    bool operator==(const ScharfetterLifetimeParameters&) const = default;
};

// Auger recombination coefficients (recombination.hpp).
struct AugerParameters {
    double Cn;  // [cm^6/s]
    double Cp;  // [cm^6/s]

    bool operator==(const AugerParameters&) const = default;
};

// Canali velocity saturation for one carrier (field_mobility.hpp).
struct CanaliParameters {
    double v_sat_cm_s;  // saturation velocity [cm/s], temperature independent (legacy)
    double beta;        // exponent

    bool operator==(const CanaliParameters&) const = default;
};

// Dopant ionization levels (ionization.hpp; Unit 15): the donor level's depth below E_c and the
// acceptor level's height above E_v, and their degeneracy factors (legacy 2 and 4). A depth of 0
// leaves that dopant completely ionized: no level is modelled for it.
struct IonizationParameters {
    double donor_eV;
    double acceptor_eV;
    double donor_degeneracy;
    double acceptor_degeneracy;

    bool operator==(const IonizationParameters&) const = default;
};

// Effective Richardson constants for thermionic emission (thermionic_emission.hpp; Unit 15)
// [A/(cm^2 K^2)]; 0 derives the emission velocity from the band's own density of states.
struct RichardsonParameters {
    double electron;
    double hole;

    bool operator==(const RichardsonParameters&) const = default;
};

// Slotboom band-gap narrowing (bandgap_narrowing.hpp).
struct SlotboomParameters {
    double E0_eV;  // [eV]
    double N0;     // [cm^-3]; no narrowing at or below N0

    bool operator==(const SlotboomParameters&) const = default;
};

// Impact ionization of one carrier, van Overstraeten-de Man (impact_ionization.hpp; Unit 19):
// alpha(E) = A exp(-B / E), the low branch (A_low, B_low) below switch_V_per_cm and the high branch
// from it on. A of 0 on both branches: the carrier does not ionize.
struct ImpactIonizationCoefficients {
    double A_low_per_cm;
    double B_low_V_per_cm;
    double A_high_per_cm;
    double B_high_V_per_cm;
    double switch_V_per_cm;
    bool operator==(const ImpactIonizationCoefficients&) const = default;
};

// Both carriers and the optical-phonon energy of the temperature factor (0: no temperature
// dependence). All zero (the default): the material does not ionize.
struct ImpactIonizationParameters {
    ImpactIonizationCoefficients electron;
    ImpactIonizationCoefficients hole;
    double phonon_energy_eV;
    bool operator==(const ImpactIonizationParameters&) const = default;
};

struct SemiconductorParameters {
    double eps_r;                   // relative permittivity
    double Eg0_eV;                  // band gap at 0 K
    double varshni_alpha_eV_per_K;
    double varshni_beta_K;
    double Nc300;                   // conduction-band effective density of states at 300 K [cm^-3]
    double Nv300;                   // valence-band effective density of states at 300 K [cm^-3]
    double electron_affinity_eV;    // chi: vacuum level to conduction-band edge [eV]
    CaugheyThomasParameters electron_mobility;
    CaugheyThomasParameters hole_mobility;
    ScharfetterLifetimeParameters lifetime;
    AugerParameters auger;
    SlotboomParameters bandgap_narrowing;
    CanaliParameters electron_saturation;
    CanaliParameters hole_saturation;
    double radiative_cm3_s;          // radiative (band-to-band) coefficient B [cm^3/s] (Unit 15)
    IonizationParameters ionization;  // Unit 15
    RichardsonParameters richardson;  // Unit 15
    ImpactIonizationParameters impact_ionization{};  // Unit 19; all zero: none

    bool operator==(const SemiconductorParameters&) const = default;
};

// Silicon, the legacy defaults (materials.py SILICON): band and density-of-states values and the
// electron affinity measured / from band structure; mobility a Caughey-Thomas fit; lifetimes a
// Scharfetter fit; Auger coefficients measured (Dziewior and Schmid); band-gap narrowing a Slotboom
// fit; velocity saturation the Canali fit; the legacy hydrogenic 45 meV levels for B, P and As
// (M13); impact ionization van Overstraeten and de Man, Solid-State Electron. 13, 583 (1970), the
// legacy values (one electron branch; the hole branches switching at 4e5 V/cm) with an
// optical-phonon energy of 0.063 eV for the temperature factor (Unit 19). The radiative coefficient
// is left 0 as in the legacy (silicon's is about 1e-14 cm^3/s, negligible against SRH), so
// silicon results do not change.
inline constexpr SemiconductorParameters silicon_parameters{
    .eps_r = 11.7,
    .Eg0_eV = 1.17,
    .varshni_alpha_eV_per_K = 4.73e-4,
    .varshni_beta_K = 636.0,
    .Nc300 = 2.86e19,
    .Nv300 = 3.10e19,
    .electron_affinity_eV = 4.05,
    .electron_mobility = {.mu_min = 92.0, .mu_max = 1360.0, .N_ref = 1.3e17, .alpha = 0.91,
                          .T_exponent = -2.33},
    .hole_mobility = {.mu_min = 47.7, .mu_max = 495.0, .N_ref = 6.3e16, .alpha = 0.76,
                      .T_exponent = -2.23},
    .lifetime = {.tau_n0 = 1.0e-5, .tau_p0 = 3.0e-6, .N_ref = 5.0e16},
    .auger = {.Cn = 2.8e-31, .Cp = 9.9e-32},
    .bandgap_narrowing = {.E0_eV = 6.92e-3, .N0 = 1.3e17},
    .electron_saturation = {.v_sat_cm_s = 1.07e7, .beta = 2.0},
    .hole_saturation = {.v_sat_cm_s = 8.37e6, .beta = 1.0},
    .radiative_cm3_s = 0.0,
    .ionization = {.donor_eV = 0.045, .acceptor_eV = 0.045, .donor_degeneracy = 2.0,
                   .acceptor_degeneracy = 4.0},
    .richardson = {.electron = 0.0, .hole = 0.0},
    .impact_ionization = {.electron = {.A_low_per_cm = 7.03e5, .B_low_V_per_cm = 1.231e6,
                                       .A_high_per_cm = 7.03e5, .B_high_V_per_cm = 1.231e6,
                                       .switch_V_per_cm = 5.0e5},
                          .hole = {.A_low_per_cm = 1.582e6, .B_low_V_per_cm = 2.036e6,
                                   .A_high_per_cm = 6.71e5, .B_high_V_per_cm = 1.693e6,
                                   .switch_V_per_cm = 4.0e5},
                          .phonon_energy_eV = 0.063},
};

// No band-gap narrowing: the Slotboom form and numbers are a silicon fit (Unit 15; the legacy left
// silicon's on every material).
inline constexpr SlotboomParameters no_bandgap_narrowing{.E0_eV = 0.0, .N0 = 1.3e17};

// The other legacy parameter sets (materials.py, M11-S1 and the 4H-SiC set; Unit 15). Band
// parameters are handbook values (Adachi); the legacy states its mobility, Auger and lifetime
// numbers as empirical fits with a 10-30% literature spread. Every field the legacy set left at its
// dataclass default keeps the silicon value, as there: the SRH lifetimes and, where not given,
// Auger and the Canali exponents. Added at Unit 15 (not in the legacy): no band-gap narrowing (the
// legacy applied silicon's Slotboom fit); radiative coefficients (Ioffe database; nextnano for
// GaAs); dopant levels where a shallow-level value is established (GaAs Si donor 5.8 meV and C
// acceptor 26.3 meV, Phys. Rev. B 43, 14734; Ge P 12.0 and B 10.4 meV; 4H-SiC N 70 and Al 220 meV,
// Ayalew, TU Wien thesis), elsewhere 0 (complete ionization); Richardson constants 0 (the
// density-of-states velocity, as the legacy; its Schottky table lists Si 252 / 32, Ge 143 / 41,
// GaAs 8 / 74 A/(cm^2 K^2) for electrons / holes, to be set by the user).
inline constexpr SemiconductorParameters germanium_parameters{
    .eps_r = 16.2,
    .Eg0_eV = 0.744,
    .varshni_alpha_eV_per_K = 4.774e-4,
    .varshni_beta_K = 235.0,
    .Nc300 = 1.04e19,
    .Nv300 = 6.0e18,
    .electron_affinity_eV = 4.13,
    .electron_mobility = {.mu_min = 0.0, .mu_max = 3900.0, .N_ref = 1.3e17, .alpha = 0.91,
                          .T_exponent = -2.33},
    .hole_mobility = {.mu_min = 0.0, .mu_max = 1900.0, .N_ref = 6.3e16, .alpha = 0.76,
                      .T_exponent = -2.23},
    .lifetime = silicon_parameters.lifetime,
    .auger = {.Cn = 2.8e-31, .Cp = 9.9e-32},
    .bandgap_narrowing = no_bandgap_narrowing,
    .electron_saturation = {.v_sat_cm_s = 7.0e6, .beta = 2.0},
    .hole_saturation = {.v_sat_cm_s = 6.3e6, .beta = 1.0},
    .radiative_cm3_s = 6.4e-14,
    .ionization = {.donor_eV = 0.012, .acceptor_eV = 0.0104, .donor_degeneracy = 2.0,
                   .acceptor_degeneracy = 4.0},
    .richardson = {.electron = 0.0, .hole = 0.0},
};

// The Gamma-valley gap.
inline constexpr SemiconductorParameters gallium_arsenide_parameters{
    .eps_r = 12.9,
    .Eg0_eV = 1.519,
    .varshni_alpha_eV_per_K = 5.405e-4,
    .varshni_beta_K = 204.0,
    .Nc300 = 4.7e17,
    .Nv300 = 7.0e18,
    .electron_affinity_eV = 4.07,
    .electron_mobility = {.mu_min = 0.0, .mu_max = 8500.0, .N_ref = 1.3e17, .alpha = 0.91,
                          .T_exponent = -2.33},
    .hole_mobility = {.mu_min = 0.0, .mu_max = 400.0, .N_ref = 6.3e16, .alpha = 0.76,
                      .T_exponent = -2.23},
    .lifetime = silicon_parameters.lifetime,
    .auger = {.Cn = 1.0e-30, .Cp = 1.0e-31},
    .bandgap_narrowing = no_bandgap_narrowing,
    .electron_saturation = {.v_sat_cm_s = 1.2e7, .beta = 2.0},
    .hole_saturation = {.v_sat_cm_s = 9.0e6, .beta = 1.0},
    .radiative_cm3_s = 7.2e-10,
    .ionization = {.donor_eV = 0.0058, .acceptor_eV = 0.0263, .donor_degeneracy = 2.0,
                   .acceptor_degeneracy = 4.0},
    .richardson = {.electron = 0.0, .hole = 0.0},
};

// In0.53Ga0.47As, lattice-matched to InP.
inline constexpr SemiconductorParameters indium_gallium_arsenide_parameters{
    .eps_r = 13.9,
    .Eg0_eV = 0.817,
    .varshni_alpha_eV_per_K = 5.78e-4,
    .varshni_beta_K = 296.0,
    .Nc300 = 2.1e17,
    .Nv300 = 7.7e18,
    .electron_affinity_eV = 4.55,
    .electron_mobility = {.mu_min = 0.0, .mu_max = 12000.0, .N_ref = 1.3e17, .alpha = 0.91,
                          .T_exponent = -2.33},
    .hole_mobility = {.mu_min = 0.0, .mu_max = 300.0, .N_ref = 6.3e16, .alpha = 0.76,
                      .T_exponent = -2.23},
    .lifetime = silicon_parameters.lifetime,
    .auger = {.Cn = 1.0e-30, .Cp = 1.0e-31},
    .bandgap_narrowing = no_bandgap_narrowing,
    .electron_saturation = {.v_sat_cm_s = 1.0e7, .beta = 2.0},
    .hole_saturation = {.v_sat_cm_s = 8.0e6, .beta = 1.0},
    .radiative_cm3_s = 0.96e-10,
    .ionization = {.donor_eV = 0.0, .acceptor_eV = 0.0, .donor_degeneracy = 2.0,
                   .acceptor_degeneracy = 4.0},
    .richardson = {.electron = 0.0, .hole = 0.0},
};

// 4H-SiC. The legacy marks eps_r, chi, Eg(300 K) = 3.23 eV and the mobility end points as
// established (Baliga); the Caughey-Thomas shape, the Varshni fit (tuned only to Eg(300 K)), the
// lifetimes (growth dependent) and Auger as carried-over approximations. Its radiative value is
// the bimolecular coefficient (Ioffe; the purely radiative part is about 1e-14 cm^3/s). The
// dopant levels are deep at 300 K (aluminium 220 meV): use incomplete ionization with it.
inline constexpr SemiconductorParameters silicon_carbide_4h_parameters{
    .eps_r = 9.7,
    .Eg0_eV = 3.30,
    .varshni_alpha_eV_per_K = 3.3e-4,
    .varshni_beta_K = 0.0,
    .Nc300 = 1.7e19,
    .Nv300 = 2.5e19,
    .electron_affinity_eV = 3.17,
    .electron_mobility = {.mu_min = 40.0, .mu_max = 950.0, .N_ref = 1.94e17, .alpha = 0.61,
                          .T_exponent = -2.15},
    .hole_mobility = {.mu_min = 15.0, .mu_max = 120.0, .N_ref = 1.76e19, .alpha = 0.34,
                      .T_exponent = -2.15},
    .lifetime = {.tau_n0 = 5.0e-7, .tau_p0 = 5.0e-7, .N_ref = 5.0e16},
    .auger = {.Cn = 2.8e-31, .Cp = 9.9e-32},
    .bandgap_narrowing = no_bandgap_narrowing,
    .electron_saturation = {.v_sat_cm_s = 2.0e7, .beta = 1.0},
    .hole_saturation = {.v_sat_cm_s = 2.0e7, .beta = 1.0},
    .radiative_cm3_s = 1.5e-12,
    .ionization = {.donor_eV = 0.070, .acceptor_eV = 0.220, .donor_degeneracy = 2.0,
                   .acceptor_degeneracy = 4.0},
    .richardson = {.electron = 0.0, .hole = 0.0},
};

// Al_x Ga_1-x As in the direct-gap regime, 0 <= x <= 0.45 (legacy materials.algaas): eps_r linear
// between GaAs and AlAs (12.9 - 2.6 x), Eg0 = 1.519 + 1.247 x, Nc and Nv scaled by (1 + 0.5 x) and
// (1 - 0.3 x) (the legacy calls this crude), mu_max 8500 - 5500 x and 400 - 150 x, and the
// conduction band taking `conduction_share` of the gap step: chi = 4.07 - share 1.247 x. The
// default share, 0.85 / 1.247 = 0.68, is the legacy's chi slope (its comment says about 85%);
// measured GaAs/AlGaAs offsets are nearer 0.62-0.65. Every other field is GaAs's (the legacy took
// silicon's defaults, including mu_min 92 and 47.7; OLD / NEW in ARCHITECTURE.md 5, Unit 15), with
// the radiative coefficient 1.8e-10 cm^3/s (Ioffe) and complete ionization (AlGaAs donors form DX
// centres, which a single shallow level does not describe).
// Errors: invalid_input if x is not in [0, 0.45] (the gap is indirect beyond) or the share not in
// [0, 1], with the offending value as the context value.
[[nodiscard]] std::expected<SemiconductorParameters, base::Error> algaas_parameters(
    double x, double conduction_share = 0.85 / 1.247);

class Semiconductor {
public:
    // Errors: invalid_input if a parameter is not finite, or eps_r, Eg0, Nc300, Nv300, mu_max,
    // N_ref, alpha, a lifetime, the band-gap-narrowing N0 or a dopant degeneracy is not positive,
    // or varshni_alpha, varshni_beta, the electron affinity, mu_min, an Auger coefficient, the
    // band-gap-narrowing E0, the radiative coefficient, a dopant level or a Richardson constant is
    // negative, or mu_min > mu_max, or a saturation velocity is not positive, or a Canali beta is
    // below 1, or an impact-ionization coefficient A, B, switch field or phonon energy is negative
    // (B and the switch field positive where A is). The message names the parameter; the context
    // value is it.
    [[nodiscard]] static std::expected<Semiconductor, base::Error> create(
        const SemiconductorParameters& parameters);

    [[nodiscard]] const SemiconductorParameters& parameters() const noexcept {
        return parameters_;
    }

private:
    explicit Semiconductor(const SemiconductorParameters& p) noexcept : parameters_(p) {}

    SemiconductorParameters parameters_;
};

// Semiconductor::create(silicon_parameters), which cannot fail.
[[nodiscard]] Semiconductor silicon();

// Whether the models of this material are usable at lattice temperature T. The device or solve
// layer calls this when it validates its input; the model functions below only require T to be
// finite and positive. Errors (invalid_input, the context value is the offending quantity):
// - T not finite and positive;
// - Eg(T) <= 0 (the Varshni fit has run past zero; n_i would exceed sqrt(Nc Nv));
// - for either carrier, mu_max(T) = mu_max (T / 300)^T_exponent not finite or below mu_min, so
//   the Caughey-Thomas mobility would rise with doping. Silicon: holes from about 857 K,
//   electrons from about 953 K. Create() applies the same rule at 300 K.
[[nodiscard]] std::expected<void, base::Error> check_temperature(const Semiconductor& m,
                                                                 double temperature_K);

// Temperature-dependent band quantities. Precondition (NITCAD_EXPECTS): temperature_K is finite
// and positive. The Varshni fit is stated by the legacy to hold for silicon over 0-500 K; outside
// that range the formulas are evaluated as written (see check_temperature).
[[nodiscard]] double band_gap_eV(const Semiconductor& m, double temperature_K);
[[nodiscard]] double conduction_band_dos(const Semiconductor& m, double temperature_K);
[[nodiscard]] double valence_band_dos(const Semiconductor& m, double temperature_K);
// 1.0674e10 cm^-3 for silicon at 300 K.
[[nodiscard]] double intrinsic_density(const Semiconductor& m, double temperature_K);
// The intrinsic level's depth below the vacuum level, chi + (E_c - E_i) with
// E_c - E_i = Eg / 2 + (k T / 2) ln(Nc / Nv) [eV] (Unit 15). Where two materials meet, the
// difference of this depth is the step of the potential at which each holds n = n_i: the band
// offsets in the equations' reference (assemble, the band shift).
[[nodiscard]] double intrinsic_level_depth_eV(const Semiconductor& m, double temperature_K);

}  // namespace NiTCAD::physics
