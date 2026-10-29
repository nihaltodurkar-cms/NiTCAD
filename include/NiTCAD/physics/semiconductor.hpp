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

// Band-to-band tunnelling (band_to_band.hpp; Unit 20). The Kane pair (A, B) of the local and the
// calibrated nonlocal rate G = A F^2 exp(-B / F); A = 0: none. The direct-gap WKB rate needs a
// direct gap and both tunnelling masses (conduction and valence band, in units of m0) from a cited
// source; 0: none (no material set ships them, as none is verified). All zero: no tunnelling.
struct BandToBandParameters {
    double A_per_cm3_s;
    double B_V_per_cm;
    bool direct_gap;
    double electron_mass;  // m_c / m0
    double hole_mass;      // m_v / m0
    bool operator==(const BandToBandParameters&) const = default;
};

// Lattice heat conduction and storage, and the thermopower's scattering exponents (thermal.hpp;
// Unit 23): kappa(T) = conductivity (T / 300)^-conductivity_exponent [W/(cm K)], rho c the
// volumetric heat capacity [J/(cm^3 K)], temperature independent. thermopower_exponent_n / _p is
// the exponent r of the carriers' momentum relaxation time in their energy, tau ~ E^r, which sets
// the energy a current carries (the thermopower, thermal.hpp); -1/2 is acoustic-phonon scattering
// (M. Lundstrom, Fundamentals of Carrier Transport, 2nd ed., Cambridge 2000, sec. 3.4).
// A conductivity of 0 (the default) means the material has no thermal data: the electrothermal
// model refuses it rather than invent values.
struct ThermalParameters {
    double conductivity_W_cmK = 0.0;
    double conductivity_exponent = 0.0;
    double heat_capacity_J_cm3K = 0.0;
    double thermopower_exponent_n = -0.5;
    double thermopower_exponent_p = -0.5;

    bool operator==(const ThermalParameters&) const = default;
};

// Passler's band gap (Unit 35; R. Passler, Phys. Rev. B 66, 085201 (2002), as given by Couderc et
// al., J. Appl. Phys. 115, 093705 (2014), eq (10)), with chi = 2 T / Theta:
//     Eg(T) = Eg0 - alpha Theta [ (1 - 3 Delta^2) / (exp(Theta / T) - 1)
//             + (3 Delta^2 / 2) ( (1 + pi^2 chi^2 / (3 (1 + Delta^2)) + (3 Delta^2 - 1) chi^3 / 4
//                                  + 8 chi^4 / 3 + chi^6)^(1/6) - 1 ) ]
// alpha of 0 (the default): no Passler gap, the Varshni form instead.
struct PasslerGapParameters {
    double alpha_eV_per_K = 0.0;  // slope of the high-temperature asymptote
    double Theta_K = 0.0;         // mean phonon temperature
    double Delta = 0.0;           // phonon dispersion coefficient, in [0, 1)

    bool operator==(const PasslerGapParameters&) const = default;
};

// Temperature-dependent density-of-states effective masses (Unit 35; Couderc et al. (2014), eqs
// (4), (5), (18), (19)): (m_d / m0)^(3/2) = c[0] T^3 + c[1] T^2 + c[2] T + c[3] per band, and
//     Nc(T) = 4.83e15 (m_dc / m0)^(3/2) T^(3/2) [cm^-3], Nv likewise,
// 4.83e15 = 2 (2 pi m0 k / h^2)^(3/2) in cm^-3 K^-3/2 as printed there. All zero (the default):
// no mass fit, Nc300 (T / 300)^(3/2) instead.
struct DensityOfStatesMassParameters {
    double conduction[4] = {0.0, 0.0, 0.0, 0.0};
    double valence[4] = {0.0, 0.0, 0.0, 0.0};

    bool operator==(const DensityOfStatesMassParameters&) const = default;
};

// Energy relaxation times of the energy-balance model (Unit 52; energy_balance.hpp; DECISIONS.md
// ED5) [s], constant. 0 (the default): not given; the model then takes silicon's cited values in
// a region with silicon's mobility and refuses any other.
struct EnergyRelaxationParameters {
    double tau_n_s = 0.0;
    double tau_p_s = 0.0;

    bool operator==(const EnergyRelaxationParameters&) const = default;
};

struct SemiconductorParameters {
    double eps_r;                   // relative permittivity
    double Eg0_eV;                  // band gap at 0 K
    double varshni_alpha_eV_per_K;  // 0 with a Passler gap
    double varshni_beta_K;          // 0 with a Passler gap
    double Nc300;  // conduction-band effective density of states at 300 K [cm^-3]; 0 with a mass fit
    double Nv300;  // valence-band effective density of states at 300 K [cm^-3]; 0 with a mass fit
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
    BandToBandParameters band_to_band{};              // Unit 20; all zero: none
    ThermalParameters thermal{};                      // Unit 23; conductivity 0: none
    PasslerGapParameters passler_gap{};               // Unit 35; alpha 0: Varshni
    DensityOfStatesMassParameters dos_mass{};         // Unit 35; all 0: (T / 300)^(3/2)
    EnergyRelaxationParameters energy_relaxation{};   // Unit 52; 0: not given

    bool operator==(const SemiconductorParameters&) const = default;
};

// Whether a parameter set uses Passler's gap, and a density-of-states mass fit (Unit 35).
[[nodiscard]] constexpr bool has_passler_gap(const SemiconductorParameters& p) noexcept {
    return p.passler_gap.alpha_eV_per_K != 0.0;
}
[[nodiscard]] constexpr bool has_dos_mass_fit(const SemiconductorParameters& p) noexcept {
    return !(p.dos_mass == DensityOfStatesMassParameters{});
}

// Silicon as the legacy had it (materials.py SILICON; the default silicon until Unit 35, kept
// selectable, deck material `silicon_legacy`; n_i(300 K) = 1.0674e10 cm^-3): band and
// density-of-states values and the
// electron affinity measured / from band structure; mobility a Caughey-Thomas fit; lifetimes a
// Scharfetter fit; Auger coefficients measured (Dziewior and Schmid); band-gap narrowing a Slotboom
// fit; velocity saturation the Canali fit; the legacy hydrogenic 45 meV levels for B, P and As
// (M13); impact ionization van Overstraeten and de Man, Solid-State Electron. 13, 583 (1970), the
// legacy values (one electron branch; the hole branches switching at 4e5 V/cm) with an
// optical-phonon energy of 0.063 eV for the temperature factor (Unit 19); band-to-band tunnelling
// the Kane pair A = 3.5e21 cm^-3 s^-1, B = 1.03e8 V/cm of Hurkx, Klaassen and Knuvers, IEEE Trans.
// Electron Devices 39, 331 (1992), Table I, the legacy values (Unit 20; indirect gap, no WKB
// masses). The radiative coefficient is left 0 as in the legacy (silicon's is about 1e-14 cm^3/s,
// negligible against SRH), so silicon results do not change. Thermal data (Unit 23): the legacy
// kappa(T) = 1.48 (T / 300)^-1.33 W/(cm K) (materials.py kappa_th); rho c = 2.329 g/cm^3 times
// 0.7 J/(g K) (Ioffe NSM archive, silicon thermal properties); acoustic-phonon thermopower
// exponents.
inline constexpr SemiconductorParameters legacy_silicon_parameters{
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
    .band_to_band = {.A_per_cm3_s = 3.5e21, .B_V_per_cm = 1.03e8, .direct_gap = false,
                     .electron_mass = 0.0, .hole_mass = 0.0},
    .thermal = {.conductivity_W_cmK = 1.48, .conductivity_exponent = 1.33,
                .heat_capacity_J_cm3K = 1.6303, .thermopower_exponent_n = -0.5,
                .thermopower_exponent_p = -0.5},
};

// Silicon, the default (Unit 35, DECISIONS.md S1-S2): the legacy set with the published band set of
// R. Couderc, M. Amara, M. Lemiti, J. Appl. Phys. 115, 093705 (2014): Passler's gap (Table I:
// Eg(0) = 1.17 eV, alpha = 3.23e-4 eV/K, Theta = 446 K, Delta = 0.51) and the density-of-states
// mass fits of eqs (18) and (19), the valence one fitted to the paper's n_i(T) (eq (15)).
// n_i(300 K) = 9.673e9 cm^-3 (Misiakos and Tsamakis measured (9.7 +/- 0.1)e9; Unit 33 V4); Eg(300 K)
// = 1.1242 eV, Nc = 2.857e19, Nv = 2.514e19 cm^-3. Every other field is the legacy's.
inline constexpr SemiconductorParameters silicon_parameters = [] {
    SemiconductorParameters p = legacy_silicon_parameters;
    p.Eg0_eV = 1.17;
    p.varshni_alpha_eV_per_K = 0.0;
    p.varshni_beta_K = 0.0;
    p.Nc300 = 0.0;
    p.Nv300 = 0.0;
    p.passler_gap = {.alpha_eV_per_K = 3.23e-4, .Theta_K = 446.0, .Delta = 0.51};
    p.dos_mass = {.conduction = {-4.609e-10, 6.753e-7, -1.312e-5, 1.094},
                  .valence = {2.525e-9, -4.689e-6, 3.376e-3, 3.426e-1}};
    return p;
}();

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

// The Gamma-valley gap. A direct gap (band_to_band.direct_gap), with no tunnelling masses: the
// direct-gap WKB rate needs cited ones from the user (Unit 20).
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
    .band_to_band = {.A_per_cm3_s = 0.0, .B_V_per_cm = 0.0, .direct_gap = true,
                     .electron_mass = 0.0, .hole_mass = 0.0},
};

// In0.53Ga0.47As, lattice-matched to InP; a direct gap, as GaAs (Unit 20).
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
    .band_to_band = {.A_per_cm3_s = 0.0, .B_V_per_cm = 0.0, .direct_gap = true,
                     .electron_mass = 0.0, .hole_mass = 0.0},
};

// 4H-SiC. The legacy marks eps_r, chi, Eg(300 K) = 3.23 eV and the mobility end points as
// established (Baliga); the Caughey-Thomas shape, the Varshni fit (tuned only to Eg(300 K)), the
// lifetimes (growth dependent) and Auger as carried-over approximations. Its radiative value is
// the bimolecular coefficient (Ioffe; the purely radiative part is about 1e-14 cm^3/s). The
// dopant levels are deep at 300 K (aluminium 220 meV): use incomplete ionization with it.
// Thermal data (Unit 23): the legacy conductivity 3.7 W/(cm K) at 300 K (materials.py) with the
// legacy exponent, which it applies to every material (silicon's 1.33; not checked against a
// source for SiC); rho c = 3.211 g/cm^3 times 0.69 J/(g K) (Ioffe NSM archive, 4H-SiC thermal
// properties). The other sets carry no thermal data (none of their values is in the legacy).
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
    .thermal = {.conductivity_W_cmK = 3.7, .conductivity_exponent = 1.33,
                .heat_capacity_J_cm3K = 2.2156, .thermopower_exponent_n = -0.5,
                .thermopower_exponent_p = -0.5},
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
    // (B and the switch field positive where A is), or a band-to-band A, B or tunnelling mass is
    // negative (B positive where A is; the masses both 0 or both positive, on a direct gap, with a
    // reduced mass below m0 / 2), or a band form is mixed or incomplete (Unit 35): a Passler gap
    // (alpha above 0) needs Theta above 0, Delta in [0, 1) and the Varshni alpha and beta 0, and
    // without one Theta and Delta are 0; a density-of-states mass fit needs Nc300 and Nv300 0 and
    // a nonzero polynomial for each band. The message names the parameter; the context value is it.
    [[nodiscard]] static std::expected<Semiconductor, base::Error> create(
        const SemiconductorParameters& parameters);

    [[nodiscard]] const SemiconductorParameters& parameters() const noexcept {
        return parameters_;
    }

private:
    explicit Semiconductor(const SemiconductorParameters& p) noexcept : parameters_(p) {}

    SemiconductorParameters parameters_;
};

// Semiconductor::create(silicon_parameters) and (legacy_silicon_parameters), which cannot fail.
[[nodiscard]] Semiconductor silicon();
[[nodiscard]] Semiconductor legacy_silicon();

// Whether the models of this material are usable at lattice temperature T. The device or solve
// layer calls this when it validates its input; the model functions below only require T to be
// finite and positive. Errors (invalid_input, the context value is the offending quantity):
// - T not finite and positive;
// - Eg(T) <= 0 (the Varshni fit has run past zero; n_i would exceed sqrt(Nc Nv));
// - with a density-of-states mass fit, either band's (m_d / m0)^(3/2) not positive (Unit 35;
//   silicon's conduction-band cubic turns negative near 2030 K);
// - for either carrier, mu_max(T) = mu_max (T / 300)^T_exponent not finite or below mu_min, so
//   the Caughey-Thomas mobility would rise with doping. Silicon: holes from about 857 K,
//   electrons from about 953 K. Create() applies the same rule at 300 K.
[[nodiscard]] std::expected<void, base::Error> check_temperature(const Semiconductor& m,
                                                                 double temperature_K);

// Temperature-dependent band quantities. Precondition (NITCAD_EXPECTS): temperature_K is finite
// and positive. The legacy Varshni fit is stated to hold for silicon over 0-500 K, Couderc's mass
// fits over 0-400 K (n_i within 1.5% of their eq (15) over 50-400 K); outside those ranges the
// formulas are evaluated as written (see check_temperature).
[[nodiscard]] double band_gap_eV(const Semiconductor& m, double temperature_K);
[[nodiscard]] double conduction_band_dos(const Semiconductor& m, double temperature_K);
[[nodiscard]] double valence_band_dos(const Semiconductor& m, double temperature_K);
// 9.673e9 cm^-3 for silicon at 300 K (1.0674e10 for legacy silicon).
[[nodiscard]] double intrinsic_density(const Semiconductor& m, double temperature_K);
// The intrinsic level's depth below the vacuum level, chi + (E_c - E_i) with
// E_c - E_i = Eg / 2 + (k T / 2) ln(Nc / Nv) [eV] (Unit 15). Where two materials meet, the
// difference of this depth is the step of the potential at which each holds n = n_i: the band
// offsets in the equations' reference (assemble, the band shift).
[[nodiscard]] double intrinsic_level_depth_eV(const Semiconductor& m, double temperature_K);

}  // namespace NiTCAD::physics
