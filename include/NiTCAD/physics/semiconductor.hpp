// Bulk semiconductor parameter set and its temperature-dependent band quantities.
//
// Formulas and silicon values are the legacy ones (materials.py, Semiconductor dataclass;
// core/include/tcad/physics/materials.hpp):
//     Eg(T)  = Eg0 - alpha T^2 / (T + beta)        Varshni [eV]
//     Nc(T)  = Nc300 (T / 300)^1.5                 [cm^-3], Nv likewise
//     n_i(T) = sqrt(Nc Nv) exp(-Eg / (2 k T))      non-degenerate (Boltzmann) limit [cm^-3]
// Carried: permittivity, band structure, Caughey-Thomas mobility and Scharfetter SRH lifetimes
// (Unit 5), Auger coefficients and Slotboom band-gap narrowing (Unit 11), electron affinity (Unit
// 12, for gate work functions), Canali velocity saturation (Unit 13). Effective masses come with
// the unit that uses them.
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

// Slotboom band-gap narrowing (bandgap_narrowing.hpp).
struct SlotboomParameters {
    double E0_eV;  // [eV]
    double N0;     // [cm^-3]; no narrowing at or below N0

    bool operator==(const SlotboomParameters&) const = default;
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

    bool operator==(const SemiconductorParameters&) const = default;
};

// Silicon, the legacy defaults (materials.py SILICON): band and density-of-states values and the
// electron affinity measured / from band structure; mobility a Caughey-Thomas fit; lifetimes a
// Scharfetter fit; Auger coefficients measured (Dziewior and Schmid); band-gap narrowing a Slotboom
// fit; velocity saturation the Canali fit.
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
};

class Semiconductor {
public:
    // Errors: invalid_input if a parameter is not finite, or eps_r, Eg0, Nc300, Nv300, mu_max,
    // N_ref, alpha, a lifetime or the band-gap-narrowing N0 is not positive, or varshni_alpha,
    // varshni_beta, the electron affinity, mu_min, an Auger coefficient or the band-gap-narrowing
    // E0 is negative, or mu_min > mu_max, or a saturation velocity is not positive, or a Canali
    // beta is below 1. The message names the parameter; the context value is it.
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

}  // namespace NiTCAD::physics
