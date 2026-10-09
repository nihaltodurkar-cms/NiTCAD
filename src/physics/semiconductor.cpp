#include "NiTCAD/physics/semiconductor.hpp"

#include <cmath>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "temperature.hpp"

namespace NiTCAD::physics {

namespace {

enum class Bound { positive, non_negative, finite };

// An invalid_input error naming the parameter as group + field (group is "" or ends in '.').
// The message is built only here, on the failure path.
base::Error parameter_error(std::string_view group, std::string_view field,
                            std::string_view requirement, double value) {
    std::string message{"semiconductor parameter "};
    message += group;
    message += field;
    message += requirement;
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = std::nullopt, .value = value}};
}

std::optional<base::Error> check(std::string_view group, std::string_view field, double value,
                                 Bound bound) {
    switch (bound) {
        case Bound::positive:
            if (std::isfinite(value) && value > 0.0) return std::nullopt;
            return parameter_error(group, field, " must be finite and positive", value);
        case Bound::non_negative:
            if (std::isfinite(value) && value >= 0.0) return std::nullopt;
            return parameter_error(group, field, " must be finite and non-negative", value);
        case Bound::finite:
            if (std::isfinite(value)) return std::nullopt;
            return parameter_error(group, field, " must be finite", value);
    }
    return std::nullopt;
}

// A >= 0 on both branches; a branch with A > 0 needs B > 0 (alpha would not vanish at low field);
// the switch field positive.
std::optional<base::Error> check_impact(std::string_view group,
                                        const ImpactIonizationCoefficients& c) {
    if (auto e = check(group, "A_low_per_cm", c.A_low_per_cm, Bound::non_negative)) return e;
    if (auto e = check(group, "A_high_per_cm", c.A_high_per_cm, Bound::non_negative)) return e;
    if (auto e = check(group, "switch_V_per_cm", c.switch_V_per_cm,
                       c.A_low_per_cm > 0.0 || c.A_high_per_cm > 0.0 ? Bound::positive
                                                                     : Bound::non_negative)) {
        return e;
    }
    if (auto e = check(group, "B_low_V_per_cm", c.B_low_V_per_cm,
                       c.A_low_per_cm > 0.0 ? Bound::positive : Bound::non_negative)) {
        return e;
    }
    if (auto e = check(group, "B_high_V_per_cm", c.B_high_V_per_cm,
                       c.A_high_per_cm > 0.0 ? Bound::positive : Bound::non_negative)) {
        return e;
    }
    return std::nullopt;
}

// A >= 0, B > 0 where A is (G would not vanish at low field); the tunnelling masses both 0 or
// both positive, only on a direct gap, with a reduced mass below m0 / 2 (the WKB kappa vanishes at
// both band edges only then, band_to_band.hpp).
std::optional<base::Error> check_band_to_band(const BandToBandParameters& b) {
    constexpr std::string_view group = "band_to_band.";
    if (auto e = check(group, "A_per_cm3_s", b.A_per_cm3_s, Bound::non_negative)) return e;
    if (auto e = check(group, "B_V_per_cm", b.B_V_per_cm,
                       b.A_per_cm3_s > 0.0 ? Bound::positive : Bound::non_negative)) {
        return e;
    }
    if (auto e = check(group, "electron_mass", b.electron_mass, Bound::non_negative)) return e;
    if (auto e = check(group, "hole_mass", b.hole_mass, Bound::non_negative)) return e;
    if ((b.electron_mass > 0.0) != (b.hole_mass > 0.0)) {
        return parameter_error(group, "hole_mass", " must be given with electron_mass",
                               b.hole_mass);
    }
    if (b.electron_mass > 0.0) {
        if (!b.direct_gap) {
            return parameter_error(group, "electron_mass",
                                   " needs a direct gap (tunnelling masses are for the direct-gap "
                                   "WKB rate)",
                                   b.electron_mass);
        }
        const double mr = b.electron_mass * b.hole_mass / (b.electron_mass + b.hole_mass);
        if (!(mr < 0.5)) {
            return parameter_error(group, "electron_mass",
                                   ": the reduced mass must be below m0 / 2", mr);
        }
    }
    return std::nullopt;
}

// beta >= 1 keeps dmu/dE finite at E = 0 (the legacy only required beta > 0; its silicon values
// are 2 and 1).
std::optional<base::Error> check_saturation(std::string_view group, const CanaliParameters& c) {
    if (auto e = check(group, "v_sat_cm_s", c.v_sat_cm_s, Bound::positive)) return e;
    if (auto e = check(group, "beta", c.beta, Bound::finite)) return e;
    if (c.beta < 1.0) return parameter_error(group, "beta", " must be at least 1", c.beta);
    return std::nullopt;
}

std::optional<base::Error> check_mobility(std::string_view group,
                                          const CaugheyThomasParameters& ct) {
    if (auto e = check(group, "mu_min", ct.mu_min, Bound::non_negative)) return e;
    if (auto e = check(group, "mu_max", ct.mu_max, Bound::positive)) return e;
    if (auto e = check(group, "N_ref", ct.N_ref, Bound::positive)) return e;
    if (auto e = check(group, "alpha", ct.alpha, Bound::positive)) return e;
    if (auto e = check(group, "T_exponent", ct.T_exponent, Bound::finite)) return e;
    if (ct.mu_min > ct.mu_max) {
        // Mobility would rise with doping. check_temperature applies the same rule at T.
        return parameter_error(group, "mu_min", " exceeds mu_max", ct.mu_min);
    }
    return std::nullopt;
}

double band_gap_unchecked(const SemiconductorParameters& p, double T) {
    return p.Eg0_eV - p.varshni_alpha_eV_per_K * (T * T) / (T + p.varshni_beta_K);
}

// (T / 300)^1.5, the temperature factor of both effective densities of states.
double dos_factor(double T) { return std::pow(T / 300.0, 1.5); }

}  // namespace

std::expected<Semiconductor, base::Error> Semiconductor::create(
    const SemiconductorParameters& p) {
    using enum Bound;
    if (auto e = check("", "eps_r", p.eps_r, positive)) return std::unexpected(std::move(*e));
    if (auto e = check("", "Eg0_eV", p.Eg0_eV, positive)) return std::unexpected(std::move(*e));
    if (auto e = check("", "varshni_alpha_eV_per_K", p.varshni_alpha_eV_per_K, non_negative)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("", "varshni_beta_K", p.varshni_beta_K, non_negative)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("", "Nc300", p.Nc300, positive)) return std::unexpected(std::move(*e));
    if (auto e = check("", "Nv300", p.Nv300, positive)) return std::unexpected(std::move(*e));
    if (auto e = check("", "electron_affinity_eV", p.electron_affinity_eV, non_negative)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check_mobility("electron_mobility.", p.electron_mobility)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check_mobility("hole_mobility.", p.hole_mobility)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("lifetime.", "tau_n0", p.lifetime.tau_n0, positive)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("lifetime.", "tau_p0", p.lifetime.tau_p0, positive)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("lifetime.", "N_ref", p.lifetime.N_ref, positive)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("auger.", "Cn", p.auger.Cn, non_negative)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("auger.", "Cp", p.auger.Cp, non_negative)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("bandgap_narrowing.", "E0_eV", p.bandgap_narrowing.E0_eV, non_negative)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("bandgap_narrowing.", "N0", p.bandgap_narrowing.N0, positive)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check_saturation("electron_saturation.", p.electron_saturation)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check_saturation("hole_saturation.", p.hole_saturation)) {
        return std::unexpected(std::move(*e));
    }
    const struct {
        std::string_view group, field;
        double value;
        Bound bound;
    } unit15[] = {
        {"", "radiative_cm3_s", p.radiative_cm3_s, non_negative},
        {"ionization.", "donor_eV", p.ionization.donor_eV, non_negative},
        {"ionization.", "acceptor_eV", p.ionization.acceptor_eV, non_negative},
        {"ionization.", "donor_degeneracy", p.ionization.donor_degeneracy, positive},
        {"ionization.", "acceptor_degeneracy", p.ionization.acceptor_degeneracy, positive},
        {"richardson.", "electron", p.richardson.electron, non_negative},
        {"richardson.", "hole", p.richardson.hole, non_negative},
    };
    for (const auto& c : unit15) {
        if (auto e = check(c.group, c.field, c.value, c.bound)) {
            return std::unexpected(std::move(*e));
        }
    }
    if (auto e = check_impact("impact_ionization.electron.", p.impact_ionization.electron)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check_impact("impact_ionization.hole.", p.impact_ionization.hole)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check("impact_ionization.", "phonon_energy_eV",
                       p.impact_ionization.phonon_energy_eV, non_negative)) {
        return std::unexpected(std::move(*e));
    }
    if (auto e = check_band_to_band(p.band_to_band)) return std::unexpected(std::move(*e));
    if (auto ok = check_thermal_parameters(p.thermal, "semiconductor"); !ok) {
        return std::unexpected(std::move(ok.error()));
    }
    return Semiconductor{p};
}

Semiconductor silicon() {
    auto m = Semiconductor::create(silicon_parameters);
    NITCAD_EXPECTS(m.has_value());
    return *m;
}

std::expected<void, base::Error> check_temperature(const Semiconductor& m,
                                                   double temperature_K) {
    const double T = temperature_K;
    if (!(std::isfinite(T) && T > 0.0)) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input, "temperature must be finite and positive",
            base::ErrorContext{.index = std::nullopt, .value = T}});
    }
    const SemiconductorParameters& p = m.parameters();
    const double Eg = band_gap_unchecked(p, T);
    if (!(Eg > 0.0)) {  // also rejects NaN
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input, "band gap Eg(T) is not positive at this temperature",
            base::ErrorContext{.index = std::nullopt, .value = Eg}});
    }
    const struct {
        std::string_view message;
        const CaugheyThomasParameters& ct;
    } carriers[] = {
        {"electron_mobility: mu_max(T) is not finite or is below mu_min at this temperature",
         p.electron_mobility},
        {"hole_mobility: mu_max(T) is not finite or is below mu_min at this temperature",
         p.hole_mobility},
    };
    for (const auto& c : carriers) {
        const double mu_max = detail::mu_max_at(c.ct, T);
        if (!(std::isfinite(mu_max) && mu_max >= c.ct.mu_min)) {
            return std::unexpected(base::Error{
                base::ErrorCode::invalid_input, std::string{c.message},
                base::ErrorContext{.index = std::nullopt, .value = mu_max}});
        }
    }
    return {};
}

double band_gap_eV(const Semiconductor& m, double temperature_K) {
    detail::expect_temperature(temperature_K);
    return band_gap_unchecked(m.parameters(), temperature_K);
}

double conduction_band_dos(const Semiconductor& m, double temperature_K) {
    detail::expect_temperature(temperature_K);
    return m.parameters().Nc300 * dos_factor(temperature_K);
}

double valence_band_dos(const Semiconductor& m, double temperature_K) {
    detail::expect_temperature(temperature_K);
    return m.parameters().Nv300 * dos_factor(temperature_K);
}

std::expected<SemiconductorParameters, base::Error> algaas_parameters(double x,
                                                                     double conduction_share) {
    if (!(std::isfinite(x) && x >= 0.0 && x <= 0.45)) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input,
            "AlGaAs mole fraction must be in [0, 0.45] (direct-gap regime)",
            base::ErrorContext{.index = std::nullopt, .value = x}});
    }
    if (!(std::isfinite(conduction_share) && conduction_share >= 0.0 && conduction_share <= 1.0)) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input, "AlGaAs conduction-band share must be in [0, 1]",
            base::ErrorContext{.index = std::nullopt, .value = conduction_share}});
    }
    SemiconductorParameters p = gallium_arsenide_parameters;
    p.eps_r = 12.9 - 2.6 * x;
    p.electron_affinity_eV = 4.07 - conduction_share * 1.247 * x;
    p.Eg0_eV = 1.519 + 1.247 * x;
    p.varshni_alpha_eV_per_K = 5.405e-4;
    p.varshni_beta_K = 204.0;
    p.Nc300 = 4.7e17 * (1.0 + 0.5 * x);
    p.Nv300 = 7.0e18 * (1.0 - 0.3 * x);
    p.electron_mobility.mu_max = 8500.0 - 5500.0 * x;
    p.hole_mobility.mu_max = 400.0 - 150.0 * x;
    p.radiative_cm3_s = 1.8e-10;
    p.ionization.donor_eV = 0.0;
    p.ionization.acceptor_eV = 0.0;
    return p;
}

double intrinsic_level_depth_eV(const Semiconductor& m, double temperature_K) {
    detail::expect_temperature(temperature_K);
    const double T = temperature_K;
    const double ln_dos = std::log(conduction_band_dos(m, T) / valence_band_dos(m, T));
    const double conduction_to_intrinsic =
        0.5 * band_gap_eV(m, T) + 0.5 * base::thermal_voltage(T) * ln_dos;
    return m.parameters().electron_affinity_eV + conduction_to_intrinsic;
}

double intrinsic_density(const Semiconductor& m, double temperature_K) {
    detail::expect_temperature(temperature_K);
    const SemiconductorParameters& p = m.parameters();
    const double T = temperature_K;
    // kT in eV is numerically V_T = kT/q in V; use the one definition of it (base).
    const double kT_eV = base::thermal_voltage(T);
    return std::sqrt(p.Nc300 * p.Nv300) * dos_factor(T) *
           std::exp(-band_gap_unchecked(p, T) / (2.0 * kT_eV));
}

}  // namespace NiTCAD::physics
