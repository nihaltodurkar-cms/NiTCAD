#include "NiTCAD/physics/semiconductor.hpp"

#include <cmath>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
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
