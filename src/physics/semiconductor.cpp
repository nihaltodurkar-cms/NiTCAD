#include "NiTCAD/physics/semiconductor.hpp"

#include <cmath>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

namespace {

enum class Bound { positive, non_negative };

// The first failing check, as an invalid_input error naming the parameter.
std::optional<base::Error> check(std::string_view name, double value, Bound bound) {
    const bool ok = std::isfinite(value) && (bound == Bound::positive ? value > 0.0 : value >= 0.0);
    if (ok) return std::nullopt;
    std::string message{"semiconductor parameter "};
    message += name;
    message += bound == Bound::positive ? " must be finite and positive"
                                        : " must be finite and non-negative";
    return base::Error{base::ErrorCode::invalid_input, std::move(message),
                       base::ErrorContext{.index = std::nullopt, .value = value}};
}

std::optional<base::Error> check_mobility(std::string_view carrier,
                                          const CaugheyThomasParameters& ct) {
    const std::string prefix = std::string{carrier} + "_mobility.";
    if (auto e = check(prefix + "mu_min", ct.mu_min, Bound::non_negative)) return e;
    if (auto e = check(prefix + "mu_max", ct.mu_max, Bound::positive)) return e;
    if (auto e = check(prefix + "N_ref", ct.N_ref, Bound::positive)) return e;
    if (auto e = check(prefix + "alpha", ct.alpha, Bound::positive)) return e;
    if (!std::isfinite(ct.T_exponent)) {
        return base::Error{base::ErrorCode::invalid_input,
                           "semiconductor parameter " + prefix + "T_exponent must be finite",
                           base::ErrorContext{.index = std::nullopt, .value = ct.T_exponent}};
    }
    if (ct.mu_min > ct.mu_max) {
        // Mobility would rise with doping.
        return base::Error{base::ErrorCode::invalid_input,
                           "semiconductor parameter " + prefix + "mu_min exceeds mu_max",
                           base::ErrorContext{.index = std::nullopt, .value = ct.mu_min}};
    }
    return std::nullopt;
}

void expect_temperature(double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(temperature_K) && temperature_K > 0.0);
}

}  // namespace

std::expected<Semiconductor, base::Error> Semiconductor::create(
    const SemiconductorParameters& p) {
    const std::optional<base::Error> failures[] = {
        check("eps_r", p.eps_r, Bound::positive),
        check("Eg0_eV", p.Eg0_eV, Bound::positive),
        check("varshni_alpha_eV_per_K", p.varshni_alpha_eV_per_K, Bound::non_negative),
        check("varshni_beta_K", p.varshni_beta_K, Bound::non_negative),
        check("Nc300", p.Nc300, Bound::positive),
        check("Nv300", p.Nv300, Bound::positive),
        check_mobility("electron", p.electron_mobility),
        check_mobility("hole", p.hole_mobility),
        check("lifetime.tau_n0", p.lifetime.tau_n0, Bound::positive),
        check("lifetime.tau_p0", p.lifetime.tau_p0, Bound::positive),
        check("lifetime.N_ref", p.lifetime.N_ref, Bound::positive),
    };
    for (const auto& failure : failures) {
        if (failure) return std::unexpected(*failure);
    }
    return Semiconductor{p};
}

Semiconductor silicon() {
    auto m = Semiconductor::create(silicon_parameters);
    NITCAD_EXPECTS(m.has_value());
    return *m;
}

double band_gap_eV(const Semiconductor& m, double temperature_K) {
    expect_temperature(temperature_K);
    const SemiconductorParameters& p = m.parameters();
    const double T = temperature_K;
    return p.Eg0_eV - p.varshni_alpha_eV_per_K * (T * T) / (T + p.varshni_beta_K);
}

double conduction_band_dos(const Semiconductor& m, double temperature_K) {
    expect_temperature(temperature_K);
    return m.parameters().Nc300 * std::pow(temperature_K / 300.0, 1.5);
}

double valence_band_dos(const Semiconductor& m, double temperature_K) {
    expect_temperature(temperature_K);
    return m.parameters().Nv300 * std::pow(temperature_K / 300.0, 1.5);
}

double intrinsic_density(const Semiconductor& m, double temperature_K) {
    const double kT_eV = base::k_B_eV_per_K * temperature_K;
    return std::sqrt(conduction_band_dos(m, temperature_K) * valence_band_dos(m, temperature_K)) *
           std::exp(-band_gap_eV(m, temperature_K) / (2.0 * kT_eV));
}

}  // namespace NiTCAD::physics
