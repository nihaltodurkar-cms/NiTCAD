#include "NiTCAD/physics/interface_traps.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

namespace {

base::Error invalid(std::string message, std::size_t index, double value) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = value}};
}

bool positive(double v) { return std::isfinite(v) && v > 0.0; }
bool non_negative(double v) { return std::isfinite(v) && v >= 0.0; }

// 6-point Gauss-Legendre rule on [-1, 1] (abscissae and weights to 25 digits).
constexpr double gl_x[3] = {0.2386191860831969086305017, 0.6612093864662645136613996,
                            0.9324695142031520278123016};
constexpr double gl_w[3] = {0.4679139345726910473898703, 0.3607615730481386075698335,
                            0.1713244923791703450402961};

}  // namespace

std::expected<void, base::Error> check_interface_traps(const InterfaceTraps& traps,
                                                       const Semiconductor& m,
                                                       double temperature_K) {
    if (!positive(traps.thermal_velocity_n_cm_s)) {
        return std::unexpected(invalid("trap electron thermal velocity is not finite and positive",
                                       0, traps.thermal_velocity_n_cm_s));
    }
    if (!positive(traps.thermal_velocity_p_cm_s)) {
        return std::unexpected(invalid("trap hole thermal velocity is not finite and positive", 0,
                                       traps.thermal_velocity_p_cm_s));
    }
    // E_c - E_i and E_v - E_i of the semiconductor.
    const double ec =
        intrinsic_level_depth_eV(m, temperature_K) - m.parameters().electron_affinity_eV;
    const double ev = ec - band_gap_eV(m, temperature_K);
    const auto in_gap = [&](double e) { return std::isfinite(e) && e >= ev && e <= ec; };
    const auto cross_sections = [&](double sn, double sp, std::size_t k,
                                    const char* what) -> std::optional<base::Error> {
        const std::string where = std::string("trap ") + what;
        if (!positive(sn)) {
            return invalid(where + " electron cross-section is not finite and positive", k, sn);
        }
        if (!positive(sp)) {
            return invalid(where + " hole cross-section is not finite and positive", k, sp);
        }
        return std::nullopt;
    };
    for (std::size_t k = 0; k < traps.levels.size(); ++k) {
        const TrapLevel& t = traps.levels[k];
        if (t.type != TrapType::donor && t.type != TrapType::acceptor) {
            return std::unexpected(invalid("trap level has an unknown type", k, 0.0));
        }
        if (!non_negative(t.density_cm2)) {
            return std::unexpected(
                invalid("trap level density is not finite and >= 0", k, t.density_cm2));
        }
        if (!in_gap(t.energy_eV)) {
            return std::unexpected(
                invalid("trap level energy is not inside the gap", k, t.energy_eV));
        }
        if (auto e = cross_sections(t.sigma_n_cm2, t.sigma_p_cm2, k, "level")) {
            return std::unexpected(std::move(*e));
        }
    }
    for (std::size_t k = 0; k < traps.bands.size(); ++k) {
        const TrapBand& b = traps.bands[k];
        if (b.type != TrapType::donor && b.type != TrapType::acceptor) {
            return std::unexpected(invalid("trap band has an unknown type", k, 0.0));
        }
        if (!non_negative(b.density_cm2_eV)) {
            return std::unexpected(
                invalid("trap band density is not finite and >= 0", k, b.density_cm2_eV));
        }
        if (!in_gap(b.energy_low_eV)) {
            return std::unexpected(
                invalid("trap band lower edge is not inside the gap", k, b.energy_low_eV));
        }
        if (!in_gap(b.energy_high_eV) || !(b.energy_high_eV > b.energy_low_eV)) {
            return std::unexpected(invalid(
                "trap band upper edge is not inside the gap and above the lower edge", k,
                b.energy_high_eV));
        }
        if (auto e = cross_sections(b.sigma_n_cm2, b.sigma_p_cm2, k, "band")) {
            return std::unexpected(std::move(*e));
        }
    }
    return {};
}

std::vector<TrapLevel> trap_band_levels(const TrapBand& band, double temperature_K) {
    NITCAD_EXPECTS(std::isfinite(temperature_K) && temperature_K > 0.0);
    NITCAD_EXPECTS(band.energy_high_eV > band.energy_low_eV);
    const double kT = base::thermal_voltage(temperature_K);  // [eV]
    const double width = band.energy_high_eV - band.energy_low_eV;
    const auto panels = static_cast<std::size_t>(std::ceil(width / kT));
    const double h = width / static_cast<double>(panels);
    std::vector<TrapLevel> levels;
    levels.reserve(6 * panels);
    for (std::size_t k = 0; k < panels; ++k) {
        const double mid = band.energy_low_eV + (static_cast<double>(k) + 0.5) * h;
        for (std::size_t q = 0; q < 3; ++q) {
            for (const double sign : {-1.0, 1.0}) {
                levels.push_back({.type = band.type,
                                  .density_cm2 = band.density_cm2_eV * 0.5 * h * gl_w[q],
                                  .energy_eV = mid + sign * 0.5 * h * gl_x[q],
                                  .sigma_n_cm2 = band.sigma_n_cm2,
                                  .sigma_p_cm2 = band.sigma_p_cm2});
            }
        }
    }
    return levels;
}

FermiOccupancy fermi_occupancy(double x) noexcept {
    // Each of f and 1 - f from the exponential that cannot overflow.
    double f, g;
    if (x > 0.0) {
        const double e = std::exp(-x);
        f = e / (1.0 + e);
        g = 1.0 / (1.0 + e);
    } else {
        const double e = std::exp(x);
        f = 1.0 / (1.0 + e);
        g = e / (1.0 + e);
    }
    return {f, g, f * g};
}

}  // namespace NiTCAD::physics
