// Private to the solve layer: the electrothermal model's input checks and the converged
// temperature's range (Unit 23; DECISIONS.md T7, T10, T12), shared by the solve paths.
#pragma once

#include <cmath>
#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <string>
#include <variant>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::solve::detail {

// Whether temperature T [K] is usable by every semiconductor region's models; the first
// region's error otherwise.
inline std::optional<base::Error> temperature_range(const device::Device& device, double T) {
    for (const device::Region& r : device.regions()) {
        if (device::is_insulator(r)) continue;
        auto ok = physics::check_temperature(std::get<physics::Semiconductor>(r.material), T);
        if (!ok) {
            return base::Error{base::ErrorCode::invalid_input,
                               "region '" + r.name + "': " + ok.error().message, {}};
        }
    }
    return std::nullopt;
}

// The thermal contacts' temperatures of each of `points` bias points: `given` (one list per point,
// one temperature per thermal contact, each finite, positive and in every semiconductor's range)
// or, empty, the device's own at every point. Errors (invalid_input) name the point.
inline std::expected<std::vector<std::vector<double>>, base::Error> thermal_points(
    const device::Device& device, std::span<const std::vector<double>> given,
    std::size_t points) {
    const auto contacts = device.thermal_contacts();
    if (given.empty()) {
        std::vector<double> own;
        for (const device::ThermalContact& c : contacts) own.push_back(c.temperature_K);
        return std::vector<std::vector<double>>(points, own);
    }
    const auto invalid = [](std::string message, std::size_t k, std::optional<double> value) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input, std::move(message),
                                           base::ErrorContext{.index = k, .value = value}});
    };
    if (given.size() != points) {
        return invalid("thermal_bias_K needs one list per bias point", 0,
                       static_cast<double>(given.size()));
    }
    for (std::size_t k = 0; k < points; ++k) {
        const std::string where = "bias point " + std::to_string(k);
        if (given[k].size() != contacts.size()) {
            return invalid(where + ": one temperature per thermal contact is needed", k,
                           static_cast<double>(given[k].size()));
        }
        for (std::size_t c = 0; c < contacts.size(); ++c) {
            const double T = given[k][c];
            const std::string what = where + ", thermal contact '" + contacts[c].name + "': ";
            if (!(std::isfinite(T) && T > 0.0)) {
                return invalid(what + "temperature must be finite and positive", k, T);
            }
            if (auto e = temperature_range(device, T)) return invalid(what + e->message, k, T);
        }
    }
    return std::vector<std::vector<double>>(given.begin(), given.end());
}

// An initial state's temperature: empty, or one finite, positive value per node.
inline bool valid_temperature(const results::NodeFields& s, std::size_t nodes) {
    if (s.temperature_K.empty()) return true;
    if (s.temperature_K.size() != nodes) return false;
    for (const double T : s.temperature_K) {
        if (!(std::isfinite(T) && T > 0.0)) return false;
    }
    return true;
}

// T10: a converged temperature outside a node's material range (non_convergence, context: the
// node and its temperature [K]). Insulator nodes are not checked (their models have no range).
inline std::optional<base::Error> converged_temperature(const device::Device& device,
                                                        std::span<const double> x, double T0) {
    for (std::size_t i = 0; i < device.mesh().node_count(); ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        if (device.is_insulator(node)) continue;
        const double T = T0 + T0 * x[4 * i + 3];  // the rise tau = (T - T0) / T0
        if (auto ok = physics::check_temperature(device.material(node), T); !ok) {
            return base::Error{base::ErrorCode::non_convergence,
                               "converged lattice temperature at node " + std::to_string(i) +
                                   " is outside its material's range: " + ok.error().message,
                               base::ErrorContext{.index = i, .value = T}};
        }
    }
    return std::nullopt;
}

}  // namespace NiTCAD::solve::detail
