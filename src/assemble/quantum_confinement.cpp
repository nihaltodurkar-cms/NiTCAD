#include "NiTCAD/assemble/quantum_confinement.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <variant>

#include "NiTCAD/physics/quantum_confinement.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::assemble {

namespace {

base::Error refuse(std::string message) {
    return {base::ErrorCode::invalid_input, "quantum_confinement mlda: " + std::move(message), {}};
}

}  // namespace

std::expected<ConfinementFactors, base::Error> mlda_factors(const device::Device& device,
                                                            const PhysicsModels& models) {
    for (const double m : {models.mlda_mass_n, models.mlda_mass_p}) {
        if (!(std::isfinite(m) && m > 0.0)) {
            return std::unexpected(refuse("the masses must be finite and positive"));
        }
    }
    if (models.electrothermal) {
        return std::unexpected(refuse("electrothermal is not supported (lambda is taken at the "
                                      "device temperature)"));
    }
    for (const device::Region& r : device.regions()) {
        const auto* s = std::get_if<physics::Semiconductor>(&r.material);
        if (s && !(s->parameters().electron_mobility ==
                       physics::silicon_parameters.electron_mobility &&
                   s->parameters().hole_mobility == physics::silicon_parameters.hole_mobility)) {
            return std::unexpected(refuse("region '" + r.name +
                                          "' is not silicon (the masses are silicon's)"));
        }
    }
    const mesh::Mesh& m = device.mesh();
    const std::size_t n = m.node_count();
    const auto node_region = device.node_region();
    std::vector<char> insulator(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        insulator[i] = device.is_insulator(static_cast<mesh::NodeId>(i)) ? 1 : 0;
    }

    ConfinementFactors c;
    const double T = device.temperature_K();
    c.wavelength_n_cm = physics::mlda_wavelength_cm(models.mlda_mass_n, T);
    c.wavelength_p_cm = physics::mlda_wavelength_cm(models.mlda_mass_p, T);
    c.distance_cm.assign(n, 0.0);
    c.half_cell.assign(n, 0);
    c.f_n.assign(n, 1.0);
    c.f_p.assign(n, 1.0);
    for (const device::Contact& contact : device.contacts()) {
        if (contact.kind != device::ContactKind::gate) continue;
        for (const mesh::NodeId v : contact.nodes) c.half_cell[static_cast<std::size_t>(v)] = 1;
    }

    // The interface points, and each lumped-gate node's shortest semiconductor edge.
    std::vector<mesh::Point> points;
    std::vector<double> shortest(n, std::numeric_limits<double>::infinity());
    for (const mesh::Edge& e : m.edges()) {
        const auto i = static_cast<std::size_t>(e.first);
        const auto j = static_cast<std::size_t>(e.second);
        const device::RegionId ra = node_region[i], rb = node_region[j];
        if (insulator[i] == 0 && insulator[j] == 0) {
            if (ra != rb &&
                device.transport(ra, rb) == device::InterfaceTransport::thermionic_emission) {
                return std::unexpected(refuse("a thermionic-emission interface is not supported"));
            }
            shortest[i] = std::min(shortest[i], e.length);
            shortest[j] = std::min(shortest[j], e.length);
        } else if (insulator[i] != insulator[j]) {
            const mesh::Point& a = m.points()[i];
            const mesh::Point& b = m.points()[j];
            points.push_back({0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])});
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (c.half_cell[i] != 0) points.push_back(m.points()[i]);
    }
    if (points.empty()) {
        return std::unexpected(refuse("the device has no insulator and no lumped gate"));
    }

    const int D = m.dimension();
    for (std::size_t i = 0; i < n; ++i) {
        if (insulator[i] != 0) continue;
        if (c.half_cell[i] != 0) {
            if (!std::isfinite(shortest[i])) {
                return std::unexpected(refuse("a lumped-gate node has no semiconductor edge"));
            }
            c.f_n[i] = physics::mlda_mean_factor(0.5 * shortest[i], c.wavelength_n_cm);
            c.f_p[i] = physics::mlda_mean_factor(0.5 * shortest[i], c.wavelength_p_cm);
            continue;
        }
        const mesh::Point& x = m.points()[i];
        double best = std::numeric_limits<double>::infinity();
        for (const mesh::Point& p : points) {
            double d2 = 0.0;
            for (int d = 0; d < D; ++d) d2 += (x[d] - p[d]) * (x[d] - p[d]);
            best = std::min(best, d2);
        }
        c.distance_cm[i] = std::sqrt(best);
        c.f_n[i] = physics::mlda_factor(c.distance_cm[i], c.wavelength_n_cm);
        c.f_p[i] = physics::mlda_factor(c.distance_cm[i], c.wavelength_p_cm);
    }
    return c;
}

}  // namespace NiTCAD::assemble
