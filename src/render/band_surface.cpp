// The band surface (band_surface.hpp).
#include "NiTCAD/render/band_surface.hpp"

#include <variant>

namespace NiTCAD::render {

std::vector<std::int32_t> band_material_keys(const device::Device& device) {
    // A key per region: insulators -1; a semiconductor the key of the first earlier region with equal
    // parameters (the solver's heterointerface test), else a new one.
    const auto regions = device.regions();
    std::vector<std::int32_t> region_key(regions.size(), -1);
    std::int32_t next = 0;
    for (std::size_t r = 0; r < regions.size(); ++r) {
        if (device::is_insulator(regions[r])) continue;
        const auto& params = std::get<physics::Semiconductor>(regions[r].material).parameters();
        for (std::size_t q = 0; q < r && region_key[r] < 0; ++q) {
            if (region_key[q] >= 0 && std::get<physics::Semiconductor>(regions[q].material).parameters() == params) {
                region_key[r] = region_key[q];
            }
        }
        if (region_key[r] < 0) region_key[r] = next++;
    }
    std::vector<std::int32_t> keys;
    keys.reserve(device.node_region().size());
    for (const auto r : device.node_region()) keys.push_back(region_key[static_cast<std::size_t>(r)]);
    return keys;
}

double band_value(const Geometry& geometry, std::size_t element, const Point& x, NodeId owner,
                  std::span<const std::int32_t> keys, std::span<const double> level) {
    const Interpolant at = geometry.interpolant(element, x);
    const std::int32_t key = keys[static_cast<std::size_t>(owner)];
    bool one = true;
    double weight = 0.0, sum = 0.0;
    for (std::size_t k = 0; k < at.nodes.size(); ++k) {
        if (at.weights[k] == 0.0) continue;
        const auto n = static_cast<std::size_t>(at.nodes[k]);
        if (keys[n] != key) {
            one = false;
            continue;
        }
        weight += at.weights[k];
        sum += at.weights[k] * level[n];
    }
    if (one) return interpolate(at, level);
    if (weight > 0.0) return sum / weight;
    return level[static_cast<std::size_t>(owner)];
}

}  // namespace NiTCAD::render
