// Shared fixtures of the deck tests: a 1D diode as a deck and built by hand.
#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"

namespace deck_test {

inline std::filesystem::path scratch(const std::string& name) {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "nitcad_deck_test" / name;
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    return dir;
}

inline std::expected<NiTCAD::deck::Plan, NiTCAD::base::Error> plan_of(const std::string& text,
                                                                      const std::filesystem::path& dir = {}) {
    const auto d = NiTCAD::deck::parse(text, "t.deck", dir);
    if (!d) return std::unexpected(d.error());
    return NiTCAD::deck::build(*d);
}

// The plan of `text`, which must build.
inline NiTCAD::deck::Plan must(const std::string& text, const std::filesystem::path& dir = {}) {
    auto p = plan_of(text, dir);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

// linspace as build.hpp states it: a + (b - a) k / (n - 1), the last value b.
inline std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) {
        v.push_back(k == n - 1 ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(n - 1));
    }
    return v;
}

// The diode without its analysis section.
inline const std::string diode_device = R"(
# The 1D diode of the deck tests.
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 2um, 201)

[region silicon]
material = silicon

[doping]
profile = acceptor constant 1e17 x=-1:1um
profile = donor constant 1e17 x=1um:3um

[contact anode]
kind = ohmic
patch = x_min

[contact cathode]
kind = ohmic
patch = x_max
)";

inline const std::string diode_deck = diode_device + R"(
[sweep]
contact anode = 0 : 0.2 : 0.1
)";

inline NiTCAD::device::Device hand_diode() {
    using namespace NiTCAD;
    auto m = *mesh::make_tensor_grid(linspace(0.0, 2e-4, 201));
    const std::size_t n = m.node_count();
    std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    const std::vector<device::RegionId> node_region(n, 0);
    device::DopingProfile a{.dopant = device::Dopant::acceptor, .peak_cm3 = 1e17};
    a.lo[0] = -1.0;
    a.hi[0] = 1e-4;
    device::DopingProfile d{.dopant = device::Dopant::donor, .peak_cm3 = 1e17};
    d.lo[0] = 1e-4;
    d.hi[0] = 3e-4;
    const std::vector<device::DopingProfile> profiles{a, d};
    auto doping = *device::sample(profiles, m, node_region, regions);
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    auto cells = *mesh::TensorCells::create(m, linspace(0.0, 2e-4, 201));
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = std::move(regions),
                                    .node_region = node_region,
                                    .donors = std::move(doping.donors),
                                    .acceptors = std::move(doping.acceptors),
                                    .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                 {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
                                    .cells = std::move(cells)});
}

}  // namespace deck_test
