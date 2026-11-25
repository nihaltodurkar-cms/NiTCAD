// Klaassen mobility in a solved device (Unit 66; DECISIONS.md KC6): the Unit 9 diode swept to
// 0.8 V with klaassen_equilibrium converges at every bias; its currents and Newton steps are
// printed beside Masetti's for the record.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

// The Unit 9 diode (default_physics_test.cpp): legacy graded mesh, 1e17 / 1e17, silicon.
device::Device diode() {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

}  // namespace

TEST_CASE("klaassen: the Unit 9 diode IV converges with klaassen_equilibrium (KC6)",
          "[klaassen]") {
    const device::Device d = diode();
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 16; ++k) points.push_back({0.05 * k, 0.0});
    const solve::BiasOptions masetti{
        .models = {.hole_mobility = assemble::HoleMobility::masetti,
                   .electron_mobility = assemble::ElectronMobility::masetti}};
    const solve::BiasOptions klaassen{
        .models = {.hole_mobility = assemble::HoleMobility::klaassen_equilibrium,
                   .electron_mobility = assemble::ElectronMobility::klaassen_equilibrium}};
    const auto m = solve::sweep_bias(d, points, masetti);
    const auto k = solve::sweep_bias(d, points, klaassen);
    REQUIRE(m.has_value());
    REQUIRE(k.has_value());
    REQUIRE_FALSE(m->stopped.has_value());
    REQUIRE_FALSE(k->stopped.has_value());
    REQUIRE(k->points.size() == points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& pm = m->points[i];
        const auto& pk = k->points[i];
        std::printf("KC6 V=%.2f  J masetti %+.5e  klaassen_equilibrium %+.5e A/cm^2 (%+.2f%%)  "
                    "Newton %zu / %zu\n",
                    points[i][0], pm.terminal_current[0], pk.terminal_current[0],
                    pm.terminal_current[0] != 0.0
                        ? 100.0 * (pk.terminal_current[0] / pm.terminal_current[0] - 1.0)
                        : 0.0,
                    pm.convergence.iterations.size(), pk.convergence.iterations.size());
        CHECK(pk.convergence.converged);
        CHECK(std::isfinite(pk.terminal_current[0]));
    }
}
