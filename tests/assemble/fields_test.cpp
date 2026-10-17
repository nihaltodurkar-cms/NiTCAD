// Fields and the unknown map (Unit 40 slice 1; DECISIONS.md FE1, FG7): the drift-diffusion
// system's map reproduces the node-interleaved layout (psi, n, p, then the temperature rise with
// the electrothermal model), and every unknown's field and node can be asked for.
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::FieldId;

namespace {

device::Device resistor(bool thermal) {
    std::vector<double> x(11);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = 1e-5 * static_cast<double>(i);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    std::vector<device::ThermalContact> sinks;
    if (thermal) {
        sinks = {{"left", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
                 {"right", "x_max", device::ThermalContactKind::isothermal, 300.0, 0.0}};
    }
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::vector<double>(n, 1e16),
                                    .acceptors = std::vector<double>(n, 0.0),
                                    .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                 {"right", device::ContactKind::ohmic, std::move(right)}},
                                    .thermal_contacts = std::move(sinks)});
}

}  // namespace

TEST_CASE("fields: the drift-diffusion map is the node-interleaved layout") {
    for (const bool thermal : {false, true}) {
        CAPTURE(thermal);
        const device::Device d = resistor(thermal);
        const auto scaling = *assemble::make_scaling(d);
        const auto system =
            *assemble::DriftDiffusion::create(d, scaling, {.electrothermal = thermal});
        const assemble::DofMap& dofs = system.dofs();
        const std::size_t m = thermal ? 4 : 3;
        REQUIRE(dofs.block() == m);
        REQUIRE(system.stride() == m);
        REQUIRE(dofs.node_count() == system.node_count());
        REQUIRE(dofs.unknowns() == system.unknowns());
        REQUIRE(dofs.has(FieldId::temperature) == thermal);
        const FieldId order[] = {FieldId::potential, FieldId::electrons, FieldId::holes,
                                 FieldId::temperature};
        for (std::size_t k = 0; k < m; ++k) {
            REQUIRE(dofs.slot(order[k]) == k);
            REQUIRE(dofs.fields()[k].id == order[k]);
            REQUIRE(dofs.fields()[k].unknown);
            REQUIRE(dofs.fields()[k].location == assemble::FieldLocation::node);
        }
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            for (std::size_t k = 0; k < m; ++k) {
                const std::size_t u = dofs.index(order[k], i);
                REQUIRE(u == m * i + k);
                REQUIRE(dofs.field_of(u).id == order[k]);
                REQUIRE(dofs.node_of(u) == i);
            }
        }
    }
}

TEST_CASE("fields: the descriptors") {
    REQUIRE(assemble::potential_field.domain == assemble::FieldDomain::all);
    REQUIRE(assemble::electron_field.domain == assemble::FieldDomain::semiconductor);
    REQUIRE(assemble::hole_field.domain == assemble::FieldDomain::semiconductor);
    REQUIRE(assemble::temperature_field.domain == assemble::FieldDomain::all);
    const assemble::Field two[] = {assemble::hole_field, assemble::potential_field};
    const auto map = assemble::DofMap::node_interleaved(two, 5);
    REQUIRE(map.unknowns() == 10);
    REQUIRE(map.slot(FieldId::holes) == 0);
    REQUIRE(map.slot(FieldId::potential) == 1);
    REQUIRE(!map.slot(FieldId::electrons).has_value());
    REQUIRE(map.index(FieldId::potential, 3) == 7);
    REQUIRE(map.field_of(7).id == FieldId::potential);
    REQUIRE(map.node_of(7) == 3);
}
