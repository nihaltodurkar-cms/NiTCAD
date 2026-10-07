// Insulator regions, electrode contacts and semiconductor-insulator interfaces in the device
// description (ARCHITECTURE.md section 11, Unit 15b): each validation rule of Device::create, the
// topology rules with insulators, and the queries.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;

namespace {

// 1D: oxide nodes 0-2 (x < 0), silicon nodes 3-9; electrode on x_min, substrate on x_max.
std::vector<double> axis() {
    return {-3e-7, -2e-7, -1e-7, 1e-7, 2e-7, 3e-7, 4e-7, 5e-7, 6e-7, 7e-7};
}

device::DeviceDescription mos() {
    mesh::Mesh m = *mesh::make_tensor_grid(axis());
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 1);
    for (std::size_t i = 0; i < 3; ++i) region[i] = 0;
    std::vector<double> acceptors(n, 1e17);
    for (std::size_t i = 0; i < 3; ++i) acceptors[i] = 0.0;
    device::Contact gate{"gate", device::ContactKind::electrode, {0}};
    return {.mesh = std::move(m),
            .temperature_K = 300.0,
            .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
            .node_region = std::move(region),
            .donors = std::vector<double>(n, 0.0),
            .acceptors = std::move(acceptors),
            .contacts = {std::move(gate), {"substrate", device::ContactKind::ohmic, {9}}},
            .interfaces = {{"oxide", "silicon"}}};
}

bool rejected(device::DeviceDescription d, const std::string& fragment) {
    const auto r = device::Device::create(std::move(d));
    if (r) return false;
    CAPTURE(r.error().message);
    return r.error().code == base::ErrorCode::invalid_input &&
           r.error().message.find(fragment) != std::string::npos;
}

}  // namespace

TEST_CASE("insulator device: a meshed MOS stack is accepted and queried") {
    auto d = mos();
    d.interfaces[0].fixed_charge_cm2 = 5e10;
    d.interfaces[0].traps.levels = {{.density_cm2 = 1e10}};
    d.interfaces[0].traps.bands = {{.density_cm2_eV = 1e11, .energy_low_eV = -0.5,
                                    .energy_high_eV = 0.5}};
    d.interfaces[0].recombination_velocity_n_cm_s = 100.0;
    d.interfaces[0].recombination_velocity_p_cm_s = 100.0;
    const auto r = device::Device::create(std::move(d));
    REQUIRE(r.has_value());
    REQUIRE(r->is_insulator(0));
    REQUIRE(r->is_insulator(2));
    REQUIRE_FALSE(r->is_insulator(3));
    REQUIRE(r->relative_permittivity(1) == 3.9);
    REQUIRE(r->relative_permittivity(5) == 11.7);
    REQUIRE(r->reference_node() == 3);
    REQUIRE(r->find_interface(0, 1) == 0);
    REQUIRE(r->find_interface(1, 0) == 0);
    REQUIRE(device::has_interface_charge_or_recombination(r->interfaces()[0]));
    REQUIRE(r->material(3).parameters() == physics::silicon_parameters);
    // Without insulators the reference node is node 0.
    auto plain = mos();
    plain.regions[0] = {"oxide", physics::silicon()};
    plain.contacts[0].kind = device::ContactKind::ohmic;
    plain.interfaces.clear();
    REQUIRE(device::Device::create(std::move(plain))->reference_node() == 0);
}

TEST_CASE("insulator device: regions and doping") {
    {
        auto d = mos();
        d.regions[1] = {"silicon", physics::silicon_dioxide()};
        d.contacts[1].kind = device::ContactKind::electrode;
        d.acceptors.assign(d.acceptors.size(), 0.0);
        d.interfaces.clear();
        REQUIRE(rejected(std::move(d), "at least one semiconductor region"));
    }
    {
        auto d = mos();
        d.donors[1] = 1e15;
        REQUIRE(rejected(std::move(d), "insulator node has a nonzero donor or acceptor"));
    }
}

TEST_CASE("insulator device: contacts and insulators") {
    {
        auto d = mos();  // an ohmic contact in the oxide
        d.contacts[0].kind = device::ContactKind::ohmic;
        REQUIRE(rejected(std::move(d), "node is in an insulator"));
    }
    {
        auto d = mos();  // an electrode on the semiconductor
        d.contacts[1].kind = device::ContactKind::electrode;
        REQUIRE(rejected(std::move(d), "electrode node is not in an insulator"));
    }
    {
        auto d = mos();  // a lumped gate in the oxide
        d.contacts[0].kind = device::ContactKind::gate;
        d.contacts[0].gate = {.boundary = "x_min", .oxide_thickness_cm = 5e-7};
        REQUIRE(rejected(std::move(d), "node is in an insulator"));
    }
    {
        auto d = mos();
        d.contacts[0].electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 0.0};
        REQUIRE(rejected(std::move(d), "work function is not finite and positive"));
    }
    {
        auto d = mos();
        d.contacts[0].electrode.kind = static_cast<device::GateElectrode>(9);
        REQUIRE(rejected(std::move(d), "unknown electrode kind"));
    }
    {
        auto d = mos();
        d.contacts[0].kind = static_cast<device::ContactKind>(9);
        REQUIRE(rejected(std::move(d), "unknown contact kind"));
    }
    auto metal = mos();
    metal.contacts[0].electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.6};
    REQUIRE(device::Device::create(std::move(metal)).has_value());
}

TEST_CASE("insulator device: topology with insulators") {
    {
        // Oxide on both sides of a silicon island that has no ohmic contact: x = silicon (0-2),
        // oxide (3-5), silicon (6-9) with the substrate contact. The island floats.
        auto d = mos();
        d.node_region = {1, 1, 1, 0, 0, 0, 1, 1, 1, 1};
        d.acceptors = {1e17, 1e17, 1e17, 0, 0, 0, 1e17, 1e17, 1e17, 1e17};
        d.contacts = {{"substrate", device::ContactKind::ohmic, {9}}};
        const auto r = device::Device::create(std::move(d));
        REQUIRE_FALSE(r.has_value());
        REQUIRE(r.error().message.find("no ohmic contact") != std::string::npos);
        REQUIRE(r.error().context->index == 0u);
    }
    {
        // The same island with its own ohmic contact on x_min is accepted: carriers do not cross
        // the oxide, the potential is continuous through it.
        auto d = mos();
        d.node_region = {1, 1, 1, 0, 0, 0, 1, 1, 1, 1};
        d.acceptors = {1e17, 1e17, 1e17, 0, 0, 0, 1e17, 1e17, 1e17, 1e17};
        d.contacts = {{"left", device::ContactKind::ohmic, {0}},
                      {"substrate", device::ContactKind::ohmic, {9}}};
        REQUIRE(device::Device::create(std::move(d)).has_value());
    }
    {
        // An oxide joined to the rest by no edge: two meshes in one, the oxide part has neither an
        // ohmic contact nor an electrode. Built from parts: nodes 0-1 oxide (one edge), 2-3
        // silicon (one edge).
        mesh::Mesh m = *mesh::Mesh::from_parts(
            1, {{0, 0, 0}, {1e-7, 0, 0}, {5e-7, 0, 0}, {6e-7, 0, 0}},
            {0.5e-7, 0.5e-7, 0.5e-7, 0.5e-7},
            {{0, 1, 1e-7, 1.0}, {2, 3, 1e-7, 1.0}},
            {{"left", {0}, {1.0}}, {"right", {3}, {1.0}}});
        device::DeviceDescription d{
            .mesh = std::move(m),
            .temperature_K = 300.0,
            .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
            .node_region = {0, 0, 1, 1},
            .donors = {0, 0, 0, 0},
            .acceptors = {0, 0, 1e17, 1e17},
            .contacts = {{"substrate", device::ContactKind::ohmic, {3}}}};
        const auto r = device::Device::create(d);
        REQUIRE_FALSE(r.has_value());
        REQUIRE(r.error().message.find("no ohmic contact or electrode") != std::string::npos);
        REQUIRE(r.error().context->index == 0u);
        d.contacts.push_back({"gate", device::ContactKind::electrode, {0}});
        REQUIRE(device::Device::create(std::move(d)).has_value());
    }
}

TEST_CASE("insulator device: interface rules") {
    {
        auto d = mos();
        d.interfaces[0].transport = device::InterfaceTransport::thermionic_emission;
        REQUIRE(rejected(std::move(d), "thermionic emission needs two semiconductor regions"));
    }
    {
        // Interface charge between two semiconductors.
        auto d = mos();
        d.regions[0] = {"oxide", physics::silicon()};
        d.contacts[0].kind = device::ContactKind::ohmic;
        d.interfaces[0].fixed_charge_cm2 = 1e10;
        REQUIRE(rejected(std::move(d), "need a semiconductor-insulator interface"));
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    {
        auto d = mos();
        d.interfaces[0].fixed_charge_cm2 = nan;
        REQUIRE(rejected(std::move(d), "fixed charge is not finite"));
    }
    {
        auto d = mos();
        d.interfaces[0].recombination_velocity_n_cm_s = -1.0;
        REQUIRE(rejected(std::move(d), "recombination velocity is not finite and >= 0"));
    }
    {
        auto d = mos();
        d.interfaces[0].recombination_velocity_p_cm_s = nan;
        REQUIRE(rejected(std::move(d), "recombination velocity is not finite and >= 0"));
    }
    {
        auto d = mos();
        d.interfaces[0].traps.levels = {{.density_cm2 = 1e10, .energy_eV = 0.7}};
        REQUIRE(rejected(std::move(d), "trap level energy is not inside the gap"));
    }
    {
        // The gap checked is the semiconductor side's, in either declaration order.
        auto d = mos();
        std::swap(d.interfaces[0].region_a, d.interfaces[0].region_b);
        d.interfaces[0].traps.bands = {{.density_cm2_eV = 1e11, .energy_low_eV = -0.2,
                                        .energy_high_eV = 0.6}};
        REQUIRE(rejected(std::move(d), "trap band upper edge"));
    }
    {
        // A negative fixed charge is allowed.
        auto d = mos();
        d.interfaces[0].fixed_charge_cm2 = -3e12;
        REQUIRE(device::Device::create(std::move(d)).has_value());
    }
}
