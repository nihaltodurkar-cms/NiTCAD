// Device construction and validation (ARCHITECTURE.md section 11, Unit 6): a valid description is
// read back unchanged, and every invalid input returns an error naming what is wrong.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::device;
using NiTCAD::base::ErrorCode;
using NiTCAD::mesh::Mesh;
using NiTCAD::mesh::NodeId;
using NiTCAD::physics::silicon;

namespace {

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();
constexpr double infinity = std::numeric_limits<double>::infinity();

std::vector<double> axis(double length, int nodes) {
    std::vector<double> x(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) x[static_cast<std::size_t>(i)] = length * i / (nodes - 1);
    return x;
}

std::vector<NodeId> patch_nodes(const Mesh& m, const char* name) {
    return m.find_boundary(name)->nodes;
}

// A silicon pn diode, 2 um long, 1e17 / 1e17, junction at 1 um, contacts on x_min (anode, p side)
// and x_max (cathode, n side). 1D by default, 2D with ny > 0.
DeviceDescription diode(int ny = 0) {
    const auto x = axis(2e-4, 21);
    Mesh m = ny > 0 ? *NiTCAD::mesh::make_tensor_grid(x, axis(1e-4, ny))
                    : *NiTCAD::mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0);
    std::vector<double> acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    }
    std::vector<Contact> contacts{{"anode", ContactKind::ohmic, patch_nodes(m, "x_min")},
                                  {"cathode", ContactKind::ohmic, patch_nodes(m, "x_max")}};
    return DeviceDescription{.mesh = std::move(m),
                             .temperature_K = 300.0,
                             .regions = {{"silicon", silicon()}},
                             .node_region = std::vector<RegionId>(n, 0),
                             .donors = std::move(donors),
                             .acceptors = std::move(acceptors),
                             .contacts = std::move(contacts)};
}

// Applies spoil to a valid diode and returns the error, which must be invalid_input.
NiTCAD::base::Error rejected(const std::function<void(DeviceDescription&)>& spoil, int ny = 0) {
    DeviceDescription d = diode(ny);
    spoil(d);
    auto device = Device::create(std::move(d));
    REQUIRE_FALSE(device.has_value());
    REQUIRE(device.error().code == ErrorCode::invalid_input);
    return device.error();
}

bool mentions(const NiTCAD::base::Error& e, const std::string& text) {
    return e.message.find(text) != std::string::npos;
}

}  // namespace

TEST_CASE("device: a valid 1D diode is accepted and read back unchanged") {
    const auto d = Device::create(diode());
    REQUIRE(d.has_value());
    REQUIRE(d->mesh().node_count() == 21);
    REQUIRE(d->temperature_K() == 300.0);
    REQUIRE(d->regions().size() == 1);
    REQUIRE(d->regions()[0].name == "silicon");
    REQUIRE(d->contacts().size() == 2);
    REQUIRE(d->find_contact("anode")->nodes == std::vector<NodeId>{0});
    REQUIRE(d->find_contact("cathode")->nodes == std::vector<NodeId>{20});
    REQUIRE(d->find_contact("gate") == nullptr);
    REQUIRE(d->net_doping(0) == -1e17);
    REQUIRE(d->net_doping(20) == 1e17);
    REQUIRE(d->total_impurity(0) == 1e17);
    REQUIRE(d->material(5).parameters().eps_r == 11.7);
}

TEST_CASE("device: a 2D diode with line contacts is accepted") {
    const auto d = Device::create(diode(6));
    REQUIRE(d.has_value());
    REQUIRE(d->mesh().dimension() == 2);
    REQUIRE(d->find_contact("anode")->nodes.size() == 6);
    REQUIRE(d->find_contact("cathode")->nodes.size() == 6);
}

TEST_CASE("device: compensated doping gives net and total separately") {
    DeviceDescription desc = diode();
    desc.donors[3] = 3e17;  // the node already has 1e17 acceptors
    const auto d = Device::create(std::move(desc));
    REQUIRE(d.has_value());
    REQUIRE(d->net_doping(3) == 2e17);
    REQUIRE(d->total_impurity(3) == 4e17);
}

TEST_CASE("device: several regions, each with its own material") {
    DeviceDescription desc = diode();
    NiTCAD::physics::SemiconductorParameters p = NiTCAD::physics::silicon_parameters;
    p.eps_r = 12.0;
    desc.regions.push_back({"other", *NiTCAD::physics::Semiconductor::create(p)});
    for (std::size_t i = 10; i < desc.node_region.size(); ++i) desc.node_region[i] = 1;
    const auto d = Device::create(std::move(desc));
    REQUIRE(d.has_value());
    REQUIRE(d->material(9).parameters().eps_r == 11.7);
    REQUIRE(d->material(10).parameters().eps_r == 12.0);
}

TEST_CASE("device: regions are validated") {
    REQUIRE(mentions(rejected([](DeviceDescription& d) { d.regions.clear(); }), "region"));
    const auto unnamed = rejected([](DeviceDescription& d) { d.regions[0].name.clear(); });
    REQUIRE(unnamed.context->index == 0);
    const auto repeated = rejected([](DeviceDescription& d) {
        d.regions.push_back({"silicon", silicon()});
        d.node_region.back() = 1;
    });
    REQUIRE(mentions(repeated, "repeated"));
    REQUIRE(repeated.context->index == 1);
}

TEST_CASE("device: the temperature is checked against every region's material") {
    for (const double T : {0.0, -1.0, not_a_number, infinity, 3500.0, 1000.0}) {
        CAPTURE(T);
        const auto e = rejected([T](DeviceDescription& d) { d.temperature_K = T; });
        REQUIRE(mentions(e, "region 'silicon'"));
        REQUIRE(e.context->index == 0);
    }
    REQUIRE(mentions(rejected([](DeviceDescription& d) { d.temperature_K = 3500.0; }), "band gap"));
    REQUIRE(mentions(rejected([](DeviceDescription& d) { d.temperature_K = 1000.0; }), "mobility"));
}

TEST_CASE("device: per-node arrays are validated") {
    REQUIRE(mentions(rejected([](DeviceDescription& d) { d.node_region.pop_back(); }),
                     "one entry per mesh node"));
    rejected([](DeviceDescription& d) { d.donors.push_back(0.0); });
    rejected([](DeviceDescription& d) { d.acceptors.clear(); });
    for (const RegionId r : {-1, 1}) {
        const auto e = rejected([r](DeviceDescription& d) { d.node_region[7] = r; });
        REQUIRE(mentions(e, "out of range"));
        REQUIRE(e.context->index == 7);
    }
    const auto unused =
        rejected([](DeviceDescription& d) { d.regions.push_back({"spare", silicon()}); });
    REQUIRE(mentions(unused, "no nodes"));
    REQUIRE(unused.context->index == 1);
}

TEST_CASE("device: doping must be finite and non-negative") {
    for (const double c : {-1.0, not_a_number, infinity}) {
        CAPTURE(c);
        const auto donor = rejected([c](DeviceDescription& d) { d.donors[12] = c; });
        REQUIRE(donor.context->index == 12);
        const auto acceptor = rejected([c](DeviceDescription& d) { d.acceptors[4] = c; });
        REQUIRE(acceptor.context->index == 4);
    }
    const auto value = rejected([](DeviceDescription& d) { d.acceptors[4] = -2.5; });
    REQUIRE(value.context->value == -2.5);
}

TEST_CASE("device: contacts are validated") {
    const auto unnamed = rejected([](DeviceDescription& d) { d.contacts[1].name.clear(); });
    REQUIRE(unnamed.context->index == 1);
    REQUIRE(mentions(rejected([](DeviceDescription& d) { d.contacts[1].name = "anode"; }),
                     "repeated"));
    REQUIRE(mentions(rejected([](DeviceDescription& d) {
                         d.contacts[0].kind = static_cast<ContactKind>(7);
                     }),
                     "kind"));
    REQUIRE(mentions(rejected([](DeviceDescription& d) { d.contacts[0].nodes.clear(); }),
                     "no nodes"));
    for (const NodeId bad : {-1, 21}) {
        const auto e = rejected([bad](DeviceDescription& d) { d.contacts[1].nodes = {bad}; });
        REQUIRE(mentions(e, "out of range"));
        REQUIRE(e.context->index == 1);
    }
}

TEST_CASE("device: contact nodes are strictly increasing, on the boundary and not shared") {
    // 2D: x_min has nodes 0, 21, 42, ... (x varies fastest).
    REQUIRE(mentions(rejected([](DeviceDescription& d) { d.contacts[0].nodes = {21, 0}; }, 4),
                     "strictly increasing"));
    REQUIRE(mentions(rejected([](DeviceDescription& d) { d.contacts[0].nodes = {0, 0}; }, 4),
                     "strictly increasing"));
    const auto interior = rejected([](DeviceDescription& d) { d.contacts[0].nodes = {0, 22}; }, 4);
    REQUIRE(mentions(interior, "boundary"));
    REQUIRE(interior.context->value == 22.0);
    // Node 0 is in both x_min and y_min; giving it to a second contact is a conflict.
    const auto shared = rejected([](DeviceDescription& d) {
        d.contacts.push_back({"substrate", ContactKind::ohmic, {0, 1, 2}});
    }, 4);
    REQUIRE(mentions(shared, "another contact"));
    REQUIRE(shared.context->index == 2);
    REQUIRE(shared.context->value == 0.0);
}

namespace {

// A gate on part of the 2D diode's y_min face, between the contacts (nodes 5..15).
Contact gate_contact() {
    Contact g{"gate", ContactKind::gate, {}};
    for (NodeId v = 5; v <= 15; ++v) g.nodes.push_back(v);
    g.gate = {.boundary = "y_min", .oxide_thickness_cm = 5e-7};
    return g;
}

NiTCAD::base::Error gate_rejected(const std::function<void(Contact&)>& spoil) {
    return rejected(
        [&](DeviceDescription& d) {
            d.contacts.push_back(gate_contact());
            spoil(d.contacts.back());
        },
        4);
}

}  // namespace

TEST_CASE("device: a gate contact is accepted and read back") {
    DeviceDescription desc = diode(4);
    desc.contacts.push_back(gate_contact());
    desc.contacts.back().gate.electrode = GateElectrode::metal;
    desc.contacts.back().gate.work_function_eV = 4.1;
    desc.contacts.back().gate.fixed_charge_cm2 = -1e11;  // a negative fixed charge is allowed
    const auto d = Device::create(std::move(desc));
    REQUIRE(d.has_value());
    const Contact* g = d->find_contact("gate");
    REQUIRE(g != nullptr);
    REQUIRE(g->kind == ContactKind::gate);
    REQUIRE(g->nodes.size() == 11);
    REQUIRE(g->gate.boundary == "y_min");
    REQUIRE(g->gate.oxide_thickness_cm == 5e-7);
    REQUIRE(g->gate.oxide_relative_permittivity == 3.9);  // SiO2 default
    REQUIRE(g->gate.work_function_eV == 4.1);
    REQUIRE(g->gate.fixed_charge_cm2 == -1e11);
    // The default stack of an ohmic contact is not read.
    REQUIRE(d->find_contact("anode")->gate.oxide_thickness_cm == 0.0);
}

TEST_CASE("device: a gate stack is validated") {
    for (const double t : {0.0, -5e-7, not_a_number, infinity}) {
        CAPTURE(t);
        const auto e = gate_rejected([t](Contact& g) { g.gate.oxide_thickness_cm = t; });
        REQUIRE(mentions(e, "gate 'gate': oxide thickness"));
        REQUIRE(e.context->index == 2);
    }
    for (const double eps : {0.0, -3.9, not_a_number}) {
        CAPTURE(eps);
        REQUIRE(mentions(
            gate_rejected([eps](Contact& g) { g.gate.oxide_relative_permittivity = eps; }),
            "oxide permittivity"));
    }
    REQUIRE(mentions(gate_rejected([](Contact& g) {
                         g.gate.electrode = static_cast<GateElectrode>(9);
                     }),
                     "unknown gate electrode"));
    for (const double phi : {0.0, -4.1, not_a_number}) {
        CAPTURE(phi);
        REQUIRE(mentions(gate_rejected([phi](Contact& g) {
                             g.gate.electrode = GateElectrode::metal;
                             g.gate.work_function_eV = phi;
                         }),
                         "work function"));
    }
    // A polysilicon electrode takes its work function from the semiconductor; the field is unused.
    DeviceDescription poly = diode(4);
    poly.contacts.push_back(gate_contact());
    poly.contacts.back().gate.work_function_eV = not_a_number;
    REQUIRE(Device::create(std::move(poly)).has_value());
    for (const double q : {not_a_number, -infinity}) {
        CAPTURE(q);
        REQUIRE(mentions(gate_rejected([q](Contact& g) { g.gate.fixed_charge_cm2 = q; }),
                         "fixed oxide charge"));
    }
}

TEST_CASE("device: a gate lies on its boundary patch") {
    const auto missing = gate_rejected([](Contact& g) { g.gate.boundary = "top"; });
    REQUIRE(mentions(missing, "no mesh boundary patch 'top'"));
    REQUIRE(missing.context->index == 2);
    // Node 70 is on y_max (the last row starts at 63), not y_min.
    const auto off = gate_rejected([](Contact& g) { g.nodes.push_back(70); });
    REQUIRE(mentions(off, "not on boundary patch 'y_min'"));
    REQUIRE(off.context->value == 70.0);
    // A gate node is still a contact node: it may not be shared with an ohmic contact.
    REQUIRE(mentions(gate_rejected([](Contact& g) { g.nodes.insert(g.nodes.begin(), 0); }),
                     "another contact"));
}

// Unit 15: declared interfaces.

namespace {

// The diode with its n side (nodes 10 on) in a second region "other", plus a third region "far"
// on the last node only, so "silicon" and "far" share no edge.
DeviceDescription three_regions() {
    DeviceDescription d = diode();
    d.regions.push_back({"other", silicon()});
    d.regions.push_back({"far", silicon()});
    for (std::size_t i = 10; i < d.node_region.size(); ++i) d.node_region[i] = 1;
    d.node_region.back() = 2;
    return d;
}

}  // namespace

TEST_CASE("device: interfaces are validated and looked up in either order") {
    DeviceDescription ok = three_regions();
    ok.interfaces = {{"other", "silicon", InterfaceTransport::thermionic_emission}};
    const auto d = Device::create(std::move(ok));
    REQUIRE(d.has_value());
    REQUIRE(d->transport(0, 1) == InterfaceTransport::thermionic_emission);
    REQUIRE(d->transport(1, 0) == InterfaceTransport::thermionic_emission);
    REQUIRE(d->transport(1, 2) == InterfaceTransport::drift_diffusion);  // undeclared
    const auto rejected_with = [](std::vector<Interface> interfaces) {
        DeviceDescription desc = three_regions();
        desc.interfaces = std::move(interfaces);
        auto r = Device::create(std::move(desc));
        REQUIRE_FALSE(r.has_value());
        REQUIRE(r.error().code == ErrorCode::invalid_input);
        return r.error();
    };
    REQUIRE(rejected_with({{"silicon", "nowhere"}}).context->index == 0u);
    REQUIRE(mentions(rejected_with({{"silicon", "silicon"}}), "itself"));
    REQUIRE(mentions(rejected_with({{"silicon", "other"}, {"other", "silicon"}}), "twice"));
    REQUIRE(rejected_with({{"silicon", "other"}, {"other", "silicon"}}).context->index == 1u);
    REQUIRE(mentions(rejected_with({{"silicon", "far"}}), "no mesh edge"));
    REQUIRE(mentions(rejected_with({{"silicon", "other", static_cast<InterfaceTransport>(7)}}),
                     "unknown transport"));
}

TEST_CASE("device: thermal contacts (Unit 23) are validated and read back") {
    DeviceDescription ok = diode();
    ok.thermal_contacts = {
        {"sink", "x_min", ThermalContactKind::isothermal, 300.0, 0.0},
        {"package", "x_max", ThermalContactKind::resistance, 320.0, 2.5}};
    const auto d = Device::create(std::move(ok));
    REQUIRE(d.has_value());
    REQUIRE(d->thermal_contacts().size() == 2);
    REQUIRE(d->thermal_contacts()[1].name == "package");
    REQUIRE(d->thermal_contacts()[1].kind == ThermalContactKind::resistance);
    REQUIRE(d->thermal_contacts()[1].temperature_K == 320.0);
    REQUIRE(d->thermal_contacts()[1].resistance_K_cm2_W == 2.5);
    REQUIRE(Device::create(diode())->thermal_contacts().empty());

    const auto spoiled = [](std::vector<ThermalContact> contacts) {
        return rejected([&](DeviceDescription& desc) { desc.thermal_contacts = contacts; });
    };
    const ThermalContact sink{"sink", "x_min", ThermalContactKind::isothermal, 300.0, 0.0};
    REQUIRE(mentions(spoiled({{"", "x_min"}}), "empty or repeated"));
    const auto twice = spoiled({sink, {"sink", "x_max"}});
    REQUIRE(mentions(twice, "empty or repeated"));
    REQUIRE(twice.context->index == 1u);
    REQUIRE(mentions(spoiled({{"t", "x_min", static_cast<ThermalContactKind>(9)}}),
                     "unknown kind"));
    REQUIRE(mentions(spoiled({{"t", "nowhere"}}), "no boundary patch 'nowhere'"));
    const auto shared = spoiled({sink, {"other", "x_min"}});
    REQUIRE(mentions(shared, "shares a node"));
    REQUIRE(shared.context->index == 1u);
    for (const double T : {0.0, -5.0, not_a_number, infinity}) {
        REQUIRE(mentions(spoiled({{"t", "x_min", ThermalContactKind::isothermal, T}}),
                         "finite and positive"));
    }
    // A temperature silicon's models refuse (mu_max below mu_min for holes above about 857 K).
    const auto hot = spoiled({{"t", "x_min", ThermalContactKind::isothermal, 1000.0}});
    REQUIRE(mentions(hot, "region 'silicon'"));
    REQUIRE(hot.context->value == 1000.0);
    for (const double R : {0.0, -1.0, not_a_number, infinity}) {
        REQUIRE(mentions(spoiled({{"t", "x_min", ThermalContactKind::resistance, 300.0, R}}),
                         "thermal resistance"));
    }
    // An isothermal contact does not read its resistance.
    DeviceDescription iso = diode();
    iso.thermal_contacts = {{"t", "x_min", ThermalContactKind::isothermal, 300.0, not_a_number}};
    REQUIRE(Device::create(std::move(iso)).has_value());
}
