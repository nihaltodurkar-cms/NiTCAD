#include "NiTCAD/device/device.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::device {

namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt,
                    std::optional<double> value = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = value}};
}

// Disjoint sets over node ids, for the connected parts of the mesh graph.
class Components {
public:
    explicit Components(std::size_t n) : parent_(n) {
        std::iota(parent_.begin(), parent_.end(), mesh::NodeId{0});
    }
    mesh::NodeId root(mesh::NodeId v) {
        while (parent_[static_cast<std::size_t>(v)] != v) {
            auto& p = parent_[static_cast<std::size_t>(v)];
            p = parent_[static_cast<std::size_t>(p)];  // path halving
            v = p;
        }
        return v;
    }
    void join(mesh::NodeId a, mesh::NodeId b) {
        a = root(a);
        b = root(b);
        if (a != b) parent_[static_cast<std::size_t>(a < b ? b : a)] = a < b ? a : b;
    }

private:
    std::vector<mesh::NodeId> parent_;
};

std::optional<base::Error> check_regions(const DeviceDescription& d) {
    if (d.regions.empty()) return invalid("device needs at least one region");
    std::set<std::string_view> names;
    for (std::size_t r = 0; r < d.regions.size(); ++r) {
        const Region& region = d.regions[r];
        if (region.name.empty() || !names.insert(region.name).second) {
            return invalid("region name is empty or repeated", r);
        }
        if (auto ok = physics::check_temperature(region.material, d.temperature_K); !ok) {
            return invalid("region '" + region.name + "': " + ok.error().message, r,
                           ok.error().context ? ok.error().context->value : std::nullopt);
        }
    }
    return std::nullopt;
}

std::optional<base::Error> check_nodes(const DeviceDescription& d) {
    const std::size_t n = d.mesh.node_count();
    if (d.node_region.size() != n || d.donors.size() != n || d.acceptors.size() != n) {
        return invalid("node_region, donors and acceptors need one entry per mesh node");
    }
    std::vector<bool> used(d.regions.size(), false);
    for (std::size_t i = 0; i < n; ++i) {
        const RegionId r = d.node_region[i];
        if (r < 0 || static_cast<std::size_t>(r) >= d.regions.size()) {
            return invalid("node region out of range", i);
        }
        used[static_cast<std::size_t>(r)] = true;
        for (const double c : {d.donors[i], d.acceptors[i]}) {
            if (!(std::isfinite(c) && c >= 0.0)) {
                return invalid("donor or acceptor concentration is not finite and >= 0", i, c);
            }
        }
    }
    for (std::size_t r = 0; r < used.size(); ++r) {
        if (!used[r]) return invalid("region has no nodes", r);
    }
    return std::nullopt;
}

// The gate stack of a gate contact; the nodes are checked to be on its boundary patch.
std::optional<base::Error> check_gate(const DeviceDescription& d, std::size_t c) {
    const Contact& contact = d.contacts[c];
    const GateStack& g = contact.gate;
    const std::string where = "gate '" + contact.name + "': ";
    const auto positive = [](double v) { return std::isfinite(v) && v > 0.0; };
    if (!positive(g.oxide_thickness_cm)) {
        return invalid(where + "oxide thickness is not finite and positive", c,
                       g.oxide_thickness_cm);
    }
    if (!positive(g.oxide_relative_permittivity)) {
        return invalid(where + "oxide permittivity is not finite and positive", c,
                       g.oxide_relative_permittivity);
    }
    switch (g.electrode) {
        case GateElectrode::n_poly:
        case GateElectrode::p_poly:
            break;
        case GateElectrode::metal:
            if (!positive(g.work_function_eV)) {
                return invalid(where + "work function is not finite and positive", c,
                               g.work_function_eV);
            }
            break;
        default:
            return invalid(where + "unknown gate electrode", c);
    }
    if (!std::isfinite(g.fixed_charge_cm2)) {
        return invalid(where + "fixed oxide charge is not finite", c, g.fixed_charge_cm2);
    }
    const mesh::BoundaryPatch* patch = d.mesh.find_boundary(g.boundary);
    if (patch == nullptr) {
        return invalid(where + "no mesh boundary patch '" + g.boundary + "'", c);
    }
    // Both lists are strictly increasing.
    auto it = patch->nodes.begin();
    for (const mesh::NodeId v : contact.nodes) {
        while (it != patch->nodes.end() && *it < v) ++it;
        if (it == patch->nodes.end() || *it != v) {
            return invalid(where + "node is not on boundary patch '" + g.boundary + "'", c,
                           static_cast<double>(v));
        }
    }
    return std::nullopt;
}

std::optional<base::Error> check_contacts(const DeviceDescription& d) {
    const std::size_t n = d.mesh.node_count();
    std::vector<bool> on_boundary(n, false);
    for (const mesh::BoundaryPatch& patch : d.mesh.boundary()) {
        for (const mesh::NodeId v : patch.nodes) on_boundary[static_cast<std::size_t>(v)] = true;
    }
    std::vector<bool> in_contact(n, false);
    std::set<std::string_view> names;
    for (std::size_t c = 0; c < d.contacts.size(); ++c) {
        const Contact& contact = d.contacts[c];
        if (contact.name.empty() || !names.insert(contact.name).second) {
            return invalid("contact name is empty or repeated", c);
        }
        if (contact.kind != ContactKind::ohmic && contact.kind != ContactKind::gate) {
            return invalid("unknown contact kind", c);
        }
        if (contact.nodes.empty()) return invalid("contact '" + contact.name + "' has no nodes", c);
        for (std::size_t k = 0; k < contact.nodes.size(); ++k) {
            const mesh::NodeId v = contact.nodes[k];
            if (v < 0 || static_cast<std::size_t>(v) >= n ||
                (k > 0 && v <= contact.nodes[k - 1])) {
                return invalid("contact '" + contact.name +
                                   "': nodes out of range or not strictly increasing",
                               c);
            }
            const auto i = static_cast<std::size_t>(v);
            if (!on_boundary[i]) {
                return invalid("contact '" + contact.name + "': node is on no mesh boundary patch",
                               c, static_cast<double>(v));
            }
            if (in_contact[i]) {
                return invalid("contact '" + contact.name + "': node belongs to another contact",
                               c, static_cast<double>(v));
            }
            in_contact[i] = true;
        }
        if (contact.kind == ContactKind::gate) {
            if (auto e = check_gate(d, c)) return e;
        }
    }
    return std::nullopt;
}

// Every connected part of the mesh graph must contain an ohmic contact node (ARCHITECTURE.md 6.10,
// the topology gate for Unit 6). A gate fixes no carrier density, so on its own it leaves the
// continuity equations of its part singular.
std::optional<base::Error> check_topology(const DeviceDescription& d) {
    const std::size_t n = d.mesh.node_count();
    Components parts(n);
    for (const mesh::Edge& e : d.mesh.edges()) parts.join(e.first, e.second);
    std::vector<bool> anchored(n, false);
    for (const Contact& contact : d.contacts) {
        if (contact.kind != ContactKind::ohmic) continue;
        for (const mesh::NodeId v : contact.nodes) {
            anchored[static_cast<std::size_t>(parts.root(v))] = true;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        // Nodes are visited in increasing order, so the first one found is the lowest node of a
        // floating part.
        if (!anchored[static_cast<std::size_t>(parts.root(static_cast<mesh::NodeId>(i)))]) {
            return invalid("floating region: part of the mesh is connected to no ohmic contact",
                           i);
        }
    }
    return std::nullopt;
}

std::optional<RegionId> region_named(const DeviceDescription& d, const std::string& name) {
    for (std::size_t r = 0; r < d.regions.size(); ++r) {
        if (d.regions[r].name == name) return static_cast<RegionId>(r);
    }
    return std::nullopt;
}

std::optional<base::Error> check_interfaces(const DeviceDescription& d) {
    std::set<std::pair<RegionId, RegionId>> pairs;
    for (std::size_t k = 0; k < d.interfaces.size(); ++k) {
        const Interface& f = d.interfaces[k];
        const auto a = region_named(d, f.region_a), b = region_named(d, f.region_b);
        if (!a || !b) return invalid("interface names an unknown region", k);
        if (*a == *b) return invalid("interface joins region '" + f.region_a + "' to itself", k);
        if (f.transport != InterfaceTransport::drift_diffusion &&
            f.transport != InterfaceTransport::thermionic_emission) {
            return invalid("interface has an unknown transport", k);
        }
        if (!pairs.insert(std::minmax(*a, *b)).second) {
            return invalid("interface between '" + f.region_a + "' and '" + f.region_b +
                               "' is declared twice",
                           k);
        }
        bool joined = false;
        for (const mesh::Edge& e : d.mesh.edges()) {
            const RegionId ra = d.node_region[static_cast<std::size_t>(e.first)];
            const RegionId rb = d.node_region[static_cast<std::size_t>(e.second)];
            if (std::minmax(ra, rb) == std::minmax(*a, *b)) {
                joined = true;
                break;
            }
        }
        if (!joined) {
            return invalid("no mesh edge joins regions '" + f.region_a + "' and '" + f.region_b +
                               "'",
                           k);
        }
    }
    return std::nullopt;
}

}  // namespace

std::expected<Device, base::Error> Device::create(DeviceDescription description) {
    for (auto check :
         {check_regions, check_nodes, check_contacts, check_topology, check_interfaces}) {
        if (auto e = check(description)) return std::unexpected(std::move(*e));
    }
    return Device{std::move(description)};
}

InterfaceTransport Device::transport(RegionId a, RegionId b) const {
    NITCAD_EXPECTS(a >= 0 && static_cast<std::size_t>(a) < d_.regions.size());
    NITCAD_EXPECTS(b >= 0 && static_cast<std::size_t>(b) < d_.regions.size());
    for (const Interface& f : d_.interfaces) {
        const std::string& na = d_.regions[static_cast<std::size_t>(a)].name;
        const std::string& nb = d_.regions[static_cast<std::size_t>(b)].name;
        if ((f.region_a == na && f.region_b == nb) || (f.region_a == nb && f.region_b == na)) {
            return f.transport;
        }
    }
    return InterfaceTransport::drift_diffusion;
}

const Contact* Device::find_contact(std::string_view name) const noexcept {
    for (const Contact& c : d_.contacts) {
        if (c.name == name) return &c;
    }
    return nullptr;
}

double Device::net_doping(mesh::NodeId node) const {
    NITCAD_EXPECTS(node >= 0 && static_cast<std::size_t>(node) < d_.donors.size());
    const auto i = static_cast<std::size_t>(node);
    return d_.donors[i] - d_.acceptors[i];
}

double Device::total_impurity(mesh::NodeId node) const {
    NITCAD_EXPECTS(node >= 0 && static_cast<std::size_t>(node) < d_.donors.size());
    const auto i = static_cast<std::size_t>(node);
    return d_.donors[i] + d_.acceptors[i];
}

const physics::Semiconductor& Device::material(mesh::NodeId node) const {
    NITCAD_EXPECTS(node >= 0 && static_cast<std::size_t>(node) < d_.node_region.size());
    return d_.regions[static_cast<std::size_t>(d_.node_region[static_cast<std::size_t>(node)])]
        .material;
}

}  // namespace NiTCAD::device
