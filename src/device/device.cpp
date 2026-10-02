#include "NiTCAD/device/device.hpp"

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
        if (contact.kind != ContactKind::ohmic) return invalid("unknown contact kind", c);
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
    }
    return std::nullopt;
}

// Every connected part of the mesh graph must contain a contact node (ARCHITECTURE.md 6.10, the
// topology gate for Unit 6).
std::optional<base::Error> check_topology(const DeviceDescription& d) {
    const std::size_t n = d.mesh.node_count();
    Components parts(n);
    for (const mesh::Edge& e : d.mesh.edges()) parts.join(e.first, e.second);
    std::vector<bool> anchored(n, false);
    for (const Contact& contact : d.contacts) {
        for (const mesh::NodeId v : contact.nodes) {
            anchored[static_cast<std::size_t>(parts.root(v))] = true;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        // Nodes are visited in increasing order, so the first one found is the lowest node of a
        // floating part.
        if (!anchored[static_cast<std::size_t>(parts.root(static_cast<mesh::NodeId>(i)))]) {
            return invalid("floating region: part of the mesh is connected to no contact", i);
        }
    }
    return std::nullopt;
}

}  // namespace

std::expected<Device, base::Error> Device::create(DeviceDescription description) {
    for (auto check : {check_regions, check_nodes, check_contacts, check_topology}) {
        if (auto e = check(description)) return std::unexpected(std::move(*e));
    }
    return Device{std::move(description)};
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
