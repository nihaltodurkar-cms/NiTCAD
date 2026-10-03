#include "scaled_device.hpp"

#include <cmath>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::assemble::detail {

namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

}  // namespace

std::expected<ScaledDevice, base::Error> make_scaled_device(const device::Device& device,
                                                            const Scaling& scaling,
                                                            const PhysicsModels& models) {
    if (scaling.temperature_K != device.temperature_K()) {
        return std::unexpected(invalid("scaling temperature differs from the device's"));
    }
    const mesh::Mesh& m = device.mesh();
    const auto regions = device.regions();
    const auto node_region = device.node_region();
    const auto region_of = [&](mesh::NodeId v) {
        return static_cast<std::size_t>(node_region[static_cast<std::size_t>(v)]);
    };
    for (std::size_t e = 0; e < m.edges().size(); ++e) {
        const std::size_t ra = region_of(m.edges()[e].first);
        const std::size_t rb = region_of(m.edges()[e].second);
        if (ra != rb &&
            !(regions[ra].material.parameters() == regions[rb].material.parameters())) {
            return std::unexpected(invalid(
                "heterojunction between regions '" + regions[ra].name + "' and '" +
                    regions[rb].name + "': band offsets and permittivity steps are deferred",
                e));
        }
    }

    const std::size_t n = m.node_count();
    const int D = m.dimension();
    const double volume_scale = std::pow(scaling.L_D, D);
    const double coupling_scale = std::pow(scaling.L_D, D - 2);

    ScaledDevice s;
    s.volume.resize(n);
    s.doping.resize(n);
    s.n_ie.resize(n);
    s.contact.assign(n, -1);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        s.volume[i] = m.volumes()[i] / volume_scale;
        s.doping[i] = device.net_doping(node) / scaling.Ns;
        const physics::Semiconductor& material = device.material(node);
        const double n_ie =
            models.bgn ? physics::effective_intrinsic_density(
                             material, device.total_impurity(node), scaling.temperature_K)
                       : physics::intrinsic_density(material, scaling.temperature_K);
        s.n_ie[i] = n_ie / scaling.Ns;
    }
    std::vector<std::int32_t> gate(n, -1);
    std::vector<GateTerm> gate_terms(n, GateTerm{0.0, 0.0, 0.0});
    const auto contacts = device.contacts();
    for (std::size_t c = 0; c < contacts.size(); ++c) {
        const device::Contact& contact = contacts[c];
        s.kinds.push_back(contact.kind);
        if (contact.kind == device::ContactKind::ohmic) {
            for (const mesh::NodeId v : contact.nodes) {
                s.contact[static_cast<std::size_t>(v)] = static_cast<std::int32_t>(c);
            }
            continue;
        }
        // A gate: the device has checked that every node is on the patch (both lists increase).
        const mesh::BoundaryPatch* patch = m.find_boundary(contact.gate.boundary);
        NITCAD_EXPECTS(patch != nullptr);
        std::size_t k = 0;
        for (const mesh::NodeId v : contact.nodes) {
            while (patch->nodes[k] < v) ++k;
            const auto i = static_cast<std::size_t>(v);
            gate[i] = static_cast<std::int32_t>(c);
            gate_terms[i] =
                gate_term(contact.gate, device.material(v), patch->areas[k], D, scaling);
        }
    }
    s.gates = GateNodes(std::move(gate), std::move(gate_terms), scaling.V_T);
    s.edges.reserve(m.edges().size());
    for (const mesh::Edge& edge : m.edges()) {
        // Both ends are in the same material (heterojunctions are rejected above).
        const double eps_r = device.material(edge.first).parameters().eps_r;
        s.edges.push_back({static_cast<std::size_t>(edge.first),
                           static_cast<std::size_t>(edge.second),
                           edge.coupling_area / edge.length / coupling_scale, edge.length,
                           permittivity_ratio(eps_r, scaling)});
    }
    return s;
}

std::size_t position(const linalg::SparseMatrix& m, std::size_t row, std::size_t col) {
    const auto offsets = m.row_offsets();
    const auto cols = m.col_indices();
    for (auto k = static_cast<std::size_t>(offsets[row]);
         k < static_cast<std::size_t>(offsets[row + 1]); ++k) {
        if (static_cast<std::size_t>(cols[k]) == col) return k;
    }
    NITCAD_EXPECTS(false);
    return 0;
}

}  // namespace NiTCAD::assemble::detail
