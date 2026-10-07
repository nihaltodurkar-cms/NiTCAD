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

    const std::size_t n = m.node_count();
    const int D = m.dimension();
    const double volume_scale = std::pow(scaling.L_D, D);
    const double coupling_scale = std::pow(scaling.L_D, D - 2);

    ScaledDevice s;
    s.volume.resize(n);
    s.doping.resize(n);
    s.donors.resize(n);
    s.acceptors.resize(n);
    s.levels.resize(n);
    s.radiative.resize(n);
    s.n_ie.resize(n);
    s.log_dos_n.resize(n);
    s.log_dos_p.resize(n);
    s.band_shift.resize(n);
    s.contact.assign(n, -1);
    const double T = scaling.temperature_K;
    const double reference_depth = physics::intrinsic_level_depth_eV(device.material(0), T);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        s.volume[i] = m.volumes()[i] / volume_scale;
        s.doping[i] = device.net_doping(node) / scaling.Ns;
        s.donors[i] = device.donors()[i] / scaling.Ns;
        s.acceptors[i] = device.acceptors()[i] / scaling.Ns;
        const physics::Semiconductor& material = device.material(node);
        const double n_ie =
            models.bgn ? physics::effective_intrinsic_density(
                             material, device.total_impurity(node), scaling.temperature_K)
                       : physics::intrinsic_density(material, scaling.temperature_K);
        s.n_ie[i] = n_ie / scaling.Ns;
        s.log_dos_n[i] = std::log(physics::conduction_band_dos(material, T) / n_ie);
        s.log_dos_p[i] = std::log(physics::valence_band_dos(material, T) / n_ie);
        const physics::IonizationParameters& ion = material.parameters().ionization;
        s.levels[i] = {ion.donor_eV / scaling.V_T, ion.acceptor_eV / scaling.V_T,
                       ion.donor_degeneracy, ion.acceptor_degeneracy};
        s.radiative[i] = material.parameters().radiative_cm3_s;
        const double depth = physics::intrinsic_level_depth_eV(material, T);
        s.band_shift[i] = (depth - reference_depth) / scaling.V_T;  // exactly 0 for one material
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
            // The electrode potential in the equations' reference: psi_G - s (exact for s = 0).
            gate_terms[i].offset += s.band_shift[i];
        }
    }
    s.gates = GateNodes(std::move(gate), std::move(gate_terms), scaling.V_T);
    s.edges.reserve(m.edges().size());
    const auto node_region = device.node_region();
    for (const mesh::Edge& edge : m.edges()) {
        const physics::SemiconductorParameters& a = device.material(edge.first).parameters();
        const physics::SemiconductorParameters& b = device.material(edge.second).parameters();
        const device::RegionId ra = node_region[static_cast<std::size_t>(edge.first)];
        const device::RegionId rb = node_region[static_cast<std::size_t>(edge.second)];
        const bool thermionic =
            ra != rb &&
            device.transport(ra, rb) == device::InterfaceTransport::thermionic_emission;
        const double eta = permittivity_ratio(a.eps_r, scaling);
        const double etb = permittivity_ratio(b.eps_r, scaling);
        s.edges.push_back({static_cast<std::size_t>(edge.first),
                           static_cast<std::size_t>(edge.second),
                           edge.coupling_area / edge.length / coupling_scale, edge.length,
                           eta == etb ? eta : 2.0 * eta * etb / (eta + etb), !(a == b),
                           thermionic});
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
