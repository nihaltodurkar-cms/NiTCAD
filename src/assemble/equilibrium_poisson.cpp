#include "NiTCAD/assemble/equilibrium_poisson.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/assemble/ohmic.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::assemble {

namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

// Position of (row, col) in the CSR values; the entry must exist.
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

}  // namespace

std::expected<EquilibriumPoisson, base::Error> EquilibriumPoisson::create(
    const device::Device& device, const Scaling& scaling) {
    if (scaling.temperature_K != device.temperature_K()) {
        return std::unexpected(invalid("scaling temperature differs from the device's"));
    }
    const mesh::Mesh& m = device.mesh();
    const auto regions = device.regions();
    const auto node_region = device.node_region();
    for (std::size_t e = 0; e < m.edges().size(); ++e) {
        const mesh::Edge& edge = m.edges()[e];
        const auto region_of = [&](mesh::NodeId v) {
            return static_cast<std::size_t>(node_region[static_cast<std::size_t>(v)]);
        };
        const std::size_t ra = region_of(edge.first);
        const std::size_t rb = region_of(edge.second);
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
    const double eps_r_ref = scaling.eps_F_per_cm / base::eps0_F_per_cm;

    EquilibriumPoisson p;
    p.volume_.resize(n);
    p.doping_.resize(n);
    p.n_ie_.resize(n);
    p.psi0_.assign(n, 0.0);
    p.contact_.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        p.volume_[i] = m.volumes()[i] / volume_scale;
        p.doping_[i] = device.net_doping(node) / scaling.Ns;
        p.n_ie_[i] = physics::intrinsic_density(device.material(node), scaling.temperature_K) /
                     scaling.Ns;
    }
    for (const device::Contact& contact : device.contacts()) {
        for (const mesh::NodeId v : contact.nodes) {
            const auto i = static_cast<std::size_t>(v);
            p.contact_[i] = 1;
            p.psi0_[i] = ohmic_contact_value(p.doping_[i], p.n_ie_[i], 0.0).psi;
        }
    }

    std::vector<linalg::Triplet> triplets;
    triplets.reserve(n + 2 * m.edges().size());
    for (std::size_t i = 0; i < n; ++i) {
        const auto k = static_cast<linalg::Index>(i);
        triplets.push_back({k, k, 0.0});
    }
    for (const mesh::Edge& edge : m.edges()) {
        triplets.push_back({edge.first, edge.second, 0.0});
        triplets.push_back({edge.second, edge.first, 0.0});
    }
    const auto size = static_cast<linalg::Index>(n);
    auto pattern = linalg::SparseMatrix::from_triplets(size, size, triplets);
    if (!pattern) return std::unexpected(std::move(pattern.error()));
    p.pattern_ = std::move(*pattern);

    p.diag_.resize(n);
    for (std::size_t i = 0; i < n; ++i) p.diag_[i] = position(p.pattern_, i, i);
    p.edges_.reserve(m.edges().size());
    for (const mesh::Edge& edge : m.edges()) {
        const auto i = static_cast<std::size_t>(edge.first);
        const auto j = static_cast<std::size_t>(edge.second);
        // Both ends are in the same material (heterojunctions are rejected above); its permittivity
        // relative to the reference one is the legacy et.
        const double et = device.material(edge.first).parameters().eps_r / eps_r_ref;
        p.edges_.push_back({i, j, et * edge.coupling_area / edge.length / coupling_scale,
                            position(p.pattern_, i, j), position(p.pattern_, j, i)});
    }
    return p;
}

linalg::SparseMatrix EquilibriumPoisson::make_jacobian() const { return pattern_; }

void EquilibriumPoisson::residual_into(std::span<const double> psi, std::span<double> residual,
                                       std::span<double> jacobian_values) const {
    const std::size_t n = unknowns();
    NITCAD_EXPECTS(psi.size() == n && residual.size() == n);
    for (std::size_t i = 0; i < n; ++i) {
        if (contact_[i] != 0) {
            residual[i] = psi[i] - psi0_[i];
            continue;
        }
        const double carriers_n = physics::boltzmann_density(n_ie_[i], psi[i]).density;
        const double carriers_p = physics::boltzmann_density(n_ie_[i], -psi[i]).density;
        residual[i] = -volume_[i] * (carriers_n - carriers_p - doping_[i]);
        // d/dpsi of the charge term, onto the diagonal (n' = n, p' = -p).
        if (!jacobian_values.empty()) {
            jacobian_values[diag_[i]] = -volume_[i] * (carriers_n + carriers_p);
        }
    }
    for (const EdgeTerm& e : edges_) {
        const double flux = e.c * (psi[e.j] - psi[e.i]);
        if (contact_[e.i] == 0) residual[e.i] += flux;
        if (contact_[e.j] == 0) residual[e.j] -= flux;
    }
}

void EquilibriumPoisson::residual(std::span<const double> psi, std::span<double> residual) const {
    residual_into(psi, residual, {});
}

void EquilibriumPoisson::evaluate(std::span<const double> psi, std::span<double> residual,
                                  linalg::SparseMatrix& jacobian) const {
    NITCAD_EXPECTS(jacobian.has_same_pattern(pattern_));
    const std::span<double> values = jacobian.values();
    std::fill(values.begin(), values.end(), 0.0);
    residual_into(psi, residual, values);
    for (std::size_t i = 0; i < unknowns(); ++i) {
        if (contact_[i] != 0) values[diag_[i]] = 1.0;
    }
    for (const EdgeTerm& e : edges_) {
        if (contact_[e.i] == 0) {
            values[e.ij] += e.c;
            values[diag_[e.i]] -= e.c;
        }
        if (contact_[e.j] == 0) {
            values[e.ji] += e.c;
            values[diag_[e.j]] -= e.c;
        }
    }
}

std::vector<double> EquilibriumPoisson::charge_neutral_potential() const {
    std::vector<double> psi(unknowns());
    for (std::size_t i = 0; i < psi.size(); ++i) {
        psi[i] = contact_[i] != 0
                     ? psi0_[i]
                     : physics::boltzmann_neutral_equilibrium(doping_[i], n_ie_[i]).eta;
    }
    return psi;
}

}  // namespace NiTCAD::assemble
