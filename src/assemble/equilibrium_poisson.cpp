#include "NiTCAD/assemble/equilibrium_poisson.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/ohmic.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "scaled_device.hpp"

namespace NiTCAD::assemble {

std::expected<EquilibriumPoisson, base::Error> EquilibriumPoisson::create(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models) {
    auto scaled = detail::make_scaled_device(device, scaling, models);
    if (!scaled) return std::unexpected(std::move(scaled.error()));
    const std::size_t n = scaled->volume.size();

    EquilibriumPoisson p;
    p.volume_ = std::move(scaled->volume);
    p.doping_ = std::move(scaled->doping);
    p.n_ie_ = std::move(scaled->n_ie);
    p.psi0_.assign(n, 0.0);
    p.contact_.assign(n, 0);
    p.gates_ = std::move(scaled->gates);
    p.kinds_ = std::move(scaled->kinds);
    for (std::size_t i = 0; i < n; ++i) {
        if (scaled->contact[i] < 0) continue;
        p.contact_[i] = 1;
        p.psi0_[i] = ohmic_contact_value(p.doping_[i], p.n_ie_[i], 0.0).psi;
    }

    std::vector<linalg::Triplet> triplets;
    triplets.reserve(n + 2 * scaled->edges.size());
    for (std::size_t i = 0; i < n; ++i) {
        const auto k = static_cast<linalg::Index>(i);
        triplets.push_back({k, k, 0.0});
    }
    for (const detail::ScaledEdge& e : scaled->edges) {
        const auto i = static_cast<linalg::Index>(e.i);
        const auto j = static_cast<linalg::Index>(e.j);
        triplets.push_back({i, j, 0.0});
        triplets.push_back({j, i, 0.0});
    }
    const auto size = static_cast<linalg::Index>(n);
    auto pattern = linalg::SparseMatrix::from_triplets(size, size, triplets);
    if (!pattern) return std::unexpected(std::move(pattern.error()));
    p.pattern_ = std::move(*pattern);

    p.diag_.resize(n);
    for (std::size_t i = 0; i < n; ++i) p.diag_[i] = detail::position(p.pattern_, i, i);
    p.edges_.reserve(scaled->edges.size());
    for (const detail::ScaledEdge& e : scaled->edges) {
        p.edges_.push_back({e.i, e.j, e.et * e.geometry, detail::position(p.pattern_, e.i, e.j),
                            detail::position(p.pattern_, e.j, e.i)});
    }
    return p;
}

std::expected<void, base::Error> EquilibriumPoisson::set_bias(std::span<const double> bias_V) {
    if (auto ok = check_contact_bias(kinds_, bias_V, true); !ok) return ok;
    gates_.set_bias(bias_V);
    return {};
}

std::vector<double> EquilibriumPoisson::gate_charges(std::span<const double> psi) const {
    NITCAD_EXPECTS(psi.size() == unknowns());
    return gates_.charges(psi, 1, kinds_.size());
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
        if (gates_.on_gate(i)) {
            residual[i] += gates_.row_term(i, psi[i]);
            if (!jacobian_values.empty()) jacobian_values[diag_[i]] -= gates_.term(i).coupling;
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

void EquilibriumPoisson::carriers(std::span<const double> psi, std::span<double> n,
                                  std::span<double> p) const {
    NITCAD_EXPECTS(psi.size() == unknowns() && n.size() == unknowns() && p.size() == unknowns());
    for (std::size_t i = 0; i < unknowns(); ++i) {
        n[i] = physics::boltzmann_density(n_ie_[i], psi[i]).density;
        p[i] = physics::boltzmann_density(n_ie_[i], -psi[i]).density;
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
