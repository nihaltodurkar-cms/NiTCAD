#include "NiTCAD/assemble/drift_diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/ohmic.hpp"
#include "NiTCAD/assemble/sg_flux.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "scaled_device.hpp"

namespace NiTCAD::assemble {

namespace {

// The (row, column) offsets inside a node's three unknowns that an edge couples, in EdgeTerm order:
// Poisson to psi, electron continuity to psi and n, hole continuity to psi and p.
constexpr std::size_t edge_rows[5] = {0, 1, 1, 2, 2};
constexpr std::size_t edge_cols[5] = {0, 0, 1, 0, 2};

double harmonic_mean(double a, double b) { return 2.0 * a * b / (a + b); }

}  // namespace

std::expected<DriftDiffusion, base::Error> DriftDiffusion::create(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models) {
    auto scaled = detail::make_scaled_device(device, scaling, models);
    if (!scaled) return std::unexpected(std::move(scaled.error()));
    const std::size_t n = scaled->volume.size();
    const double T = scaling.temperature_K;

    DriftDiffusion s;
    s.volume_ = std::move(scaled->volume);
    s.doping_ = std::move(scaled->doping);
    s.n_ie_ = std::move(scaled->n_ie);
    s.contact_ = std::move(scaled->contact);
    s.gates_ = std::move(scaled->gates);
    s.kinds_ = std::move(scaled->kinds);
    s.contact_count_ = device.contacts().size();
    s.rate_scale_ = scaling.Ns / scaling.R0;
    s.Ns_ = scaling.Ns;
    s.auger_ = models.auger;
    s.V_T_ = scaling.V_T;
    s.srh_ = models.srh;
    s.field_mobility_ = models.field_mobility;
    s.tau_n_.resize(n);
    s.tau_p_.resize(n);
    s.auger_n_.resize(n);
    s.auger_p_.resize(n);
    std::vector<double> mu_n(n), mu_p(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        const physics::Semiconductor& m = device.material(node);
        const double N = device.total_impurity(node);
        const double N_mobility = models.doping_mobility ? N : 0.0;
        mu_n[i] = physics::caughey_thomas_mobility(m, physics::Carrier::electron, N_mobility, T);
        mu_p[i] = physics::caughey_thomas_mobility(m, physics::Carrier::hole, N_mobility, T);
        s.tau_n_[i] = physics::scharfetter_lifetime(m, physics::Carrier::electron, N);
        s.tau_p_[i] = physics::scharfetter_lifetime(m, physics::Carrier::hole, N);
        s.auger_n_[i] = m.parameters().auger.Cn;
        s.auger_p_[i] = m.parameters().auger.Cp;
    }

    std::vector<linalg::Triplet> triplets;
    triplets.reserve(9 * n + 10 * scaled->edges.size());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                triplets.push_back({static_cast<linalg::Index>(3 * i + r),
                                    static_cast<linalg::Index>(3 * i + c), 0.0});
            }
        }
    }
    for (const detail::ScaledEdge& e : scaled->edges) {
        for (std::size_t k = 0; k < 5; ++k) {
            triplets.push_back({static_cast<linalg::Index>(3 * e.i + edge_rows[k]),
                                static_cast<linalg::Index>(3 * e.j + edge_cols[k]), 0.0});
            triplets.push_back({static_cast<linalg::Index>(3 * e.j + edge_rows[k]),
                                static_cast<linalg::Index>(3 * e.i + edge_cols[k]), 0.0});
        }
    }
    const auto size = static_cast<linalg::Index>(3 * n);
    auto pattern = linalg::SparseMatrix::from_triplets(size, size, triplets);
    if (!pattern) return std::unexpected(std::move(pattern.error()));
    s.pattern_ = std::move(*pattern);

    s.block_.resize(9 * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                s.block_[9 * i + 3 * r + c] = detail::position(s.pattern_, 3 * i + r, 3 * i + c);
            }
        }
    }
    s.edges_.reserve(scaled->edges.size());
    for (const detail::ScaledEdge& e : scaled->edges) {
        EdgeTerm t{};
        t.a = e.i;
        t.b = e.j;
        t.c = e.et * e.geometry;
        t.mu_n = harmonic_mean(mu_n[e.i], mu_n[e.j]);
        t.mu_p = harmonic_mean(mu_p[e.i], mu_p[e.j]);
        t.an = t.mu_n * scaling.V_T / scaling.D0 * e.geometry;
        t.ap = t.mu_p * scaling.V_T / scaling.D0 * e.geometry;
        t.field = scaling.V_T / e.length_cm;
        const physics::Semiconductor& m = device.material(static_cast<mesh::NodeId>(e.i));
        t.sat_n = physics::saturation(m, physics::Carrier::electron);
        t.sat_p = physics::saturation(m, physics::Carrier::hole);
        t.dln = std::log(s.n_ie_[e.j] / s.n_ie_[e.i]);
        for (std::size_t k = 0; k < 5; ++k) {
            t.ab[k] = detail::position(s.pattern_, 3 * e.i + edge_rows[k], 3 * e.j + edge_cols[k]);
            t.ba[k] = detail::position(s.pattern_, 3 * e.j + edge_rows[k], 3 * e.i + edge_cols[k]);
        }
        s.edges_.push_back(t);
    }

    s.psi0_.assign(n, 0.0);
    s.n0_.assign(n, 0.0);
    s.p0_.assign(n, 0.0);
    const std::vector<double> zero(s.contact_count_, 0.0);
    if (auto ok = s.set_bias(zero); !ok) return std::unexpected(std::move(ok.error()));
    return s;
}

std::expected<void, base::Error> DriftDiffusion::set_bias(std::span<const double> bias_V) {
    if (auto ok = check_contact_bias(kinds_, bias_V, false); !ok) return ok;
    gates_.set_bias(bias_V);
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (contact_[i] < 0) continue;
        const double bias = bias_V[static_cast<std::size_t>(contact_[i])] / V_T_;
        const OhmicValue v = ohmic_contact_value(doping_[i], n_ie_[i], bias);
        psi0_[i] = v.psi;
        n0_[i] = v.n;
        p0_[i] = v.p;
    }
    return {};
}

std::vector<double> DriftDiffusion::state_from_potential(std::span<const double> psi) const {
    NITCAD_EXPECTS(psi.size() == node_count());
    std::vector<double> x(unknowns());
    for (std::size_t i = 0; i < node_count(); ++i) {
        x[3 * i] = psi[i];
        x[3 * i + 1] = physics::boltzmann_density(n_ie_[i], psi[i]).density;
        x[3 * i + 2] = physics::boltzmann_density(n_ie_[i], -psi[i]).density;
    }
    stamp_contacts(x);
    return x;
}

void DriftDiffusion::stamp_contacts(std::span<double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (contact_[i] < 0) continue;
        x[3 * i] = psi0_[i];
        x[3 * i + 1] = n0_[i];
        x[3 * i + 2] = p0_[i];
    }
}

linalg::SparseMatrix DriftDiffusion::make_jacobian() const { return pattern_; }

std::pair<EdgeFlux, EdgeFlux> DriftDiffusion::edge_fluxes(const EdgeTerm& e,
                                                          std::span<const double> x) const {
    // The flux functions take the driving term as psi2 - psi1; its partials are those with
    // respect to psi_a and psi_b, since ln(n_ie) does not depend on the state.
    const double dpsi = x[3 * e.b] - x[3 * e.a];
    const double na = x[3 * e.a + 1], nb = x[3 * e.b + 1];
    const double pa = x[3 * e.a + 2], pb = x[3 * e.b + 2];
    if (!field_mobility_) {
        return {sg_electron_flux(e.an, 0.0, dpsi + e.dln, na, nb),
                sg_hole_flux(e.ap, 0.0, dpsi - e.dln, pa, pb)};
    }
    // Edge factor a(E) = a0 mu_C(mu0, E) / mu0, E = field |dpsi|; the flux is linear in a, so
    // d flux / d dpsi gains (flux / a) da/d dpsi.
    const double E = e.field * std::abs(dpsi);
    const double dE = dpsi > 0.0 ? e.field : dpsi < 0.0 ? -e.field : 0.0;  // dE / d dpsi
    const auto scaled = [&](EdgeFlux unit, double a0, double mu0,
                            const physics::CanaliParameters& sat) {
        const physics::FieldMobility mu = physics::canali_mobility(mu0, E, sat);
        const double a = a0 * (mu.mobility / mu0);
        const double da = a0 * (mu.d_dE / mu0) * dE;
        return EdgeFlux{a * unit.flux, a * unit.d_psi1 - unit.flux * da,
                        a * unit.d_psi2 + unit.flux * da, a * unit.d_c1, a * unit.d_c2};
    };
    return {scaled(sg_electron_flux(1.0, 0.0, dpsi + e.dln, na, nb), e.an, e.mu_n, e.sat_n),
            scaled(sg_hole_flux(1.0, 0.0, dpsi - e.dln, pa, pb), e.ap, e.mu_p, e.sat_p)};
}

void DriftDiffusion::assemble(std::span<const double> x, std::span<double> f,
                              std::span<double> values) const {
    NITCAD_EXPECTS(x.size() == unknowns() && f.size() == unknowns());
    const bool jacobian = !values.empty();
    const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
        return values[block_[9 * node + 3 * r + c]];
    };
    for (std::size_t i = 0; i < node_count(); ++i) {
        const double psi = x[3 * i], n = x[3 * i + 1], p = x[3 * i + 2];
        if (contact_[i] >= 0) {
            f[3 * i] = psi - psi0_[i];
            f[3 * i + 1] = n - n0_[i];
            f[3 * i + 2] = p - p0_[i];
            if (jacobian) at(i, 0, 0) = at(i, 1, 1) = at(i, 2, 2) = 1.0;
            continue;
        }
        const double V = volume_[i];
        f[3 * i] = -V * (n - p - doping_[i]);
        f[3 * i + 1] = 0.0;
        f[3 * i + 2] = 0.0;
        if (jacobian) {
            at(i, 0, 1) = -V;
            at(i, 0, 2) = V;
        }
        if (gates_.on_gate(i)) {
            f[3 * i] += gates_.row_term(i, psi);
            if (jacobian) at(i, 0, 0) -= gates_.term(i).coupling;
        }
        if (srh_) {
            const physics::RecombinationRate r = physics::srh_recombination(
                n, p, physics::boltzmann_equilibrium_product(n_ie_[i]), n_ie_[i], tau_n_[i],
                tau_p_[i]);
            const double k = V * rate_scale_;
            f[3 * i + 1] -= k * r.rate;
            f[3 * i + 2] += k * r.rate;
            if (jacobian) {
                at(i, 1, 1) -= k * r.d_dn;
                at(i, 1, 2) -= k * r.d_dp;
                at(i, 2, 1) += k * r.d_dn;
                at(i, 2, 2) += k * r.d_dp;
            }
        }
        if (auger_) {
            // Physical densities; R / R0 = rate_scale (R / Ns), d(R / R0)/dn' = rate_scale dR/dn.
            const double nie = n_ie_[i] * Ns_;
            const physics::RecombinationRate r = physics::auger_recombination(
                n * Ns_, p * Ns_, physics::boltzmann_equilibrium_product(nie), auger_n_[i],
                auger_p_[i]);
            const double k = V * rate_scale_;
            f[3 * i + 1] -= k * (r.rate / Ns_);
            f[3 * i + 2] += k * (r.rate / Ns_);
            if (jacobian) {
                at(i, 1, 1) -= k * r.d_dn;
                at(i, 1, 2) -= k * r.d_dp;
                at(i, 2, 1) += k * r.d_dn;
                at(i, 2, 2) += k * r.d_dp;
            }
        }
    }
    for (const EdgeTerm& e : edges_) {
        const std::size_t a = e.a, b = e.b;
        const double psi_a = x[3 * a], psi_b = x[3 * b];
        const double poisson = e.c * (psi_b - psi_a);
        const auto [fn, fp] = edge_fluxes(e, x);
        if (contact_[a] < 0) {
            f[3 * a] += poisson;
            f[3 * a + 1] += fn.flux;
            f[3 * a + 2] += fp.flux;
            if (jacobian) {
                at(a, 0, 0) -= e.c;
                values[e.ab[0]] += e.c;
                at(a, 1, 0) += fn.d_psi1;
                values[e.ab[1]] += fn.d_psi2;
                at(a, 1, 1) += fn.d_c1;
                values[e.ab[2]] += fn.d_c2;
                at(a, 2, 0) += fp.d_psi1;
                values[e.ab[3]] += fp.d_psi2;
                at(a, 2, 2) += fp.d_c1;
                values[e.ab[4]] += fp.d_c2;
            }
        }
        if (contact_[b] < 0) {
            f[3 * b] -= poisson;
            f[3 * b + 1] -= fn.flux;
            f[3 * b + 2] -= fp.flux;
            if (jacobian) {
                at(b, 0, 0) -= e.c;
                values[e.ba[0]] += e.c;
                values[e.ba[1]] -= fn.d_psi1;
                at(b, 1, 0) -= fn.d_psi2;
                values[e.ba[2]] -= fn.d_c1;
                at(b, 1, 1) -= fn.d_c2;
                values[e.ba[3]] -= fp.d_psi1;
                at(b, 2, 0) -= fp.d_psi2;
                values[e.ba[4]] -= fp.d_c1;
                at(b, 2, 2) -= fp.d_c2;
            }
        }
    }
}

void DriftDiffusion::evaluate(std::span<const double> x, std::span<double> residual,
                              linalg::SparseMatrix& jacobian) const {
    NITCAD_EXPECTS(jacobian.has_same_pattern(pattern_));
    const std::span<double> values = jacobian.values();
    std::fill(values.begin(), values.end(), 0.0);
    assemble(x, residual, values);
}

void DriftDiffusion::residual(std::span<const double> x, std::span<double> residual) const {
    assemble(x, residual, {});
}

double DriftDiffusion::update_size(std::span<const double> x, std::span<const double> dx) const {
    NITCAD_EXPECTS(x.size() == unknowns() && dx.size() == unknowns());
    constexpr double infinity = std::numeric_limits<double>::infinity();
    double size = 0.0;
    for (std::size_t i = 0; i < node_count(); ++i) {
        const double n = x[3 * i + 1], p = x[3 * i + 2];
        size = std::max({size, std::abs(dx[3 * i]),
                         n > 0.0 ? std::abs(dx[3 * i + 1]) / n : infinity,
                         p > 0.0 ? std::abs(dx[3 * i + 2]) / p : infinity});
    }
    return size;
}

void DriftDiffusion::apply_update(std::span<double> x, std::span<const double> dx,
                                  double max_update) const {
    NITCAD_EXPECTS(x.size() == unknowns() && dx.size() == unknowns());
    for (std::size_t i = 0; i < node_count(); ++i) {
        x[3 * i] += std::clamp(dx[3 * i], -max_update, max_update);
        for (const std::size_t k : {3 * i + 1, 3 * i + 2}) {
            x[k] = std::clamp(x[k] + dx[k], 0.1 * x[k], 10.0 * x[k]);
        }
    }
}

std::vector<std::pair<double, double>> DriftDiffusion::edge_currents(
    std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<std::pair<double, double>> currents;
    currents.reserve(edges_.size());
    for (const EdgeTerm& e : edges_) {
        const auto [fn, fp] = edge_fluxes(e, x);
        currents.emplace_back(fn.flux, fp.flux);
    }
    return currents;
}

std::vector<double> DriftDiffusion::terminal_currents(std::span<const double> x) const {
    const auto currents = edge_currents(x);
    std::vector<double> terminal(contact_count_, 0.0);
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const std::int32_t ca = contact_[edges_[k].a], cb = contact_[edges_[k].b];
        if (ca == cb) continue;  // both inside one contact, or both in the device
        const double total = currents[k].first + currents[k].second;
        if (ca >= 0) terminal[static_cast<std::size_t>(ca)] += total;  // leaves contact a
        if (cb >= 0) terminal[static_cast<std::size_t>(cb)] -= total;  // flows into contact b
    }
    return terminal;
}

std::vector<double> DriftDiffusion::gate_charges(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    return gates_.charges(x, 3, contact_count_);
}

}  // namespace NiTCAD::assemble
