#include "terms.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <string_view>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/band_to_band.hpp"
#include "composition.hpp"
#include "scaled_device.hpp"

namespace NiTCAD::assemble {

namespace {

constexpr FieldId psi_only[] = {FieldId::potential};
constexpr FieldId densities[] = {FieldId::electrons, FieldId::holes};
constexpr FieldId all3[] = {FieldId::potential, FieldId::electrons, FieldId::holes};

}  // namespace

using Evaluation = TermAccess::Evaluation;

void DriftDiffusion::Term::pattern(const DriftDiffusion&, std::vector<linalg::Triplet>&) const {}

physics::EquilibriumProduct DriftDiffusion::Terms::product(const DriftDiffusion& s,
                                                           const Evaluation& ev, std::size_t i) {
    return s.fermi_dirac_
               ? physics::fermi_dirac_equilibrium_product(s.n_ie_[i], ev.g[i].n, ev.g[i].p)
               : physics::boltzmann_equilibrium_product(s.n_ie_[i]);
}

// R / R0 = rate_scale (R / Ns), d(R / R0)/dn' = rate_scale dR/dn. The scaled product
// E' = E / Ns^2 has dE/dn = Ns dE'/dn'.
physics::EquilibriumProduct DriftDiffusion::Terms::physical_product(const DriftDiffusion& s,
                                                                    const Evaluation& ev,
                                                                    std::size_t i) {
    const double nie = s.n_ie_[i] * s.Ns_;
    physics::EquilibriumProduct E = physics::boltzmann_equilibrium_product(nie);
    if (s.fermi_dirac_) {
        const physics::EquilibriumProduct p = product(s, ev, i);
        E = {p.value * s.Ns_ * s.Ns_, p.d_dn * s.Ns_, p.d_dp * s.Ns_};
    }
    return E;
}

void DriftDiffusion::Terms::add_physical(const DriftDiffusion& s, const Evaluation& ev,
                                         std::size_t i, const physics::RecombinationRate& r) {
    const std::size_t m = s.m_;
    const double k = s.volume_[i] * s.rate_scale_;
    ev.f[m * i + 1] -= k * (r.rate / s.Ns_);
    ev.f[m * i + 2] += k * (r.rate / s.Ns_);
    if (ev.jacobian()) {
        at(s, ev, i, 1, 1) -= k * r.d_dn;
        at(s, ev, i, 1, 2) -= k * r.d_dp;
        at(s, ev, i, 2, 1) += k * r.d_dn;
        at(s, ev, i, 2, 2) += k * r.d_dp;
    }
}

// An ohmic contact node: its rows are psi - psi0, n - n0, p - p0 (assigned).
class DriftDiffusion::Terms::ContactRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"ohmic contact rows", TermDomain::nodes, all3, all3, true};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.contact_[i] < 0) continue;
            const double psi = ev.x[m * i], n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            ev.f[m * i] = psi - s.psi0_[i];
            ev.f[m * i + 1] = n - s.n0_[i];
            ev.f[m * i + 2] = p - s.p0_[i];
            if (ev.jacobian()) at(s, ev, i, 0, 0) = at(s, ev, i, 1, 1) = at(s, ev, i, 2, 2) = 1.0;
        }
    }
};

// An insulator node: no charge, no carriers (n and p held at 0); an electrode holds psi
// (assigned).
class DriftDiffusion::Terms::InsulatorRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"insulator rows", TermDomain::nodes, all3, all3, true};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.contact_[i] >= 0 || s.insulator_[i] == 0) continue;
            const double psi = ev.x[m * i], n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            const bool electrode = s.electrode_[i] >= 0;
            ev.f[m * i] = electrode ? psi - s.psi0_[i] : 0.0;
            ev.f[m * i + 1] = n;
            ev.f[m * i + 2] = p;
            if (ev.jacobian()) {
                if (electrode) at(s, ev, i, 0, 0) = 1.0;
                at(s, ev, i, 1, 1) = at(s, ev, i, 2, 2) = 1.0;
            }
        }
    }
};

// A semiconductor node's space charge -V (n - p - N) in its Poisson row; its continuity rows start
// at 0 (assigned).
class DriftDiffusion::Terms::PoissonCharge final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"Poisson charge", TermDomain::nodes, all3, densities, true};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            const double V = s.volume_[i];
            ev.f[m * i] = -V * (n - p - s.doping_[i]);
            ev.f[m * i + 1] = 0.0;
            ev.f[m * i + 2] = 0.0;
            if (ev.jacobian()) {
                at(s, ev, i, 0, 1) = -V;
                at(s, ev, i, 0, 2) = V;
            }
        }
    }
};

// Incomplete ionization: the Poisson row's charge with the ionized dopants, -V (n - p - N+ + N-),
// in place of the Poisson charge's (assigned).
class DriftDiffusion::Terms::Ionization final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"incomplete ionization", TermDomain::nodes, psi_only, densities, true};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        const std::span<const NodeDegeneracy> g = ev.g;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            const double V = s.volume_[i];
            // eta_c = ln(n / n_ie) - ln gamma_n - g_n, d eta_c / dn = 1/n - d ln gamma_n / dn;
            // eta_v likewise from p.
            const double ln = s.fermi_dirac_ ? g[i].n.log_gamma : 0.0;
            const double lp = s.fermi_dirac_ ? g[i].p.log_gamma : 0.0;
            const double wn = s.fermi_dirac_ ? g[i].n.d_density : 0.0;
            const double wp = s.fermi_dirac_ ? g[i].p.d_density : 0.0;
            const double eta_c = s.donors_[i] > 0.0
                                     ? std::log(n / s.n_ie_[i]) - ln - s.log_dos_n_[i]
                                     : 0.0;
            const double eta_v = s.acceptors_[i] > 0.0
                                     ? std::log(p / s.n_ie_[i]) - lp - s.log_dos_p_[i]
                                     : 0.0;
            const detail::IonizedCharge c = detail::ionized_charge(
                s.donors_[i], s.acceptors_[i], eta_c, eta_v, s.levels_[i]);
            ev.f[m * i] = -V * (n - p - c.value);
            if (ev.jacobian()) {
                at(s, ev, i, 0, 1) = -V * (1.0 - c.d_eta_c * (1.0 / n - wn));
                at(s, ev, i, 0, 2) = V * (1.0 + c.d_eta_v * (1.0 / p - wp));
            }
        }
    }
};

// A gate node's coupling to its gate electrode in its Poisson row (gate.hpp).
class DriftDiffusion::Terms::Gate final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"gate", TermDomain::nodes, psi_only, psi_only, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i) || !s.gates_.on_gate(i)) continue;
            const double psi = ev.x[m * i];
            ev.f[m * i] += s.gates_.row_term(i, psi);
            if (ev.jacobian()) at(s, ev, i, 0, 0) -= s.gates_.term(i).coupling;
        }
    }
};

// Shockley-Read-Hall recombination: -V R in the electron row, +V R in the hole row.
class DriftDiffusion::Terms::Srh final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"SRH recombination", TermDomain::nodes, densities, densities, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            const double V = s.volume_[i];
            const physics::RecombinationRate r = physics::srh_recombination(
                n, p, product(s, ev, i), s.n_ie_[i], s.tau_n_[i], s.tau_p_[i]);
            const double k = V * s.rate_scale_;
            ev.f[m * i + 1] -= k * r.rate;
            ev.f[m * i + 2] += k * r.rate;
            if (ev.jacobian()) {
                at(s, ev, i, 1, 1) -= k * r.d_dn;
                at(s, ev, i, 1, 2) -= k * r.d_dp;
                at(s, ev, i, 2, 1) += k * r.d_dn;
                at(s, ev, i, 2, 2) += k * r.d_dp;
            }
        }
    }
};

// Auger recombination (physical densities, add_physical).
class DriftDiffusion::Terms::Auger final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"Auger recombination", TermDomain::nodes, densities, densities, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            add_physical(s, ev, i,
                         physics::auger_recombination(n * s.Ns_, p * s.Ns_,
                                                      physical_product(s, ev, i), s.auger_n_[i],
                                                      s.auger_p_[i]));
        }
    }
};

// Radiative recombination on the nodes whose material has a coefficient (physical densities,
// add_physical).
class DriftDiffusion::Terms::Radiative final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"radiative recombination", TermDomain::nodes, densities, densities, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i) || !(s.radiative_[i] > 0.0)) continue;
            const double n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            add_physical(s, ev, i,
                         physics::radiative_recombination(n * s.Ns_, p * s.Ns_,
                                                          physical_product(s, ev, i),
                                                          s.radiative_[i]));
        }
    }
};

// The BDF storage of a time step, -V rate (S_n - S_n,old) in the electron row and
// +V rate (S_p - S_p,old) in the hole row; nothing without a step.
class DriftDiffusion::Terms::Storage final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"storage", TermDomain::nodes, densities, densities, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        if (ev.step == nullptr) return;
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double V = s.volume_[i];
            const NodeStorage st = s.node_storage(i, ev.x, ev.g);
            const double k = V * ev.step->rate;
            ev.f[m * i + 1] -= k * (st.n - ev.step->storage[2 * i]);
            ev.f[m * i + 2] += k * (st.p - ev.step->storage[2 * i + 1]);
            if (ev.jacobian()) {
                at(s, ev, i, 1, 1) -= k * st.d_n;
                at(s, ev, i, 2, 2) += k * st.d_p;
            }
        }
    }
};

// The declared interfaces' half-edge fluxes, fixed and trapped charge and surface recombination
// (interface_edges.hpp): the insulator node's Poisson row, then the semiconductor node's rows.
class DriftDiffusion::Terms::Interfaces final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"interfaces", TermDomain::interface_edges, all3, all3, false};
    }
    // The insulator's Poisson row reads the semiconductor's n and p.
    void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t) const override {
        interface_pattern(s, t);
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        for (std::size_t k = 0; k < s.interfaces_.edges().size(); ++k) {
            const InterfaceEdge& v = s.interfaces_.edges()[k];
            const std::size_t o = v.insulator, sn = v.semiconductor;
            const InterfaceDrift t = s.interface_at(k, ev.x, ev.step);
            // Edge positions: row of o with column psi of sn, and the rows of sn with column psi
            // of o.
            const EdgeTerm& e = s.edges_[v.edge];
            const bool o_first = e.a == o;
            const std::size_t os = o_first ? e.ab[0] : e.ba[0];
            const std::size_t so[3] = {o_first ? e.ba[0] : e.ab[0], o_first ? e.ba[1] : e.ab[1],
                                       o_first ? e.ba[3] : e.ab[3]};
            if (s.electrode_[o] < 0) {  // the insulator node's Poisson row
                f[m * o] += t.flux_insulator;
                if (jacobian) {
                    at(s, ev, o, 0, 0) += t.d_flux_insulator[0];
                    values[os] += t.d_flux_insulator[1];
                    values[s.interface_np_[k].first] += t.d_flux_insulator[2];
                    values[s.interface_np_[k].second] += t.d_flux_insulator[3];
                }
            }
            if (s.contact_[sn] >= 0) continue;  // a Dirichlet semiconductor node
            f[m * sn] += t.flux_semiconductor;
            f[m * sn + 1] -= t.rate;
            f[m * sn + 2] += t.rate_p;
            if (jacobian) {
                values[so[0]] += t.d_flux_semiconductor[0];
                values[so[1]] -= t.d_rate[0];
                values[so[2]] += t.d_rate_p[0];
                for (std::size_t c = 1; c < 4; ++c) {
                    at(s, ev, sn, 0, c - 1) += t.d_flux_semiconductor[c];
                    at(s, ev, sn, 1, c - 1) -= t.d_rate[c];
                    at(s, ev, sn, 2, c - 1) += t.d_rate_p[c];
                }
            }
        }
    }
};

// The displacement flux c_e (psi_b - psi_a) of every mesh edge but an interface's, into the
// Poisson rows that are not Dirichlet rows (+ at a, - at b).
class DriftDiffusion::Terms::PoissonFlux final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"Poisson flux", TermDomain::mesh_edges, psi_only, psi_only, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        for (const EdgeTerm& e : s.edges_) {
            const std::size_t a = e.a, b = e.b;
            const double psi_a = ev.x[m * a], psi_b = ev.x[m * b];
            const double poisson = e.c * (psi_b - psi_a);
            if (e.charged) continue;  // the interface's half-edge fluxes (Interfaces)
            if (!e.carriers) {  // into the rows that are not Dirichlet
                if (s.contact_[a] < 0 && s.electrode_[a] < 0) {
                    f[m * a] += poisson;
                    if (jacobian) {
                        at(s, ev, a, 0, 0) -= e.c;
                        values[e.ab[0]] += e.c;
                    }
                }
                if (s.contact_[b] < 0 && s.electrode_[b] < 0) {
                    f[m * b] -= poisson;
                    if (jacobian) {
                        at(s, ev, b, 0, 0) -= e.c;
                        values[e.ba[0]] += e.c;
                    }
                }
                continue;
            }
            if (s.contact_[a] < 0) {
                f[m * a] += poisson;
                if (jacobian) {
                    at(s, ev, a, 0, 0) -= e.c;
                    values[e.ab[0]] += e.c;
                }
            }
            if (s.contact_[b] < 0) {
                f[m * b] -= poisson;
                if (jacobian) {
                    at(s, ev, b, 0, 0) -= e.c;
                    values[e.ba[0]] += e.c;
                }
            }
        }
    }
};

// The Scharfetter-Gummel electron and hole fluxes of every carrier edge but an interface's
// (edge_fluxes: field mobility, Fermi-Dirac, thermionic emission), into the continuity rows that
// are not Dirichlet rows (+ at a, - at b). One term for both carriers: edge_fluxes computes them
// together, and they write disjoint rows.
class DriftDiffusion::Terms::CarrierFluxes final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"carrier fluxes", TermDomain::mesh_edges, densities, all3, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        for (const EdgeTerm& e : s.edges_) {
            if (e.charged || !e.carriers) continue;
            const std::size_t a = e.a, b = e.b;
            const auto [fn, fp] = s.edge_fluxes(e, ev.x, ev.g);
            if (s.contact_[a] < 0) {
                f[m * a + 1] += fn.flux;
                f[m * a + 2] += fp.flux;
                if (jacobian) {
                    at(s, ev, a, 1, 0) += fn.d_psi1;
                    values[e.ab[1]] += fn.d_psi2;
                    at(s, ev, a, 1, 1) += fn.d_c1;
                    values[e.ab[2]] += fn.d_c2;
                    at(s, ev, a, 2, 0) += fp.d_psi1;
                    values[e.ab[3]] += fp.d_psi2;
                    at(s, ev, a, 2, 2) += fp.d_c1;
                    values[e.ab[4]] += fp.d_c2;
                }
            }
            if (s.contact_[b] < 0) {
                f[m * b + 1] -= fn.flux;
                f[m * b + 2] -= fp.flux;
                if (jacobian) {
                    values[e.ba[1]] -= fn.d_psi1;
                    at(s, ev, b, 1, 0) -= fn.d_psi2;
                    values[e.ba[2]] -= fn.d_c1;
                    at(s, ev, b, 1, 1) -= fn.d_c2;
                    values[e.ba[3]] -= fp.d_psi1;
                    at(s, ev, b, 2, 0) -= fp.d_psi2;
                    values[e.ba[4]] -= fp.d_c1;
                    at(s, ev, b, 2, 2) -= fp.d_c2;
                }
            }
        }
    }
};

// Impact ionization (Unit 19; see drift_diffusion.hpp's header comment): the pair generation G at
// each node that ionizes, +V G in its electron row and -V G in its hole row.
class DriftDiffusion::Terms::ImpactGeneration final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"impact ionization", TermDomain::global, densities, all3, false};
    }
    // A node's electron generation reads its neighbours' p (through the hole current) and its
    // hole generation their n.
    void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t) const override {
        impact_pattern(s, t);
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        if (s.impact_nodes_.empty()) return;
        const std::size_t stride = s.m_;
        const std::span<const double> x = ev.x;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const int D = s.dimension_;
        // The carrier edges' fluxes, once.
        std::vector<std::pair<EdgeFlux, EdgeFlux>> flux(s.edges_.size());
        for (std::size_t k = 0; k < s.edges_.size(); ++k) {
            if (s.edges_[k].carriers) flux[k] = s.edge_fluxes(s.edges_[k], x, ev.g);
        }
        constexpr double q = base::q_C;
        constexpr double floor2 = impact_current_floor * impact_current_floor;
        const double rho2 = s.impact_resolution_ * s.impact_resolution_;
        for (const ImpactNode& node : s.impact_nodes_) {
            const std::size_t i = node.node;
            const double count = static_cast<double>(node.last - node.first);
            // The field and the two current densities at the node, reconstructed from their edge
            // components by least squares: v = P sum_k t_k v_k, P the pseudo-inverse of
            // sum t t^T. Edge components: E_k = -field (psi_b - psi_a), j_k = current * flux.
            // And the size of the opposing flux terms whose difference each current is, averaged
            // over the edges: R_c = mean_k current_k a_k (c_a + c_b) sqrt(1 + Delta_k^2), a the
            // edge's low-field factor, c the carrier's density, Delta = psi_b - psi_a.
            double E[3] = {}, jn[3] = {}, jp[3] = {}, R[2] = {0.0, 0.0};
            for (std::size_t k = node.first; k < node.last; ++k) {
                const ImpactEdge& ie = s.impact_edges_[k];
                const EdgeTerm& e = s.edges_[ie.edge];
                const double delta = x[stride * e.b] - x[stride * e.a];
                const double Ek = -ie.field * delta;
                const double jnk = ie.current * flux[ie.edge].first.flux;
                const double jpk = ie.current * flux[ie.edge].second.flux;
                const double root = std::sqrt(1.0 + delta * delta);
                R[0] += ie.current * e.an * (x[stride * e.a + 1] + x[stride * e.b + 1]) *
                        root / count;
                R[1] += ie.current * e.ap * (x[stride * e.a + 2] + x[stride * e.b + 2]) *
                        root / count;
                for (int r = 0; r < D; ++r) {
                    double w = 0.0;  // (P t_k)_r
                    for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * ie.t[c];
                    E[r] += w * Ek;
                    jn[r] += w * jnk;
                    jp[r] += w * jpk;
                }
            }
            // Per carrier, with m = sqrt(|j|^2 + eps^2), eps^2 = floor^2 + (rho R)^2:
            //     E_par = |E . j| / m,   G = alpha(E_par) |j|^2 / m / q.
            // E_par is the field along the current (|E . j^| for |j| >> eps), at most |E|;
            // |j|^2 / m is |j| for |j| >> eps, at most |j|, and smooth through j = 0, where it and
            // its gradient vanish. eps is the current's resolution at the node: rho =
            // models.impact_current_resolution of the opposing flux terms, above
            // impact_current_floor; a current at or below it (rounding, in the dense spill-over
            // layer at a junction, where it is a difference of fluxes 1e14 times larger and flips
            // sign between Newton iterates) does not ionize and has no kink.
            // Partials, s = sign(E . j), w = |j|^2 / m:
            //     dG/dE   = alpha' w s j / m / q,
            //     dG/dj   = (alpha' w (s E / m - |E . j| j / m^3)
            //                + alpha j (|j|^2 + 2 eps^2) / m^3) / q,
            //     dG/deps = -(alpha' w |E . j| + alpha |j|^2) eps / m^3 / q,
            //     deps/dR = rho^2 R / eps.
            double G = 0.0, dE[3] = {}, dj[2][3] = {}, dR[2] = {0.0, 0.0};
            for (int c = 0; c < 2; ++c) {
                const double* j = c == 0 ? jn : jp;
                double mag2 = 0.0, Ej = 0.0;
                for (int r = 0; r < D; ++r) {
                    mag2 += j[r] * j[r];
                    Ej += E[r] * j[r];
                }
                const double eps2 = floor2 + rho2 * R[c] * R[c];
                const double m = std::sqrt(mag2 + eps2);
                const double Epar = std::abs(Ej) / m;
                const physics::ImpactIonizationRate a = physics::impact_ionization_coefficient(
                    c == 0 ? node.n : node.p, node.gamma, Epar);
                if (a.alpha == 0.0 && a.d_dE == 0.0) continue;
                const double w = mag2 / m;
                G += a.alpha * w / q;
                const double sign = Ej > 0.0 ? 1.0 : Ej < 0.0 ? -1.0 : 0.0;
                const double m3 = m * m * m;
                for (int r = 0; r < D; ++r) {
                    dE[r] += a.d_dE * w * sign * j[r] / m / q;
                    dj[c][r] = (a.d_dE * w * (sign * E[r] / m - std::abs(Ej) * j[r] / m3) +
                                a.alpha * j[r] * (mag2 + 2.0 * eps2) / m3) /
                               q;
                }
                // dG/deps deps/dR = -(alpha' w |E . j| + alpha |j|^2) / m^3 / q * rho^2 R.
                dR[c] = -(a.d_dE * w * std::abs(Ej) + a.alpha * mag2) / m3 / q * rho2 * R[c];
            }
            const double k_row = s.volume_[i] * s.generation_scale_;
            f[stride * i + 1] += k_row * G;
            f[stride * i + 2] -= k_row * G;
            if (!jacobian) continue;
            const auto add_entry = [&](std::size_t node_j, int comp, std::size_t k, double v) {
                // Row n of node i gets +v, row p -v, at column comp of node_j through edge k.
                const ImpactEdge& ie = s.impact_edges_[k];
                const EdgeTerm& e = s.edges_[ie.edge];
                const bool own = node_j == i;
                const bool i_is_a = e.a == i;
                for (const int row : {1, 2}) {
                    const double sign = row == 1 ? k_row : -k_row;
                    std::size_t pos;
                    if (own) {
                        pos = s.block_[stride * stride * i +
                                       stride * static_cast<std::size_t>(row) +
                                       static_cast<std::size_t>(comp)];
                    } else {
                        const std::size_t* p = i_is_a ? e.ab : e.ba;
                        // Edge positions: (psi,psi), (n,psi), (n,n), (p,psi), (p,p).
                        if (row == 1) {
                            pos = comp == 0 ? p[1] : comp == 1 ? p[2] : ie.n_p;
                        } else {
                            pos = comp == 0 ? p[3] : comp == 1 ? ie.p_n : p[4];
                        }
                    }
                    values[pos] += sign * v;
                }
            };
            for (std::size_t k = node.first; k < node.last; ++k) {
                const ImpactEdge& ie = s.impact_edges_[k];
                const EdgeTerm& e = s.edges_[ie.edge];
                const EdgeFlux& fn = flux[ie.edge].first;
                const EdgeFlux& fp = flux[ie.edge].second;
                double wE = 0.0, wn = 0.0, wp = 0.0;  // dG / d(edge component)
                for (int r = 0; r < D; ++r) {
                    double w = 0.0;
                    for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * ie.t[c];
                    wE += w * dE[r];
                    wn += w * dj[0][r];
                    wp += w * dj[1][r];
                }
                wn *= ie.current;
                wp *= ie.current;
                // R's own partials through this edge: d/dc_a = d/dc_b = current a root / count,
                // d/dpsi_b = -d/dpsi_a = current a (c_a + c_b) Delta / root / count.
                const double delta = x[stride * e.b] - x[stride * e.a];
                const double root = std::sqrt(1.0 + delta * delta);
                const double rn = ie.current * e.an / count, rp = ie.current * e.ap / count;
                const double sn = x[stride * e.a + 1] + x[stride * e.b + 1];
                const double sp = x[stride * e.a + 2] + x[stride * e.b + 2];
                const double dpsi = (dR[0] * rn * sn + dR[1] * rp * sp) * delta / root;
                // E_k = -field (psi_b - psi_a); the fluxes' own partials.
                add_entry(e.a, 0, k, wE * ie.field + wn * fn.d_psi1 + wp * fp.d_psi1 - dpsi);
                add_entry(e.b, 0, k, -wE * ie.field + wn * fn.d_psi2 + wp * fp.d_psi2 + dpsi);
                add_entry(e.a, 1, k, wn * fn.d_c1 + dR[0] * rn * root);
                add_entry(e.b, 1, k, wn * fn.d_c2 + dR[0] * rn * root);
                add_entry(e.a, 2, k, wp * fp.d_c1 + dR[1] * rp * root);
                add_entry(e.b, 2, k, wp * fp.d_c2 + dR[1] * rp * root);
            }
        }
    }
};

// Local band-to-band tunnelling (Unit 20): Kane's G(|E|) at each node with a Kane pair, the field
// by least squares from its edges, +V G in its electron row and -V G in its hole row.
class DriftDiffusion::Terms::LocalTunnelling final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"local band-to-band tunnelling", TermDomain::global, densities, psi_only, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        if (s.tunnel_nodes_.empty()) return;
        const std::size_t stride = s.m_;
        const std::span<const double> x = ev.x;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const int D = s.dimension_;
        for (const TunnelNode& node : s.tunnel_nodes_) {
            const std::size_t i = node.node;
            // The field at the node by least squares from its edges' components
            // E_k = -field (psi_b - psi_a): E = P sum_k t_k E_k.
            double E[3] = {};
            for (std::size_t k = node.first; k < node.last; ++k) {
                const FieldEdge& fe = s.field_edges_[k];
                const EdgeTerm& e = s.edges_[fe.edge];
                const double Ek = -fe.field * (x[stride * e.b] - x[stride * e.a]);
                for (int r = 0; r < D; ++r) {
                    double w = 0.0;  // (P t_k)_r
                    for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * fe.t[c];
                    E[r] += w * Ek;
                }
            }
            double F2 = 0.0;
            for (int r = 0; r < D; ++r) F2 += E[r] * E[r];
            const double F = std::sqrt(F2);
            const physics::KaneRate g = physics::kane_generation(node.A, node.B, F);
            if (g.rate == 0.0 && g.d_dF == 0.0) continue;
            const double k_row = s.volume_[i] * s.generation_scale_;
            f[stride * i + 1] += k_row * g.rate;
            f[stride * i + 2] -= k_row * g.rate;
            if (!jacobian) continue;
            // dG/dE_r = G' E_r / F (F > 0 here: G vanishes at F = 0); an edge's component enters
            // through (P t_k), and E_k through its ends' psi.
            for (std::size_t k = node.first; k < node.last; ++k) {
                const FieldEdge& fe = s.field_edges_[k];
                const EdgeTerm& e = s.edges_[fe.edge];
                double dEk = 0.0;  // dG / dE_k
                for (int r = 0; r < D; ++r) {
                    double w = 0.0;
                    for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * fe.t[c];
                    dEk += w * g.d_dF * E[r] / F;
                }
                const double d_b = -dEk * fe.field;  // dG / dpsi_b; dG / dpsi_a = -d_b
                const bool i_is_a = e.a == i;
                const std::size_t* p = i_is_a ? e.ab : e.ba;  // row i, the other end's psi column
                const double own = i_is_a ? -d_b : d_b, other = i_is_a ? d_b : -d_b;
                values[s.block_[stride * stride * i + stride]] += k_row * own;      // (n, psi) of i
                values[s.block_[stride * stride * i + 2 * stride]] -= k_row * own;  // (p, psi) of i
                values[p[1]] += k_row * other;                          // (n, psi) of the other end
                values[p[3]] -= k_row * other;                          // (p, psi) of the other end
            }
        }
    }
};

// Nonlocal band-to-band tunnelling along the frozen paths (Unit 20): each path's pairs P leave its
// start's hole row and arrive, shared, in its deposit nodes' electron rows.
class DriftDiffusion::Terms::PathTunnelling final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"nonlocal band-to-band tunnelling", TermDomain::global, densities, psi_only,
                false};
    }
    // The start's hole row and the deposit nodes' electron rows read the psi of every column.
    void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t) const override {
        const std::size_t stride = s.m_;
        for (std::size_t k = 0; k < s.paths_.paths.size(); ++k) {
            const TunnelPath& p = s.paths_.paths[k];
            const PathLayout& l = s.layouts_[k];
            for (const std::size_t c : l.columns) {
                t.push_back({static_cast<linalg::Index>(stride * p.start + 2),
                             static_cast<linalg::Index>(stride * c), 0.0});
                for (const std::size_t d : l.deposit) {
                    t.push_back({static_cast<linalg::Index>(stride * d + 1),
                                 static_cast<linalg::Index>(stride * c), 0.0});
                }
            }
        }
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        if (s.paths_.paths.empty()) return;
        const std::size_t stride = s.m_;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        PathTerms t;
        for (std::size_t k = 0; k < s.paths_.paths.size(); ++k) {
            const PathLayout& l = s.layouts_[k];
            s.path_terms(s.paths_.paths[k], l, ev.x, t);
            if (t.pairs == 0.0 &&
                std::ranges::all_of(t.d_pairs, [](double v) { return v == 0.0; })) {
                continue;
            }
            const std::size_t nc = l.columns.size();
            const std::size_t i = s.paths_.paths[k].start;
            f[stride * i + 2] -= t.pairs;
            for (std::size_t d = 0; d < l.deposit.size(); ++d) {
                f[stride * l.deposit[d] + 1] += t.share[d] * t.pairs;
            }
            if (!jacobian) continue;
            for (std::size_t c = 0; c < nc; ++c) {
                values[l.hole[c]] -= t.d_pairs[c];
                for (std::size_t d = 0; d < l.deposit.size(); ++d) {
                    values[l.electron[d * nc + c]] +=
                        t.share[d] * t.d_pairs[c] + t.pairs * t.d_share[d * nc + c];
                }
            }
        }
    }
};

std::shared_ptr<const DriftDiffusion::Term> DriftDiffusion::Terms::poisson_charge() {
    return std::make_shared<PoissonCharge>();
}

namespace {

template <class T>
std::shared_ptr<const TermAccess::Term> make_term() {
    return std::make_shared<T>();
}

}  // namespace

std::span<const CompositionEntry> DriftDiffusion::Terms::isothermal_table() {
    using E = EquationSet;
    using M = const PhysicsModels&;
    const auto always = [](const DriftDiffusion&, M) { return true; };
    const auto contacts = [](const DriftDiffusion& s, M) {
        return std::ranges::any_of(s.contact_, [](std::int32_t c) { return c >= 0; });
    };
    const auto insulators = [](const DriftDiffusion& s, M) {
        return std::ranges::any_of(s.insulator_, [](char c) { return c != 0; });
    };
    const auto gated = [](const DriftDiffusion& s, M) {
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.gates_.on_gate(i)) return true;
        }
        return false;
    };
    const auto interfaces = [](const DriftDiffusion& s, M) {
        return !s.interfaces_.edges().empty();
    };
    static const CompositionEntry table[] = {
        {+E::contacts, contacts, &make_term<ContactRows>},
        {+E::oxide, insulators, &make_term<InsulatorRows>},
        {+E::poisson, always, &make_term<PoissonCharge>},
        {+E::poisson, [](const DriftDiffusion&, M m) { return m.incomplete_ionization; },
         &make_term<Ionization>},
        {+E::gate, gated, &make_term<Gate>},
        {+E::continuity, [](const DriftDiffusion&, M m) { return m.srh; }, &make_term<Srh>},
        {+E::continuity, [](const DriftDiffusion&, M m) { return m.auger; },
         &make_term<Auger>},
        {+E::continuity, [](const DriftDiffusion&, M m) { return m.radiative; },
         &make_term<Radiative>},
        {+E::continuity, always, &make_term<Storage>},
        {+E::interfaces, interfaces, &make_term<Interfaces>},
        {+E::poisson, always, &make_term<PoissonFlux>},
        {+E::continuity, always, &make_term<CarrierFluxes>},
        {+E::generation, [](const DriftDiffusion&, M m) { return m.impact_ionization; },
         &make_term<ImpactGeneration>},
        {+E::generation, [](const DriftDiffusion&, M m) { return m.btbt_local; },
         &make_term<LocalTunnelling>},
        {+E::generation,
         [](const DriftDiffusion&, M m) { return m.btbt_nonlocal != NonlocalTunnelling::off; },
         &make_term<PathTunnelling>},
    };
    return table;
}

std::span<const std::shared_ptr<const DriftDiffusion::Term>> DriftDiffusion::Terms::list(
    const DriftDiffusion& s) {
    return s.composition_->terms();
}

void DriftDiffusion::Terms::interface_pattern(const DriftDiffusion& s,
                                              std::vector<linalg::Triplet>& t) {
    const std::size_t stride = s.m_;
    for (const InterfaceEdge& f : s.interfaces_.edges()) {
        for (const std::size_t c : {std::size_t{1}, std::size_t{2}}) {
            t.push_back({static_cast<linalg::Index>(stride * f.insulator),
                         static_cast<linalg::Index>(stride * f.semiconductor + c), 0.0});
        }
    }
}

void DriftDiffusion::Terms::impact_pattern(const DriftDiffusion& s,
                                           std::vector<linalg::Triplet>& t) {
    const std::size_t stride = s.m_;
    for (const EdgeTerm& e : s.edges_) {
        if (!e.carriers) continue;
        for (const auto [a, b] : {std::pair{e.a, e.b}, std::pair{e.b, e.a}}) {
            t.push_back({static_cast<linalg::Index>(stride * a + 1),
                         static_cast<linalg::Index>(stride * b + 2), 0.0});
            t.push_back({static_cast<linalg::Index>(stride * a + 2),
                         static_cast<linalg::Index>(stride * b + 1), 0.0});
        }
    }
}

void DriftDiffusion::Terms::pattern(const DriftDiffusion& s,
                                    std::vector<linalg::Triplet>& triplets) {
    s.composition_->pattern(s, triplets);
}

void DriftDiffusion::Terms::evaluate(const DriftDiffusion& s, const Term& t,
                                     std::span<const double> x, const TimeStep* step,
                                     std::span<double> f, linalg::SparseMatrix& jacobian) {
    NITCAD_EXPECTS(x.size() == s.unknowns() && f.size() == s.unknowns());
    NITCAD_EXPECTS(jacobian.has_same_pattern(s.pattern_));
    const std::size_t per_node = s.electrothermal() ? 3 : 2;  // storage entries
    NITCAD_EXPECTS(step == nullptr || (step->storage.size() == per_node * s.node_count() &&
                                       step->traps.size() == s.trap_slots()));
    const std::span<double> values = jacobian.values();
    std::ranges::fill(f, 0.0);
    std::ranges::fill(values, 0.0);
    if (s.electrothermal()) {
        const std::vector<ThermalState> thermal = s.thermal_states(x);
        const std::vector<NodeLevels> levels = s.thermal_levels(x, thermal);
        t.add(s, {x, f, values, step, {}, thermal, levels, true});
        return;
    }
    const std::vector<NodeDegeneracy> g = s.degeneracies(x);
    t.add(s, {x, f, values, step, g});
}

// --- Unit 45: rates read from the terms --------------------------------------------------------

namespace {

// The mechanism of a term, by its name (terms.cpp and thermal_terms.cpp), or none.
std::optional<DriftDiffusion::RateMechanism> mechanism_of(std::string_view name) noexcept {
    using M = DriftDiffusion::RateMechanism;
    if (name.starts_with("SRH recombination")) return M::srh;
    if (name.starts_with("Auger recombination")) return M::auger;
    if (name.starts_with("radiative recombination")) return M::radiative;
    if (name.starts_with("interfaces")) return M::interfaces;
    if (name.starts_with("impact ionization")) return M::impact;
    if (name == "local band-to-band tunnelling" || name == "nonlocal band-to-band tunnelling") {
        return M::band_to_band;
    }
    return std::nullopt;
}

}  // namespace

std::vector<DriftDiffusion::MechanismRate> DriftDiffusion::rates(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    const std::size_t n = node_count();
    // A row's scaled residual r is V R / R0 (rate_scale_ = Ns / R0, generation_scale_ = 1 / R0).
    const double R0 = Ns_ / rate_scale_;
    std::vector<double> f(unknowns());
    std::vector<NodeDegeneracy> g;
    std::vector<ThermalState> thermal;
    std::vector<NodeLevels> levels;
    if (electrothermal()) {
        thermal = thermal_states(x);
        levels = thermal_levels(x, thermal);
    } else {
        g = degeneracies(x);
    }
    std::vector<MechanismRate> out;
    for (const std::shared_ptr<const Term>& t : Terms::list(*this)) {
        const auto m = mechanism_of(t->info().name);
        if (!m) continue;
        std::ranges::fill(f, 0.0);
        const Terms::Evaluation ev{x, f, {}, nullptr, g, thermal, levels, true};
        t->add(*this, ev);
        auto it = std::ranges::find(out, *m, &MechanismRate::mechanism);
        if (it == out.end()) {
            out.push_back({*m, std::vector<double>(n, 0.0), std::vector<double>(n, 0.0)});
            it = out.end() - 1;
        }
        for (std::size_t i = 0; i < n; ++i) {
            const double k = R0 / volume_[i];
            it->electron[i] -= f[m_ * i + 1] * k;  // the electron row subtracts V R / R0
            it->hole[i] += f[m_ * i + 2] * k;      // the hole row adds it
        }
    }
    return out;
}

}  // namespace NiTCAD::assemble
