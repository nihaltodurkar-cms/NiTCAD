#include "terms.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <string_view>
#include <utility>

#include "NiTCAD/assemble/observer.hpp"
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

using Evaluation = SystemData::Terms::Evaluation;

void SystemData::Term::add_frequency(const SystemData&, const SystemState&,
                                     const Terms::FrequencyEvaluation&) const {}

void SystemData::Term::pattern(const SystemData&, const SystemState&,
                               std::vector<linalg::Triplet>&) const {}

physics::EquilibriumProduct SystemData::Terms::product(const SystemData& s,
                                                           const Evaluation& ev, std::size_t i) {
    return s.fermi_dirac_
               ? physics::fermi_dirac_equilibrium_product(s.n_ie_[i], ev.g[i].n, ev.g[i].p)
               : physics::boltzmann_equilibrium_product(s.n_ie_[i]);
}

// R / R0 = rate_scale (R / Ns), d(R / R0)/dn' = rate_scale dR/dn. The scaled product
// E' = E / Ns^2 has dE/dn = Ns dE'/dn'.
physics::EquilibriumProduct SystemData::Terms::physical_product(const SystemData& s,
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

void SystemData::Terms::add_physical(const SystemData& s, const Evaluation& ev,
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
class SystemData::Terms::ContactRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"ohmic contact rows", TermDomain::nodes, all3, all3, true};
    }
    void add(const SystemData& s, const SystemState& st, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.contact_[i] < 0) continue;
            const double psi = ev.x[m * i], n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            ev.f[m * i] = psi - st.psi0_[i];
            ev.f[m * i + 1] = n - st.n0_[i];
            ev.f[m * i + 2] = p - st.p0_[i];
            if (ev.jacobian()) at(s, ev, i, 0, 0) = at(s, ev, i, 1, 1) = at(s, ev, i, 2, 2) = 1.0;
        }
    }
};

// A Schottky contact node (Unit 48, DECISIONS.md SD3), last in the table: its Poisson row becomes
// psi - psi_S, replacing what the Poisson terms added (assigned: the row's Jacobian values are
// cleared first); its continuity rows keep every box term and lose the carriers emitted into the
// metal, -K_n (n - n_B) in the electron row and +K_p (p - p_B) in the hole row (the sign of a
// recombination: a sink of electrons and holes).
class SystemData::Terms::SchottkyRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"Schottky contact rows", TermDomain::nodes, all3, all3, false};
    }
    void add(const SystemData& s, const SystemState& st, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        const std::span<const linalg::Index> rows = st.pattern_.row_offsets();
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.schottky_[i] < 0) continue;
            const double psi = ev.x[m * i], n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            const double kn = s.schottky_kn_[i], kp = s.schottky_kp_[i];
            ev.f[m * i] = psi - st.psi0_[i];
            ev.f[m * i + 1] -= kn * (n - st.n0_[i]);
            ev.f[m * i + 2] += kp * (p - st.p0_[i]);
            if (!ev.jacobian()) continue;
            const auto r = static_cast<std::size_t>(rows[m * i]);
            const auto r1 = static_cast<std::size_t>(rows[m * i + 1]);
            std::fill(ev.values.begin() + static_cast<std::ptrdiff_t>(r),
                      ev.values.begin() + static_cast<std::ptrdiff_t>(r1), 0.0);
            at(s, ev, i, 0, 0) = 1.0;
            at(s, ev, i, 1, 1) -= kn;
            at(s, ev, i, 2, 2) += kp;
        }
    }
};

// An insulator node: no charge, no carriers (n and p held at 0); an electrode holds psi
// (assigned).
class SystemData::Terms::InsulatorRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"insulator rows", TermDomain::nodes, all3, all3, true};
    }
    void add(const SystemData& s, const SystemState& st, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.contact_[i] >= 0 || s.insulator_[i] == 0) continue;
            const double psi = ev.x[m * i], n = ev.x[m * i + 1], p = ev.x[m * i + 2];
            const bool electrode = s.electrode_[i] >= 0;
            ev.f[m * i] = electrode ? psi - st.psi0_[i] : 0.0;
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
class SystemData::Terms::PoissonCharge final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"Poisson charge", TermDomain::nodes, all3, densities, true};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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
class SystemData::Terms::Ionization final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"incomplete ionization", TermDomain::nodes, psi_only, densities, true};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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
class SystemData::Terms::Gate final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"gate", TermDomain::nodes, psi_only, psi_only, false};
    }
    void add(const SystemData& s, const SystemState& st, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i) || !st.gates_.on_gate(i)) continue;
            const double psi = ev.x[m * i];
            ev.f[m * i] += st.gates_.row_term(i, psi);
            if (ev.jacobian()) at(s, ev, i, 0, 0) -= st.gates_.term(i).coupling;
        }
    }
};

// Shockley-Read-Hall recombination: -V R in the electron row, +V R in the hole row.
class SystemData::Terms::Srh final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"SRH recombination", TermDomain::nodes, densities, densities, false};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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
class SystemData::Terms::Auger final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"Auger recombination", TermDomain::nodes, densities, densities, false};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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
class SystemData::Terms::Radiative final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"radiative recombination", TermDomain::nodes, densities, densities, false};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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
class SystemData::Terms::Storage final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"storage", TermDomain::nodes, densities, densities, false};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
        if (ev.step == nullptr) return;
        const std::size_t m = s.m_, w = s.storage_width();  // w 4 with energy balance
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double V = s.volume_[i];
            const NodeStorage ns = s.node_storage(i, ev.x, ev.g);
            const double k = V * ev.step->rate;
            ev.f[m * i + 1] -= k * (ns.n - ev.step->storage[w * i]);
            ev.f[m * i + 2] += k * (ns.p - ev.step->storage[w * i + 1]);
            if (ev.jacobian()) {
                at(s, ev, i, 1, 1) -= k * ns.d_n;
                at(s, ev, i, 2, 2) += k * ns.d_p;
            }
        }
    }
};

// Manufactured-solution sources (verification only; models.hpp ManufacturedSources, TVQ3): +V S/Ns
// in the Poisson row (extra positive charge, as the doping enters PoissonCharge), +V G_n/R0 in the
// electron row and -V G_p/R0 in the hole row (generation, as ImpactGeneration). Independent of
// the state: no Jacobian entries. Contact and electrode rows are boundary rows and are skipped.
// Each source is scaled by the time profile's w at the step's end time (at t = 0 without a step;
// models.hpp, DT-Q3 A); without a profile w = 1 and w y is y exactly.
class SystemData::Terms::Manufactured final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"manufactured source", TermDomain::nodes, all3, {}, false};  // no partials
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
        const ManufacturedSources& ms = *s.manufactured_;
        const std::size_t m = s.m_;
        const double to_rate = s.rate_scale_ / s.Ns_;  // 1 / R0
        const double w = manufactured_weight(ms, ev.step != nullptr ? ev.step->time_s : 0.0);
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.contact_[i] >= 0) continue;
            const double V = s.volume_[i];
            if (s.insulator_[i] != 0) {
                if (s.electrode_[i] < 0 && !ms.charge_cm3.empty()) {
                    ev.f[m * i] += V * ((w * ms.charge_cm3[i]) / s.Ns_);
                }
                continue;
            }
            if (!ms.charge_cm3.empty()) ev.f[m * i] += V * ((w * ms.charge_cm3[i]) / s.Ns_);
            if (!ms.electron_generation.empty()) {
                ev.f[m * i + 1] += V * ((w * ms.electron_generation[i]) * to_rate);
            }
            if (!ms.hole_generation.empty()) ev.f[m * i + 2] -= V * ((w * ms.hole_generation[i]) * to_rate);
        }
    }
};

// The declared interfaces' half-edge fluxes, fixed and trapped charge and surface recombination
// (interface_edges.hpp): the insulator node's Poisson row, then the semiconductor node's rows.
class SystemData::Terms::Interfaces final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"interfaces", TermDomain::interface_edges, all3, all3, false};
    }
    // The insulator's Poisson row reads the semiconductor's n and p.
    void pattern(const SystemData& s, const SystemState&,
                 std::vector<linalg::Triplet>& t) const override {
        interface_pattern(s, t);
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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

// The interface edges' traps in the frequency domain (Unit 22's small-signal terms; Unit 57
// GE3): with their occupancy's dynamics, in place of the steady Interfaces term. Moved from
// small_signal_matrix with its expressions and order.
class SystemData::Terms::InterfacesFrequency final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"interface traps, frequency domain", TermDomain::interface_edges, all3, all3,
                false};
    }
    void add(const SystemData&, const SystemState&, const Evaluation&) const override {
        NITCAD_EXPECTS(false);  // a frequency-domain term has no steady rows
    }
    void add_frequency(const SystemData& d, const SystemState&,
                       const FrequencyEvaluation& ev) const override {
        using Complex = std::complex<double>;
        const std::span<const double> x = ev.x;
        const std::span<Complex> values = ev.values;
        const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> Complex& {
            return values[d.block_[d.m_ * d.m_ * node + d.m_ * r + c]];
        };
        // The traps' weight Ns / s_physical, s_physical = s / t0.
        const Complex inv_weight = ev.s / (d.rate_scale_ * d.Ns_);
        for (std::size_t k = 0; k < d.interfaces_.edges().size(); ++k) {
            const InterfaceEdge& v = d.interfaces_.edges()[k];
            const std::size_t o = v.insulator, sn = v.semiconductor;
            const InterfaceSmallSignal t = d.interfaces_.small_signal(
                v, x[d.m_ * o], x[d.m_ * sn], x[d.m_ * sn + 1], x[d.m_ * sn + 2],
                d.statistics(sn), inv_weight);
            // As assemble's interface terms.
            const EdgeTerm& e = d.edges_[v.edge];
            const bool o_first = e.a == o;
            const std::size_t os = o_first ? e.ab[0] : e.ba[0];
            const std::size_t so[3] = {o_first ? e.ba[0] : e.ab[0], o_first ? e.ba[1] : e.ab[1],
                                       o_first ? e.ba[3] : e.ab[3]};
            if (d.electrode_[o] < 0) {
                at(o, 0, 0) += t.d_flux_insulator[0];
                values[os] += t.d_flux_insulator[1];
                values[d.interface_np_[k].first] += t.d_flux_insulator[2];
                values[d.interface_np_[k].second] += t.d_flux_insulator[3];
            }
            if (d.contact_[sn] >= 0) continue;
            const bool poisson =
                d.schottky_[sn] < 0;  // a Schottky node's Poisson row is Dirichlet
            if (poisson) values[so[0]] += t.d_flux_semiconductor[0];
            values[so[1]] -= t.d_rate[0];
            values[so[2]] += t.d_rate_p[0];
            for (std::size_t c = 1; c < 4; ++c) {
                if (poisson) at(sn, 0, c - 1) += t.d_flux_semiconductor[c];
                at(sn, 1, c - 1) -= t.d_rate[c];
                at(sn, 2, c - 1) += t.d_rate_p[c];
            }
        }
    }
};

// The displacement flux c_e (psi_b - psi_a) of every mesh edge but an interface's, into the
// Poisson rows that are not Dirichlet rows (+ at a, - at b).
class SystemData::Terms::PoissonFlux final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"Poisson flux", TermDomain::mesh_edges, psi_only, psi_only, false};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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
class SystemData::Terms::CarrierFluxes final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"carrier fluxes", TermDomain::mesh_edges, densities, all3, false};
    }
    // Lombardi (Unit 46): the continuity rows of both ends with the psi columns of the edge's
    // E_perp stencil; none without the model.
    void pattern(const SystemData& s, const SystemState&,
                 std::vector<linalg::Triplet>& t) const override {
        surface_pattern(s, t, 2);
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        for (const EdgeTerm& e : s.edges_) {
            if (e.charged || !e.carriers) continue;
            const std::size_t a = e.a, b = e.b;
            SurfaceSlope slope{};
            const auto [fn, fp] = s.edge_fluxes(e, ev.x, ev.g, jacobian ? &slope : nullptr);
            if (jacobian && e.surface) {
                // d flux / d psi_j = d flux / d E_perp * d E_perp / d psi_j (+ at a, - at b).
                for (std::size_t k = e.surface_first; k < e.surface_last; ++k) {
                    const SurfaceStencil& j = s.surface_stencil_[k];
                    const double w = s.surface_gradient(slope.field, j);
                    if (s.contact_[a] < 0) {
                        values[j.pos[0][0]] += slope.n * w;
                        values[j.pos[0][1]] += slope.p * w;
                    }
                    if (s.contact_[b] < 0) {
                        values[j.pos[1][0]] -= slope.n * w;
                        values[j.pos[1][1]] -= slope.p * w;
                    }
                }
            }
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

// Impact ionization (Unit 19; see system.hpp's header comment): the pair generation G at
// each node that ionizes, +V G in its electron row and -V G in its hole row.
class SystemData::Terms::ImpactGeneration final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"impact ionization", TermDomain::global, densities, all3, false};
    }
    // A node's electron generation reads its neighbours' p (through the hole current) and its
    // hole generation their n.
    void pattern(const SystemData& s, const SystemState&,
                 std::vector<linalg::Triplet>& t) const override {
        impact_pattern(s, t);
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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
class SystemData::Terms::LocalTunnelling final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"local band-to-band tunnelling", TermDomain::global, densities, psi_only, false};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
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
class SystemData::Terms::PathTunnelling final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"nonlocal band-to-band tunnelling", TermDomain::global, densities, psi_only,
                false};
    }
    // The start's hole row and the deposit nodes' electron rows read the psi of every column.
    void pattern(const SystemData& s, const SystemState& st,
                 std::vector<linalg::Triplet>& t) const override {
        const std::size_t stride = s.m_;
        for (std::size_t k = 0; k < st.paths_.paths.size(); ++k) {
            const TunnelPath& p = st.paths_.paths[k];
            const PathLayout& l = st.layouts_[k];
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
    void add(const SystemData& s, const SystemState& st, const Evaluation& ev) const override {
        if (st.paths_.paths.empty()) return;
        const std::size_t stride = s.m_;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        PathTerms t;
        for (std::size_t k = 0; k < st.paths_.paths.size(); ++k) {
            const PathLayout& l = st.layouts_[k];
            s.path_terms(st.paths_.paths[k], l, ev.x, t);
            if (t.pairs == 0.0 &&
                std::ranges::all_of(t.d_pairs, [](double v) { return v == 0.0; })) {
                continue;
            }
            const std::size_t nc = l.columns.size();
            const std::size_t i = st.paths_.paths[k].start;
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

// Gate tunnelling (Unit 51, DECISIONS.md GD6): each tunnel's sinks, -T_n in its node's electron
// row and +T_p in its hole row (the signs of the Schottky emission), reading the node's psi, n and
// p and the psi of its path's nodes.
class SystemData::Terms::GateTunnelling final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"gate tunnelling", TermDomain::global, densities, all3, false};
    }
    void pattern(const SystemData& s, const SystemState&,
                 std::vector<linalg::Triplet>& t) const override {
        const std::size_t stride = s.m_;
        for (const GateTunnel& g : s.gate_tunnels_) {
            for (const std::size_t v : g.nodes) {
                for (const std::size_t r : {std::size_t{1}, std::size_t{2}}) {
                    t.push_back({static_cast<linalg::Index>(stride * g.node + r),
                                 static_cast<linalg::Index>(stride * v), 0.0});
                }
            }
        }
    }
    void add(const SystemData& s, const SystemState& st, const Evaluation& ev) const override {
        const std::size_t m = s.m_;
        std::vector<GateTunnelFlux> own;
        if (ev.tunnels == nullptr) s.gate_tunnel_fluxes(st, ev.x, ev.g, own);
        const std::vector<GateTunnelFlux>& fluxes = ev.tunnels != nullptr ? ev.tunnels->finish() : own;
        for (std::size_t t = 0; t < s.gate_tunnels_.size(); ++t) {
            const GateTunnel& g = s.gate_tunnels_[t];
            const GateTunnelFlux& flux = fluxes[t];
            const std::size_t i = g.node;
            ev.f[m * i + 1] -= flux.T_n;
            ev.f[m * i + 2] += flux.T_p;
            if (!ev.jacobian()) continue;
            for (std::size_t c = 0; c < 3; ++c) {
                at(s, ev, i, 1, c) -= flux.d_n[c];
                at(s, ev, i, 2, c) += flux.d_p[c];
            }
            for (std::size_t k = 0; k < g.nodes.size(); ++k) {
                ev.values[g.pos_n[k]] -= flux.d_n[3 + k];
                ev.values[g.pos_p[k]] += flux.d_p[3 + k];
            }
        }
    }
};

std::shared_ptr<const SystemData::Term> SystemData::Terms::poisson_charge() {
    return std::make_shared<PoissonCharge>();
}
std::shared_ptr<const SystemData::Term> SystemData::Terms::gate() {
    return std::make_shared<Gate>();
}
std::shared_ptr<const SystemData::Term> SystemData::Terms::poisson_flux() {
    return std::make_shared<PoissonFlux>();
}

namespace {

template <class T>
std::shared_ptr<const SystemData::Term> make_term() {
    return std::make_shared<T>();
}

}  // namespace

std::span<const CompositionEntry> SystemData::Terms::isothermal_table() {
    using E = EquationSet;
    using M = const PhysicsModels&;
    const auto always = [](const SystemData&, const SystemState&, M) { return true; };
    const auto contacts = [](const SystemData& s, const SystemState&, M) {
        return std::ranges::any_of(s.contact_, [](std::int32_t c) { return c >= 0; });
    };
    const auto insulators = [](const SystemData& s, const SystemState&, M) {
        return std::ranges::any_of(s.insulator_, [](char c) { return c != 0; });
    };
    const auto gated = [](const SystemData& s, const SystemState& st, M) {
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (st.gates_.on_gate(i)) return true;
        }
        return false;
    };
    const auto interfaces = [](const SystemData& s, const SystemState&, M) {
        return !s.interfaces_.edges().empty();
    };
    const auto schottky = [](const SystemData& s, const SystemState&, M) {
        return std::ranges::any_of(s.schottky_, [](std::int32_t c) { return c >= 0; });
    };
    static const CompositionEntry table[] = {
        {+E::contacts, contacts, &make_term<ContactRows>},
        {+E::oxide, insulators, &make_term<InsulatorRows>},
        {+E::poisson, always, &make_term<PoissonCharge>},
        {+E::poisson,
         [](const SystemData&, const SystemState&, M m) { return m.incomplete_ionization; },
         &make_term<Ionization>},
        {+E::gate, gated, &make_term<Gate>},
        {+E::continuity, [](const SystemData&, const SystemState&, M m) { return m.srh; },
         &make_term<Srh>},
        {+E::continuity, [](const SystemData&, const SystemState&, M m) { return m.auger; },
         &make_term<Auger>},
        {+E::continuity, [](const SystemData&, const SystemState&, M m) { return m.radiative; },
         &make_term<Radiative>},
        {+E::continuity, always, &make_term<Storage>},
        {E::poisson | E::continuity,
         [](const SystemData&, const SystemState&, M m) { return m.manufactured != nullptr; },
         &make_term<Manufactured>},
        {+E::energy, [](const SystemData&, const SystemState&, M m) { return m.energy_balance; },
         &energy_rows},
        {+E::interfaces, interfaces, &make_term<Interfaces>},
        {+E::interfaces, interfaces, &make_term<InterfacesFrequency>, true},
        {+E::poisson, always, &make_term<PoissonFlux>},
        {+E::continuity,
         [](const SystemData&, const SystemState&, M m) { return !m.energy_balance; },
         &make_term<CarrierFluxes>},
        {E::continuity | E::energy,
         [](const SystemData&, const SystemState&, M m) { return m.energy_balance; },
         &energy_fluxes},
        {+E::generation,
         [](const SystemData&, const SystemState&, M m) {
             return m.impact_ionization && !m.energy_balance;
         },
         &make_term<ImpactGeneration>},
        {E::generation | E::energy,
         [](const SystemData&, const SystemState&, M m) {
             return m.impact_ionization && m.energy_balance;
         },
         &energy_impact},
        {+E::generation, [](const SystemData&, const SystemState&, M m) { return m.btbt_local; },
         &make_term<LocalTunnelling>},
        {+E::generation,
         [](const SystemData&, const SystemState&, M m) {
             return m.btbt_nonlocal != NonlocalTunnelling::off;
         },
         &make_term<PathTunnelling>},
        {+E::continuity,
         [](const SystemData& s, const SystemState&, M) { return !s.gate_tunnels_.empty(); },
         &make_term<GateTunnelling>},
        {+E::contacts, schottky, &make_term<SchottkyRows>},
    };
    return table;
}

std::span<const std::shared_ptr<const SystemData::Term>> SystemData::Terms::list(
    const System& s) {
    return of(s).composition_->terms();
}

std::span<const std::shared_ptr<const SystemData::Term>> SystemData::Terms::list(
    const SystemData& s) {
    return s.composition_->terms();
}

void SystemData::Terms::interface_pattern(const SystemData& s,
                                              std::vector<linalg::Triplet>& t) {
    const std::size_t stride = s.m_;
    for (const InterfaceEdge& f : s.interfaces_.edges()) {
        for (const std::size_t c : {std::size_t{1}, std::size_t{2}}) {
            t.push_back({static_cast<linalg::Index>(stride * f.insulator),
                         static_cast<linalg::Index>(stride * f.semiconductor + c), 0.0});
        }
    }
}

void SystemData::Terms::impact_pattern(const SystemData& s,
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

void SystemData::Terms::surface_pattern(const SystemData& s,
                                            std::vector<linalg::Triplet>& t, std::size_t rows) {
    const std::size_t stride = s.m_;
    for (const EdgeTerm& e : s.edges_) {
        for (std::size_t k = e.surface_first; k < e.surface_last; ++k) {
            const std::size_t column = stride * s.surface_stencil_[k].node;
            for (const std::size_t node : {e.a, e.b}) {
                for (std::size_t r = 1; r <= rows; ++r) {
                    t.push_back({static_cast<linalg::Index>(stride * node + r),
                                 static_cast<linalg::Index>(column), 0.0});
                }
            }
        }
    }
}

void SystemData::Terms::pattern(const SystemData& s, const SystemState& st,
                                std::vector<linalg::Triplet>& triplets) {
    s.composition_->pattern(s, st, triplets);
}

void SystemData::Terms::evaluate(const System& dd, const Term& t,
                                 std::span<const double> x, const TimeStep* step,
                                 std::span<double> f, linalg::SparseMatrix& jacobian) {
    const SystemData& s = of(dd);
    const SystemState& st = state_of(dd);
    NITCAD_EXPECTS(x.size() == s.unknowns() && f.size() == s.unknowns());
    NITCAD_EXPECTS(jacobian.has_same_pattern(st.pattern_));
    NITCAD_EXPECTS(step == nullptr || (step->storage.size() == s.storage_width() * s.node_count() &&
                                       step->traps.size() == s.trap_slots()));
    const std::span<double> values = jacobian.values();
    std::ranges::fill(f, 0.0);
    std::ranges::fill(values, 0.0);
    if (s.electrothermal()) {
        const std::vector<ThermalState> thermal = s.thermal_states(x);
        const std::vector<NodeLevels> levels = s.thermal_levels(x, thermal);
        t.add(s, st, {x, f, values, step, {}, thermal, levels, true});
        return;
    }
    const std::vector<NodeDegeneracy> g = s.degeneracies(x);
    t.add(s, st, {x, f, values, step, g});
}

// --- Unit 45: rates read from the terms --------------------------------------------------------

namespace {

// The mechanism of a term, by its name (terms.cpp and thermal_terms.cpp), or none.
std::optional<RateMechanism> mechanism_of(std::string_view name) noexcept {
    using M = RateMechanism;
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

std::vector<MechanismRate> SystemData::rates(const SystemState& st, std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    const std::size_t n = node_count();
    // A row's scaled residual r is V R / R0 (rate_scale_ = Ns / R0, generation_scale_ = 1 / R0).
    const double R0 = Ns_ / rate_scale_;
    std::vector<double> f(unknowns());
    std::vector<SystemData::NodeDegeneracy> g;
    std::vector<SystemData::ThermalState> thermal;
    std::vector<SystemData::NodeLevels> levels;
    if (electrothermal()) {
        thermal = thermal_states(x);
        levels = thermal_levels(x, thermal);
    } else {
        g = degeneracies(x);
    }
    std::vector<MechanismRate> out;
    for (const std::shared_ptr<const SystemData::Term>& t : SystemData::Terms::list(*this)) {
        const auto m = mechanism_of(t->info().name);
        if (!m) continue;
        std::ranges::fill(f, 0.0);
        const SystemData::Terms::Evaluation ev{x, f, {}, nullptr, g, thermal, levels, true};
        t->add(*this, st, ev);
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
    for (const MechanismRate& r : out) {  // the trace (observer.hpp; Unit 57, GE1 ii)
        detail::observe("rates.mechanism", static_cast<double>(r.mechanism));
        detail::observe("rates.electron", r.electron);
        detail::observe("rates.hole", r.hole);
    }
    return out;
}

}  // namespace NiTCAD::assemble
