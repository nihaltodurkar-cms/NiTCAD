// Thermal equilibrium as a composition (Unit 57, DECISIONS.md GE5; equilibrium.hpp): the data,
// the bias, equilibrium's own term variants and its composition table, and its queries. The
// variants keep the expressions and the order of additions of the assembler before them
// (EquilibriumPoisson): per entry the charge, the gate, the interfaces, the edge flux, then the
// Dirichlet rows.
#include "NiTCAD/assemble/equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/ohmic.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "composition.hpp"
#include "scaled_device.hpp"
#include "terms.hpp"

namespace NiTCAD::assemble {

namespace {

constexpr FieldId psi_only[] = {FieldId::potential};

using Evaluation = SystemData::Terms::Evaluation;

template <class T>
std::shared_ptr<const SystemData::Term> make_term() {
    return std::make_shared<T>();
}

}  // namespace

// --- the data and the bias ---------------------------------------------------------------------

std::expected<SystemData::Made, base::Error> SystemData::make_equilibrium(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models) {
    auto scaled = detail::make_scaled_device(device, scaling, models);
    if (!scaled) return std::unexpected(std::move(scaled.error()));
    const std::size_t n = scaled->volume.size();

    Made made{std::make_shared<SystemData>(), {}};
    SystemData& s = *made.data;
    SystemState& st = made.state;
    s.equilibrium_ = true;
    s.volume_ = std::move(scaled->volume);
    s.doping_ = std::move(scaled->doping);
    s.n_ie_ = std::move(scaled->n_ie);
    s.log_dos_n_ = std::move(scaled->log_dos_n);
    s.log_dos_p_ = std::move(scaled->log_dos_p);
    s.band_shift_ = std::move(scaled->band_shift);
    s.donors_ = std::move(scaled->donors);
    s.acceptors_ = std::move(scaled->acceptors);
    s.levels_ = std::move(scaled->levels);
    s.fermi_dirac_ = models.fermi_dirac;
    s.ionization_ = models.incomplete_ionization;
    s.contact_ = std::move(scaled->contact);
    s.schottky_ = std::move(scaled->schottky);
    s.schottky_eta_ = std::move(scaled->schottky_eta);
    s.insulator_ = std::move(scaled->insulator);
    s.electrode_ = std::move(scaled->electrode);
    s.electrode_potential_ = std::move(scaled->electrode_potential);
    s.interfaces_ = std::move(scaled->interfaces);
    s.kinds_ = std::move(scaled->kinds);
    s.contact_count_ = device.contacts().size();
    s.temperature_K_ = scaling.temperature_K;
    s.Ns_ = scaling.Ns;
    s.V_T_ = scaling.V_T;
    st.gates_ = std::move(scaled->gates);

    // The Dirichlet values: an electrode's at zero bias, a Schottky contact's psi_S at zero bias
    // (Unit 48, SD2, SD6), an ohmic contact's at zero bias.
    st.psi0_.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (s.electrode_[i] >= 0) st.psi0_[i] = s.electrode_potential_[i];
        if (s.schottky_[i] >= 0) {
            st.psi0_[i] = s.schottky_eta_[i] - s.band_shift_[i];
            continue;
        }
        if (s.contact_[i] < 0) continue;
        const physics::NeutralEquilibrium e = detail::neutral_equilibrium(
            s.fermi_dirac_, s.ionization_, s.doping_[i], s.donors_[i], s.acceptors_[i], s.n_ie_[i],
            s.log_dos_n_[i], s.log_dos_p_[i], s.levels_[i]);
        st.psi0_[i] = ohmic_contact_value(e, 0.0).psi - s.band_shift_[i];
    }

    s.composition_ = Composition::compose(s, st, models);
    s.m_ = s.composition_->dofs().block();  // 1

    // The pattern: the diagonal, both directions of every edge, and the terms' own entries.
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
    Terms::pattern(s, st, triplets);
    const auto size = static_cast<linalg::Index>(n);
    auto pattern = linalg::SparseMatrix::from_triplets(size, size, triplets);
    if (!pattern) return std::unexpected(std::move(pattern.error()));
    st.pattern_ = std::move(*pattern);

    s.block_.resize(n);
    for (std::size_t i = 0; i < n; ++i) s.block_[i] = detail::position(st.pattern_, i, i);
    s.edges_.reserve(scaled->edges.size());
    for (const detail::ScaledEdge& e : scaled->edges) {
        EdgeTerm t{};
        t.a = e.i;
        t.b = e.j;
        t.c = e.et * e.geometry;
        t.carriers = e.carriers;
        t.charged = e.charged;
        t.ab[0] = detail::position(st.pattern_, e.i, e.j);
        t.ba[0] = detail::position(st.pattern_, e.j, e.i);
        s.edges_.push_back(t);
    }
    return made;
}

std::expected<void, base::Error> SystemData::set_equilibrium_bias(
    SystemState& st, std::span<const double> bias_V) const {
    if (auto ok = check_contact_bias(kinds_, bias_V, true); !ok) return ok;
    st.gates_.set_bias(bias_V);
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (electrode_[i] < 0) continue;
        st.psi0_[i] = bias_V[static_cast<std::size_t>(electrode_[i])] / V_T_ + electrode_potential_[i];
    }
    return {};
}

// --- the terms -----------------------------------------------------------------------------------

// An insulator node off the electrodes: no charge, its row 0 and its diagonal 0 (assigned).
class SystemData::Terms::EquilibriumInsulatorRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"equilibrium insulator rows", TermDomain::nodes, psi_only, psi_only, true};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (dirichlet(s, i) || s.insulator_[i] == 0) continue;
            ev.f[i] = 0.0;
            if (ev.jacobian()) at(s, ev, i, 0, 0) = 0.0;
        }
    }
};

// A semiconductor node off the contacts: the space charge -V (n - p - C) with the carriers slaved
// to psi (Boltzmann or Fermi-Dirac), with incomplete ionization the ionized doping, and its psi
// derivative on the diagonal (assigned).
class SystemData::Terms::EquilibriumCharge final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"equilibrium charge", TermDomain::nodes, psi_only, psi_only, true};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (dirichlet(s, i) || s.insulator_[i] != 0) continue;
            const double eta = ev.x[i] + s.band_shift_[i];
            const physics::DensityResult carriers_n =
                detail::density(s.fermi_dirac_, s.n_ie_[i], s.log_dos_n_[i], eta);
            const physics::DensityResult carriers_p =
                detail::density(s.fermi_dirac_, s.n_ie_[i], s.log_dos_p_[i], -eta);
            // d/dpsi of the charge term, onto the diagonal (dn/dpsi = n.d_eta, dp/dpsi = -p.d_eta;
            // Boltzmann: n and -p; the ionized doping's: d_eta_c - d_eta_v).
            if (s.ionization_) {
                const detail::IonizedCharge c =
                    detail::ionized_charge(s.donors_[i], s.acceptors_[i], eta - s.log_dos_n_[i],
                                           -eta - s.log_dos_p_[i], s.levels_[i]);
                ev.f[i] = -s.volume_[i] * (carriers_n.density - carriers_p.density - c.value);
                if (ev.jacobian()) {
                    at(s, ev, i, 0, 0) =
                        -s.volume_[i] *
                        (carriers_n.d_eta + carriers_p.d_eta - c.d_eta_c + c.d_eta_v);
                }
            } else {
                ev.f[i] = -s.volume_[i] * (carriers_n.density - carriers_p.density - s.doping_[i]);
                if (ev.jacobian()) {
                    at(s, ev, i, 0, 0) = -s.volume_[i] * (carriers_n.d_eta + carriers_p.d_eta);
                }
            }
        }
    }
};

// The interface edges' half-edge fluxes with the traps at Fermi occupancy
// (InterfaceEdges::equilibrium), into the rows that are not Dirichlet rows.
class SystemData::Terms::EquilibriumInterfaces final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"equilibrium interfaces", TermDomain::interface_edges, psi_only, psi_only, false};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
        const std::span<double> values = ev.values;
        for (const InterfaceEdge& f : s.interfaces_.edges()) {
            const std::size_t a = f.insulator, b = f.semiconductor;
            const InterfaceEquilibrium q =
                s.interfaces_.equilibrium(f, ev.x[a], ev.x[b], s.statistics(b));
            // The insulator row a and the semiconductor row b, columns (a, b).
            const EdgeTerm& e = s.edges_[f.edge];
            const std::size_t ab = e.a == a ? e.ab[0] : e.ba[0];
            const std::size_t ba = e.a == a ? e.ba[0] : e.ab[0];
            if (!dirichlet(s, a)) {
                ev.f[a] += q.flux_insulator;
                if (ev.jacobian()) {
                    at(s, ev, a, 0, 0) += q.d_flux_insulator[0];
                    values[ab] += q.d_flux_insulator[1];
                }
            }
            if (!dirichlet(s, b)) {
                ev.f[b] += q.flux_semiconductor;
                if (ev.jacobian()) {
                    values[ba] += q.d_flux_semiconductor[0];
                    at(s, ev, b, 0, 0) += q.d_flux_semiconductor[1];
                }
            }
        }
    }
};

// The Dirichlet rows, last: psi - psi0 on an ohmic or Schottky contact node and on an electrode
// node, the row's Jacobian values cleared and its diagonal 1 (assigned; the shared Poisson flux
// adds into a Schottky node's row, which this replaces).
class SystemData::Terms::EquilibriumContacts final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"equilibrium contact rows", TermDomain::nodes, psi_only, psi_only, true};
    }
    void add(const SystemData& s, const SystemState& st, const Evaluation& ev) const override {
        const std::span<const linalg::Index> rows = st.pattern_.row_offsets();
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!dirichlet(s, i)) continue;
            ev.f[i] = ev.x[i] - st.psi0_[i];
            if (!ev.jacobian()) continue;
            std::fill(ev.values.begin() + static_cast<std::ptrdiff_t>(rows[i]),
                      ev.values.begin() + static_cast<std::ptrdiff_t>(rows[i + 1]), 0.0);
            at(s, ev, i, 0, 0) = 1.0;
        }
    }
};

std::span<const CompositionEntry> SystemData::Terms::equilibrium_table() {
    using E = EquationSet;
    using M = const PhysicsModels&;
    const auto always = [](const SystemData&, const SystemState&, M) { return true; };
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
    const auto contacts = [](const SystemData& s, const SystemState&, M) {
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (dirichlet(s, i)) return true;
        }
        return false;
    };
    static const CompositionEntry table[] = {
        {+E::oxide, insulators, &make_term<EquilibriumInsulatorRows>},
        {+E::poisson, always, &make_term<EquilibriumCharge>},
        {+E::gate, gated, &gate},
        {+E::interfaces, interfaces, &make_term<EquilibriumInterfaces>},
        {+E::poisson, always, &poisson_flux},
        {+E::contacts, contacts, &make_term<EquilibriumContacts>},
    };
    return table;
}

// --- the queries ---------------------------------------------------------------------------------

void SystemData::equilibrium_carriers(std::span<const double> psi, std::span<double> n,
                                      std::span<double> p) const {
    NITCAD_EXPECTS(equilibrium_);
    NITCAD_EXPECTS(psi.size() == node_count() && n.size() == node_count() &&
                   p.size() == node_count());
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) {
            n[i] = p[i] = 0.0;
            continue;
        }
        const double eta = psi[i] + band_shift_[i];
        n[i] = detail::density(fermi_dirac_, n_ie_[i], log_dos_n_[i], eta).density;
        p[i] = detail::density(fermi_dirac_, n_ie_[i], log_dos_p_[i], -eta).density;
    }
}

std::vector<double> SystemData::charge_neutral_potential(const SystemState& st) const {
    NITCAD_EXPECTS(equilibrium_);
    std::vector<double> psi(node_count());
    for (std::size_t i = 0; i < psi.size(); ++i) {
        psi[i] = Terms::dirichlet(*this, i) ? st.psi0_[i]
                 : insulator_[i] != 0       ? 0.0
                     : detail::neutral_equilibrium(fermi_dirac_, ionization_, doping_[i],
                                                   donors_[i], acceptors_[i], n_ie_[i],
                                                   log_dos_n_[i], log_dos_p_[i], levels_[i])
                               .eta -
                           band_shift_[i];
    }
    return psi;
}

BandEdges SystemData::equilibrium_band_edges(std::span<const double> psi) const {
    NITCAD_EXPECTS(psi.size() == node_count());
    const std::size_t n = node_count();
    std::vector<double> carriers_n(n), carriers_p(n);
    equilibrium_carriers(psi, carriers_n, carriers_p);
    BandEdges b{std::vector<double>(n), std::vector<double>(n), std::vector<double>(n),
                std::vector<double>(n)};
    for (std::size_t i = 0; i < n; ++i) {
        if (insulator_[i] != 0) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            b.conduction[i] = b.valence[i] = b.electron_fermi[i] = b.hole_fermi[i] = nan;
            continue;
        }
        const detail::NodeBands e =
            detail::node_bands(fermi_dirac_, n_ie_[i], log_dos_n_[i], log_dos_p_[i],
                               psi[i] + band_shift_[i], carriers_n[i], carriers_p[i]);
        b.conduction[i] = e.conduction;
        b.valence[i] = e.valence;
        b.electron_fermi[i] = e.electron_fermi;
        b.hole_fermi[i] = e.hole_fermi;
    }
    return b;
}

std::vector<double> SystemData::equilibrium_gate_charges(const SystemState& st,
                                                         std::span<const double> psi) const {
    NITCAD_EXPECTS(psi.size() == node_count());
    std::vector<double> charge = st.gates_.charges(psi, 1, kinds_.size());
    // An electrode's charge is the displacement flux leaving its nodes (Gauss over their boxes,
    // which hold no other charge).
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ci = electrode_[e.a], cj = electrode_[e.b];
        if (ci == cj || e.charged) continue;
        const double flux = e.c * (psi[e.a] - psi[e.b]);  // from a to b
        if (ci >= 0) charge[static_cast<std::size_t>(ci)] += flux;
        if (cj >= 0) charge[static_cast<std::size_t>(cj)] -= flux;
    }
    // An electrode node on an interface edge: the flux leaving it is minus its half-edge flux.
    for (const InterfaceEdge& f : interfaces_.edges()) {
        const std::int32_t c = electrode_[f.insulator];
        if (c < 0) continue;
        charge[static_cast<std::size_t>(c)] -=
            interfaces_
                .equilibrium(f, psi[f.insulator], psi[f.semiconductor],
                             statistics(f.semiconductor))
                .flux_insulator;
    }
    return charge;
}

std::vector<double> SystemData::equilibrium_trap_charges(std::span<const double> psi) const {
    NITCAD_EXPECTS(psi.size() == node_count());
    std::vector<double> q(interfaces_.interface_count(), 0.0);
    for (const InterfaceEdge& e : interfaces_.edges()) {
        q[e.interface] +=
            interfaces_
                .equilibrium(e, psi[e.insulator], psi[e.semiconductor], statistics(e.semiconductor))
                .trapped;
    }
    return q;
}

// --- equilibrium.hpp -----------------------------------------------------------------------------

std::expected<System, base::Error> create_equilibrium(const device::Device& device,
                                                      const Scaling& scaling,
                                                      const PhysicsModels& models) {
    return System::create(device, scaling, models, Formulation::equilibrium);
}

void equilibrium_carriers(const System& system, std::span<const double> psi, std::span<double> n,
                          std::span<double> p) {
    SystemData::of(system).equilibrium_carriers(psi, n, p);
}

std::vector<double> charge_neutral_potential(const System& system) {
    return SystemData::of(system).charge_neutral_potential(SystemData::state_of(system));
}

std::span<const double> contact_potential(const System& system) {
    NITCAD_EXPECTS(SystemData::of(system).equilibrium_);
    return SystemData::state_of(system).psi0_;
}

bool is_contact(const System& system, std::size_t node) {
    const SystemData& s = SystemData::of(system);
    NITCAD_EXPECTS(s.equilibrium_ && node < s.node_count());
    return SystemData::Terms::dirichlet(s, node);
}

}  // namespace NiTCAD::assemble
