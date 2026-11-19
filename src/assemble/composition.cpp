#include "composition.hpp"

#include <expected>
#include <optional>
#include <utility>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::assemble {

std::shared_ptr<const SystemData::Composition> SystemData::Composition::compose(
    const SystemData& s, const SystemState& st, const PhysicsModels& models,
    const CompositionExtension* extension) {
    auto c = std::make_shared<Composition>();
    // The fields the equations solve for: psi, n, p at every node, and with the heat equation
    // the temperature rise last (node-interleaved; insulator nodes keep their n and p, held at
    // 0, DECISIONS.md FD5).
    const Field fields[] = {potential_field, electron_field, hole_field, temperature_field};
    // With energy balance (Unit 52) the carrier temperature rises instead, tau_n then tau_p.
    const Field energy[] = {potential_field, electron_field, hole_field,
                            electron_temperature_field, hole_temperature_field};
    // Thermal equilibrium (Unit 57 GE5): psi alone, its own table whatever the models.
    // With density gradient (Unit 62, QP2) the quantum potentials after them: psi, n, p, Lambda_n,
    // Lambda_p, and at equilibrium psi, Lambda_n, Lambda_p.
    const Field quantum[] = {potential_field, electron_field, hole_field, quantum_potential_n_field,
                             quantum_potential_p_field};
    const Field quantum_equilibrium[] = {potential_field, quantum_potential_n_field,
                                         quantum_potential_p_field};
    const std::span<const Field> built_in =
        s.density_gradient_     ? (s.equilibrium_ ? std::span<const Field>(quantum_equilibrium)
                                                  : std::span<const Field>(quantum))
        : s.equilibrium_        ? std::span(fields).first(1)
        : models.energy_balance ? std::span<const Field>(energy)
                                : std::span(fields).first(models.electrothermal ? 4 : 3);
    if (extension == nullptr) {
        c->dofs_ = DofMap::node_interleaved(built_in, s.node_count());
    } else {  // Unit 57 GG8: the extension's fields after the built-in ones
        NITCAD_EXPECTS(!s.equilibrium_);
        std::vector<Field> all(built_in.begin(), built_in.end());
        all.insert(all.end(), extension->fields.begin(), extension->fields.end());
        c->dofs_ = DofMap::node_interleaved(all, s.node_count());
    }
    const std::span<const CompositionEntry> table =
        s.equilibrium_          ? Terms::equilibrium_table()
        : models.electrothermal ? Terms::thermal_table()
                                : Terms::isothermal_table();
    for (const CompositionEntry& e : table) {
        if (!e.active(s, st, models)) continue;
        if (e.frequency) {  // Unit 57 GE3: not in the steady list, its own sets already are
            c->frequency_terms_.push_back(e.make());
            const auto domain = static_cast<unsigned>(c->frequency_terms_.back()->info().domain);
            c->frequency_domains_ |= 1u << domain;
            continue;
        }
        c->terms_.push_back(e.make());
        c->sets_.push_back(e.sets);
        c->all_ |= e.sets;
    }
    if (extension != nullptr) {  // after the table's terms, under no equation set
        for (const std::shared_ptr<const Term>& t : extension->terms) {
            NITCAD_EXPECTS(t != nullptr);
            c->terms_.push_back(t);
            c->sets_.push_back(0);
        }
    }
    return c;
}

void SystemData::Composition::pattern(const SystemData& s, const SystemState& st,
                                          std::vector<linalg::Triplet>& triplets) const {
    for (const std::shared_ptr<const Term>& t : terms_) t->pattern(s, st, triplets);
}

void SystemData::Composition::assemble(const SystemData& s, const SystemState& st,
                                       const Terms::Evaluation& e,
                                       bool frequency) const {
    // Gate tunnelling's fluxes start on worker threads now and the other terms add theirs
    // meanwhile; the gate tunnelling term collects them (GR11).
    std::optional<GateTunnelJob> tunnels;
    Terms::Evaluation ev = e;
    if (!s.gate_tunnels_.empty() && ev.tunnels == nullptr) {
        tunnels.emplace(s, st, e.x, e.g);
        ev.tunnels = &*tunnels;
    }
    for (const std::shared_ptr<const Term>& t : terms_) {
        if (frequency &&
            (frequency_domains_ & (1u << static_cast<unsigned>(t->info().domain))) != 0) {
            continue;
        }
        t->add(s, st, ev);
    }
}

void SystemData::Composition::assemble_frequency(const SystemData& s, const SystemState& st,
                                                 const Terms::FrequencyEvaluation& e) const {
    for (const std::shared_ptr<const Term>& t : frequency_terms_) t->add_frequency(s, st, e);
}


std::expected<System, base::Error> SystemData::create_extended(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models,
    const CompositionExtension& extension) {
    NITCAD_EXPECTS(!models.electrothermal && !models.energy_balance && !extension.fields.empty());
    for (const Field& f : extension.fields) {
        NITCAD_EXPECTS(f.kind == FieldKind::unknown && f.location == FieldLocation::node);
    }
    // System's constructor is its own: a system created without the extension is the shell whose
    // data and state are replaced (System::replace_data, state; the passkey).
    auto out = System::create(device, scaling, models);
    if (!out) return out;
    auto made = make_drift_diffusion(device, scaling, models, &extension);
    if (!made) return std::unexpected(std::move(made.error()));
    const std::vector<double> zero(made->data->contact_count_, 0.0);
    replace(*out, std::move(made->data));
    state_of(*out) = std::move(made->state);
    if (auto ok = out->set_bias(zero); !ok) return std::unexpected(std::move(ok.error()));
    return out;
}

}  // namespace NiTCAD::assemble
