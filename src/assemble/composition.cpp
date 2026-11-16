#include "composition.hpp"

#include <optional>
#include <utility>

namespace NiTCAD::assemble {

std::shared_ptr<const SystemData::Composition> SystemData::Composition::compose(
    const SystemData& s, const SystemState& st, const PhysicsModels& models) {
    auto c = std::make_shared<Composition>();
    // The fields the equations solve for: psi, n, p at every node, and with the heat equation
    // the temperature rise last (node-interleaved; insulator nodes keep their n and p, held at
    // 0, DECISIONS.md FD5).
    const Field fields[] = {potential_field, electron_field, hole_field, temperature_field};
    // With energy balance (Unit 52) the carrier temperature rises instead, tau_n then tau_p.
    const Field energy[] = {potential_field, electron_field, hole_field,
                            electron_temperature_field, hole_temperature_field};
    // Thermal equilibrium (Unit 57 GE5): psi alone, its own table whatever the models.
    c->dofs_ = s.equilibrium_ ? DofMap::node_interleaved(std::span(fields).first(1), s.node_count())
               : models.energy_balance
                   ? DofMap::node_interleaved(energy, s.node_count())
                   : DofMap::node_interleaved(std::span(fields).first(models.electrothermal ? 4 : 3),
                                              s.node_count());
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

}  // namespace NiTCAD::assemble
