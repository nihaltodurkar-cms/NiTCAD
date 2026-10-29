#include "composition.hpp"

#include <optional>
#include <utility>

namespace NiTCAD::assemble {

std::shared_ptr<const DriftDiffusion::Composition> DriftDiffusion::Composition::compose(
    const DriftDiffusion& s, const PhysicsModels& models) {
    auto c = std::make_shared<Composition>();
    // The fields the equations solve for: psi, n, p at every node, and with the heat equation
    // the temperature rise last (node-interleaved; insulator nodes keep their n and p, held at
    // 0, DECISIONS.md FD5).
    const Field fields[] = {potential_field, electron_field, hole_field, temperature_field};
    // With energy balance (Unit 52) the carrier temperature rises instead, tau_n then tau_p.
    const Field energy[] = {potential_field, electron_field, hole_field,
                            electron_temperature_field, hole_temperature_field};
    c->dofs_ = models.energy_balance
                   ? DofMap::node_interleaved(energy, s.node_count())
                   : DofMap::node_interleaved(std::span(fields).first(models.electrothermal ? 4 : 3),
                                              s.node_count());
    for (const CompositionEntry& e :
         models.electrothermal ? Terms::thermal_table() : Terms::isothermal_table()) {
        if (!e.active(s, models)) continue;
        c->terms_.push_back(e.make());
        c->sets_.push_back(e.sets);
        c->all_ |= e.sets;
    }
    return c;
}

void DriftDiffusion::Composition::pattern(const DriftDiffusion& s,
                                          std::vector<linalg::Triplet>& triplets) const {
    for (const std::shared_ptr<const Term>& t : terms_) t->pattern(s, triplets);
}

void DriftDiffusion::Composition::assemble(const DriftDiffusion& s, const Terms::Evaluation& e,
                                           bool interfaces) const {
    // Gate tunnelling's fluxes start on worker threads now and the other terms add theirs
    // meanwhile; the gate tunnelling term collects them (GR11).
    std::optional<GateTunnelJob> tunnels;
    Terms::Evaluation ev = e;
    if (!s.gate_tunnels_.empty() && ev.tunnels == nullptr) {
        tunnels.emplace(s, e.x, e.g);
        ev.tunnels = &*tunnels;
    }
    for (const std::shared_ptr<const Term>& t : terms_) {
        if (!interfaces && t->info().domain == TermDomain::interface_edges) continue;
        t->add(s, ev);
    }
}

}  // namespace NiTCAD::assemble
