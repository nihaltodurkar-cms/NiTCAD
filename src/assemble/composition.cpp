#include "composition.hpp"

#include <utility>

namespace NiTCAD::assemble {

std::shared_ptr<const DriftDiffusion::Composition> DriftDiffusion::Composition::compose(
    const DriftDiffusion& s, const PhysicsModels& models) {
    auto c = std::make_shared<Composition>();
    // The fields the equations solve for: psi, n, p at every node, and with the heat equation
    // the temperature rise last (node-interleaved; insulator nodes keep their n and p, held at
    // 0, DECISIONS.md FD5).
    const Field fields[] = {potential_field, electron_field, hole_field, temperature_field};
    c->dofs_ = DofMap::node_interleaved(std::span(fields).first(models.electrothermal ? 4 : 3),
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
    for (const std::shared_ptr<const Term>& t : terms_) {
        if (!interfaces && t->info().domain == TermDomain::interface_edges) continue;
        t->add(s, e);
    }
}

}  // namespace NiTCAD::assemble
