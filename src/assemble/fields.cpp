#include "NiTCAD/assemble/fields.hpp"

#include <algorithm>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::assemble {

DofMap DofMap::node_interleaved(std::span<const Field> fields, std::size_t nodes) {
    DofMap m;
    for (const Field& f : fields) {
        NITCAD_EXPECTS(f.location == FieldLocation::node && f.unknown);
        NITCAD_EXPECTS(std::ranges::none_of(m.fields_, [&](const Field& g) { return g.id == f.id; }));
        m.fields_.push_back(f);
    }
    m.nodes_ = nodes;
    return m;
}

std::optional<std::size_t> DofMap::slot(FieldId id) const noexcept {
    for (std::size_t k = 0; k < fields_.size(); ++k) {
        if (fields_[k].id == id) return k;
    }
    return std::nullopt;
}

std::size_t DofMap::index(FieldId id, std::size_t node) const {
    const std::optional<std::size_t> k = slot(id);
    NITCAD_EXPECTS(k.has_value() && node < nodes_);
    return block() * node + *k;
}

const Field& DofMap::field_of(std::size_t unknown) const {
    NITCAD_EXPECTS(unknown < unknowns());
    return fields_[unknown % block()];
}

std::size_t DofMap::node_of(std::size_t unknown) const {
    NITCAD_EXPECTS(unknown < unknowns());
    return unknown / block();
}

}  // namespace NiTCAD::assemble
