// Derived fields (Unit 57 slice 3, DECISIONS.md GE4, FD7; ARCHITECTURE.md 6.12.1).
//
// The quantities a system computes from a state rather than solves for, as one table: each row
// is a diagnostic, the fields it gives (kind derived, fields.hpp) and the function that evaluates
// them at a state x. A row's fields are computed together because they come from one pass (the
// four band levels from one band_edges, every mechanism's rate from one rates). The evaluators
// call the topic functions (band_edges.hpp, rates.hpp, currents.hpp, charges.hpp), so the values
// are theirs operand for operand; the solve layer fills results::BiasPoint, RateFields,
// BandDiagram and TransientSnapshot from them.
#pragma once

#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

#include "NiTCAD/assemble/fields.hpp"

namespace NiTCAD::assemble {

class System;

enum class DiagnosticId {
    bands,             // E_c, E_v, E_Fn, E_Fp per node, in V_T (NaN on insulator nodes)
    rates,             // per mechanism and net, per node, in cm^-3 s^-1
    edge_currents,     // Jn, Jp per mesh edge, scaled as currents.hpp edge_currents
    edge_flux_terms,   // per mesh edge, as currents.hpp edge_flux_terms
    trap_occupancies,  // f per trap slot, the steady state
    temperatures,      // T, T_n, T_p per node in K
};

// The values of a diagnostic's fields at a state, one vector per field in the row's order. A
// vector is empty when the system lacks the field: a rate mechanism that is off, the lattice
// temperature without the electrothermal model, the carrier temperatures without energy
// balance, the occupancies without traps.
using DerivedValues = std::vector<std::vector<double>>;

struct Diagnostic {
    DiagnosticId id;
    std::string_view name;
    std::span<const Field> fields;
    // Precondition (NITCAD_EXPECTS, of the topic functions): x has system.unknowns() entries
    // and positive densities off the insulators.
    DerivedValues (*evaluate)(const System& system, std::span<const double> x);
};

// Every diagnostic, in DiagnosticId order.
[[nodiscard]] std::span<const Diagnostic> diagnostics() noexcept;
[[nodiscard]] const Diagnostic& diagnostic(DiagnosticId id) noexcept;

// The position of a field in a diagnostic's row. Precondition (NITCAD_EXPECTS): the row has it.
[[nodiscard]] std::size_t position(const Diagnostic& diagnostic, FieldId field);

// One field's values, moved out of a diagnostic's evaluated values (the other fields' are left
// as they were). Precondition (NITCAD_EXPECTS): the row has the field and `values` one vector
// per field of the row.
[[nodiscard]] std::vector<double> take(DerivedValues& values, const Diagnostic& diagnostic,
                                       FieldId field);

}  // namespace NiTCAD::assemble
