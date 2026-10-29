// Private to the assemble layer: the composition of a drift-diffusion system from equation terms
// (Unit 40 slice 3; DECISIONS.md FE4; ARCHITECTURE.md 6.12.2).
//
// A system is a set of equations activated from the device and the models, not a separate solver:
//     diode   = Poisson + electron continuity + hole continuity (+ the contacts' conditions)
//     MOS     = the diode's + the oxide's Poisson + gate and interface conditions
//     thermal = the semiconductor equations at T + the heat equation and its sinks
// with generation (impact ionization, band-to-band tunnelling) and recombination models adding
// terms. A composition table lists, in assembly order (DECISIONS.md FE3), every term a system can
// have: the equation sets it belongs to, whether it is active for a device and models, and its
// factory. The isothermal table is in terms.cpp, the electrothermal one in thermal_terms.cpp.
// Composition picks the table, keeps the active rows and builds the unknowns' map from the
// fields the equations solve for (fields.hpp).
//
// The model data the terms read (per node, per edge, per contact) stays in DriftDiffusion: its
// diagnostics read it too (FE5). DriftDiffusion is the compatibility wrapper (FE4) and that data's
// store.
#pragma once

#include <memory>
#include <span>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "terms.hpp"

namespace NiTCAD::assemble {

// The equation sets of ARCHITECTURE.md 6.12.2 a term belongs to (a bit set: the electrothermal
// edge term is the Poisson, continuity and heat equations' at once).
enum class EquationSet : unsigned {
    poisson = 1u << 0,     // the electrostatic potential's equation, semiconductor and oxide
    continuity = 1u << 1,  // electron and hole continuity: transport, recombination, storage
    contacts = 1u << 2,    // the ohmic contacts' Dirichlet rows; the Schottky contacts' rows
    oxide = 1u << 3,       // insulator nodes: no carriers, an electrode's potential
    gate = 1u << 4,        // lumped gates' oxide coupling
    interfaces = 1u << 5,  // declared interfaces: charge, traps, surface recombination
    generation = 1u << 6,  // impact ionization, band-to-band tunnelling
    heat = 1u << 7,        // the lattice heat equation and its thermal contacts
    energy = 1u << 8,      // the carriers' energy-balance equations (Unit 52)
};
using EquationSets = unsigned;
// One set as a bit set.
[[nodiscard]] constexpr EquationSets operator+(EquationSet s) noexcept {
    return static_cast<EquationSets>(s);
}
[[nodiscard]] constexpr EquationSets operator|(EquationSet a, EquationSet b) noexcept {
    return static_cast<EquationSets>(a) | static_cast<EquationSets>(b);
}
[[nodiscard]] constexpr EquationSets operator|(EquationSets a, EquationSet b) noexcept {
    return a | static_cast<EquationSets>(b);
}
[[nodiscard]] constexpr bool contains(EquationSets sets, EquationSet s) noexcept {
    return (sets & static_cast<EquationSets>(s)) != 0;
}

// One row of a composition table.
struct CompositionEntry {
    EquationSets sets;
    // Whether the term is in a system with these models; reads the device data create has set
    // before the unknowns' map (contacts, insulators, gates, interfaces).
    bool (*active)(const DriftDiffusion& s, const PhysicsModels& models);
    std::shared_ptr<const TermAccess::Term> (*make)();
};

class DriftDiffusion::Composition {
public:
    // The unknowns and the active terms of a system's device data and models.
    [[nodiscard]] static std::shared_ptr<const Composition> compose(const DriftDiffusion& s,
                                                                    const PhysicsModels& models);

    [[nodiscard]] const DofMap& dofs() const noexcept { return dofs_; }
    [[nodiscard]] std::span<const std::shared_ptr<const Term>> terms() const noexcept {
        return terms_;
    }
    // The equation sets of term k, and of the whole system (their union).
    [[nodiscard]] EquationSets sets(std::size_t k) const { return sets_[k]; }
    [[nodiscard]] EquationSets sets() const noexcept { return all_; }

    // Appends every term's own pattern entries (Term::pattern).
    void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& triplets) const;
    // Every term's add, in order; without `interfaces`, the interface-edge terms are left out
    // (the isothermal small_signal_matrix adds its own). assemble_heat always passes true: the
    // electrothermal small_signal_matrix adds none.
    void assemble(const DriftDiffusion& s, const Terms::Evaluation& e, bool interfaces) const;

private:
    DofMap dofs_;
    std::vector<std::shared_ptr<const Term>> terms_;
    std::vector<EquationSets> sets_;
    EquationSets all_ = 0;
};

}  // namespace NiTCAD::assemble
