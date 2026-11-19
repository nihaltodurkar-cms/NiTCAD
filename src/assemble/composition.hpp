// Private to the assemble layer: the composition of a drift-diffusion system from equation terms
// (Unit 40 slice 3; DECISIONS.md FE4; ARCHITECTURE.md 6.12.2).
//
// A system is a set of equations activated from the device and the models, not a separate solver:
//     diode   = Poisson + electron continuity + hole continuity (+ the contacts' conditions)
//     MOS     = the diode's + the oxide's Poisson + gate and interface conditions
//     thermal = the semiconductor equations at T + the heat equation and its sinks
// with generation (impact ionization, band-to-band tunnelling) and recombination models adding
// terms; thermal equilibrium (Unit 57 GE5) = Poisson alone over psi, the carriers slaved to it, with
// the oxide, gate, interface and contact conditions. A composition table lists, in assembly order (DECISIONS.md FE3), every term a system can
// have: the equation sets it belongs to, whether it is active for a device and models, and its
// factory. The isothermal table is in terms.cpp, the electrothermal one in thermal_terms.cpp, the
// equilibrium one in equilibrium.cpp.
// Composition picks the table, keeps the active rows and builds the unknowns' map from the
// fields the equations solve for (fields.hpp).
//
// The model data the terms read (per node, per edge, per contact) is SystemData's and what changes
// after creation SystemState's (system_data.hpp, Unit 57 GE2); System is the compatibility
// wrapper (FE4) holding them.
#pragma once

#include <memory>
#include <span>
#include <vector>

#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/system.hpp"
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
    quantum = 1u << 9,     // the quantum potentials of density gradient (Unit 62)
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
    bool (*active)(const SystemData& s, const SystemState& st, const PhysicsModels& models);
    std::shared_ptr<const SystemData::Term> (*make)();
    // A term of the frequency domain alone (Term::add_frequency; Unit 57 GE3).
    bool frequency = false;
};

// An equation added to the drift-diffusion composition from outside its tables (Unit 57 GG8,
// SystemData::create_extended): its fields, node unknowns placed after the built-in fields in
// every node block, and its terms, added after every term of the table, in order, under no
// equation set. Its terms write their own rows and any others: the rows of its fields are theirs
// to assign (no table term writes them), and the Jacobian entries beyond every node's full block
// are theirs to add (Term::pattern), found by position at each add (they hold no positions in the
// system's data). The Newton hooks move each added unknown as the potential: clipped to
// +-max_update, measured by |dx| (System::update_size, apply_update). Steady evaluations only (no
// time step, no small-signal matrix).
struct SystemData::CompositionExtension {
    std::vector<Field> fields;
    std::vector<std::shared_ptr<const Term>> terms;
};

class SystemData::Composition {
public:
    // The unknowns and the active terms of a system's device data and models, and an extension's
    // (nullptr: none).
    [[nodiscard]] static std::shared_ptr<const Composition> compose(
        const SystemData& s, const SystemState& st, const PhysicsModels& models,
        const CompositionExtension* extension);

    [[nodiscard]] const DofMap& dofs() const noexcept { return dofs_; }
    [[nodiscard]] std::span<const std::shared_ptr<const Term>> terms() const noexcept {
        return terms_;
    }
    // The equation sets of term k, and of the whole system (their union).
    [[nodiscard]] EquationSets sets(std::size_t k) const { return sets_[k]; }
    [[nodiscard]] EquationSets sets() const noexcept { return all_; }

    // Appends every term's own pattern entries (Term::pattern).
    void pattern(const SystemData& s, const SystemState& st,
                 std::vector<linalg::Triplet>& triplets) const;
    // Every steady term's add, in order. In the frequency domain (small_signal_matrix), the steady
    // terms of a domain that a frequency-domain term has are left out: the frequency-domain term
    // takes their place (assemble_frequency).
    void assemble(const SystemData& s, const SystemState& st, const Terms::Evaluation& e,
                  bool frequency) const;
    // Every frequency-domain term's add_frequency, in order.
    void assemble_frequency(const SystemData& s, const SystemState& st,
                            const Terms::FrequencyEvaluation& e) const;
    [[nodiscard]] std::span<const std::shared_ptr<const Term>> frequency_terms() const noexcept {
        return frequency_terms_;
    }

private:
    DofMap dofs_;
    std::vector<std::shared_ptr<const Term>> terms_;
    std::vector<std::shared_ptr<const Term>> frequency_terms_;
    // The domains whose steady terms the frequency-domain terms replace (a bit per TermDomain).
    unsigned frequency_domains_ = 0;
    std::vector<EquationSets> sets_;
    EquationSets all_ = 0;
};

}  // namespace NiTCAD::assemble
