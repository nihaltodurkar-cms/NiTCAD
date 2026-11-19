// Fields and the unknown map (Unit 40 slice 1; ARCHITECTURE.md 6.12.1; DECISIONS.md FE1).
//
// A Field describes a named physical quantity of a system: its units, where it lives (nodes only
// so far), whether it is a Newton unknown, and on which nodes it means something. A DofMap places
// the unknowns of a system's fields: the unknown of (field, node) and the field and node of an
// unknown. The drift-diffusion system's map is node-interleaved, the layout every result, file and
// checkpoint already uses: node i's unknowns are block() i ... block() i + block() - 1, in the order
// of the fields (psi, n, p, then the temperature rise with the electrothermal model). Every field
// has an unknown at every node: on an insulator node the densities are unknowns held at 0
// (DECISIONS.md FD5).
#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace NiTCAD::assemble {

// The quantities a system may carry: its unknowns, then the derived quantities of derived.hpp
// (Unit 57 GE4), each computed from a state.
enum class FieldId {
    potential,    // psi: the electrostatic potential
    electrons,    // n: the electron density
    holes,        // p: the hole density
    temperature,  // tau = (T - T0) / T0: the lattice temperature rise (electrothermal, Unit 23)
    electron_temperature,  // tau_n = (T_n - T_L) / T_L (energy balance, Unit 52)
    hole_temperature,      // tau_p = (T_p - T_L) / T_L
    // Density gradient (Unit 62, DECISIONS.md QP2): the quantum potentials Lambda_n / kT and
    // Lambda_p / kT, after psi, n, p (stride 5); after psi alone at thermal equilibrium (stride 3).
    quantum_potential_n,
    quantum_potential_p,
    // Derived (derived.hpp).
    conduction_band,     // E_c, band_edges.hpp
    valence_band,        // E_v
    electron_fermi,      // E_Fn
    hole_fermi,          // E_Fp
    srh_rate,            // per mechanism, R - G (rates.hpp)
    auger_rate,
    radiative_rate,
    impact_rate,
    interface_rate_n,    // per carrier where the mechanism moves the carriers apart
    interface_rate_p,
    band_to_band_rate_n,
    band_to_band_rate_p,
    net_rate_n,          // the sum of every mechanism, as each carrier's rows see it
    net_rate_p,
    electron_edge_current,  // Jn, Jp of each mesh edge (currents.hpp edge_currents)
    hole_edge_current,
    edge_flux_terms,        // currents.hpp edge_flux_terms
    trap_occupancy,         // f of each trap slot (charges.hpp trap_occupancies)
    lattice_temperature,     // T in K (electrothermal)
    electron_temperature_K,  // T_n, T_p in K (energy balance)
    hole_temperature_K,
};

// Where a field's values are: one per mesh node, per mesh edge, or per trap slot of the interface
// edges (each interface edge's trap levels, in the assembler's slot order).
enum class FieldLocation { node, edge, interface_edge };

// A Newton unknown; a value fixed by the device (none listed yet); or derived from a state.
enum class FieldKind { unknown, fixed, derived };

// The nodes on which a field is physical: every node, or the semiconductor nodes (elsewhere its
// unknown is held at 0, a derived node field NaN or 0 as its function states).
enum class FieldDomain { all, semiconductor };

struct Field {
    FieldId id;
    std::string_view name;
    std::string_view units;  // of the scaled unknown, or of the derived values as evaluated
    FieldLocation location;
    FieldKind kind;
    FieldDomain domain;
};

// The drift-diffusion system's fields (scaled by de Mari, 6.1).
inline constexpr Field potential_field{FieldId::potential, "potential", "V_T", FieldLocation::node,
                                       FieldKind::unknown, FieldDomain::all};
inline constexpr Field electron_field{FieldId::electrons, "electrons", "Ns", FieldLocation::node,
                                      FieldKind::unknown, FieldDomain::semiconductor};
inline constexpr Field hole_field{FieldId::holes, "holes", "Ns", FieldLocation::node,
                                  FieldKind::unknown, FieldDomain::semiconductor};
inline constexpr Field temperature_field{FieldId::temperature, "temperature rise", "T0",
                                         FieldLocation::node, FieldKind::unknown,
                                         FieldDomain::all};
inline constexpr Field electron_temperature_field{FieldId::electron_temperature,
                                                  "electron temperature rise", "T_L",
                                                  FieldLocation::node, FieldKind::unknown,
                                                  FieldDomain::semiconductor};
inline constexpr Field hole_temperature_field{FieldId::hole_temperature, "hole temperature rise",
                                              "T_L", FieldLocation::node, FieldKind::unknown,
                                              FieldDomain::semiconductor};
inline constexpr Field quantum_potential_n_field{FieldId::quantum_potential_n,
                                                 "electron quantum potential", "kT",
                                                 FieldLocation::node, FieldKind::unknown,
                                                 FieldDomain::semiconductor};
inline constexpr Field quantum_potential_p_field{FieldId::quantum_potential_p,
                                                 "hole quantum potential", "kT",
                                                 FieldLocation::node, FieldKind::unknown,
                                                 FieldDomain::semiconductor};

class DofMap {
public:
    DofMap() = default;
    // Node-interleaved: every field an unknown at every node, in the given order. Precondition
    // (NITCAD_EXPECTS): the fields are node unknowns, none twice.
    [[nodiscard]] static DofMap node_interleaved(std::span<const Field> fields, std::size_t nodes);

    [[nodiscard]] std::span<const Field> fields() const noexcept { return fields_; }
    [[nodiscard]] std::size_t node_count() const noexcept { return nodes_; }
    // Unknowns per node.
    [[nodiscard]] std::size_t block() const noexcept { return fields_.size(); }
    [[nodiscard]] std::size_t unknowns() const noexcept { return block() * nodes_; }
    [[nodiscard]] bool has(FieldId id) const noexcept { return slot(id).has_value(); }
    // The position of a field in each node's block, if the system has it.
    [[nodiscard]] std::optional<std::size_t> slot(FieldId id) const noexcept;
    // The unknown of a field at a node. Preconditions (NITCAD_EXPECTS): the system has the field;
    // node < node_count().
    [[nodiscard]] std::size_t index(FieldId id, std::size_t node) const;
    // The field and the node of an unknown. Precondition (NITCAD_EXPECTS): unknown < unknowns().
    [[nodiscard]] const Field& field_of(std::size_t unknown) const;
    [[nodiscard]] std::size_t node_of(std::size_t unknown) const;

private:
    std::vector<Field> fields_;
    std::size_t nodes_ = 0;
};

}  // namespace NiTCAD::assemble
