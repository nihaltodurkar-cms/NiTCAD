// Private to the assemble layer: the equation terms of the drift-diffusion rows (Unit 40,
// DECISIONS.md FE2, FE3).
//
// A term is one piece of the residual F and its exact Jacobian: it names the rows (fields) it
// writes and the fields its partials are taken in, the domain it loops over, the pattern entries
// it needs beyond every node's full block and the edges' stencil, and an add() that adds its
// residual and partials at the Jacobian positions the system located at creation (block_,
// EdgeTerm::ab / ba, interface_np_, ImpactEdge, PathLayout). It reads the system's data at each
// call and holds none of its own, so a system stays copyable and a term is shared between copies.
//
// The order of the additions into every entry of F and J is the assembler's before Unit 40
// (DECISIONS.md, "The current accumulation order"): the node terms in the order ohmic contact
// rows, insulator rows, Poisson charge, incomplete ionization, gate, SRH, Auger, radiative,
// storage, each over the nodes in index order (a node term writes only its own node's rows, so one
// loop per term adds into every entry in the same order as one fused loop); then the interface
// edges, the mesh edges' Poisson flux and carrier fluxes (disjoint rows), impact ionization, local
// and nonlocal band-to-band tunnelling. The contact rows, the insulator rows, the Poisson charge
// and the ionization override assign their rows, as before; every other term adds. The terms keep
// the expressions of the assembler before them operand for operand.
//
// The electrothermal rows are not terms yet (assemble_heat; slice 2b).
#pragma once

#include <memory>
#include <span>
#include <string_view>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::assemble {

// What a term loops over, once per evaluation (DECISIONS.md FD3: one call per term per pass).
enum class TermDomain { nodes, interface_edges, mesh_edges, global };

struct TermInfo {
    std::string_view name;
    TermDomain domain;
    std::span<const FieldId> rows;     // the fields whose rows it writes
    std::span<const FieldId> columns;  // the fields its partials are taken in
    bool assigns;                      // it sets its rows rather than adding to them
};

struct DriftDiffusion::Terms {
    // One evaluation: the state, F (every row is assigned by a node term), the Jacobian's values
    // (empty: the residual alone; else zero on entry), the time step (nullptr: steady) and the
    // degeneracy factors (degeneracies(x)).
    struct Evaluation {
        std::span<const double> x;
        std::span<double> f;
        std::span<double> values;
        const TimeStep* step;
        std::span<const NodeDegeneracy> g;
        [[nodiscard]] bool jacobian() const noexcept { return !values.empty(); }
    };

    // The terms of a system's models in assembly order. Inclusion follows the models alone, since
    // the system's tables (impact nodes, tunnel nodes) are built after the pattern the terms add
    // to; a term whose table is empty adds nothing.
    [[nodiscard]] static std::vector<std::shared_ptr<const Term>> make(const DriftDiffusion& s,
                                                                       const PhysicsModels& models);
    // Appends every term's own pattern entries (Term::pattern) to triplets.
    static void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& triplets);
    // A system's terms (for its tests).
    [[nodiscard]] static std::span<const std::shared_ptr<const Term>> list(
        const DriftDiffusion& s) noexcept {
        return s.terms_;
    }
    // One term alone at x: F and the values set to 0, then the term's add (for its tests: each
    // term's partials against finite differences of its own residual). Preconditions
    // (NITCAD_EXPECTS): as evaluate's, step's spans of evaluate's sizes when not nullptr.
    static void evaluate(const DriftDiffusion& s, const Term& t, std::span<const double> x,
                         const TimeStep* step, std::span<double> f,
                         linalg::SparseMatrix& jacobian);

    // Whether node i takes the semiconductor node terms: neither an ohmic contact node nor an
    // insulator node.
    [[nodiscard]] static bool semiconductor(const DriftDiffusion& s, std::size_t i) noexcept {
        return s.contact_[i] < 0 && s.insulator_[i] == 0;
    }
    // The Jacobian value of (row r, column c) of node i's block.
    [[nodiscard]] static double& at(const DriftDiffusion& s, const Evaluation& e, std::size_t node,
                                    std::size_t r, std::size_t c) noexcept {
        return e.values[s.block_[s.m_ * s.m_ * node + s.m_ * r + c]];
    }
    // The equilibrium product of node i, scaled.
    [[nodiscard]] static physics::EquilibriumProduct product(const DriftDiffusion& s,
                                                             const Evaluation& e, std::size_t i);
    // Auger and radiative recombination take physical densities: node i's physical equilibrium
    // product, and the addition of a physical rate R into its continuity rows.
    [[nodiscard]] static physics::EquilibriumProduct physical_product(const DriftDiffusion& s,
                                                                      const Evaluation& e,
                                                                      std::size_t i);
    static void add_physical(const DriftDiffusion& s, const Evaluation& e, std::size_t i,
                             const physics::RecombinationRate& r);

    class ContactRows;
    class InsulatorRows;
    class PoissonCharge;
    class Ionization;
    class Gate;
    class Srh;
    class Auger;
    class Radiative;
    class Storage;
    class Interfaces;
    class PoissonFlux;
    class CarrierFluxes;
    class ImpactGeneration;
    class LocalTunnelling;
    class PathTunnelling;
};

class DriftDiffusion::Term {
public:
    Term() = default;
    Term(const Term&) = delete;
    Term& operator=(const Term&) = delete;
    virtual ~Term() = default;

    [[nodiscard]] virtual TermInfo info() const noexcept = 0;
    // Appends the Jacobian entries the term writes beyond every node's full block and the mesh
    // edges' stencil (an interface's cross-node density columns, impact ionization's other
    // carrier, the tunnel paths'); none by default.
    virtual void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& triplets) const;
    virtual void add(const DriftDiffusion& s, const Terms::Evaluation& e) const = 0;
};

}  // namespace NiTCAD::assemble
