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
// With the electrothermal model the rows are the electrothermal terms (thermal_terms.cpp), in
// assemble_heat's order: per node the heat row (sinks, lattice storage) before the electrical
// rows, the carriers' storage last (it adds to the heat row too); then the interfaces, the edges
// (Poisson, carriers and heat in one term: the heat row takes the carrier fluxes), impact
// ionization.
#pragma once

#include <atomic>
#include <cstdint>
#include <memory>
#include <span>
#include <string_view>
#include <thread>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::assemble {

struct CompositionEntry;  // composition.hpp

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
    // degeneracy factors (degeneracies(x); isothermal). With the electrothermal model
    // (assemble_heat): the nodes' thermal states and levels (thermal_states, thermal_levels) and
    // whether the thermal contacts' terms are in (sinks).
    struct Evaluation {
        std::span<const double> x;
        std::span<double> f;
        std::span<double> values;
        const TimeStep* step;
        std::span<const NodeDegeneracy> g;
        std::span<const ThermalState> thermal = {};
        std::span<const NodeLevels> levels = {};
        bool sinks = true;
        // The gate tunnels' fluxes at x under way (Composition::assemble); nullptr: the gate
        // tunnelling term computes them itself.
        GateTunnelJob* tunnels = nullptr;
        [[nodiscard]] bool jacobian() const noexcept { return !values.empty(); }
    };

    // The composition tables (composition.hpp): every term a system can have, in assembly order,
    // with its equation sets and activation. Activation follows the models and the device data set
    // before the unknowns' map, since the system's other tables (impact nodes, tunnel nodes) are
    // built after the pattern the terms add to; a term whose table is empty adds nothing.
    // Isothermal in terms.cpp; with the electrothermal model thermal_terms.cpp.
    [[nodiscard]] static std::span<const CompositionEntry> isothermal_table();
    [[nodiscard]] static std::span<const CompositionEntry> thermal_table();
    // The Poisson charge term, in both tables.
    [[nodiscard]] static std::shared_ptr<const Term> poisson_charge();
    // Appends every term's own pattern entries (Term::pattern) to triplets.
    static void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& triplets);
    // A system's terms and its composition (for its tests).
    [[nodiscard]] static std::span<const std::shared_ptr<const Term>> list(const DriftDiffusion& s);
    [[nodiscard]] static const Composition& composition(const DriftDiffusion& s) noexcept {
        return *s.composition_;
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
    // The pattern entries of the interface edges (the insulator's Poisson row reads the
    // semiconductor's n and p) and of impact ionization (a node's electron generation reads its
    // neighbours' p, its hole generation their n), shared by the isothermal and electrothermal
    // terms.
    static void interface_pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t);
    static void impact_pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t);
    // Lombardi (Unit 46): rows 1 .. rows of both ends of every surface edge with the psi
    // columns of its E_perp stencil (rows 2: the continuity rows; 3: and the heat row).
    static void surface_pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t,
                                std::size_t rows);

    // Electrothermal: whether node i's heat row is an isothermal sink's Dirichlet row.
    [[nodiscard]] static bool held(const DriftDiffusion& s, const Evaluation& e,
                                   std::size_t i) noexcept {
        const std::int32_t c = s.thermal_nodes_[i].contact;
        return e.sinks && c >= 0 && s.thermal_isothermal_[static_cast<std::size_t>(c)] != 0;
    }
    // Electrothermal: node i's equilibrium product at its T (scaled), n_ie^2 or, under
    // Fermi-Dirac statistics, n_ie^2 gamma_n gamma_p with its density partials, and its theta
    // derivative dE_t.
    static void thermal_product(const DriftDiffusion& s, const Evaluation& e, std::size_t i,
                                physics::EquilibriumProduct& E, double& dE_t);
    // Electrothermal: k times a recombination rate R / R0 into node i's continuity rows, with its
    // density partials and d_t its theta derivative.
    static void add_thermal_rate(const DriftDiffusion& s, const Evaluation& e, std::size_t i,
                                 double k, double rate, double d_n, double d_p, double d_t);

    class ContactRows;
    class SchottkyRows;
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
    class GateTunnelling;
    // Electrothermal (thermal_terms.cpp; PoissonCharge is shared).
    class HeatSinks;
    class LatticeStorage;
    class ThermalContactRows;
    class ThermalInsulatorRows;
    class ThermalGate;
    class ThermalSrh;
    class ThermalAuger;
    class ThermalRadiative;
    class ThermalCarrierStorage;
    class ThermalInterfaces;
    class ThermalEdges;
    class ThermalImpact;
};

// The gate tunnels' fluxes at one state (GR11; gate_tunnels.cpp): worker threads start on them at
// construction, claiming chunks of tunnels in turn; finish() has the calling thread take chunks too
// until none is left, waits for the workers and returns the fluxes, one per tunnel. Each flux is
// computed whole by one thread, so they do not depend on the thread count. Destruction waits for
// the workers when finish() was not called.
class DriftDiffusion::GateTunnelJob {
public:
    GateTunnelJob(const DriftDiffusion& s, std::span<const double> x,
                  std::span<const NodeDegeneracy> g);
    GateTunnelJob(const GateTunnelJob&) = delete;
    GateTunnelJob& operator=(const GateTunnelJob&) = delete;
    ~GateTunnelJob() = default;
    [[nodiscard]] std::vector<GateTunnelFlux>& finish();

private:
    void work();
    const DriftDiffusion& s_;
    std::span<const double> x_;
    std::span<const NodeDegeneracy> g_;
    std::vector<GateTunnelFlux> out_;
    std::size_t chunks_;
    std::atomic<std::size_t> next_{0};
    std::vector<std::jthread> workers_;  // last: joined before the rest is destroyed
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

// The names of DriftDiffusion's private term types, for the assemble library and its tests
// (this header is internal: not under include/).
struct TermAccess {
    using Term = DriftDiffusion::Term;
    using Terms = DriftDiffusion::Terms;
    using Composition = DriftDiffusion::Composition;
    using Evaluation = DriftDiffusion::Terms::Evaluation;
};

}  // namespace NiTCAD::assemble
