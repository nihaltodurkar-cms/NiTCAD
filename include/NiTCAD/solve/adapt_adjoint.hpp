// Goal adjoints and the adjoint error estimate (Unit 17, slice 4; DECISIONS.md Unit 17 AQ4, AD3,
// AD6, AR3).
//
// The goal G of contact c at a converged state x of mesh M is its weighted current (Unit 32 W1, when
// the options compute it) or its terminal current T_c. With w the weighted current's weights (1 on
// c's nodes, 0 on the other ohmic contacts', the conductance-weighted Laplace solution elsewhere),
//     G(x) = sum over edges of (Jn + Jp)_e (w_a - w_b) = T_c(x) + sum_i w_i B_i(x),
// the sum over the semiconductor nodes off the ohmic contacts, B_i = F_n,i + F_p,i the node's
// total-current balance (its electron and hole rows: generation and recombination cancel). Its
// gradient is dT_c/dx (the assembler's linearized conduction current, small_signal_currents at s =
// 0, with every model's chain rule) plus w_i times those two Jacobian rows. The adjoint z solves
//     J^T z = dG/dx,
// one transposed solve with J's factorization (linalg solve_transposed; gmres_ilut falls back to
// SparseLU for it). It is computed as z = y + w on the electron and hole rows, J^T y = dT_c/dx: the
// w part is exact (J^T w = sum w_i (rows)), and y carries no cancellation of large terms.
//
// The estimate (AR3; the dual-weighted residual): M' is M refined uniformly (AD11's mesh), u the
// state of M carried onto M' by the Scharfetter-Gummel prolongation (adapt_transfer.hpp
// exponential_fit), R' = F'(u) the residual of M''s equations there (both meshes scaled alike: M's
// concentration scale on M'). The adjoint z' of M''s goal is computed at u on M' (one factorization
// of J'(u), one transposed solve per goal), and
//     D = G'(u) - G_M - z'^T R'
// is the goal's change over Newton's first correction on M', linearized (G'(u) - G_M, the
// prolongated state's own change, is at rounding level for the prolongation: it reproduces M's
// fluxes). D estimates G_M' - G_M without solving on M'. The correction itself is not applied: from
// a prolongated state it is large and clamped, and G' after it measured orders off. With p = 1
// (AD11's default, the box method's order at abrupt features), the error of M is
// 2^p / (2^p - 1) |D| = 2 |D|: the estimate AD6 and AD13 use. The local contribution of a fine
// node, z'_q^T R'_q (the dual-weighted residual), is summed onto the coarse node of its largest
// interpolation weight in its own region; the magnitude of a coarse node's sum is AD3's weight
// before normalization.
// Why on M', not M's adjoint interpolated: M''s residual of a prolongated state moves the
// recombination between old and new nodes (each of order h R), which cancels only to O(h^p)
// against a smooth weight; an interpolated coarse adjoint is not smooth to that order, and
// measured on the legacy diode its estimate grew with refinement at -1 V (to 370 times the actual
// change at 321 nodes), where the adjoint on M' tracked it.
//
// An estimate, not a bound: the second-order terms beyond one Newton correction are left out
// (AC7 measures it against the actual error).
#pragma once

#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/locate.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/bias.hpp"

namespace NiTCAD::solve {

// The adjoints of goals at one converged point of mesh M (scaled units; block layout of the
// drift-diffusion unknowns, three per node).
struct PointAdjoint {
    std::vector<double> state;                   // x
    std::vector<std::vector<double>> weights;    // per goal: w (weighted current), empty otherwise
    std::vector<double> goal;                    // per goal: G(x)
    std::vector<std::vector<double>> gradient;   // per goal: dG/dx
    std::vector<std::vector<double>> adjoint;    // per goal: z
    std::vector<std::vector<double>> terminal;   // per goal: y (the terminal current's adjoint)
    double current_scale = 0.0;                  // public units per scaled current
    bool fell_back = false;                      // gmres_ilut: the transposed solve by SparseLU
};

// Refusals (invalid_input) shared by the adjoint functions and adapt()'s weighting = adjoint: not
// the drift-diffusion equations; electrothermal, energy balance or nonlocal tunnelling; a linear
// solver that cannot transpose (gmres_ilut without its fallback).
[[nodiscard]] std::optional<base::Error> check_adjoint(const BiasOptions& options);

// Errors: check_adjoint's; invalid_input for a contact out of range or not ohmic, or a point whose
// fields are not one per node; those of the assembly's bias, the factorization and the solves.
[[nodiscard]] std::expected<PointAdjoint, base::Error> point_adjoint(const device::Device& device,
                                                                     const BiasOptions& options,
                                                                     const results::BiasPoint& point,
                                                                     std::span<const std::size_t> contacts);

// G of a contact at state x with the given weights (empty: the terminal current), scaled; for
// tests of the gradient. Preconditions as point_adjoint's.
[[nodiscard]] std::expected<double, base::Error> goal_at(const device::Device& device, const BiasOptions& options,
                                                         std::span<const double> bias_V, std::span<const double> x,
                                                         std::size_t contact, std::span<const double> weights);

struct GoalEstimate {
    std::vector<double> difference;                 // per goal: D (public units)
    std::vector<double> estimate;                   // per goal: 2 |D|
    std::vector<std::vector<double>> node_weight;   // per goal: per node of M, |sum of its contributions|
    std::vector<double> prolongated_part;           // per goal: G'(u) - G_M (public units)
    bool fell_back = false;
};

// The estimate on `fine` (M') of the goals at `point` of `coarse` (M), whose locator is given.
// Errors: point_adjoint's and the transfer's.
[[nodiscard]] std::expected<GoalEstimate, base::Error> estimate_goals(const device::Device& coarse,
                                                                      const mesh::PointLocator& coarse_locator,
                                                                      const results::BiasPoint& point,
                                                                      const device::Device& fine,
                                                                      const BiasOptions& options,
                                                                      std::span<const std::size_t> contacts);

}  // namespace NiTCAD::solve
