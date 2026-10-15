// Private to the solve layer: the weighted terminal current (Unit 32; DECISIONS.md W1-W8).
//
// In a steady state the total current J = J_n + J_p is conserved at every node that is not an
// ohmic contact node (generation and recombination create pairs, so they cancel in the total).
// For any weight w on the nodes that is 1 on contact c's nodes and 0 on the other ohmic contacts'
// nodes, summing w_i times each node's total-current balance gives, by summation by parts,
//     I_c = sum over edges of J_e (w_a - w_b),
// exactly the contact's current when every balance vanishes. The choice of w only decides where
// the current is read. Here w solves the conductance-weighted Laplace problem
//     sum over the edges at i of g_e (w_i - w_j) = 0
// at every other semiconductor node, with g_e half the edge's one-sided flux terms (its
// conductance to a quasi-Fermi difference near equilibrium; assemble edge_flux_terms): w is flat
// where the conductance is large (heavily doped regions, whose edge currents are small
// differences of huge fluxes) and varies where it is small (depletion regions, a channel), where
// the current is resolved. A floor of 1e-30 of the largest conductance keeps the problem regular;
// any w is exact, so the floor costs nothing.
//
// The resolution: the adjoint (dual-weighted residual) estimate of the Newton residual's effect on
// I_c is grad(I_c) J^-1 F, which is the change of I_c over the Newton correction dx = -J^-1 F, so it
// is measured as |I_c(x + dx) - I_c(x)| with w held fixed (one factorization at the state, no
// transposed solve and no gradient of the fluxes); plus the rounding of the weighted edge
// currents, 8 eps sum over edges of t_e |w_a - w_b|. If the correction cannot be computed, the
// conventional terminal_current_resolution stands in for the estimate.
#pragma once

#include <cstdint>
#include <expected>
#include <optional>
#include <utility>
#include <span>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/bias.hpp"

namespace NiTCAD::solve::detail {

// The option's own rule (invalid_input), checked before anything is solved: drift-diffusion only,
// and not with nonlocal tunnelling, whose pairs are created at different nodes (holes at a path's
// start, electrons at its crossing), so the total current is not carried by the edges alone.
inline std::optional<base::Error> check_weighted_current(const BiasOptions& options) {
    if (!options.weighted_current) return std::nullopt;
    const auto invalid = [](const char* m) {
        return base::Error{base::ErrorCode::invalid_input, m, std::nullopt};
    };
    if (options.equations != Equations::drift_diffusion) {
        return invalid("weighted_current needs the drift-diffusion equations: a quasi-static sweep "
                       "carries no current");
    }
    if (options.models.btbt_nonlocal != assemble::NonlocalTunnelling::off) {
        return invalid("weighted_current does not support nonlocal tunnelling: its pairs are not "
                       "carried along the mesh edges (DECISIONS.md W6)");
    }
    return std::nullopt;
}

class WeightedCurrents {
public:
    // Errors: those of linalg::LinearSolver::create.
    [[nodiscard]] static std::expected<WeightedCurrents, base::Error> create(
        const device::Device& device, const linalg::SolverConfig& linear);

    // The weights of each contact at state x (the inner vector per node; empty for a contact that
    // is not ohmic). Errors: those of the linear solve (singular_system, inaccurate_solve).
    [[nodiscard]] std::expected<std::vector<std::vector<double>>, base::Error> weights(
        const assemble::DriftDiffusion& system, std::span<const double> x);

    // Scaled weighted currents of x with these weights (0 for a contact without weights).
    [[nodiscard]] std::vector<double> currents(const assemble::DriftDiffusion& system,
                                               std::span<const double> x,
                                               const std::vector<std::vector<double>>& w) const;

    // point.weighted_current and point.weighted_current_resolution of the converged state x,
    // in the public units; point.terminal_current_resolution must be filled. Errors: those of
    // weights().
    [[nodiscard]] std::expected<void, base::Error> fill(const assemble::DriftDiffusion& system,
                                                        std::span<const double> x,
                                                        double max_update,
                                                        const assemble::Scaling& scaling,
                                                        int dimension, results::BiasPoint& point);

private:
    WeightedCurrents(linalg::LinearSolver laplace, linalg::LinearSolver newton)
        : laplace_(std::move(laplace)), newton_(std::move(newton)) {}

    const device::Device* device_ = nullptr;
    std::vector<std::int32_t> ohmic_;   // per node: its ohmic contact, or -1
    std::vector<std::int32_t> row_;     // per node: its unknown in the Laplace problem, or -1
    std::vector<std::uint8_t> carrier_;  // per node: a semiconductor node
    std::size_t rows_ = 0;
    linalg::LinearSolver laplace_;
    linalg::LinearSolver newton_;
};

}  // namespace NiTCAD::solve::detail
