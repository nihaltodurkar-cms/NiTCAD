// The intelligent decision layer of adaptive refinement (Unit 17, slice 5; DECISIONS.md Unit 17
// ID1-ID7, ID9, AD12, AR5, AR7).
//
// The layer decides only where nodes are spent (ID1, OA6): the estimate, the stop rule, AD11's
// verification and the certificate are computed from solves alone and read nothing here.
//
// A priori prediction (ID2), before any solve, from the device and the analysis's biases:
// - junctions: every edge between semiconductor nodes of opposite net doping gives a junction
//   point (where linear interpolation of N_D - N_A crosses zero) with the edge's direction as its
//   normal. N_1, N_2: the largest |N| of each side's sign within the depletion width of the point
//   (a fixed point from the end nodes' values, three iterations); N_eff = N_1 N_2 / (N_1 + N_2);
//   V_bi = V_T ln(N_1 N_2 / n_i^2) (>= 0); W = sqrt(2 eps (V_bi + V) / (q N_eff)) with V the
//   largest difference between non-gate contact biases over the analysis; the target size is
//   L_D / 2, L_D = sqrt(eps V_T / (q N_eff)), for nodes within W of the point.
// - insulator interfaces under gates: a gate contact's nodes (its oxide thickness and
//   permittivity), and semiconductor nodes next to an insulator when the device has electrode
//   contacts (t_ox the distance to the nearest electrode node, the insulator's permittivity): the
//   surface field E_s = eps_ox V / (eps_s t_ox), V the largest difference between the gate (or an
//   electrode) and the other contacts over the analysis, gives the classical inversion length
//   lambda = V_T / E_s (with MLDA on, the larger of it and MLDA's electron wavelength); the target
//   size is lambda for nodes within 5 lambda; and the depletion region under it (AR5): L_D / 2,
//   L_D at the node's |N|, for nodes within W_d = sqrt(2 eps_s psi / (q N)), psi = min(2 phi_F, V).
//   The normal points into the semiconductor (the mean direction to the node's semiconductor
//   neighbours off the interface).
// - contacts: none. Impact ionization and band-to-band tunnelling add nothing beyond their
//   junction's target (AR5).
// A node's target is the smallest over the features covering it, with that feature's normal. On a
// tensor or swept axis an edge needs refinement when its length projected on the normal exceeds the
// target at either end (a junction or an interface is resolved across, not along, itself); on a
// generated 2D mesh or cross-section the target enters the size field as is (isotropic).
//
// Bias prediction (ID3): from one coarse sweep on the start mesh, per point the united indicator
// excess; the points are the last one and those whose field changed most from the previous point
// (the weighted 2-norm of the difference, weights 1), up to max_points, in sweep order.
//
// Benefit, cost and allocation (ID4-ID6, AR7; in adapt.cpp, with adjoint weighting only: the benefit
// is the adjoint estimate's share): candidates are axis intervals (tensor axes, the swept axis; up
// to three levels of halving per pass) and nodes (a generated 2D size field, a swept cross-section;
// one level). A candidate's benefit is E (S_c / S) 2^-(l-1)p (1 - 2^-p) for its level l, E the
// relative estimate, S_c its sum of weight times excess, p its fitted order; the allocation ranks
// by benefit over a deterministic work, N_newton ((W + dW)^1.5 - W^1.5) with W = nnz + fill, and
// takes candidates until the predicted estimate is at most min(1/2, E/2), the pass's DOF growth
// cap (2.5) or the budget, both held on the built mesh (a proposal over one is allocated again with
// grading's measured inflation). The timing-fitted cost model below predicts each candidate's and
// pass's seconds (reported, AC15) and enforces a time budget; it never ranks, so the mesh is
// deterministic (AD9).
#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <functional>
#include <limits>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/solve/adapt_indicators.hpp"
#include "NiTCAD/solve/bias.hpp"

namespace NiTCAD::solve {

enum class Policy : unsigned char { basic, intelligent };

[[nodiscard]] constexpr const char* to_string(Policy p) noexcept {
    return p == Policy::basic ? "basic" : "intelligent";
}

// Per node: a target size [cm] (infinity where nothing is predicted) and the unit normal of the
// feature that gave it (all zero: isotropic).
struct SizeTargets {
    std::vector<double> size;
    std::vector<std::array<double, 3>> normal;
};

struct Prediction {
    SizeTargets targets;
    std::size_t junction_points = 0;
    std::size_t interface_points = 0;
    std::size_t predicted_nodes = 0;  // nodes with a finite target
};

// ID2's feature points, before they become per-node targets (exposed for Unit 17b's learning
// features, DECISIONS.md LB3). Junction points: span W, size L_D / 2; interface points: span 5 lambda,
// size lambda; their depletion region (AR5): span W_d, size L_D / 2. v_junction and v_gate: the
// analysis's largest difference between non-gate contacts and between a gate-like contact and the
// rest (V). Errors (invalid_input): points not one bias per contact.
struct FeaturePoint {
    enum class Kind : unsigned char { junction, interface, depletion };
    Kind kind = Kind::junction;
    mesh::Point at{};
    std::array<double, 3> normal{};
    double span = 0.0;
    double size = 0.0;
};

struct FeaturePoints {
    std::vector<FeaturePoint> points;
    double v_junction = 0.0;
    double v_gate = 0.0;
};

[[nodiscard]] std::expected<FeaturePoints, base::Error> feature_points(const device::Device& device,
                                                                       std::span<const std::vector<double>> points,
                                                                       const BiasOptions& options);

// ID2. Errors (invalid_input): points not one bias per contact.
[[nodiscard]] std::expected<Prediction, base::Error> predict_sizes(const device::Device& device,
                                                                   std::span<const std::vector<double>> points,
                                                                   const BiasOptions& options);

// Whether edge e of the device's mesh is longer than the targets ask (projected on the normal).
[[nodiscard]] bool edge_exceeds(const device::Device& device, const SizeTargets& targets, std::size_t e);

// ID3: the adaptation points from per-point indicator fields (one per sweep point, same size),
// among the eligible points (empty: all; the loop excludes points where a goal is at or below its
// current resolution, which cannot carry a relative tolerance, AR5): the last eligible point and
// the eligible points of largest change. Errors (invalid_input): no fields, fields of different
// sizes, max_points 0, eligibility not one flag per point, no eligible point.
[[nodiscard]] std::expected<std::vector<std::size_t>, base::Error> choose_points(
    std::span<const IndicatorField> fields, std::size_t max_points, std::span<const char> eligible = {});

// ID5's cost model (OA4), per pass: T_pass = N_newton (a nnz^b + c fill + d DOF) + T_mesh + T_transfer,
// T_mesh and T_transfer per node. One sample per solved pass.
struct CostSample {
    double newton_iterations = 0.0;  // summed over the sweep's points
    double matrix_nonzeros = 0.0;    // the Jacobian's
    double fill = 0.0;               // nnz(L + U); 0 when the backend does not report it
    double dofs = 0.0;
    double nodes = 0.0;
    double seconds_solve = 0.0;
    double seconds_mesh = 0.0;
    double seconds_transfer = 0.0;
};

struct CostModel {
    // Per Newton iteration: a nnz^b + c fill + d DOF (a, c, d >= 0; b in [1, 2]).
    double a = 0.0, b = 1.5, c = 0.0, d = 0.0;
    double mesh_per_node = 0.0, transfer_per_node = 0.0;
    std::size_t samples = 0;
    double residual = 0.0;  // root mean square relative error of the per-iteration fit
    // The node-count-only model fitted to the same passes (AC15's comparison): T_pass = alpha N^beta.
    double alpha = 0.0, beta = 1.0;
    double residual_nodes = 0.0;
    [[nodiscard]] double per_iteration(double nnz, double fill, double dofs) const;
    // A pass's predicted time (NaN without samples).
    [[nodiscard]] double pass_seconds(double newton_iterations, double nnz, double fill, double dofs, double nodes) const;
    [[nodiscard]] double pass_seconds_nodes(double nodes) const;
};

// The fit (ID5, AR7): with one sample a alone (b = 1.5, c = d = 0); with two, a and b through both
// (b clamped); with three or more, b on a grid of 0.05 in [1, 2] and a, c, d by non-negative least
// squares of the relative error, plus a ridge (0.01 / n) (b - 1.5)^2 fading with the n samples. T_mesh and T_transfer per node: their
// sums over the nodes' (transfer: of the passes that transferred). The node-only model: beta from a
// least-squares line in log-log (1 with one sample; clamped to [0.5, 3]), alpha by relative least squares.
[[nodiscard]] CostModel fit_cost_model(std::span<const CostSample> samples);

// ID5's local order: a refined candidate's realized ratio of its local error after / before (AR7).
struct OrderSample {
    Indicator kind = Indicator::potential;
    std::size_t region = 0;
    double ratio = 1.0;
};

struct OrderFit {
    Indicator kind = Indicator::potential;
    std::size_t region = 0;
    double order = 1.0;
    std::size_t samples = 0;
};

// Per (kind, region), the least squares of -log2 ratio with a ridge of one sample toward 1, clamped to
// [0.5, 3]; samples with a ratio not finite and positive are skipped. Sorted by kind, then region.
[[nodiscard]] std::vector<OrderFit> fit_orders(std::span<const OrderSample> samples);
// The fitted order of (kind, region), 1 when none.
[[nodiscard]] double order_of(std::span<const OrderFit> fits, Indicator kind, std::size_t region);

// The layer's options (adapt.hpp's AdaptOptions::intelligent).
struct PolicyOptions {
    // ID3: with AdaptOptions::at empty, up to this many points chosen from a coarse sweep; 0 keeps
    // AdaptOptions' default points (the first and last).
    std::size_t max_points = 0;
    // A relative fall the estimate must make for a proposal to count as lowering it (ID1, AR5).
    double miss_margin = 0.1;
    // A hand-written prediction in place of the layer's own (ID2 at the first mesh, and every later
    // proposal): per-node targets on the given device. Used for AC17's wrong prediction; empty: the
    // layer's.
    std::function<SizeTargets(const device::Device&)> prediction;
    // Unit 17b (LB4): a learned prediction in place of ID2's for the start mesh only (the prediction
    // rounds before pass 0); every later pass is the layer's own (or `prediction`'s, when set: AC21's
    // fixture combines the two).
    std::function<SizeTargets(const device::Device&)> first_prediction;
    // ID6's budget mode: the run's DOFs (final mesh) and its time (measured passes plus the next one
    // predicted) at most these; 0: tolerance mode. A pass whose allocation fits nothing stops the run
    // (max_nodes, unverified).
    std::size_t dof_budget = 0;
    double seconds_budget = 0.0;
};

}  // namespace NiTCAD::solve
