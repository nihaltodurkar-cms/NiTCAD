// The adaptive refinement loop (Unit 17, slice 3; DECISIONS.md Unit 17 AD1, AD4-AD6, AD9, AD11,
// AD13; AQ1-AQ3).
//
// adapt() solves the analysis (a bias sweep, one mesh for the whole sweep, AQ3) on the device the
// builder makes from the deck's mesh, then repeats: indicators at the adaptation points (AD2),
// united, marked (AD3), a refinement of the mesh (AD4: a graded size field for a generated 2D mesh
// and a swept mesh's cross-section, midpoint lines for tensor axes and the swept axis), the device
// rebuilt by the builder with it (AQ1: regenerated, so every certificate comes from the generators'
// checked path), the previous state transferred as Newton's initial state (AD5), and the sweep
// solved again. The physics, the solver and the analysis are not changed: the builder and the
// sweep are passed in (AD1).
//
// Goals (AD6): terminal currents of named contacts at the adaptation points (the weighted current,
// Unit 32, when the options compute it; else the terminal current). The loop proposes to stop when
// every goal's change from the previous pass is within tolerance, |G_k - G_k-1| <= tol max(|G_k|,
// r) with r the point's current resolution, and otherwise ends at max_passes, max_nodes (the next
// mesh would exceed it; the current one is kept) or when nothing is marked. The estimate of a goal
// is its last change with `weighting = none` (feature-only marks, w = 1). With `weighting =
// adjoint` (slice 4; adapt_adjoint.hpp, DECISIONS.md AR3) every pass builds M' (AD11's mesh), and
// per adaptation point computes the adjoint estimate 2 |D| on it and the dual-weighted residual's
// node weights: the stop then also needs every estimate within tolerance (AD6), the marks are
// e_i w_i with w the normalized weights (the largest over goals and points; AD3), and AD11 reuses
// that M'.
//
// Verification (AD11): after a proposed stop, the final mesh M is refined uniformly once (every
// size halved, every tensor and swept-axis interval outside the straddles split) to M', the sweep
// solved on M' from the transferred state, and the verified error is |G_M' - G_M| / (2^p - 1), p
// the order observed from the last two passes and M' (log2 of the ratio of the changes, clamped to
// [0.5, 3]) when the loop made at least two passes, else 1. It reads no indicator, weight or
// prediction. The result stays on M. A stop is accepted only when, for every goal, the estimate and
// the verified error are both within tol max(|G|, r); a rejected stop (DECISIONS.md AR2) is recorded
// and the loop goes on refining (by its marks, or to M' when nothing is marked) within max_passes
// and max_nodes. Two consecutive small changes can be stagnation, not convergence: this is the
// guard.
//
// The certificate (AD13), per goal and adaptation point, from the final mesh's verification:
// certified when the estimate and the verified error are both within tol max(|G|, r); not certified
// when either exceeds it (a stop rejected at the last pass, or a budget stop right after one);
// unverified when the final mesh was not verified (verification off, or a budget stop: max_passes,
// max_nodes, nothing marked).
//
// Determinism (AD9): every step is deterministic, so two runs on one machine and thread count are
// bit-identical.
#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/size_field.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/solve/adapt_indicators.hpp"
#include "NiTCAD/solve/adapt_policy.hpp"
#include "NiTCAD/solve/adapt_transfer.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/control.hpp"

namespace NiTCAD::solve {

// What the builder is asked for, relative to the deck's own mesh; empty members leave it as given.
struct MeshRefinement {
    std::shared_ptr<const mesh::SizeField> size_field;  // generated 2D, or a swept mesh's cross-section
    std::array<std::vector<double>, 3> axes;            // tensor: the whole axes, straddles included
    std::vector<double> sweep_lines;                    // swept: the given lines (the straddles' added by the generator)
};

enum class AdaptMeshKind : unsigned char { tensor, generated2d, swept };

// A device as the builder made it, with what the loop needs to refine its mesh.
struct AdaptableMesh {
    device::Device device;
    AdaptMeshKind kind = AdaptMeshKind::tensor;
    std::optional<mesh::GeneratedMesh> generated;  // generated 2D and swept
    std::optional<mesh::SimplexMesh> section;      // swept: the cross-section's 2D mesh
    std::array<std::vector<mesh::AxisRange>, 3> axis_protected;  // tensor: the straddles per axis
    int sweep_axis = 2;                                 // swept: the axis (0, 1, 2)
    std::vector<mesh::AxisRange> sweep_protected;       // swept: the straddles along the axis
    double interface_resolution = 0.0;                  // generated 2D and swept: s (0 without interfaces)
    double max_size = 0.0;                              // generated 2D and swept: the deck's h
};

using DeviceBuilder = std::function<std::expected<AdaptableMesh, base::Error>(const MeshRefinement&)>;

struct AdaptGoal {
    std::size_t contact = 0;  // a terminal current, device contact order
};

enum class Weighting : unsigned char { none, adjoint };

struct AdaptOptions {
    std::vector<AdaptGoal> goals;           // at least one
    double tolerance = 1e-3;                // relative, in (0, 1)
    std::vector<std::size_t> at;            // adaptation points, indices into the sweep; empty: its first and last
    IndicatorOptions indicators;
    Weighting weighting = Weighting::none;  // adjoint: AD3's weights and AD6's estimate (AR3)
    double theta = 0.3;                     // AD3
    double grading = 1.5;                   // AD4, g >= 1
    double min_size = 0.0;                  // cm; 0: s / 4 with interfaces, else h / 64 (axes: per axis, see adapt.cpp)
    int max_passes = 6;                     // solves on adapted meshes, the first mesh included; >= 1
    std::size_t max_nodes = 0;              // 0: no limit
    bool verify = true;                     // AD11
    Policy policy = Policy::basic;          // the decision layer (adapt_policy.hpp); basic: AD3 alone
    PolicyOptions intelligent;              // read only with policy = intelligent
};

// How a pass's mesh was made (ID7).
enum class MeshSource : unsigned char {
    start,       // the builder's mesh as given
    prediction,  // the start refined toward the a priori prediction (ID2)
    proposal,    // the decision layer's proposal (ID1)
    marks,       // AD3's marks
    difference,  // AD3 on the verified error's distribution (AR2, AD12 c)
    uniform,     // M' after a rejected stop with nothing marked
};

[[nodiscard]] constexpr const char* to_string(MeshSource s) noexcept {
    switch (s) {
        case MeshSource::start: return "start";
        case MeshSource::prediction: return "prediction";
        case MeshSource::proposal: return "proposal";
        case MeshSource::marks: return "marks";
        case MeshSource::difference: return "difference";
        case MeshSource::uniform: return "uniform";
    }
    return "?";
}

// AD12's triggers.
enum class FallbackTrigger : unsigned char { none, two_misses, no_progress, rejected_stop };

[[nodiscard]] constexpr const char* to_string(FallbackTrigger t) noexcept {
    switch (t) {
        case FallbackTrigger::none: return "none";
        case FallbackTrigger::two_misses: return "two_misses";
        case FallbackTrigger::no_progress: return "no_progress";
        case FallbackTrigger::rejected_stop: return "rejected_stop";
    }
    return "?";
}

enum class AdaptStop : unsigned char { converged, max_passes, max_nodes, nothing_marked };

[[nodiscard]] constexpr const char* to_string(AdaptStop s) noexcept {
    switch (s) {
        case AdaptStop::converged: return "converged";
        case AdaptStop::max_passes: return "max_passes";
        case AdaptStop::max_nodes: return "max_nodes";
        case AdaptStop::nothing_marked: return "nothing_marked";
    }
    return "?";
}

struct AdaptPass {
    std::size_t nodes = 0;
    std::size_t edges = 0;
    // Per goal and adaptation point, goal-major (goal g, point a at g * at.size() + a).
    std::vector<double> goal;
    std::vector<double> resolution;  // the point's current resolution, same layout
    std::vector<double> change;      // |G_k - G_k-1|; empty on the first pass
    // weighting = adjoint (AR3): the adjoint estimate 2 |D| per goal and point (same layout; empty
    // without weighting), M''s node count it was computed on, its time (building M' included) and
    // whether a transposed solve fell back to SparseLU (gmres_ilut).
    std::vector<double> estimate;
    std::size_t estimate_nodes = 0;
    double seconds_estimate = 0.0;
    bool estimate_fell_back = false;
    Marking marking;                 // the marks this pass made for the next (empty count on the last)
    std::array<std::size_t, indicator_count> nodes_with_excess{};
    int newton_iterations = 0;       // summed over the sweep's points
    bool rejected = false;           // a proposed stop on this mesh that AD11 rejected (AR2)
    bool marked_by_difference = false;  // the marks came from the verified error's distribution
    // The decision layer (ID1, ID7): how this pass's mesh was made; the largest over goals and points
    // of the estimate over its allowance (infinity when there is none); whether this pass's mesh was
    // a proposal that missed (did not lower it by the miss margin, or refined nothing).
    MeshSource source = MeshSource::start;
    double relative_estimate = 0.0;
    bool miss = false;
    bool empty_proposal = false;  // the proposal made on this pass refined nothing (a miss; AD3's marks used)
    TransferReport transfer;         // onto this pass's mesh (nodes 0 on the first)
    double seconds_mesh = 0.0;       // building this pass's device
    double seconds_transfer = 0.0;
    double seconds_solve = 0.0;
};

enum class CertificateStatus : unsigned char { certified, not_certified, unverified };

[[nodiscard]] constexpr const char* to_string(CertificateStatus s) noexcept {
    switch (s) {
        case CertificateStatus::certified: return "certified";
        case CertificateStatus::not_certified: return "not_certified";
        case CertificateStatus::unverified: return "unverified";
    }
    return "?";
}

struct GoalCertificate {
    std::size_t goal = 0;
    std::size_t point = 0;    // index into the sweep
    double value = 0.0;       // on the final mesh
    double estimate = 0.0;    // AD6's: the adjoint's 2 |D|, or the last change (infinity after a single pass)
    double verified_error = 0.0;  // AD11 (NaN when unverified)
    double order = 1.0;           // p used
    double allowed = 0.0;         // tol max(|G|, r)
    CertificateStatus status = CertificateStatus::unverified;
};

struct AdaptResult {
    std::vector<AdaptPass> passes;
    AdaptStop stop = AdaptStop::converged;
    std::size_t rejected_stops = 0;     // proposed stops AD11 rejected
    std::optional<AdaptableMesh> mesh;  // the final mesh
    results::Sweep sweep;               // the sweep on the final mesh
    // AD11's solve on M' (absent when not run): its node count, time and goals (same layout).
    std::size_t verification_nodes = 0;
    double verification_seconds = 0.0;
    std::vector<double> verification_goal;
    std::vector<GoalCertificate> certificates;
    // The decision layer (ID2, ID3, ID7, AD12).
    std::vector<std::size_t> points;    // the adaptation points used
    Policy policy = Policy::basic;      // as asked
    FallbackTrigger fallback = FallbackTrigger::none;
    int fallback_pass = -1;             // the pass whose outcome switched the run to basic
    std::size_t misses = 0;
    std::size_t coarse_nodes = 0;       // ID3's coarse sweep (0 when not run)
    double coarse_seconds = 0.0;
    std::size_t junction_points = 0;    // ID2 on the start mesh
    std::size_t interface_points = 0;
    std::size_t predicted_nodes = 0;    // start-mesh nodes with a finite target
    int prediction_rounds = 0;          // refinements toward the prediction before pass 0
};

// The refinement of `current` for the next pass (AD4), from per-node marks; exposed for tests.
// `uniform` asks for AD11's M' instead (every size halved, every unprotected interval split).
// `previous` is the refinement `current` was built with. Errors: those of the size field and axis
// refinement; invalid_input for marks not one per node.
[[nodiscard]] std::expected<MeshRefinement, base::Error> next_refinement(
    const AdaptableMesh& current, const MeshRefinement& previous, const std::vector<char>& marked,
    const EdgeIndicators* edges, const IndicatorOptions& indicators, const AdaptOptions& options, bool uniform);

// Errors (invalid_input, before solving): no points; no goals, a goal contact out of range or a
// gate; tolerance not in (0, 1); an adaptation point out of range; theta not in (0, 1]; grading
// below 1; min_size negative; max_passes below 1; weighting adjoint with check_adjoint's refusals
// (adapt_adjoint.hpp) or a goal contact that is not ohmic. Then those
// of the builder, the sweep (a sweep that stops early is an error: non_convergence and the rest),
// the transfer and the refinement.
[[nodiscard]] std::expected<AdaptResult, base::Error> adapt(const DeviceBuilder& build,
                                                            std::span<const std::vector<double>> points,
                                                            const BiasOptions& bias, const AdaptOptions& options,
                                                            const RunControl& control = {});

}  // namespace NiTCAD::solve
