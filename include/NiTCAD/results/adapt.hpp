// The adaptation report (Unit 17, slice 6; DECISIONS.md Unit 17 AD8, AD13, ID7, AR8): what an
// adaptive run did, pass by pass, and the accuracy certificate of its result.
//
// Plain data, filled by solve::adapt_report from solve::AdaptResult, saved by io::save_adapt_report
// and shown by `nitcad run`, `nitcad plot --view passes|adapt` and the app. Text fields hold the
// solve layer's names (to_string of its enums), so this layer knows nothing of solve.
//
// The certificate (AD13), per goal and adaptation point: the value, the estimate (the adjoint's
// 2 |D|, or the last change without weighting), the verified error of AD11 and the order it used, the
// allowance tol max(|G|, r), and the status: "certified" (estimate and verified error both within the
// allowance), "not_certified" (either exceeds it; the margin is the larger of the two over the
// allowance) or "unverified" (verification off, or a budget stop). The report's status is the
// worst over its certificates (not_certified, then unverified, then certified). A certificate is
// evidence from two independent estimates, not a proof of the error bound.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "NiTCAD/results/run.hpp"

namespace NiTCAD::results {

// One candidate the allocation took (ID4, ID7).
struct AdaptCandidateRow {
    std::int64_t axis = -1;          // -1: a node; else the interval's axis
    std::uint64_t index = 0;
    double x = 0.0, y = 0.0, z = 0.0;  // cm
    std::uint64_t region = 0;
    std::string indicator;
    double order = 1.0;
    std::int64_t levels = 1;
    double benefit = 0.0;            // predicted, relative to the allowance
    double cost = 0.0;               // predicted seconds (reported only, AR7)
    double work = 0.0;               // the allocation's deterministic cost
    std::uint64_t added_nodes = 0, added_dofs = 0, added_nonzeros = 0;
    double realized_ratio = 0.0;     // NaN until measured
    double realized = 0.0;           // NaN until measured
    bool operator==(const AdaptCandidateRow&) const = default;
};

// A fitted local order (ID5).
struct AdaptOrderRow {
    std::string indicator;
    std::uint64_t region = 0;
    double order = 1.0;
    std::uint64_t samples = 0;
    bool operator==(const AdaptOrderRow&) const = default;
};

struct AdaptPassRow {
    std::uint64_t nodes = 0, edges = 0, dofs = 0;
    std::uint64_t mesh_digest = 0;               // mesh::graph_digest of the pass's mesh
    std::string source;                          // how the mesh was made ("start", "proposal", ...)
    // Per goal and adaptation point, goal-major (as AdaptReport::goal_contacts x points).
    std::vector<double> goal, resolution, change, estimate;
    double relative_estimate = 0.0;              // the largest estimate over its allowance
    std::uint64_t estimate_nodes = 0;            // M''s nodes (adjoint weighting)
    // The marks this pass made for the next mesh: their count, by dominant indicator (indexed as
    // AdaptReport::indicators) and by device region.
    std::uint64_t marked = 0;
    std::vector<std::uint64_t> marked_by_indicator, marked_by_region;
    std::vector<std::uint64_t> nodes_with_excess;  // per indicator
    std::int64_t newton_iterations = 0;
    std::uint64_t matrix_nonzeros = 0, factor_fill = 0;
    bool rejected = false;                       // a proposed stop AD11 rejected
    bool miss = false;                           // a proposal that did not lower the estimate (ID1)
    // The decision layer (ID4-ID7): the predicted time by the cost model and the node-only model,
    // the allocation's counts and predictions, the candidates listed (at most 20) and the orders.
    double predicted_seconds = 0.0, predicted_seconds_nodes = 0.0;  // NaN on the first pass
    std::uint64_t candidate_count = 0, candidates_taken = 0;
    double predicted_relative = 0.0;             // NaN without an allocation
    std::uint64_t predicted_dofs = 0;
    std::vector<AdaptCandidateRow> candidates;
    std::vector<AdaptOrderRow> orders;
    // The cost model fitted up to this pass (ID5): a nnz^b + c fill + d DOF per Newton iteration,
    // mesh and transfer per node, its sample count and residual; the node-only alpha N^beta.
    double cost_a = 0.0, cost_b = 1.5, cost_c = 0.0, cost_d = 0.0;
    double cost_mesh_per_node = 0.0, cost_transfer_per_node = 0.0, cost_residual = 0.0;
    double cost_alpha = 0.0, cost_beta = 1.0, cost_residual_nodes = 0.0;
    std::uint64_t cost_samples = 0;
    double seconds_mesh = 0.0, seconds_transfer = 0.0, seconds_solve = 0.0, seconds_estimate = 0.0;
    bool operator==(const AdaptPassRow&) const = default;
};

struct AdaptCertificateRow {
    std::uint64_t goal = 0;      // index into AdaptReport::goal_contacts
    std::uint64_t point = 0;     // index into the sweep
    double value = 0.0;
    double estimate = 0.0;
    double verified_error = 0.0;  // NaN when unverified
    double order = 1.0;
    double allowed = 0.0;
    std::string status;           // "certified", "not_certified", "unverified"
    bool operator==(const AdaptCertificateRow&) const = default;
};

struct AdaptReport {
    RunRecord run;  // the final sweep's record (identity stamped with the settings, AR8)
    // The settings (AD7): the goals' contacts (names, device contact order), the tolerance, the
    // adaptation points, the policy and weighting asked, verification on, certification required.
    std::vector<std::string> goal_contacts;
    double tolerance = 0.0;
    std::vector<std::uint64_t> points;
    std::string policy, weighting;
    bool verify = true;
    bool certify_required = false;
    std::vector<std::string> indicators;  // the indicator kinds, in the rows' order
    // The outcome: the stopping reason, rejected stops, the fallback (AD12) and its pass, misses.
    std::string stop;
    std::uint64_t rejected_stops = 0;
    std::string fallback;                 // "none" without one
    std::int64_t fallback_pass = -1;
    std::uint64_t misses = 0;
    // ID2 and ID3 before pass 0.
    std::uint64_t junction_points = 0, interface_points = 0, predicted_nodes = 0;
    std::int64_t prediction_rounds = 0;
    std::uint64_t coarse_nodes = 0;
    double coarse_seconds = 0.0;
    // AD11's verification solve on M' (nodes 0 when not run), and the goals it gave.
    std::uint64_t verification_nodes = 0;
    double verification_seconds = 0.0;
    std::vector<double> verification_goal;
    std::vector<AdaptPassRow> passes;
    std::vector<AdaptCertificateRow> certificates;
    std::string status;                   // the worst certificate's
    std::uint64_t final_mesh_digest = 0;  // the last pass's mesh: the result's
    bool operator==(const AdaptReport&) const = default;
};

}  // namespace NiTCAD::results
