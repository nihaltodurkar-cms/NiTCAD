// Cross-run learning data (Unit 17b; DECISIONS.md Unit 17 ID8, ID10, Unit 17b LB2, LB3, LB6, LB11).
//
// Plain data: an adaptation history (what one adaptive run teaches: per solved pass, one sample per
// edge, its features and its realized refinement against the run's final mesh), and the manifest
// and training set built from a list of histories. Filled by solve/adapt_learn.hpp, saved and read
// by io/learn_files.hpp.
//
// Schema 1's features (LB3), per edge, index: name.
//   0 region_kind (category: 0 semiconductor, 1 insulator: both ends insulator)
//   1 material (category: the index of the material's key in the list it is read with)
//   2 log10 max(|N_D - N_A|, n_i) (the mean over the semiconductor ends; 0 on an insulator edge)
//   3 |grad log10 max(|N|, n_i)| along the edge, per um (0 unless both ends are semiconductor)
//   4 log10 distance to the nearest junction point [cm] (0, i.e. 1 cm, when there is none)
//   5 log10 distance to the nearest insulator-interface point [cm] (0 when none)
//   6 log10 distance to the nearest contact node [cm] (0 when none)
//   7 log10 W of the nearest junction point [cm] (0 when none)
//   8 log10 L_D of the nearest junction point [cm] (0 when none)
//   9 log10 lambda of the nearest interface point [cm] (0 when none)
//  10 the largest junction bias difference at one of the analysis's points [V] (ID2's)
//  11 the largest gate bias difference at one of its points [V] (ID2's)
//  12 the goal kind (0: a terminal current)
//  13 log10 of the edge's length [cm]
//  14 |cos| between the edge and the nearest junction or interface point's normal (0 when none)
//  15 log10 (the nearest junction's distance / its W) (6 when none)
//  16 log10 (the nearest interface point's distance / its lambda) (6 when none)
// Distances are from the edge's midpoint. The label (LB2): l = max(0, log2(L_e / h_f)), h_f the final
// mesh's extent along the edge's direction in the cell holding its midpoint; NaN when the midpoint
// lies outside every final cell.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/results/run.hpp"

namespace NiTCAD::results {

inline constexpr std::uint64_t learn_schema_version = 1;
inline constexpr std::size_t learn_feature_count = 17;

inline constexpr std::array<std::string_view, learn_feature_count> learn_feature_names{
    "region_kind",         "material",           "log_doping",         "log_doping_gradient",
    "log_junction_distance", "log_interface_distance", "log_contact_distance", "log_junction_width",
    "log_junction_debye",  "log_interface_lambda", "junction_bias",     "gate_bias",
    "goal_kind",           "log_edge_length",    "normal_cosine",      "log_junction_scaled",
    "log_interface_scaled"};

// Whether feature k is a category (LB3, LB5, LB8).
[[nodiscard]] constexpr bool learn_feature_is_category(std::size_t k) noexcept { return k == 0 || k == 1; }

// One solved pass of an adaptive run.
struct HistoryPass {
    std::uint64_t device_digest = 0;  // solve::device_digest of the pass's device
    std::uint64_t mesh_digest = 0;    // io::mesh_digest of its mesh
    std::uint64_t rows = 0;           // its edges
    std::vector<double> features;     // rows x learn_feature_count, row-major
    std::vector<double> labels;       // per edge; NaN: not labelled
    bool operator==(const HistoryPass&) const = default;
};

struct AdaptHistory {
    std::uint64_t schema = learn_schema_version;
    std::string deck_text;               // the deck the run read
    RunRecord run;                       // the final sweep's run record
    std::vector<std::string> materials;  // sorted material keys; feature 1 indexes this list
    std::vector<HistoryPass> passes;     // in pass order
};

// A manifest (LB6): the histories of a training set, by file digest, sorted by digest.
struct ManifestEntry {
    std::uint64_t digest = 0;  // FNV-1a 64 of the whole history file
    std::string path;
    bool operator==(const ManifestEntry&) const = default;
};

struct Manifest {
    std::uint64_t schema = learn_schema_version;
    std::vector<ManifestEntry> entries;
    bool operator==(const Manifest&) const = default;
};

// The samples of a manifest (LB6): its histories' labelled edges in manifest order, within a history
// in pass then edge order, at most `max_per_device` per history by a fixed stride; feature 1
// re-indexed to `materials`, the sorted union of the histories' keys.
struct TrainingSet {
    std::uint64_t schema = learn_schema_version;
    std::vector<std::string> materials;
    std::vector<ManifestEntry> devices;  // the manifest's entries; a sample's device indexes them
    std::size_t rows = 0;
    std::vector<double> features;        // rows x learn_feature_count
    std::vector<double> labels;
    std::vector<std::uint32_t> device;
    std::vector<std::uint32_t> pass;
    bool operator==(const TrainingSet&) const = default;
};

// The model file's format (ID10); a model of another format or schema is refused.
inline constexpr std::uint64_t learn_model_format = 1;

// A regression tree node: an internal node sends x[feature] <= threshold left; a leaf (feature -1)
// holds its value (shrinkage applied).
struct TreeNode {
    std::int64_t feature = -1;
    double threshold = 0.0;
    std::int64_t left = -1, right = -1;
    double value = 0.0;
    bool operator==(const TreeNode&) const = default;
};

// The held-out metrics (ID10, LB10), on the held-out devices' samples.
struct HeldOutMetrics {
    std::uint64_t samples = 0;           // all held-out samples (rank correlations)
    std::uint64_t first_pass_samples = 0; // their first pass's (recall, precision)
    double rank_correlation = 0.0;       // Spearman, trees' prediction against the label
    double ridge_rank_correlation = 0.0; // the ridge baseline's (ID9)
    double recall = 0.0;                 // of l >= 0.5 by l^ >= 0.5, first pass
    double precision = 0.0;
    bool operator==(const HeldOutMetrics&) const = default;
};

// A trained cross-run predictor (ID8, ID10, LB5): gradient-boosted trees predicting l from schema
// `schema`'s features, with what it was trained on and how it did on held-out devices.
struct LearnedModel {
    std::uint64_t format = learn_model_format;
    std::uint64_t schema = learn_schema_version;
    std::string model_class = "gradient_boosted_trees";
    // LB5's hyperparameters, as trained.
    std::uint64_t tree_count = 0, depth = 0, min_leaf = 0, cut_points = 0;
    double shrinkage = 0.0;
    double base = 0.0;                                 // the training labels' mean
    std::vector<std::vector<TreeNode>> trees;          // node 0 the root
    std::vector<std::string> materials;                // feature 1's categories index this list
    std::vector<double> feature_min, feature_max;      // per feature, over the training samples
    std::vector<std::vector<double>> categories;       // per feature: the values seen (categories only)
    std::vector<double> importance;                    // per feature: summed squared-error decrease
    std::uint64_t manifest_digest = 0;                 // FNV-1a 64 over the manifest's entry digests
    std::vector<ManifestEntry> training, held_out;     // the devices of each side
    std::uint64_t training_samples = 0;
    Provenance provenance;                             // the build that trained it
    HeldOutMetrics metrics;
    bool operator==(const LearnedModel&) const = default;
};

}  // namespace NiTCAD::results
