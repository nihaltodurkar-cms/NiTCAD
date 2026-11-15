// Cross-run learning for adaptive refinement: features, labels and the history of a run (Unit 17b,
// slice 1; DECISIONS.md Unit 17 ID8, Unit 17b LB2, LB3, LB11). The features are schema 1's
// (results/learn.hpp), per edge of a device's mesh, from ID2's feature points (feature_points) and
// the device's doping, materials and contacts; the labels the realized refinement against a final
// mesh. Deterministic: every search takes the lowest index on ties.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/locate.hpp"
#include "NiTCAD/results/learn.hpp"
#include "NiTCAD/solve/adapt.hpp"
#include "NiTCAD/solve/bias.hpp"

namespace NiTCAD::solve {

// A material's category key: its kind and the parameters that tell materials apart (eps_r; for a
// semiconductor also Eg0, Nc300 and Nv300), as exact hexadecimal floats.
[[nodiscard]] std::string material_key(const device::Region& region);
// The device's material keys, sorted, without duplicates.
[[nodiscard]] std::vector<std::string> material_keys(const device::Device& device);

// FNV-1a 64 over the device's mesh digest, node regions, donors, acceptors and contacts' nodes and
// kinds (doubles by bit pattern).
[[nodiscard]] std::uint64_t device_digest(const device::Device& device);

// Schema 1's features of every edge, with the edges' geometry for labelling.
struct EdgeFeatures {
    std::size_t rows = 0;
    std::vector<double> features;                 // rows x results::learn_feature_count
    std::vector<mesh::Point> midpoint;            // cm
    std::vector<std::array<double, 3>> direction; // unit
    std::vector<double> length;                   // cm
};

// Feature 1 is the index of the edge's material key in `materials` (-1 when absent: an unseen
// category, LB8). An edge's material: its semiconductor end's (the lower node's when both or
// neither are). Errors: those of feature_points.
[[nodiscard]] std::expected<EdgeFeatures, base::Error> edge_features(const device::Device& device,
                                                                     std::span<const std::vector<double>> points,
                                                                     const BiasOptions& options,
                                                                     std::span<const std::string> materials);

// The realized refinement of each edge against a final mesh (LB2): max(0, log2(L_e / h_f)), NaN
// where the midpoint is outside every final cell.
[[nodiscard]] std::vector<double> realized_levels(const EdgeFeatures& edges, const mesh::PointLocator& final_mesh,
                                                  const mesh::Mesh& final_points);

// The locator of an adaptable mesh (tensor cells, or the generated simplices).
[[nodiscard]] std::expected<mesh::PointLocator, base::Error> adaptable_locator(const AdaptableMesh& mesh);

// Gathers a run's history (LB11): record() each solved pass as AdaptOptions::on_mesh gives it, then
// finish() with the final mesh. The materials are the first pass's device's keys plus any later
// pass's, sorted; features are recomputed against them at finish.
class HistoryRecorder {
public:
    HistoryRecorder(std::vector<std::vector<double>> points, BiasOptions options);
    [[nodiscard]] std::expected<void, base::Error> record(std::size_t pass, const AdaptableMesh& mesh);
    [[nodiscard]] std::expected<results::AdaptHistory, base::Error> finish(const AdaptableMesh& final_mesh) const;

private:
    struct Pass {
        std::uint64_t device_digest = 0, mesh_digest = 0;
        std::vector<std::string> materials;
        EdgeFeatures edges;
    };
    std::vector<std::vector<double>> points_;
    BiasOptions options_;
    std::vector<Pass> passes_;
};

// --- Slice 2: the model (LB5, LB7, LB10; ID9, ID10) ---------------------------------------------

// LB5's fixed constants.
struct TreeOptions {
    std::size_t trees = 200;
    std::size_t depth = 4;
    double shrinkage = 0.1;
    std::size_t min_leaf = 20;
    std::size_t cut_points = 64;
};

// LB5's split candidates of one feature: the values sorted, the cut at rank floor(j n / (m + 1)) for
// j = 1..m (m = max_cuts), duplicates removed, a cut equal to the largest value dropped (nothing
// would go right); with at most m distinct values, every distinct value but the largest.
[[nodiscard]] std::vector<double> cut_points(std::span<const double> values, std::size_t max_cuts);

// Squared-loss gradient boosting (LB5) on the rows `rows` of a training set's features (row-major,
// results::learn_feature_count columns). Returns the trees and the base; each node's split the one
// of largest squared-error decrease with both sides at least min_leaf, ties to the lower feature,
// then the lower cut; a node without such a split, or at the depth, is a leaf of the mean residual
// times the shrinkage. `importance` (per feature) receives the summed decreases.
struct BoostedTrees {
    double base = 0.0;
    std::vector<std::vector<results::TreeNode>> trees;
    std::vector<double> importance;
};
[[nodiscard]] BoostedTrees fit_boosted_trees(std::span<const double> features, std::span<const double> labels,
                                             std::span<const std::size_t> rows, const TreeOptions& options = {});
[[nodiscard]] double predict_trees(double base, std::span<const std::vector<results::TreeNode>> trees, const double* row);

// The ridge baseline (ID9, LB5): continuous features standardized by the rows' mean and deviation
// (a feature of zero deviation dropped), categories one-hot over the values in the rows, an
// unpenalized intercept, lambda = 1, by Cholesky.
struct Ridge {
    std::vector<double> mean, scale;                 // per feature (continuous; scale 0: dropped)
    std::vector<std::vector<double>> categories;     // per feature (categories)
    std::vector<double> coefficients;                // intercept first, then the columns in feature order
};
[[nodiscard]] std::expected<Ridge, base::Error> fit_ridge(std::span<const double> features, std::span<const double> labels,
                                                          std::span<const std::size_t> rows, double lambda = 1.0);
[[nodiscard]] double predict_ridge(const Ridge& ridge, const double* row);

// Spearman's rank correlation (ranks averaged over ties; 0 when either side is constant). Errors: none;
// sizes must match (NITCAD_EXPECTS).
[[nodiscard]] double rank_correlation(std::span<const double> a, std::span<const double> b);

// LB7's default held-out devices of n: every fourth from the fourth (indices 3, 7, ...); with fewer
// than four, the last. Empty for n < 2.
[[nodiscard]] std::vector<std::size_t> default_held_out(std::size_t devices);

// Trains the predictor (LB5, LB7, LB10): trees on the training devices' samples, the ridge baseline on
// the same, both measured on the held-out devices' samples; ranges, categories and the manifest's
// digest recorded; the provenance this build's. Errors (invalid_input): a held-out index out of range
// or repeated, no training or no held-out device, no training samples or no held-out samples, another
// schema.
[[nodiscard]] std::expected<results::LearnedModel, base::Error> train_predictor(const results::TrainingSet& set,
                                                                                std::span<const std::size_t> held_out,
                                                                                const TreeOptions& options = {});

// The model's l^ for one schema-1 row.
[[nodiscard]] double predict_levels(const results::LearnedModel& model, const double* row);

// --- Slice 3: the predictor in the loop (LB4, LB8) -----------------------------------------------

// LB8: whether a row is outside the model's training range (a continuous feature outside [min - 5%
// of the range, max + 5%], or a category not seen), and the fraction of a device's edges that are.
[[nodiscard]] bool out_of_range(const results::LearnedModel& model, const double* row);
[[nodiscard]] double out_of_range_fraction(const results::LearnedModel& model, const EdgeFeatures& edges);
inline constexpr double out_of_range_limit = 0.05;  // LB8: a device above it falls back to ID2

// The run identity with a predictor (ID10, LB9): input_identity continued by FNV-1a 64 over
// "adapt.predictor", the model file's digest, its format and schema and whether it was used; the
// settings gain adapt.predictor.digest_high and _low (its 32-bit halves), .format, .schema and .used.
void stamp_predictor(results::RunRecord& run, std::uint64_t file_digest, const results::LearnedModel& model, bool used);

// LB4's targets on a device: per edge L_e 2^-l^ where l^ >= 0.5; per node the smallest over its edges,
// with that edge's direction as the normal; infinity (no target) elsewhere. Errors: those of
// edge_features.
[[nodiscard]] std::expected<SizeTargets, base::Error> learned_targets(const results::LearnedModel& model,
                                                                      const device::Device& device,
                                                                      std::span<const std::vector<double>> points,
                                                                      const BiasOptions& options);

}  // namespace NiTCAD::solve
