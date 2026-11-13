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

}  // namespace NiTCAD::solve
