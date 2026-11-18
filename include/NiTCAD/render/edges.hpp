// Feature edges and the enhanced 3D look (Unit 26f; DECISIONS.md KD1, KD3, KQ2-KQ4).
//
// feature_edges finds the edges of a triangle batch worth drawing as lines: an edge two faces
// share whose normals differ by more than an angle (either winding, so slices and interfaces wound
// either way count once), and an edge of only one face (the batch's open boundary) or of more than
// two. Vertices are matched by their exact positions; degenerate triangles are ignored. The result
// is sorted, each segment with its lower end first, so drawing it is deterministic whatever the
// triangle order.
//
// set_look sets a scene's look; on a 3D scene made enhanced it appends one LineBatch, "feature
// edges", built once from the region solids, field surfaces, slices and electrodes (not arrows,
// streamline tubes, isosurfaces or the outline shell), so a camera move redraws without rebuilding
// them. Setting the look again leaves a single such batch.
#pragma once

#include <array>
#include <vector>

#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::render {

using Segment = std::array<Float3, 2>;

inline constexpr float feature_angle_deg = 30.0f;  // KQ2

[[nodiscard]] std::vector<Segment> feature_edges(const TriangleBatch& batch, float angle_deg = feature_angle_deg);

// Whether a batch takes feature edges in the enhanced look (KQ3).
[[nodiscard]] bool takes_feature_edges(const TriangleBatch& batch) noexcept;

void set_look(FieldScene& scene, Look look);

}  // namespace NiTCAD::render
