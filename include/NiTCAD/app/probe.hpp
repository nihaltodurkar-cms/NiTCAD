// Probe and measure (Unit 26, DECISIONS.md VA8): the mesh node under a point of the view, its position
// and field value, and the distance between two probed nodes.
//
// The nodes that can be picked are those the view shows: every node of a 2D map; in 3D the nodes of
// the boundary surface (the tetrahedra's outer faces), of which the front-most within the pick radius
// wins. A 1D device has none (its views are plots).
#pragma once

#include <array>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::app {

struct ProbeHit {
    mesh::NodeId node = 0;
    std::array<double, 3> position_um{};
    double value = 0.0;          // the field's value at the node (NaN when no values are given)
    float distance_px = 0.0f;    // from the pick point
};

class Prober {
public:
    explicit Prober(const Workspace& w);

    [[nodiscard]] std::span<const mesh::NodeId> candidates() const noexcept { return candidates_; }
    // The node drawn nearest the pixel (x, y) of the scene drawn at width x height, within the radius:
    // the nearest in 2D, the front-most in 3D. `values` has one entry per node, or is empty.
    [[nodiscard]] std::optional<ProbeHit> pick(const render::FieldScene& scene, int width, int height, float x,
                                               float y, std::span<const double> values = {},
                                               float radius_px = 12.0f) const;
    // A point of the device under the pixel, um (Unit 26b, a cut line's ends): in 2D the point itself, anywhere
    // in the view's plane; in 3D the front-most surface node within the radius (as pick).
    [[nodiscard]] std::optional<std::array<double, 3>> point_at(const render::FieldScene& scene, int width, int height,
                                                                float x, float y, float radius_px = 12.0f) const;

private:
    const mesh::Mesh* mesh_ = nullptr;
    int dimension_ = 0;
    std::vector<mesh::NodeId> candidates_;
};

[[nodiscard]] double distance_um(const ProbeHit& a, const ProbeHit& b) noexcept;

}  // namespace NiTCAD::app
