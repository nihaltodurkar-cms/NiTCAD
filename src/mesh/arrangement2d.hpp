// The planar arrangement behind layered_partition2d's slanted path, labelled by the caller (Unit 26l slice 4, SD6).
// Internal to the mesh layer.
//
// The shapes' edges and the domain's sides are split where they meet and each split edge is labelled side by side,
// as the slanted 2D path does, but the label of a side is `label(holding, at)`: `holding` the indices of the shapes
// whose inside that side is (ascending; a shape's own edge counts by its direction), `at` the edge's midpoint. A
// label is >= 0, or -1 for a point left uncovered (an error). Each label's parts become simple polygons as in 2D
// (loops traced, holes bridged, collinear split points dropped), the polygons' region being the label. `pins` are
// points on the domain's sides kept as polygon vertices (a side split there). The shapes' regions and the
// background are not used; the shapes are taken as given (the caller has validated them).
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <functional>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/generate.hpp"

namespace NiTCAD::mesh::detail {

struct ArrangementLabels {
    std::function<std::int32_t(const std::vector<std::size_t>& holding, const Point& at)> label;
    std::vector<Point> pins;
};

[[nodiscard]] std::expected<Geometry2D, base::Error> labelled_partition2d(const LayeredGeometry2D& shapes, const ArrangementLabels& how);

}  // namespace NiTCAD::mesh::detail
