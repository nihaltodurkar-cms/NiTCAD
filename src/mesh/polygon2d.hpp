// Exact 2D polygon tests on double coordinates (orient2d only), shared by the complex, its
// generator and the interface rule (Unit 53). Internal to the mesh layer. Points carry their two
// coordinates in [0] and [1]; [2] is ignored.
#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh::detail {

// The point in the plane of the two coordinates other than `drop`, cyclically.
[[nodiscard]] Point project(const Point& p, int drop);

// Collinear in 3D, exactly: collinear in all three coordinate projections.
[[nodiscard]] bool collinear(const Point& a, const Point& b, const Point& c);

// Whether x lies in the closed 2D triangle abc (non-degenerate).
[[nodiscard]] bool in_triangle_2d(const Point& a, const Point& b, const Point& c, const Point& x);

// Whether x lies on the closed 2D segment ab (a != b).
[[nodiscard]] bool on_segment_2d(const Point& a, const Point& b, const Point& x);

// Whether the open 2D segments ab and cd cross at one point interior to both (not collinear).
[[nodiscard]] bool cross_properly_2d(const Point& a, const Point& b, const Point& c, const Point& d);

// Whether the closed 2D segments ab and cd intersect at all.
[[nodiscard]] bool segments_meet_2d(const Point& a, const Point& b, const Point& c, const Point& d);

// Whether a simple 2D polygon is traversed counterclockwise: the turn at its lexicographically
// smallest vertex (never zero for a simple polygon without reversals).
[[nodiscard]] int polygon_turn(const std::vector<Point>& q);

// Whether a 2D polygon is simple: consecutive edges meet only at their shared vertex (collinear
// continuation allowed, reversal not), other edges not at all.
[[nodiscard]] bool simple_2d(const std::vector<Point>& q);

// Whether p lies strictly inside a simple 2D polygon (crossing number with the half-open rule;
// points on the boundary are not inside).
[[nodiscard]] bool strictly_inside_2d(const std::vector<Point>& q, const Point& p);

// Whether p lies in the closed simple 2D polygon (inside or on its boundary).
[[nodiscard]] bool inside_closed_2d(const std::vector<Point>& q, const Point& p);

// Ear clipping of a simple 2D polygon with `turn` its orientation: an ear is a strictly convex
// vertex whose triangle (closed) holds no other remaining vertex, so no triangle is degenerate.
// Returns local index triples in the polygon's orientation, or nothing if no ear is found.
[[nodiscard]] std::optional<std::vector<std::array<std::size_t, 3>>> ear_clip(const std::vector<Point>& q, int turn);

}  // namespace NiTCAD::mesh::detail
