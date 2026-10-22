// A field along a straight line through the device (Unit 26b; DECISIONS.md XD5).
//
// The segment from a to b is clipped by every element it crosses. Barycentric weights are affine
// along it, so each element holds one interval of the line, between the points where a weight
// reaches 0; the samples are the ends of these intervals, the crossings with element faces (one
// sample where neighbouring elements meet). A nodal field read at the samples by their elements'
// weights is exactly the piecewise-linear field along the line: it is linear between two samples,
// with kinks at the samples only. Where the line runs outside the mesh there is a gap: the next
// sample is marked. A line along a face shared by two elements reads either (the field is
// continuous there).
//
// Coordinates are the mesh's (cm); only the mesh's dimensions of a and b are used.
#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/render/geometry.hpp"

namespace NiTCAD::render {

struct CutSample {
    double distance = 0.0;  // from a, cm
    Interpolant at;         // weights summing to 1, zero for the nodes off the sample's face
    std::size_t element = 0;
    bool gap_before = false;  // the line ran outside the mesh between the previous sample and this
};

struct CutLine {
    Point a{};
    Point b{};
    std::vector<CutSample> samples;  // by distance; empty when the line misses the mesh
};

// Errors (invalid_input): a geometry without elements (a 2D or 3D graph); a and b equal in the
// mesh's dimensions, or not finite.
[[nodiscard]] std::expected<CutLine, base::Error> cut_line(const Geometry& geometry, const Point& a,
                                                           const Point& b);
// A nodal field at the samples (one value per sample). Nodes of zero weight are not read, so a NaN
// node (bands on an insulator) leaves the samples on the faces away from it defined.
[[nodiscard]] std::vector<double> cut_values(const CutLine& cut, std::span<const double> field);

}  // namespace NiTCAD::render
