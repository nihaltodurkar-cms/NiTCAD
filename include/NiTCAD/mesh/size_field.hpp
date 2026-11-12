// Size fields for adaptive refinement (Unit 17, slice 1; DECISIONS.md Unit 17 AD4, AD10, AQ1).
//
// A SizeField is a nodal element size h on a background triangulation of the cross-section or 2D
// domain (normally the previous pass's generated mesh). The 2D generator, and through it the swept
// generator's cross-section, take it as MeshOptions2D::size_field: the size at a point is then the
// smallest of max_size, the region's size, the boxes' and the field's, so a field only refines. The
// field is evaluated by exact point location on the background (orientation predicates; a point on
// a shared edge or vertex takes the lowest-indexed containing triangle) and linear (barycentric)
// interpolation of the nodal sizes. A point outside every triangle (rounding at an oblique
// boundary) takes the value at its nearest point of the nearest triangle (lowest index on ties).
// With no field, every generator, digest and result is unchanged (AD10).
//
// Tensor and swept axes (AD4): refine_axis splits marked intervals at their midpoints and then
// grades the axis so neighbouring intervals differ by at most the grading factor; intervals inside
// a protected range (an interface straddle) are never split.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <memory>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"

namespace NiTCAD::mesh {

class SizeField {
public:
    // Errors (invalid_input, the context index naming the entry): no points or no triangles; sizes
    // not one per point; a point not finite or with z != 0; a triangle index out of range or
    // repeated; a flat triangle (exact orientation); a size not finite and positive.
    [[nodiscard]] static std::expected<std::shared_ptr<const SizeField>, base::Error> create(
        std::vector<Point> points, std::vector<NodeId> triangles, std::vector<double> sizes);

    // The size at p (cm).
    [[nodiscard]] double at(const Point& p) const;

    [[nodiscard]] const std::vector<Point>& points() const { return points_; }
    [[nodiscard]] const std::vector<NodeId>& triangles() const { return triangles_; }
    [[nodiscard]] const std::vector<double>& sizes() const { return sizes_; }
    [[nodiscard]] double min_size() const;

    // A 64-bit FNV-1a digest of the points, triangles and sizes (doubles by bit pattern).
    [[nodiscard]] std::uint64_t digest() const;

private:
    SizeField() = default;
    [[nodiscard]] double interpolate(std::size_t t, const Point& p) const;

    std::vector<Point> points_;
    std::vector<NodeId> triangles_;  // counterclockwise after create
    std::vector<double> sizes_;
    // A uniform bucket grid over the bounding box; each bucket lists the triangles whose bounding
    // box meets it, in increasing index.
    Point lo_{};
    double cell_ = 1.0;
    std::size_t nx_ = 1;
    std::size_t ny_ = 1;
    std::vector<std::uint32_t> bucket_start_;
    std::vector<std::uint32_t> bucket_items_;
};

// The local size of every node of a 2D simplex mesh: its longest incident edge (cm).
[[nodiscard]] std::vector<double> local_sizes(const SimplexMesh& simplices);

// Grades nodal sizes on a 2D or 3D simplex mesh's edges in place, from the seed nodes (all nodes
// when `seeds` is empty): each node takes the smallest of its own size and g^k h_s over the seeds s
// k edges away, so sizes only decrease and, with every node a seed, h_j <= g h_i on every edge.
// With seeds (the nodes a pass refined, AD4), the grading follows the refinement only: ratios the
// mesh already had (interface pairs beside coarse elements) do not spread. Errors (invalid_input):
// g < 1 or not finite; sizes (or non-empty seeds) not one per point, or a size not finite and
// positive.
[[nodiscard]] std::expected<void, base::Error> grade_sizes(const SimplexMesh& simplices, std::vector<double>& sizes,
                                                           double g, const std::vector<char>& seeds = {});

// A closed range of an axis where no line may be inserted (an interface straddle).
struct AxisRange {
    double lo = 0.0;
    double hi = 0.0;
};

struct AxisRefinement {
    std::vector<double> lines;
    std::size_t split = 0;      // intervals split because marked
    std::size_t graded = 0;     // intervals split by grading
    std::size_t protected_skips = 0;  // marked or graded intervals left whole: inside a protected range
};

// Splits each marked interval [lines[i], lines[i+1]] (marked.size() == lines.size() - 1) at its
// midpoint, then repeatedly splits any interval longer than max(g, 2) times a neighbour made by a
// split in this call (fixed point; ratios of the given axis are kept; a midpoint split always leaves
// a ratio of 2, so a smaller factor would cascade over the axis), never
// below min_interval and never inside a protected range (an interval meeting one in more than an
// end point). Deterministic. Errors (invalid_input): lines not finite or not strictly increasing,
// fewer than two; marked of the wrong size; g < 1 or not finite; min_interval negative.
[[nodiscard]] std::expected<AxisRefinement, base::Error> refine_axis(const std::vector<double>& lines,
                                                                     const std::vector<char>& marked, double g,
                                                                     double min_interval,
                                                                     const std::vector<AxisRange>& protected_ranges = {});

}  // namespace NiTCAD::mesh
