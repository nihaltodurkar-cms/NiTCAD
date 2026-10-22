// Mesh quality (Unit 26b; DECISIONS.md XD7): the shape of each linear element as an angle, its summary
// and histogram.
//
// A triangle's quality is its smallest interior angle; a tetrahedron's its smallest dihedral angle
// (between the planes of two faces along their common edge). Both are 60 degrees or less in 2D and
// arccos(1/3) = 70.5288 degrees or less in 3D (the equilateral triangle and the regular tetrahedron),
// and near 0 for a sliver. Angles are measured by atan2 of the sine and cosine parts, accurate to
// rounding at every angle.
#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/render/geometry.hpp"

namespace NiTCAD::render {

// The smallest angle, degrees, of a triangle (3 points) or the smallest dihedral angle of a
// tetrahedron (4 points). 0 for a degenerate element (coincident or collinear points, a flat
// tetrahedron).
[[nodiscard]] double min_angle_deg(std::span<const Point> vertices) noexcept;

// Per element of a 2D or 3D geometry, in element order. Errors (invalid_input): a 1D geometry or one
// without elements (a graph).
[[nodiscard]] std::expected<std::vector<double>, base::Error> element_quality(const Geometry& geometry);

// The angle below which an element counts as poor: 20 degrees in 2D, 10 in 3D.
[[nodiscard]] double quality_threshold_deg(int dimension) noexcept;
// The largest quality an element can have: 60 degrees in 2D, arccos(1/3) in 3D.
[[nodiscard]] double best_quality_deg(int dimension) noexcept;

struct QualitySummary {
    std::size_t elements = 0;
    double min_deg = 0.0;
    double median_deg = 0.0;  // the mean of the two middle values for an even count
    double threshold_deg = 0.0;
    std::size_t below = 0;    // elements strictly below the threshold
};

// The summary of per-element angles (none: every number 0 but the threshold).
[[nodiscard]] QualitySummary quality_summary(std::span<const double> angles, int dimension);

// Counts in 5-degree bins from 0 to 90 degrees: bin k holds [5k, 5k + 5); an angle of 90 or more is
// counted in the last bin, so every element is counted once.
struct QualityHistogram {
    static constexpr double bin_deg = 5.0;
    std::array<std::size_t, 18> counts{};
};
[[nodiscard]] QualityHistogram quality_histogram(std::span<const double> angles) noexcept;

}  // namespace NiTCAD::render
