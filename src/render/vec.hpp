// 3D vector helpers on render::Point and the barycentric weights of a point in a simplex, shared by the
// render sources (internal to src/render; not an installed header).
#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <span>

#include "NiTCAD/render/geometry.hpp"

namespace NiTCAD::render::vec {

inline Point sub(const Point& a, const Point& b) noexcept { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Point add(const Point& a, const Point& b) noexcept { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Point scale(const Point& a, double s) noexcept { return {a[0] * s, a[1] * s, a[2] * s}; }
inline double dot(const Point& a, const Point& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Point cross(const Point& a, const Point& b) noexcept {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline double norm(const Point& a) noexcept { return std::sqrt(dot(a, a)); }
inline Point lerp(const Point& a, const Point& b, double t) noexcept { return add(a, scale(sub(b, a), t)); }

// The barycentric weights of x in a simplex of 2, 3 or 4 nodes (a segment along x, a triangle in the xy plane, a
// tetrahedron), as affine functions extended beyond it (Cramer's rule): weights[k] for element[k], unused slots 0,
// summing to 1. Nothing for a degenerate simplex (zero length, area or volume).
inline std::optional<std::array<double, 4>> barycentric(std::span<const Point> points, std::span<const mesh::NodeId> element,
                                                        const Point& x) noexcept {
    const std::size_t k = element.size();
    const Point& x0 = points[static_cast<std::size_t>(element[0])];
    std::array<Point, 3> col{};
    for (std::size_t i = 1; i < k; ++i) col[i - 1] = sub(points[static_cast<std::size_t>(element[i])], x0);
    const Point r = sub(x, x0);
    std::array<double, 4> w{};
    if (k == 2) {
        if (col[0][0] == 0.0) return std::nullopt;
        w[1] = r[0] / col[0][0];
    } else if (k == 3) {
        const double det = col[0][0] * col[1][1] - col[1][0] * col[0][1];
        if (det == 0.0) return std::nullopt;
        w[1] = (r[0] * col[1][1] - col[1][0] * r[1]) / det;
        w[2] = (col[0][0] * r[1] - r[0] * col[0][1]) / det;
    } else {
        const double det = dot(col[0], cross(col[1], col[2]));
        if (det == 0.0) return std::nullopt;
        w[1] = dot(r, cross(col[1], col[2])) / det;
        w[2] = dot(col[0], cross(r, col[2])) / det;
        w[3] = dot(col[0], cross(col[1], r)) / det;
    }
    w[0] = 1.0;
    for (std::size_t i = 1; i < k; ++i) w[0] -= w[i];
    return w;
}

}  // namespace NiTCAD::render::vec
