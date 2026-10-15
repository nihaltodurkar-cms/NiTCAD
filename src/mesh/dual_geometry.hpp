// Circumcentric dual pieces of single simplices (Unit 16; DECISIONS.md M1). Internal to the mesh
// layer; make_dual_mesh (simplex_mesh.hpp) sums them.
//
// Every piece is signed: for an obtuse simplex some are negative (the circumcentre lies outside),
// and only sums over a Delaunay mesh are guaranteed non-negative (M2, M3). Determinants, dot
// products and circumcentre numerators are computed exactly (exact_integer.hpp) and rounded once,
// so a piece is accurate to a few ulps of its own terms however flat the simplex; edge lengths are
// std::hypot of the coordinate differences, as Mesh::from_parts checks them.
//
// 2D (per unit depth): the coupling piece of edge ij in triangle ijk is the signed distance from
// the edge midpoint to the circumcentre, L_ij cot(theta_k) / 2.
// 3D: the coupling piece of edge ij in tetrahedron ijkl is the signed area of the dual polygon
// (midpoint, circumcentre of face ijk, circumcentre of the tetrahedron, circumcentre of face ijl),
// (h_k d_k + h_l d_l) / 2, with h_f the in-face distance from the midpoint to face f's circumcentre
// and d_f the signed distance from face f to the tetrahedron's circumcentre (toward the fourth
// vertex). A vertex's volume piece is the sum over its edges of L piece / (2 D); its boundary piece
// on a facet is its share of the facet's own circumcentric partition (L / 2 per endpoint in 2D).
#pragma once

#include <array>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh::detail {

// Local edges of a triangle and of a tetrahedron, in this order.
inline constexpr std::array<std::array<int, 2>, 3> triangle_edges{{{0, 1}, {0, 2}, {1, 2}}};
inline constexpr std::array<std::array<int, 2>, 6> tetrahedron_edges{
    {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};

// Preconditions (NITCAD_EXPECTS): the simplex is positively oriented (orient2d / orient3d > 0).
[[nodiscard]] std::array<double, 3> triangle_couplings(const std::array<Point, 3>& p);
[[nodiscard]] std::array<double, 6> tetrahedron_couplings(const std::array<Point, 4>& p);

// The boundary piece of each vertex of a 3D boundary triangle (area, cm^2). Precondition: not
// collinear.
[[nodiscard]] std::array<double, 3> facet_pieces(const std::array<Point, 3>& p);

// The circumcentre of a 3D triangle (for neighbour searches and P4; rounded). Precondition: not
// collinear. In 2D (z = 0) it is the triangle's circumcentre in the plane.
[[nodiscard]] Point triangle_circumcentre(const std::array<Point, 3>& p);
// The circumcentre of a tetrahedron (rounded once from exact numerators). Precondition: not flat.
[[nodiscard]] Point tetrahedron_circumcentre(const std::array<Point, 4>& p);

[[nodiscard]] double distance(const Point& a, const Point& b) noexcept;

}  // namespace NiTCAD::mesh::detail
