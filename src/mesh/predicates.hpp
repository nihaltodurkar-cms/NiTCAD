// Exact geometric predicates and their deterministic tie-break (Unit 16; DECISIONS.md M4).
// Internal to the mesh layer.
//
// Each predicate returns the exact sign (-1, 0 or +1) of a polynomial in the input coordinates,
// for any finite doubles (no range restriction). Two stages:
// - a filter: the polynomial in double precision, with a forward error bound computed alongside it
//   (|computed - exact| <= gamma_d * M, with M the polynomial evaluated on absolute values and d
//   accumulated per operation, max + 1 for sums and a + b + 1 for products; the standard
//   rounding-error induction, as in Higham, Accuracy and Stability of Numerical Algorithms, ch. 3).
//   The sign is returned when |computed| exceeds 2 d 2^-53 M, unless a product underflowed or a
//   value overflowed (the rounding model then fails);
// - otherwise the same polynomial in exact integers (exact_integer.hpp), over the inputs scaled by
//   a common power of two.
// Both stages evaluate one shared formula, so they cannot disagree on what is computed.
//
// Conventions (2D predicates read x and y only):
// - orient2d(a, b, c) > 0 when a, b, c turn counterclockwise (twice the signed area);
// - orient3d(a, b, c, d) = det[b - a; c - a; d - a] > 0 when d is on the side of the plane abc
//   toward which (b - a) x (c - a) points (positive volume, right-handed);
// - incircle(a, b, c, d) > 0 when d is strictly inside the circle through a, b, c and
//   orient2d(a, b, c) > 0; the sign flips with the orientation; 0 when cocircular;
// - insphere(a, b, c, d, e) > 0 when e is strictly inside the sphere through a, b, c, d and
//   orient3d(a, b, c, d) > 0; the sign flips with the orientation; 0 when cospherical;
// - in_diametral_circle(a, b, p) > 0 when p is strictly inside the circle with diameter ab;
// - in_diametral_sphere(a, b, c, p) > 0 when p is strictly inside the smallest sphere through the
//   triangle abc (centre at its circumcentre). Precondition: a, b, c not collinear.
// incircle and insphere are the lifted determinants of their 4 or 5 points, so they are
// antisymmetric in their arguments: an odd permutation flips the sign.
//
// Tie-break (symbolic perturbation of the lifting; Devillers and Teillaud, "Perturbations for
// Delaunay and weighted Delaunay 3D triangulations", 2011, [cited, to verify]): each point's lifted
// value |p|^2 is raised by eps^r, r its rank in decreasing node id, eps infinitesimal. The lifted
// determinant is linear in the lifted values, so the perturbed sign is the sign of the first nonzero
// cofactor in decreasing id order. The orientation predicates are never perturbed, so a tie never
// creates a flat element. The perturbed predicates never return 0; they require distinct ids and
// a non-degenerate base simplex (NITCAD_EXPECTS).
#pragma once

#include <array>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh::predicates {

[[nodiscard]] int orient2d(const Point& a, const Point& b, const Point& c);
[[nodiscard]] int orient3d(const Point& a, const Point& b, const Point& c, const Point& d);
[[nodiscard]] int incircle(const Point& a, const Point& b, const Point& c, const Point& d);
[[nodiscard]] int insphere(const Point& a, const Point& b, const Point& c, const Point& d,
                           const Point& e);
[[nodiscard]] int in_diametral_circle(const Point& a, const Point& b, const Point& p);
[[nodiscard]] int in_diametral_sphere(const Point& a, const Point& b, const Point& c,
                                      const Point& p);

// Precondition: orient2d(p[0], p[1], p[2]) != 0; ids distinct.
[[nodiscard]] int incircle_perturbed(const std::array<Point, 4>& p,
                                     const std::array<NodeId, 4>& id);
// Precondition: orient3d(p[0], p[1], p[2], p[3]) != 0; ids distinct.
[[nodiscard]] int insphere_perturbed(const std::array<Point, 5>& p,
                                     const std::array<NodeId, 5>& id);

// The exact stage alone, without the filter; for tests that check the filter against it.
namespace exact {
[[nodiscard]] int orient2d(const Point& a, const Point& b, const Point& c);
[[nodiscard]] int orient3d(const Point& a, const Point& b, const Point& c, const Point& d);
[[nodiscard]] int incircle(const Point& a, const Point& b, const Point& c, const Point& d);
[[nodiscard]] int insphere(const Point& a, const Point& b, const Point& c, const Point& d,
                           const Point& e);
[[nodiscard]] int in_diametral_sphere(const Point& a, const Point& b, const Point& c,
                                      const Point& p);
}  // namespace exact

}  // namespace NiTCAD::mesh::predicates
