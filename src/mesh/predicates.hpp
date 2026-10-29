// Exact geometric predicates and their deterministic tie-break (Unit 16; DECISIONS.md M4).
// Internal to the mesh layer.
//
// Each predicate returns the exact sign (-1, 0 or +1) of a polynomial in the input coordinates,
// for any finite doubles (no range restriction). Two stages:
// - a filter: the polynomial in double precision, with a forward error bound carried alongside;
//   the sign is returned when |computed| > 2 d u m (u = 2^-53, m and d below), otherwise it defers;
// - the same polynomial in exact integers (exact_integer.hpp), over the inputs scaled by a common
//   power of two: exact by construction, for any finite inputs.
// Both stages evaluate one shared formula, so they cannot evaluate different polynomials.
//
// Filter soundness. Arithmetic model: IEEE 754 binary64, round to nearest even, gradual underflow
// (no flush-to-zero or denormals-are-zero), no extended precision and no contraction across the
// filter's operations (MSVC x64 with /fp:precise; nothing in the project changes the
// floating-point environment). Each node of the formula's tree carries v (its computed value), m
// (the computed value of the same tree on absolute values; m >= |v| by monotone rounding) and d.
// Claim: if no product's v or m falls below the normal range and nothing overflows, then
// |v - exact| <= gamma_d M, with M the exact absolute-value tree and gamma_k = k u / (1 - k u).
// Induction:
// - a leaf (an input or a constant 0, +-1): v exact, d = 0;
// - a sum or difference z = x +- y: fl(x' +- y') = (x' +- y')(1 + delta), |delta| <= u (delta = 0
//   when the result is subnormal: such a sum is exact under gradual underflow), so
//   |z' - z| <= gamma_a Mx + gamma_b My + u (1 + gamma_max)(Mx + My) <= gamma_(max(a,b)+1) Mz:
//   d = max(a, b) + 1. Cancellation is covered: the bound is relative to M, not to |z|;
// - a product: |x'y' - xy| <= (gamma_a + gamma_b + gamma_a gamma_b) Mx My, and the rounding adds
//   u (1 + gamma_a)(1 + gamma_b) Mx My; in all (1 + gamma_a)(1 + gamma_b)(1 + u) - 1 <=
//   gamma_(a+b+1): d = a + b + 1. A product below the normal range has an absolute error up to
//   2^-1075 that this relative model does not cover, so it marks the evaluation `underflow`;
// - unary minus is exact.
// The same induction on the absolute-value tree gives m >= M (1 - gamma_d), so gamma_d M <=
// gamma_d m / (1 - gamma_d) <= d u m (1 + 3 d u) for d u <= 2^-40 (d <= 19 here). The bound
// 2 d u m is computed with one rounding. In the normal range that costs u relative. A subnormal
// bound needs no special case: every formula's root is a sum of products, and a product's m is 0 or
// at least DBL_MIN (else it is flagged), so a nonzero m is at least DBL_MIN, the exact bound at
// least 2 d u DBL_MIN = d 2^-1073, and its rounding error (at most 2^-1075) at most 1 / (4 d) <=
// 1/16 relative; the factor 2 still covers d u m (1 + 3 d u). Overflow: an infinity or NaN anywhere
// propagates to v or m (inf - inf and inf * 0 give NaN), so finite v and m mean no overflow. Then
// |v| > bound implies sign(v) = sign(exact).
// Which underflow checks are load-bearing: the value check catches a nonzero product that rounds
// below the normal range or to zero (a factor's value may be small by cancellation while its
// magnitude is not); the magnitude check catches m below the normal range when a factor's value
// is exactly zero, which keeps a product's m at 0 or at least DBL_MIN, as the subnormal-bound
// argument above needs. The depths of the formulas (pinned by tests): orient2d 4, orient3d 8, incircle 11,
// insphere 17, diametral circle 4, diametral sphere 19, incircle cofactors 7, insphere cofactors 12.
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
#include <cstdint>
#include <span>

#include "NiTCAD/mesh/mesh.hpp"
#include "exact_integer.hpp"

namespace NiTCAD::mesh::predicates {

[[nodiscard]] int orient2d(const Point& a, const Point& b, const Point& c);
[[nodiscard]] int orient3d(const Point& a, const Point& b, const Point& c, const Point& d);
[[nodiscard]] int incircle(const Point& a, const Point& b, const Point& c, const Point& d);
[[nodiscard]] int insphere(const Point& a, const Point& b, const Point& c, const Point& d,
                           const Point& e);
[[nodiscard]] int in_diametral_circle(const Point& a, const Point& b, const Point& p);
// 2D angle tests at o between the rays to a and b (both distinct from o): the sign of
// (a - o) . (b - o) (0 exactly at 90 degrees), and whether the angle (in [0, 180]) is below 60
// degrees, cos > 1/2, decided exactly as (a - o).(b - o) > 0 and 4 ((a - o).(b - o))^2 >
// |a - o|^2 |b - o|^2.
[[nodiscard]] int dot_sign(const Point& o, const Point& a, const Point& b);
[[nodiscard]] bool angle_below_60(const Point& o, const Point& a, const Point& b);
// The same in 3D: the sign of (a - o) . (b - o), and whether the angle at o between the rays to a
// and b is below 60 degrees.
[[nodiscard]] int dot_sign_3d(const Point& o, const Point& a, const Point& b);
// The sign of the cosine of the dihedral angle between the half-planes bounded by the line ab and
// containing p and q (both off the line): + below 90 degrees, 0 at 90, - above.
[[nodiscard]] int dihedral_dot_sign(const Point& a, const Point& b, const Point& p, const Point& q);
[[nodiscard]] bool angle_below_60_3d(const Point& o, const Point& a, const Point& b);
// Whether the dihedral angle (in [0, 180]) between the half-planes bounded by the line ab and
// containing p and q (both off the line) is below 60 degrees, exactly (Unit 53, FD7).
[[nodiscard]] bool dihedral_below_60(const Point& a, const Point& b, const Point& p, const Point& q);
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

// The filter's internals for one evaluation, with the exact value beside it: for tests of the
// soundness argument above (never used by the predicates themselves).
namespace probe {

enum class Kind : std::uint8_t {
    orient2d,           // 3 points
    orient3d,           // 4 points
    incircle,           // 4 points
    insphere,           // 5 points
    diametral_circle,   // 3 points
    diametral_sphere,   // 4 points
    incircle_cofactor,  // 4 points; column 0-3: the lifted-value cofactor of that point
    insphere_cofactor,  // 5 points; column 0-4
};

struct Result {
    double value = 0.0;      // v
    double magnitude = 0.0;  // m
    int depth = 0;           // d
    bool underflow = false;
    double bound = 0.0;      // 2 d u m (0 when the filter could not form one)
    int filter_sign = 0;     // the filter's decision; 0 when it defers
    int exact_sign = 0;
    detail::ExactInteger exact;  // the exact value is exact * 2^exponent
    int exponent = 0;
};

// Precondition (NITCAD_EXPECTS): the point count and column fit the kind.
[[nodiscard]] Result evaluate(Kind kind, std::span<const Point> points, int column = 0);

}  // namespace probe

}  // namespace NiTCAD::mesh::predicates
