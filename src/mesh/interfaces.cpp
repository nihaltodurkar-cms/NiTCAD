#include "NiTCAD/mesh/interfaces.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "dual_geometry.hpp"
#include "polygon2d.hpp"
#include "predicates.hpp"

namespace NiTCAD::mesh {

namespace {

// P4a: deviation tolerance outside the corner zones, relative to the edge and to the coordinates'
// magnitude (interfaces.hpp). Zone radius in units of the resolution.
constexpr double flat_relative = 1e-12;
constexpr double flat_rounding = 1e-13;
constexpr double zone_radius = 2.0;

using Vec3 = std::array<double, 3>;

Vec3 minus(const Point& a, const Point& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::hypot(a[0], a[1], a[2]); }

double segment_distance(const Point& x, const Point& a, const Point& b) {
    const Vec3 ab = minus(b, a);
    const double len2 = dot(ab, ab);
    const double t = len2 > 0.0 ? std::clamp(dot(minus(x, a), ab) / len2, 0.0, 1.0) : 0.0;
    const Point foot{a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]};
    return detail::distance(x, foot);
}

// A piece prepared for distance queries: a segment, or a convex polygon with its unit normal, or
// (Unit 53, FD9) a planar simple non-convex polygon with its convex decomposition.
struct Piece {
    std::vector<Point> v;
    Vec3 normal{0.0, 0.0, 0.0};
    std::int32_t a = 0;
    std::int32_t b = 0;
    // Non-convex polygons only (empty `parts` for segments and convex polygons).
    std::vector<Piece> parts;  // the ear-clipped triangles, each a convex piece
    int axis = -1;             // the coordinate constant on every vertex, or -1 (oblique)
    std::vector<Point> flat;   // axis >= 0: the vertices with that coordinate dropped (exact)
    std::vector<Point> grown;  // axis >= 0: `flat` grown outward by a mitred offset (rounded coordinates)
    double grown_reach = 0.0;  // the largest distance from a vertex of `grown` to the polygon

    // FD9: max over V of the distance to the polygon. For an axis-aligned polygon whose exact 2D
    // projection contains the convex hull of V's projection, the orthogonal projection of every
    // point of that hull lies in the polygon, so the nearest point of the polygon to each is on
    // the plane and the distance is the plane distance, linear in the point: its max over the
    // hull is attained at a point of V. When the hull lies only in the grown polygon, every point
    // of it is within grown_reach of the polygon in the plane (the grown polygon's farthest
    // points from it are its vertices), so the distance is at most the plane distance plus
    // grown_reach: this covers face vertices rounded just outside a feature line, where the exact
    // test fails by an ulp (they lie in corner zones, judged against s). Otherwise the
    // decomposition's min over parts of max over V (an upper bound, as for several pieces of one
    // pair). Oblique non-convex polygons always take the decomposition: their orthogonal
    // projection is not exact.
    [[nodiscard]] double nonconvex_deviation(const std::vector<Point>& V) const {
        if (axis >= 0) {
            const auto ax = static_cast<std::size_t>(axis);
            double worst = 0.0;
            for (const Point& x : V) worst = std::max(worst, std::abs(x[ax] - v[0][ax]));
            if (hull_inside(V, flat)) return worst;
            if (!grown.empty() && hull_inside(V, grown)) return worst + grown_reach;
        }
        double best = std::numeric_limits<double>::infinity();
        for (const Piece& t : parts) {
            double worst = 0.0;
            for (const Point& x : V) worst = std::max(worst, t.distance(x));
            best = std::min(best, worst);
        }
        return best;
    }
    // Whether the convex hull of V's projection lies in the closed projected polygon `poly`,
    // exactly: every hull vertex in it and no polygon edge meeting the hull's open interior (if
    // the hull had an interior point outside, the segment from a hull vertex to it would cross
    // the polygon's boundary inside the hull). A hull with no interior is not judged (false).
    [[nodiscard]] bool hull_inside(const std::vector<Point>& V, const std::vector<Point>& poly) const {
        std::vector<Point> q;
        for (const Point& x : V) q.push_back(detail::project(x, axis));
        std::ranges::sort(q, [](const Point& x, const Point& y) { return x[0] < y[0] || (x[0] == y[0] && x[1] < y[1]); });
        q.erase(std::unique(q.begin(), q.end()), q.end());
        if (q.size() < 3) return false;
        // Andrew's monotone chain, counterclockwise, collinear points dropped (exact orient2d).
        std::vector<Point> h(2 * q.size());
        std::size_t k = 0;
        for (std::size_t i = 0; i < q.size(); ++i) {
            while (k >= 2 && predicates::orient2d(h[k - 2], h[k - 1], q[i]) <= 0) --k;
            h[k++] = q[i];
        }
        for (std::size_t i = q.size() - 1, t = k + 1; i-- > 0;) {
            while (k >= t && predicates::orient2d(h[k - 2], h[k - 1], q[i]) <= 0) --k;
            h[k++] = q[i];
        }
        h.resize(k - 1);
        if (h.size() < 3) return false;
        for (const Point& x : h) {
            if (!detail::inside_closed_2d(poly, x)) return false;
        }
        const std::size_t n = poly.size();
        for (std::size_t i = 0; i < n; ++i) {
            const Point& p = poly[i];
            const Point& r = poly[(i + 1) % n];
            // Separated from the open hull by a hull edge's line, or by the edge's own line.
            bool separated = false;
            for (std::size_t j = 0; j < h.size() && !separated; ++j) {
                const Point& c = h[j];
                const Point& d = h[(j + 1) % h.size()];
                separated = predicates::orient2d(c, d, p) <= 0 && predicates::orient2d(c, d, r) <= 0;
            }
            if (!separated) {
                bool left = true;
                bool right = true;
                for (const Point& c : h) {
                    const int o = predicates::orient2d(p, r, c);
                    left = left && o >= 0;
                    right = right && o <= 0;
                }
                separated = left || right;
            }
            if (!separated) return false;
        }
        return true;
    }

    [[nodiscard]] bool separates(std::int32_t x, std::int32_t y) const noexcept {
        return (a == x && b == y) || (a == y && b == x);
    }
    [[nodiscard]] double distance(const Point& x) const {
        if (v.size() == 2) return segment_distance(x, v[0], v[1]);
        const double h = dot(minus(x, v[0]), normal);
        const Point foot{x[0] - h * normal[0], x[1] - h * normal[1], x[2] - h * normal[2]};
        bool inside = true;
        for (std::size_t k = 0; k < v.size() && inside; ++k) {
            const Point& p = v[k];
            const Point& q = v[(k + 1) % v.size()];
            inside = dot(cross(minus(q, p), minus(foot, p)), normal) >= 0.0;
        }
        if (inside) return std::abs(h);
        return boundary_distance(x);
    }
    // Distance to the relative boundary: the end points of a segment, the edges of a polygon.
    [[nodiscard]] double boundary_distance(const Point& x) const {
        if (v.size() == 2) return std::min(detail::distance(x, v[0]), detail::distance(x, v[1]));
        double d = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < v.size(); ++k) d = std::min(d, segment_distance(x, v[k], v[(k + 1) % v.size()]));
        return d;
    }
};

base::Error error(base::ErrorCode code, std::string message, std::optional<std::size_t> index,
                  std::optional<double> value = std::nullopt) {
    return {code, std::move(message), base::ErrorContext{.index = index, .value = value}};
}

std::string point_text(const Point& p) {
    char at[96];
    std::snprintf(at, sizeof at, " at (%.6g, %.6g, %.6g) cm", p[0], p[1], p[2]);
    return at;
}

std::string pair_text(NodeId a, NodeId b) {
    return "(" + std::to_string(a) + ", " + std::to_string(b) + ")";
}

}  // namespace

std::expected<InterfaceReport, base::Error> check_interfaces(const SimplexMesh& s, const Mesh& dual,
                                                             std::span<const std::int32_t> node_region,
                                                             const InterfaceSet& interfaces) {
    using base::ErrorCode;
    const int dim = s.dimension;
    const auto D = static_cast<std::size_t>(dim);
    const std::size_t per_cell = D + 1;
    if (dual.dimension() != dim || dual.node_count() != s.points.size() ||
        !std::ranges::equal(dual.points(), s.points) || (dim != 2 && dim != 3) ||
        s.cells.size() % per_cell != 0) {
        return std::unexpected(error(ErrorCode::invalid_input,
                                     "the dual mesh is not built from these simplices", std::nullopt));
    }
    const std::size_t n = s.points.size();
    if (node_region.size() != n) {
        return std::unexpected(error(ErrorCode::invalid_input, "node_region needs one region per node",
                                     std::nullopt));
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (node_region[i] < 0) {
            return std::unexpected(error(ErrorCode::invalid_input, "node region is negative", i));
        }
    }
    const double res = interfaces.resolution;
    if (!(std::isfinite(res) && res > 0.0)) {
        return std::unexpected(error(ErrorCode::invalid_input,
                                     "interface resolution must be finite and positive",
                                     std::nullopt, res));
    }

    // Pieces.
    std::vector<Piece> pieces;
    for (std::size_t k = 0; k < interfaces.pieces.size(); ++k) {
        const InterfacePiece& in = interfaces.pieces[k];
        const bool count_ok = dim == 2 ? in.vertices.size() == 2 : in.vertices.size() >= 3;
        bool finite = true;
        for (const Point& p : in.vertices) {
            for (std::size_t c = 0; c < 3; ++c) finite = finite && std::isfinite(p[c]) && (c < D || p[c] == 0.0);
        }
        if (!count_ok || !finite || in.region_a < 0 || in.region_b < 0 || in.region_a == in.region_b) {
            return std::unexpected(error(ErrorCode::invalid_input,
                                         "interface piece: wrong vertex count, a vertex not finite, or "
                                         "regions equal or negative",
                                         k));
        }
        Piece p{in.vertices, {0.0, 0.0, 0.0}, in.region_a, in.region_b};
        if (dim == 2) {
            if (!(detail::distance(p.v[0], p.v[1]) > 0.0)) {
                return std::unexpected(error(ErrorCode::invalid_input, "interface segment is degenerate", k));
            }
        } else {
            // Newell's normal; planarity and convexity relative to the polygon's extent.
            Vec3 normal{0.0, 0.0, 0.0};
            double extent = 0.0;
            for (std::size_t i = 0; i < p.v.size(); ++i) {
                const Point& a = p.v[i];
                const Point& b = p.v[(i + 1) % p.v.size()];
                normal[0] += (a[1] - b[1]) * (a[2] + b[2]);
                normal[1] += (a[2] - b[2]) * (a[0] + b[0]);
                normal[2] += (a[0] - b[0]) * (a[1] + b[1]);
                extent = std::max(extent, detail::distance(a, p.v[0]));
            }
            const double len = norm(normal);
            bool ok = len > 0.0 && extent > 0.0;
            bool planar = ok;
            if (ok) {
                for (double& c : normal) c /= len;
                for (const Point& a : p.v) ok = ok && std::abs(dot(minus(a, p.v[0]), normal)) <= 1e-12 * extent;
                planar = ok;
                for (std::size_t i = 0; i < p.v.size() && ok; ++i) {
                    const Point& a = p.v[i];
                    const Point& b = p.v[(i + 1) % p.v.size()];
                    const Point& c = p.v[(i + 2) % p.v.size()];
                    ok = dot(cross(minus(b, a), minus(c, b)), normal) > 0.0;
                }
            }
            p.normal = normal;
            if (!ok && planar) {
                // FD9: not convex by the rule above (which keeps every polygon it accepts on the
                // original path), but planar and simple in the exact projection dropping the
                // normal's largest component; decomposed by ear clipping.
                int drop = 0;
                for (int c = 1; c < 3; ++c) {
                    if (std::abs(normal[static_cast<std::size_t>(c)]) > std::abs(normal[static_cast<std::size_t>(drop)])) drop = c;
                }
                std::vector<Point> q;
                for (const Point& a : p.v) q.push_back(detail::project(a, drop));
                std::optional<std::vector<std::array<std::size_t, 3>>> ears;
                if (detail::simple_2d(q)) ears = detail::ear_clip(q, detail::polygon_turn(q));
                if (ears) {
                    ok = true;
                    for (const auto& t : *ears) {
                        Piece part{{p.v[t[0]], p.v[t[1]], p.v[t[2]]}, {0.0, 0.0, 0.0}, p.a, p.b};
                        const Vec3 tn = cross(minus(part.v[1], part.v[0]), minus(part.v[2], part.v[0]));
                        const double sign = dot(tn, normal) < 0.0 ? -1.0 : 1.0;
                        const double tl = norm(tn);
                        ok = ok && tl > 0.0;
                        if (!ok) break;
                        // Oriented with the polygon's normal (distance() needs counterclockwise).
                        if (sign < 0.0) std::swap(part.v[1], part.v[2]);
                        part.normal = normal;
                        p.parts.push_back(std::move(part));
                    }
                    for (int c = 0; c < 3 && ok; ++c) {
                        const auto cc = static_cast<std::size_t>(c);
                        if (std::ranges::all_of(p.v, [&](const Point& a) { return a[cc] == p.v[0][cc]; })) {
                            p.axis = c;
                            for (const Point& a : p.v) p.flat.push_back(detail::project(a, c));
                        }
                    }
                    if (p.axis >= 0) {
                        // Grown outward by 1e-6 s: each edge's line moved out, consecutive lines
                        // met (a collinear vertex moved along its edge's normal).
                        const double grow = 1e-6 * res;
                        const int turn = detail::polygon_turn(p.flat);
                        const std::size_t m = p.flat.size();
                        const auto outward = [&](const Point& a, const Point& b) {
                            const double dx = b[0] - a[0];
                            const double dy = b[1] - a[1];
                            const double len = std::hypot(dx, dy);
                            return turn > 0 ? Point{dy / len, -dx / len, 0.0} : Point{-dy / len, dx / len, 0.0};
                        };
                        for (std::size_t i = 0; i < m; ++i) {
                            const Point& pa = p.flat[(i + m - 1) % m];
                            const Point& pb = p.flat[i];
                            const Point& pc = p.flat[(i + 1) % m];
                            const Point n1 = outward(pa, pb);
                            const Point n2 = outward(pb, pc);
                            const Point d1{pb[0] - pa[0], pb[1] - pa[1], 0.0};
                            const Point d2{pc[0] - pb[0], pc[1] - pb[1], 0.0};
                            const double den = d1[0] * d2[1] - d1[1] * d2[0];
                            Point g{pb[0] + grow * n2[0], pb[1] + grow * n2[1], 0.0};
                            if (predicates::orient2d(pa, pb, pc) != 0) {
                                const Point q1{pa[0] + grow * n1[0], pa[1] + grow * n1[1], 0.0};
                                const Point q2{pb[0] + grow * n2[0], pb[1] + grow * n2[1], 0.0};
                                const double t = ((q2[0] - q1[0]) * d2[1] - (q2[1] - q1[1]) * d2[0]) / den;
                                g = {q1[0] + t * d1[0], q1[1] + t * d1[1], 0.0};
                            }
                            p.grown.push_back(g);
                            p.grown_reach = std::max(p.grown_reach, detail::distance(g, pb));
                        }
                        if (!detail::simple_2d(p.grown)) p.grown.clear();  // then the path is not taken
                    }
                }
            }
            if (!ok) {
                return std::unexpected(error(ErrorCode::invalid_input,
                                             "interface polygon is not planar, simple and non-degenerate", k));
            }
        }
        pieces.push_back(std::move(p));
    }

    // The kept dual edges.
    std::vector<std::pair<NodeId, NodeId>> kept;
    kept.reserve(dual.edges().size());
    for (const Edge& e : dual.edges()) kept.emplace_back(e.first, e.second);
    std::ranges::sort(kept);

    // Boundary facets (faces of exactly one cell), as sorted keys.
    const std::size_t cell_count = s.cells.size() / per_cell;
    const auto node = [&](std::size_t c, std::size_t k) { return s.cells[c * per_cell + k]; };
    std::vector<std::array<NodeId, 3>> faces;
    for (std::size_t c = 0; c < cell_count; ++c) {
        for (std::size_t l = 0; l < per_cell; ++l) {
            std::array<NodeId, 3> key{-1, -1, -1};
            for (std::size_t q = 0, k = 0; q < per_cell; ++q) {
                if (q != l) key[k++] = node(c, q);
            }
            std::sort(key.begin(), key.begin() + static_cast<std::ptrdiff_t>(D));
            faces.push_back(key);
        }
    }
    std::ranges::sort(faces);
    std::vector<std::array<NodeId, 3>> boundary;
    for (std::size_t g = 0; g < faces.size();) {
        std::size_t e = g + 1;
        while (e < faces.size() && faces[e] == faces[g]) ++e;
        if (e - g == 1) boundary.push_back(faces[g]);
        g = e;
    }

    // Each boundary facet's edges (in 2D the facet is the edge), sorted, for the face vertices.
    struct BoundaryEdge {
        std::pair<NodeId, NodeId> edge;
        std::size_t facet;
    };
    std::vector<BoundaryEdge> boundary_edges;
    for (std::size_t k = 0; k < boundary.size(); ++k) {
        for (std::size_t i = 0; i < D; ++i) {
            for (std::size_t j = i + 1; j < D; ++j) boundary_edges.push_back({{boundary[k][i], boundary[k][j]}, k});
        }
    }
    std::ranges::sort(boundary_edges, [](const BoundaryEdge& x, const BoundaryEdge& y) {
        return std::tie(x.edge, x.facet) < std::tie(y.edge, y.facet);
    });

    // Cross-region edges with the cells around each.
    struct Incidence {
        NodeId a;
        NodeId b;
        std::size_t cell;
    };
    std::vector<Incidence> cross_edges;
    for (std::size_t c = 0; c < cell_count; ++c) {
        for (std::size_t i = 0; i < per_cell; ++i) {
            for (std::size_t j = i + 1; j < per_cell; ++j) {
                const NodeId a = std::min(node(c, i), node(c, j));
                const NodeId b = std::max(node(c, i), node(c, j));
                if (node_region[static_cast<std::size_t>(a)] != node_region[static_cast<std::size_t>(b)]) {
                    cross_edges.push_back({a, b, c});
                }
            }
        }
    }
    std::ranges::sort(cross_edges, [](const Incidence& x, const Incidence& y) {
        return std::tie(x.a, x.b, x.cell) < std::tie(y.a, y.b, y.cell);
    });

    const auto P = [&](NodeId v) -> const Point& { return s.points[static_cast<std::size_t>(v)]; };
    InterfaceReport report;
    for (std::size_t g = 0; g < cross_edges.size();) {
        const NodeId a = cross_edges[g].a;
        const NodeId b = cross_edges[g].b;
        std::vector<Point> vertices;
        for (; g < cross_edges.size() && cross_edges[g].a == a && cross_edges[g].b == b; ++g) {
            const std::size_t c = cross_edges[g].cell;
            if (dim == 2) {
                vertices.push_back(detail::triangle_circumcentre({P(node(c, 0)), P(node(c, 1)), P(node(c, 2))}));
            } else {
                vertices.push_back(
                    detail::tetrahedron_circumcentre({P(node(c, 0)), P(node(c, 1)), P(node(c, 2)), P(node(c, 3))}));
            }
        }
        if (!std::ranges::binary_search(kept, std::pair{a, b})) continue;  // a zero coupling: no face

        const Point& pa = P(a);
        const Point& pb = P(b);
        const Point mid{0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])};
        const auto [lo, hi] = std::ranges::equal_range(boundary_edges, std::pair{a, b}, {},
                                                       &BoundaryEdge::edge);
        for (auto it = lo; it != hi; ++it) {
            const auto& f = boundary[it->facet];
            if (dim == 3) vertices.push_back(detail::triangle_circumcentre({P(f[0]), P(f[1]), P(f[2])}));
        }
        if (lo != hi) vertices.push_back(mid);

        const std::int32_t ra = node_region[static_cast<std::size_t>(a)];
        const std::int32_t rb = node_region[static_cast<std::size_t>(b)];
        double deviation = std::numeric_limits<double>::infinity();
        bool declared = false;
        for (const Piece& p : pieces) {
            if (!p.separates(ra, rb)) continue;
            declared = true;
            if (!p.parts.empty()) {
                deviation = std::min(deviation, p.nonconvex_deviation(vertices));
                continue;
            }
            double worst = 0.0;
            for (const Point& v : vertices) worst = std::max(worst, p.distance(v));
            deviation = std::min(deviation, worst);
        }
        if (!declared) {
            return std::unexpected(error(ErrorCode::degenerate_mesh,
                                         "interface rule (P4): regions " + std::to_string(ra) + " and " +
                                             std::to_string(rb) + " meet at edge " + pair_text(a, b) +
                                             " with no declared interface",
                                         static_cast<std::size_t>(a)));
        }
        ++report.cross_faces;
        bool in_zone = false;
        for (const Piece& p : pieces) {
            in_zone = in_zone || p.boundary_distance(mid) < zone_radius * res;
            for (const Point& v : vertices) in_zone = in_zone || p.boundary_distance(v) < zone_radius * res;
        }
        if (in_zone) {
            ++report.corner_faces;
            report.max_corner_deviation = std::max(report.max_corner_deviation, deviation);
            if (!(deviation <= res)) {
                return std::unexpected(error(ErrorCode::degenerate_mesh,
                                             "interface rule (P4b): the dual face of edge " + pair_text(a, b) +
                                                 " in a corner zone deviates beyond the resolution" + point_text(mid),
                                             static_cast<std::size_t>(a), deviation));
            }
        } else {
            double magnitude = 0.0;
            for (const Point& v : vertices) magnitude = std::max({magnitude, std::abs(v[0]), std::abs(v[1]), std::abs(v[2])});
            for (const Point& v : {pa, pb}) magnitude = std::max({magnitude, std::abs(v[0]), std::abs(v[1]), std::abs(v[2])});
            const double tolerance = flat_relative * detail::distance(pa, pb) + flat_rounding * magnitude;
            report.max_flat_deviation = std::max(report.max_flat_deviation, deviation / tolerance);
            if (!(deviation <= tolerance)) {
                return std::unexpected(error(ErrorCode::degenerate_mesh,
                                             "interface rule (P4a): the dual face of edge " + pair_text(a, b) +
                                                 " is not on the interface",
                                             static_cast<std::size_t>(a), deviation));
            }
        }
    }
    return report;
}

namespace {

void digest_word(std::uint64_t& h, std::uint64_t v) {
    for (int i = 0; i < 8; ++i) {
        h ^= (v >> (8 * i)) & 0xffU;
        h *= 0x100000001b3ULL;
    }
}

void digest_point(std::uint64_t& h, const Point& p) {
    for (const double c : p) digest_word(h, std::bit_cast<std::uint64_t>(c));
}

// Without a patch set, the words of Unit 16 (existing digests unchanged).
std::uint64_t interfaces_digest(const InterfaceSet& set, const PatchBoundarySet* patches) {
    std::uint64_t h = 0xcbf29ce484222325ULL;
    digest_word(h, std::bit_cast<std::uint64_t>(set.resolution));
    digest_word(h, set.pieces.size());
    for (const InterfacePiece& p : set.pieces) {
        digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.region_a)));
        digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.region_b)));
        digest_word(h, p.vertices.size());
        for (const Point& v : p.vertices) digest_point(h, v);
    }
    if (patches != nullptr) {
        digest_word(h, std::bit_cast<std::uint64_t>(patches->resolution));
        digest_word(h, patches->boundaries.size());
        for (const PatchBoundary& b : patches->boundaries) {
            digest_point(h, b.a);
            digest_point(h, b.b);
            digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(b.patch_a)));
            digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(b.patch_b)));
        }
        digest_word(h, patches->zone_centres.size());
        for (const Point& z : patches->zone_centres) digest_point(h, z);
        digest_word(h, patches->patch_area.size());
        for (const double a : patches->patch_area) digest_word(h, std::bit_cast<std::uint64_t>(a));
    }
    return h;
}

}  // namespace

std::expected<PatchBoundaryReport, base::Error> check_patch_boundaries(const SimplexMesh& s, const Mesh& dual,
                                                                       const PatchBoundarySet& patches) {
    using base::ErrorCode;
    if (s.dimension != 3 || dual.dimension() != 3 || dual.node_count() != s.points.size() ||
        !std::ranges::equal(dual.points(), s.points) || s.facets.size() % 3 != 0 ||
        s.facet_patch.size() != s.facets.size() / 3 ||
        (!s.facet_corner_patch.empty() && s.facet_corner_patch.size() != s.facets.size())) {
        return std::unexpected(error(ErrorCode::invalid_input, "the dual mesh is not built from these 3D simplices", std::nullopt));
    }
    const double res = patches.resolution;
    if (!(std::isfinite(res) && res > 0.0)) {
        return std::unexpected(error(ErrorCode::invalid_input, "patch resolution must be finite and positive", std::nullopt, res));
    }
    const std::size_t np = s.patch_names.size();
    if (patches.patch_area.size() != np ||
        !std::ranges::all_of(patches.patch_area, [](double a) { return std::isfinite(a) && a >= 0.0; })) {
        return std::unexpected(error(ErrorCode::invalid_input, "patch_area needs one finite non-negative area per patch", std::nullopt));
    }
    const auto finite = [](const Point& p) { return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]); };
    for (std::size_t k = 0; k < patches.boundaries.size(); ++k) {
        const PatchBoundary& b = patches.boundaries[k];
        if (!finite(b.a) || !finite(b.b) || !(detail::distance(b.a, b.b) > 0.0) || b.patch_a < 0 || b.patch_b < 0 ||
            static_cast<std::size_t>(b.patch_a) >= np || static_cast<std::size_t>(b.patch_b) >= np || b.patch_a == b.patch_b) {
            return std::unexpected(error(ErrorCode::invalid_input,
                                         "patch boundary: an end not finite, degenerate, or patches out of range or equal", k));
        }
    }
    for (std::size_t k = 0; k < patches.zone_centres.size(); ++k) {
        if (!finite(patches.zone_centres[k])) {
            return std::unexpected(error(ErrorCode::invalid_input, "patch zone centre not finite", k));
        }
    }
    const double zone = zone_radius * res;
    const auto in_zone = [&](const Point& x) {
        return std::ranges::any_of(patches.zone_centres, [&](const Point& z) { return detail::distance(x, z) <= zone; });
    };
    const auto magnitude = [](const Point& p) { return std::max({std::abs(p[0]), std::abs(p[1]), std::abs(p[2])}); };

    PatchBoundaryReport report;
    const std::size_t nf = s.facets.size() / 3;
    // Where a failing piece is: its edge's midpoint and the two patches (for the message).
    const auto where = [&](const Point& m, std::int32_t pa, std::int32_t pb) {
        return point_text(m) + " between '" + s.patch_names[static_cast<std::size_t>(pa)] + "' and '" +
               s.patch_names[static_cast<std::size_t>(pb)] + "'";
    };
    for (std::size_t f = 0; f < nf; ++f) {
        std::array<std::int32_t, 3> cp{};
        for (std::size_t k = 0; k < 3; ++k) cp[k] = s.facet_corner_patch.empty() ? s.facet_patch[f] : s.facet_corner_patch[3 * f + k];
        if (cp[0] == cp[1] && cp[1] == cp[2]) continue;
        const std::array<Point, 3> q{s.points[static_cast<std::size_t>(s.facets[3 * f])],
                                     s.points[static_cast<std::size_t>(s.facets[3 * f + 1])],
                                     s.points[static_cast<std::size_t>(s.facets[3 * f + 2])]};
        const Point c = detail::triangle_circumcentre(q);
        for (std::size_t i = 0; i < 3; ++i) {
            const std::size_t j = (i + 1) % 3;
            if (cp[i] == cp[j]) continue;
            const Point m{0.5 * (q[i][0] + q[j][0]), 0.5 * (q[i][1] + q[j][1]), 0.5 * (q[i][2] + q[j][2])};
            double deviation = std::numeric_limits<double>::infinity();
            for (const PatchBoundary& b : patches.boundaries) {
                if (!((b.patch_a == cp[i] && b.patch_b == cp[j]) || (b.patch_a == cp[j] && b.patch_b == cp[i]))) continue;
                deviation = std::min(deviation, std::max(segment_distance(m, b.a, b.b), segment_distance(c, b.a, b.b)));
            }
            if (!std::isfinite(deviation)) {
                return std::unexpected(error(ErrorCode::degenerate_mesh,
                                             "patches '" + s.patch_names[static_cast<std::size_t>(cp[i])] + "' and '" +
                                                 s.patch_names[static_cast<std::size_t>(cp[j])] +
                                                 "' meet with no declared boundary between them (P6)",
                                             f));
            }
            ++report.cross_pieces;
            if (in_zone(m) || in_zone(c)) {
                ++report.corner_pieces;
                report.max_corner_deviation = std::max(report.max_corner_deviation, deviation);
                if (!(deviation <= res)) {
                    return std::unexpected(error(ErrorCode::degenerate_mesh,
                                                 "a boundary dual piece in a patch zone deviates beyond s from its patch boundary (P6)" + where(m, cp[i], cp[j]),
                                                 f, deviation));
                }
                continue;
            }
            const double tolerance = flat_relative * detail::distance(q[i], q[j]) +
                                     flat_rounding * std::max({magnitude(m), magnitude(c), magnitude(q[i]), magnitude(q[j])});
            report.max_flat_deviation = std::max(report.max_flat_deviation, deviation / tolerance);
            if (!(deviation <= tolerance)) {
                return std::unexpected(error(ErrorCode::degenerate_mesh,
                                             "a boundary dual piece deviates from its patch boundary (P6)" + where(m, cp[i], cp[j]), f, deviation));
            }
        }
    }

    // Areas, and the declared boundary inside the zones (each boundary sampled finely enough
    // for a report: 64 steps per zone radius).
    report.area_error.assign(np, 0.0);
    report.zone_length.assign(np, 0.0);
    for (std::size_t p = 0; p < np; ++p) {
        double area = 0.0;
        if (const BoundaryPatch* bp = dual.find_boundary(s.patch_names[p])) {
            for (const double a : bp->areas) area += a;
        }
        report.area_error[p] = area - patches.patch_area[p];
    }
    for (const PatchBoundary& b : patches.boundaries) {
        const double len = detail::distance(b.a, b.b);
        const auto steps = static_cast<std::size_t>(std::ceil(64.0 * len / zone));
        double inside = 0.0;
        for (std::size_t k = 0; k < steps; ++k) {
            const double t = (static_cast<double>(k) + 0.5) / static_cast<double>(steps);
            const Point x{b.a[0] + t * (b.b[0] - b.a[0]), b.a[1] + t * (b.b[1] - b.a[1]), b.a[2] + t * (b.b[2] - b.a[2])};
            if (in_zone(x)) inside += len / static_cast<double>(steps);
        }
        report.zone_length[static_cast<std::size_t>(b.patch_a)] += inside;
        report.zone_length[static_cast<std::size_t>(b.patch_b)] += inside;
    }
    return report;
}

std::uint64_t region_digest(std::span<const std::int32_t> node_region) {
    std::uint64_t h = 0xcbf29ce484222325ULL;
    digest_word(h, node_region.size());
    for (const std::int32_t r : node_region) {
        digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(r)));
    }
    return h;
}

std::expected<InterfaceCertificate, base::Error> certify_interfaces(
    const SimplexMesh& simplices, const Mesh& dual, std::span<const std::int32_t> node_region,
    const InterfaceSet& interfaces, const PatchBoundarySet* patches) {
    if (!simplices.facet_corner_patch.empty() && patches == nullptr) {
        return std::unexpected(error(base::ErrorCode::invalid_input,
                                     "simplices with per-corner patches need a patch-boundary set (P6, FD10)", std::nullopt));
    }
    auto report = check_interfaces(simplices, dual, node_region, interfaces);
    if (!report) return std::unexpected(std::move(report.error()));
    std::optional<PatchBoundaryReport> p6;
    if (patches != nullptr) {
        auto checked = check_patch_boundaries(simplices, dual, *patches);
        if (!checked) return std::unexpected(std::move(checked.error()));
        p6 = std::move(*checked);
    }
    return InterfaceCertificate{.mesh_digest = graph_digest(dual),
                                .region_digest = region_digest(node_region),
                                .interfaces_digest = interfaces_digest(interfaces, patches),
                                .report = *report,
                                .patches = std::move(p6)};
}

bool certifies(const InterfaceCertificate& c, const Mesh& mesh,
               std::span<const std::int32_t> node_region) {
    return c.mesh_digest == graph_digest(mesh) && c.region_digest == region_digest(node_region);
}


std::expected<InterfaceCertificate, base::Error> certify_interfaces_1d(
    const Mesh& mesh, std::span<const std::int32_t> node_region,
    std::span<const InterfacePoint1D> interfaces) {
    using base::ErrorCode;
    if (mesh.dimension() != 1) {
        return std::unexpected(error(ErrorCode::invalid_input, "certify_interfaces_1d needs a 1D mesh", std::nullopt));
    }
    if (node_region.size() != mesh.node_count()) {
        return std::unexpected(error(ErrorCode::invalid_input, "node_region needs one region per node", std::nullopt));
    }
    for (std::size_t i = 0; i < node_region.size(); ++i) {
        if (node_region[i] < 0) return std::unexpected(error(ErrorCode::invalid_input, "node region is negative", i));
    }
    for (std::size_t k = 0; k < interfaces.size(); ++k) {
        const InterfacePoint1D& p = interfaces[k];
        if (!std::isfinite(p.x) || p.region_a < 0 || p.region_b < 0 || p.region_a == p.region_b) {
            return std::unexpected(error(ErrorCode::invalid_input,
                                         "interface point not finite, or regions equal or negative", k));
        }
    }
    InterfaceReport report;
    for (const Edge& e : mesh.edges()) {
        const std::int32_t ra = node_region[static_cast<std::size_t>(e.first)];
        const std::int32_t rb = node_region[static_cast<std::size_t>(e.second)];
        if (ra == rb) continue;
        const double xa = mesh.points()[static_cast<std::size_t>(e.first)][0];
        const double xb = mesh.points()[static_cast<std::size_t>(e.second)][0];
        const double mid = 0.5 * (xa + xb);
        const double tolerance = flat_relative * e.length + flat_rounding * std::max({std::abs(xa), std::abs(xb)});
        bool declared = false;
        double best = std::numeric_limits<double>::infinity();
        for (const InterfacePoint1D& p : interfaces) {
            if (!((p.region_a == ra && p.region_b == rb) || (p.region_a == rb && p.region_b == ra))) continue;
            declared = true;
            best = std::min(best, std::abs(mid - p.x));
        }
        if (!declared) {
            return std::unexpected(error(ErrorCode::degenerate_mesh,
                                         "interface rule (P4, 1D): regions " + std::to_string(ra) + " and " +
                                             std::to_string(rb) + " meet at edge " + pair_text(e.first, e.second) +
                                             " with no declared interface",
                                         static_cast<std::size_t>(e.first)));
        }
        ++report.cross_faces;
        report.max_flat_deviation = std::max(report.max_flat_deviation, best / tolerance);
        if (!(best <= tolerance)) {
            return std::unexpected(error(ErrorCode::degenerate_mesh,
                                         "interface rule (P4, 1D): the dual face of edge " + pair_text(e.first, e.second) +
                                             " (its midpoint) is not on a declared interface",
                                         static_cast<std::size_t>(e.first), best));
        }
    }
    std::uint64_t h = 0xcbf29ce484222325ULL;
    digest_word(h, interfaces.size());
    for (const InterfacePoint1D& p : interfaces) {
        digest_word(h, std::bit_cast<std::uint64_t>(p.x));
        digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.region_a)));
        digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.region_b)));
    }
    return InterfaceCertificate{.mesh_digest = graph_digest(mesh),
                                .region_digest = region_digest(node_region),
                                .interfaces_digest = h,
                                .report = report};
}

}  // namespace NiTCAD::mesh
