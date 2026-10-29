#include "NiTCAD/mesh/interfaces.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
                                                 " in a corner zone deviates beyond the resolution",
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

std::uint64_t interfaces_digest(const InterfaceSet& set) {
    std::uint64_t h = 0xcbf29ce484222325ULL;
    digest_word(h, std::bit_cast<std::uint64_t>(set.resolution));
    digest_word(h, set.pieces.size());
    for (const InterfacePiece& p : set.pieces) {
        digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.region_a)));
        digest_word(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.region_b)));
        digest_word(h, p.vertices.size());
        for (const Point& v : p.vertices) {
            for (const double c : v) digest_word(h, std::bit_cast<std::uint64_t>(c));
        }
    }
    return h;
}

}  // namespace

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
    const InterfaceSet& interfaces) {
    auto report = check_interfaces(simplices, dual, node_region, interfaces);
    if (!report) return std::unexpected(std::move(report.error()));
    return InterfaceCertificate{.mesh_digest = graph_digest(dual),
                                .region_digest = region_digest(node_region),
                                .interfaces_digest = interfaces_digest(interfaces),
                                .report = *report};
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
