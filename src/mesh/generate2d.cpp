#include "NiTCAD/mesh/generate.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "dual_geometry.hpp"
#include "predicates.hpp"
#include "triangulation2d.hpp"

namespace NiTCAD::mesh {

namespace {

using detail::Triangulation2D;
using ErrorCode = base::ErrorCode;

constexpr double min_feature_angle_text = 60.0;  // U12, for messages; the test is exact

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}

std::string where(const Point& p) {
    return "(" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ")";
}

bool same(const Point& a, const Point& b) { return a[0] == b[0] && a[1] == b[1]; }
bool less_point(const Point& a, const Point& b) { return std::tie(a[0], a[1]) < std::tie(b[0], b[1]); }
double dist(const Point& a, const Point& b) { return std::hypot(b[0] - a[0], b[1] - a[1]); }

// The circumradius in double: for quality decisions only.
double circumradius(const Point& a, const Point& b, const Point& c) {
    const double ab = dist(a, b), bc = dist(b, c), ca = dist(c, a);
    const double twice_area = std::abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
    return ab * bc * ca / (2.0 * twice_area);
}

// The circumcentre in double: for reports and band tests on existing triangles only.
Point circumcentre_fast(const Point& a, const Point& b, const Point& c) {
    const double bx = b[0] - a[0], by = b[1] - a[1], cx = c[0] - a[0], cy = c[1] - a[1];
    const double d = 2.0 * (bx * cy - by * cx);
    const double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    return {a[0] + (cy * b2 - by * c2) / d, a[1] + (bx * c2 - cx * b2) / d, 0.0};
}

// Whether p (collinear with a, b) lies in the closed segment [a, b]: exact comparisons.
bool within(const Point& p, const Point& a, const Point& b) {
    return std::min(a[0], b[0]) <= p[0] && p[0] <= std::max(a[0], b[0]) && std::min(a[1], b[1]) <= p[1] &&
           p[1] <= std::max(a[1], b[1]);
}

// Whether the closed segments ab and cd meet anywhere other than at a shared end point (where a
// collinear overlap beyond it also counts). Exact.
bool segments_conflict(const Point& a, const Point& b, const Point& c, const Point& d) {
    const int o1 = predicates::orient2d(a, b, c);
    const int o2 = predicates::orient2d(a, b, d);
    const int o3 = predicates::orient2d(c, d, a);
    const int o4 = predicates::orient2d(c, d, b);
    const bool ac = same(a, c), ad = same(a, d), bc = same(b, c), bd = same(b, d);
    if ((ac && bd) || (ad && bc)) return true;  // the same segment twice
    if (ac || ad || bc || bd) {
        // One shared end point: a conflict only if they are collinear and overlap.
        const Point& shared = ac || ad ? a : b;
        const Point& own = ac || ad ? b : a;
        const Point& other = ac || bc ? d : c;
        if (predicates::orient2d(shared, own, other) != 0) return false;
        return predicates::dot_sign(shared, own, other) > 0;
    }
    if (o1 * o2 < 0 && o3 * o4 < 0) return true;
    if (o1 == 0 && within(c, a, b)) return true;
    if (o2 == 0 && within(d, a, b)) return true;
    if (o3 == 0 && within(a, c, d)) return true;
    if (o4 == 0 && within(b, c, d)) return true;
    return false;
}

// Point in a counterclockwise polygon: 1 inside, 0 on its boundary, -1 outside. Exact.
int point_in_polygon(const Point& p, const std::vector<Point>& v) {
    bool inside = false;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const Point& a = v[i];
        const Point& b = v[(i + 1) % v.size()];
        const int o = predicates::orient2d(a, b, p);
        if (o == 0 && within(p, a, b)) return 0;
        if ((a[1] <= p[1]) != (b[1] <= p[1])) {
            // The edge straddles the horizontal line through p; it crosses the ray to +x when p
            // is on its left (going up) or on its right (going down).
            if ((a[1] <= p[1] && o > 0) || (b[1] <= p[1] && o < 0)) inside = !inside;
        }
    }
    return inside ? 1 : -1;
}

double segment_distance(const Point& x, const Point& a, const Point& b) {
    const double ux = b[0] - a[0], uy = b[1] - a[1];
    const double len2 = ux * ux + uy * uy;
    const double t = len2 > 0.0 ? std::clamp(((x[0] - a[0]) * ux + (x[1] - a[1]) * uy) / len2, 0.0, 1.0) : 0.0;
    return std::hypot(x[0] - (a[0] + t * ux), x[1] - (a[1] + t * uy));
}

double segment_segment_distance(const Point& a, const Point& b, const Point& c, const Point& d) {
    if (segments_conflict(a, b, c, d)) return 0.0;
    return std::min({segment_distance(a, c, d), segment_distance(b, c, d), segment_distance(c, a, b),
                     segment_distance(d, a, b)});
}

// The polygons, counterclockwise.
struct Poly {
    std::vector<Point> v;
    std::int32_t region;
};

struct PolyEdge {
    Point a;  // counterclockwise in its polygon: the interior is on the left
    Point b;
    std::size_t poly;
};

// A straight run of interface edges between two features.
struct Chain {
    Point start;
    Point end;
    std::int32_t left;   // the region on the left of start -> end
    std::int32_t right;
    double band = 0.0;   // the protected band's half-width
    std::vector<double> positions;  // pair positions along it, from start
};

enum class NodeKind : std::uint8_t { free, boundary, interface_node };

struct Generator {
    const Geometry2D& g;
    const MeshOptions2D& o;
    std::vector<Poly> polys;
    std::vector<PolyEdge> boundary_edges;
    std::vector<PolyEdge> interface_edges;  // one per shared edge, from the left region's polygon
    std::vector<std::int32_t> interface_right;
    std::vector<Chain> chains;
    std::vector<Point> domain_vertices;  // polygon vertices on the boundary
    std::vector<Point> junctions;        // interface ends on the boundary
    double s = 0.0;
    std::size_t region_count = 0;

    // Initial nodes.
    std::vector<Point> nodes;
    std::vector<NodeKind> kind;
    std::vector<std::int32_t> hint;  // a boundary node's region (from its edge), else -1
    std::vector<std::vector<Point>> edge_nodes;  // per boundary edge, its nodes in order along it

    // The region at p (a polygon containing it, on its boundary counted), or nullopt outside.
    [[nodiscard]] std::optional<std::int32_t> region_at(const Point& p) const {
        std::optional<std::int32_t> r;
        for (const Poly& q : polys) {
            if (point_in_polygon(p, q.v) >= 0) {
                if (r && *r != q.region) return std::nullopt;  // on an interface: never a node
                r = q.region;
            }
        }
        return r;
    }

    [[nodiscard]] double size_at(const Point& p, std::optional<std::int32_t> region) const {
        double h = o.max_size;
        if (region && static_cast<std::size_t>(*region) < o.region_max_size.size() &&
            o.region_max_size[static_cast<std::size_t>(*region)] > 0.0) {
            h = std::min(h, o.region_max_size[static_cast<std::size_t>(*region)]);
        }
        for (const SizeBox& box : o.boxes) {
            if (p[0] >= box.lo[0] && p[0] <= box.hi[0] && p[1] >= box.lo[1] && p[1] <= box.hi[1]) {
                h = std::min(h, box.max_size);
            }
        }
        return h;
    }

    [[nodiscard]] double region_size(std::int32_t r) const {
        double h = o.max_size;
        if (static_cast<std::size_t>(r) < o.region_max_size.size() && o.region_max_size[static_cast<std::size_t>(r)] > 0.0) {
            h = std::min(h, o.region_max_size[static_cast<std::size_t>(r)]);
        }
        return h;
    }

    std::optional<base::Error> validate_options() {
        if (!(std::isfinite(o.max_size) && o.max_size > 0.0)) return invalid("max_size must be finite and positive");
        for (const double h : o.region_max_size) {
            if (!(std::isfinite(h) && h >= 0.0)) return invalid("region_max_size entries must be finite and >= 0");
        }
        for (std::size_t k = 0; k < o.boxes.size(); ++k) {
            const SizeBox& b = o.boxes[k];
            if (!(std::isfinite(b.max_size) && b.max_size > 0.0) || !(b.lo[0] <= b.hi[0] && b.lo[1] <= b.hi[1])) {
                return invalid("size box: max_size must be finite and positive, lo <= hi", k);
            }
        }
        if (!(o.min_angle_deg > 0.0 && o.min_angle_deg <= 33.8)) return invalid("min_angle_deg must be in (0, 33.8]");
        if (o.max_nodes < 3) return invalid("max_nodes must be at least 3");
        return std::nullopt;
    }

    std::optional<base::Error> read_polygons() {
        if (g.polygons.empty()) return invalid("the geometry has no polygons");
        for (std::size_t k = 0; k < g.polygons.size(); ++k) {
            const RegionPolygon& in = g.polygons[k];
            if (in.region < 0 || in.vertices.size() < 3) return invalid("polygon: region < 0 or fewer than three vertices", k);
            for (const Point& p : in.vertices) {
                if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || p[2] != 0.0) {
                    return invalid("polygon vertex not finite or z != 0", k);
                }
            }
            for (std::size_t i = 0; i < in.vertices.size(); ++i) {
                if (same(in.vertices[i], in.vertices[(i + 1) % in.vertices.size()])) {
                    return invalid("polygon has a repeated consecutive vertex", k);
                }
            }
            double area = 0.0;
            double extent = 0.0;
            for (std::size_t i = 0; i < in.vertices.size(); ++i) {
                const Point& a = in.vertices[i];
                const Point& b = in.vertices[(i + 1) % in.vertices.size()];
                area += a[0] * b[1] - b[0] * a[1];
                extent = std::max(extent, dist(a, in.vertices[0]));
            }
            if (!(std::abs(area) > 1e-14 * extent * extent)) return invalid("polygon has no area", k);
            Poly p{in.vertices, in.region};
            if (area < 0.0) std::ranges::reverse(p.v);
            polys.push_back(std::move(p));
            region_count = std::max(region_count, static_cast<std::size_t>(in.region) + 1);
        }
        // Simple, non-overlapping, conforming: no two edges meet except as one shared edge (of
        // two polygons) or at shared end points.
        std::vector<PolyEdge> edges;
        for (std::size_t k = 0; k < polys.size(); ++k) {
            for (std::size_t i = 0; i < polys[k].v.size(); ++i) {
                edges.push_back({polys[k].v[i], polys[k].v[(i + 1) % polys[k].v.size()], k});
            }
        }
        std::vector<int> partner(edges.size(), -1);
        for (std::size_t i = 0; i < edges.size(); ++i) {
            for (std::size_t j = i + 1; j < edges.size(); ++j) {
                const PolyEdge& e = edges[i];
                const PolyEdge& f = edges[j];
                if (same(e.a, f.b) && same(e.b, f.a) && e.poly != f.poly) {
                    if (partner[i] >= 0 || partner[j] >= 0) return invalid("an edge is shared by more than two polygons", e.poly);
                    partner[i] = static_cast<int>(j);
                    partner[j] = static_cast<int>(i);
                    continue;
                }
                if (segments_conflict(e.a, e.b, f.a, f.b)) {
                    return invalid("polygons overlap, are not simple, or meet without sharing whole edges and "
                                   "vertices, near " + where(e.a),
                                   e.poly);
                }
            }
        }
        for (std::size_t k = 0; k < polys.size(); ++k) {
            for (std::size_t q = 0; q < polys.size(); ++q) {
                if (q == k) continue;
                for (const Point& p : polys[k].v) {
                    if (point_in_polygon(p, polys[q].v) > 0) return invalid("a polygon lies inside another", k);
                }
            }
        }
        for (std::size_t i = 0; i < edges.size(); ++i) {
            if (partner[i] < 0) {
                boundary_edges.push_back(edges[i]);
            } else if (partner[i] > static_cast<int>(i)) {
                const PolyEdge& f = edges[static_cast<std::size_t>(partner[i])];
                if (polys[edges[i].poly].region != polys[f.poly].region) {
                    interface_edges.push_back(edges[i]);
                    interface_right.push_back(polys[f.poly].region);
                }
            }
        }
        for (const PolyEdge& e : boundary_edges) {
            if (std::ranges::find_if(domain_vertices, [&](const Point& p) { return same(p, e.a); }) == domain_vertices.end()) {
                domain_vertices.push_back(e.a);
            }
        }
        return std::nullopt;
    }

    std::optional<base::Error> read_patches() {
        std::set<std::string, std::less<>> names;
        for (std::size_t k = 0; k < g.patches.size(); ++k) {
            const PatchSegment& p = g.patches[k];
            if (p.name.empty() || !names.insert(p.name).second) return invalid("patch name empty or repeated", k);
            if (same(p.a, p.b)) return invalid("patch has zero length", k);
            bool on = false;
            for (const PolyEdge& e : boundary_edges) {
                on = on || (predicates::orient2d(e.a, e.b, p.a) == 0 && predicates::orient2d(e.a, e.b, p.b) == 0 &&
                            within(p.a, e.a, e.b) && within(p.b, e.a, e.b));
            }
            if (!on) return invalid("patch is not inside one boundary edge", k);
            for (std::size_t q = 0; q < k; ++q) {
                const PatchSegment& r = g.patches[q];
                if (predicates::orient2d(p.a, p.b, r.a) == 0 && predicates::orient2d(p.a, p.b, r.b) == 0) {
                    const bool overlap = (within(r.a, p.a, p.b) && !same(r.a, p.a) && !same(r.a, p.b)) ||
                                         (within(r.b, p.a, p.b) && !same(r.b, p.a) && !same(r.b, p.b)) ||
                                         (within(p.a, r.a, r.b) && !same(p.a, r.a) && !same(p.a, r.b)) ||
                                         (same(p.a, r.a) && same(p.b, r.b)) || (same(p.a, r.b) && same(p.b, r.a));
                    if (overlap) return invalid("patches overlap", k);
                }
            }
        }
        return std::nullopt;
    }

    // Features, the corner rules (U12) and the chains.
    std::optional<base::Error> read_interfaces() {
        if (interface_edges.empty()) return std::nullopt;
        s = o.interface_resolution;
        if (!(std::isfinite(s) && s > 0.0)) return invalid("interface_resolution must be finite and positive");
        // Vertices of the interface graph.
        std::vector<Point> vertices;
        for (const PolyEdge& e : interface_edges) {
            for (const Point& p : {e.a, e.b}) {
                if (std::ranges::find_if(vertices, [&](const Point& q) { return same(p, q); }) == vertices.end()) {
                    vertices.push_back(p);
                }
            }
        }
        std::ranges::sort(vertices, less_point);
        const auto incident = [&](const Point& v, const std::vector<PolyEdge>& list) {
            std::vector<Point> far;
            for (const PolyEdge& e : list) {
                if (same(e.a, v)) far.push_back(e.b);
                if (same(e.b, v)) far.push_back(e.a);
            }
            return far;
        };
        std::vector<char> feature(vertices.size(), 0);
        for (std::size_t k = 0; k < vertices.size(); ++k) {
            const Point& v = vertices[k];
            const std::vector<Point> arms = incident(v, interface_edges);
            const std::vector<Point> walls = incident(v, boundary_edges);
            if (!walls.empty()) {
                if (walls.size() != 2 || predicates::orient2d(walls[0], v, walls[1]) != 0) {
                    return invalid("an interface meets the boundary at a domain corner (unsupported), at " + where(v));
                }
                if (arms.size() != 1) {
                    return invalid("two or more interfaces meet the boundary at one point (unsupported), at " + where(v));
                }
                feature[k] = 1;
                junctions.push_back(v);
            } else if (arms.size() != 2) {
                feature[k] = 1;
            } else {
                const bool straight = predicates::orient2d(arms[0], v, arms[1]) == 0 &&
                                      predicates::dot_sign(v, arms[0], arms[1]) < 0;
                feature[k] = straight ? 0 : 1;
            }
            if (!feature[k]) continue;
            // Sectors between consecutive walls, counterclockwise; at the boundary the exterior one
            // (from the wall to the edge entering v to the wall to the edge leaving v) is skipped.
            std::vector<Point> rays = arms;
            rays.insert(rays.end(), walls.begin(), walls.end());
            const auto half = [&](const Point& p) { return p[1] > v[1] || (p[1] == v[1] && p[0] > v[0]) ? 0 : 1; };
            std::ranges::sort(rays, [&](const Point& p, const Point& q) {
                if (half(p) != half(q)) return half(p) < half(q);
                return predicates::orient2d(v, p, q) > 0;
            });
            Point entering{}, leaving{};
            for (const PolyEdge& e : boundary_edges) {
                if (same(e.b, v)) entering = e.a;
                if (same(e.a, v)) leaving = e.b;
            }
            for (std::size_t i = 0; i < rays.size(); ++i) {
                const Point& a = rays[i];
                const Point& b = rays[(i + 1) % rays.size()];
                if (!walls.empty() && same(a, entering) && same(b, leaving)) continue;  // the exterior
                const int turn = predicates::orient2d(v, a, b);
                if (rays.size() > 1 && turn > 0 && predicates::angle_below_60(v, a, b)) {
                    return invalid("an interface sector below " + std::to_string(static_cast<int>(min_feature_angle_text)) +
                                   " degrees (refused until justified, DECISIONS.md U12), at " + where(v));
                }
            }
            // Corner nodes (U12): none at a 90-degree bend; else one per sector on its bisector.
            if (walls.empty()) {
                const bool right_bend = arms.size() == 2 && predicates::dot_sign(v, arms[0], arms[1]) == 0;
                if (!right_bend) {
                    for (std::size_t i = 0; i < rays.size(); ++i) {
                        const Point& a = rays[i];
                        const Point& b = rays[(i + 1) % rays.size()];
                        const double ax = a[0] - v[0], ay = a[1] - v[1], bx = b[0] - v[0], by = b[1] - v[1];
                        const double la = std::hypot(ax, ay), lb = std::hypot(bx, by);
                        double dx = ax / la + bx / lb, dy = ay / la + by / lb;
                        const int turn = predicates::orient2d(v, a, b);
                        if (turn < 0 || (rays.size() == 1)) {
                            dx = -dx;
                            dy = -dy;
                        }
                        if (turn == 0) {  // a straight sector: the left normal of a
                            dx = -ay / la;
                            dy = ax / la;
                        }
                        const double ld = std::hypot(dx, dy);
                        add(Point{v[0] + 0.5 * s * dx / ld, v[1] + 0.5 * s * dy / ld, 0.0}, NodeKind::interface_node);
                    }
                }
            }
        }
        // Chains: straight runs through non-feature vertices.
        std::vector<char> used(interface_edges.size(), 0);
        const auto is_feature = [&](const Point& p) {
            const auto it = std::ranges::lower_bound(vertices, p, less_point);
            return it != vertices.end() && same(*it, p) && feature[static_cast<std::size_t>(it - vertices.begin())] != 0;
        };
        for (std::size_t e0 = 0; e0 < interface_edges.size(); ++e0) {
            if (used[e0]) continue;
            // Walk back to a feature, then forward to the next one.
            std::size_t e = e0;
            bool forward = true;  // traversal along a -> b of edge e
            for (std::size_t guard = 0; guard < interface_edges.size(); ++guard) {
                const Point& tail = forward ? interface_edges[e].a : interface_edges[e].b;
                if (is_feature(tail)) break;
                std::optional<std::pair<std::size_t, bool>> prev;
                for (std::size_t f = 0; f < interface_edges.size(); ++f) {
                    if (f == e) continue;
                    if (same(interface_edges[f].b, tail)) prev = std::pair{f, true};
                    if (same(interface_edges[f].a, tail)) prev = std::pair{f, false};
                }
                if (!prev) break;
                e = prev->first;
                forward = prev->second;
            }
            Chain c;
            c.start = forward ? interface_edges[e].a : interface_edges[e].b;
            c.left = forward ? polys[interface_edges[e].poly].region : interface_right[e];
            c.right = forward ? interface_right[e] : polys[interface_edges[e].poly].region;
            for (std::size_t guard = 0; guard <= interface_edges.size(); ++guard) {
                used[e] = 1;
                const Point& head = forward ? interface_edges[e].b : interface_edges[e].a;
                c.end = head;
                if (is_feature(head)) break;
                std::optional<std::pair<std::size_t, bool>> next;
                for (std::size_t f = 0; f < interface_edges.size(); ++f) {
                    if (f == e || used[f]) continue;
                    if (same(interface_edges[f].a, head)) next = std::pair{f, true};
                    if (same(interface_edges[f].b, head)) next = std::pair{f, false};
                }
                if (!next) break;
                e = next->first;
                forward = next->second;
            }
            chains.push_back(c);
        }
        // Pair positions and bands.
        for (Chain& c : chains) {
            const double len = dist(c.start, c.end);
            if (!(len >= 6.0 * s)) {
                return invalid("an interface segment is shorter than 6 s (lower interface_resolution), from " +
                               where(c.start) + " to " + where(c.end));
            }
            const double spacing = std::min(4.0 * s, std::min(region_size(c.left), region_size(c.right)));
            std::vector<double>& t = c.positions;
            t = {1.5 * s, 3.0 * s};
            const double middle = len - 6.0 * s;
            const int steps = std::max(1, static_cast<int>(std::ceil(middle / spacing)));
            for (int k = 1; k < steps; ++k) t.push_back(3.0 * s + middle * k / steps);
            t.push_back(len - 3.0 * s);
            t.push_back(len - 1.5 * s);
            double gap = 1.5 * s;
            for (std::size_t k = 1; k < t.size(); ++k) gap = std::max(gap, t[k] - t[k - 1]);
            c.band = std::hypot(0.5 * gap, 0.5 * s);
            const double ux = (c.end[0] - c.start[0]) / len, uy = (c.end[1] - c.start[1]) / len;
            for (const double x : t) {
                const Point p{c.start[0] + x * ux, c.start[1] + x * uy, 0.0};
                add(Point{p[0] - 0.5 * s * uy, p[1] + 0.5 * s * ux, 0.0}, NodeKind::interface_node);  // left
                add(Point{p[0] + 0.5 * s * uy, p[1] - 0.5 * s * ux, 0.0}, NodeKind::interface_node);  // right
            }
        }
        // Separation (4 s) between unrelated chains, and between chains and boundary edges not at
        // their junctions.
        for (std::size_t i = 0; i < chains.size(); ++i) {
            const Chain& a = chains[i];
            const auto shares = [&](const Point& p, const Point& q) {
                return same(p, a.start) || same(p, a.end) || same(q, a.start) || same(q, a.end);
            };
            for (std::size_t j = i + 1; j < chains.size(); ++j) {
                const Chain& b = chains[j];
                if (shares(b.start, b.end)) continue;
                if (segment_segment_distance(a.start, a.end, b.start, b.end) < 4.0 * s) {
                    return invalid("two interfaces closer than 4 s (lower interface_resolution), near " + where(a.start));
                }
            }
            for (const PolyEdge& e : boundary_edges) {
                if (shares(e.a, e.b)) continue;
                if (segment_segment_distance(a.start, a.end, e.a, e.b) < 4.0 * s) {
                    return invalid("an interface closer than 4 s to the boundary (lower interface_resolution), near " +
                                   where(a.start));
                }
            }
        }
        for (std::size_t k = 0; k < g.patches.size(); ++k) {
            for (const Point& j : junctions) {
                if (dist(g.patches[k].a, j) < 3.0 * s || dist(g.patches[k].b, j) < 3.0 * s) {
                    return invalid("a patch end within 3 s of an interface junction (unsupported)", k);
                }
            }
        }
        return std::nullopt;
    }

    void add(const Point& p, NodeKind k, std::int32_t region_hint = -1) {
        for (const Point& q : nodes) {
            if (same(p, q)) return;
        }
        nodes.push_back(p);
        kind.push_back(k);
        hint.push_back(region_hint);
    }

    // Boundary nodes along each boundary edge: its end points (not junctions), patch ends, the
    // junction straddle (s / 2 on either side), then at most h apart.
    void place_boundary() {
        for (const PolyEdge& e : boundary_edges) {
            edge_nodes.emplace_back();
            const double len = dist(e.a, e.b);
            const double ux = (e.b[0] - e.a[0]) / len, uy = (e.b[1] - e.a[1]) / len;
            std::vector<std::pair<double, NodeKind>> t;
            for (const auto& [end, from_start] : {std::pair{e.a, true}, std::pair{e.b, false}}) {
                const bool junction = std::ranges::find_if(junctions, [&](const Point& j) { return same(j, end); }) != junctions.end();
                if (!junction) {
                    t.emplace_back(from_start ? 0.0 : len, NodeKind::boundary);
                    continue;
                }
                for (const double d : {0.5 * s}) {  // the straddle; refinement splits as needed
                    t.emplace_back(from_start ? d : len - d, NodeKind::interface_node);
                }
            }
            for (const PatchSegment& p : g.patches) {
                for (const Point& q : {p.a, p.b}) {
                    if (predicates::orient2d(e.a, e.b, q) == 0 && within(q, e.a, e.b)) {
                        t.emplace_back(dist(e.a, q), NodeKind::boundary);
                    }
                }
            }
            std::ranges::sort(t);
            std::vector<std::pair<double, NodeKind>> filled;
            const double h = region_size(polys[e.poly].region);
            for (std::size_t k = 0; k < t.size(); ++k) {
                if (k > 0 && t[k].first == t[k - 1].first) continue;
                if (k > 0) {
                    const double gap = t[k].first - filled.back().first;
                    const int steps = static_cast<int>(std::ceil(gap / h));
                    const double base = filled.back().first;
                    for (int m = 1; m < steps; ++m) filled.emplace_back(base + gap * m / steps, NodeKind::boundary);
                }
                filled.push_back(t[k]);
            }
            for (const auto& [x, k] : filled) {
                Point p = x == 0.0 ? e.a : (x == len ? e.b : Point{e.a[0] + x * ux, e.a[1] + x * uy, 0.0});
                // A patch end lies exactly on the edge: use it as given.
                for (const PatchSegment& ps : g.patches) {
                    for (const Point& q : {ps.a, ps.b}) {
                        if (predicates::orient2d(e.a, e.b, q) == 0 && within(q, e.a, e.b) && dist(e.a, q) == x) p = q;
                    }
                }
                add(p, k, polys[e.poly].region);  // on its edge, though rounding may put it a hair outside
                edge_nodes.back().push_back(p);
            }
        }
    }
};

// Hilbert-curve order of the points in their bounding box (2^16 cells a side): the insertion
// order, so that each point is located from a nearby triangle.
std::vector<std::size_t> hilbert_order(const std::vector<Point>& p) {
    double x0 = p[0][0], x1 = p[0][0], y0 = p[0][1], y1 = p[0][1];
    for (const Point& q : p) {
        x0 = std::min(x0, q[0]);
        x1 = std::max(x1, q[0]);
        y0 = std::min(y0, q[1]);
        y1 = std::max(y1, q[1]);
    }
    const double span = std::max({x1 - x0, y1 - y0, 1e-300});
    constexpr std::uint32_t n = 1u << 16;
    std::vector<std::pair<std::uint64_t, std::size_t>> key(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) {
        auto x = static_cast<std::uint32_t>(std::min((p[i][0] - x0) / span * (n - 1), double(n - 1)));
        auto y = static_cast<std::uint32_t>(std::min((p[i][1] - y0) / span * (n - 1), double(n - 1)));
        std::uint64_t d = 0;
        for (std::uint32_t sq = n / 2; sq > 0; sq /= 2) {
            const std::uint32_t rx = (x & sq) > 0 ? 1u : 0u;
            const std::uint32_t ry = (y & sq) > 0 ? 1u : 0u;
            d += static_cast<std::uint64_t>(sq) * sq * ((3u * rx) ^ ry);
            if (ry == 0) {
                if (rx == 1) {
                    x = n - 1 - x;
                    y = n - 1 - y;
                }
                std::swap(x, y);
            }
        }
        key[i] = {d, i};
    }
    std::ranges::sort(key);
    std::vector<std::size_t> order(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) order[i] = key[i].second;
    return order;
}

}  // namespace

std::expected<GeneratedMesh, base::Error> generate_mesh(const Geometry2D& geometry, const MeshOptions2D& options) {
    Generator gen{geometry, options};
    if (auto e = gen.validate_options()) return std::unexpected(*e);
    if (auto e = gen.read_polygons()) return std::unexpected(*e);
    if (auto e = gen.read_patches()) return std::unexpected(*e);
    if (auto e = gen.read_interfaces()) return std::unexpected(*e);
    gen.place_boundary();

    // Triangulate the initial nodes in Hilbert order.
    Point lo = gen.nodes[0], hi = gen.nodes[0];
    for (const Poly& p : gen.polys) {
        for (const Point& q : p.v) {
            lo = {std::min(lo[0], q[0]), std::min(lo[1], q[1]), 0.0};
            hi = {std::max(hi[0], q[0]), std::max(hi[1], q[1]), 0.0};
        }
    }
    Triangulation2D tri(lo, hi);
    std::vector<NodeKind> kind{NodeKind::free, NodeKind::free, NodeKind::free};  // super vertices
    std::vector<std::int32_t> region{-1, -1, -1};
    const auto record = [&](int v, NodeKind k, const Point& p, std::int32_t hinted) {
        if (static_cast<std::size_t>(v) >= kind.size()) {
            kind.resize(static_cast<std::size_t>(v) + 1, NodeKind::free);
            region.resize(static_cast<std::size_t>(v) + 1, -1);
        }
        kind[static_cast<std::size_t>(v)] = k;
        const auto r = hinted >= 0 ? std::optional<std::int32_t>(hinted) : gen.region_at(p);
        region[static_cast<std::size_t>(v)] = r ? *r : -1;
        return r.has_value();
    };
    for (const std::size_t i : hilbert_order(gen.nodes)) {
        const int v = tri.insert(gen.nodes[i]);
        if (!record(v, gen.kind[i], gen.nodes[i], gen.hint[i])) {
            return std::unexpected(base::Error{ErrorCode::degenerate_mesh,
                                               "a generated node is outside the domain or on an interface at " +
                                                   where(gen.nodes[i]),
                                               std::nullopt});
        }
    }

    // Boundary segments between consecutive boundary nodes along each boundary edge.
    struct Segment {
        int a;
        int b;
        std::int32_t patch;   // -1: none
        std::int32_t region;  // of the polygon whose edge it is
        bool alive;
    };
    std::vector<Segment> segments;
    std::map<std::pair<int, int>, std::size_t> segment_of;
    const auto patch_of = [&](const Point& a, const Point& b) -> std::int32_t {
        const Point m{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.0};
        for (std::size_t k = 0; k < geometry.patches.size(); ++k) {
            const PatchSegment& ps = geometry.patches[k];
            const double len = dist(ps.a, ps.b);
            if (segment_distance(m, ps.a, ps.b) <= 1e-12 * len && segment_distance(a, ps.a, ps.b) <= 1e-12 * len &&
                segment_distance(b, ps.a, ps.b) <= 1e-12 * len) {
                return static_cast<std::int32_t>(k);
            }
        }
        return -1;
    };
    const auto add_segment = [&](int a, int b, std::int32_t patch, std::int32_t region) {
        segment_of[{std::min(a, b), std::max(a, b)}] = segments.size();
        segments.push_back({a, b, patch, region, true});
    };
    std::map<std::pair<double, double>, int> vertex_at;
    for (std::size_t v = 3; v < tri.vertex_count(); ++v) {
        const Point& p = tri.point(static_cast<int>(v));
        vertex_at[{p[0], p[1]}] = static_cast<int>(v);
    }
    for (std::size_t edge = 0; edge < gen.edge_nodes.size(); ++edge) {
        const std::vector<Point>& along = gen.edge_nodes[edge];
        const std::int32_t owner = gen.polys[gen.boundary_edges[edge].poly].region;
        for (std::size_t k = 0; k + 1 < along.size(); ++k) {
            const int a = vertex_at.at({along[k][0], along[k][1]});
            const int b = vertex_at.at({along[k + 1][0], along[k + 1][1]});
            add_segment(a, b, patch_of(tri.point(a), tri.point(b)), owner);
        }
    }
    std::set<Point, decltype([](const Point& a, const Point& b) { return less_point(a, b); })> corners(
        gen.domain_vertices.begin(), gen.domain_vertices.end());

    MeshReport report;
    for (std::size_t v = 3; v < kind.size(); ++v) report.interface_nodes += kind[v] == NodeKind::interface_node ? 1 : 0;
    const double sin_min = std::sin(options.min_angle_deg * 3.141592653589793 / 180.0);
    std::deque<std::size_t> segment_queue;
    for (std::size_t k = 0; k < segments.size(); ++k) segment_queue.push_back(k);
    std::deque<int> triangle_queue;
    for (std::size_t t = 0; t < tri.triangle_slots(); ++t) {
        if (tri.triangle(static_cast<int>(t)).alive) triangle_queue.push_back(static_cast<int>(t));
    }

    const auto inserted_triangles = [&](std::size_t before_slots, const Triangulation2D::Cavity& c) {
        // Triangles created by an insertion: the reused cavity slots and any new slots.
        for (const int t : c.triangles) {
            if (tri.triangle(t).alive) triangle_queue.push_back(t);
        }
        for (std::size_t t = before_slots; t < tri.triangle_slots(); ++t) triangle_queue.push_back(static_cast<int>(t));
    };

    const auto node_cap = [&]() { return tri.vertex_count() - 3 >= options.max_nodes; };

    const auto encroached = [&](const Segment& sg) {
        const auto t1 = tri.triangle_with_edge(sg.a, sg.b);
        const auto t2 = tri.triangle_with_edge(sg.b, sg.a);
        if (!t1 && !t2) return true;
        for (const auto& t : {t1, t2}) {
            if (!t) continue;
            const auto& tr = tri.triangle(*t);
            for (const int apex : tr.v) {
                if (apex == sg.a || apex == sg.b || tri.is_super(apex)) continue;
                if (predicates::in_diametral_circle(tri.point(sg.a), tri.point(sg.b), tri.point(apex)) > 0) return true;
            }
        }
        return false;
    };

    const auto split = [&](std::size_t k) -> std::optional<base::Error> {
        Segment sg = segments[k];
        const Point& pa = tri.point(sg.a);
        const Point& pb = tri.point(sg.b);
        const double len = dist(pa, pb);
        double t = 0.5;
        const bool ca = corners.contains(pa), cb = corners.contains(pb);
        if (ca != cb) {
            // Concentric shells around a domain corner: a power of two from it, near the middle.
            const double d = std::exp2(std::round(std::log2(0.5 * len)));
            const double f = d / len;
            if (f > 0.25 && f < 0.75) t = ca ? f : 1.0 - f;
        }
        Point m{pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1]), 0.0};
        if (pa[0] == pb[0]) m[0] = pa[0];  // keep axis-aligned edges exactly straight
        if (pa[1] == pb[1]) m[1] = pa[1];
        if (node_cap()) return base::Error{ErrorCode::resource_exhausted, "max_nodes reached", std::nullopt};
        const std::size_t before = tri.triangle_slots();
        const int start = tri.locate(m, tri.last_triangle());
        const auto cav = tri.cavity(m, start);
        const int v = tri.insert(m, cav);
        if (!record(v, NodeKind::boundary, m, sg.region)) {
            return base::Error{ErrorCode::degenerate_mesh, "a boundary split fell outside the domain at " + where(m), std::nullopt};
        }
        inserted_triangles(before, cav);
        segments[k].alive = false;
        segment_of.erase({std::min(sg.a, sg.b), std::max(sg.a, sg.b)});
        add_segment(sg.a, v, sg.patch, sg.region);
        add_segment(v, sg.b, sg.patch, sg.region);
        segment_queue.push_back(segments.size() - 2);
        segment_queue.push_back(segments.size() - 1);
        ++report.segment_splits;
        return std::nullopt;
    };

    const auto in_band = [&](const Point& c) {
        for (const Chain& ch : gen.chains) {
            if (segment_distance(c, ch.start, ch.end) < ch.band) return true;
        }
        return false;
    };

    const std::size_t iteration_cap = 64 * options.max_nodes + 1024;
    for (std::size_t iteration = 0;; ++iteration) {
        if (iteration > iteration_cap) {
            return std::unexpected(base::Error{ErrorCode::resource_exhausted,
                                               "refinement did not settle within its iteration cap", std::nullopt});
        }
        if (!segment_queue.empty()) {
            const std::size_t k = segment_queue.front();
            segment_queue.pop_front();
            if (!segments[k].alive || !encroached(segments[k])) continue;
            if (auto e = split(k)) return std::unexpected(*e);
            continue;
        }
        if (triangle_queue.empty()) break;
        const int t = triangle_queue.front();
        triangle_queue.pop_front();
        const auto& tr = tri.triangle(t);
        if (!tr.alive || tri.is_super(tr.v[0]) || tri.is_super(tr.v[1]) || tri.is_super(tr.v[2])) continue;
        const Point& a = tri.point(tr.v[0]);
        const Point& b = tri.point(tr.v[1]);
        const Point& c = tri.point(tr.v[2]);
        const Point centroid{(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, 0.0};
        const auto r = gen.region_at(centroid);
        if (!r) continue;  // outside the domain
        // Quality in double (a decision, not a construction); the exact circumcentre only for a
        // point that may be inserted.
        const double radius = circumradius(a, b, c);
        const double shortest = std::min({dist(a, b), dist(b, c), dist(c, a)});
        const bool skinny = shortest / (2.0 * radius) < sin_min;
        const bool large = radius > gen.size_at(centroid, r) / std::sqrt(3.0);
        if (!skinny && !large) continue;
        const bool all_interface = kind[static_cast<std::size_t>(tr.v[0])] == NodeKind::interface_node &&
                                   kind[static_cast<std::size_t>(tr.v[1])] == NodeKind::interface_node &&
                                   kind[static_cast<std::size_t>(tr.v[2])] == NodeKind::interface_node;
        if (all_interface) {
            ++report.exempt_triangles;
            continue;
        }
        const Point cc = detail::triangle_circumcentre({a, b, c});
        if (in_band(cc)) {
            ++report.exempt_triangles;
            continue;
        }
        const int start = tri.locate(cc, t);
        const auto cav = tri.cavity(cc, start);
        // Ruppert: a candidate that encroaches boundary segments is not inserted; those segments
        // are split instead (whether or not an existing vertex encroaches them), and the triangle
        // is tried again.
        std::vector<std::size_t> hit;
        for (const int ct : cav.triangles) {
            const auto& ctr = tri.triangle(ct);
            for (std::size_t i = 0; i < 3; ++i) {
                const int x = ctr.v[(i + 1) % 3], y = ctr.v[(i + 2) % 3];
                const auto it = segment_of.find({std::min(x, y), std::max(x, y)});
                if (it == segment_of.end()) continue;
                if (predicates::in_diametral_circle(tri.point(x), tri.point(y), cc) > 0 &&
                    std::ranges::find(hit, it->second) == hit.end()) {
                    hit.push_back(it->second);
                }
            }
        }
        if (!hit.empty()) {
            for (const std::size_t k : hit) {
                if (!segments[k].alive) continue;
                if (auto e = split(k)) return std::unexpected(*e);
            }
            triangle_queue.push_back(t);
            continue;
        }
        if (!gen.region_at(cc)) {
            ++report.exempt_triangles;
            continue;
        }
        if (node_cap()) {
            return std::unexpected(base::Error{ErrorCode::resource_exhausted, "max_nodes reached", std::nullopt});
        }
        const std::size_t before = tri.triangle_slots();
        const int v = tri.insert(cc, cav);
        record(v, NodeKind::free, cc, -1);
        inserted_triangles(before, cav);
        ++report.refinement_points;
    }

    // Every boundary segment is an edge.
    for (const Segment& sg : segments) {
        if (sg.alive && !tri.triangle_with_edge(sg.a, sg.b) && !tri.triangle_with_edge(sg.b, sg.a)) {
            return std::unexpected(base::Error{ErrorCode::degenerate_mesh, "a boundary segment is missing", std::nullopt});
        }
    }

    // Extract the inside triangles.
    SimplexMesh simplices;
    SimplexMesh& sm = simplices;
    sm.dimension = 2;
    for (std::size_t v = 3; v < tri.vertex_count(); ++v) sm.points.push_back(tri.point(static_cast<int>(v)));
    double min_angle = 180.0;
    for (std::size_t t = 0; t < tri.triangle_slots(); ++t) {
        const auto& tr = tri.triangle(static_cast<int>(t));
        if (!tr.alive || tri.is_super(tr.v[0]) || tri.is_super(tr.v[1]) || tri.is_super(tr.v[2])) continue;
        const Point& a = tri.point(tr.v[0]);
        const Point& b = tri.point(tr.v[1]);
        const Point& c = tri.point(tr.v[2]);
        const Point centroid{(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, 0.0};
        if (!gen.region_at(centroid)) continue;
        for (const int v : tr.v) sm.cells.push_back(static_cast<NodeId>(v - 3));
        const Point cc = circumcentre_fast(a, b, c);
        const bool protected_tri = kind[static_cast<std::size_t>(tr.v[0])] == NodeKind::interface_node ||
                                   kind[static_cast<std::size_t>(tr.v[1])] == NodeKind::interface_node ||
                                   kind[static_cast<std::size_t>(tr.v[2])] == NodeKind::interface_node || in_band(cc);
        if (!protected_tri) {
            const double shortest = std::min({dist(a, b), dist(b, c), dist(c, a)});
            min_angle = std::min(min_angle, std::asin(std::min(1.0, shortest / (2.0 * dist(cc, a)))) * 180.0 / 3.141592653589793);
        }
    }
    sm.patch_names.clear();
    for (const PatchSegment& p : geometry.patches) sm.patch_names.push_back(p.name);
    for (const Segment& sg : segments) {
        if (!sg.alive || sg.patch < 0) continue;
        sm.facets.push_back(static_cast<NodeId>(sg.a - 3));
        sm.facets.push_back(static_cast<NodeId>(sg.b - 3));
        sm.facet_patch.push_back(sg.patch);
    }
    std::vector<std::int32_t> node_region(region.begin() + 3, region.end());
    auto dual = make_dual_mesh(sm);
    if (!dual) {
        return std::unexpected(base::Error{ErrorCode::degenerate_mesh,
                                           "the generated mesh failed an import check: " + dual.error().message,
                                           dual.error().context});
    }
    std::optional<InterfaceCertificate> certificate;
    std::set<std::int32_t> distinct(node_region.begin(), node_region.end());
    if (distinct.size() > 1) {
        // With no interfaces (regions that do not touch) any positive resolution serves.
        InterfaceSet set{{}, gen.s > 0.0 ? gen.s : options.max_size};
        for (std::size_t k = 0; k < gen.interface_edges.size(); ++k) {
            const PolyEdge& e = gen.interface_edges[k];
            set.pieces.push_back({{e.a, e.b}, gen.polys[e.poly].region, gen.interface_right[k]});
        }
        auto cert = certify_interfaces(sm, *dual, node_region, set);
        if (!cert) {
            return std::unexpected(base::Error{ErrorCode::degenerate_mesh,
                                               "the generated mesh did not meet P4 (the corner rule, DECISIONS.md "
                                               "U12): " + cert.error().message,
                                               cert.error().context});
        }
        certificate = std::move(*cert);
    }
    report.nodes = dual->node_count();
    report.triangles = sm.cells.size() / 3;
    report.min_angle_deg = min_angle;
    return GeneratedMesh{std::move(*dual), std::move(simplices), std::move(node_region), std::move(certificate),
                         report};
}

}  // namespace NiTCAD::mesh
