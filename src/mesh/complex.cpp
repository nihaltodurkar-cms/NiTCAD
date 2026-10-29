// The polyhedral complex, its exact validation and point location, and the sectioned extrusion
// (Unit 53; DECISIONS.md FD1, FD2, FD6).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "complex_index.hpp"
#include "predicates.hpp"

namespace NiTCAD::mesh {

namespace {

using detail::ComplexIndex;
using detail::FacetKind;
using detail::IndexedFacet;
using predicates::orient2d;
using predicates::orient3d;

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}

bool finite(const Point& p) { return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]); }

// The point in the plane of the two coordinates other than `drop`, cyclically.
Point project(const Point& p, int drop) {
    return {p[static_cast<std::size_t>((drop + 1) % 3)], p[static_cast<std::size_t>((drop + 2) % 3)], 0.0};
}

// Collinear in 3D, exactly: collinear in all three coordinate projections.
bool collinear(const Point& a, const Point& b, const Point& c) {
    for (int k = 0; k < 3; ++k) {
        if (orient2d(project(a, k), project(b, k), project(c, k)) != 0) return false;
    }
    return true;
}

// Whether x lies in the closed 2D triangle abc (non-degenerate).
bool in_triangle_2d(const Point& a, const Point& b, const Point& c, const Point& x) {
    const int o1 = orient2d(a, b, x);
    const int o2 = orient2d(b, c, x);
    const int o3 = orient2d(c, a, x);
    return (o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0);
}

// Whether x lies on the closed 2D segment ab (a != b).
bool on_segment_2d(const Point& a, const Point& b, const Point& x) {
    if (orient2d(a, b, x) != 0) return false;
    for (std::size_t k = 0; k < 2; ++k) {
        if (x[k] < std::min(a[k], b[k]) || x[k] > std::max(a[k], b[k])) return false;
    }
    return true;
}

// Whether the open 2D segments ab and cd cross at one point interior to both (not collinear).
bool cross_properly_2d(const Point& a, const Point& b, const Point& c, const Point& d) {
    const int o1 = orient2d(a, b, c);
    const int o2 = orient2d(a, b, d);
    const int o3 = orient2d(c, d, a);
    const int o4 = orient2d(c, d, b);
    return o1 * o2 < 0 && o3 * o4 < 0;
}

// Whether the closed 2D segments ab and cd intersect at all.
bool segments_meet_2d(const Point& a, const Point& b, const Point& c, const Point& d) {
    return cross_properly_2d(a, b, c, d) || on_segment_2d(a, b, c) || on_segment_2d(a, b, d) ||
           on_segment_2d(c, d, a) || on_segment_2d(c, d, b);
}

// Whether a simple 2D polygon is traversed counterclockwise: the turn at its lexicographically
// smallest vertex (never zero for a simple polygon without reversals).
int polygon_turn(const std::vector<Point>& q) {
    std::size_t m = 0;
    for (std::size_t i = 1; i < q.size(); ++i) {
        if (q[i][0] < q[m][0] || (q[i][0] == q[m][0] && q[i][1] < q[m][1])) m = i;
    }
    const std::size_t n = q.size();
    return orient2d(q[(m + n - 1) % n], q[m], q[(m + 1) % n]);
}

// Whether a 2D polygon is simple: consecutive edges meet only at their shared vertex (collinear
// continuation allowed, reversal not), other edges not at all.
bool simple_2d(const std::vector<Point>& q) {
    const std::size_t n = q.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = q[i];
        const Point& b = q[(i + 1) % n];
        const Point& c = q[(i + 2) % n];
        if (orient2d(a, b, c) == 0 && predicates::dot_sign(b, a, c) >= 0) return false;  // a reversal
        for (std::size_t j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) continue;  // consecutive through vertex 0
            if (segments_meet_2d(a, b, q[j], q[(j + 1) % n])) return false;
        }
    }
    return true;
}

// Whether p lies strictly inside a simple 2D polygon (crossing number with the half-open rule;
// points on the boundary are not inside).
bool strictly_inside_2d(const std::vector<Point>& q, const Point& p) {
    const std::size_t n = q.size();
    bool in = false;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = q[i];
        const Point& b = q[(i + 1) % n];
        if (on_segment_2d(a, b, p)) return false;
        if ((a[1] > p[1]) != (b[1] > p[1])) {
            const int o = orient2d(a, b, p);
            if (b[1] > a[1] ? o > 0 : o < 0) in = !in;
        }
    }
    return in;
}

// Ear clipping of a simple 2D polygon with `turn` its orientation: an ear is a strictly convex
// vertex whose triangle (closed) holds no other remaining vertex, so no triangle is degenerate.
// Returns local index triples in the polygon's orientation, or nothing if no ear is found.
std::optional<std::vector<std::array<std::size_t, 3>>> ear_clip(const std::vector<Point>& q, int turn) {
    std::vector<std::size_t> left(q.size());
    for (std::size_t i = 0; i < q.size(); ++i) left[i] = i;
    std::vector<std::array<std::size_t, 3>> out;
    while (left.size() > 3) {
        bool found = false;
        const std::size_t n = left.size();
        for (std::size_t k = 0; k < n && !found; ++k) {
            const std::size_t a = left[(k + n - 1) % n];
            const std::size_t b = left[k];
            const std::size_t c = left[(k + 1) % n];
            if (orient2d(q[a], q[b], q[c]) != turn) continue;
            bool empty = true;
            for (const std::size_t v : left) {
                if (v != a && v != b && v != c && in_triangle_2d(q[a], q[b], q[c], q[v])) {
                    empty = false;
                    break;
                }
            }
            if (!empty) continue;
            out.push_back({a, b, c});
            left.erase(left.begin() + static_cast<std::ptrdiff_t>(k));
            found = true;
        }
        if (!found) return std::nullopt;
    }
    if (orient2d(q[left[0]], q[left[1]], q[left[2]]) != turn) return std::nullopt;
    out.push_back({left[0], left[1], left[2]});
    return out;
}

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

struct Tri {
    std::array<std::size_t, 3> v;
    std::size_t facet;
    Point lo, hi;
};

// Whether the closed segment pq meets the closed triangle t only inside conv(S), S the vertices
// shared by the two triangles (by index). Not coplanar: the meeting is at most one point, allowed
// only at a shared end of pq. Coplanar (in t's projection): no end of pq outside S in t, no proper
// crossing of an edge of t, no collinear overlap with an edge of t except the shared edge itself,
// no vertex of t outside S on pq; together these cover every way the convex intersection can
// leave conv(S) (DECISIONS.md FD1).
bool segment_meets_only_shared(const std::vector<Point>& P, std::size_t pi, std::size_t qi, const Tri& t, int drop,
                               const std::vector<std::size_t>& S) {
    const auto shared = [&](std::size_t v) { return std::ranges::find(S, v) != S.end(); };
    const Point& a = P[t.v[0]];
    const Point& b = P[t.v[1]];
    const Point& c = P[t.v[2]];
    const Point& p = P[pi];
    const Point& q = P[qi];
    const int op = orient3d(a, b, c, p);
    const int oq = orient3d(a, b, c, q);
    if (op * oq > 0) return true;
    const Point A = project(a, drop);
    const Point B = project(b, drop);
    const Point C = project(c, drop);
    if (op != 0 || oq != 0) {
        if (op == 0) return !in_triangle_2d(A, B, C, project(p, drop)) || shared(pi);
        if (oq == 0) return !in_triangle_2d(A, B, C, project(q, drop)) || shared(qi);
        const int s1 = orient3d(p, q, a, b);
        const int s2 = orient3d(p, q, b, c);
        const int s3 = orient3d(p, q, c, a);
        const bool inside = (s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0);
        return !inside;
    }
    const Point Pp = project(p, drop);
    const Point Pq = project(q, drop);
    if (!shared(pi) && in_triangle_2d(A, B, C, Pp)) return false;
    if (!shared(qi) && in_triangle_2d(A, B, C, Pq)) return false;
    const std::array<std::size_t, 3> tv = t.v;
    for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t xi = tv[k];
        const std::size_t yi = tv[(k + 1) % 3];
        const Point X = project(P[xi], drop);
        const Point Y = project(P[yi], drop);
        if (cross_properly_2d(Pp, Pq, X, Y)) return false;
        if (orient2d(Pp, Pq, X) == 0 && orient2d(Pp, Pq, Y) == 0) {
            // Collinear: the overlap along a coordinate on which the line is not constant.
            const std::size_t m = Pp[0] != Pq[0] ? 0 : 1;
            const double lo = std::max(std::min(Pp[m], Pq[m]), std::min(X[m], Y[m]));
            const double hi = std::min(std::max(Pp[m], Pq[m]), std::max(X[m], Y[m]));
            if (lo < hi) {
                const bool same = (xi == pi && yi == qi) || (xi == qi && yi == pi);
                if (!same) return false;
            } else if (lo == hi) {
                // One touching point: an end of both segments; it must be a shared vertex.
                const std::size_t at = Pp[m] == lo ? pi : qi;
                if (!shared(at)) return false;
            }
        }
        if (!shared(xi) && on_segment_2d(Pp, Pq, X)) return false;
    }
    return true;
}

}  // namespace

namespace detail {

std::expected<ComplexIndex, base::Error> ComplexIndex::create(PolyhedralGeometry geometry) {
    ComplexIndex out;
    out.g_ = std::move(geometry);
    const PolyhedralGeometry& g = out.g_;
    const std::vector<Point>& P = g.vertices;
    if (P.empty() || g.facets.empty() || g.volumes.empty()) return std::unexpected(invalid("the complex is empty"));
    for (std::size_t i = 0; i < P.size(); ++i) {
        if (!finite(P[i])) return std::unexpected(invalid("a vertex is not finite", i));
    }
    {
        std::vector<std::size_t> order(P.size());
        for (std::size_t i = 0; i < P.size(); ++i) order[i] = i;
        std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return P[a] < P[b]; });
        for (std::size_t k = 1; k < order.size(); ++k) {
            if (P[order[k]] == P[order[k - 1]]) return std::unexpected(invalid("two vertices coincide", std::max(order[k], order[k - 1])));
        }
    }
    out.lo_ = P[0];
    out.hi_ = P[0];
    for (const Point& p : P) {
        for (std::size_t k = 0; k < 3; ++k) {
            out.lo_[k] = std::min(out.lo_[k], p[k]);
            out.hi_[k] = std::max(out.hi_[k], p[k]);
        }
    }

    // Facets: indices, planarity, simplicity, triangulation.
    out.facets_.resize(g.facets.size());
    for (std::size_t f = 0; f < g.facets.size(); ++f) {
        const std::vector<std::size_t>& v = g.facets[f].vertices;
        IndexedFacet& F = out.facets_[f];
        if (v.size() < 3) return std::unexpected(invalid("a facet has fewer than three vertices", f));
        std::set<std::size_t> seen;
        for (const std::size_t i : v) {
            if (i >= P.size()) return std::unexpected(invalid("a facet's vertex index is out of range", f));
            if (!seen.insert(i).second) return std::unexpected(invalid("a facet repeats a vertex", f));
        }
        std::optional<std::size_t> third;
        for (std::size_t k = 2; k < v.size() && !third; ++k) {
            if (!collinear(P[v[0]], P[v[1]], P[v[k]])) third = k;
        }
        if (!third) return std::unexpected(invalid("a facet's vertices are collinear", f));
        F.plane = {v[0], v[1], v[*third]};
        for (const std::size_t i : v) {
            if (orient3d(P[F.plane[0]], P[F.plane[1]], P[F.plane[2]], P[i]) != 0) {
                return std::unexpected(invalid("a facet is not planar", f));
            }
        }
        F.drop = -1;
        for (int k = 2; k >= 0 && F.drop < 0; --k) {
            if (orient2d(project(P[F.plane[0]], k), project(P[F.plane[1]], k), project(P[F.plane[2]], k)) != 0) F.drop = k;
        }
        std::vector<Point> q;
        q.reserve(v.size());
        for (const std::size_t i : v) q.push_back(project(P[i], F.drop));
        if (!simple_2d(q)) return std::unexpected(invalid("a facet is not a simple polygon", f));
        const auto ears = ear_clip(q, polygon_turn(q));
        if (!ears) return std::unexpected(invalid("a facet could not be triangulated", f));
        for (const auto& e : *ears) F.triangles.push_back({v[e[0]], v[e[1]], v[e[2]]});
    }

    // Volumes: closed, oriented, manifold, non-zero; each facet on one or two volumes.
    std::vector<std::vector<int>> on(g.facets.size());
    std::vector<std::vector<int>> outward(g.volumes.size());  // per volume and listed facet: +1 or -1
    for (std::size_t c = 0; c < g.volumes.size(); ++c) {
        const ComplexVolume& V = g.volumes[c];
        if (V.region < 0) return std::unexpected(invalid("a volume's region is negative", c));
        if (V.facets.empty()) return std::unexpected(invalid("a volume has no facets", c));
        std::set<std::size_t> distinct;
        for (const std::size_t f : V.facets) {
            if (f >= g.facets.size()) return std::unexpected(invalid("a volume's facet index is out of range", c));
            if (!distinct.insert(f).second) return std::unexpected(invalid("a volume lists a facet twice", c));
            on[f].push_back(static_cast<int>(c));
        }
        // Each undirected edge: the listed facets using it, with the direction they traverse it.
        std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<std::size_t, int>>> edges;
        for (std::size_t k = 0; k < V.facets.size(); ++k) {
            const auto& v = g.facets[V.facets[k]].vertices;
            for (std::size_t i = 0; i < v.size(); ++i) {
                const std::size_t a = v[i];
                const std::size_t b = v[(i + 1) % v.size()];
                edges[{std::min(a, b), std::max(a, b)}].emplace_back(k, a < b ? 1 : -1);
            }
        }
        for (const auto& [e, users] : edges) {
            if (users.size() != 2) return std::unexpected(invalid("a volume is not closed edge to edge (an edge not in exactly two of its facets)", c));
        }
        // Orientation by breadth-first search: neighbours traverse a shared edge in opposite directions.
        std::vector<int> flip(V.facets.size(), 0);
        std::vector<std::size_t> queue{0};
        flip[0] = 1;
        std::map<std::size_t, std::vector<std::pair<std::size_t, std::pair<std::size_t, std::size_t>>>> by_facet;
        for (const auto& [e, users] : edges) {
            by_facet[users[0].first].push_back({users[1].first, e});
            by_facet[users[1].first].push_back({users[0].first, e});
        }
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const std::size_t k = queue[head];
            for (const auto& [m, e] : by_facet[k]) {
                const auto& users = edges[e];
                const int dk = users[0].first == k ? users[0].second : users[1].second;
                const int dm = users[0].first == m ? users[0].second : users[1].second;
                const int want = -dk * flip[k] * dm;
                if (flip[m] == 0) {
                    flip[m] = want;
                    queue.push_back(m);
                } else if (flip[m] != want) {
                    return std::unexpected(invalid("a volume is not orientable", c));
                }
            }
        }
        if (queue.size() != V.facets.size()) return std::unexpected(invalid("a volume's facets are not one connected surface", c));
        // No pinched vertices: around each vertex the facets form one cycle.
        std::map<std::size_t, std::vector<std::size_t>> at_vertex;
        for (std::size_t k = 0; k < V.facets.size(); ++k) {
            for (const std::size_t i : g.facets[V.facets[k]].vertices) at_vertex[i].push_back(k);
        }
        for (const auto& [x, around] : at_vertex) {
            std::set<std::size_t> reached{around[0]};
            std::vector<std::size_t> stack{around[0]};
            while (!stack.empty()) {
                const std::size_t k = stack.back();
                stack.pop_back();
                for (const auto& [m, e] : by_facet[k]) {
                    if ((e.first == x || e.second == x) && reached.insert(m).second) stack.push_back(m);
                }
            }
            if (reached.size() != around.size()) return std::unexpected(invalid("a volume is pinched at a vertex", c));
        }
        double vol = 0.0;
        for (std::size_t k = 0; k < V.facets.size(); ++k) {
            for (const auto& t : out.facets_[V.facets[k]].triangles) {
                const Point& a = P[t[0]];
                const Point& b = P[t[1]];
                const Point& d = P[t[2]];
                const double det = a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
                                   a[2] * (b[0] * d[1] - b[1] * d[0]);
                vol += flip[k] * det / 6.0;
            }
        }
        if (vol < 0.0) {
            vol = -vol;
            for (int& s : flip) s = -s;
        }
        double extent = 0.0;
        for (std::size_t k = 0; k < 3; ++k) extent = std::max(extent, out.hi_[k] - out.lo_[k]);
        if (!(vol > 1e-12 * extent * extent * extent)) return std::unexpected(invalid("a volume has zero volume", c));
        outward[c] = flip;
        out.report_.volume.push_back(vol);
    }
    for (std::size_t f = 0; f < g.facets.size(); ++f) {
        IndexedFacet& F = out.facets_[f];
        if (on[f].empty()) return std::unexpected(invalid("a facet is on no volume", f));
        if (on[f].size() > 2) return std::unexpected(invalid("a facet is on more than two volumes", f));
        F.volumes[0] = on[f][0];
        if (on[f].size() == 1) {
            F.kind = FacetKind::boundary;
            if (g.facets[f].patch.empty()) return std::unexpected(invalid("a boundary facet has no patch name", f));
            auto it = std::ranges::find(out.report_.patches, g.facets[f].patch);
            if (it == out.report_.patches.end()) {
                out.report_.patches.push_back(g.facets[f].patch);
                it = out.report_.patches.end() - 1;
            }
            F.patch = static_cast<int>(it - out.report_.patches.begin());
            ++out.report_.boundary_facets;
            continue;
        }
        F.volumes[1] = on[f][1];
        // The two volumes must lie on opposite sides: their outward orientations of f are opposite.
        std::array<int, 2> o{};
        for (std::size_t s = 0; s < 2; ++s) {
            const auto c = static_cast<std::size_t>(F.volumes[s]);
            const auto& list = g.volumes[c].facets;
            o[s] = outward[c][static_cast<std::size_t>(std::ranges::find(list, f) - list.begin())];
        }
        if (o[0] == o[1]) return std::unexpected(invalid("two volumes lie on the same side of a facet (they overlap)", f));
        const bool same = g.volumes[static_cast<std::size_t>(F.volumes[0])].region == g.volumes[static_cast<std::size_t>(F.volumes[1])].region;
        F.kind = same ? FacetKind::internal : FacetKind::interface;
        ++(same ? out.report_.internal_facets : out.report_.interface_facets);
    }

    // Facets meet only in shared vertices and edges: every pair of triangles from different facets.
    std::vector<Tri> tris;
    for (std::size_t f = 0; f < out.facets_.size(); ++f) {
        for (const auto& t : out.facets_[f].triangles) {
            Tri x{t, f, P[t[0]], P[t[0]]};
            for (const std::size_t i : t) {
                for (std::size_t k = 0; k < 3; ++k) {
                    x.lo[k] = std::min(x.lo[k], P[i][k]);
                    x.hi[k] = std::max(x.hi[k], P[i][k]);
                }
            }
            tris.push_back(x);
        }
    }
    std::vector<std::size_t> order(tris.size());
    for (std::size_t i = 0; i < tris.size(); ++i) order[i] = i;
    std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return tris[a].lo[0] < tris[b].lo[0]; });
    for (std::size_t ii = 0; ii < order.size(); ++ii) {
        const Tri& t1 = tris[order[ii]];
        for (std::size_t jj = ii + 1; jj < order.size(); ++jj) {
            const Tri& t2 = tris[order[jj]];
            if (t2.lo[0] > t1.hi[0]) break;
            if (t1.facet == t2.facet) continue;
            if (t2.lo[1] > t1.hi[1] || t1.lo[1] > t2.hi[1] || t2.lo[2] > t1.hi[2] || t1.lo[2] > t2.hi[2]) continue;
            std::vector<std::size_t> S;
            for (const std::size_t a : t1.v) {
                if (std::ranges::find(t2.v, a) != t2.v.end()) S.push_back(a);
            }
            bool ok = S.size() < 3;
            for (std::size_t k = 0; k < 3 && ok; ++k) {
                ok = segment_meets_only_shared(P, t1.v[k], t1.v[(k + 1) % 3], t2, out.facets_[t2.facet].drop, S) &&
                     segment_meets_only_shared(P, t2.v[k], t2.v[(k + 1) % 3], t1, out.facets_[t1.facet].drop, S);
            }
            if (!ok) return std::unexpected(invalid("two facets intersect other than along shared edges and vertices", std::max(t1.facet, t2.facet)));
        }
    }

    // Volumes with disjoint surfaces must not nest: a vertex of one not on the other is outside it.
    std::vector<std::pair<Point, Point>> vbox(g.volumes.size(), {out.hi_, out.lo_});
    std::vector<std::set<std::size_t>> vverts(g.volumes.size());
    for (std::size_t c = 0; c < g.volumes.size(); ++c) {
        for (const std::size_t f : g.volumes[c].facets) {
            for (const std::size_t i : g.facets[f].vertices) {
                vverts[c].insert(i);
                for (std::size_t k = 0; k < 3; ++k) {
                    vbox[c].first[k] = std::min(vbox[c].first[k], P[i][k]);
                    vbox[c].second[k] = std::max(vbox[c].second[k], P[i][k]);
                }
            }
        }
    }
    for (std::size_t a = 0; a < g.volumes.size(); ++a) {
        for (std::size_t b = 0; b < g.volumes.size(); ++b) {
            if (a == b) continue;
            bool overlap = true;
            for (std::size_t k = 0; k < 3; ++k) {
                overlap = overlap && vbox[a].first[k] <= vbox[b].second[k] && vbox[b].first[k] <= vbox[a].second[k];
            }
            if (!overlap) continue;
            const auto it = std::ranges::find_if(vverts[a], [&](std::size_t i) { return !vverts[b].contains(i); });
            if (it == vverts[a].end()) continue;
            ComplexIndex only_b;
            only_b.g_.vertices = P;
            only_b.lo_ = out.lo_;
            only_b.hi_ = out.hi_;
            for (const std::size_t f : g.volumes[b].facets) {
                IndexedFacet F = out.facets_[f];
                F.volumes = {0, -1};
                only_b.facets_.push_back(std::move(F));
            }
            const Location l = only_b.locate(P[*it]);
            if (l.kind != Location::Kind::outside) return std::unexpected(invalid("two volumes overlap", std::max(a, b)));
        }
    }

    // The volumes connected through shared facets.
    std::vector<bool> reached(g.volumes.size(), false);
    std::vector<int> stack{0};
    reached[0] = true;
    while (!stack.empty()) {
        const int c = stack.back();
        stack.pop_back();
        for (const std::size_t f : g.volumes[static_cast<std::size_t>(c)].facets) {
            for (const int d : out.facets_[f].volumes) {
                if (d >= 0 && !reached[static_cast<std::size_t>(d)]) {
                    reached[static_cast<std::size_t>(d)] = true;
                    stack.push_back(d);
                }
            }
        }
    }
    const auto lone = std::ranges::find(reached, false);
    if (lone != reached.end()) return std::unexpected(invalid("the volumes are not connected through shared facets", static_cast<std::size_t>(lone - reached.begin())));
    return out;
}

ComplexIndex::Location ComplexIndex::locate(const Point& p) const {
    const std::vector<Point>& P = g_.vertices;
    for (std::size_t f = 0; f < facets_.size(); ++f) {
        const IndexedFacet& F = facets_[f];
        if (orient3d(P[F.plane[0]], P[F.plane[1]], P[F.plane[2]], p) != 0) continue;
        const Point x = project(p, F.drop);
        for (const auto& t : F.triangles) {
            if (in_triangle_2d(project(P[t[0]], F.drop), project(P[t[1]], F.drop), project(P[t[2]], F.drop), x)) {
                return {Location::Kind::on_facet, static_cast<int>(f)};
            }
        }
    }
    double extent = 0.0;
    for (std::size_t k = 0; k < 3; ++k) extent = std::max(extent, hi_[k] - lo_[k]);
    extent = std::max(extent, 1.0e-30);
    for (std::uint64_t attempt = 0; attempt < 64; ++attempt) {
        // A far point beyond the box in +x, at a deterministic pseudo-random height and depth.
        const auto unit = [&](std::uint64_t salt) {
            return static_cast<double>(splitmix64(0x436F6D706C6578ull ^ (attempt * 2 + salt)) >> 11) * 0x1.0p-53;
        };
        const Point q{hi_[0] + 2.0 * extent, lo_[1] + extent * (3.0 * unit(0) - 1.0), lo_[2] + extent * (3.0 * unit(1) - 1.0)};
        Point slo{}, shi{};
        for (std::size_t k = 0; k < 3; ++k) {
            slo[k] = std::min(p[k], q[k]);
            shi[k] = std::max(p[k], q[k]);
        }
        std::vector<int> parity;
        bool degenerate = false;
        for (std::size_t f = 0; f < facets_.size() && !degenerate; ++f) {
            const IndexedFacet& F = facets_[f];
            for (const auto& t : F.triangles) {
                const Point& a = P[t[0]];
                const Point& b = P[t[1]];
                const Point& c = P[t[2]];
                bool apart = false;
                for (std::size_t k = 0; k < 3; ++k) {
                    apart = apart || std::max({a[k], b[k], c[k]}) < slo[k] || std::min({a[k], b[k], c[k]}) > shi[k];
                }
                if (apart) continue;
                const int op = orient3d(a, b, c, p);
                const int oq = orient3d(a, b, c, q);
                if (op * oq >= 0) continue;  // same side, or touching the plane at p or q only (p is on no facet; q is outside the box)
                const int s1 = orient3d(p, q, a, b);
                const int s2 = orient3d(p, q, b, c);
                const int s3 = orient3d(p, q, c, a);
                if (s1 == 0 || s2 == 0 || s3 == 0) {
                    if ((s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0)) {
                        degenerate = true;  // through an edge or a vertex
                        break;
                    }
                    continue;
                }
                if (s1 != s2 || s2 != s3) continue;
                for (const int v : F.volumes) {
                    if (v < 0) continue;
                    if (parity.size() <= static_cast<std::size_t>(v)) parity.resize(static_cast<std::size_t>(v) + 1, 0);
                    parity[static_cast<std::size_t>(v)] ^= 1;
                }
            }
        }
        if (degenerate) continue;
        for (std::size_t v = 0; v < parity.size(); ++v) {
            if (parity[v] != 0) return {Location::Kind::inside, static_cast<int>(v)};
        }
        return {Location::Kind::outside, -1};
    }
    return {Location::Kind::undetermined, -1};
}

}  // namespace detail

std::expected<ComplexReport, base::Error> check_complex(const PolyhedralGeometry& geometry) {
    auto index = ComplexIndex::create(geometry);
    if (!index) return std::unexpected(std::move(index.error()));
    return index->report();
}

std::expected<PolyhedralGeometry, base::Error> extrude(const ExtrudedGeometry& x) {
    if (x.axis < 0 || x.axis > 2) return std::unexpected(invalid("the extrusion axis must be 0, 1 or 2"));
    const std::vector<double>& t = x.stations;
    if (t.size() < 2) return std::unexpected(invalid("an extrusion needs at least two stations"));
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (!std::isfinite(t[i])) return std::unexpected(invalid("a station is not finite", i));
        if (i > 0 && !(t[i] > t[i - 1])) return std::unexpected(invalid("stations must increase strictly", i));
    }
    const std::size_t sections = t.size() - 1;
    if (x.cells.empty()) return std::unexpected(invalid("an extrusion needs at least one cell"));
    if (x.regions.size() != x.cells.size()) return std::unexpected(invalid("regions needs one row per cell"));
    // Cells: simple, non-zero, uniquely named; stored counterclockwise.
    std::vector<std::vector<Point>> poly(x.cells.size());
    std::set<std::string> names;
    for (std::size_t c = 0; c < x.cells.size(); ++c) {
        const ExtrusionCell& cell = x.cells[c];
        if (cell.name.empty() || !names.insert(cell.name).second) return std::unexpected(invalid("cell names must be non-empty and unique", c));
        if (cell.polygon.size() < 3) return std::unexpected(invalid("a cell has fewer than three vertices", c));
        for (const Point& p : cell.polygon) {
            if (!std::isfinite(p[0]) || !std::isfinite(p[1])) return std::unexpected(invalid("a cell vertex is not finite", c));
            poly[c].push_back({p[0], p[1], 0.0});
        }
        std::set<std::pair<double, double>> distinct;
        for (const Point& p : poly[c]) distinct.insert({p[0], p[1]});
        if (distinct.size() != poly[c].size()) return std::unexpected(invalid("a cell repeats a vertex", c));
        if (!simple_2d(poly[c])) return std::unexpected(invalid("a cell is not a simple polygon", c));
        const int turn = polygon_turn(poly[c]);
        if (turn == 0) return std::unexpected(invalid("a cell has zero area", c));
        if (turn < 0) std::ranges::reverse(poly[c]);
        if (x.regions[c].size() != sections) return std::unexpected(invalid("a cell's regions need one entry per section", c));
        for (const std::int32_t r : x.regions[c]) {
            if (r < 0 && r != ExtrudedGeometry::void_region) return std::unexpected(invalid("a region index must be >= 0 or void", c));
        }
    }
    // The tiling: cells meet only in shared vertices and whole shared edges (opposite directions);
    // no vertex of one strictly inside another.
    std::map<std::pair<std::pair<double, double>, std::pair<double, double>>, std::size_t> directed;
    for (std::size_t c = 0; c < poly.size(); ++c) {
        const auto& q = poly[c];
        for (std::size_t i = 0; i < q.size(); ++i) {
            const auto& a = q[i];
            const auto& b = q[(i + 1) % q.size()];
            if (!directed.emplace(std::pair{std::pair{a[0], a[1]}, std::pair{b[0], b[1]}}, c).second) {
                return std::unexpected(invalid("two cells overlap (an edge traversed the same way twice)", c));
            }
        }
    }
    for (std::size_t a = 0; a < poly.size(); ++a) {
        for (std::size_t b = a + 1; b < poly.size(); ++b) {
            for (const Point& v : poly[a]) {
                if (strictly_inside_2d(poly[b], v)) return std::unexpected(invalid("two cells overlap", b));
            }
            for (const Point& v : poly[b]) {
                if (strictly_inside_2d(poly[a], v)) return std::unexpected(invalid("two cells overlap", b));
            }
            for (std::size_t i = 0; i < poly[a].size(); ++i) {
                const Point& p = poly[a][i];
                const Point& q = poly[a][(i + 1) % poly[a].size()];
                for (std::size_t j = 0; j < poly[b].size(); ++j) {
                    const Point& r = poly[b][j];
                    const Point& s = poly[b][(j + 1) % poly[b].size()];
                    const bool same_edge = (p == r && q == s) || (p == s && q == r);
                    if (same_edge) continue;
                    // Allowed contact: a shared end point only.
                    if (cross_properly_2d(p, q, r, s)) return std::unexpected(invalid("two cells' edges cross", b));
                    for (const auto& [e0, e1, w] : {std::tuple{p, q, r}, std::tuple{p, q, s}, std::tuple{r, s, p}, std::tuple{r, s, q}}) {
                        if (w != e0 && w != e1 && on_segment_2d(e0, e1, w)) {
                            return std::unexpected(invalid("a cell vertex lies inside another cell's edge (not conforming)", b));
                        }
                    }
                    if (orient2d(p, q, r) == 0 && orient2d(p, q, s) == 0) {
                        const std::size_t m = p[0] != q[0] ? 0 : 1;
                        const double lo = std::max(std::min(p[m], q[m]), std::min(r[m], s[m]));
                        const double hi = std::min(std::max(p[m], q[m]), std::max(r[m], s[m]));
                        if (lo < hi) return std::unexpected(invalid("two cells' edges overlap without being one shared edge", b));
                    }
                }
            }
        }
    }
    bool any = false;
    for (const auto& row : x.regions) {
        for (const std::int32_t r : row) any = any || r >= 0;
    }
    if (!any) return std::unexpected(invalid("every cell is void in every section"));
    for (std::size_t k = 0; k < x.patches.size(); ++k) {
        const PatchBox& b = x.patches[k];
        if (b.name.empty() || b.name == x.default_patch) return std::unexpected(invalid("a patch box needs a name other than the default patch", k));
        if (!finite(b.lo) || !finite(b.hi)) return std::unexpected(invalid("a patch box is not finite", k));
        for (std::size_t c = 0; c < 3; ++c) {
            if (b.lo[c] > b.hi[c]) return std::unexpected(invalid("a patch box has lo above hi", k));
        }
    }
    if (x.default_patch.empty()) return std::unexpected(invalid("the default patch name is empty"));

    // The complex.
    PolyhedralGeometry g;
    const auto ax = static_cast<std::size_t>(x.axis);
    const auto uaxis = static_cast<std::size_t>((x.axis + 1) % 3);
    const auto vaxis = static_cast<std::size_t>((x.axis + 2) % 3);
    std::map<std::pair<std::pair<double, double>, std::size_t>, std::size_t> vertex_id;
    const auto vertex = [&](const Point& q, std::size_t station) {
        const auto key = std::pair{std::pair{q[0], q[1]}, station};
        const auto it = vertex_id.find(key);
        if (it != vertex_id.end()) return it->second;
        Point p{};
        p[ax] = t[station];
        p[uaxis] = q[0];
        p[vaxis] = q[1];
        g.vertices.push_back(p);
        vertex_id.emplace(key, g.vertices.size() - 1);
        return g.vertices.size() - 1;
    };
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> station_facet;  // (cell, station)
    std::map<std::tuple<std::pair<double, double>, std::pair<double, double>, std::size_t>, std::size_t> lateral_facet;
    std::vector<std::vector<std::size_t>> facet_volumes;
    const auto facet = [&](std::vector<std::size_t> verts, std::size_t volume) {
        g.facets.push_back({std::move(verts), {}});
        facet_volumes.push_back({volume});
        return g.facets.size() - 1;
    };
    for (std::size_t c = 0; c < poly.size(); ++c) {
        const auto& q = poly[c];
        for (std::size_t i = 0; i < sections; ++i) {
            if (x.regions[c][i] < 0) continue;
            const std::size_t vol = g.volumes.size();
            g.volumes.push_back({{}, x.regions[c][i]});
            for (const std::size_t s : {i, i + 1}) {
                const auto it = station_facet.find({c, s});
                if (it != station_facet.end()) {
                    facet_volumes[it->second].push_back(vol);
                    g.volumes[vol].facets.push_back(it->second);
                    continue;
                }
                std::vector<std::size_t> verts;
                for (const Point& p : q) verts.push_back(vertex(p, s));
                const std::size_t f = facet(std::move(verts), vol);
                station_facet.emplace(std::pair{c, s}, f);
                g.volumes[vol].facets.push_back(f);
            }
            for (std::size_t k = 0; k < q.size(); ++k) {
                const Point& a = q[k];
                const Point& b = q[(k + 1) % q.size()];
                auto ka = std::pair{a[0], a[1]};
                auto kb = std::pair{b[0], b[1]};
                if (kb < ka) std::swap(ka, kb);
                const auto key = std::tuple{ka, kb, i};
                const auto it = lateral_facet.find(key);
                if (it != lateral_facet.end()) {
                    facet_volumes[it->second].push_back(vol);
                    g.volumes[vol].facets.push_back(it->second);
                    continue;
                }
                const std::size_t f = facet({vertex(a, i), vertex(b, i), vertex(b, i + 1), vertex(a, i + 1)}, vol);
                lateral_facet.emplace(key, f);
                g.volumes[vol].facets.push_back(f);
            }
        }
    }
    // Patches on the boundary facets.
    std::vector<bool> used(x.patches.size(), false);
    for (std::size_t f = 0; f < g.facets.size(); ++f) {
        if (facet_volumes[f].size() != 1) continue;
        g.facets[f].patch = x.default_patch;
        for (std::size_t k = 0; k < x.patches.size(); ++k) {
            const PatchBox& b = x.patches[k];
            const bool inside = std::ranges::all_of(g.facets[f].vertices, [&](std::size_t i) {
                const Point& p = g.vertices[i];
                for (std::size_t c = 0; c < 3; ++c) {
                    if (p[c] < b.lo[c] || p[c] > b.hi[c]) return false;
                }
                return true;
            });
            if (inside) {
                g.facets[f].patch = b.name;
                used[k] = true;
                break;
            }
        }
    }
    const auto unused = std::ranges::find(used, false);
    if (unused != used.end()) return std::unexpected(invalid("a patch box selects no boundary facet (none inside it, or all taken by an earlier box)", static_cast<std::size_t>(unused - used.begin())));
    auto checked = check_complex(g);
    if (!checked) return std::unexpected(std::move(checked.error()));
    return g;
}

}  // namespace NiTCAD::mesh
