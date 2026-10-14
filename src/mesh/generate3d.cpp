#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "dual_geometry.hpp"
#include "predicates.hpp"
#include "triangulation2d.hpp"
#include "triangulation3d.hpp"
#include "vec_ops.hpp"

namespace NiTCAD::mesh {

namespace {

using detail::Triangulation2D;
using detail::Triangulation3D;
using ErrorCode = base::ErrorCode;

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}

base::Error degenerate(std::string message) { return {ErrorCode::degenerate_mesh, std::move(message), std::nullopt}; }

double dist3(const Point& a, const Point& b) { return std::sqrt(detail::dot3(detail::sub(a, b), detail::sub(a, b))); }

// Collinear in 3D, exactly: collinear in all three coordinate projections.
bool collinear(const Point& a, const Point& b, const Point& c) {
    const auto proj = [](const Point& p, int i, int j) { return Point{p[static_cast<std::size_t>(i)], p[static_cast<std::size_t>(j)], 0.0}; };
    for (const auto& [i, j] : {std::pair{0, 1}, std::pair{1, 2}, std::pair{0, 2}}) {
        if (predicates::orient2d(proj(a, i, j), proj(b, i, j), proj(c, i, j)) != 0) return false;
    }
    return true;
}

// The circumcentre of a tetrahedron in double: for quality decisions only.
Point circumcentre_fast(const Point& a, const Point& b, const Point& c, const Point& d) {
    const Point u = detail::sub(b, a);
    const Point v = detail::sub(c, a);
    const Point w = detail::sub(d, a);
    const double det = 2.0 * detail::dot3(u, detail::cross(v, w));
    const Point num = detail::add(detail::add(detail::scale(detail::dot3(w, w), detail::cross(u, v)),
                                              detail::scale(detail::dot3(v, v), detail::cross(w, u))),
                                  detail::scale(detail::dot3(u, u), detail::cross(v, w)));
    return detail::add(a, detail::scale(1.0 / det, num));
}

// The smallest dihedral angle of a tetrahedron [degrees].
double min_dihedral(const std::array<Point, 4>& p) {
    double best = 180.0;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = i + 1; j < 4; ++j) {
            std::size_t k = 0;
            while (k == i || k == j) ++k;
            std::size_t l = k + 1;
            while (l == i || l == j) ++l;
            const Point e = detail::sub(p[j], p[i]);
            const double ee = detail::dot3(e, e);
            const auto perp = [&](const Point& q) {
                const Point r = detail::sub(q, p[i]);
                return detail::sub(r, detail::scale(detail::dot3(r, e) / ee, e));
            };
            const Point u = perp(p[k]);
            const Point v = perp(p[l]);
            const double c = detail::dot3(u, v) / std::sqrt(detail::dot3(u, u) * detail::dot3(v, v));
            best = std::min(best, std::acos(std::clamp(c, -1.0, 1.0)) * 180.0 / 3.141592653589793);
        }
    }
    return best;
}

struct Face {
    std::vector<std::size_t> vid;  // polyhedron vertex indices, in order
    std::array<Point, 3> plane{};  // three non-collinear vertices
    int inside = 0;                // the sign of orient3d(plane, x) for x inside the polyhedron
    int axis = -1;                 // the coordinate constant on the face, or -1 (oblique)
    Point origin{}, u{}, w{};      // an orthonormal frame in the plane (oblique faces)
    std::vector<Point> local_polygon;  // counterclockwise in face coordinates
    std::unique_ptr<Triangulation2D> tri;
    std::vector<int> global;  // per face-triangulation vertex: the 3D vertex (-1: super)
    std::map<int, int> local;  // 3D vertex -> face-triangulation vertex

    [[nodiscard]] Point to_local(const Point& p) const {
        if (axis >= 0) {
            const auto a = static_cast<std::size_t>((axis + 1) % 3);
            const auto b = static_cast<std::size_t>((axis + 2) % 3);
            return {p[a], p[b], 0.0};
        }
        const Point r = detail::sub(p, origin);
        return {detail::dot3(r, u), detail::dot3(r, w), 0.0};
    }
    // Whether a face-coordinate point lies in the closed polygon.
    [[nodiscard]] bool contains_local(const Point& q) const {
        for (std::size_t i = 0; i < local_polygon.size(); ++i) {
            if (predicates::orient2d(local_polygon[i], local_polygon[(i + 1) % local_polygon.size()], q) < 0) return false;
        }
        return true;
    }
};

struct DomainEdge {
    std::size_t a = 0;
    std::size_t b = 0;
    std::array<int, 2> faces{-1, -1};
};

struct Subsegment {
    int a = 0;  // 3D vertices
    int b = 0;
    int edge = 0;
    bool alive = true;
};

struct Generator3D {
    const Geometry3D& g;
    const MeshOptions3D& o;
    std::vector<Face> faces;
    std::vector<DomainEdge> edges;

    std::optional<base::Error> validate_options() const {
        if (!(std::isfinite(o.max_size) && o.max_size > 0.0)) return invalid("max_size must be finite and positive");
        for (std::size_t k = 0; k < o.boxes.size(); ++k) {
            const SizeBox& b = o.boxes[k];
            if (!(std::isfinite(b.max_size) && b.max_size > 0.0) ||
                !(b.lo[0] <= b.hi[0] && b.lo[1] <= b.hi[1] && b.lo[2] <= b.hi[2])) {
                return invalid("size box: max_size must be finite and positive, lo <= hi", k);
            }
        }
        if (!(std::isfinite(o.radius_edge) && o.radius_edge >= 1.414)) return invalid("radius_edge must be finite and >= 1.414");
        if (g.region < 0) return invalid("region must be >= 0");
        if (o.max_nodes < 4) return invalid("max_nodes must be at least 4");
        return std::nullopt;
    }

    std::optional<base::Error> read_polyhedron() {
        const ConvexPolyhedron& d = g.domain;
        if (d.vertices.size() < 4 || d.faces.size() < 4) return invalid("a polyhedron needs at least four vertices and faces");
        for (std::size_t i = 0; i < d.vertices.size(); ++i) {
            for (const double c : d.vertices[i]) {
                if (!std::isfinite(c)) return invalid("polyhedron vertex not finite", i);
            }
            for (std::size_t j = 0; j < i; ++j) {
                if (d.vertices[i] == d.vertices[j]) return invalid("polyhedron vertices repeated", i);
            }
        }
        std::set<std::string> names;
        std::map<std::pair<std::size_t, std::size_t>, std::vector<int>> edge_faces;
        std::vector<int> used(d.vertices.size(), 0);
        for (std::size_t f = 0; f < d.faces.size(); ++f) {
            const PolyhedronFace& in = d.faces[f];
            if (in.name.empty() || !names.insert(in.name).second) return invalid("face names must be non-empty and unique", f);
            if (in.vertices.size() < 3) return invalid("a face needs at least three vertices", f);
            std::set<std::size_t> distinct;
            for (const std::size_t v : in.vertices) {
                if (v >= d.vertices.size()) return invalid("face vertex index out of range", f);
                if (!distinct.insert(v).second) return invalid("face vertex repeated", f);
                used[v] = 1;
            }
            Face face;
            face.vid = in.vertices;
            const auto P = [&](std::size_t k) -> const Point& { return d.vertices[in.vertices[k % in.vertices.size()]]; };
            // A non-collinear triple of consecutive vertices defines the plane.
            bool found = false;
            for (std::size_t k = 0; k < in.vertices.size() && !found; ++k) {
                if (!collinear(P(k), P(k + 1), P(k + 2))) {
                    face.plane = {P(k), P(k + 1), P(k + 2)};
                    found = true;
                }
            }
            if (!found) return invalid("face has no area", f);
            for (std::size_t k = 0; k < in.vertices.size(); ++k) {
                if (predicates::orient3d(face.plane[0], face.plane[1], face.plane[2], P(k)) != 0) {
                    return invalid("face is not planar", f);
                }
            }
            // The other vertices: strictly on one side (convexity); a vertex on the plane would
            // lie on this face without being one of its vertices.
            for (std::size_t v = 0; v < d.vertices.size(); ++v) {
                if (distinct.contains(v)) continue;
                const int s = predicates::orient3d(face.plane[0], face.plane[1], face.plane[2], d.vertices[v]);
                if (s == 0) return invalid("a polyhedron vertex lies on the plane of a face it is not a vertex of", f);
                if (face.inside == 0) face.inside = s;
                if (s != face.inside) return invalid("the polyhedron is not convex (DECISIONS.md M6, U7)", f);
            }
            // Convex and simple in its plane: every vertex strictly on the same side of every
            // edge, seen from a point off the plane.
            Point off{};
            for (std::size_t v = 0; v < d.vertices.size(); ++v) {
                if (!distinct.contains(v)) {
                    off = d.vertices[v];
                    break;
                }
            }
            int turn = 0;
            for (std::size_t k = 0; k < in.vertices.size(); ++k) {
                for (std::size_t m = 0; m < in.vertices.size(); ++m) {
                    if (m == k || m == (k + 1) % in.vertices.size()) continue;
                    const int s = predicates::orient3d(P(k), P(k + 1), P(m), off);
                    if (s == 0 || (turn != 0 && s != turn)) return invalid("face is not strictly convex and simple", f);
                    turn = s;
                }
                const std::size_t a = in.vertices[k];
                const std::size_t b = in.vertices[(k + 1) % in.vertices.size()];
                edge_faces[{std::min(a, b), std::max(a, b)}].push_back(static_cast<int>(f));
            }
            for (int axis = 0; axis < 3 && face.axis < 0; ++axis) {
                const auto ax = static_cast<std::size_t>(axis);
                if (std::ranges::all_of(in.vertices, [&](std::size_t v) { return d.vertices[v][ax] == P(0)[ax]; })) {
                    face.axis = axis;
                }
            }
            if (face.axis < 0) {
                face.origin = face.plane[0];
                const Point e = detail::sub(face.plane[1], face.plane[0]);
                face.u = detail::scale(1.0 / std::sqrt(detail::dot3(e, e)), e);
                Point n = detail::cross(e, detail::sub(face.plane[2], face.plane[0]));
                n = detail::scale(1.0 / std::sqrt(detail::dot3(n, n)), n);
                face.w = detail::cross(n, face.u);
            }
            for (std::size_t k = 0; k < in.vertices.size(); ++k) face.local_polygon.push_back(face.to_local(P(k)));
            double area = 0.0;
            for (std::size_t k = 0; k < face.local_polygon.size(); ++k) {
                const Point& a = face.local_polygon[k];
                const Point& b = face.local_polygon[(k + 1) % face.local_polygon.size()];
                area += a[0] * b[1] - b[0] * a[1];
            }
            if (area < 0.0) std::ranges::reverse(face.local_polygon);
            faces.push_back(std::move(face));
        }
        if (std::ranges::find(used, 0) != used.end()) return invalid("a polyhedron vertex is in no face");
        for (const auto& [key, fs] : edge_faces) {
            if (fs.size() != 2) return invalid("the polyhedron is not closed: an edge is not in exactly two faces", key.first);
            edges.push_back({key.first, key.second, {fs[0], fs[1]}});
        }
        const auto euler = static_cast<long long>(d.vertices.size()) - static_cast<long long>(edges.size()) +
                           static_cast<long long>(d.faces.size());
        if (euler != 2) return invalid("the polyhedron is not a closed convex surface (V - E + F != 2)");
        return std::nullopt;
    }

    [[nodiscard]] bool inside(const Point& p) const {
        for (const Face& f : faces) {
            const int s = predicates::orient3d(f.plane[0], f.plane[1], f.plane[2], p);
            if (s != 0 && s != f.inside) return false;
        }
        return true;
    }

    [[nodiscard]] double size_at(const Point& p) const {
        double h = o.max_size;
        for (const SizeBox& b : o.boxes) {
            if (p[0] >= b.lo[0] && p[0] <= b.hi[0] && p[1] >= b.lo[1] && p[1] <= b.hi[1] && p[2] >= b.lo[2] &&
                p[2] <= b.hi[2]) {
                h = std::min(h, b.max_size);
            }
        }
        return h;
    }

    [[nodiscard]] double volume() const {
        // Sum over faces of the fan tetrahedra from the first vertex of the polyhedron.
        const Point& apex = g.domain.vertices[0];
        double v = 0.0;
        for (const Face& f : faces) {
            const auto P = [&](std::size_t k) -> const Point& { return g.domain.vertices[f.vid[k]]; };
            for (std::size_t k = 1; k + 1 < f.vid.size(); ++k) {
                v += std::abs(detail::dot3(detail::sub(P(0), apex),
                                           detail::cross(detail::sub(P(k), apex), detail::sub(P(k + 1), apex)))) / 6.0;
            }
        }
        return v;
    }
};

}  // namespace

std::expected<GeneratedMesh, base::Error> generate_mesh(const Geometry3D& geometry, const MeshOptions3D& options) {
    Generator3D gen{geometry, options, {}, {}};
    if (auto e = gen.validate_options()) return std::unexpected(*e);
    if (auto e = gen.read_polyhedron()) return std::unexpected(*e);
    const ConvexPolyhedron& d = geometry.domain;
    if (d.vertices.size() > options.max_nodes) {
        return std::unexpected(base::Error{ErrorCode::resource_exhausted, "max_nodes reached", std::nullopt});
    }

    Point lo = d.vertices[0];
    Point hi = d.vertices[0];
    for (const Point& p : d.vertices) {
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    }
    Triangulation3D tri(lo, hi);
    // Per 3D vertex (super vertices included): the faces and domain edges it lies on.
    std::vector<std::vector<int>> faces_of(4);
    std::vector<std::vector<int>> edges_of(4);
    std::vector<std::vector<int>> segments_of(4);  // the subsegments it ends

    std::deque<int> tet_queue;
    std::deque<std::pair<int, int>> facet_queue;  // (face, face-triangulation triangle)
    std::deque<int> segment_queue;
    std::vector<Subsegment> segments;
    std::map<std::pair<int, int>, int> segment_at;
    MeshReport report;

    for (Face& f : gen.faces) {
        Point flo = f.local_polygon[0];
        Point fhi = f.local_polygon[0];
        for (const Point& q : f.local_polygon) {
            flo = {std::min(flo[0], q[0]), std::min(flo[1], q[1]), 0.0};
            fhi = {std::max(fhi[0], q[0]), std::max(fhi[1], q[1]), 0.0};
        }
        f.tri = std::make_unique<Triangulation2D>(flo, fhi);
        f.global = {-1, -1, -1};
    }

    const auto queue_face_triangles = [&](int fi, const std::vector<int>& ts) {
        for (const int t : ts) facet_queue.emplace_back(fi, t);
    };

    // Inserts a point into the tetrahedralization and into the triangulations of its faces.
    const auto insert = [&](const Point& p, const std::vector<int>& on_faces, const std::vector<int>& on_edges,
                            int hint) -> int {
        const int start = tri.locate(p, hint >= 0 ? hint : tri.last_tet());
        const auto cav = tri.cavity(p, start);
        const int v = tri.insert(p, cav);
        faces_of.push_back(on_faces);
        edges_of.push_back(on_edges);
        segments_of.emplace_back();
        for (const int t : tri.created()) tet_queue.push_back(t);
        for (const int fi : on_faces) {
            Face& f = gen.faces[static_cast<std::size_t>(fi)];
            const std::size_t before = f.tri->triangle_slots();
            const Point q = f.to_local(p);
            const int at = f.tri->locate(q, f.tri->last_triangle());
            const auto fc = f.tri->cavity(q, at);
            const int lv = f.tri->insert(q, fc);
            if (static_cast<std::size_t>(lv) >= f.global.size()) f.global.resize(static_cast<std::size_t>(lv) + 1, -1);
            f.global[static_cast<std::size_t>(lv)] = v;
            f.local[v] = lv;
            std::vector<int> fresh;
            for (const int t : fc.triangles) {
                if (f.tri->triangle(t).alive) fresh.push_back(t);
            }
            for (std::size_t t = before; t < f.tri->triangle_slots(); ++t) fresh.push_back(static_cast<int>(t));
            queue_face_triangles(fi, fresh);
        }
        // The subfacets and subsegments around the cavity may now be encroached by p.
        std::set<int> around;
        for (const int t : tri.created()) {
            for (const int x : tri.tet(t).v) {
                if (x != v && !tri.is_super(x)) around.insert(x);
            }
        }
        for (const int x : around) {
            for (const int fi : faces_of[static_cast<std::size_t>(x)]) {
                const Face& f = gen.faces[static_cast<std::size_t>(fi)];
                queue_face_triangles(fi, f.tri->incident(f.local.at(x)));
            }
            for (const int s : segments_of[static_cast<std::size_t>(x)]) segment_queue.push_back(s);
        }
        return v;
    };

    const auto add_segment = [&](int a, int b, int edge) {
        const int k = static_cast<int>(segments.size());
        segments.push_back({a, b, edge, true});
        segment_at[{std::min(a, b), std::max(a, b)}] = k;
        segments_of[static_cast<std::size_t>(a)].push_back(k);
        segments_of[static_cast<std::size_t>(b)].push_back(k);
        segment_queue.push_back(k);
    };

    // The polyhedron's vertices, in BRIO order.
    std::vector<int> vertex_of(d.vertices.size(), -1);
    for (const std::size_t i : detail::brio_order(d.vertices)) {
        std::vector<int> fs;
        for (std::size_t f = 0; f < gen.faces.size(); ++f) {
            if (std::ranges::find(gen.faces[f].vid, i) != gen.faces[f].vid.end()) fs.push_back(static_cast<int>(f));
        }
        std::vector<int> es;
        for (std::size_t e = 0; e < gen.edges.size(); ++e) {
            if (gen.edges[e].a == i || gen.edges[e].b == i) es.push_back(static_cast<int>(e));
        }
        vertex_of[i] = insert(d.vertices[i], fs, es, -1);
    }
    for (std::size_t e = 0; e < gen.edges.size(); ++e) {
        add_segment(vertex_of[gen.edges[e].a], vertex_of[gen.edges[e].b], static_cast<int>(e));
    }
    const auto is_corner = [&](int v) { return v >= 4 && static_cast<std::size_t>(v) < 4 + d.vertices.size(); };

    const auto common_face = [&](std::initializer_list<int> vs) -> int {
        const std::vector<int>& first = faces_of[static_cast<std::size_t>(*vs.begin())];
        for (const int f : first) {
            if (std::ranges::all_of(vs, [&](int x) {
                    return std::ranges::find(faces_of[static_cast<std::size_t>(x)], f) != faces_of[static_cast<std::size_t>(x)].end();
                })) {
                return f;
            }
        }
        return -1;
    };
    const auto on_one_edge = [&](std::initializer_list<int> vs) {
        const std::vector<int>& first = edges_of[static_cast<std::size_t>(*vs.begin())];
        return std::ranges::any_of(first, [&](int e) {
            return std::ranges::all_of(vs, [&](int x) {
                return std::ranges::find(edges_of[static_cast<std::size_t>(x)], e) != edges_of[static_cast<std::size_t>(x)].end();
            });
        });
    };

    // The 3D tetrahedra holding vertices a and b (and c, when c >= 0).
    const auto tets_with = [&](int a, int b, int c) {
        std::vector<int> out;
        for (const int t : tri.incident(a)) {
            const auto& v = tri.tet(t).v;
            if (std::ranges::find(v, b) == v.end()) continue;
            if (c >= 0 && std::ranges::find(v, c) == v.end()) continue;
            out.push_back(t);
        }
        return out;
    };

    const auto segment_encroached = [&](const Subsegment& s) {
        const std::vector<int> around = tets_with(s.a, s.b, -1);
        if (around.empty()) return true;
        const Point& A = tri.point(s.a);
        const Point& B = tri.point(s.b);
        for (const int t : around) {
            for (const int x : tri.tet(t).v) {
                if (x == s.a || x == s.b || tri.is_super(x)) continue;
                const Point& X = tri.point(x);
                if (detail::dot3(detail::sub(A, X), detail::sub(B, X)) < 0.0) return true;
            }
        }
        return false;
    };

    // A face-triangulation triangle that is a subfacet: alive, no super vertex, inside the face,
    // and not three vertices of one domain edge. Its 3D vertices, or nullopt.
    const auto subfacet = [&](int fi, int t) -> std::optional<std::array<int, 3>> {
        const Face& f = gen.faces[static_cast<std::size_t>(fi)];
        const auto& tr = f.tri->triangle(t);
        if (!tr.alive || f.tri->is_super(tr.v[0]) || f.tri->is_super(tr.v[1]) || f.tri->is_super(tr.v[2])) return std::nullopt;
        const std::array<int, 3> v{f.global[static_cast<std::size_t>(tr.v[0])], f.global[static_cast<std::size_t>(tr.v[1])],
                                   f.global[static_cast<std::size_t>(tr.v[2])]};
        if (on_one_edge({v[0], v[1], v[2]})) return std::nullopt;
        const Point& a = f.tri->point(tr.v[0]);
        const Point& b = f.tri->point(tr.v[1]);
        const Point& c = f.tri->point(tr.v[2]);
        if (!f.contains_local({(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, 0.0})) return std::nullopt;
        return v;
    };

    const auto facet_encroached = [&](const std::array<int, 3>& v) {
        const std::vector<int> around = tets_with(v[0], v[1], v[2]);
        if (around.empty()) return true;
        for (const int t : around) {
            for (const int x : tri.tet(t).v) {
                if (x == v[0] || x == v[1] || x == v[2] || tri.is_super(x)) continue;
                if (predicates::in_diametral_sphere(tri.point(v[0]), tri.point(v[1]), tri.point(v[2]), tri.point(x)) > 0) {
                    return true;
                }
            }
        }
        return false;
    };

    const auto node_cap = [&]() { return tri.vertex_count() - 4 >= options.max_nodes; };
    const auto exhausted = []() { return base::Error{ErrorCode::resource_exhausted, "max_nodes reached", std::nullopt}; };

    const auto split_segment = [&](int k) -> std::optional<base::Error> {
        const Subsegment s = segments[static_cast<std::size_t>(k)];
        const Point& A = tri.point(s.a);
        const Point& B = tri.point(s.b);
        const double len = dist3(A, B);
        double t = 0.5;
        const bool ca = is_corner(s.a);
        const bool cb = is_corner(s.b);
        if (ca != cb) {
            // Concentric shells around a polyhedron vertex: a power of two from it, near the middle.
            const double dd = std::exp2(std::round(std::log2(0.5 * len)));
            const double f = dd / len;
            if (f > 0.25 && f < 0.75) t = ca ? f : 1.0 - f;
        }
        Point m{};
        for (std::size_t c = 0; c < 3; ++c) m[c] = A[c] == B[c] ? A[c] : A[c] + t * (B[c] - A[c]);
        if (node_cap()) return exhausted();
        const DomainEdge& e = gen.edges[static_cast<std::size_t>(s.edge)];
        const int v = insert(m, {e.faces[0], e.faces[1]}, {s.edge}, tri.vertex_tet(s.a));
        segments[static_cast<std::size_t>(k)].alive = false;
        segment_at.erase({std::min(s.a, s.b), std::max(s.a, s.b)});
        add_segment(s.a, v, s.edge);
        add_segment(v, s.b, s.edge);
        ++report.segment_splits;
        return std::nullopt;
    };

    // Splits a subfacet at its circumcentre, unless that encroaches subsegments of its face (they
    // are split instead) or lies outside the face (the nearest subsegment beyond it is split).
    const auto split_facet = [&](int fi, const std::array<int, 3>& v) -> std::optional<base::Error> {
        const Face& f = gen.faces[static_cast<std::size_t>(fi)];
        // Still a triangle of the face (an earlier split may have removed it, and cocircular
        // subfacets share a circumcentre that must be inserted once).
        {
            const auto la = f.local.find(v[0]);
            const auto lb = f.local.find(v[1]);
            const auto lc = f.local.find(v[2]);
            if (la == f.local.end() || lb == f.local.end() || lc == f.local.end()) return std::nullopt;
            bool present = false;
            for (const auto& [x, y] : {std::pair{la->second, lb->second}, std::pair{lb->second, la->second}}) {
                const auto t = f.tri->triangle_with_edge(x, y);
                if (t && std::ranges::find(f.tri->triangle(*t).v, lc->second) != f.tri->triangle(*t).v.end()) present = true;
            }
            if (!present) return std::nullopt;
        }
        Point c = detail::triangle_circumcentre({tri.point(v[0]), tri.point(v[1]), tri.point(v[2])});
        if (f.axis >= 0) c[static_cast<std::size_t>(f.axis)] = f.plane[0][static_cast<std::size_t>(f.axis)];
        std::vector<int> hit;
        int nearest = -1;
        double nearest_d = 0.0;
        for (std::size_t k = 0; k < segments.size(); ++k) {
            const Subsegment& s = segments[k];
            if (!s.alive) continue;
            const DomainEdge& e = gen.edges[static_cast<std::size_t>(s.edge)];
            if (e.faces[0] != fi && e.faces[1] != fi) continue;
            const Point& A = tri.point(s.a);
            const Point& B = tri.point(s.b);
            if (detail::dot3(detail::sub(A, c), detail::sub(B, c)) < 0.0) hit.push_back(static_cast<int>(k));
            const Point mid = detail::scale(0.5, detail::add(A, B));
            const double dm = dist3(mid, c);
            if (nearest < 0 || dm < nearest_d) {
                nearest = static_cast<int>(k);
                nearest_d = dm;
            }
        }
        if (hit.empty() && !f.contains_local(f.to_local(c))) hit.push_back(nearest);
        if (!hit.empty()) {
            for (const int k : hit) {
                if (!segments[static_cast<std::size_t>(k)].alive) continue;
                if (auto e = split_segment(k)) return e;
            }
            return std::nullopt;
        }
        if (node_cap()) return exhausted();
        insert(c, {fi}, {}, tri.vertex_tet(v[0]));
        ++report.facet_splits;
        return std::nullopt;
    };

    for (std::size_t t = 0; t < tri.tet_slots(); ++t) {
        if (tri.tet(static_cast<int>(t)).alive) tet_queue.push_back(static_cast<int>(t));
    }
    const double bound = options.radius_edge;
    const double size_factor = std::sqrt(3.0 / 8.0);
    const std::size_t iteration_cap = 256 * options.max_nodes + 4096;
    for (std::size_t iteration = 0;; ++iteration) {
        if (iteration > iteration_cap) {
            return std::unexpected(base::Error{ErrorCode::resource_exhausted,
                                               "refinement did not settle within its iteration cap", std::nullopt});
        }
        if (!segment_queue.empty()) {
            const int k = segment_queue.front();
            segment_queue.pop_front();
            if (!segments[static_cast<std::size_t>(k)].alive || !segment_encroached(segments[static_cast<std::size_t>(k)])) continue;
            if (auto e = split_segment(k)) return std::unexpected(*e);
            continue;
        }
        if (!facet_queue.empty()) {
            const auto [fi, t] = facet_queue.front();
            facet_queue.pop_front();
            const auto v = subfacet(fi, t);
            if (!v || !facet_encroached(*v)) continue;
            if (auto e = split_facet(fi, *v)) return std::unexpected(*e);
            facet_queue.emplace_back(fi, t);  // tried again (it may have been split by a segment)
            continue;
        }
        if (tet_queue.empty()) break;
        const int t = tet_queue.front();
        tet_queue.pop_front();
        const auto& tt = tri.tet(t);
        if (!tt.alive || std::ranges::any_of(tt.v, [&](int x) { return tri.is_super(x); })) continue;
        if (common_face({tt.v[0], tt.v[1], tt.v[2], tt.v[3]}) >= 0) continue;  // flat on a face: dropped
        const std::array<Point, 4> p{tri.point(tt.v[0]), tri.point(tt.v[1]), tri.point(tt.v[2]), tri.point(tt.v[3])};
        const Point cf = circumcentre_fast(p[0], p[1], p[2], p[3]);
        const double radius = dist3(cf, p[0]);
        double shortest = dist3(p[0], p[1]);
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) shortest = std::min(shortest, dist3(p[i], p[j]));
        }
        const Point centroid = detail::scale(0.25, detail::add(detail::add(p[0], p[1]), detail::add(p[2], p[3])));
        if (!(radius > bound * shortest) && !(radius > size_factor * gen.size_at(centroid))) continue;
        const Point cc = detail::tetrahedron_circumcentre(p);
        const int start = tri.locate(cc, t);
        const auto cav = tri.cavity(cc, start);
        // A candidate that encroaches subsegments or subfacets is not inserted; they are split
        // instead, and the tetrahedron is tried again.
        std::vector<int> seg_hit;
        std::vector<std::pair<int, std::array<int, 3>>> facet_hit;
        for (const int ct : cav.tetrahedra) {
            const auto& cv = tri.tet(ct).v;
            for (std::size_t i = 0; i < 4; ++i) {
                for (std::size_t j = i + 1; j < 4; ++j) {
                    const auto it = segment_at.find({std::min(cv[i], cv[j]), std::max(cv[i], cv[j])});
                    if (it == segment_at.end()) continue;
                    const Point& A = tri.point(cv[i]);
                    const Point& B = tri.point(cv[j]);
                    if (detail::dot3(detail::sub(A, cc), detail::sub(B, cc)) < 0.0 &&
                        std::ranges::find(seg_hit, it->second) == seg_hit.end()) {
                        seg_hit.push_back(it->second);
                    }
                }
            }
            for (std::size_t i = 0; i < 4; ++i) {
                std::array<int, 3> fv{};
                for (std::size_t a = 0, m = 0; a < 4; ++a) {
                    if (a != i) fv[m++] = cv[a];
                }
                if (std::ranges::any_of(fv, [&](int x) { return tri.is_super(x); })) continue;
                const int fi = common_face({fv[0], fv[1], fv[2]});
                if (fi < 0) continue;
                if (predicates::in_diametral_sphere(tri.point(fv[0]), tri.point(fv[1]), tri.point(fv[2]), cc) > 0) {
                    std::ranges::sort(fv);
                    if (std::ranges::find(facet_hit, std::pair{fi, fv}) == facet_hit.end()) facet_hit.emplace_back(fi, fv);
                }
            }
        }
        if (!seg_hit.empty()) {
            for (const int k : seg_hit) {
                if (!segments[static_cast<std::size_t>(k)].alive) continue;
                if (auto e = split_segment(k)) return std::unexpected(*e);
            }
            tet_queue.push_back(t);
            continue;
        }
        if (!facet_hit.empty()) {
            for (const auto& [fi, fv] : facet_hit) {
                if (auto e = split_facet(fi, fv)) return std::unexpected(*e);
            }
            tet_queue.push_back(t);
            continue;
        }
        if (!gen.inside(cc)) {
            ++report.exempt_tetrahedra;
            continue;
        }
        if (node_cap()) return std::unexpected(exhausted());
        insert(cc, {}, {}, t);
        ++report.refinement_points;
    }

    // Extraction: drop the super tetrahedra and those flat on one face.
    SimplexMesh sm;
    sm.dimension = 3;
    for (std::size_t v = 4; v < tri.vertex_count(); ++v) sm.points.push_back(tri.point(static_cast<int>(v)));
    std::vector<char> kept(tri.tet_slots(), 0);
    double kept_volume = 0.0;
    double min_dihedral_all = 180.0;
    double max_ratio = 0.0;
    for (std::size_t t = 0; t < tri.tet_slots(); ++t) {
        const auto& tt = tri.tet(static_cast<int>(t));
        if (!tt.alive || std::ranges::any_of(tt.v, [&](int x) { return tri.is_super(x); })) continue;
        if (common_face({tt.v[0], tt.v[1], tt.v[2], tt.v[3]}) >= 0) {
            ++report.flat_dropped;
            continue;
        }
        kept[t] = 1;
        for (const int v : tt.v) sm.cells.push_back(static_cast<NodeId>(v - 4));
        const std::array<Point, 4> p{tri.point(tt.v[0]), tri.point(tt.v[1]), tri.point(tt.v[2]), tri.point(tt.v[3])};
        kept_volume += detail::dot3(detail::sub(p[1], p[0]), detail::cross(detail::sub(p[2], p[0]), detail::sub(p[3], p[0]))) / 6.0;
        const double dih = min_dihedral(p);
        min_dihedral_all = std::min(min_dihedral_all, dih);
        report.slivers += dih < 5.0 ? 1 : 0;
        double shortest = dist3(p[0], p[1]);
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) shortest = std::min(shortest, dist3(p[i], p[j]));
        }
        max_ratio = std::max(max_ratio, dist3(circumcentre_fast(p[0], p[1], p[2], p[3]), p[0]) / shortest);
    }
    // The boundary check: every exposed triangle on one face; the volume of the polyhedron.
    for (std::size_t f = 0; f < d.faces.size(); ++f) sm.patch_names.push_back(d.faces[f].name);
    for (std::size_t t = 0; t < tri.tet_slots(); ++t) {
        if (kept[t] == 0) continue;
        const auto& tt = tri.tet(static_cast<int>(t));
        for (std::size_t i = 0; i < 4; ++i) {
            const int nb = tt.n[i];
            if (nb >= 0 && kept[static_cast<std::size_t>(nb)] != 0) continue;
            std::array<int, 3> fv{};
            for (std::size_t a = 0, m = 0; a < 4; ++a) {
                if (a != i) fv[m++] = tt.v[a];
            }
            const int fi = common_face({fv[0], fv[1], fv[2]});
            if (fi < 0) return std::unexpected(degenerate("a boundary triangle of the generated mesh is not on one face"));
            for (const int v : fv) sm.facets.push_back(static_cast<NodeId>(v - 4));
            sm.facet_patch.push_back(fi);
        }
    }
    const double volume = gen.volume();
    if (!(std::abs(kept_volume - volume) <= 1e-12 * volume)) {
        return std::unexpected(degenerate("the generated mesh does not fill the polyhedron (volume " +
                                          std::to_string(kept_volume) + " of " + std::to_string(volume) + ")"));
    }
    DualMeshReport dual_report;
    auto dual = make_dual_mesh(sm, &dual_report);
    if (!dual) {
        return std::unexpected(base::Error{ErrorCode::degenerate_mesh,
                                           "the generated mesh failed an import check: " + dual.error().message,
                                           dual.error().context});
    }
    report.nodes = dual->node_count();
    report.tetrahedra = sm.cells.size() / 4;
    report.min_dihedral_deg = min_dihedral_all;
    report.max_radius_edge = max_ratio;
    std::vector<std::int32_t> node_region(sm.points.size(), geometry.region);
    return GeneratedMesh{std::move(*dual), std::move(sm), std::move(node_region), std::nullopt, report};
}

}  // namespace NiTCAD::mesh
