// The 3D generator over a polyhedral complex (Unit 53; DECISIONS.md FD3-FD5, FD12): Shewchuk's
// conforming Delaunay refinement of a non-convex domain, extending generate3d.cpp's convex
// generator, which stays as it is (FD12).
//
// The boundary is grouped into faces: maximal connected sets of coplanar boundary facets. A face
// may be non-convex and have holes; each keeps a 2D Delaunay triangulation of its nodes in face
// coordinates, as in the convex generator, and its subfacets are the triangles of that
// triangulation inside the face. Domain edges are the chains of facet edges where two faces meet;
// patch edges are the facet edges inside a face between facets of different patches. Both are
// subsegments. A subsegment is encroached when it is missing from the tetrahedralization or from
// the triangulation of a face it bounds, or a vertex lies strictly inside its diametral ball, so
// a missing boundary is recovered by splitting (Steiner points), never by constrained faces.
// Extraction: a flood fill from the bounding tetrahedra across triangles that are not subfacets;
// what it reaches is outside the domain.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "complex_index.hpp"
#include "complex_interfaces.hpp"
#include "dual_geometry.hpp"
#include "predicates.hpp"
#include "triangulation2d.hpp"
#include "triangulation3d.hpp"
#include "vec_ops.hpp"

namespace NiTCAD::mesh {

namespace {

using detail::ComplexIndex;
using detail::FacetKind;
using detail::Triangulation2D;
using detail::Triangulation3D;
using ErrorCode = base::ErrorCode;

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}
base::Error degenerate(std::string message) { return {ErrorCode::degenerate_mesh, std::move(message), std::nullopt}; }
base::Error exhausted() { return {ErrorCode::resource_exhausted, "max_nodes reached", std::nullopt}; }

double dist3(const Point& a, const Point& b) { return std::sqrt(detail::dot3(detail::sub(a, b), detail::sub(a, b))); }

bool collinear(const Point& a, const Point& b, const Point& c) {
    for (int k = 0; k < 3; ++k) {
        const auto i = static_cast<std::size_t>((k + 1) % 3);
        const auto j = static_cast<std::size_t>((k + 2) % 3);
        if (predicates::orient2d({a[i], a[j], 0.0}, {b[i], b[j], 0.0}, {c[i], c[j], 0.0}) != 0) return false;
    }
    return true;
}

bool in_triangle_2d(const Point& a, const Point& b, const Point& c, const Point& x) {
    const int o1 = predicates::orient2d(a, b, x);
    const int o2 = predicates::orient2d(b, c, x);
    const int o3 = predicates::orient2d(c, a, x);
    return (o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0);
}

Point circumcentre_fast(const Point& a, const Point& b, const Point& c, const Point& d) {
    const Point u = detail::sub(b, a);
    const Point v = detail::sub(c, a);
    const Point w = detail::sub(d, a);
    const double uu = detail::dot3(u, u);
    const double vv = detail::dot3(v, v);
    const double ww = detail::dot3(w, w);
    const Point vw = detail::cross(v, w);
    const Point wu = detail::cross(w, u);
    const Point uv = detail::cross(u, v);
    const double den = 2.0 * detail::dot3(u, vw);
    return detail::add(a, detail::scale(1.0 / den, detail::add(detail::add(detail::scale(uu, vw), detail::scale(vv, wu)), detail::scale(ww, uv))));
}

double min_dihedral(const std::array<Point, 4>& p) {
    double best = 180.0;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = i + 1; j < 4; ++j) {
            std::array<std::size_t, 2> o{};
            for (std::size_t k = 0, m = 0; k < 4; ++k) {
                if (k != i && k != j) o[m++] = k;
            }
            const Point e = detail::sub(p[j], p[i]);
            const Point n1 = detail::cross(e, detail::sub(p[o[0]], p[i]));
            const Point n2 = detail::cross(e, detail::sub(p[o[1]], p[i]));
            const double c = detail::dot3(n1, n2) / std::sqrt(detail::dot3(n1, n1) * detail::dot3(n2, n2));
            best = std::min(best, std::acos(std::clamp(c, -1.0, 1.0)) * 180.0 / 3.14159265358979323846);
        }
    }
    return best;
}

struct Face {
    std::vector<std::size_t> facets;
    std::array<Point, 3> plane{};
    int axis = -1;               // the coordinate constant on the face, or -1 (oblique)
    Point origin{}, u{}, w{};    // an orthonormal frame in the plane (oblique faces)
    std::vector<std::array<Point, 3>> local_tris;  // its facets' triangles in face coordinates
    std::vector<int> local_patch;                  // the patch of each
    double area = 0.0;
    std::unique_ptr<Triangulation2D> tri;
    std::vector<int> global;
    std::map<int, int> local;

    [[nodiscard]] Point to_local(const Point& p) const {
        if (axis >= 0) {
            const auto a = static_cast<std::size_t>((axis + 1) % 3);
            const auto b = static_cast<std::size_t>((axis + 2) % 3);
            return {p[a], p[b], 0.0};
        }
        const Point r = detail::sub(p, origin);
        return {detail::dot3(r, u), detail::dot3(r, w), 0.0};
    }
    // The facet triangle holding a face-coordinate point (closed), or -1.
    [[nodiscard]] int triangle_at(const Point& q) const {
        for (std::size_t k = 0; k < local_tris.size(); ++k) {
            if (in_triangle_2d(local_tris[k][0], local_tris[k][1], local_tris[k][2], q)) return static_cast<int>(k);
        }
        return -1;
    }
    [[nodiscard]] bool contains_local(const Point& q) const { return triangle_at(q) >= 0; }
};

struct DomainEdge {
    std::vector<std::size_t> chain;  // complex vertices along the edge, in order
    std::array<int, 2> faces{-1, -1};
    int face = -1;  // a patch edge: its face (both entries of `faces` are it)
};

struct Subsegment {
    int a = 0;
    int b = 0;
    int edge = 0;
    bool alive = true;
    int face = -1;        // a patch edge: its face, judged in face coordinates
    bool locked = false;  // across an interface crossing: never split (its midpoint is on the interface)
};

struct Trace {
    detail::InterfaceTrace geometry;
    std::vector<double> params;  // the pair positions along it
};

}  // namespace

std::expected<GeneratedMesh, base::Error> generate_mesh(const PolyhedralGeometry& geometry, const MeshOptions3D& options) {
    if (!(std::isfinite(options.max_size) && options.max_size > 0.0)) return std::unexpected(invalid("max_size must be finite and positive"));
    for (std::size_t k = 0; k < options.boxes.size(); ++k) {
        const SizeBox& b = options.boxes[k];
        if (!(std::isfinite(b.max_size) && b.max_size > 0.0) || !(b.lo[0] <= b.hi[0] && b.lo[1] <= b.hi[1] && b.lo[2] <= b.hi[2])) {
            return std::unexpected(invalid("size box: max_size must be finite and positive, lo <= hi", k));
        }
    }
    if (!(std::isfinite(options.radius_edge) && options.radius_edge >= 1.414)) return std::unexpected(invalid("radius_edge must be finite and >= 1.414"));
    for (const double h : options.region_max_size) {
        if (!(std::isfinite(h) && h >= 0.0)) return std::unexpected(invalid("region_max_size entries must be finite and >= 0"));
    }
    if (options.max_nodes < 4) return std::unexpected(invalid("max_nodes must be at least 4"));
    auto made = ComplexIndex::create(geometry);
    if (!made) return std::unexpected(std::move(made.error()));
    const ComplexIndex& cx = *made;
    const PolyhedralGeometry& g = cx.geometry();
    const std::vector<Point>& P = g.vertices;

    // The boundary's faces: coplanar boundary facets joined across shared edges.
    std::vector<std::size_t> boundary;
    for (std::size_t f = 0; f < cx.facets().size(); ++f) {
        if (cx.facets()[f].kind == FacetKind::boundary) boundary.push_back(f);
    }
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> edge_facets;
    for (const std::size_t f : boundary) {
        const auto& v = g.facets[f].vertices;
        for (std::size_t i = 0; i < v.size(); ++i) {
            const std::size_t a = v[i];
            const std::size_t b = v[(i + 1) % v.size()];
            edge_facets[{std::min(a, b), std::max(a, b)}].push_back(f);
        }
    }
    for (const auto& [e, fs] : edge_facets) {
        if (fs.size() != 2) return std::unexpected(invalid("the domain's boundary is not a manifold surface (an edge in more than two boundary facets)", fs[0]));
    }
    std::vector<std::size_t> parent(g.facets.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto find = [&](std::size_t x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    const auto coplanar = [&](std::size_t f, std::size_t h) {
        const auto& pl = cx.facets()[f].plane;
        return std::ranges::all_of(g.facets[h].vertices, [&](std::size_t i) { return predicates::orient3d(P[pl[0]], P[pl[1]], P[pl[2]], P[i]) == 0; });
    };
    for (const auto& [e, fs] : edge_facets) {
        if (coplanar(fs[0], fs[1])) parent[find(fs[0])] = find(fs[1]);
    }
    std::vector<Face> faces;
    std::vector<int> face_of(g.facets.size(), -1);
    {
        std::map<std::size_t, int> root_face;
        for (const std::size_t f : boundary) {
            const auto [it, fresh] = root_face.emplace(find(f), static_cast<int>(faces.size()));
            if (fresh) faces.emplace_back();
            faces[static_cast<std::size_t>(it->second)].facets.push_back(f);
            face_of[f] = it->second;
        }
    }
    for (Face& face : faces) {
        const auto& pl = cx.facets()[face.facets[0]].plane;
        face.plane = {P[pl[0]], P[pl[1]], P[pl[2]]};
        std::set<std::size_t> verts;
        for (const std::size_t f : face.facets) verts.insert(g.facets[f].vertices.begin(), g.facets[f].vertices.end());
        for (int axis = 0; axis < 3 && face.axis < 0; ++axis) {
            const auto ax = static_cast<std::size_t>(axis);
            if (std::ranges::all_of(verts, [&](std::size_t i) { return P[i][ax] == face.plane[0][ax]; })) face.axis = axis;
        }
        if (face.axis < 0) {
            face.origin = face.plane[0];
            const Point e = detail::sub(face.plane[1], face.plane[0]);
            face.u = detail::scale(1.0 / std::sqrt(detail::dot3(e, e)), e);
            Point n = detail::cross(e, detail::sub(face.plane[2], face.plane[0]));
            n = detail::scale(1.0 / std::sqrt(detail::dot3(n, n)), n);
            face.w = detail::cross(n, face.u);
        }
        Point flo{}, fhi{};
        bool first = true;
        for (const std::size_t f : face.facets) {
            for (const auto& t : cx.facets()[f].triangles) {
                std::array<Point, 3> q{face.to_local(P[t[0]]), face.to_local(P[t[1]]), face.to_local(P[t[2]])};
                for (const Point& x : q) {
                    if (first) {
                        flo = x;
                        fhi = x;
                        first = false;
                    }
                    flo = {std::min(flo[0], x[0]), std::min(flo[1], x[1]), 0.0};
                    fhi = {std::max(fhi[0], x[0]), std::max(fhi[1], x[1]), 0.0};
                }
                face.local_tris.push_back(q);
                face.local_patch.push_back(cx.facets()[f].patch);
                face.area += 0.5 * std::sqrt(detail::dot3(detail::cross(detail::sub(P[t[1]], P[t[0]]), detail::sub(P[t[2]], P[t[0]])),
                                                          detail::cross(detail::sub(P[t[1]], P[t[0]]), detail::sub(P[t[2]], P[t[0]]))));
            }
        }
        face.tri = std::make_unique<Triangulation2D>(flo, fhi);
        face.global = {-1, -1, -1};
    }

    // Domain edges (faces meeting) and patch edges (patches meeting inside a face), each a chain
    // of collinear facet edges with one key; a chain continues through a vertex only when exactly
    // two pieces of its key meet there and they are collinear.
    // A patch edge on an interface trace is straddled, not a chain (FD8 a).
    std::set<std::pair<std::size_t, std::size_t>> interface_edges;
    for (std::size_t f = 0; f < cx.facets().size(); ++f) {
        if (cx.facets()[f].kind != FacetKind::interface) continue;
        const auto& v = g.facets[f].vertices;
        for (std::size_t i = 0; i < v.size(); ++i) interface_edges.insert({std::min(v[i], v[(i + 1) % v.size()]), std::max(v[i], v[(i + 1) % v.size()])});
    }
    struct PatchEdge {
        std::pair<std::size_t, std::size_t> e;
        std::array<int, 2> patches;
        int face;
        bool straddled;
    };
    std::vector<PatchEdge> patch_edges;
    std::vector<DomainEdge> edges;
    {
        std::map<std::array<int, 3>, std::vector<std::pair<std::size_t, std::size_t>>> pieces;  // key: faces, patch-edge face
        for (const auto& [e, fs] : edge_facets) {
            const int fa = face_of[fs[0]];
            const int fb = face_of[fs[1]];
            if (fa != fb) {
                pieces[{std::min(fa, fb), std::max(fa, fb), -1}].push_back(e);
            } else if (cx.facets()[fs[0]].patch != cx.facets()[fs[1]].patch) {
                const bool straddled = interface_edges.contains(e);
                patch_edges.push_back({e, {cx.facets()[fs[0]].patch, cx.facets()[fs[1]].patch}, fa, straddled});
                if (!straddled) pieces[{fa, fa, fa}].push_back(e);
            }
        }
        for (const auto& [key, list] : pieces) {
            std::map<std::size_t, std::vector<std::size_t>> at;
            for (std::size_t k = 0; k < list.size(); ++k) {
                at[list[k].first].push_back(k);
                at[list[k].second].push_back(k);
            }
            const auto other = [&](std::size_t k, std::size_t v) { return list[k].first == v ? list[k].second : list[k].first; };
            const auto continues = [&](std::size_t v) {
                const auto& ks = at[v];
                return ks.size() == 2 && collinear(P[other(ks[0], v)], P[v], P[other(ks[1], v)]);
            };
            std::vector<bool> used(list.size(), false);
            for (std::size_t k0 = 0; k0 < list.size(); ++k0) {
                if (used[k0]) continue;
                // Walk back to an end of the chain, then forward.
                std::size_t k = k0;
                std::size_t v = list[k].first;
                for (std::size_t guard = 0; continues(v) && guard < list.size(); ++guard) {
                    const auto& ks = at[v];
                    const std::size_t next = ks[0] == k ? ks[1] : ks[0];
                    if (next == k0) break;  // a closed collinear loop is impossible; guard
                    k = next;
                    v = other(k, v);
                }
                DomainEdge de;
                de.faces = {key[0], key[1]};
                de.face = key[2];
                std::size_t from = v;
                de.chain.push_back(from);
                for (;;) {
                    used[k] = true;
                    const std::size_t to = other(k, from);
                    de.chain.push_back(to);
                    if (!continues(to)) break;
                    const auto& ks = at[to];
                    const std::size_t next = ks[0] == k ? ks[1] : ks[0];
                    if (used[next]) break;
                    k = next;
                    from = to;
                }
                edges.push_back(std::move(de));
            }
        }
    }

    // Regions and their sizes.
    std::int32_t max_region = 0;
    for (const ComplexVolume& v : g.volumes) max_region = std::max(max_region, v.region);
    std::vector<double> region_h(static_cast<std::size_t>(max_region) + 1, options.max_size);
    for (std::size_t r = 0; r < region_h.size() && r < options.region_max_size.size(); ++r) {
        if (options.region_max_size[r] > 0.0) region_h[r] = std::min(region_h[r], options.region_max_size[r]);
    }
    double smallest_region_h = std::numeric_limits<double>::infinity();
    for (const ComplexVolume& v : g.volumes) smallest_region_h = std::min(smallest_region_h, region_h[static_cast<std::size_t>(v.region)]);

    // Interfaces (FD7): pieces, traces, crossings and the rules.
    detail::InterfaceLayout layout;
    {
        detail::InterfaceInput in;
        in.complex = &cx;
        in.face_of = face_of;
        for (const Face& f : faces) {
            in.face_plane.push_back(f.plane);
            auto& tris = in.face_triangles.emplace_back();
            for (const std::size_t fct : f.facets) {
                for (const auto& t : cx.facets()[fct].triangles) tris.push_back({P[t[0]], P[t[1]], P[t[2]]});
            }
        }
        for (const DomainEdge& e : edges) in.chains.push_back({e.chain, e.faces, e.face});
        in.region_size = region_h;
        in.resolution = options.interface_resolution;
        auto read = detail::read_interfaces(in);
        if (!read) return std::unexpected(std::move(read.error()));
        layout = std::move(*read);
    }
    const double res = options.interface_resolution;
    const bool has_interfaces = !layout.pieces.empty();

    // The complex's vertices that bound a face: those on a domain edge or a patch edge.
    std::vector<std::vector<int>> vertex_edges(P.size());
    for (std::size_t e = 0; e < edges.size(); ++e) {
        for (const std::size_t v : edges[e].chain) {
            auto& list = vertex_edges[v];
            if (std::ranges::find(list, static_cast<int>(e)) == list.end()) list.push_back(static_cast<int>(e));
        }
    }
    std::vector<std::vector<int>> vertex_faces(P.size());
    for (std::size_t fi = 0; fi < faces.size(); ++fi) {
        for (const std::size_t f : faces[fi].facets) {
            for (const std::size_t v : g.facets[f].vertices) {
                auto& list = vertex_faces[v];
                if (std::ranges::find(list, static_cast<int>(fi)) == list.end()) list.push_back(static_cast<int>(fi));
            }
        }
    }
    std::vector<Point> corner_points;
    std::vector<std::size_t> corner_ids;
    for (std::size_t v = 0; v < P.size(); ++v) {
        // A vertex on an interface (a crossing) is never a node: nodes at +- s/2 replace it.
        if (!vertex_edges[v].empty() && layout.on_interface[v] == 0) {
            corner_points.push_back(P[v]);
            corner_ids.push_back(v);
        }
    }
    if (corner_points.size() > options.max_nodes) return std::unexpected(exhausted());

    Point lo = corner_points[0];
    Point hi = corner_points[0];
    for (const Point& p : corner_points) {
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    }
    Triangulation3D tri(lo, hi);
    std::vector<std::vector<int>> faces_of(4);
    std::vector<std::vector<int>> edges_of(4);
    std::vector<std::vector<int>> segments_of(4);
    std::deque<int> tet_queue;
    std::deque<std::pair<int, int>> facet_queue;
    std::deque<int> segment_queue;
    std::vector<Subsegment> segments;
    std::map<std::pair<int, int>, int> segment_at;
    MeshReport report;

    const auto insert = [&](const Point& p, const std::vector<int>& on_faces, const std::vector<int>& on_edges, int hint) -> int {
        const int start = tri.locate(p, hint >= 0 ? hint : tri.last_tet());
        const auto cav = tri.cavity(p, start);
        const int v = tri.insert(p, cav);
        faces_of.push_back(on_faces);
        edges_of.push_back(on_edges);
        segments_of.emplace_back();
        for (const int t : tri.created()) tet_queue.push_back(t);
        for (const int fi : on_faces) {
            Face& f = faces[static_cast<std::size_t>(fi)];
            const std::size_t before = f.tri->triangle_slots();
            const Point q = f.to_local(p);
            const int at = f.tri->locate(q, f.tri->last_triangle());
            const auto fc = f.tri->cavity(q, at);
            const int lv = f.tri->insert(q, fc);
            if (static_cast<std::size_t>(lv) >= f.global.size()) f.global.resize(static_cast<std::size_t>(lv) + 1, -1);
            f.global[static_cast<std::size_t>(lv)] = v;
            f.local[v] = lv;
            for (const int t : fc.triangles) {
                if (f.tri->triangle(t).alive) facet_queue.emplace_back(fi, t);
            }
            for (std::size_t t = before; t < f.tri->triangle_slots(); ++t) facet_queue.emplace_back(fi, static_cast<int>(t));
        }
        std::set<int> around;
        for (const int t : tri.created()) {
            for (const int x : tri.tet(t).v) {
                if (x != v && !tri.is_super(x)) around.insert(x);
            }
        }
        for (const int x : around) {
            for (const int fi : faces_of[static_cast<std::size_t>(x)]) {
                const Face& f = faces[static_cast<std::size_t>(fi)];
                for (const int t : f.tri->incident(f.local.at(x))) facet_queue.emplace_back(fi, t);
            }
            for (const int s : segments_of[static_cast<std::size_t>(x)]) segment_queue.push_back(s);
        }
        // A face's triangulation changes the presence of its subsegments in it.
        for (const int fi : on_faces) {
            const Face& f = faces[static_cast<std::size_t>(fi)];
            for (const int t : f.tri->incident(f.local.at(v))) {
                for (const int lx : f.tri->triangle(t).v) {
                    if (f.tri->is_super(lx)) continue;
                    for (const int s : segments_of[static_cast<std::size_t>(f.global[static_cast<std::size_t>(lx)])]) segment_queue.push_back(s);
                }
            }
        }
        return v;
    };
    const auto add_segment = [&](int a, int b, int edge, int face) {
        const int k = static_cast<int>(segments.size());
        segments.push_back({a, b, edge, true, face});
        segment_at[{std::min(a, b), std::max(a, b)}] = k;
        segments_of[static_cast<std::size_t>(a)].push_back(k);
        segments_of[static_cast<std::size_t>(b)].push_back(k);
        segment_queue.push_back(k);
    };

    const auto node_cap = [&]() { return tri.vertex_count() - 4 >= options.max_nodes; };
    std::vector<int> vertex_of(P.size(), -1);
    for (const std::size_t i : detail::brio_order(corner_points)) {
        const std::size_t v = corner_ids[i];
        vertex_of[v] = insert(P[v], vertex_faces[v], vertex_edges[v], -1);
    }
    const int corner_end = static_cast<int>(tri.vertex_count());

    // Interface nodes (M7, FD7): no node on an interface and none at its feature lines or points.
    std::vector<char> interface_node;
    const auto mark_interface = [&](int v) {
        if (static_cast<std::size_t>(v) >= interface_node.size()) interface_node.resize(static_cast<std::size_t>(v) + 1, 0);
        interface_node[static_cast<std::size_t>(v)] = 1;
    };
    const auto snap = [&](const Face& f, Point p) {
        if (f.axis >= 0) p[static_cast<std::size_t>(f.axis)] = f.plane[0][static_cast<std::size_t>(f.axis)];
        return p;
    };
    // On each crossed domain edge or patch edge: nodes at +- s/2 from the crossing, toward the
    // chain's start and end (a patch edge ending on a trace, FD8 b: only the side it runs to, the
    // other -1).
    std::map<std::pair<int, std::size_t>, std::array<int, 2>> crossing_nodes;
    for (const detail::InterfaceCrossing& c : layout.crossings) {
        const DomainEdge& de = edges[static_cast<std::size_t>(c.chain)];
        const Point& X = P[de.chain[c.position]];
        const Point& A = P[de.chain.front()];
        const Point& B = P[de.chain.back()];
        const Point ed = detail::scale(1.0 / dist3(A, B), detail::sub(B, A));
        const std::vector<int> on_faces = de.face >= 0 ? std::vector<int>{de.face} : std::vector<int>{de.faces[0], de.faces[1]};
        std::array<int, 2> ends{-1, -1};
        for (std::size_t side = 0; side < 2; ++side) {
            if ((side == 0 && c.position == 0) || (side == 1 && c.position + 1 == de.chain.size())) continue;
            Point x = detail::add(X, detail::scale(side == 0 ? -0.5 * res : 0.5 * res, ed));
            for (std::size_t a = 0; a < 3; ++a) {
                if (A[a] == B[a]) x[a] = A[a];
            }
            if (node_cap()) return std::unexpected(exhausted());
            ends[side] = insert(x, on_faces, {c.chain}, -1);
            mark_interface(ends[side]);
        }
        crossing_nodes[{c.chain, c.position}] = ends;
    }
    // Along each trace: pairs at +- s/2 within the face, graded from each end. Where a patch edge
    // meets the trace (FD8 b), its nodes stand in for a pair: no lattice pair within 0.75 s of it,
    // and its position is kept so that refinement adds none within s / 2.
    std::vector<Trace> traces;
    for (const detail::InterfaceTrace& tg : layout.traces) {
        Trace tr{tg, detail::trace_positions(tg.start, tg.dir, tg.length, res)};
        std::vector<double> met;
        for (const detail::InterfaceCrossing& c : layout.crossings) {
            const DomainEdge& de = edges[static_cast<std::size_t>(c.chain)];
            if (de.face != tg.face) continue;
            const Point r = detail::sub(P[de.chain[c.position]], tg.start);
            const double t = detail::dot3(r, tg.dir);
            if (t > 0.0 && t < tg.length && dist3(detail::add(tg.start, detail::scale(t, tg.dir)), P[de.chain[c.position]]) <= 1e-9 * tg.length) {
                met.push_back(t);
            }
        }
        std::erase_if(tr.params, [&](double p) { return std::ranges::any_of(met, [&](double t) { return std::abs(p - t) < 0.75 * res; }); });
        const Face& f = faces[static_cast<std::size_t>(tg.face)];
        for (const double t : tr.params) {
            const Point c = detail::add(tg.start, detail::scale(t, tg.dir));
            for (const double sign : {-0.5, 0.5}) {
                if (node_cap()) return std::unexpected(exhausted());
                mark_interface(insert(snap(f, detail::add(c, detail::scale(sign * res, tg.normal))), {tg.face}, {}, -1));
            }
        }
        tr.params.insert(tr.params.end(), met.begin(), met.end());
        traces.push_back(std::move(tr));
    }
    // The pair set P of each piece, mirrored at +- s/2.
    for (const detail::InterfacePieceLayout& pc : layout.pieces) {
        const auto pairs = detail::pair_set(pc, res, options.max_nodes);
        if (!pairs) return std::unexpected(invalid("the pair set of an interface could not be meshed: " + pairs.error().message, pc.facets[0]));
        for (const Point& a : *pairs) {
            for (const double sign : {-0.5, 0.5}) {
                if (node_cap()) return std::unexpected(exhausted());
                mark_interface(insert(pc.to_world(a, sign * res), {}, {}, -1));
            }
        }
    }
    // The subsegments along each chain; the one across a crossing is locked.
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto& c = edges[e].chain;
        std::vector<int> nodes;
        std::vector<std::size_t> locked_at;
        for (std::size_t k = 0; k < c.size(); ++k) {
            const auto it = crossing_nodes.find({static_cast<int>(e), k});
            if (it == crossing_nodes.end()) {
                nodes.push_back(vertex_of[c[k]]);
                continue;
            }
            if (it->second[0] >= 0 && it->second[1] >= 0) locked_at.push_back(nodes.size());
            for (const int v : it->second) {
                if (v >= 0) nodes.push_back(v);
            }
        }
        const std::size_t first = segments.size();
        for (std::size_t k = 0; k + 1 < nodes.size(); ++k) add_segment(nodes[k], nodes[k + 1], static_cast<int>(e), edges[e].face);
        for (const std::size_t k : locked_at) segments[first + k].locked = true;
    }
    const auto is_interface = [&](int v) {
        return static_cast<std::size_t>(v) < interface_node.size() && interface_node[static_cast<std::size_t>(v)] != 0;
    };
    const auto in_band = [&](const Point& x) {
        return std::ranges::any_of(layout.pieces, [&](const detail::InterfacePieceLayout& pc) { return pc.distance(x) < layout.band; });
    };
    const auto is_corner = [&](int v) { return v >= 4 && v < corner_end; };
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
    // Whether the face triangulation holds the edge ab (3D vertices) with an apex inside the face
    // on at least one side; for a patch edge, also no apex strictly inside its diametral circle.
    const auto in_face_triangulation = [&](const Face& f, int a, int b, bool judge_circle) -> bool {
        const auto la = f.local.find(a);
        const auto lb = f.local.find(b);
        if (la == f.local.end() || lb == f.local.end()) return false;
        bool present = false;
        for (const auto& [x, y] : {std::pair{la->second, lb->second}, std::pair{lb->second, la->second}}) {
            const auto t = f.tri->triangle_with_edge(x, y);
            if (!t) continue;
            for (const int apex : f.tri->triangle(*t).v) {
                if (apex == la->second || apex == lb->second || f.tri->is_super(apex)) continue;
                const Point& A = f.tri->point(la->second);
                const Point& B = f.tri->point(lb->second);
                const Point& C = f.tri->point(apex);
                if (f.contains_local({(A[0] + B[0] + C[0]) / 3.0, (A[1] + B[1] + C[1]) / 3.0, 0.0})) present = true;
                if (judge_circle && predicates::in_diametral_circle(A, B, C) > 0) return false;
            }
        }
        return present;
    };
    const auto segment_encroached = [&](const Subsegment& s) {
        if (s.locked) return false;
        if (s.face >= 0) return !in_face_triangulation(faces[static_cast<std::size_t>(s.face)], s.a, s.b, true);
        const DomainEdge& de = edges[static_cast<std::size_t>(s.edge)];
        for (const int fi : de.faces) {
            if (!in_face_triangulation(faces[static_cast<std::size_t>(fi)], s.a, s.b, false)) return true;
        }
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
    const auto subfacet = [&](int fi, int t) -> std::optional<std::array<int, 3>> {
        const Face& f = faces[static_cast<std::size_t>(fi)];
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
                if (predicates::in_diametral_sphere(tri.point(v[0]), tri.point(v[1]), tri.point(v[2]), tri.point(x)) > 0) return true;
            }
        }
        return false;
    };

    const auto split_segment = [&](int k) -> std::optional<base::Error> {
        const Subsegment s = segments[static_cast<std::size_t>(k)];
        if (s.locked) return std::nullopt;
        const Point& A = tri.point(s.a);
        const Point& B = tri.point(s.b);
        const double len = dist3(A, B);
        double t = 0.5;
        const bool ca = is_corner(s.a);
        const bool cb = is_corner(s.b);
        if (ca != cb) {
            const double dd = std::exp2(std::round(std::log2(0.5 * len)));
            const double f = dd / len;
            if (f > 0.25 && f < 0.75) t = ca ? f : 1.0 - f;
        }
        Point m{};
        for (std::size_t c = 0; c < 3; ++c) m[c] = A[c] == B[c] ? A[c] : A[c] + t * (B[c] - A[c]);
        if (node_cap()) return exhausted();
        const DomainEdge& de = edges[static_cast<std::size_t>(s.edge)];
        std::vector<int> on_faces{de.faces[0]};
        if (de.faces[1] != de.faces[0]) on_faces.push_back(de.faces[1]);
        const int v = insert(m, on_faces, {s.edge}, tri.vertex_tet(s.a));
        segments[static_cast<std::size_t>(k)].alive = false;
        segment_at.erase({std::min(s.a, s.b), std::max(s.a, s.b)});
        add_segment(s.a, v, s.edge, s.face);
        add_segment(v, s.b, s.edge, s.face);
        ++report.segment_splits;
        return std::nullopt;
    };
    const auto split_facet = [&](int fi, const std::array<int, 3>& v) -> std::optional<base::Error> {
        const Face& f = faces[static_cast<std::size_t>(fi)];
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
        if (has_interfaces && in_band(c)) {
            // Near a trace: a new trace pair at the projection instead (M11), unless one is already
            // within s / 2 along the trace; then c itself, judged by P4.
            Trace* best = nullptr;
            double best_d = 0.0;
            for (Trace& tr : traces) {
                if (tr.geometry.face != fi) continue;
                const Point r = detail::sub(c, tr.geometry.start);
                const double along = std::clamp(detail::dot3(r, tr.geometry.dir), 0.0, tr.geometry.length);
                const double dd = dist3(c, detail::add(tr.geometry.start, detail::scale(along, tr.geometry.dir)));
                if (best == nullptr || dd < best_d) {
                    best = &tr;
                    best_d = dd;
                }
            }
            if (best != nullptr && best_d < layout.band) {
                const auto& tg = best->geometry;
                const double t = std::clamp(detail::dot3(detail::sub(c, tg.start), tg.dir), 1.5 * res, tg.length - 1.5 * res);
                const bool near = std::ranges::any_of(best->params, [&](double x) { return std::abs(x - t) < 0.5 * res; });
                if (!near) {
                    if (node_cap()) return exhausted();
                    const Point centre = detail::add(tg.start, detail::scale(t, tg.dir));
                    for (const double sign : {-0.5, 0.5}) {
                        mark_interface(insert(snap(f, detail::add(centre, detail::scale(sign * res, tg.normal))), {fi}, {}, tri.vertex_tet(v[0])));
                    }
                    best->params.push_back(t);
                    ++report.facet_splits;
                    return std::nullopt;
                }
            }
            ++report.band_facet_points;
        }
        const Point cl = f.to_local(c);
        std::vector<int> hit;
        int nearest = -1;
        double nearest_d = 0.0;
        for (std::size_t k = 0; k < segments.size(); ++k) {
            const Subsegment& s = segments[k];
            if (!s.alive || s.locked) continue;
            const DomainEdge& de = edges[static_cast<std::size_t>(s.edge)];
            if (de.faces[0] != fi && de.faces[1] != fi) continue;
            const Point& A = tri.point(s.a);
            const Point& B = tri.point(s.b);
            if (s.face >= 0) {
                if (predicates::in_diametral_circle(f.tri->point(f.local.at(s.a)), f.tri->point(f.local.at(s.b)), cl) > 0) {
                    hit.push_back(static_cast<int>(k));
                }
            } else if (detail::dot3(detail::sub(A, c), detail::sub(B, c)) < 0.0) {
                hit.push_back(static_cast<int>(k));
            }
            const double dm = dist3(detail::scale(0.5, detail::add(A, B)), c);
            if (nearest < 0 || dm < nearest_d) {
                nearest = static_cast<int>(k);
                nearest_d = dm;
            }
        }
        if (hit.empty() && !f.contains_local(cl)) hit.push_back(nearest);
        if (!hit.empty()) {
            for (const int k : hit) {
                if (k < 0 || !segments[static_cast<std::size_t>(k)].alive) continue;
                if (auto e = split_segment(k)) return e;
            }
            return std::nullopt;
        }
        if (node_cap()) return exhausted();
        insert(c, {fi}, {}, tri.vertex_tet(v[0]));
        ++report.facet_splits;
        return std::nullopt;
    };
    // h at p for a region size h_r (the boxes may lower it).
    const auto size_at = [&](const Point& p, double h_r) {
        double h = h_r;
        for (const SizeBox& b : options.boxes) {
            if (p[0] >= b.lo[0] && p[0] <= b.hi[0] && p[1] >= b.lo[1] && p[1] <= b.hi[1] && p[2] >= b.lo[2] && p[2] <= b.hi[2]) {
                h = std::min(h, b.max_size);
            }
        }
        return h;
    };
    // The region size at a point strictly inside the domain or on an internal facet (the smaller
    // of the two regions on an interface), or nothing outside or on the boundary.
    const auto domain_size = [&](const Point& p) -> std::optional<double> {
        const auto l = cx.locate(p);
        const auto h_of_volume = [&](int vol) { return region_h[static_cast<std::size_t>(g.volumes[static_cast<std::size_t>(vol)].region)]; };
        if (l.kind == ComplexIndex::Location::Kind::inside) return h_of_volume(l.index);
        if (l.kind != ComplexIndex::Location::Kind::on_facet) return std::nullopt;
        const auto& fc = cx.facets()[static_cast<std::size_t>(l.index)];
        if (fc.kind == FacetKind::boundary) return std::nullopt;
        return std::min(h_of_volume(fc.volumes[0]), h_of_volume(fc.volumes[1]));
    };

    for (std::size_t t = 0; t < tri.tet_slots(); ++t) {
        if (tri.tet(static_cast<int>(t)).alive) tet_queue.push_back(static_cast<int>(t));
    }
    const double bound = options.radius_edge;
    const double size_factor = std::sqrt(3.0 / 8.0);
    const std::size_t iteration_cap = 256 * options.max_nodes + 4096;
    for (std::size_t iteration = 0;; ++iteration) {
        if (iteration > iteration_cap) {
            return std::unexpected(base::Error{ErrorCode::resource_exhausted, "refinement did not settle within its iteration cap", std::nullopt});
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
            facet_queue.emplace_back(fi, t);
            continue;
        }
        if (tet_queue.empty()) break;
        const int t = tet_queue.front();
        tet_queue.pop_front();
        const auto& tt = tri.tet(t);
        if (!tt.alive || std::ranges::any_of(tt.v, [&](int x) { return tri.is_super(x); })) continue;
        if (common_face({tt.v[0], tt.v[1], tt.v[2], tt.v[3]}) >= 0) continue;
        const std::array<Point, 4> p{tri.point(tt.v[0]), tri.point(tt.v[1]), tri.point(tt.v[2]), tri.point(tt.v[3])};
        const Point cf = circumcentre_fast(p[0], p[1], p[2], p[3]);
        const double radius = dist3(cf, p[0]);
        double shortest = dist3(p[0], p[1]);
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) shortest = std::min(shortest, dist3(p[i], p[j]));
        }
        const Point centroid = detail::scale(0.25, detail::add(detail::add(p[0], p[1]), detail::add(p[2], p[3])));
        const bool poor = radius > bound * shortest;
        if (!poor && !(radius > size_factor * size_at(centroid, smallest_region_h))) continue;
        // Only tetrahedra inside the domain are refined (outside ones are dropped at extraction).
        const auto h_here = domain_size(centroid);
        if (!h_here) continue;
        if (!poor && !(radius > size_factor * size_at(centroid, *h_here))) continue;
        if (std::ranges::all_of(tt.v, [&](int x) { return is_interface(x); })) {
            ++report.exempt_tetrahedra;  // pair nodes only (M11)
            continue;
        }
        const Point cc = detail::tetrahedron_circumcentre(p);
        if (has_interfaces && in_band(cc)) {
            ++report.exempt_tetrahedra;  // no circumcentre in a protected band (M7)
            continue;
        }
        const int start = tri.locate(cc, t);
        const auto cav = tri.cavity(cc, start);
        std::vector<int> seg_hit;
        std::vector<std::pair<int, std::array<int, 3>>> facet_hit;
        for (const int ct : cav.tetrahedra) {
            const auto& cv = tri.tet(ct).v;
            for (std::size_t i = 0; i < 4; ++i) {
                for (std::size_t j = i + 1; j < 4; ++j) {
                    const auto it = segment_at.find({std::min(cv[i], cv[j]), std::max(cv[i], cv[j])});
                    if (it == segment_at.end() || segments[static_cast<std::size_t>(it->second)].locked) continue;
                    const Point& A = tri.point(cv[i]);
                    const Point& B = tri.point(cv[j]);
                    if (detail::dot3(detail::sub(A, cc), detail::sub(B, cc)) < 0.0 && std::ranges::find(seg_hit, it->second) == seg_hit.end()) {
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
        if (!domain_size(cc)) {
            ++report.exempt_tetrahedra;
            continue;
        }
        if (node_cap()) return std::unexpected(exhausted());
        insert(cc, {}, {}, t);
        ++report.refinement_points;
    }

    // Extraction: flood from the bounding tetrahedra across triangles that are not subfacets.
    const auto is_subfacet = [&](const std::array<int, 3>& fv) {
        if (std::ranges::any_of(fv, [&](int x) { return tri.is_super(x); })) return false;
        const int fi = common_face({fv[0], fv[1], fv[2]});
        if (fi < 0 || on_one_edge({fv[0], fv[1], fv[2]})) return false;
        const Face& f = faces[static_cast<std::size_t>(fi)];
        Point c{0, 0, 0};
        for (const int v : fv) {
            const Point q = f.to_local(tri.point(v));
            c = {c[0] + q[0] / 3.0, c[1] + q[1] / 3.0, 0.0};
        }
        return f.contains_local(c);
    };
    std::vector<char> outside(tri.tet_slots(), 0);
    std::vector<char> flat(tri.tet_slots(), 0);
    std::vector<int> stack;
    for (std::size_t t = 0; t < tri.tet_slots(); ++t) {
        const auto& tt = tri.tet(static_cast<int>(t));
        if (!tt.alive) continue;
        if (std::ranges::any_of(tt.v, [&](int x) { return tri.is_super(x); })) {
            outside[t] = 1;
            stack.push_back(static_cast<int>(t));
        } else if (common_face({tt.v[0], tt.v[1], tt.v[2], tt.v[3]}) >= 0) {
            flat[t] = 1;
        }
    }
    while (!stack.empty()) {
        const int t = stack.back();
        stack.pop_back();
        const auto& tt = tri.tet(t);
        for (std::size_t i = 0; i < 4; ++i) {
            const int nb = tt.n[i];
            if (nb < 0 || outside[static_cast<std::size_t>(nb)] != 0 || flat[static_cast<std::size_t>(nb)] != 0) continue;
            std::array<int, 3> fv{};
            for (std::size_t a = 0, m = 0; a < 4; ++a) {
                if (a != i) fv[m++] = tt.v[a];
            }
            if (is_subfacet(fv)) continue;
            outside[static_cast<std::size_t>(nb)] = 1;
            stack.push_back(nb);
        }
    }
    SimplexMesh sm;
    sm.dimension = 3;
    sm.patch_names = cx.report().patches;
    std::vector<char> kept(tri.tet_slots(), 0);
    double kept_volume = 0.0;
    double min_dihedral_all = 180.0;
    double max_ratio = 0.0;
    std::vector<int> renumber(tri.vertex_count(), -1);
    for (std::size_t t = 0; t < tri.tet_slots(); ++t) {
        const auto& tt = tri.tet(static_cast<int>(t));
        if (!tt.alive || outside[t] != 0) continue;
        if (flat[t] != 0) {
            ++report.flat_dropped;
            continue;
        }
        kept[t] = 1;
        for (const int v : tt.v) renumber[static_cast<std::size_t>(v)] = 0;
    }
    for (std::size_t v = 4; v < tri.vertex_count(); ++v) {
        if (renumber[v] < 0) return std::unexpected(degenerate("a generated node is in no tetrahedron inside the domain"));
        renumber[v] = static_cast<int>(sm.points.size());
        sm.points.push_back(tri.point(static_cast<int>(v)));
    }
    for (std::size_t t = 0; t < tri.tet_slots(); ++t) {
        if (kept[t] == 0) continue;
        const auto& tt = tri.tet(static_cast<int>(t));
        for (const int v : tt.v) sm.cells.push_back(static_cast<NodeId>(renumber[static_cast<std::size_t>(v)]));
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
    std::vector<double> face_area(faces.size(), 0.0);
    // FD8 (a): on a face with a straddled patch edge, each corner takes the patch of the facets
    // holding its node, unless those differ (a node on a patch edge): then the triangle's.
    std::vector<char> straddled_face(faces.size(), 0);
    for (const PatchEdge& pe : patch_edges) {
        if (pe.straddled) straddled_face[static_cast<std::size_t>(pe.face)] = 1;
    }
    std::vector<std::int32_t> corner_patch;
    bool corners_differ = false;
    const auto patch_of_node = [&](const Face& f, const Point& p, int fallback) {
        const Point q = f.to_local(p);
        int found = -1;
        for (std::size_t k = 0; k < f.local_tris.size(); ++k) {
            if (!in_triangle_2d(f.local_tris[k][0], f.local_tris[k][1], f.local_tris[k][2], q)) continue;
            if (found >= 0 && found != f.local_patch[k]) return fallback;
            found = f.local_patch[k];
        }
        return found >= 0 ? found : fallback;
    };
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
            if (!is_subfacet(fv)) return std::unexpected(degenerate("a boundary triangle of the generated mesh is not inside one face"));
            const int fi = common_face({fv[0], fv[1], fv[2]});
            const Face& f = faces[static_cast<std::size_t>(fi)];
            Point c{0, 0, 0};
            for (const int v : fv) {
                const Point q = f.to_local(tri.point(v));
                c = {c[0] + q[0] / 3.0, c[1] + q[1] / 3.0, 0.0};
            }
            const int at = f.triangle_at(c);
            const int patch = f.local_patch[static_cast<std::size_t>(at)];
            for (const int v : fv) sm.facets.push_back(static_cast<NodeId>(renumber[static_cast<std::size_t>(v)]));
            sm.facet_patch.push_back(patch);
            for (const int v : fv) {
                const int cp = straddled_face[static_cast<std::size_t>(fi)] != 0 ? patch_of_node(f, tri.point(v), patch) : patch;
                corners_differ = corners_differ || cp != patch;
                corner_patch.push_back(cp);
            }
            const Point& A = tri.point(fv[0]);
            const Point n = detail::cross(detail::sub(tri.point(fv[1]), A), detail::sub(tri.point(fv[2]), A));
            face_area[static_cast<std::size_t>(fi)] += 0.5 * std::sqrt(detail::dot3(n, n));
        }
    }
    if (corners_differ) sm.facet_corner_patch = std::move(corner_patch);
    double volume = 0.0;
    for (const double v : cx.report().volume) volume += v;
    if (!(std::abs(kept_volume - volume) <= 1e-12 * volume)) {
        return std::unexpected(degenerate("the generated mesh does not fill the domain (volume " + std::to_string(kept_volume) + " of " +
                                          std::to_string(volume) + ")"));
    }
    for (std::size_t fi = 0; fi < faces.size(); ++fi) {
        if (!(std::abs(face_area[fi] - faces[fi].area) <= 1e-12 * faces[fi].area)) {
            return std::unexpected(degenerate("the generated mesh does not cover a boundary face (area " + std::to_string(face_area[fi]) +
                                              " of " + std::to_string(faces[fi].area) + ")"));
        }
    }
    DualMeshReport dual_report;
    auto dual = make_dual_mesh(sm, &dual_report);
    if (!dual) {
        return std::unexpected(base::Error{ErrorCode::degenerate_mesh, "the generated mesh failed an import check: " + dual.error().message,
                                           dual.error().context});
    }
    report.nodes = dual->node_count();
    report.tetrahedra = sm.cells.size() / 4;
    report.min_dihedral_deg = min_dihedral_all;
    report.max_radius_edge = max_ratio;
    // Regions (FD6). Without interfaces, the one region (slice 2). Otherwise each node located
    // exactly; a node on an interface facet is a generator bug. A node on an oblique boundary face
    // is rounded off it and may locate just outside: it takes the region beside the nearest
    // boundary facet, within 1e-9 of the domain's extent (pairs lie s / 2 from any interface).
    std::vector<std::int32_t> node_region(sm.points.size(), g.volumes[0].region);
    if (has_interfaces) {
        const double extent = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
        for (std::size_t i = 0; i < sm.points.size(); ++i) {
            const Point& p = sm.points[i];
            const auto l = cx.locate(p);
            int vol = -1;
            if (l.kind == ComplexIndex::Location::Kind::inside) {
                vol = l.index;
            } else if (l.kind == ComplexIndex::Location::Kind::on_facet) {
                const auto& fc = cx.facets()[static_cast<std::size_t>(l.index)];
                if (fc.kind == FacetKind::interface) return std::unexpected(degenerate("a generated node lies on an interface"));
                vol = fc.volumes[0];
            } else {
                double nearest = 1e-9 * extent;
                for (const detail::IndexedFacet& fc : cx.facets()) {
                    if (fc.kind != FacetKind::boundary) continue;
                    for (const auto& t : fc.triangles) {
                        const double d = detail::triangle_distance(p, P[t[0]], P[t[1]], P[t[2]]);
                        if (d <= nearest) {
                            nearest = d;
                            vol = fc.volumes[0];
                        }
                    }
                }
            }
            if (vol < 0) return std::unexpected(degenerate("a generated node could not be located in the domain"));
            node_region[i] = g.volumes[static_cast<std::size_t>(vol)].region;
        }
    }
    for (std::size_t v = 4; v < tri.vertex_count(); ++v) report.interface_nodes += is_interface(static_cast<int>(v)) ? 1 : 0;
    std::optional<InterfaceCertificate> certificate;
    if (has_interfaces) {
        InterfaceSet set{{}, res};
        for (const detail::InterfacePieceLayout& pc : layout.pieces) set.pieces.push_back({pc.polygon, pc.region_a, pc.region_b});
        // P6 (FD10) wherever a patch edge meets an interface (FD8 a or b): every patch edge
        // declared; zones where a straddled run ends or meets a chain, and where a patch edge
        // meets a trace.
        std::optional<PatchBoundarySet> patch_set;
        const bool patch_crossing = std::ranges::any_of(layout.crossings, [&](const detail::InterfaceCrossing& c) {
            return edges[static_cast<std::size_t>(c.chain)].face >= 0;
        });
        if (patch_crossing || std::ranges::any_of(patch_edges, &PatchEdge::straddled)) {
            PatchBoundarySet& ps = patch_set.emplace();
            ps.resolution = res;
            std::map<std::size_t, std::vector<std::size_t>> straddled_ends;  // vertex -> the other ends of its straddled edges
            for (const PatchEdge& pe : patch_edges) {
                ps.boundaries.push_back({P[pe.e.first], P[pe.e.second], pe.patches[0], pe.patches[1]});
                if (pe.straddled) {
                    straddled_ends[pe.e.first].push_back(pe.e.second);
                    straddled_ends[pe.e.second].push_back(pe.e.first);
                }
            }
            for (const auto& [v, ends] : straddled_ends) {
                // A run's end, a chain meeting it, or a corner of the run (two traces crossing on the face).
                bool zone = ends.size() != 2 || !vertex_edges[v].empty();
                if (!zone) {
                    const Point n = detail::cross(detail::sub(P[ends[0]], P[v]), detail::sub(P[ends[1]], P[v]));
                    const double a = std::sqrt(detail::dot3(detail::sub(P[ends[0]], P[v]), detail::sub(P[ends[0]], P[v])));
                    const double b = std::sqrt(detail::dot3(detail::sub(P[ends[1]], P[v]), detail::sub(P[ends[1]], P[v])));
                    zone = std::sqrt(detail::dot3(n, n)) > 1e-9 * a * b;
                }
                if (zone) ps.zone_centres.push_back(P[v]);
            }
            for (const detail::InterfaceCrossing& c : layout.crossings) {
                const DomainEdge& de = edges[static_cast<std::size_t>(c.chain)];
                if (de.face >= 0) ps.zone_centres.push_back(P[de.chain[c.position]]);
            }
            ps.patch_area.assign(sm.patch_names.size(), 0.0);
            for (const std::size_t f : boundary) {
                for (const auto& t : cx.facets()[f].triangles) {
                    const Point n = detail::cross(detail::sub(P[t[1]], P[t[0]]), detail::sub(P[t[2]], P[t[0]]));
                    ps.patch_area[static_cast<std::size_t>(cx.facets()[f].patch)] += 0.5 * std::sqrt(detail::dot3(n, n));
                }
            }
        }
        auto cert = certify_interfaces(sm, *dual, node_region, set, patch_set ? &*patch_set : nullptr);
        if (!cert) {
            return std::unexpected(base::Error{ErrorCode::degenerate_mesh, "the generated mesh did not meet P4: " + cert.error().message,
                                               cert.error().context});
        }
        certificate = std::move(*cert);
    }
    return GeneratedMesh{std::move(*dual), std::move(sm), std::move(node_region), std::move(certificate), report};
}

}  // namespace NiTCAD::mesh
