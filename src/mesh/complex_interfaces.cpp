#include "complex_interfaces.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>

#include "NiTCAD/mesh/generate.hpp"
#include "polygon2d.hpp"
#include "predicates.hpp"
#include "vec_ops.hpp"

namespace NiTCAD::mesh::detail {

namespace {

base::Error invalid(std::string message, std::size_t index) {
    return {base::ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}

double norm3(const Point& a) { return std::sqrt(dot3(a, a)); }
Point unit(const Point& a) { return scale(1.0 / norm3(a), a); }
double dist(const Point& a, const Point& b) { return norm3(sub(a, b)); }

// Closest distance between the segments pq and rs (Ericson, Real-Time Collision Detection, 5.1.9).
double segment_segment(const Point& p, const Point& q, const Point& r, const Point& s) {
    const Point d1 = sub(q, p);
    const Point d2 = sub(s, r);
    const Point w = sub(p, r);
    const double a = dot3(d1, d1);
    const double e = dot3(d2, d2);
    const double f = dot3(d2, w);
    double t1 = 0.0;
    double t2 = 0.0;
    if (a <= 0.0 && e <= 0.0) return dist(p, r);
    if (a <= 0.0) {
        t2 = std::clamp(f / e, 0.0, 1.0);
    } else {
        const double c = dot3(d1, w);
        if (e <= 0.0) {
            t1 = std::clamp(-c / a, 0.0, 1.0);
        } else {
            const double b = dot3(d1, d2);
            const double den = a * e - b * b;
            t1 = den > 0.0 ? std::clamp((b * f - c * e) / den, 0.0, 1.0) : 0.0;
            t2 = (b * t1 + f) / e;
            if (t2 < 0.0) {
                t2 = 0.0;
                t1 = std::clamp(-c / a, 0.0, 1.0);
            } else if (t2 > 1.0) {
                t2 = 1.0;
                t1 = std::clamp((b - c) / a, 0.0, 1.0);
            }
        }
    }
    return dist(add(p, scale(t1, d1)), add(r, scale(t2, d2)));
}

// Closest distance from x to the triangle abc (Ericson 5.1.5).
double point_triangle(const Point& x, const Point& a, const Point& b, const Point& c) {
    const Point ab = sub(b, a);
    const Point ac = sub(c, a);
    const Point ap = sub(x, a);
    const double d1 = dot3(ab, ap);
    const double d2 = dot3(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return dist(x, a);
    const Point bp = sub(x, b);
    const double d3 = dot3(ab, bp);
    const double d4 = dot3(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return dist(x, b);
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) return dist(x, add(a, scale(d1 / (d1 - d3), ab)));
    const Point cp = sub(x, c);
    const double d5 = dot3(ab, cp);
    const double d6 = dot3(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return dist(x, c);
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) return dist(x, add(a, scale(d2 / (d2 - d6), ac)));
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        return dist(x, add(b, scale((d4 - d3) / ((d4 - d3) + (d5 - d6)), sub(c, b))));
    }
    const double den = 1.0 / (va + vb + vc);
    return dist(x, add(a, add(scale(vb * den, ab), scale(vc * den, ac))));
}

// Sorts half-planes (each given by a point off the line ab) cyclically about the line, exactly.
void sort_about_line(const Point& a, const Point& b, std::vector<std::pair<Point, std::size_t>>& h) {
    if (h.empty()) return;
    const Point r = h[0].first;
    const auto half = [&](const Point& p) {
        const int o = predicates::orient3d(a, b, r, p);
        if (o != 0) return o > 0 ? 0 : 1;
        return predicates::dihedral_dot_sign(a, b, r, p) > 0 ? 0 : 1;
    };
    std::ranges::stable_sort(h, [&](const auto& x, const auto& y) {
        const int hx = half(x.first);
        const int hy = half(y.first);
        if (hx != hy) return hx < hy;
        return predicates::orient3d(a, b, x.first, y.first) > 0;
    });
}

// Sorts rays from o in a plane (q strictly off it) cyclically, exactly.
void sort_about_point(const Point& o, const Point& q, std::vector<Point>& rays) {
    if (rays.empty()) return;
    const Point r = rays[0];
    const auto half = [&](const Point& p) {
        const int s = predicates::orient3d(o, r, p, q);
        if (s != 0) return s > 0 ? 0 : 1;
        return predicates::dot_sign_3d(o, r, p) > 0 ? 0 : 1;
    };
    std::ranges::stable_sort(rays, [&](const Point& x, const Point& y) {
        const int hx = half(x);
        const int hy = half(y);
        if (hx != hy) return hx < hy;
        return predicates::orient3d(o, x, y, q) > 0;
    });
}

}  // namespace

double triangle_distance(const Point& x, const Point& a, const Point& b, const Point& c) { return point_triangle(x, a, b, c); }

Point InterfacePieceLayout::to_local(const Point& p) const {
    if (axis >= 0) return project(p, axis);
    const Point r = sub(p, origin);
    return {dot3(r, u), dot3(r, w), 0.0};
}

Point InterfacePieceLayout::to_world(const Point& l, double offset) const {
    if (axis >= 0) {
        Point x{};
        x[static_cast<std::size_t>(axis)] = polygon[0][static_cast<std::size_t>(axis)] + offset;
        x[static_cast<std::size_t>((axis + 1) % 3)] = l[0];
        x[static_cast<std::size_t>((axis + 2) % 3)] = l[1];
        return x;
    }
    return add(origin, add(add(scale(l[0], u), scale(l[1], w)), scale(offset, normal)));
}

double InterfacePieceLayout::distance(const Point& x) const {
    const double h = axis >= 0 ? x[static_cast<std::size_t>(axis)] - polygon[0][static_cast<std::size_t>(axis)]
                               : dot3(sub(x, origin), normal);
    if (inside_closed_2d(local, to_local(x))) return std::abs(h);
    double d = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < polygon.size(); ++k) {
        d = std::min(d, segment_segment(x, x, polygon[k], polygon[(k + 1) % polygon.size()]));
    }
    return d;
}

namespace {

// Positions along a segment (from start, along unit dir) at which its line's coordinate is a multiple of
// 1.5 s, at least `margin` from both ends: every piece and trace beside one feature line places its
// row of pairs at the same points along it, so near the line the nodes are an extrusion of a 2D
// cross-section (the configuration U15 measured).
std::vector<double> line_lattice(const Point& start, const Point& dir, double length, double s, double margin) {
    // The line's own coordinate: the projection on its direction with the sign fixed (the first
    // nonzero component positive), so every segment of one line agrees on it.
    Point d = dir;
    for (std::size_t c = 0; c < 3; ++c) {
        if (std::abs(d[c]) > 1e-12) {
            if (d[c] < 0.0) d = scale(-1.0, d);
            break;
        }
    }
    const double step = 1.5 * s;
    const double along = dot3(dir, d);  // +1 or -1
    const double t0 = dot3(start, d);
    const double t1 = t0 + along * length;
    const double lo = std::min(t0, t1);
    const double hi = std::max(t0, t1);
    std::vector<double> at;
    for (auto k = static_cast<long long>(std::ceil((lo + margin) / step)); static_cast<double>(k) * step <= hi - margin; ++k) {
        at.push_back((static_cast<double>(k) * step - t0) * along);
    }
    std::ranges::sort(at);
    return at;
}

}  // namespace

std::vector<double> trace_positions(const Point& start, const Point& dir, double length, double s) {
    std::vector<double> at{1.5 * s, length - 1.5 * s};
    for (const double t : line_lattice(start, dir, length, s, 2.0 * s)) at.push_back(t);
    std::ranges::sort(at);
    at.erase(std::unique(at.begin(), at.end(), [&](double x, double y) { return y - x < 0.25 * s; }), at.end());
    return at;
}

std::expected<std::vector<Point>, base::Error> pair_set(const InterfacePieceLayout& piece, double s, std::size_t max_nodes) {
    // The inset polygon with its edges' graded points as extra (near-collinear) vertices.
    std::vector<Point> poly;
    const std::size_t n = piece.inset.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = piece.inset[i];
        const Point& b = piece.inset[(i + 1) % n];
        poly.push_back(a);
        const double len = std::hypot(b[0] - a[0], b[1] - a[1]);
        // The feature line's lattice, carried to the inset edge parallel to it.
        const Point& A = piece.polygon[i];
        const Point& B = piece.polygon[(i + 1) % piece.polygon.size()];
        const Point dir = unit(sub(B, A));
        const Point aw = piece.to_world(a, 0.0);
        for (const double t : line_lattice(aw, dir, len, s, 0.5 * s)) {
            Point x{a[0] + (b[0] - a[0]) * (t / len), a[1] + (b[1] - a[1]) * (t / len), 0.0};
            for (std::size_t c = 0; c < 2; ++c) {
                if (a[c] == b[c]) x[c] = a[c];
            }
            poly.push_back(x);
        }
    }
    MeshOptions2D o2;
    o2.max_size = piece.spacing;
    o2.max_nodes = max_nodes;
    auto m = generate_mesh(Geometry2D{{RegionPolygon{poly, 0}}, {}}, o2);
    if (!m) return std::unexpected(std::move(m.error()));
    return std::move(m->simplices.points);
}

std::expected<InterfaceLayout, base::Error> read_interfaces(const InterfaceInput& in) {
    const ComplexIndex& cx = *in.complex;
    const PolyhedralGeometry& g = cx.geometry();
    const std::vector<Point>& P = g.vertices;
    const auto& F = cx.facets();
    InterfaceLayout out;
    out.on_interface.assign(P.size(), 0);
    std::vector<std::size_t> interface_facets;
    for (std::size_t f = 0; f < F.size(); ++f) {
        if (F[f].kind != FacetKind::interface) continue;
        interface_facets.push_back(f);
        for (const std::size_t v : g.facets[f].vertices) out.on_interface[v] = 1;
    }
    if (interface_facets.empty()) return out;
    const double s = in.resolution;
    if (!(std::isfinite(s) && s > 0.0)) return std::unexpected(invalid("interface_resolution must be finite and positive with interfaces", interface_facets[0]));

    using EdgeKey = std::pair<std::size_t, std::size_t>;
    const auto key = [](std::size_t a, std::size_t b) { return EdgeKey{std::min(a, b), std::max(a, b)}; };
    std::map<EdgeKey, std::vector<std::size_t>> edge_facets;
    for (std::size_t f = 0; f < F.size(); ++f) {
        const auto& v = g.facets[f].vertices;
        for (std::size_t i = 0; i < v.size(); ++i) edge_facets[key(v[i], v[(i + 1) % v.size()])].push_back(f);
    }
    const auto regions_of = [&](std::size_t f) {
        const std::int32_t a = g.volumes[static_cast<std::size_t>(F[f].volumes[0])].region;
        const std::int32_t b = g.volumes[static_cast<std::size_t>(F[f].volumes[1])].region;
        return std::pair{std::min(a, b), std::max(a, b)};
    };
    const auto coplanar = [&](std::size_t f, std::size_t h) {
        const auto& pl = F[f].plane;
        return std::ranges::all_of(g.facets[h].vertices, [&](std::size_t i) { return predicates::orient3d(P[pl[0]], P[pl[1]], P[pl[2]], P[i]) == 0; });
    };
    const auto h_of = [&](std::int32_t r) {
        return static_cast<std::size_t>(r) < in.region_size.size() ? in.region_size[static_cast<std::size_t>(r)] : in.region_size.back();
    };

    // Pieces: interface facets of one pair joined across shared edges when coplanar.
    std::vector<std::size_t> parent(F.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto find = [&](std::size_t x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    for (const auto& [e, fs] : edge_facets) {
        for (std::size_t i = 0; i < fs.size(); ++i) {
            for (std::size_t j = i + 1; j < fs.size(); ++j) {
                const std::size_t f = fs[i];
                const std::size_t h = fs[j];
                if (F[f].kind != FacetKind::interface || F[h].kind != FacetKind::interface) continue;
                if (regions_of(f) != regions_of(h) || !coplanar(f, h)) continue;
                parent[find(f)] = find(h);
            }
        }
    }
    std::vector<int> piece_of(F.size(), -1);
    {
        std::map<std::size_t, int> root_piece;
        for (const std::size_t f : interface_facets) {
            const auto [it, fresh] = root_piece.emplace(find(f), static_cast<int>(out.pieces.size()));
            if (fresh) out.pieces.emplace_back();
            out.pieces[static_cast<std::size_t>(it->second)].facets.push_back(f);
            piece_of[f] = it->second;
        }
    }

    // Each piece's boundary loop, corners and frame.
    for (InterfacePieceLayout& pc : out.pieces) {
        const std::size_t f0 = pc.facets[0];
        std::tie(pc.region_a, pc.region_b) = regions_of(f0);
        std::map<EdgeKey, int> count;
        for (const std::size_t f : pc.facets) {
            const auto& v = g.facets[f].vertices;
            for (std::size_t i = 0; i < v.size(); ++i) ++count[key(v[i], v[(i + 1) % v.size()])];
        }
        std::map<std::size_t, std::vector<std::size_t>> next;
        std::size_t edges = 0;
        for (const auto& [e, c] : count) {
            if (c != 1) continue;
            next[e.first].push_back(e.second);
            next[e.second].push_back(e.first);
            ++edges;
        }
        for (const auto& [v, list] : next) {
            if (list.size() != 2) return std::unexpected(invalid("an interface piece is not a simple polygon (a pinched vertex; unsupported)", f0));
        }
        std::size_t prev = next.begin()->first;
        std::size_t at = next.begin()->second[0];
        pc.loop.push_back(prev);
        while (at != pc.loop[0]) {
            pc.loop.push_back(at);
            const auto& nb = next[at];
            const std::size_t to = nb[0] == prev ? nb[1] : nb[0];
            prev = at;
            at = to;
            if (pc.loop.size() > edges) break;
        }
        if (pc.loop.size() != edges) return std::unexpected(invalid("an interface piece has holes (unsupported)", f0));
        const std::size_t m = pc.loop.size();
        std::vector<std::size_t> corners;
        for (std::size_t i = 0; i < m; ++i) {
            if (!collinear(P[pc.loop[(i + m - 1) % m]], P[pc.loop[i]], P[pc.loop[(i + 1) % m]])) corners.push_back(pc.loop[i]);
        }
        for (int a = 0; a < 3 && pc.axis < 0; ++a) {
            const auto aa = static_cast<std::size_t>(a);
            if (std::ranges::all_of(pc.loop, [&](std::size_t v) { return P[v][aa] == P[pc.loop[0]][aa]; })) pc.axis = a;
        }
        for (const std::size_t v : corners) pc.polygon.push_back(P[v]);
        if (pc.axis >= 0) {
            pc.normal = {0.0, 0.0, 0.0};
            pc.normal[static_cast<std::size_t>(pc.axis)] = 1.0;
        } else {
            const auto& pl = F[f0].plane;
            pc.normal = unit(cross(sub(P[pl[1]], P[pl[0]]), sub(P[pl[2]], P[pl[0]])));
            pc.origin = pc.polygon[0];
            pc.u = unit(sub(pc.polygon[1], pc.polygon[0]));
            pc.w = cross(pc.normal, pc.u);
        }
        for (const Point& x : pc.polygon) pc.local.push_back(pc.to_local(x));
        if (polygon_turn(pc.local) < 0) {
            std::ranges::reverse(pc.polygon);
            std::ranges::reverse(pc.local);
            std::ranges::reverse(pc.loop);
        }
        pc.spacing = std::min({h_of(pc.region_a), h_of(pc.region_b), 4.0 * s});
        out.band = std::max(out.band, std::sqrt(pc.spacing * pc.spacing / 3.0 + s * s / 4.0));
    }

    // Feature lines: the sectors about every piece edge, and the traces.
    std::set<EdgeKey> checked;
    std::vector<std::vector<int>> loop_trace_face(out.pieces.size());  // per loop edge: its face if a trace, else -1
    for (std::size_t k = 0; k < out.pieces.size(); ++k) {
        const InterfacePieceLayout& pc = out.pieces[k];
        const std::size_t f0 = pc.facets[0];
        const std::size_t m = pc.loop.size();
        loop_trace_face[k].assign(m, -1);
        for (std::size_t i = 0; i < m; ++i) {
            const std::size_t a = pc.loop[i];
            const std::size_t b = pc.loop[(i + 1) % m];
            const auto& fs = edge_facets.at(key(a, b));
            std::vector<std::size_t> bnd;
            for (const std::size_t f : fs) {
                if (F[f].kind == FacetKind::boundary) bnd.push_back(f);
            }
            if (!bnd.empty()) {
                if (bnd.size() != 2 || in.face_of[bnd[0]] != in.face_of[bnd[1]]) {
                    return std::unexpected(invalid("an interface ends on a domain edge (the boundary must be planar across a trace; unsupported)", f0));
                }
                for (const std::size_t f : fs) {
                    if (F[f].kind == FacetKind::interface && piece_of[f] != static_cast<int>(k)) {
                        return std::unexpected(invalid("two or more interfaces end on the boundary along one line (unsupported)", f0));
                    }
                }
                loop_trace_face[k][i] = in.face_of[bnd[0]];
            }
            if (!checked.insert(key(a, b)).second) continue;
            // Sectors between consecutive interface and boundary facets about the edge; the
            // exterior one (between the two boundary facets) is skipped.
            std::vector<std::pair<Point, std::size_t>> h;
            for (const std::size_t f : fs) {
                if (F[f].kind == FacetKind::internal) continue;
                for (const auto& t : F[f].triangles) {
                    if (std::ranges::find(t, a) == t.end() || std::ranges::find(t, b) == t.end()) continue;
                    for (const std::size_t v : t) {
                        if (v != a && v != b) h.emplace_back(P[v], f);
                    }
                    break;
                }
            }
            sort_about_line(P[a], P[b], h);
            for (std::size_t j = 0; j < h.size(); ++j) {
                const auto& [p, fp] = h[j];
                const auto& [q, fq] = h[(j + 1) % h.size()];
                if (fp == fq) continue;
                if (F[fp].kind == FacetKind::boundary && F[fq].kind == FacetKind::boundary) continue;
                if (predicates::orient3d(P[a], P[b], p, q) > 0 && predicates::dihedral_below_60(P[a], P[b], p, q)) {
                    return std::unexpected(invalid("an interface sector below 60 degrees (refused until justified, DECISIONS.md U12, FD7)", f0));
                }
            }
        }
    }

    // Traces: maximal straight runs of trace edges in one face.
    std::vector<std::vector<int>> vertex_chains(P.size());
    for (std::size_t c = 0; c < in.chains.size(); ++c) {
        for (const std::size_t v : in.chains[c].chain) vertex_chains[v].push_back(static_cast<int>(c));
    }
    for (std::size_t k = 0; k < out.pieces.size(); ++k) {
        const InterfacePieceLayout& pc = out.pieces[k];
        const std::size_t m = pc.loop.size();
        const auto& tf = loop_trace_face[k];
        const auto is_corner = [&](std::size_t i) {
            return !collinear(P[pc.loop[(i + m - 1) % m]], P[pc.loop[i]], P[pc.loop[(i + 1) % m]]);
        };
        // A break before loop edge i: a corner at its start, or a change of trace face.
        const auto breaks = [&](std::size_t i) { return is_corner(i) || tf[i] != tf[(i + m - 1) % m]; };
        std::size_t i0 = 0;
        while (i0 < m && !breaks(i0)) ++i0;
        if (i0 == m) continue;  // impossible for a polygon; defensive
        for (std::size_t n = 0; n < m;) {
            const std::size_t i = (i0 + n) % m;
            std::size_t len = 1;
            while (n + len < m && !breaks((i + len) % m)) ++len;
            if (tf[i] >= 0) {
                InterfaceTrace tr;
                tr.piece = static_cast<int>(k);
                tr.face = tf[i];
                const Point& a = P[pc.loop[i]];
                const Point& b = P[pc.loop[(i + len) % m]];
                tr.start = a;
                tr.length = dist(a, b);
                tr.dir = scale(1.0 / tr.length, sub(b, a));
                const auto& pl = in.face_plane[static_cast<std::size_t>(tr.face)];
                tr.normal = unit(cross(cross(sub(pl[1], pl[0]), sub(pl[2], pl[0])), tr.dir));
                if (!(tr.length >= 6.0 * s)) return std::unexpected(invalid("an interface trace on the boundary is shorter than 6 s", pc.facets[0]));
                out.traces.push_back(tr);
            }
            n += len;
        }
    }

    // Feature points on the boundary: crossings of domain edges, and points in a face where traces meet.
    std::map<std::size_t, std::vector<std::size_t>> loops_at;  // complex vertex -> pieces whose loop holds it
    for (std::size_t k = 0; k < out.pieces.size(); ++k) {
        for (const std::size_t v : out.pieces[k].loop) loops_at[v].push_back(k);
    }
    for (const auto& [x, pieces] : loops_at) {
        const std::size_t f0 = out.pieces[pieces[0]].facets[0];
        // The trace edges at x: (piece, neighbour vertex, face).
        struct Arm {
            std::size_t piece;
            std::size_t to;
            int face;
        };
        std::vector<Arm> arms;
        for (const std::size_t k : pieces) {
            const auto& loop = out.pieces[k].loop;
            const std::size_t m = loop.size();
            const auto i = static_cast<std::size_t>(std::ranges::find(loop, x) - loop.begin());
            const int before = loop_trace_face[k][(i + m - 1) % m];
            const int after = loop_trace_face[k][i];
            if (before >= 0) arms.push_back({k, loop[(i + m - 1) % m], before});
            if (after >= 0) arms.push_back({k, loop[(i + 1) % m], after});
        }
        const auto& chains = vertex_chains[x];
        if (chains.empty()) {
            if (arms.size() < 2) continue;
            // In one face's interior: sectors between consecutive traces.
            const int fi = arms[0].face;
            const auto& pl = in.face_plane[static_cast<std::size_t>(fi)];
            const Point nf = cross(sub(pl[1], pl[0]), sub(pl[2], pl[0]));
            const Point q = add(P[x], unit(nf));
            std::vector<Point> rays;
            for (const Arm& r : arms) rays.push_back(P[r.to]);
            sort_about_point(P[x], q, rays);
            for (std::size_t j = 0; j < rays.size(); ++j) {
                const Point& a = rays[j];
                const Point& b = rays[(j + 1) % rays.size()];
                if (predicates::orient3d(P[x], a, b, q) > 0 && predicates::angle_below_60_3d(P[x], a, b)) {
                    return std::unexpected(invalid("interface traces meet in a boundary face at under 60 degrees", f0));
                }
            }
            continue;
        }
        // FD8 (b), DECISIONS.md SQ5: a constrained patch edge meeting a trace inside its face,
        // crossing it or ending on it; the rules of a domain-edge crossing.
        const auto patch_chain = std::ranges::find_if(chains, [&](int c) { return in.chains[static_cast<std::size_t>(c)].face >= 0; });
        if (patch_chain != chains.end()) {
            const BoundaryChain& pch = in.chains[static_cast<std::size_t>(*patch_chain)];
            if (chains.size() != 1) {
                return std::unexpected(invalid("a patch edge meets an interface trace at a domain edge or another patch edge (unsupported)", f0));
            }
            if (pieces.size() != 1 || arms.size() != 2 || arms[0].face != pch.face || arms[1].face != pch.face) {
                return std::unexpected(invalid("a patch edge meets an interface other than at a single trace in its face (unsupported)", f0));
            }
            const auto pos = static_cast<std::size_t>(std::ranges::find(pch.chain, x) - pch.chain.begin());
            std::vector<std::size_t> along;  // the chain's vertices beside x
            if (pos > 0) along.push_back(pch.chain[pos - 1]);
            if (pos + 1 < pch.chain.size()) along.push_back(pch.chain[pos + 1]);
            for (const Arm& r : arms) {
                for (const std::size_t e : along) {
                    if (predicates::angle_below_60_3d(P[x], P[r.to], P[e])) {
                        return std::unexpected(invalid("a patch edge meets an interface trace at under 60 degrees", f0));
                    }
                }
            }
            for (const std::size_t e : along) {
                if (!(dist(P[x], P[e]) >= 4.0 * s)) {
                    return std::unexpected(invalid("a patch-edge vertex within 4 s of where a patch edge meets an interface trace", f0));
                }
            }
            out.crossings.push_back({*patch_chain, pos});
            continue;
        }
        const BoundaryChain& ch = in.chains[static_cast<std::size_t>(chains[0])];
        const auto pos = static_cast<std::size_t>(std::ranges::find(ch.chain, x) - ch.chain.begin());
        if (chains.size() != 1 || pos == 0 || pos + 1 == ch.chain.size()) {
            return std::unexpected(invalid("an interface meets the boundary at a domain corner (unsupported)", f0));
        }
        if (pieces.size() != 1 || arms.size() != 2 || arms[0].face == arms[1].face ||
            (arms[0].face != ch.faces[0] && arms[0].face != ch.faces[1]) || (arms[1].face != ch.faces[0] && arms[1].face != ch.faces[1])) {
            return std::unexpected(invalid("an interface feature line meets a domain edge (only a single interface crossing it; unsupported)", f0));
        }
        const Point& e1 = P[ch.chain[pos - 1]];
        const Point& e2 = P[ch.chain[pos + 1]];
        for (const Arm& r : arms) {
            for (const Point& e : {e1, e2}) {
                if (predicates::angle_below_60_3d(P[x], P[r.to], e)) {
                    return std::unexpected(invalid("an interface trace meets a domain edge at under 60 degrees", f0));
                }
            }
        }
        if (!(dist(P[x], e1) >= 4.0 * s && dist(P[x], e2) >= 4.0 * s)) {
            return std::unexpected(invalid("a domain-edge vertex within 4 s of an interface crossing", f0));
        }
        out.crossings.push_back({chains[0], pos});
    }

    // Lengths and separations (double; FD7 and the thin-layer reading, DECISIONS.md slice 3).
    struct Side {
        std::size_t piece;
        Point a, b;
        std::vector<std::size_t> vertices;  // the loop vertices along it
    };
    std::vector<Side> sides;
    for (std::size_t k = 0; k < out.pieces.size(); ++k) {
        const auto& loop = out.pieces[k].loop;
        const std::size_t m = loop.size();
        std::size_t i0 = 0;
        while (collinear(P[loop[(i0 + m - 1) % m]], P[loop[i0]], P[loop[(i0 + 1) % m]]) == true) ++i0;
        Side cur{k, P[loop[i0]], {}, {loop[i0]}};
        for (std::size_t n = 1; n <= m; ++n) {
            const std::size_t i = (i0 + n) % m;
            cur.vertices.push_back(loop[i]);
            if (n == m || !collinear(P[loop[(i + m - 1) % m]], P[loop[i]], P[loop[(i + 1) % m]])) {
                cur.b = P[loop[i]];
                if (!(dist(cur.a, cur.b) >= 6.0 * s)) return std::unexpected(invalid("an interface feature line is shorter than 6 s", out.pieces[k].facets[0]));
                sides.push_back(cur);
                cur = Side{k, P[loop[i]], {}, {loop[i]}};
            }
        }
    }
    const auto share = [](const std::vector<std::size_t>& x, const std::vector<std::size_t>& y) {
        return std::ranges::any_of(x, [&](std::size_t v) { return std::ranges::find(y, v) != y.end(); });
    };
    for (std::size_t i = 0; i < sides.size(); ++i) {
        for (std::size_t j = i + 1; j < sides.size(); ++j) {
            if (share(sides[i].vertices, sides[j].vertices)) continue;
            if (!(segment_segment(sides[i].a, sides[i].b, sides[j].a, sides[j].b) >= 4.0 * s)) {
                return std::unexpected(invalid("two unrelated interface features within 4 s", out.pieces[sides[i].piece].facets[0]));
            }
        }
        for (const BoundaryChain& ch : in.chains) {
            if (share(sides[i].vertices, ch.chain)) continue;
            if (!(segment_segment(sides[i].a, sides[i].b, P[ch.chain.front()], P[ch.chain.back()]) >= 4.0 * s)) {
                // A patch edge (slice 4) as a domain edge: its subsegments must stay clear of the pairs.
                return std::unexpected(invalid(ch.face >= 0 ? "a patch edge within 4 s of an interface feature it does not meet"
                                                            : "an interface feature within 4 s of a domain edge",
                                               out.pieces[sides[i].piece].facets[0]));
            }
        }
    }
    // Pieces against pieces and boundary faces that they do not touch, and parallel ones.
    std::vector<std::set<std::size_t>> piece_vertices(out.pieces.size());
    for (std::size_t k = 0; k < out.pieces.size(); ++k) {
        for (const std::size_t f : out.pieces[k].facets) piece_vertices[k].insert(g.facets[f].vertices.begin(), g.facets[f].vertices.end());
    }
    const auto parallel = [](const Point& n1, const Point& n2) { return norm3(cross(n1, n2)) <= 1e-12; };
    for (std::size_t i = 0; i < out.pieces.size(); ++i) {
        const InterfacePieceLayout& A = out.pieces[i];
        for (std::size_t j = i + 1; j < out.pieces.size(); ++j) {
            const InterfacePieceLayout& B = out.pieces[j];
            const bool touching = std::ranges::any_of(piece_vertices[i], [&](std::size_t v) { return piece_vertices[j].contains(v); });
            const bool par = parallel(A.normal, B.normal);
            const bool coplanar_pieces = std::ranges::all_of(B.polygon, [&](const Point& x) {
                return predicates::orient3d(A.polygon[0], A.polygon[1], A.polygon[2], x) == 0;
            });
            if (coplanar_pieces) continue;  // neighbours in one plane (a junction at 180 degrees)
            if (touching && !par) continue;
            double d = std::numeric_limits<double>::infinity();
            for (const Point& x : A.polygon) d = std::min(d, B.distance(x));
            for (const Point& x : B.polygon) d = std::min(d, A.distance(x));
            for (std::size_t a = 0; a < A.polygon.size(); ++a) {
                for (std::size_t b = 0; b < B.polygon.size(); ++b) {
                    d = std::min(d, segment_segment(A.polygon[a], A.polygon[(a + 1) % A.polygon.size()], B.polygon[b],
                                                    B.polygon[(b + 1) % B.polygon.size()]));
                }
            }
            if (!(d >= 4.0 * s)) return std::unexpected(invalid("two interfaces within 4 s (a layer thinner than 4 s, or unrelated interfaces)", B.facets[0]));
        }
        for (std::size_t fi = 0; fi < in.face_plane.size(); ++fi) {
            bool touching = false;
            for (std::size_t f = 0; f < F.size() && !touching; ++f) {
                if (in.face_of[f] != static_cast<int>(fi)) continue;
                touching = std::ranges::any_of(g.facets[f].vertices, [&](std::size_t v) { return piece_vertices[i].contains(v); });
            }
            if (touching) continue;
            double d = std::numeric_limits<double>::infinity();
            for (const auto& t : in.face_triangles[fi]) {
                for (const Point& x : A.polygon) d = std::min(d, point_triangle(x, t[0], t[1], t[2]));
                for (const Point& x : t) d = std::min(d, A.distance(x));
                for (std::size_t a = 0; a < A.polygon.size(); ++a) {
                    for (std::size_t e = 0; e < 3; ++e) {
                        d = std::min(d, segment_segment(A.polygon[a], A.polygon[(a + 1) % A.polygon.size()], t[e], t[(e + 1) % 3]));
                    }
                }
            }
            if (!(d >= 4.0 * s)) return std::unexpected(invalid("an interface within 4 s of a boundary face it does not meet", A.facets[0]));
        }
    }

    // The inset polygons (1.5 s): each offset edge keeps its direction and the result is simple.
    for (InterfacePieceLayout& pc : out.pieces) {
        const std::vector<Point>& L = pc.local;
        const std::size_t nv = L.size();
        for (std::size_t i = 0; i < nv; ++i) {
            const auto offset_line = [&](std::size_t a, std::size_t b) {
                const Point dd = sub(L[b], L[a]);
                const double len = std::hypot(dd[0], dd[1]);
                const Point inward{-dd[1] / len, dd[0] / len, 0.0};
                return std::pair{add(L[a], scale(1.5 * s, inward)), dd};
            };
            const auto [p1, d1] = offset_line((i + nv - 1) % nv, i);
            const auto [p2, d2] = offset_line(i, (i + 1) % nv);
            const double den = d1[0] * d2[1] - d1[1] * d2[0];
            const double t = ((p2[0] - p1[0]) * d2[1] - (p2[1] - p1[1]) * d2[0]) / den;
            pc.inset.push_back({p1[0] + t * d1[0], p1[1] + t * d1[1], 0.0});
        }
        bool ok = simple_2d(pc.inset) && polygon_turn(pc.inset) > 0;
        for (std::size_t i = 0; i < nv && ok; ++i) {
            const Point d0 = sub(L[(i + 1) % nv], L[i]);
            const Point d1 = sub(pc.inset[(i + 1) % nv], pc.inset[i]);
            ok = d0[0] * d1[0] + d0[1] * d1[1] > 0.0;
        }
        if (!ok) return std::unexpected(invalid("an interface piece is too narrow for its 1.5 s inset", pc.facets[0]));
    }
    return out;
}

}  // namespace NiTCAD::mesh::detail
