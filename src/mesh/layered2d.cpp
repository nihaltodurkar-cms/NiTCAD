// Layered 2D geometry resolved into a tiling (Unit 26k, GD1-GD2; generate.hpp).
//
// The overlay: every shape's x and y breakpoints (cut to the domain) make a grid of cells; each cell takes the
// region of the last shape covering it (an even-odd scan of the shape's vertical edges per row, exact: only
// comparisons of input coordinates), else the background; with `fill` shapes (Unit 26m MT4), the first in
// paint_order's order. Same-region cells sharing an edge form parts; each part's
// boundary is traced with the part on its left (counter-clockwise). A part whose boundary touches itself at a grid
// vertex (two outgoing edges there) or has more than one loop (a hole) is cut along the vertical breakpoint line
// through that vertex or the hole's leftmost line, and its pieces are traced again: each cut strictly divides the
// part, so this ends. Last, every polygon takes on its edges the corners of all polygons (T-junctions).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "arrangement2d.hpp"
#include "polygon2d.hpp"
#include "predicates.hpp"

namespace NiTCAD::mesh {

namespace {

std::unexpected<base::Error> invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return std::unexpected(base::Error{base::ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}});
}

// The order a point asks the shapes in, the first holding it taking it (Unit 26m MT4): the replacing shapes top
// first, then the filling ones bottom first. With no `fill` shape, top first as before.
std::vector<std::size_t> paint_order(const LayeredGeometry2D& g) {
    std::vector<std::size_t> order;
    order.reserve(g.shapes.size());
    for (std::size_t s = g.shapes.size(); s-- > 0;) {
        if (!g.shapes[s].fill) order.push_back(s);
    }
    for (std::size_t s = 0; s < g.shapes.size(); ++s) {
        if (g.shapes[s].fill) order.push_back(s);
    }
    return order;
}

// The breakpoints along one axis: the domain's ends and every vertex coordinate inside them.
std::vector<double> breakpoints(const LayeredGeometry2D& g, int axis) {
    std::vector<double> v{g.lo[static_cast<std::size_t>(axis)], g.hi[static_cast<std::size_t>(axis)]};
    for (const LayeredShape& s : g.shapes) {
        for (const Point& p : s.vertices) {
            const double c = p[static_cast<std::size_t>(axis)];
            if (c > v[0] && c < v[1]) v.push_back(c);
        }
    }
    std::ranges::sort(v);
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

// The index of the breakpoint at or nearest inside the domain for a coordinate (cut to the domain).
std::size_t index_of(const std::vector<double>& b, double c) {
    if (c <= b.front()) return 0;
    if (c >= b.back()) return b.size() - 1;
    return static_cast<std::size_t>(std::ranges::lower_bound(b, c) - b.begin());  // c is a breakpoint
}


// --- slanted edges (Unit 26k slice 4, GD10): the layers as a planar arrangement -------------------------------------
//
// Every shape's edges and the domain's sides, split where they meet: a touch or a collinear overlap at the input's own
// points (exact, orient2d), a proper crossing at a point computed once, rounded, and shared by every segment through
// it (computed crossings within a tolerance of one another, or of an input point, are one point). Each split edge
// knows, side by side, the region there: the topmost shape it is inside (its own shape by the edge's direction, others
// by an exact point-in-polygon test of its midpoint), else the background, or nothing outside the domain. A region's
// boundary is its edges with another region (or nothing) across; its loops are traced with the region on the left,
// turning at each vertex to the first boundary edge clockwise, and split where they pass a vertex twice. A hole (a
// clockwise loop) is cut by two bridges between existing corners (inside the part, crossing nothing; on the domain's
// boundary off every patch's inside first, then the shorter), and the region traced again. A polygon keeps the
// vertices where three or more boundary edges meet or where its two edges come from lines that differ: the points
// along a straight edge that only split it are dropped, the rounded ones too.

int orient(const Point& a, const Point& b, const Point& c) { return predicates::orient2d(a, b, c); }
bool same_pt(const Point& a, const Point& b) { return a[0] == b[0] && a[1] == b[1]; }
bool lower(const Point& a, const Point& b) { return a[1] < b[1] || (a[1] == b[1] && a[0] < b[0]); }  // lowest, leftmost

// The side of the line ab that `p` is on (exact, orient2d), 0 when on it or within `tol` of it.
int side(const Point& a, const Point& b, const Point& p, double tol) {
    const int o = orient(a, b, p);
    if (o == 0) return 0;
    const double bx = b[0] - a[0], by = b[1] - a[1], px = p[0] - a[0], py = p[1] - a[1];
    return std::abs(bx * py - by * px) <= tol * std::hypot(bx, by) ? 0 : o;
}
// `p`, on the line ab or within `tol` of it (side 0), touches the segment: exactly on it, or projecting inside it.
bool touches(const Point& a, const Point& b, const Point& p, double tol) {
    if (orient(a, b, p) == 0) return detail::on_segment_2d(a, b, p);
    if (same_pt(a, p) || same_pt(b, p)) return true;
    const double bx = b[0] - a[0], by = b[1] - a[1], px = p[0] - a[0], py = p[1] - a[1];
    const double u = (px * bx + py * by) / (bx * bx + by * by);
    return u > 0.0 && u < 1.0 && side(a, b, p, tol) == 0;
}

std::expected<Geometry2D, base::Error> general_partition(const LayeredGeometry2D& g, const detail::ArrangementLabels* how = nullptr) {
    struct Orig {
        Point a, b;
        std::int32_t shape;  // -1: the domain's side
    };
    const Point c0{g.lo[0], g.lo[1], 0}, c1{g.hi[0], g.lo[1], 0}, c2{g.hi[0], g.hi[1], 0}, c3{g.lo[0], g.hi[1], 0};
    std::vector<Orig> orig{{c0, c1, -1}, {c1, c2, -1}, {c2, c3, -1}, {c3, c0, -1}};
    std::vector<std::vector<Point>> ccw(g.shapes.size());
    std::vector<std::array<double, 4>> box(g.shapes.size());
    for (std::size_t s = 0; s < g.shapes.size(); ++s) {
        std::vector<Point> v;
        for (const Point& p : g.shapes[s].vertices) v.push_back({p[0], p[1], 0});
        if (detail::polygon_turn(v) < 0) std::ranges::reverse(v);
        box[s] = {v[0][0], v[0][0], v[0][1], v[0][1]};
        for (std::size_t k = 0; k < v.size(); ++k) {
            orig.push_back({v[k], v[(k + 1) % v.size()], static_cast<std::int32_t>(s)});
            box[s] = {std::min(box[s][0], v[k][0]), std::max(box[s][1], v[k][0]), std::min(box[s][2], v[k][1]),
                      std::max(box[s][3], v[k][1])};
        }
        ccw[s] = std::move(v);
    }
    const double tol = 1e-12 * std::max(g.hi[0] - g.lo[0], g.hi[1] - g.lo[1]);

    // The split points of every segment.
    const std::size_t n = orig.size();
    std::vector<std::vector<Point>> splits(n);
    for (std::size_t i = 0; i < n; ++i) splits[i] = {orig[i].a, orig[i].b};
    // The points so far by x: every input vertex, then each crossing made. A crossing within the tolerance of one
    // (any shape's vertex, not only its two segments' ends) is that point.
    std::multimap<double, Point> made;
    for (const Orig& o : orig) made.emplace(o.a[0], o.a);
    if (how != nullptr) {
        for (const Point& p : how->pins) made.emplace(p[0], Point{p[0], p[1], 0});
    }
    const auto share = [&](Point p, const Orig& s, const Orig& t) {
        const auto near = [&](const Point& q) { return std::abs(q[0] - p[0]) <= tol && std::abs(q[1] - p[1]) <= tol; };
        for (const Point* q : {&s.a, &s.b, &t.a, &t.b}) {
            if (near(*q)) return *q;
        }
        for (auto it = made.lower_bound(p[0] - tol); it != made.end() && it->first <= p[0] + tol; ++it) {
            if (near(it->second)) return it->second;
        }
        made.emplace(p[0], p);
        return p;
    };
    // The constant coordinates of the horizontal and vertical segments: a computed crossing within the tolerance of
    // one takes it exactly, so two crossings merged by `share` both lie on the line (one on it, one rounded off it).
    std::array<std::vector<double>, 2> axes;
    for (const Orig& o : orig) {
        for (std::size_t k = 0; k < 2; ++k) {
            if (o.a[k] == o.b[k]) axes[k].push_back(o.a[k]);
        }
    }
    for (auto& a : axes) {
        std::ranges::sort(a);
        a.erase(std::unique(a.begin(), a.end()), a.end());
    }
    const auto on_axes = [&](Point p) {
        for (std::size_t k = 0; k < 2; ++k) {
            const auto it = std::ranges::lower_bound(axes[k], p[k]);
            if (it != axes[k].end() && *it - p[k] <= tol) p[k] = *it;
            else if (it != axes[k].begin() && p[k] - *(it - 1) <= tol) p[k] = *(it - 1);
        }
        return p;
    };
    const auto crossing = [&on_axes](const Orig& s, const Orig& t) {
        // Along an axis-parallel segment the crossing keeps its constant coordinate exactly.
        const auto axis_hit = [](const Orig& h, const Orig& o, std::size_t k) {
            const std::size_t m = 1 - k;
            const double c = h.a[k], u = (c - o.a[k]) / (o.b[k] - o.a[k]);
            Point r{0, 0, 0};
            r[k] = c;
            r[m] = o.a[m] + u * (o.b[m] - o.a[m]);
            return r;
        };
        Point p{0, 0, 0};
        if (s.a[1] == s.b[1]) p = axis_hit(s, t, 1);
        else if (s.a[0] == s.b[0]) p = axis_hit(s, t, 0);
        else if (t.a[1] == t.b[1]) p = axis_hit(t, s, 1);
        else if (t.a[0] == t.b[0]) p = axis_hit(t, s, 0);
        else {
            const double ax = s.b[0] - s.a[0], ay = s.b[1] - s.a[1], bx = t.b[0] - t.a[0], by = t.b[1] - t.a[1];
            const double u = ((t.a[0] - s.a[0]) * by - (t.a[1] - s.a[1]) * bx) / (ax * by - ay * bx);
            p = {s.a[0] + u * ax, s.a[1] + u * ay, 0};
        }
        for (std::size_t k = 0; k < 2; ++k) {
            const double lo = std::max(std::min(s.a[k], s.b[k]), std::min(t.a[k], t.b[k]));
            const double hi = std::min(std::max(s.a[k], s.b[k]), std::max(t.a[k], t.b[k]));
            p[k] = std::clamp(p[k], lo, hi);
        }
        return on_axes(p);
    };
    std::vector<std::size_t> by_x(n);
    for (std::size_t i = 0; i < n; ++i) by_x[i] = i;
    const auto xmin = [&](std::size_t i) { return std::min(orig[i].a[0], orig[i].b[0]); };
    const auto xmax = [&](std::size_t i) { return std::max(orig[i].a[0], orig[i].b[0]); };
    std::ranges::sort(by_x, [&](std::size_t a, std::size_t b) { return xmin(a) < xmin(b) || (xmin(a) == xmin(b) && a < b); });
    for (std::size_t p = 0; p < n; ++p) {
        const std::size_t i = by_x[p];
        const Orig& s = orig[i];
        for (std::size_t q = p + 1; q < n && xmin(by_x[q]) <= xmax(i) + tol; ++q) {
            const std::size_t j = by_x[q];
            const Orig& t = orig[j];
            if (std::max(s.a[1], s.b[1]) + tol < std::min(t.a[1], t.b[1]) || std::max(t.a[1], t.b[1]) + tol < std::min(s.a[1], s.b[1])) continue;
            // Sides within the tolerance of the line count as on it, touches within it of the segment's inside as
            // touches: a vertex, or a whole edge, meant to lie on a slanted edge but rounded off it (decimal
            // coordinates) splits it there instead of crossing it at a computed point, so no sliver is left between.
            const int o1 = side(s.a, s.b, t.a, tol), o2 = side(s.a, s.b, t.b, tol);
            const int o3 = side(t.a, t.b, s.a, tol), o4 = side(t.a, t.b, s.b, tol);
            if (o1 * o2 < 0 && o3 * o4 < 0) {
                const Point x = share(crossing(s, t), s, t);
                splits[i].push_back(x);
                splits[j].push_back(x);
                continue;
            }
            if (o1 == 0 && touches(s.a, s.b, t.a, tol)) splits[i].push_back(t.a);
            if (o2 == 0 && touches(s.a, s.b, t.b, tol)) splits[i].push_back(t.b);
            if (o3 == 0 && touches(t.a, t.b, s.a, tol)) splits[j].push_back(s.a);
            if (o4 == 0 && touches(t.a, t.b, s.b, tol)) splits[j].push_back(s.b);
        }
    }

    // Pins (labelled use): every segment through one, exactly, is split there.
    std::vector<Point> pins;
    if (how != nullptr) {
        for (const Point& p0 : how->pins) {
            const Point p{p0[0], p0[1], 0};
            pins.push_back(p);
            for (std::size_t i = 0; i < n; ++i) {
                if (orient(orig[i].a, orig[i].b, p) == 0 && touches(orig[i].a, orig[i].b, p, 0.0)) splits[i].push_back(p);
            }
        }
    }

    // The vertices, and the split edges (one per pair of vertices, every segment along it its source).
    std::vector<Point> pts;
    for (const auto& sp : splits) pts.insert(pts.end(), sp.begin(), sp.end());
    const auto by_xy = [](const Point& a, const Point& b) { return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]); };
    std::ranges::sort(pts, by_xy);
    pts.erase(std::unique(pts.begin(), pts.end(), same_pt), pts.end());
    const auto id = [&](const Point& p) { return static_cast<std::size_t>(std::ranges::lower_bound(pts, p, by_xy) - pts.begin()); };
    struct Edge {
        std::size_t u, v;
        std::vector<std::size_t> src;  // the segments along it (a bridge: none)
        bool cut = false;
    };
    std::vector<Edge> edges;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge_of;
    for (std::size_t i = 0; i < n; ++i) {
        auto& sp = splits[i];
        const Point a = orig[i].a;
        const double dx = orig[i].b[0] - a[0], dy = orig[i].b[1] - a[1];
        std::ranges::sort(sp, [&](const Point& p, const Point& q) {
            return (p[0] - a[0]) * dx + (p[1] - a[1]) * dy < (q[0] - a[0]) * dx + (q[1] - a[1]) * dy;
        });
        sp.erase(std::unique(sp.begin(), sp.end(), same_pt), sp.end());
        for (std::size_t k = 0; k + 1 < sp.size(); ++k) {
            const std::size_t u = id(sp[k]), v = id(sp[k + 1]);
            if (u == v) continue;
            const auto key = std::minmax(u, v);
            const auto [it, fresh] = edge_of.try_emplace({key.first, key.second}, edges.size());
            if (fresh) edges.push_back({u, v, {i}, false});
            else edges[it->second].src.push_back(i);
        }
    }

    // Half-edges 2e (u to v) and 2e + 1, each vertex's outgoing ones counter-clockwise by direction (exact).
    std::vector<std::vector<std::size_t>> out(pts.size());
    std::vector<std::size_t> pos;
    const auto from = [&](std::size_t h) { return h % 2 == 0 ? edges[h / 2].u : edges[h / 2].v; };
    const auto to = [&](std::size_t h) { return h % 2 == 0 ? edges[h / 2].v : edges[h / 2].u; };
    const auto sort_out = [&](std::size_t v) {
        const Point& o = pts[v];
        const auto half = [&](std::size_t h) {
            const Point& d = pts[to(h)];
            return d[1] > o[1] || (d[1] == o[1] && d[0] > o[0]) ? 0 : 1;
        };
        std::ranges::sort(out[v], [&](std::size_t a, std::size_t b) {
            const int ha = half(a), hb = half(b);
            return ha != hb ? ha < hb : orient(o, pts[to(a)], pts[to(b)]) > 0;
        });
        for (std::size_t k = 0; k < out[v].size(); ++k) pos[out[v][k]] = k;
    };
    pos.assign(2 * edges.size(), 0);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        out[edges[e].u].push_back(2 * e);
        out[edges[e].v].push_back(2 * e + 1);
    }
    for (std::size_t v = 0; v < pts.size(); ++v) sort_out(v);

    // The region on each half-edge's left: -2 outside the domain.
    constexpr std::int32_t outside = -2;
    const std::vector<std::size_t> order = paint_order(g);
    std::vector<std::int32_t> left(2 * edges.size(), outside);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const Point& u = pts[edges[e].u];
        const Point& v = pts[edges[e].v];
        const Point m{0.5 * (u[0] + v[0]), 0.5 * (u[1] + v[1]), 0};
        const auto along = [&](std::size_t i) {  // segment i runs from u to v
            return (v[0] - u[0]) * (orig[i].b[0] - orig[i].a[0]) + (v[1] - u[1]) * (orig[i].b[1] - orig[i].a[1]) > 0;
        };
        // In the domain: by its side's direction, else the midpoint.
        std::array<bool, 2> in{false, false};  // left, right
        bool on_side = false;
        for (const std::size_t i : edges[e].src) {
            if (orig[i].shape != -1) continue;
            on_side = true;
            in = {along(i), !along(i)};
        }
        if (!on_side) {
            const bool inside = m[0] > g.lo[0] && m[0] < g.hi[0] && m[1] > g.lo[1] && m[1] < g.hi[1];
            in = {inside, inside};
        }
        std::array<std::int32_t, 2> r{in[0] ? g.background : outside, in[1] ? g.background : outside};
        std::array<bool, 2> done{!in[0], !in[1]};
        if (how != nullptr) {  // labelled: every shape holding each side, then the caller's label
            std::array<std::vector<std::size_t>, 2> hold;
            for (std::size_t s = 0; s < g.shapes.size(); ++s) {
                std::optional<std::size_t> own;
                for (const std::size_t i : edges[e].src) {
                    if (orig[i].shape == static_cast<std::int32_t>(s)) own = i;
                }
                std::array<bool, 2> inside{false, false};
                if (own) {
                    inside = {along(*own), !along(*own)};
                } else {
                    if (m[0] < box[s][0] || m[0] > box[s][1] || m[1] < box[s][2] || m[1] > box[s][3]) continue;
                    const bool k = detail::strictly_inside_2d(ccw[s], m);
                    inside = {k, k};
                }
                for (std::size_t side = 0; side < 2; ++side) {
                    if (inside[side]) hold[side].push_back(s);
                }
            }
            for (std::size_t side = 0; side < 2; ++side) {
                if (in[side]) r[side] = how->label(hold[side], m);
            }
            done = {true, true};
        }
        for (std::size_t o = 0; o < order.size() && !(done[0] && done[1]); ++o) {
            const std::size_t s = order[o];
            std::optional<std::size_t> own;
            for (const std::size_t i : edges[e].src) {
                if (orig[i].shape == static_cast<std::int32_t>(s)) own = i;
            }
            std::array<bool, 2> inside{false, false};
            if (own) {
                inside = {along(*own), !along(*own)};
            } else {
                if (m[0] < box[s][0] || m[0] > box[s][1] || m[1] < box[s][2] || m[1] > box[s][3]) continue;
                const bool k = detail::strictly_inside_2d(ccw[s], m);
                inside = {k, k};
            }
            for (std::size_t side = 0; side < 2; ++side) {
                if (!done[side] && inside[side]) {
                    r[side] = g.shapes[s].region;
                    done[side] = true;
                }
            }
        }
        for (std::size_t side = 0; side < 2; ++side) {
            if (in[side] && r[side] == -1) {
                return invalid("the layers leave part of the domain uncovered (near (" + std::to_string(m[0]) + ", " +
                               std::to_string(m[1]) + ")) and no region is the background");
            }
        }
        left[2 * e] = r[0];
        left[2 * e + 1] = r[1];
    }
    const auto boundary = [&](std::size_t h) { return left[h] >= 0 && (edges[h / 2].cut || left[h] != left[h ^ 1]); };

    // A region's loops: half-edges, the region on the left, each passing a vertex once.
    const auto trace = [&](std::int32_t region) {
        std::vector<std::vector<std::size_t>> found;
        std::vector<char> seen(left.size(), 0);
        for (std::size_t h0 = 0; h0 < left.size(); ++h0) {
            if (seen[h0] != 0 || left[h0] != region || !boundary(h0)) continue;
            std::vector<std::size_t> walk;
            for (std::size_t h = h0; seen[h] == 0;) {
                seen[h] = 1;
                walk.push_back(h);
                const std::size_t v = to(h), t = h ^ 1, k = out[v].size();
                std::size_t next = h0;
                for (std::size_t step = 1; step <= k; ++step) {
                    const std::size_t c = out[v][(pos[t] + k - step) % k];
                    if (left[c] == region && boundary(c)) {
                        next = c;
                        break;
                    }
                }
                h = next;
            }
            // Split where it passes a vertex twice.
            std::vector<std::size_t> path;
            std::map<std::size_t, std::size_t> at;
            for (const std::size_t h : walk) {
                const std::size_t o = from(h);
                if (const auto it = at.find(o); it != at.end()) {
                    const std::size_t k = it->second;
                    found.emplace_back(path.begin() + static_cast<std::ptrdiff_t>(k), path.end());
                    for (std::size_t q = k; q < path.size(); ++q) at.erase(from(path[q]));
                    path.resize(k);
                }
                at[o] = path.size();
                path.push_back(h);
            }
            found.push_back(std::move(path));
        }
        return found;
    };
    const auto points_of = [&](const std::vector<std::size_t>& loop) {
        std::vector<Point> p;
        for (const std::size_t h : loop) p.push_back(pts[from(h)]);
        return p;
    };
    const auto turns_at = [&](const std::vector<std::size_t>& loop, std::size_t k) {
        const std::size_t m = loop.size();
        return orient(pts[from(loop[(k + m - 1) % m])], pts[from(loop[k])], pts[to(loop[k])]) != 0;
    };
    const auto in_patch = [&](const Point& p) {
        for (const PatchSegment& s : g.patches) {
            const Point a{s.a[0], s.a[1], 0}, b{s.b[0], s.b[1], 0};
            if (detail::on_segment_2d(a, b, p) && !same_pt(a, p) && !same_pt(b, p)) return true;
        }
        return false;
    };
    const auto add_bridge = [&](std::size_t a, std::size_t b, std::int32_t region) {
        const std::size_t e = edges.size();
        edges.push_back({a, b, {}, true});
        left.push_back(region);
        left.push_back(region);
        pos.push_back(0);
        pos.push_back(0);
        out[a].push_back(2 * e);
        out[b].push_back(2 * e + 1);
        sort_out(a);
        sort_out(b);
    };
    const auto contains = [](const std::vector<Point>& loop, const Point& p) {
        return std::ranges::find_if(loop, [&](const Point& q) { return same_pt(p, q); }) != loop.end();
    };

    std::vector<std::int32_t> regions;
    for (const std::int32_t r : left) {
        if (r >= 0 && std::ranges::find(regions, r) == regions.end()) regions.push_back(r);
    }
    std::ranges::sort(regions);
    std::vector<std::pair<std::int32_t, std::vector<std::size_t>>> loops;
    for (const std::int32_t region : regions) {
        for (int round = 0;; ++round) {
            auto mine = trace(region);
            std::vector<std::size_t> holes;
            for (std::size_t k = 0; k < mine.size(); ++k) {
                if (detail::polygon_turn(points_of(mine[k])) < 0) holes.push_back(k);
            }
            if (holes.empty()) {
                for (auto& l : mine) loops.emplace_back(region, std::move(l));
                break;
            }
            if (round > static_cast<int>(mine.size()) + 64) return invalid("internal: a slanted layout's hole could not be cut");
            // The first hole (by its lowest vertex), the outer loop around it, and the part's other holes.
            const auto low = [&](std::size_t k) {
                const auto p = points_of(mine[k]);
                return *std::ranges::min_element(p, lower);
            };
            const std::size_t hk = *std::ranges::min_element(holes, [&](std::size_t a, std::size_t b) { return lower(low(a), low(b)); });
            const auto inside_of = [&](const std::vector<Point>& o, const std::vector<Point>& l) {
                for (const Point& p : l) {
                    if (!contains(o, p)) return detail::strictly_inside_2d(o, p);
                }
                return false;
            };
            const auto hole = points_of(mine[hk]);
            std::optional<std::size_t> outer;
            for (std::size_t k = 0; k < mine.size() && !outer; ++k) {
                if (std::ranges::find(holes, k) == holes.end() && inside_of(points_of(mine[k]), hole)) outer = k;
            }
            if (!outer) return invalid("internal: a slanted layout's hole has no outer loop");
            const auto outer_pts = points_of(mine[*outer]);  // part_pts[0] below
            std::vector<std::size_t> part{*outer};
            for (const std::size_t k : holes) {
                if (inside_of(outer_pts, points_of(mine[k]))) part.push_back(k);
            }
            // Bridge ends: the hole's corners, and the corners of the part's other loops (with their loop: both
            // bridges reach one loop, so the two divide the part).
            std::vector<std::vector<Point>> part_pts;
            for (const std::size_t k : part) part_pts.push_back(points_of(mine[k]));
            std::vector<std::size_t> starts;
            std::vector<std::pair<std::size_t, std::size_t>> ends;  // vertex, loop (its index in `part`)
            for (std::size_t q = 0; q < mine[hk].size(); ++q) {
                if (turns_at(mine[hk], q)) starts.push_back(from(mine[hk][q]));
            }
            for (std::size_t pk = 0; pk < part.size(); ++pk) {
                if (part[pk] == hk) continue;
                const auto& l = mine[part[pk]];
                for (std::size_t q = 0; q < l.size(); ++q) {
                    if (turns_at(l, q)) ends.emplace_back(from(l[q]), pk);
                }
            }
            const auto clear = [&](std::size_t a, std::size_t b, std::optional<std::pair<std::size_t, std::size_t>> extra) {
                const Point& pa = pts[a];
                const Point& pb = pts[b];
                const auto meets = [&](std::size_t x, std::size_t y) {
                    const bool xa = x == a || x == b, ya = y == a || y == b;
                    if (xa && ya) return true;  // the same edge
                    if (xa || ya) {
                        const std::size_t s = xa ? x : y, w = xa ? y : x, other = s == a ? b : a;
                        // Along the bridge, or so nearly (rounded crossings on one line): a sliver, refused.
                        const Point& q = pts[w];
                        const double bx = pts[other][0] - pts[s][0], by = pts[other][1] - pts[s][1];
                        const double qx = q[0] - pts[s][0], qy = q[1] - pts[s][1];
                        const double extent = std::max(g.hi[0] - g.lo[0], g.hi[1] - g.lo[1]);
                        if (bx * qx + by * qy > 0 && std::abs(bx * qy - by * qx) <= 1e-9 * std::hypot(bx, by) * extent) return true;
                        return orient(pa, pb, q) == 0 &&
                               (detail::on_segment_2d(pa, pb, q) || detail::on_segment_2d(pts[s], q, pts[other]));
                    }
                    return detail::segments_meet_2d(pa, pb, pts[x], pts[y]);
                };
                if (extra && meets(extra->first, extra->second)) return false;
                const Point m{0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0};
                if (!detail::strictly_inside_2d(part_pts[0], m)) return false;
                for (std::size_t pk = 1; pk < part.size(); ++pk) {
                    if (detail::strictly_inside_2d(part_pts[pk], m)) return false;
                }
                for (const std::size_t k : part) {
                    for (const std::size_t h : mine[k]) {
                        if (meets(from(h), to(h))) return false;
                    }
                }
                return true;
            };
            // Every candidate ranked first (on the domain's boundary off every patch's inside, then the shorter),
            // tested in that order: the first clear one, then the first clear to the same loop from another corner.
            struct Bridge {
                int score;
                double length;
                std::size_t a, b, loop;
            };
            std::vector<Bridge> ranked;
            for (const std::size_t a : starts) {
                for (const auto& [b, pk] : ends) {
                    if (a == b) continue;
                    const Point& p = pts[b];
                    const bool on_x = p[0] == g.lo[0] || p[0] == g.hi[0], on_y = p[1] == g.lo[1] || p[1] == g.hi[1];
                    const int score = in_patch(p) ? 16 : (on_x || on_y) ? (on_x && on_y ? 0 : 1) : 8;
                    const double dx = p[0] - pts[a][0], dy = p[1] - pts[a][1];
                    ranked.push_back({score, dx * dx + dy * dy, a, b, pk});
                }
            }
            std::ranges::sort(ranked, [](const Bridge& x, const Bridge& y) {
                return std::tie(x.score, x.length, x.a, x.b, x.loop) < std::tie(y.score, y.length, y.a, y.b, y.loop);
            });
            std::optional<Bridge> one, two;
            for (const Bridge& r : ranked) {
                if (!clear(r.a, r.b, std::nullopt)) continue;
                const auto other = std::ranges::find_if(ranked, [&](const Bridge& s) {
                    return s.loop == r.loop && s.a != r.a && s.b != r.b && clear(s.a, s.b, std::pair{r.a, r.b});
                });
                if (other == ranked.end()) continue;
                one = r;
                two = *other;
                break;
            }
            if (!one) return invalid("internal: a slanted layout's hole could not be cut");
            add_bridge(one->a, one->b, region);
            add_bridge(two->a, two->b, region);
        }
    }

    // The vertices kept: three or more boundary edges meet there, or its two come from lines that differ.
    std::vector<int> degree(pts.size(), 0);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        if (boundary(2 * e) || boundary(2 * e + 1)) {
            ++degree[edges[e].u];
            ++degree[edges[e].v];
        }
    }
    const auto line_of = [&](std::size_t h) -> std::pair<Point, Point> {
        const Edge& e = edges[h / 2];
        if (e.src.empty()) return {pts[e.u], pts[e.v]};
        return {orig[e.src[0]].a, orig[e.src[0]].b};
    };
    Geometry2D result;
    result.patches = g.patches;
    std::vector<std::pair<std::int32_t, std::vector<Point>>> polys;
    for (const auto& [region, loop] : loops) {
        std::vector<Point> p;
        const std::size_t m = loop.size();
        for (std::size_t k = 0; k < m; ++k) {
            const std::size_t v = from(loop[k]);
            bool keep = degree[v] >= 3 || std::ranges::any_of(pins, [&](const Point& q) { return same_pt(q, pts[v]); });
            if (!keep) {
                const auto [a, b] = line_of(loop[(k + m - 1) % m]);
                const auto [c, d] = line_of(loop[k]);
                keep = orient(a, b, c) != 0 || orient(a, b, d) != 0;
            }
            if (keep) p.push_back(pts[v]);
        }
        if (p.size() < 3) return invalid("internal: a slanted layout's polygon lost its corners");
        std::ranges::rotate(p, std::ranges::min_element(p, lower));
        if (!detail::simple_2d(p) || detail::polygon_turn(p) <= 0) {
            return invalid("the slanted layout has edges too close to resolve exactly (near (" + std::to_string(p[0][0]) +
                           ", " + std::to_string(p[0][1]) + "))");
        }
        polys.emplace_back(region, std::move(p));
    }
    std::ranges::sort(polys, [](const auto& x, const auto& y) {
        return lower(x.second[0], y.second[0]) || (same_pt(x.second[0], y.second[0]) && x.first < y.first);
    });
    for (auto& [region, p] : polys) result.polygons.push_back({std::move(p), region});
    return result;
}

}  // namespace

namespace detail {

std::expected<Geometry2D, base::Error> labelled_partition2d(const LayeredGeometry2D& shapes, const ArrangementLabels& how) {
    return general_partition(shapes, &how);
}

}  // namespace detail

std::expected<Geometry2D, base::Error> layered_partition2d(const LayeredGeometry2D& g) {
    for (int a = 0; a < 2; ++a) {
        const auto k = static_cast<std::size_t>(a);
        if (!(std::isfinite(g.lo[k]) && std::isfinite(g.hi[k]) && g.lo[k] < g.hi[k])) {
            return invalid("the layered domain must be finite with lo < hi");
        }
    }
    bool slanted = false;
    for (std::size_t s = 0; s < g.shapes.size(); ++s) {
        const auto& v = g.shapes[s].vertices;
        if (g.shapes[s].region < 0) return invalid("a layered shape's region must be >= 0", s);
        if (v.size() < 3) return invalid("a layered shape needs at least three vertices", s);
        for (std::size_t k = 0; k < v.size(); ++k) {
            const Point& p = v[k];
            const Point& q = v[(k + 1) % v.size()];
            if (!(std::isfinite(p[0]) && std::isfinite(p[1]))) return invalid("a layered shape's vertices must be finite", s);
            const bool horizontal = p[1] == q[1], vertical = p[0] == q[0];
            if (horizontal && vertical) return invalid("a layered shape repeats a vertex", s);
            slanted = slanted || (!horizontal && !vertical);
        }
        if (!detail::simple_2d(v) || detail::polygon_turn(v) == 0) return invalid("a layered shape must be simple, of non-zero area", s);
    }
    // Slice 4 (GD10): a slanted edge anywhere, the general arrangement; rectilinear layouts keep the grid below.
    if (slanted) return general_partition(g);

    std::vector<double> xs = breakpoints(g, 0), ys = breakpoints(g, 1);
    std::size_t nx = xs.size() - 1, ny = ys.size() - 1;
    // Each row painted from the top shape down, every cell once: `free_from` skips the cells already taken
    // (path-halving union-find over the row's columns).
    std::vector<std::int32_t> owner(nx * ny, g.background);
    std::vector<std::array<double, 2>> yrange(g.shapes.size());
    for (std::size_t s = 0; s < g.shapes.size(); ++s) {
        double lo = g.shapes[s].vertices[0][1], hi = lo;
        for (const Point& p : g.shapes[s].vertices) {
            lo = std::min(lo, p[1]);
            hi = std::max(hi, p[1]);
        }
        yrange[s] = {lo, hi};
    }
    const std::vector<std::size_t> order = paint_order(g);
    std::vector<std::size_t> cross, free_from(nx + 1);
    const auto find = [&](std::size_t i) {
        while (free_from[i] != i) {
            free_from[i] = free_from[free_from[i]];
            i = free_from[i];
        }
        return i;
    };
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) free_from[i] = i;
        for (const std::size_t s : order) {
            if (yrange[s][0] > ys[j] || yrange[s][1] < ys[j + 1]) continue;
            const auto& v = g.shapes[s].vertices;
            cross.clear();
            for (std::size_t k = 0; k < v.size(); ++k) {
                const Point& p = v[k];
                const Point& q = v[(k + 1) % v.size()];
                if (p[0] != q[0]) continue;  // horizontal
                const double lo = std::min(p[1], q[1]), hi = std::max(p[1], q[1]);
                if (lo <= ys[j] && hi >= ys[j + 1]) cross.push_back(index_of(xs, p[0]));
            }
            std::ranges::sort(cross);
            for (std::size_t c = 0; c + 1 < cross.size(); c += 2) {
                for (std::size_t i = find(cross[c]); i < cross[c + 1]; i = find(i + 1)) {
                    owner[j * nx + i] = g.shapes[s].region;
                    free_from[i] = i + 1;
                }
            }
        }
    }
    if (const auto it = std::ranges::find(owner, -1); it != owner.end()) {
        const auto c = static_cast<std::size_t>(it - owner.begin());
        return invalid("the layers leave part of the domain uncovered (near (" + std::to_string(xs[c % nx]) + ", " +
                       std::to_string(ys[c / nx]) + ")) and no region is the background");
    }
    // Breakpoint lines with the same owners on both sides along their whole length (hidden under later shapes)
    // carry no edge of any polygon: merged away, exactly, so the later passes see a smaller grid.
    {
        std::vector<char> keep_x(nx + 1, 0), keep_y(ny + 1, 0);
        keep_x[0] = keep_x[nx] = keep_y[0] = keep_y[ny] = 1;
        for (std::size_t j = 0; j < ny; ++j) {
            const std::int32_t* r = owner.data() + j * nx;
            for (std::size_t i = 1; i < nx; ++i) keep_x[i] |= static_cast<char>(r[i - 1] != r[i]);
            if (j > 0) {
                const std::int32_t* below = r - nx;
                for (std::size_t i = 0; i < nx && keep_y[j] == 0; ++i) keep_y[j] = static_cast<char>(below[i] != r[i]);
            }
        }
        std::vector<std::size_t> cx, cy;  // the first old column (row) of each new one
        std::vector<double> x2, y2;
        for (std::size_t i = 0; i <= nx; ++i) {
            if (keep_x[i] == 0) continue;
            x2.push_back(xs[i]);
            if (i < nx) cx.push_back(i);
        }
        for (std::size_t j = 0; j <= ny; ++j) {
            if (keep_y[j] == 0) continue;
            y2.push_back(ys[j]);
            if (j < ny) cy.push_back(j);
        }
        if (cx.size() < nx || cy.size() < ny) {
            std::vector<std::int32_t> o2(cx.size() * cy.size());
            for (std::size_t j = 0; j < cy.size(); ++j) {
                for (std::size_t i = 0; i < cx.size(); ++i) o2[j * cx.size() + i] = owner[cy[j] * nx + cx[i]];
            }
            owner = std::move(o2);
            xs = std::move(x2);
            ys = std::move(y2);
            nx = cx.size();
            ny = cy.size();
        }
    }
    // Each cell's column and row, looked up (a division per cell per pass costs more than the pass).
    std::vector<std::uint32_t> col(nx * ny), row(nx * ny);
    for (std::size_t j = 0, c = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i, ++c) {
            col[c] = static_cast<std::uint32_t>(i);
            row[c] = static_cast<std::uint32_t>(j);
        }
    }

    // Parts: same-region cells sharing an edge, labelled by a flood from each unlabelled cell in row-major order.
    // A part's boundary is its directed unit edges with the part on the left, found for every part in one pass.
    constexpr std::uint32_t none = UINT32_MAX;
    std::vector<std::uint32_t> label(nx * ny, none);
    std::vector<std::size_t> stack;
    struct PartInfo {
        std::int32_t region;
        std::size_t seed;                                       // a cell of the part
        std::vector<std::pair<std::size_t, std::size_t>> edges;  // grid vertex keys
    };
    std::vector<PartInfo> parts;
    const std::size_t vx = nx + 1;
    const auto key = [&](std::size_t i, std::size_t j) { return j * vx + i; };
    // Flood from c0 over its region's cells with label `from` whose column lies in [lo, hi): they take `to`. A span
    // fill: each popped seed's row run is taken whole, and one seed is pushed per run of matches above and below.
    const auto flood = [&](std::size_t c0, std::uint32_t from, std::uint32_t to, std::size_t lo, std::size_t hi,
                           std::vector<std::size_t>* cells) {
        const std::int32_t region = owner[c0];
        const auto match = [&](std::size_t n) { return label[n] == from && owner[n] == region; };
        stack.assign(1, c0);
        while (!stack.empty()) {
            std::size_t c = stack.back();
            stack.pop_back();
            if (!match(c)) continue;
            const std::size_t j = row[c], base = j * nx;
            std::size_t a = col[c];
            while (a > lo && match(base + a - 1)) --a;
            std::size_t b = a;
            while (b < hi && match(base + b)) {
                label[base + b] = to;
                if (cells) cells->push_back(base + b);
                ++b;
            }
            for (const std::size_t r : {j - 1, j + 1}) {  // j - 1 wraps past 0 when j == 0: skipped as >= ny
                if (r >= ny) continue;
                bool run = false;
                for (std::size_t x = a; x < b; ++x) {
                    const bool m = match(r * nx + x);
                    if (m && !run) stack.push_back(r * nx + x);
                    run = m;
                }
            }
        }
    };
    for (std::size_t c = 0; c < nx * ny; ++c) {
        if (label[c] != none) continue;
        parts.push_back({owner[c], c, {}});
        flood(c, none, static_cast<std::uint32_t>(parts.size() - 1), 0, nx, nullptr);
    }
    // The boundary edges of the cells given (all of them, or a cut part's).
    const auto boundary = [&](std::size_t c) {
        const std::uint32_t l = label[c];
        auto& e = parts[l].edges;
        const std::size_t i = col[c], j = row[c];
        if (j == 0 || label[c - nx] != l) e.emplace_back(key(i, j), key(i + 1, j));
        if (i + 1 == nx || label[c + 1] != l) e.emplace_back(key(i + 1, j), key(i + 1, j + 1));
        if (j + 1 == ny || label[c + nx] != l) e.emplace_back(key(i + 1, j + 1), key(i, j + 1));
        if (i == 0 || label[c - 1] != l) e.emplace_back(key(i, j + 1), key(i, j));
    };
    for (std::size_t c = 0; c < nx * ny; ++c) boundary(c);

    // A hole's cut (a ring has no simple outline, Unit 26k slice 1): two rays along grid lines from the hole's
    // corners out to the part's outer loop. A straight cut across the whole part would also meet the domain's
    // boundary, and a vertex inside a patch there breaks generate2d's rule that a patch lies inside one boundary
    // edge, and one ending on an interface makes a corner zone generate2d's P4b check can refuse; so a ray ends, best
    // first: on the domain's boundary off every patch's inside (a corner, then not), on an interface's corner, on an
    // interface, inside a patch; then the shorter. The
    // part's cells, walled along the rays, fall into pieces, each a new part (traced again: other holes, pinches).
    // 1: cut; 0: fewer than two rays (nothing changed); -1: the rays did not divide the part.
    std::vector<char> on_outer((nx + 1) * (ny + 1), 0), vwall(nx * ny, 0), hwall(nx * ny, 0);
    std::vector<std::size_t> cells;
    const auto in_patch = [&](std::size_t v) {
        const std::size_t i = v % vx, j = v / vx;
        if (i != 0 && i != nx && j != 0 && j != ny) return false;
        const double x = xs[i], y = ys[j];
        for (const PatchSegment& p : g.patches) {
            if (p.a[0] == x && p.b[0] == x && y > std::min(p.a[1], p.b[1]) && y < std::max(p.a[1], p.b[1])) return true;
            if (p.a[1] == y && p.b[1] == y && x > std::min(p.a[0], p.b[0]) && x < std::max(p.a[0], p.b[0])) return true;
        }
        return false;
    };
    const auto corner_at = [](const std::vector<std::size_t>& l, std::size_t k) {
        const std::size_t n = l.size();
        const auto d1 = static_cast<std::ptrdiff_t>(l[k]) - static_cast<std::ptrdiff_t>(l[(k + n - 1) % n]);
        const auto d2 = static_cast<std::ptrdiff_t>(l[(k + 1) % n]) - static_cast<std::ptrdiff_t>(l[k]);
        return d1 != d2;
    };
    const auto ray_cut = [&](std::size_t q, const std::vector<std::vector<std::size_t>>& ls) -> int {
        const auto lq = static_cast<std::uint32_t>(q);
        const auto mine_cell = [&](std::ptrdiff_t i, std::ptrdiff_t j) {
            return i >= 0 && j >= 0 && i < static_cast<std::ptrdiff_t>(nx) && j < static_cast<std::ptrdiff_t>(ny) &&
                   label[static_cast<std::size_t>(j) * nx + static_cast<std::size_t>(i)] == lq;
        };
        // One unit step from (i, j) in direction d (0 up, 1 down, 2 right, 3 left): whether both cells along it are
        // the part's; `wall` the wall it sets (`vertical`: a wall on the cell's left side, else on its bottom).
        const auto step = [&](std::ptrdiff_t& i, std::ptrdiff_t& j, int d, bool set) {
            std::ptrdiff_t ai = i - 1, aj = j, bi = i, bj = j, ni = i, nj = j + 1;
            if (d == 1) { aj = j - 1; bj = j - 1; nj = j - 1; }
            if (d == 2) { ai = i; aj = j - 1; nj = j; ni = i + 1; }
            if (d == 3) { ai = i - 1; aj = j - 1; bi = i - 1; bj = j; ni = i - 1; nj = j; }
            if (!mine_cell(ai, aj) || !mine_cell(bi, bj)) return false;
            if (set) {
                const std::size_t wc = static_cast<std::size_t>(bj) * nx + static_cast<std::size_t>(bi);
                (d < 2 ? vwall : hwall)[wc] = 1;
            }
            i = ni;
            j = nj;
            return true;
        };
        const auto& outer = ls[0];
        const auto& hole = ls[1];
        for (std::size_t k = 0; k < outer.size(); ++k) on_outer[outer[k]] = corner_at(outer, k) ? 2 : 1;
        struct Ray {
            std::size_t start;
            int dir, score;
            std::size_t steps;
        };
        std::vector<Ray> found;
        for (std::size_t k = 0; k < hole.size(); ++k) {
            if (!corner_at(hole, k)) continue;
            for (int d = 0; d < 4; ++d) {
                auto i = static_cast<std::ptrdiff_t>(hole[k] % vx), j = static_cast<std::ptrdiff_t>(hole[k] / vx);
                std::size_t n = 0;
                while (step(i, j, d, false)) ++n;
                const std::size_t end = static_cast<std::size_t>(j) * vx + static_cast<std::size_t>(i);
                if (n == 0 || on_outer[end] == 0) continue;
                const std::size_t ei = end % vx, ej = end / vx;
                const bool wall = ei == 0 || ei == nx || ej == 0 || ej == ny, corner = on_outer[end] == 2;
                const int score = in_patch(end) ? 16 : wall ? (corner ? 0 : 1) : (corner ? 4 : 8);
                found.push_back({hole[k], d, score, n});
            }
        }
        for (const std::size_t v : outer) on_outer[v] = 0;
        std::ranges::sort(found, [](const Ray& a, const Ray& b) { return std::tie(a.score, a.steps, a.start, a.dir) < std::tie(b.score, b.steps, b.start, b.dir); });
        // The best ray, and the best from another corner (two from one corner leave the rest touching itself there).
        const auto second = std::ranges::find_if(found, [&](const Ray& r) { return found.empty() || r.start != found.front().start; });
        if (found.empty() || second == found.end()) return 0;
        for (const Ray* r : {&found.front(), &*second}) {
            auto i = static_cast<std::ptrdiff_t>(r->start % vx), j = static_cast<std::ptrdiff_t>(r->start / vx);
            while (step(i, j, r->dir, true)) {}
        }
        // The pieces: a plain flood that stops at walls.
        const std::uint32_t work = none - 1;
        cells.clear();
        flood(parts[q].seed, lq, work, 0, nx, &cells);
        const auto [lo_c, hi_c] = std::ranges::minmax(cells);
        const std::size_t first = parts.size();
        for (std::size_t c0 = lo_c; c0 <= hi_c; ++c0) {
            if (label[c0] != work) continue;
            const auto to = static_cast<std::uint32_t>(parts.size());
            parts.push_back({parts[q].region, c0, {}});
            stack.assign(1, c0);
            label[c0] = to;
            while (!stack.empty()) {
                const std::size_t c = stack.back();
                stack.pop_back();
                const std::size_t i = col[c], j = row[c];
                const auto visit = [&](std::size_t n) {
                    if (label[n] == work) {
                        label[n] = to;
                        stack.push_back(n);
                    }
                };
                if (i > 0 && vwall[c] == 0) visit(c - 1);
                if (i + 1 < nx && vwall[c + 1] == 0) visit(c + 1);
                if (j > 0 && hwall[c] == 0) visit(c - nx);
                if (j + 1 < ny && hwall[c + nx] == 0) visit(c + nx);
            }
        }
        for (const std::size_t c : cells) vwall[c] = hwall[c] = 0;
        if (parts.size() - first < 2) return -1;
        for (const std::size_t c : cells) boundary(c);
        return 1;
    };

    // Grid vertex (i, j) as one key; a loop as its unit steps' vertices.
    std::vector<std::pair<std::int32_t, std::vector<std::size_t>>> loops;  // region, vertices (every unit step)
    constexpr std::size_t npos = SIZE_MAX;
    std::vector<std::size_t> next((nx + 1) * (ny + 1), npos);  // a boundary vertex's outgoing edge's end
    std::vector<char> seen((nx + 1) * (ny + 1), 0);
    std::vector<std::size_t> touched, cuts;
    for (std::size_t q = 0; q < parts.size(); ++q) {
        // A vertex with two outgoing edges is a pinch: a cut along its column line.
        touched.clear();
        cuts.clear();
        for (const auto& [a, b] : parts[q].edges) {
            if (next[a] != npos) {
                cuts.push_back(a % vx);
                continue;
            }
            next[a] = b;
            touched.push_back(a);
        }
        std::vector<std::vector<std::size_t>> mine;
        int rays = 0;
        if (cuts.empty()) {
            // Loops from their lowest-then-leftmost vertex (keys order by row, then column); the first is the outer
            // one, each other a hole: the first hole cut by two rays (ray_cut), or, without two, every hole by a cut
            // along its leftmost line.
            std::ranges::sort(touched);
            for (const std::size_t s : touched) {
                if (seen[s] != 0) continue;
                std::vector<std::size_t> loop;
                for (std::size_t v = s; seen[v] == 0; v = next[v]) {
                    seen[v] = 1;
                    loop.push_back(v);
                }
                mine.push_back(std::move(loop));
            }
            if (mine.size() > 1) {
                rays = ray_cut(q, mine);
                if (rays < 0) return invalid("internal: a layered part's cut did not divide it");
                for (std::size_t k = 1; rays == 0 && k < mine.size(); ++k) {
                    std::size_t left = npos;
                    for (const std::size_t v : mine[k]) left = std::min(left, v % vx);
                    cuts.push_back(left);
                }
            }
        }
        for (const std::size_t v : touched) {
            next[v] = npos;
            seen[v] = 0;
        }
        parts[q].edges.clear();
        parts[q].edges.shrink_to_fit();
        if (rays > 0) continue;
        if (cuts.empty()) {
            loops.emplace_back(parts[q].region, std::move(mine.front()));
            continue;
        }
        // Every cut at once: the part's cells split into the bands between cut lines, each band's pieces new parts
        // (cutting adds no hole or pinch), their edges found again.
        std::ranges::sort(cuts);
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        const std::uint32_t old = static_cast<std::uint32_t>(q), work = none - 1;
        cells.clear();
        flood(parts[q].seed, old, work, 0, nx, &cells);
        // In row-major order without a sort: a scan of the part's index range.
        const auto [lo_c, hi_c] = std::ranges::minmax(cells);
        const std::size_t first = parts.size();
        for (std::size_t c = lo_c; c <= hi_c; ++c) {
            if (label[c] != work) continue;
            const auto band = static_cast<std::size_t>(std::ranges::upper_bound(cuts, col[c]) - cuts.begin());
            const std::size_t lo = band == 0 ? 0 : cuts[band - 1], hi = band == cuts.size() ? nx : cuts[band];
            parts.push_back({parts[q].region, c, {}});
            flood(c, work, static_cast<std::uint32_t>(parts.size() - 1), lo, hi, nullptr);
        }
        if (parts.size() - first < 2) return invalid("internal: a layered part's cut did not divide it");
        for (const std::size_t c : cells) boundary(c);
    }
    // The corners of every loop; each polygon keeps the vertices along its boundary that are any loop's corner.
    std::vector<char> corners((nx + 1) * (ny + 1), 0);
    const auto is_corner = [&](const std::vector<std::size_t>& l, std::size_t k) {
        const std::size_t a = l[(k + l.size() - 1) % l.size()], b = l[k], c = l[(k + 1) % l.size()];
        const bool ab_h = a / vx == b / vx, bc_h = b / vx == c / vx;
        return ab_h != bc_h;
    };
    for (const auto& [region, l] : loops) {
        for (std::size_t k = 0; k < l.size(); ++k) {
            if (is_corner(l, k)) corners[l[k]] = 1;
        }
    }
    Geometry2D out;
    out.patches = g.patches;
    for (const auto& [region, l] : loops) {
        RegionPolygon poly;
        poly.region = region;
        for (const std::size_t v : l) {
            if (corners[v] != 0) poly.vertices.push_back({xs[v % vx], ys[v / vx], 0.0});
        }
        out.polygons.push_back(std::move(poly));
    }
    return out;
}

}  // namespace NiTCAD::mesh
