// Layered 3D solids resolved into a polyhedral complex (Unit 26l, SD1-SD2; complex.hpp).
//
// The solids' coordinates on each axis (and the patch boxes') make a grid. Each cell goes to the topmost solid covering
// its centre: the solids are painted top first (`fill` ones after, bottom first: Unit 26m MT4), each row's unpainted cells found by a skip list, so every cell is
// painted once. Cells of one region merge into blocks (runs along x, then y, then z), and each block is a volume. On
// each grid plane the faces of the blocks on its two sides are overlaid: every overlay cell knows the block below and
// above it (and, on the domain's boundary, its patch), and cells alike merge into rectangular facets. Each facet's
// outline then takes every vertex lying on it, so facets meet edge to edge. All of it is on grid indices, so it is
// exact; coordinates are looked up at the end.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "arrangement2d.hpp"
#include "polygon2d.hpp"

namespace NiTCAD::mesh {

namespace {

std::unexpected<base::Error> invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return std::unexpected(base::Error{base::ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}});
}

// The order a cell asks the solids in, the first holding it taking it (Unit 26m MT4): the replacing solids top
// first, then the filling ones bottom first. With no `fill` solid, top first as before.
std::vector<std::size_t> paint_order(const LayeredGeometry3D& g) {
    std::vector<std::size_t> order;
    order.reserve(g.solids.size());
    for (std::size_t s = g.solids.size(); s-- > 0;) {
        if (!g.solids[s].fill) order.push_back(s);
    }
    for (std::size_t s = 0; s < g.solids.size(); ++s) {
        if (g.solids[s].fill) order.push_back(s);
    }
    return order;
}

using Idx = std::array<int, 3>;

// The breakpoints of one axis: the domain's ends and every coordinate inside it, coordinates within `tol` of one
// another one breakpoint (the first in order, so 20nm and 0.02um are one plane).
std::vector<double> merged(std::vector<double> v, double lo, double hi, double tol) {
    std::vector<double> in{lo, hi};
    for (const double c : v) {
        if (c > lo + tol && c < hi - tol) in.push_back(c);
    }
    std::ranges::sort(in);
    std::vector<double> out;
    for (const double c : in) {
        if (out.empty() || c - out.back() > tol) out.push_back(c);
    }
    out.back() = hi;
    return out;
}

// The breakpoint nearest `c`, clamped to the axis.
int index_of(const std::vector<double>& b, double c) {
    if (c <= b.front()) return 0;
    if (c >= b.back()) return static_cast<int>(b.size()) - 1;
    const auto it = std::ranges::lower_bound(b, c);
    const auto k = static_cast<int>(it - b.begin());
    return (c - b[static_cast<std::size_t>(k - 1)] < b[static_cast<std::size_t>(k)] - c) ? k - 1 : k;
}

struct Block {
    Idx lo, hi;  // cell index ranges [lo, hi)
    std::int32_t region;
};

struct Facet {
    int axis, at;          // on the plane of breakpoint `at` across `axis`
    int u0, u1, v0, v1;    // breakpoint ranges along (axis + 1) % 3 and (axis + 2) % 3
    int below, above;      // blocks (-1: outside the domain)
    int patch;             // a boundary facet's patch box (-1: its face's name; -2: internal)
};

// --- slanted prisms (slice 4, SD6) ----------------------------------------------------------------------------------
//
// Every slanted face is parallel to one axis A (SQ4). Cut across A at its breakpoints the layout is a stack of slabs,
// each a 2D layout in (u, v), the axes after A. One planar arrangement (arrangement2d.hpp) serves every slab: each
// point is labelled by its signature (its region in every slab, and the patch box holding it on each end plane), so
// the arrangement's polygons are the pieces across which nothing changes along A, and the sectioned extrusion of those
// polygons through the breakpoints (`extrude`) is the complex. Crossings are computed once and shared by every slab.
// The lateral facets are polygon edges swept through one slab: vertical quads, exactly planar even through a rounded
// crossing (two such quads meet at a fold of about an ulp where a slanted face is split there). Coordinates within
// the tolerance of one another on u or v are first made one (20nm and 0.02um), as the rectilinear grid does.
std::expected<PolyhedralGeometry, base::Error> slanted_partition(const LayeredGeometry3D& g, int axis, const std::vector<double>& t, double tol) {
    const auto A = static_cast<std::size_t>(axis), ua = (A + 1) % 3, va = (A + 2) % 3;
    const std::size_t ns = t.size() - 1;
    struct Entry {
        std::vector<Point> polygon;  // (u, v, 0)
        int solid = -1;              // a solid's cross-section, or (solid -1) a patch box on an end plane:
        int patch = -1, end = 0;     // its index, and 0 for the low end plane, 1 for the high one
        std::vector<char> slabs;     // a cross-section: the slabs it holds in
    };
    std::vector<Entry> entries;
    const auto rect = [](double u0, double u1, double v0, double v1) {
        return std::vector<Point>{{u0, v0, 0.0}, {u1, v0, 0.0}, {u1, v1, 0.0}, {u0, v1, 0.0}};
    };
    std::map<std::pair<int, std::vector<std::pair<double, double>>>, std::size_t> entry_of;
    const auto add = [&](int solid, std::vector<Point> polygon, int k0, int k1) {
        std::vector<std::pair<double, double>> key;
        for (const Point& p : polygon) key.emplace_back(p[0], p[1]);
        const auto [it, fresh] = entry_of.try_emplace({solid, std::move(key)}, entries.size());
        if (fresh) entries.push_back({std::move(polygon), solid, -1, 0, std::vector<char>(ns, 0)});
        for (int k = k0; k < k1; ++k) entries[it->second].slabs[static_cast<std::size_t>(k)] = 1;
    };
    for (std::size_t s = 0; s < g.solids.size(); ++s) {
        const LayeredSolid& v = g.solids[s];
        const int si = static_cast<int>(s);
        if (v.kind == LayeredSolid::Kind::box) {
            const int k0 = index_of(t, v.lo[A]), k1 = index_of(t, v.hi[A]);
            if (k0 < k1) add(si, rect(v.lo[ua], v.hi[ua], v.lo[va], v.hi[va]), k0, k1);
            continue;
        }
        if (v.axis == axis) {
            const int k0 = index_of(t, v.from), k1 = index_of(t, v.to);
            std::vector<Point> polygon;
            for (const Point& p : v.polygon) polygon.push_back({p[0], p[1], 0.0});
            if (k0 < k1) add(si, std::move(polygon), k0, k1);
            continue;
        }
        // A rectilinear prism across A: in each slab, a scanline at the slab's centre across its polygon gives
        // intervals along the polygon's other coordinate w, each swept along the prism's axis b: rectangles.
        const auto b = static_cast<std::size_t>(v.axis);
        const std::size_t pa = (b + 1) % 3 == A ? 0 : 1, pw = 1 - pa;
        const std::size_t n = v.polygon.size();
        std::vector<int> at(n);
        for (std::size_t i = 0; i < n; ++i) at[i] = index_of(t, v.polygon[i][pa]);
        const auto [kmin, kmax] = std::ranges::minmax(at);
        std::vector<double> cuts;
        for (int k = kmin; k < kmax; ++k) {
            cuts.clear();
            for (std::size_t e = 0; e < n; ++e) {
                const std::size_t f = (e + 1) % n;
                if (v.polygon[e][pw] != v.polygon[f][pw]) continue;  // along w: not a crossing
                if ((at[e] <= k) != (at[f] <= k)) cuts.push_back(v.polygon[e][pw]);
            }
            std::ranges::sort(cuts);
            for (std::size_t c = 0; c + 1 < cuts.size(); c += 2) {
                if (!(cuts[c] < cuts[c + 1])) continue;
                add(si, ua == b ? rect(v.from, v.to, cuts[c], cuts[c + 1]) : rect(cuts[c], cuts[c + 1], v.from, v.to), k, k + 1);
            }
        }
    }
    // Patch boxes: on an end plane, a rectangle the arrangement splits along; on a side of the domain, pins where the
    // box's edges cross that side (its extent along A is among the breakpoints already).
    std::vector<Point> pins;
    for (std::size_t p = 0; p < g.patches.size(); ++p) {
        const PatchBox& b = g.patches[p];
        for (const int end : {0, 1}) {
            const double plane = end == 1 ? t[ns] : t[0];
            if (b.lo[A] - tol <= plane && plane <= b.hi[A] + tol && b.lo[ua] < b.hi[ua] && b.lo[va] < b.hi[va]) {
                entries.push_back({rect(b.lo[ua], b.hi[ua], b.lo[va], b.hi[va]), -1, static_cast<int>(p), end, {}});
            }
        }
        for (const std::size_t c : {ua, va}) {
            const std::size_t d = c == ua ? va : ua;
            for (const double side : {g.lo[c], g.hi[c]}) {
                if (!(b.lo[c] - tol <= side && side <= b.hi[c] + tol)) continue;
                for (const double w : {b.lo[d], b.hi[d]}) {
                    if (!(w > g.lo[d] + tol && w < g.hi[d] - tol)) continue;
                    Point q{0.0, 0.0, 0.0};
                    q[c == ua ? 0 : 1] = side;
                    q[c == ua ? 1 : 0] = w;
                    pins.push_back(q);
                }
            }
        }
    }
    // Coordinates within the tolerance of one another made one, per axis (the domain's ends kept exactly).
    for (const std::size_t c : {std::size_t{0}, std::size_t{1}}) {
        const std::size_t ax = c == 0 ? ua : va;
        std::vector<double> all{g.lo[ax], g.hi[ax]};
        for (const Entry& e : entries) {
            for (const Point& p : e.polygon) all.push_back(p[c]);
        }
        for (const Point& p : pins) all.push_back(p[c]);
        std::ranges::sort(all);
        std::vector<double> reps;
        for (const double x : all) {
            if (reps.empty() || x - reps.back() > tol) reps.push_back(x);
        }
        for (double& r : reps) {
            if (std::abs(r - g.lo[ax]) <= tol) r = g.lo[ax];
            else if (std::abs(r - g.hi[ax]) <= tol) r = g.hi[ax];
        }
        const auto snap = [&](double x) {
            const auto it = std::ranges::lower_bound(reps, x - tol);
            return it != reps.end() && std::abs(*it - x) <= tol ? *it : x;
        };
        for (Entry& e : entries) {
            for (Point& p : e.polygon) p[c] = snap(p[c]);
        }
        for (Point& p : pins) p[c] = snap(p[c]);
    }

    // The arrangement, labelled by signature.
    std::map<std::vector<std::int32_t>, std::int32_t> ids;
    std::vector<std::vector<std::int32_t>> sigs;
    std::optional<std::pair<Point, std::size_t>> uncovered;
    detail::ArrangementLabels how;
    how.pins = pins;
    how.label = [&](const std::vector<std::size_t>& hold, const Point& at) -> std::int32_t {
        std::vector<int> top(ns, -1), first_fill(ns, -1);  // the top replacing solid, the bottom filling one (MT4)
        std::array<int, 2> cap{-1, -1};
        for (const std::size_t i : hold) {
            const Entry& e = entries[i];
            if (e.solid < 0) {
                if (cap[static_cast<std::size_t>(e.end)] < 0 || e.patch < cap[static_cast<std::size_t>(e.end)]) cap[static_cast<std::size_t>(e.end)] = e.patch;
                continue;
            }
            const bool fill = g.solids[static_cast<std::size_t>(e.solid)].fill;
            for (std::size_t k = 0; k < ns; ++k) {
                if (e.slabs[k] == 0) continue;
                if (!fill && e.solid > top[k]) top[k] = e.solid;
                if (fill && (first_fill[k] < 0 || e.solid < first_fill[k])) first_fill[k] = e.solid;
            }
        }
        std::vector<std::int32_t> sig(ns + 2);
        for (std::size_t k = 0; k < ns; ++k) {
            const int s = top[k] >= 0 ? top[k] : first_fill[k];
            const std::int32_t r = s >= 0 ? g.solids[static_cast<std::size_t>(s)].region : g.background;
            if (r < 0) {
                if (!uncovered) uncovered = std::pair{at, k};
                return -1;
            }
            sig[k] = r;
        }
        sig[ns] = cap[0];
        sig[ns + 1] = cap[1];
        const auto [it, fresh] = ids.try_emplace(sig, static_cast<std::int32_t>(sigs.size()));
        if (fresh) sigs.push_back(std::move(sig));
        return it->second;
    };
    LayeredGeometry2D flat;
    flat.lo = {g.lo[ua], g.lo[va], 0.0};
    flat.hi = {g.hi[ua], g.hi[va], 0.0};
    for (const Entry& e : entries) flat.shapes.push_back({e.polygon, 0});
    auto tiling = detail::labelled_partition2d(flat, how);
    if (!tiling) {
        if (uncovered) {
            Point p{};
            p[A] = 0.5 * (t[uncovered->second] + t[uncovered->second + 1]);
            p[ua] = uncovered->first[0];
            p[va] = uncovered->first[1];
            return invalid(std::format("the layers leave part of the domain uncovered (near ({:.6f}, {:.6f}, {:.6f})) and no region "
                                       "is the background",
                                       p[0], p[1], p[2]));
        }
        return std::unexpected(std::move(tiling.error()));
    }

    // The sectioned extrusion of the pieces through the breakpoints.
    static constexpr const char* unnamed = " face";  // never a deck name (those are identifiers)
    ExtrudedGeometry x;
    x.axis = axis;
    x.stations = t;
    x.default_patch = unnamed;
    x.patches = g.patches;
    for (std::size_t i = 0; i < tiling->polygons.size(); ++i) {
        const RegionPolygon& p = tiling->polygons[i];
        x.cells.push_back({"piece_" + std::to_string(i), p.vertices});
        const auto& sig = sigs[static_cast<std::size_t>(p.region)];
        x.regions.emplace_back(sig.begin(), sig.begin() + static_cast<std::ptrdiff_t>(ns));
    }
    auto built = extrude(x);
    if (!built) return std::unexpected(std::move(built.error()));
    static constexpr std::array<std::array<const char*, 2>, 3> face_names{{{"x_min", "x_max"}, {"y_min", "y_max"}, {"z_min", "z_max"}}};
    for (std::size_t f = 0; f < built->facets.size(); ++f) {
        ComplexFacet& c = built->facets[f];
        if (c.patch != unnamed) continue;
        c.patch.clear();
        for (std::size_t k = 0; k < 3 && c.patch.empty(); ++k) {
            for (const int side : {0, 1}) {
                const double w = side == 0 ? g.lo[k] : g.hi[k];
                if (std::ranges::all_of(c.vertices, [&](std::size_t i) { return built->vertices[i][k] == w; })) {
                    c.patch = face_names[k][static_cast<std::size_t>(side)];
                    break;
                }
            }
        }
        if (c.patch.empty()) return invalid("internal: a slanted layout's boundary facet lies on no face of the domain", f);
    }
    return built;
}

}  // namespace

std::expected<PolyhedralGeometry, base::Error> layered_partition3d(const LayeredGeometry3D& g) {
    for (std::size_t k = 0; k < 3; ++k) {
        if (!std::isfinite(g.lo[k]) || !std::isfinite(g.hi[k]) || !(g.lo[k] < g.hi[k])) {
            return invalid("the layered domain needs lo < hi on each axis");
        }
    }
    if (g.background < -1) return invalid("the background region is >= 0 (or -1: none)");
    double extent = 0.0;
    for (std::size_t k = 0; k < 3; ++k) extent = std::max({extent, std::abs(g.lo[k]), std::abs(g.hi[k]), g.hi[k] - g.lo[k]});
    const double tol = 1e-12 * extent;

    // Validation, and every coordinate per axis.
    std::array<std::vector<double>, 3> coords;
    int slant_axis = -1;  // the axis of the prisms with a slanted edge (slice 4, SD6): one at most
    for (std::size_t s = 0; s < g.solids.size(); ++s) {
        const LayeredSolid& v = g.solids[s];
        if (v.region < 0) return invalid("a layered solid's region is >= 0", s);
        if (v.kind == LayeredSolid::Kind::box) {
            for (std::size_t k = 0; k < 3; ++k) {
                if (!std::isfinite(v.lo[k]) || !std::isfinite(v.hi[k]) || !(v.lo[k] < v.hi[k])) {
                    return invalid("a layered box needs lo < hi on each axis", s);
                }
                coords[k].push_back(v.lo[k]);
                coords[k].push_back(v.hi[k]);
            }
            continue;
        }
        if (v.axis < 0 || v.axis > 2) return invalid("a prism's axis is x, y or z", s);
        if (!std::isfinite(v.from) || !std::isfinite(v.to) || !(v.from < v.to)) return invalid("a prism needs from < to along its axis", s);
        const auto a = static_cast<std::size_t>(v.axis), ua = (a + 1) % 3, va = (a + 2) % 3;
        const std::size_t n = v.polygon.size();
        if (n < 3) return invalid("a prism's polygon needs at least three vertices", s);
        std::vector<Point> flat;
        for (std::size_t i = 0; i < n; ++i) {
            const Point& p = v.polygon[i];
            const Point& q = v.polygon[(i + 1) % n];
            if (!std::isfinite(p[0]) || !std::isfinite(p[1])) return invalid("a prism's vertex is not finite", s);
            for (std::size_t j = 0; j < i; ++j) {
                if (v.polygon[j][0] == p[0] && v.polygon[j][1] == p[1]) return invalid("a prism's polygon repeats a vertex", s);
            }
            if (p[0] != q[0] && p[1] != q[1]) {
                if (slant_axis >= 0 && slant_axis != v.axis) {
                    return invalid("slanted prisms along two axes (every slanted face must be parallel to one axis, the one all "
                                   "slanted prisms share)",
                                   s);
                }
                slant_axis = v.axis;
            }
            flat.push_back({p[0], p[1], 0.0});
            coords[ua].push_back(p[0]);
            coords[va].push_back(p[1]);
        }
        if (!detail::simple_2d(flat) || detail::polygon_turn(flat) == 0) return invalid("a prism's polygon is not simple", s);
        coords[a].push_back(v.from);
        coords[a].push_back(v.to);
    }
    for (std::size_t b = 0; b < g.patches.size(); ++b) {
        const PatchBox& p = g.patches[b];
        if (p.name.empty()) return invalid("a patch box needs a name", b);
        for (std::size_t k = 0; k < 3; ++k) {
            if (!(p.lo[k] <= p.hi[k])) return invalid("a patch box needs lo <= hi on each axis", b);
            coords[k].push_back(p.lo[k]);
            coords[k].push_back(p.hi[k]);
        }
    }
    std::array<std::vector<double>, 3> bp;
    for (std::size_t k = 0; k < 3; ++k) bp[k] = merged(std::move(coords[k]), g.lo[k], g.hi[k], tol);
    if (slant_axis >= 0) return slanted_partition(g, slant_axis, bp[static_cast<std::size_t>(slant_axis)], tol);
    const int nx = static_cast<int>(bp[0].size()) - 1, ny = static_cast<int>(bp[1].size()) - 1, nz = static_cast<int>(bp[2].size()) - 1;
    // Rows: along x, at each (j, k), the runs of one region, resolved top solid first. A row asks each solid for the
    // x-intervals it covers there (exact: grid indices, and for a prism the 2D cell centre's side of its polygon).
    struct Run {
        int i0, i1;
        std::int32_t region;
        bool taken = false;
    };
    struct Span {
        int i0, i1;
    };
    // A prism's flat polygon on the grid, its index box.
    struct Flat {
        std::vector<Point> polygon;
        int u0, u1, v0, v1, t0, t1;
    };
    std::vector<Flat> flats(g.solids.size());
    std::vector<std::array<Idx, 2>> boxes(g.solids.size());
    for (std::size_t s = 0; s < g.solids.size(); ++s) {
        const LayeredSolid& v = g.solids[s];
        if (v.kind == LayeredSolid::Kind::box) {
            for (std::size_t k = 0; k < 3; ++k) {
                boxes[s][0][k] = index_of(bp[k], v.lo[k]);
                boxes[s][1][k] = index_of(bp[k], v.hi[k]);
            }
            continue;
        }
        const auto a = static_cast<std::size_t>(v.axis), ua = (a + 1) % 3, va = (a + 2) % 3;
        Flat& f = flats[s];
        f.u0 = static_cast<int>(bp[ua].size());
        f.u1 = 0;
        f.v0 = static_cast<int>(bp[va].size());
        f.v1 = 0;
        for (const Point& p : v.polygon) {
            const int iu = index_of(bp[ua], p[0]), iv = index_of(bp[va], p[1]);
            f.u0 = std::min(f.u0, iu);
            f.u1 = std::max(f.u1, iu);
            f.v0 = std::min(f.v0, iv);
            f.v1 = std::max(f.v1, iv);
            f.polygon.push_back({static_cast<double>(iu), static_cast<double>(iv), 0.0});  // in grid indices: exact
        }
        f.t0 = index_of(bp[a], v.from);
        f.t1 = index_of(bp[a], v.to);
    }
    // The x-intervals solid `s` covers in row (j, k), appended to `out`.
    const auto cover = [&](std::size_t s, int j, int k, std::vector<Span>& out) {
        const LayeredSolid& v = g.solids[s];
        if (v.kind == LayeredSolid::Kind::box) {
            const auto& b = boxes[s];
            if (j >= b[0][1] && j < b[1][1] && k >= b[0][2] && k < b[1][2] && b[0][0] < b[1][0]) out.push_back({b[0][0], b[1][0]});
            return;
        }
        const Flat& f = flats[s];
        const auto a = static_cast<std::size_t>(v.axis);
        if (a == 0) {  // (u, v) = (y, z): the row's centre in the polygon, then the prism's extent along x
            if (j < f.u0 || j >= f.u1 || k < f.v0 || k >= f.v1 || f.t0 >= f.t1) return;
            if (detail::strictly_inside_2d(f.polygon, Point{j + 0.5, k + 0.5, 0.0})) out.push_back({f.t0, f.t1});
            return;
        }
        // Along y, (u, v) = (z, x): x is v at u = k; along z, (u, v) = (x, y): x is u at v = j.
        const int along = a == 1 ? j : k;
        if (along < f.t0 || along >= f.t1) return;
        const bool x_is_v = a == 1;
        const double fixed = (x_is_v ? k : j) + 0.5;  // the other coordinate's cell centre
        if (x_is_v ? (k < f.u0 || k >= f.u1) : (j < f.v0 || j >= f.v1)) return;
        // A scanline across the rectilinear polygon: the edges crossing it, at their constant coordinate.
        std::vector<int> cuts;
        const std::size_t n = f.polygon.size();
        for (std::size_t e = 0; e < n; ++e) {
            const Point& p = f.polygon[e];
            const Point& q = f.polygon[(e + 1) % n];
            const std::size_t fc = x_is_v ? 0 : 1, xc = x_is_v ? 1 : 0;  // the fixed and the x coordinate in (u, v)
            if (p[xc] != q[xc]) continue;                                // along x: not a crossing
            if ((p[fc] < fixed) != (q[fc] < fixed)) cuts.push_back(static_cast<int>(p[xc]));
        }
        std::ranges::sort(cuts);
        for (std::size_t c = 0; c + 1 < cuts.size(); c += 2) {
            if (cuts[c] < cuts[c + 1]) out.push_back({cuts[c], cuts[c + 1]});
        }
    };
    std::vector<std::vector<Run>> rows(static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz));
    const auto row_of = [&](int j, int k) -> std::vector<Run>& { return rows[static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) + static_cast<std::size_t>(j)]; };
    std::vector<Span> spans, open, next;
    const std::vector<std::size_t> order = paint_order(g);
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            std::vector<Run>& row = row_of(j, k);
            open = {{0, nx}};  // the parts not yet painted
            for (std::size_t n = 0; n < order.size() && !open.empty(); ++n) {
                const std::size_t s = order[n];
                spans.clear();
                cover(s, j, k, spans);
                for (const Span& c : spans) {
                    next.clear();
                    for (const Span& o : open) {
                        const int a0 = std::max(o.i0, c.i0), a1 = std::min(o.i1, c.i1);
                        if (a0 >= a1) {
                            next.push_back(o);
                            continue;
                        }
                        row.push_back({a0, a1, g.solids[s].region});
                        if (o.i0 < a0) next.push_back({o.i0, a0});
                        if (a1 < o.i1) next.push_back({a1, o.i1});
                    }
                    std::swap(open, next);
                }
            }
            for (const Span& o : open) {
                if (g.background < 0) {
                    const auto mid = [&](std::size_t a, int c) { return 0.5 * (bp[a][static_cast<std::size_t>(c)] + bp[a][static_cast<std::size_t>(c + 1)]); };
                    return invalid(std::format("the layers leave part of the domain uncovered (near ({:.6f}, {:.6f}, {:.6f})) and no region "
                                               "is the background",
                                               mid(0, o.i0), mid(1, j), mid(2, k)));
                }
                row.push_back({o.i0, o.i1, g.background});
            }
            std::ranges::sort(row, {}, &Run::i0);
            std::size_t w = 0;  // runs of one region side by side joined
            for (std::size_t r = 0; r < row.size(); ++r) {
                if (w > 0 && row[w - 1].region == row[r].region && row[w - 1].i1 == row[r].i0) row[w - 1].i1 = row[r].i1;
                else row[w++] = row[r];
            }
            row.resize(w);
        }
    }

    // Blocks: a run, extended along y over rows holding the same run (same ends and region, not taken), then along z
    // over slabs of such rows.
    std::vector<Block> blocks;
    const auto match = [&](int j, int k, const Run& r) -> Run* {
        std::vector<Run>& row = row_of(j, k);
        const auto it = std::ranges::lower_bound(row, r.i0, {}, &Run::i0);
        if (it == row.end() || it->i0 != r.i0 || it->i1 != r.i1 || it->region != r.region || it->taken) return nullptr;
        return &*it;
    };
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (Run& r : row_of(j, k)) {
                if (r.taken) continue;
                int j1 = j + 1;
                while (j1 < ny && match(j1, k, r) != nullptr) ++j1;
                int k1 = k + 1;
                const auto slab = [&](int z) {
                    for (int y = j; y < j1; ++y) {
                        if (match(y, z, r) == nullptr) return false;
                    }
                    return true;
                };
                while (k1 < nz && slab(k1)) ++k1;
                for (int z = k; z < k1; ++z) {
                    for (int y = j; y < j1; ++y) {
                        if (y == j && z == k) continue;
                        match(y, z, r)->taken = true;
                    }
                }
                r.taken = true;
                blocks.push_back({{r.i0, j, k}, {r.i1, j1, k1}, r.region});
            }
        }
    }

    // Facets, plane by plane: the faces of the blocks on either side overlaid.
    struct Face {
        int u0, u1, v0, v1, block;
        bool above;  // the block is on the plane's upper side
    };
    std::map<std::pair<int, int>, std::vector<Face>> planes;
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        const Block& k = blocks[b];
        for (std::size_t a = 0; a < 3; ++a) {
            const std::size_t ua = (a + 1) % 3, va = (a + 2) % 3;
            planes[{static_cast<int>(a), k.lo[a]}].push_back({k.lo[ua], k.hi[ua], k.lo[va], k.hi[va], static_cast<int>(b), true});
            planes[{static_cast<int>(a), k.hi[a]}].push_back({k.lo[ua], k.hi[ua], k.lo[va], k.hi[va], static_cast<int>(b), false});
        }
    }
    std::vector<Facet> facets;
    std::vector<bool> used(g.patches.size(), false);
    for (const auto& [key, faces] : planes) {
        const auto [axis, at] = key;
        const auto a = static_cast<std::size_t>(axis), ua = (a + 1) % 3, va = (a + 2) % 3;
        const bool boundary = at == 0 || at == static_cast<int>(bp[a].size()) - 1;
        std::vector<int> us, vs;
        for (const Face& f : faces) {
            us.insert(us.end(), {f.u0, f.u1});
            vs.insert(vs.end(), {f.v0, f.v1});
        }
        if (boundary) {
            for (const PatchBox& p : g.patches) {
                us.insert(us.end(), {index_of(bp[ua], p.lo[ua]), index_of(bp[ua], p.hi[ua])});
                vs.insert(vs.end(), {index_of(bp[va], p.lo[va]), index_of(bp[va], p.hi[va])});
            }
        }
        for (auto* w : {&us, &vs}) {
            std::ranges::sort(*w);
            w->erase(std::unique(w->begin(), w->end()), w->end());
        }
        const int mu = static_cast<int>(us.size()) - 1, mv = static_cast<int>(vs.size()) - 1;
        const auto local = [](const std::vector<int>& w, int c) { return static_cast<int>(std::ranges::lower_bound(w, c) - w.begin()); };
        struct Side {
            int below = -1, above = -1, patch = -2;
            bool operator==(const Side&) const = default;
        };
        std::vector<Side> grid(static_cast<std::size_t>(mu) * static_cast<std::size_t>(mv));
        const auto at_grid = [&](int x, int y) -> Side& { return grid[static_cast<std::size_t>(y) * static_cast<std::size_t>(mu) + static_cast<std::size_t>(x)]; };
        for (const Face& f : faces) {
            for (int y = local(vs, f.v0); y < local(vs, f.v1); ++y) {
                for (int x = local(us, f.u0); x < local(us, f.u1); ++x) (f.above ? at_grid(x, y).above : at_grid(x, y).below) = f.block;
            }
        }
        if (boundary) {
            for (int y = 0; y < mv; ++y) {
                for (int x = 0; x < mu; ++x) {
                    Side& c = at_grid(x, y);
                    if (c.below < 0 && c.above < 0) continue;
                    c.patch = -1;
                    for (std::size_t p = 0; p < g.patches.size(); ++p) {
                        const PatchBox& b = g.patches[p];
                        const auto in = [&](std::size_t ax, double c0, double c1) { return c0 >= b.lo[ax] - tol && c1 <= b.hi[ax] + tol; };
                        const double w = bp[a][static_cast<std::size_t>(at)];
                        if (in(a, w, w) && in(ua, bp[ua][static_cast<std::size_t>(us[static_cast<std::size_t>(x)])], bp[ua][static_cast<std::size_t>(us[static_cast<std::size_t>(x + 1)])]) &&
                            in(va, bp[va][static_cast<std::size_t>(vs[static_cast<std::size_t>(y)])], bp[va][static_cast<std::size_t>(vs[static_cast<std::size_t>(y + 1)])])) {
                            c.patch = static_cast<int>(p);
                            used[p] = true;
                            break;
                        }
                    }
                }
            }
        }
        std::vector<char> done(grid.size(), 0);
        for (int y = 0; y < mv; ++y) {
            for (int x = 0; x < mu; ++x) {
                const Side s = at_grid(x, y);
                if (done[static_cast<std::size_t>(y) * static_cast<std::size_t>(mu) + static_cast<std::size_t>(x)] != 0 || (s.below < 0 && s.above < 0)) continue;
                const auto alike = [&](int xx, int yy) {
                    return done[static_cast<std::size_t>(yy) * static_cast<std::size_t>(mu) + static_cast<std::size_t>(xx)] == 0 && at_grid(xx, yy) == s;
                };
                int x1 = x + 1;
                while (x1 < mu && alike(x1, y)) ++x1;
                int y1 = y + 1;
                while (y1 < mv) {
                    bool all = true;
                    for (int xx = x; xx < x1 && all; ++xx) all = alike(xx, y1);
                    if (!all) break;
                    ++y1;
                }
                for (int yy = y; yy < y1; ++yy) {
                    for (int xx = x; xx < x1; ++xx) done[static_cast<std::size_t>(yy) * static_cast<std::size_t>(mu) + static_cast<std::size_t>(xx)] = 1;
                }
                facets.push_back({axis, at, us[static_cast<std::size_t>(x)], us[static_cast<std::size_t>(x1)], vs[static_cast<std::size_t>(y)],
                                  vs[static_cast<std::size_t>(y1)], s.below, s.above, s.patch});
            }
        }
    }
    const auto unused = std::ranges::find(used, false);
    if (unused != used.end()) {
        return invalid("a patch box selects no boundary facet (none inside it, or all taken by an earlier box)",
                       static_cast<std::size_t>(unused - used.begin()));
    }

    // Vertices: every facet corner, and on each grid line the corners along it, so an outline takes every vertex on it.
    std::map<Idx, std::size_t> vertex;
    std::map<std::tuple<int, int, int>, std::vector<int>> lines;  // (axis, the two other indices) -> positions
    const auto line_key = [](int d, const Idx& p) { return std::tuple{d, p[static_cast<std::size_t>((d + 1) % 3)], p[static_cast<std::size_t>((d + 2) % 3)]}; };
    const auto corners = [](const Facet& f) {
        const auto a = static_cast<std::size_t>(f.axis), ua = (a + 1) % 3, va = (a + 2) % 3;
        std::array<Idx, 4> c{};
        const std::array<std::pair<int, int>, 4> uv{{{f.u0, f.v0}, {f.u1, f.v0}, {f.u1, f.v1}, {f.u0, f.v1}}};
        for (std::size_t q = 0; q < 4; ++q) {
            c[q][a] = f.at;
            c[q][ua] = uv[q].first;
            c[q][va] = uv[q].second;
        }
        return c;
    };
    for (const Facet& f : facets) {
        for (const Idx& p : corners(f)) {
            if (vertex.try_emplace(p, vertex.size()).second) {
                for (int d = 0; d < 3; ++d) lines[line_key(d, p)].push_back(p[static_cast<std::size_t>(d)]);
            }
        }
    }
    for (auto& [k, v] : lines) std::ranges::sort(v);

    PolyhedralGeometry out;
    out.vertices.resize(vertex.size());
    for (const auto& [p, i] : vertex) out.vertices[i] = {bp[0][static_cast<std::size_t>(p[0])], bp[1][static_cast<std::size_t>(p[1])], bp[2][static_cast<std::size_t>(p[2])]};
    out.volumes.resize(blocks.size());
    for (std::size_t b = 0; b < blocks.size(); ++b) out.volumes[b].region = blocks[b].region;
    static constexpr std::array<std::array<const char*, 2>, 3> face_names{{{"x_min", "x_max"}, {"y_min", "y_max"}, {"z_min", "z_max"}}};
    for (const Facet& f : facets) {
        ComplexFacet c;
        const auto q = corners(f);
        for (std::size_t e = 0; e < 4; ++e) {
            const Idx& p = q[e];
            const Idx& r = q[(e + 1) % 4];
            c.vertices.push_back(vertex.at(p));
            int d = 0;
            while (p[static_cast<std::size_t>(d)] == r[static_cast<std::size_t>(d)]) ++d;  // the edge's direction
            const auto& on = lines.at(line_key(d, p));
            const int from = p[static_cast<std::size_t>(d)], to = r[static_cast<std::size_t>(d)];
            const auto lo = std::ranges::upper_bound(on, std::min(from, to)), hi = std::ranges::lower_bound(on, std::max(from, to));
            std::vector<int> between(lo, hi);
            if (from > to) std::ranges::reverse(between);
            for (const int t : between) {
                Idx m = p;
                m[static_cast<std::size_t>(d)] = t;
                c.vertices.push_back(vertex.at(m));
            }
        }
        if (f.patch >= 0) c.patch = g.patches[static_cast<std::size_t>(f.patch)].name;
        else if (f.patch == -1) c.patch = face_names[static_cast<std::size_t>(f.axis)][f.at == 0 ? 0 : 1];
        const std::size_t index = out.facets.size();
        out.facets.push_back(std::move(c));
        for (const int b : {f.below, f.above}) {
            if (b >= 0) out.volumes[static_cast<std::size_t>(b)].facets.push_back(index);
        }
    }
    return out;  // valid by construction: the tests check it (SE2), and the generator checks it again
}

}  // namespace NiTCAD::mesh
