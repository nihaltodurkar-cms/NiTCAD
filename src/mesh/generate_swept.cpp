// The swept mesh of an extrusion (Unit 54; DECISIONS.md UD1-UD4, UD8; generate.hpp).
//
// The cross-section is meshed once by the 2D generator; its nodes are repeated on every axis line
// and each prism (a triangle between two lines) is split into three tetrahedra. Node (p, j) is
// j * n2 + p, p the 2D node. The split follows the global-index rule: a quad face (a, b) x (j, j+1)
// is cut on the diagonal from min(a, b) on line j to max(a, b) on line j + 1, so the two prisms on
// either side of it cut it alike. With the triangle's 2D nodes sorted v0 < v1 < v2, the three
// tetrahedra are (v0, v1, v2 | v2'), (v0, v1 | v1', v2') and (v0 | v0', v1', v2') (primes on line
// j + 1); each uses exactly the three diagonals the rule chooses.
//
// Why the result is Delaunay, with ties only inside prisms (UD1): the six vertices of a prism lie
// on the sphere centred at (c, t_mid) with radius sqrt(R^2 + d^2 / 4), c and R the triangle's
// circumcircle, d the line spacing. Any other node (q, t) has |q - c| >= R (the 2D mesh is
// Delaunay) and |t - t_mid| >= d / 2, so it is not strictly inside. make_dual_mesh then checks
// P0-P3 and P5 as for every generated mesh; the dual is the product of the 2D and 1D duals.
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "complex_index.hpp"
#include "polygon2d.hpp"
#include "predicates.hpp"

namespace NiTCAD::mesh {

namespace {

using detail::ComplexIndex;
using detail::FacetKind;

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}

base::Error degenerate(std::string message) { return {base::ErrorCode::degenerate_mesh, std::move(message), std::nullopt}; }

using Key2 = std::pair<double, double>;

Key2 key2(const Point& p) { return {p[0], p[1]}; }

// The distance from q to the 2D segment ab.
double segment_distance_2d(const Point& a, const Point& b, const Point& q) {
    const double ux = b[0] - a[0];
    const double uy = b[1] - a[1];
    const double len2 = ux * ux + uy * uy;
    double t = ((q[0] - a[0]) * ux + (q[1] - a[1]) * uy) / len2;
    t = std::clamp(t, 0.0, 1.0);
    const double dx = a[0] + t * ux - q[0];
    const double dy = a[1] + t * uy - q[1];
    return std::sqrt(dx * dx + dy * dy);
}

// A facet's vertex set as a sorted list of exact coordinates: the key under which the extrusion's
// boundary facets are found again.
using FacetKey = std::vector<std::array<double, 3>>;

FacetKey facet_key(std::vector<Point> pts) {
    FacetKey k;
    for (const Point& p : pts) k.push_back({p[0], p[1], p[2]});
    std::ranges::sort(k);
    return k;
}

// A boundary edge of the cross-section: a cell edge on no other non-void cell.
struct BoundaryEdge {
    Point a{}, b{};
    std::size_t cell = 0;
};

}  // namespace

std::expected<GeneratedMesh, base::Error> generate_swept(const ExtrudedGeometry& x, const SweptOptions& options) {
    auto complex = extrude(x);
    if (!complex) return std::unexpected(std::move(complex.error()));
    const std::vector<double>& t = x.stations;
    const std::size_t sections = t.size() - 1;
    const std::size_t cells = x.cells.size();

    // UD4: a cell is void in every section or in none.
    std::vector<char> filled(cells, 0);
    for (std::size_t c = 0; c < cells; ++c) {
        const auto& row = x.regions[c];
        const bool any_void = std::ranges::any_of(row, [](std::int32_t r) { return r < 0; });
        const bool any_filled = std::ranges::any_of(row, [](std::int32_t r) { return r >= 0; });
        if (any_void && any_filled) {
            return std::unexpected(invalid("a swept mesh needs each cell void in every section or in none (a void next to a filled section is not supported)", c));
        }
        filled[c] = any_filled ? 1 : 0;
    }

    // Lines and straddles (UD3).
    const std::vector<double>& given = options.lines;
    if (given.size() < 2) return std::unexpected(invalid("a swept mesh needs at least two axis lines"));
    for (std::size_t i = 0; i < given.size(); ++i) {
        if (!std::isfinite(given[i])) return std::unexpected(invalid("an axis line is not finite", i));
        if (i > 0 && !(given[i] > given[i - 1])) return std::unexpected(invalid("axis lines must increase strictly", i));
    }
    if (given.front() != t.front() || given.back() != t.back()) {
        return std::unexpected(invalid("the first and last axis lines must be the first and last stations exactly"));
    }
    std::vector<double> lines = given;
    std::vector<char> straddled(t.size(), 0);
    for (std::size_t k = 0; k < options.straddles.size(); ++k) {
        const AxisStraddle& s = options.straddles[k];
        const auto it = std::ranges::find(t, s.station);
        if (it == t.end() || it == t.begin() || it == t.end() - 1) return std::unexpected(invalid("a straddle must be at an interior station", k));
        const auto i = static_cast<std::size_t>(it - t.begin());
        if (straddled[i] != 0) return std::unexpected(invalid("a station is straddled twice", k));
        if (!std::isfinite(s.width) || !(s.width > 0.0)) return std::unexpected(invalid("a straddle width must be finite and positive", k));
        const double lo = s.station - 0.5 * s.width;
        const double hi = s.station + 0.5 * s.width;
        if (!(lo > t[i - 1]) || !(hi < t[i + 1])) return std::unexpected(invalid("a straddle's lines must lie inside the neighbouring sections", k));
        if (std::ranges::any_of(given, [&](double l) { return l >= lo && l <= hi; })) {
            return std::unexpected(invalid("an axis line lies within a straddle (between or on its two lines)", k));
        }
        straddled[i] = 1;
        lines.push_back(lo);
        lines.push_back(hi);
    }
    std::ranges::sort(lines);
    for (std::size_t i = 1; i < lines.size(); ++i) {
        if (!(lines[i] > lines[i - 1])) return std::unexpected(invalid("two straddles' lines coincide or cross", i));
    }
    for (std::size_t i = 1; i + 1 < t.size(); ++i) {
        bool change = false;
        for (std::size_t c = 0; c < cells; ++c) change = change || x.regions[c][i - 1] != x.regions[c][i];
        if (change && straddled[i] == 0) return std::unexpected(invalid("a station where a cell changes region needs a straddle", i));
    }

    // The boundary facets' patches, found by their exact vertex sets in the extrusion's complex.
    const auto ax = static_cast<std::size_t>(x.axis);
    const auto ua = static_cast<std::size_t>((x.axis + 1) % 3);
    const auto va = static_cast<std::size_t>((x.axis + 2) % 3);
    const auto lift = [&](const Point& q, double tv) {
        Point p{};
        p[ax] = tv;
        p[ua] = q[0];
        p[va] = q[1];
        return p;
    };
    auto index = ComplexIndex::create(*complex);
    if (!index) return std::unexpected(std::move(index.error()));
    const ComplexIndex& cx = *index;
    const PolyhedralGeometry& g = cx.geometry();
    const std::vector<std::string>& patch_names = cx.report().patches;
    std::map<FacetKey, int> facet_patch_of;
    for (std::size_t f = 0; f < g.facets.size(); ++f) {
        if (cx.facets()[f].kind != FacetKind::boundary) continue;
        std::vector<Point> pts;
        for (const std::size_t v : g.facets[f].vertices) pts.push_back(g.vertices[v]);
        facet_patch_of.emplace(facet_key(std::move(pts)), cx.facets()[f].patch);
    }
    std::vector<std::vector<Point>> poly(cells);
    for (std::size_t c = 0; c < cells; ++c) {
        for (const Point& p : x.cells[c].polygon) poly[c].push_back({p[0], p[1], 0.0});
    }
    // End patches per cell (station 0 and the last), lateral boundary edges with a patch per section.
    std::vector<std::array<int, 2>> end_patch(cells, {-1, -1});
    std::map<std::pair<Key2, Key2>, int> edge_use;  // undirected cell edge -> number of filled cells on it
    for (std::size_t c = 0; c < cells; ++c) {
        if (filled[c] == 0) continue;
        for (std::size_t i = 0; i < poly[c].size(); ++i) {
            auto ka = key2(poly[c][i]);
            auto kb = key2(poly[c][(i + 1) % poly[c].size()]);
            if (kb < ka) std::swap(ka, kb);
            ++edge_use[{ka, kb}];
        }
    }
    std::vector<BoundaryEdge> bedges;
    std::vector<std::vector<int>> bedge_patch;  // [boundary edge][section]
    for (std::size_t c = 0; c < cells; ++c) {
        if (filled[c] == 0) continue;
        for (const std::size_t s : {std::size_t{0}, t.size() - 1}) {
            std::vector<Point> pts;
            for (const Point& q : poly[c]) pts.push_back(lift(q, t[s]));
            const auto it = facet_patch_of.find(facet_key(std::move(pts)));
            if (it == facet_patch_of.end()) return std::unexpected(degenerate("an end facet of the extrusion was not found"));
            end_patch[c][s == 0 ? 0 : 1] = it->second;
        }
        for (std::size_t i = 0; i < poly[c].size(); ++i) {
            const Point& a = poly[c][i];
            const Point& b = poly[c][(i + 1) % poly[c].size()];
            auto ka = key2(a);
            auto kb = key2(b);
            if (kb < ka) std::swap(ka, kb);
            if (edge_use[{ka, kb}] != 1) continue;
            std::vector<int> row(sections, -1);
            for (std::size_t s = 0; s < sections; ++s) {
                const auto it = facet_patch_of.find(facet_key({lift(a, t[s]), lift(b, t[s]), lift(b, t[s + 1]), lift(a, t[s + 1])}));
                if (it == facet_patch_of.end()) return std::unexpected(degenerate("a lateral facet of the extrusion was not found"));
                row[s] = it->second;
            }
            bedges.push_back({a, b, c});
            bedge_patch.push_back(std::move(row));
        }
    }
    // A station where a lateral patch changes is straddled or is a line.
    for (std::size_t i = 1; i + 1 < t.size(); ++i) {
        if (straddled[i] != 0 || std::ranges::find(given, t[i]) != given.end()) continue;
        for (std::size_t e = 0; e < bedges.size(); ++e) {
            if (bedge_patch[e][i - 1] != bedge_patch[e][i]) {
                return std::unexpected(invalid("a station where a lateral boundary patch changes must be straddled or be an axis line", i));
            }
        }
    }

    // The cross-section's 2D regions (UD2): cells alike in regions, end patches and lateral patch
    // sets share one.
    struct Signature {
        std::vector<std::int32_t> regions;
        std::array<int, 2> ends{};
        std::vector<std::set<int>> lateral;
        auto operator<=>(const Signature&) const = default;
    };
    std::vector<Signature> sigs;
    std::vector<int> column(cells, -1);
    {
        std::vector<std::vector<std::set<int>>> lateral(cells, std::vector<std::set<int>>(sections));
        for (std::size_t e = 0; e < bedges.size(); ++e) {
            for (std::size_t s = 0; s < sections; ++s) lateral[bedges[e].cell][s].insert(bedge_patch[e][s]);
        }
        std::map<Signature, int> seen;
        for (std::size_t c = 0; c < cells; ++c) {
            if (filled[c] == 0) continue;
            Signature sg{x.regions[c], end_patch[c], lateral[c]};
            const auto [it, fresh] = seen.emplace(sg, static_cast<int>(sigs.size()));
            if (fresh) sigs.push_back(sg);
            column[c] = it->second;
        }
    }
    // A lateral patch change at a cross-section vertex only where the boundary turns or two 2D
    // regions meet (otherwise the 2D mesh has no node there to carry it).
    for (std::size_t e = 0; e < bedges.size(); ++e) {
        for (std::size_t f = 0; f < bedges.size(); ++f) {
            if (f == e || bedges[e].b != bedges[f].a) continue;
            const bool turns = predicates::orient2d(bedges[e].a, bedges[e].b, bedges[f].b) != 0;
            const bool meet = column[bedges[e].cell] != column[bedges[f].cell];
            if (turns || meet) continue;
            for (std::size_t s = 0; s < sections; ++s) {
                if (bedge_patch[e][s] != bedge_patch[f][s]) {
                    return std::unexpected(invalid("a lateral patch changes at a cross-section vertex where the boundary neither turns nor meets another 2D region", bedges[e].cell));
                }
            }
        }
    }

    bool has_interfaces = false;
    for (const auto& fc : cx.facets()) has_interfaces = has_interfaces || fc.kind == FacetKind::interface;
    const double s_res = options.section.interface_resolution;
    if (has_interfaces && !(std::isfinite(s_res) && s_res > 0.0)) {
        return std::unexpected(invalid("a swept mesh with interfaces needs a positive interface_resolution"));
    }

    // The 2D mesh.
    Geometry2D g2;
    for (std::size_t c = 0; c < cells; ++c) {
        if (filled[c] != 0) g2.polygons.push_back({poly[c], column[c]});
    }
    MeshOptions2D o2 = options.section;
    o2.region_max_size.assign(sigs.size(), 0.0);
    for (std::size_t k = 0; k < sigs.size(); ++k) {
        double h = 0.0;
        for (const std::int32_t r : sigs[k].regions) {
            const auto ri = static_cast<std::size_t>(r);
            const double hr = ri < options.section.region_max_size.size() ? options.section.region_max_size[ri] : 0.0;
            if (hr > 0.0) h = h > 0.0 ? std::min(h, hr) : hr;
        }
        o2.region_max_size[k] = h;
    }
    auto m2 = generate_mesh(g2, o2);
    if (!m2) {
        auto e = std::move(m2.error());
        e.message = "the cross-section: " + e.message;
        return std::unexpected(std::move(e));
    }
    const SimplexMesh& s2 = m2->simplices;
    const std::size_t n2 = s2.points.size();
    const std::size_t nl = lines.size();
    if (n2 > options.max_nodes / nl) return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "a swept mesh would exceed max_nodes", std::nullopt});

    // Nodes, tetrahedra.
    SimplexMesh sm;
    sm.dimension = 3;
    sm.patch_names = patch_names;
    sm.points.reserve(n2 * nl);
    for (std::size_t j = 0; j < nl; ++j) {
        for (std::size_t p = 0; p < n2; ++p) sm.points.push_back(lift(s2.points[p], lines[j]));
    }
    const auto id = [&](std::size_t p, std::size_t j) { return static_cast<NodeId>(j * n2 + p); };
    const std::size_t tris = s2.cells.size() / 3;
    sm.cells.reserve(tris * (nl - 1) * 12);
    for (std::size_t c = 0; c < tris; ++c) {
        std::array<std::size_t, 3> v{static_cast<std::size_t>(s2.cells[3 * c]), static_cast<std::size_t>(s2.cells[3 * c + 1]),
                                     static_cast<std::size_t>(s2.cells[3 * c + 2])};
        std::ranges::sort(v);
        for (std::size_t j = 0; j + 1 < nl; ++j) {
            for (const auto& tet : {std::array{id(v[0], j), id(v[1], j), id(v[2], j), id(v[2], j + 1)},
                                    std::array{id(v[0], j), id(v[1], j), id(v[1], j + 1), id(v[2], j + 1)},
                                    std::array{id(v[0], j), id(v[0], j + 1), id(v[1], j + 1), id(v[2], j + 1)}}) {
                sm.cells.insert(sm.cells.end(), tet.begin(), tet.end());
            }
        }
    }

    // Sections of lines, and the 2D node columns.
    const auto section_of = [&](double tv) {
        std::size_t i = 0;
        while (i + 1 < sections && tv >= t[i + 1]) ++i;
        return i;
    };
    std::vector<std::size_t> line_section(nl);
    for (std::size_t j = 0; j < nl; ++j) line_section[j] = section_of(lines[j]);
    const std::vector<std::int32_t>& col2 = m2->node_region;

    // Boundary triangles: both end stations, and each 2D boundary edge swept.
    double extent = 0.0;
    for (const Point& p : g.vertices) {
        for (std::size_t c = 0; c < 3; ++c) extent = std::max(extent, std::abs(p[c]));
    }
    const double on_tol = 1e-10 * extent;
    const auto lateral_candidates = [&](const Point& q) {
        std::vector<std::size_t> out;
        for (std::size_t e = 0; e < bedges.size(); ++e) {
            if (segment_distance_2d(bedges[e].a, bedges[e].b, q) <= on_tol) out.push_back(e);
        }
        return out;
    };
    const auto cell_at = [&](const Point& q) -> int {
        for (std::size_t c = 0; c < cells; ++c) {
            if (filled[c] != 0 && detail::inside_closed_2d(poly[c], q)) return static_cast<int>(c);
        }
        return -1;
    };
    std::vector<std::int32_t> corner_patch;
    bool corners_differ = false;
    // End stations: the triangle's patch from the cell holding its centroid; each corner's from its
    // 2D region (its column's end patch).
    for (const std::size_t j : {std::size_t{0}, nl - 1}) {
        const std::size_t end = j == 0 ? 0 : 1;
        for (std::size_t c = 0; c < tris; ++c) {
            const std::array<std::size_t, 3> v{static_cast<std::size_t>(s2.cells[3 * c]), static_cast<std::size_t>(s2.cells[3 * c + 1]),
                                               static_cast<std::size_t>(s2.cells[3 * c + 2])};
            Point cen{0, 0, 0};
            for (const std::size_t p : v) cen = {cen[0] + s2.points[p][0] / 3.0, cen[1] + s2.points[p][1] / 3.0, 0.0};
            const int cell = cell_at(cen);
            if (cell < 0) return std::unexpected(degenerate("a cross-section triangle's centroid lies in no cell"));
            const int patch = end_patch[static_cast<std::size_t>(cell)][end];
            for (const std::size_t p : v) sm.facets.push_back(id(p, j));
            sm.facet_patch.push_back(patch);
            for (const std::size_t p : v) {
                const int cp = sigs[static_cast<std::size_t>(col2[p])].ends[end];
                corners_differ = corners_differ || cp != patch;
                corner_patch.push_back(cp);
            }
        }
    }
    // The 2D boundary edges: edges of exactly one triangle.
    std::map<std::pair<std::size_t, std::size_t>, int> e2;
    for (std::size_t c = 0; c < tris; ++c) {
        for (std::size_t k = 0; k < 3; ++k) {
            auto a = static_cast<std::size_t>(s2.cells[3 * c + k]);
            auto b = static_cast<std::size_t>(s2.cells[3 * c + (k + 1) % 3]);
            if (b < a) std::swap(a, b);
            ++e2[{a, b}];
        }
    }
    const auto node_lateral_patch = [&](std::size_t p, std::size_t j, int fallback) {
        // The patches of every lateral facet holding the node; one, or the triangle's.
        const auto cand = lateral_candidates(s2.points[p]);
        std::set<int> found;
        const double tv = lines[j];
        for (std::size_t s = 0; s < sections; ++s) {
            if (tv < t[s] || tv > t[s + 1]) continue;
            for (const std::size_t e : cand) found.insert(bedge_patch[e][s]);
        }
        return found.size() == 1 ? *found.begin() : fallback;
    };
    for (const auto& [edge, uses] : e2) {
        if (uses != 1) continue;
        const auto [a, b] = edge;  // a < b
        for (std::size_t j = 0; j + 1 < nl; ++j) {
            // The quad (a, b) x (j, j + 1), cut from a on j to b on j + 1 (a < b).
            for (const auto& tri : {std::array{id(a, j), id(b, j), id(b, j + 1)}, std::array{id(a, j), id(b, j + 1), id(a, j + 1)}}) {
                // The triangle's patch: the lateral facet holding its centroid.
                double tc = 0.0;
                Point qc{0, 0, 0};
                for (const NodeId n : tri) {
                    const std::size_t p = static_cast<std::size_t>(n) % n2;
                    const std::size_t l = static_cast<std::size_t>(n) / n2;
                    tc += lines[l] / 3.0;
                    qc = {qc[0] + s2.points[p][0] / 3.0, qc[1] + s2.points[p][1] / 3.0, 0.0};
                }
                const auto cand = lateral_candidates(qc);
                if (cand.empty()) return std::unexpected(degenerate("a boundary triangle of the swept mesh lies on no lateral facet"));
                const int patch = bedge_patch[cand.front()][section_of(tc)];
                for (const NodeId n : tri) sm.facets.push_back(n);
                sm.facet_patch.push_back(patch);
                for (const NodeId n : tri) {
                    const int cp = node_lateral_patch(static_cast<std::size_t>(n) % n2, static_cast<std::size_t>(n) / n2, patch);
                    corners_differ = corners_differ || cp != patch;
                    corner_patch.push_back(cp);
                }
            }
        }
    }
    if (corners_differ) sm.facet_corner_patch = std::move(corner_patch);

    DualMeshReport dual_report;
    auto dual = make_dual_mesh(sm, &dual_report);
    if (!dual) {
        return std::unexpected(base::Error{base::ErrorCode::degenerate_mesh, "the swept mesh failed an import check: " + dual.error().message,
                                           dual.error().context});
    }
    double volume = 0.0;
    for (const double v : cx.report().volume) volume += v;
    double dual_volume = 0.0;
    for (const double v : dual->volumes()) dual_volume += v;
    if (!(std::abs(dual_volume - volume) <= 1e-12 * volume)) return std::unexpected(degenerate("the swept mesh does not fill the domain"));

    // Regions: the column's region in the line's section.
    std::vector<std::int32_t> node_region(sm.points.size());
    for (std::size_t j = 0; j < nl; ++j) {
        for (std::size_t p = 0; p < n2; ++p) {
            node_region[j * n2 + p] = sigs[static_cast<std::size_t>(col2[p])].regions[line_section[j]];
        }
    }

    MeshReport report;
    report.nodes = dual->node_count();
    report.triangles = tris;
    report.tetrahedra = sm.cells.size() / 4;
    report.interface_nodes = m2->report.interface_nodes * nl;
    report.exempt_triangles = m2->report.exempt_triangles;
    report.min_angle_deg = m2->report.min_angle_deg;

    std::optional<InterfaceCertificate> certificate;
    // Certified with interfaces, and wherever corners carry patches (FD10: such a mesh is certified
    // with P6 or not at all); those corners come from straddled 2D junctions, so s is given.
    if (has_interfaces || !sm.facet_corner_patch.empty()) {
        InterfaceSet set{{}, s_res};
        for (std::size_t f = 0; f < g.facets.size(); ++f) {
            const auto& fc = cx.facets()[f];
            if (fc.kind != FacetKind::interface) continue;
            InterfacePiece piece;
            for (const std::size_t v : g.facets[f].vertices) piece.vertices.push_back(g.vertices[v]);
            piece.region_a = g.volumes[static_cast<std::size_t>(fc.volumes[0])].region;
            piece.region_b = g.volumes[static_cast<std::size_t>(fc.volumes[1])].region;
            set.pieces.push_back(std::move(piece));
        }
        std::optional<PatchBoundarySet> patch_set;
        if (!sm.facet_corner_patch.empty()) {
            // P6: every edge between two boundary facets of different patches; zones at every end of
            // one (conservative: the product is exact up to them, UD8).
            PatchBoundarySet& ps = patch_set.emplace();
            ps.resolution = s_res;
            std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> edge_facets;
            for (std::size_t f = 0; f < g.facets.size(); ++f) {
                if (cx.facets()[f].kind != FacetKind::boundary) continue;
                const auto& vs = g.facets[f].vertices;
                for (std::size_t i = 0; i < vs.size(); ++i) {
                    auto a = vs[i];
                    auto b = vs[(i + 1) % vs.size()];
                    if (b < a) std::swap(a, b);
                    edge_facets[{a, b}].push_back(f);
                }
            }
            std::set<std::size_t> zone_vertices;
            for (const auto& [edge, fs] : edge_facets) {
                if (fs.size() != 2) continue;
                const int pa = cx.facets()[fs[0]].patch;
                const int pb = cx.facets()[fs[1]].patch;
                if (pa == pb) continue;
                ps.boundaries.push_back({g.vertices[edge.first], g.vertices[edge.second], pa, pb});
                zone_vertices.insert(edge.first);
                zone_vertices.insert(edge.second);
            }
            for (const std::size_t v : zone_vertices) ps.zone_centres.push_back(g.vertices[v]);
            ps.patch_area.assign(patch_names.size(), 0.0);
            for (std::size_t f = 0; f < g.facets.size(); ++f) {
                const auto& fc = cx.facets()[f];
                if (fc.kind != FacetKind::boundary) continue;
                for (const auto& tr : fc.triangles) {
                    const Point& A = g.vertices[tr[0]];
                    const Point& B = g.vertices[tr[1]];
                    const Point& C = g.vertices[tr[2]];
                    const Point u{B[0] - A[0], B[1] - A[1], B[2] - A[2]};
                    const Point w{C[0] - A[0], C[1] - A[1], C[2] - A[2]};
                    const Point n{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]};
                    ps.patch_area[static_cast<std::size_t>(fc.patch)] += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                }
            }
        }
        auto cert = certify_interfaces(sm, *dual, node_region, set, patch_set ? &*patch_set : nullptr);
        if (!cert) {
            return std::unexpected(base::Error{base::ErrorCode::degenerate_mesh,
                                               "the swept mesh did not meet the interface certificate (P4, P6): " + cert.error().message,
                                               cert.error().context});
        }
        certificate = std::move(*cert);
    }
    return GeneratedMesh{std::move(*dual), std::move(sm), std::move(node_region), std::move(certificate), report,
                         generation_digest(x, options)};
}

std::uint64_t generation_digest(const ExtrudedGeometry& x, const SweptOptions& options) {
    // FNV-1a 64 over little-endian 8-byte words, as the complex's digest; text byte by byte.
    std::uint64_t h = 0xcbf29ce484222325ULL;
    const auto word = [&](std::uint64_t v) {
        for (int i = 0; i < 8; ++i) {
            h ^= (v >> (8 * i)) & 0xffU;
            h *= 0x100000001b3ULL;
        }
    };
    const auto real = [&](double d) { word(std::bit_cast<std::uint64_t>(d)); };
    const auto point = [&](const Point& p) {
        for (const double c : p) real(c);
    };
    const auto text = [&](const std::string& s) {
        word(s.size());
        for (const char c : s) {
            h ^= static_cast<std::uint64_t>(static_cast<unsigned char>(c));
            h *= 0x100000001b3ULL;
        }
    };
    const auto reals = [&](const std::vector<double>& v) {
        word(v.size());
        for (const double d : v) real(d);
    };
    word(0x5357455054ULL);  // "SWEPT": distinct from the complex's digest
    word(static_cast<std::uint64_t>(x.axis));
    reals(x.stations);
    word(x.cells.size());
    for (const ExtrusionCell& c : x.cells) {
        text(c.name);
        word(c.polygon.size());
        for (const Point& p : c.polygon) point(p);
    }
    word(x.regions.size());
    for (const auto& row : x.regions) {
        word(row.size());
        for (const std::int32_t r : row) word(static_cast<std::uint64_t>(static_cast<std::int64_t>(r)));
    }
    word(x.patches.size());
    for (const PatchBox& b : x.patches) {
        text(b.name);
        point(b.lo);
        point(b.hi);
    }
    text(x.default_patch);
    const MeshOptions2D& s = options.section;
    real(s.max_size);
    reals(s.region_max_size);
    word(s.boxes.size());
    for (const SizeBox& b : s.boxes) {
        point(b.lo);
        point(b.hi);
        real(b.max_size);
    }
    word(s.layers.size());
    for (const BoundaryLayer& l : s.layers) {
        point(l.a);
        point(l.b);
        reals(l.along);
        reals(l.offsets);
    }
    real(s.interface_resolution);
    real(s.min_angle_deg);
    word(s.max_nodes);
    reals(options.lines);
    word(options.straddles.size());
    for (const AxisStraddle& st : options.straddles) {
        real(st.station);
        real(st.width);
    }
    word(options.max_nodes);
    return h;
}

}  // namespace NiTCAD::mesh
