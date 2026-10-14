#include "NiTCAD/mesh/simplex_mesh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "dual_geometry.hpp"
#include "dual_tolerance.hpp"
#include "predicates.hpp"

namespace NiTCAD::mesh {

namespace {

// P5: closure tolerance relative to the sum of magnitudes. The M2 rule is in dual_tolerance.hpp.
constexpr double closure_tolerance = 1e-11;
constexpr NodeId none = -1;

using Cell = std::array<NodeId, 4>;  // dimension + 1 used; the rest `none`
using Key = std::array<NodeId, 3>;   // a facet's nodes, increasing; the rest `none`

base::Error error(base::ErrorCode code, std::string message, std::optional<std::size_t> index,
                  std::optional<double> value = std::nullopt) {
    return {code, std::move(message), base::ErrorContext{.index = index, .value = value}};
}

std::string pair_text(NodeId a, NodeId b) {
    return "(" + std::to_string(a) + ", " + std::to_string(b) + ")";
}

// Exact orientation of up to four points of a cell, in the mesh dimension.
int orientation(int dim, const std::vector<Point>& pts, const Cell& c) {
    const auto at = [&](std::size_t k) -> const Point& {
        return pts[static_cast<std::size_t>(c[k])];
    };
    return dim == 2 ? predicates::orient2d(at(0), at(1), at(2))
                    : predicates::orient3d(at(0), at(1), at(2), at(3));
}

struct FacetRef {
    Key key;
    std::size_t cell;
    std::size_t local;  // the cell's vertex opposite the facet
};

// A uniform bucket grid over the points, for the Gabriel check's neighbour queries.
class PointGrid {
public:
    PointGrid(int dim, const std::vector<Point>& pts) : dim_(static_cast<std::size_t>(dim)) {
        const std::size_t n = pts.size();
        const auto per_axis = static_cast<std::size_t>(
            std::max(1.0, std::ceil(std::pow(static_cast<double>(n), 1.0 / static_cast<double>(dim)))));
        for (std::size_t a = 0; a < 3; ++a) {
            lo_[a] = hi_[a] = 0.0;
            cells_[a] = 1;
        }
        for (std::size_t a = 0; a < dim_; ++a) {
            lo_[a] = hi_[a] = pts[0][a];
            for (const Point& p : pts) {
                lo_[a] = std::min(lo_[a], p[a]);
                hi_[a] = std::max(hi_[a], p[a]);
            }
            cells_[a] = hi_[a] > lo_[a] ? per_axis : 1;
        }
        const std::size_t total = cells_[0] * cells_[1] * cells_[2];
        start_.assign(total + 1, 0);
        std::vector<std::size_t> bucket(n);
        for (std::size_t i = 0; i < n; ++i) {
            bucket[i] = flat(index_of(pts[i]));
            ++start_[bucket[i] + 1];
        }
        for (std::size_t b = 0; b < total; ++b) start_[b + 1] += start_[b];
        nodes_.resize(n);
        std::vector<std::size_t> fill(start_.begin(), start_.end() - 1);
        for (std::size_t i = 0; i < n; ++i) nodes_[fill[bucket[i]]++] = static_cast<NodeId>(i);
    }

    // Every node in the buckets that meet the box [lo, hi], in increasing bucket order.
    template <class Visit>
    void visit(const Point& lo, const Point& hi, Visit visit_node) const {
        const auto a = index_of(lo);
        const auto b = index_of(hi);
        for (std::size_t z = a[2]; z <= b[2]; ++z) {
            for (std::size_t y = a[1]; y <= b[1]; ++y) {
                for (std::size_t x = a[0]; x <= b[0]; ++x) {
                    const std::size_t f = flat({x, y, z});
                    for (std::size_t k = start_[f]; k < start_[f + 1]; ++k) visit_node(nodes_[k]);
                }
            }
        }
    }

private:
    std::array<std::size_t, 3> index_of(const Point& p) const {
        std::array<std::size_t, 3> r{0, 0, 0};
        for (std::size_t a = 0; a < dim_; ++a) {
            if (cells_[a] == 1) continue;
            const double t = (p[a] - lo_[a]) / (hi_[a] - lo_[a]) * static_cast<double>(cells_[a]);
            const double clamped = std::clamp(t, 0.0, static_cast<double>(cells_[a] - 1));
            r[a] = static_cast<std::size_t>(clamped);
        }
        return r;
    }
    std::size_t flat(const std::array<std::size_t, 3>& i) const {
        return (i[2] * cells_[1] + i[1]) * cells_[0] + i[0];
    }

    std::size_t dim_;
    std::array<double, 3> lo_{};
    std::array<double, 3> hi_{};
    std::array<std::size_t, 3> cells_{};
    std::vector<std::size_t> start_;
    std::vector<NodeId> nodes_;
};

struct Piece {
    NodeId a;  // edge: first < second; patch area: a = patch, b = node
    NodeId b;
    double value;
};

}  // namespace

std::expected<Mesh, base::Error> make_dual_mesh(const SimplexMesh& s, DualMeshReport* report) {
    using base::ErrorCode;
    const int dim = s.dimension;
    if (dim != 2 && dim != 3) {
        return std::unexpected(
            error(ErrorCode::invalid_input, "simplex mesh dimension must be 2 or 3", std::nullopt));
    }
    const auto D = static_cast<std::size_t>(dim);
    const std::size_t per_cell = D + 1;
    if (s.points.empty() || s.cells.empty()) {
        return std::unexpected(error(ErrorCode::invalid_input,
                                     "simplex mesh needs points and cells", std::nullopt));
    }
    if (s.points.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
        return std::unexpected(error(ErrorCode::resource_exhausted,
                                     "too many simplex mesh points for NodeId", std::nullopt));
    }
    if (s.cells.size() % per_cell != 0 || s.facets.size() % D != 0 ||
        s.facet_patch.size() != s.facets.size() / D) {
        return std::unexpected(error(ErrorCode::invalid_input,
                                     "cells, facets or facet_patch sized inconsistently with the "
                                     "dimension",
                                     std::nullopt));
    }
    const std::vector<Point>& pts = s.points;
    const std::size_t n = pts.size();
    const auto node_count = static_cast<NodeId>(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t c = 0; c < 3; ++c) {
            const double v = pts[i][c];
            if (!std::isfinite(v) || (c >= D && v != 0.0)) {
                return std::unexpected(error(ErrorCode::invalid_input,
                                             "node coordinate is not finite, or nonzero beyond "
                                             "the dimension",
                                             i, v));
            }
        }
    }

    // Cells: range, repeats, duplicates, use of every node, orientation (P0).
    const std::size_t cell_count = s.cells.size() / per_cell;
    std::vector<Cell> cells(cell_count, Cell{none, none, none, none});
    for (std::size_t c = 0; c < cell_count; ++c) {
        for (std::size_t k = 0; k < per_cell; ++k) {
            const NodeId v = s.cells[c * per_cell + k];
            if (v < 0 || v >= node_count) {
                return std::unexpected(
                    error(ErrorCode::invalid_input, "cell node out of range", c));
            }
            for (std::size_t q = 0; q < k; ++q) {
                if (cells[c][q] == v) {
                    return std::unexpected(
                        error(ErrorCode::invalid_input, "cell node repeated", c));
                }
            }
            cells[c][k] = v;
        }
    }
    {
        std::vector<std::pair<Cell, std::size_t>> sorted(cell_count);
        for (std::size_t c = 0; c < cell_count; ++c) {
            Cell key = cells[c];
            std::sort(key.begin(), key.begin() + static_cast<std::ptrdiff_t>(per_cell));
            sorted[c] = {key, c};
        }
        std::ranges::sort(sorted);
        for (std::size_t k = 1; k < sorted.size(); ++k) {
            if (sorted[k].first == sorted[k - 1].first) {
                return std::unexpected(
                    error(ErrorCode::invalid_input, "cell listed twice", sorted[k].second));
            }
        }
    }
    {
        std::vector<bool> used(n, false);
        for (const Cell& c : cells) {
            for (std::size_t k = 0; k < per_cell; ++k) used[static_cast<std::size_t>(c[k])] = true;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!used[i]) {
                return std::unexpected(error(ErrorCode::invalid_input, "node is in no cell", i));
            }
        }
    }
    for (std::size_t c = 0; c < cell_count; ++c) {
        const int o = orientation(dim, pts, cells[c]);
        if (o == 0) {
            return std::unexpected(
                error(ErrorCode::degenerate_mesh, "flat cell (zero orientation, P0)", c));
        }
        if (o < 0) std::swap(cells[c][0], cells[c][1]);
    }

    // Facets: manifold, no fold, locally Delaunay (P1); the boundary facets.
    std::vector<FacetRef> refs;
    refs.reserve(cell_count * per_cell);
    for (std::size_t c = 0; c < cell_count; ++c) {
        for (std::size_t l = 0; l < per_cell; ++l) {
            Key key{none, none, none};
            for (std::size_t q = 0, k = 0; q < per_cell; ++q) {
                if (q != l) key[k++] = cells[c][q];
            }
            std::sort(key.begin(), key.begin() + static_cast<std::ptrdiff_t>(D));
            refs.push_back({key, c, l});
        }
    }
    std::ranges::sort(refs, [](const FacetRef& x, const FacetRef& y) {
        return std::tie(x.key, x.cell, x.local) < std::tie(y.key, y.cell, y.local);
    });
    std::vector<FacetRef> boundary;
    for (std::size_t g = 0; g < refs.size();) {
        std::size_t e = g + 1;
        while (e < refs.size() && refs[e].key == refs[g].key) ++e;
        if (e - g == 1) {
            boundary.push_back(refs[g]);
        } else if (e - g > 2) {
            return std::unexpected(error(ErrorCode::degenerate_mesh,
                                         "facet shared by more than two cells", refs[g + 2].cell));
        } else {
            const FacetRef& r1 = refs[g];
            const FacetRef& r2 = refs[g + 1];
            const NodeId opposite = cells[r2.cell][r2.local];
            // The first cell with its own opposite vertex replaced by the second cell's: positive
            // orientation would put both cells on the same side of the shared facet.
            Cell swapped = cells[r1.cell];
            swapped[r1.local] = opposite;
            if (orientation(dim, pts, swapped) >= 0) {
                return std::unexpected(error(ErrorCode::degenerate_mesh,
                                             "two cells on the same side of their shared facet",
                                             r2.cell));
            }
            const Cell& c1 = cells[r1.cell];
            const auto P = [&](NodeId v) -> const Point& { return pts[static_cast<std::size_t>(v)]; };
            const int inside = dim == 2
                                   ? predicates::incircle(P(c1[0]), P(c1[1]), P(c1[2]), P(opposite))
                                   : predicates::insphere(P(c1[0]), P(c1[1]), P(c1[2]), P(c1[3]),
                                                          P(opposite));
            if (inside > 0) {
                return std::unexpected(error(ErrorCode::degenerate_mesh,
                                             "not locally Delaunay (P1): node " +
                                                 std::to_string(opposite) +
                                                 " is strictly inside the circumsphere of cell " +
                                                 std::to_string(r1.cell),
                                             r1.cell));
            }
        }
        g = e;
    }

    // Patches and listed facets.
    {
        std::set<std::string, std::less<>> names;
        for (std::size_t p = 0; p < s.patch_names.size(); ++p) {
            if (s.patch_names[p].empty() || !names.insert(s.patch_names[p]).second) {
                return std::unexpected(error(ErrorCode::invalid_input,
                                             "boundary patch name is empty or repeated", p));
            }
        }
    }
    const std::size_t listed_count = s.facets.size() / D;
    // Per topological boundary facet (in `boundary` order), its patch or -1.
    std::vector<std::int32_t> boundary_patch(boundary.size(), -1);
    {
        std::vector<std::pair<Key, std::size_t>> listed(listed_count);
        for (std::size_t f = 0; f < listed_count; ++f) {
            Key key{none, none, none};
            for (std::size_t k = 0; k < D; ++k) {
                const NodeId v = s.facets[f * D + k];
                if (v < 0 || v >= node_count) {
                    return std::unexpected(
                        error(ErrorCode::invalid_input, "facet node out of range", f));
                }
                key[k] = v;
            }
            std::sort(key.begin(), key.begin() + static_cast<std::ptrdiff_t>(D));
            for (std::size_t k = 1; k < D; ++k) {
                if (key[k] == key[k - 1]) {
                    return std::unexpected(
                        error(ErrorCode::invalid_input, "facet node repeated", f));
                }
            }
            const std::int32_t patch = s.facet_patch[f];
            if (patch < 0 || static_cast<std::size_t>(patch) >= s.patch_names.size()) {
                return std::unexpected(
                    error(ErrorCode::invalid_input, "facet patch index out of range", f));
            }
            listed[f] = {key, f};
        }
        std::ranges::sort(listed);
        for (std::size_t k = 1; k < listed.size(); ++k) {
            if (listed[k].first == listed[k - 1].first) {
                return std::unexpected(
                    error(ErrorCode::invalid_input, "facet listed twice", listed[k].second));
            }
        }
        // `boundary` is sorted by key (it is a subsequence of the sorted refs).
        for (const auto& [key, f] : listed) {
            const auto it = std::ranges::lower_bound(boundary, key, {}, &FacetRef::key);
            if (it == boundary.end() || it->key != key) {
                return std::unexpected(error(ErrorCode::invalid_input,
                                             "listed facet is not a boundary facet (a face of "
                                             "exactly one cell)",
                                             f));
            }
            boundary_patch[static_cast<std::size_t>(it - boundary.begin())] = s.facet_patch[f];
        }
    }

    // P2: every boundary facet is Gabriel.
    {
        const PointGrid grid(dim, pts);
        for (const FacetRef& b : boundary) {
            const auto P = [&](NodeId v) -> const Point& { return pts[static_cast<std::size_t>(v)]; };
            Point centre{0.0, 0.0, 0.0};
            if (dim == 2) {
                for (std::size_t a = 0; a < 2; ++a) centre[a] = 0.5 * (P(b.key[0])[a] + P(b.key[1])[a]);
            } else {
                centre = detail::triangle_circumcentre({P(b.key[0]), P(b.key[1]), P(b.key[2])});
            }
            // The rounded centre is within a few ulps; widen the box well beyond that.
            const double r = detail::distance(centre, P(b.key[0])) * (1.0 + 1e-6);
            Point lo = centre;
            Point hi = centre;
            for (std::size_t a = 0; a < D; ++a) {
                lo[a] -= r;
                hi[a] += r;
            }
            std::optional<NodeId> inside;
            grid.visit(lo, hi, [&](NodeId q) {
                if (inside || q == b.key[0] || q == b.key[1] || (dim == 3 && q == b.key[2])) return;
                const int t = dim == 2 ? predicates::in_diametral_circle(P(b.key[0]), P(b.key[1]), P(q))
                                       : predicates::in_diametral_sphere(P(b.key[0]), P(b.key[1]),
                                                                         P(b.key[2]), P(q));
                if (t > 0) inside = q;
            });
            if (inside) {
                return std::unexpected(error(ErrorCode::degenerate_mesh,
                                             "boundary facet is not Gabriel (P2): node " +
                                                 std::to_string(*inside) +
                                                 " is strictly inside its diametral sphere",
                                             b.cell, static_cast<double>(*inside)));
            }
        }
    }

    // The circumcentric pieces (M1).
    std::vector<Piece> edge_pieces;
    edge_pieces.reserve(cell_count * (dim == 2 ? 3 : 6));
    std::vector<double> volume(n, 0.0);
    std::vector<double> volume_magnitude(n, 0.0);
    for (std::size_t c = 0; c < cell_count; ++c) {
        const Cell& cell = cells[c];
        const auto P = [&](std::size_t k) { return pts[static_cast<std::size_t>(cell[k])]; };
        const auto add = [&](int i, int j, double piece) {
            const NodeId a = cell[static_cast<std::size_t>(i)];
            const NodeId b = cell[static_cast<std::size_t>(j)];
            edge_pieces.push_back({std::min(a, b), std::max(a, b), piece});
            const double v = detail::distance(pts[static_cast<std::size_t>(a)],
                                              pts[static_cast<std::size_t>(b)]) *
                             piece / static_cast<double>(2 * D);
            for (const NodeId q : {a, b}) {
                volume[static_cast<std::size_t>(q)] += v;
                volume_magnitude[static_cast<std::size_t>(q)] += std::abs(v);
            }
        };
        if (dim == 2) {
            const auto piece = detail::triangle_couplings({P(0), P(1), P(2)});
            for (std::size_t e = 0; e < 3; ++e) {
                add(detail::triangle_edges[e][0], detail::triangle_edges[e][1], piece[e]);
            }
        } else {
            const auto piece = detail::tetrahedron_couplings({P(0), P(1), P(2), P(3)});
            for (std::size_t e = 0; e < 6; ++e) {
                add(detail::tetrahedron_edges[e][0], detail::tetrahedron_edges[e][1], piece[e]);
            }
        }
    }
    std::ranges::stable_sort(edge_pieces, [](const Piece& x, const Piece& y) {
        return std::tie(x.a, x.b) < std::tie(y.a, y.b);
    });

    struct Total {
        NodeId first;
        NodeId second;
        double length;
        double area;
        double magnitude;
    };
    std::vector<Total> totals;
    for (std::size_t g = 0; g < edge_pieces.size();) {
        Total t{edge_pieces[g].a, edge_pieces[g].b, 0.0, 0.0, 0.0};
        for (; g < edge_pieces.size() && edge_pieces[g].a == t.first && edge_pieces[g].b == t.second;
             ++g) {
            t.area += edge_pieces[g].value;
            t.magnitude += std::abs(edge_pieces[g].value);
        }
        t.length = detail::distance(pts[static_cast<std::size_t>(t.first)],
                                    pts[static_cast<std::size_t>(t.second)]);
        totals.push_back(t);
    }

    // Per node: the mean length of its edges (the scale of M2), and the closure sums (P5).
    std::vector<double> length_sum(n, 0.0);
    std::vector<std::size_t> degree(n, 0);
    std::vector<std::array<double, 3>> closure(n, {0.0, 0.0, 0.0});
    std::vector<double> closure_magnitude(n, 0.0);
    for (const Total& t : totals) {
        const auto a = static_cast<std::size_t>(t.first);
        const auto b = static_cast<std::size_t>(t.second);
        length_sum[a] += t.length;
        length_sum[b] += t.length;
        ++degree[a];
        ++degree[b];
        for (std::size_t c = 0; c < 3; ++c) {
            const double f = t.area * (pts[b][c] - pts[a][c]) / t.length;
            closure[a][c] += f;
            closure[b][c] -= f;
        }
        closure_magnitude[a] += t.magnitude;
        closure_magnitude[b] += t.magnitude;
    }
    const auto mean_length = [&](std::size_t i) {
        return length_sum[i] / static_cast<double>(degree[i]);
    };

    // Boundary pieces: closure for every boundary facet; patch areas for the listed ones.
    std::vector<Piece> patch_pieces;
    for (std::size_t f = 0; f < boundary.size(); ++f) {
        const FacetRef& b = boundary[f];
        const Point& o = pts[static_cast<std::size_t>(cells[b.cell][b.local])];
        std::array<Point, 3> q{};
        for (std::size_t k = 0; k < D; ++k) q[k] = pts[static_cast<std::size_t>(b.key[k])];
        std::array<double, 3> piece{};
        std::array<double, 3> normal{0.0, 0.0, 0.0};
        if (dim == 2) {
            const double length = detail::distance(q[0], q[1]);
            piece = {0.5 * length, 0.5 * length, 0.0};
            normal = {(q[1][1] - q[0][1]) / length, -(q[1][0] - q[0][0]) / length, 0.0};
        } else {
            piece = detail::facet_pieces({q[0], q[1], q[2]});
            const std::array<double, 3> u{q[1][0] - q[0][0], q[1][1] - q[0][1], q[1][2] - q[0][2]};
            const std::array<double, 3> v{q[2][0] - q[0][0], q[2][1] - q[0][1], q[2][2] - q[0][2]};
            normal = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
            const double norm = std::hypot(normal[0], normal[1], normal[2]);
            for (double& x : normal) x /= norm;
        }
        // Outward: away from the cell's opposite vertex.
        const double toward = (o[0] - q[0][0]) * normal[0] + (o[1] - q[0][1]) * normal[1] +
                              (o[2] - q[0][2]) * normal[2];
        if (toward > 0.0) {
            for (double& x : normal) x = -x;
        }
        for (std::size_t k = 0; k < D; ++k) {
            const auto i = static_cast<std::size_t>(b.key[k]);
            for (std::size_t c = 0; c < 3; ++c) closure[i][c] += piece[k] * normal[c];
            closure_magnitude[i] += std::abs(piece[k]);
            if (boundary_patch[f] >= 0) {
                patch_pieces.push_back({boundary_patch[f], b.key[k], piece[k]});
            }
        }
    }

    // P5: closure.
    double max_residual = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double r = std::hypot(closure[i][0], closure[i][1], closure[i][2]) /
                         closure_magnitude[i];
        max_residual = std::max(max_residual, r);
        if (!(r <= closure_tolerance)) {
            return std::unexpected(error(ErrorCode::degenerate_mesh,
                                         "control volume does not close (P5)", i, r));
        }
    }

    // P3: volumes (dual_tolerance.hpp).
    for (std::size_t i = 0; i < n; ++i) {
        const double scale = detail::total_scale(mean_length(i), dim, volume_magnitude[i]);
        if (!detail::accept_volume(volume[i], scale)) {
            return std::unexpected(error(ErrorCode::degenerate_mesh,
                                         "control volume is not positive (P3)", i, volume[i]));
        }
    }

    // P3: couplings; the edges of the graph.
    std::vector<Edge> edges;
    edges.reserve(totals.size());
    std::size_t dropped_edges = 0;
    for (const Total& t : totals) {
        const double scale = detail::total_scale(t.length, dim - 1, t.magnitude);
        const detail::TotalClass c = detail::classify_total(t.area, scale);
        if (c == detail::TotalClass::refuse) {
            return std::unexpected(error(ErrorCode::degenerate_mesh,
                                         "negative coupling area (P3) on edge " +
                                             pair_text(t.first, t.second),
                                         static_cast<std::size_t>(t.first), t.area));
        }
        if (c == detail::TotalClass::zero) {
            ++dropped_edges;
            continue;
        }
        edges.push_back(Edge{.first = t.first, .second = t.second, .length = t.length,
                             .coupling_area = t.area});
    }

    // P3: patch areas; the boundary patches.
    std::ranges::stable_sort(patch_pieces, [](const Piece& x, const Piece& y) {
        return std::tie(x.a, x.b) < std::tie(y.a, y.b);
    });
    std::vector<BoundaryPatch> patches(s.patch_names.size());
    for (std::size_t p = 0; p < patches.size(); ++p) patches[p].name = s.patch_names[p];
    std::size_t dropped_patch_nodes = 0;
    for (std::size_t g = 0; g < patch_pieces.size();) {
        const NodeId patch = patch_pieces[g].a;
        const NodeId node = patch_pieces[g].b;
        double area = 0.0;
        double magnitude = 0.0;
        for (; g < patch_pieces.size() && patch_pieces[g].a == patch && patch_pieces[g].b == node;
             ++g) {
            area += patch_pieces[g].value;
            magnitude += std::abs(patch_pieces[g].value);
        }
        const auto i = static_cast<std::size_t>(node);
        const double scale = detail::total_scale(mean_length(i), dim - 1, magnitude);
        const detail::TotalClass c = detail::classify_total(area, scale);
        if (c == detail::TotalClass::refuse) {
            return std::unexpected(error(ErrorCode::degenerate_mesh,
                                         "negative boundary area (P3) in patch '" +
                                             s.patch_names[static_cast<std::size_t>(patch)] + "'",
                                         i, area));
        }
        if (c == detail::TotalClass::zero) {
            ++dropped_patch_nodes;
            continue;
        }
        BoundaryPatch& bp = patches[static_cast<std::size_t>(patch)];
        bp.nodes.push_back(node);
        bp.areas.push_back(area);
    }

    if (report != nullptr) {
        *report = DualMeshReport{.dropped_edges = dropped_edges,
                                 .dropped_patch_nodes = dropped_patch_nodes,
                                 .boundary_facets = boundary.size(),
                                 .max_closure_residual = max_residual};
    }
    return Mesh::from_parts(dim, pts, std::move(volume), std::move(edges), std::move(patches));
}

}  // namespace NiTCAD::mesh
