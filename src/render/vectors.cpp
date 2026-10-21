#include "NiTCAD/render/vectors.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>

namespace NiTCAD::render {

namespace {

constexpr double nan = std::numeric_limits<double>::quiet_NaN();

base::Error invalid(std::string message, std::size_t index) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

bool same_mesh(const device::Device& d, const Geometry& g) {
    if (&d.mesh() == &g.mesh()) return true;
    const auto a = d.mesh().points(), b = g.mesh().points();
    return a.size() == b.size() && std::ranges::equal(a, b) &&
           d.mesh().dimension() == g.mesh().dimension();
}

// Solves the D x D system m x = r (row-major in the leading D x D block) by Cramer's rule; false when
// m is singular relative to its size (|det| below 1e-10 of the product of its rows' norms).
bool solve(int D, const std::array<std::array<double, 3>, 3>& m, const Point& r, Point& x) noexcept {
    const auto row_norm = [&](int i) {
        double s = 0.0;
        for (int j = 0; j < D; ++j) s += m[i][j] * m[i][j];
        return std::sqrt(s);
    };
    double scale = 1.0;
    for (int i = 0; i < D; ++i) scale *= row_norm(i);
    x = {0.0, 0.0, 0.0};
    if (D == 1) {
        if (!(std::abs(m[0][0]) > 1e-10 * scale)) return false;
        x[0] = r[0] / m[0][0];
        return true;
    }
    if (D == 2) {
        const double det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        if (!(std::abs(det) > 1e-10 * scale)) return false;
        x[0] = (r[0] * m[1][1] - m[0][1] * r[1]) / det;
        x[1] = (m[0][0] * r[1] - r[0] * m[1][0]) / det;
        return true;
    }
    const auto det3 = [](const std::array<std::array<double, 3>, 3>& a) {
        return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
               a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
               a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    };
    const double det = det3(m);
    if (!(std::abs(det) > 1e-10 * scale)) return false;
    for (int c = 0; c < 3; ++c) {
        auto mc = m;
        for (int i = 0; i < 3; ++i) mc[i][c] = r[i];
        x[c] = det3(mc) / det;
    }
    return true;
}

std::uint64_t edge_key(NodeId a, NodeId b) noexcept {
    const auto lo = static_cast<std::uint64_t>(std::min(a, b));
    const auto hi = static_cast<std::uint64_t>(std::max(a, b));
    return (hi << 32) | lo;
}

}  // namespace

std::expected<ElementVectors, base::Error> element_vectors(const device::Device& device,
                                                           const Geometry& g,
                                                           const results::BiasPoint& point,
                                                           VectorField field) {
    if (!same_mesh(device, g)) {
        return std::unexpected(invalid("the geometry's mesh is not the device's", 0));
    }
    if (g.element_count() == 0) {
        return std::unexpected(invalid("the geometry has no elements to reconstruct vectors on", 0));
    }
    const mesh::Mesh& m = device.mesh();
    const int D = m.dimension();
    const auto points = m.points();
    const std::size_t E = g.element_count();
    ElementVectors out{std::vector<Point>(E, Point{}), std::vector<std::uint8_t>(E, 0)};

    if (field == VectorField::electric) {
        const auto& psi = point.fields.potential_V;
        if (psi.size() != m.node_count()) {
            return std::unexpected(invalid("the result's potential is not one value per mesh node", psi.size()));
        }
        for (std::size_t e = 0; e < E; ++e) {
            const auto nodes = g.element(e);
            const Point& x0 = points[static_cast<std::size_t>(nodes[0])];
            std::array<std::array<double, 3>, 3> a{};
            Point r{};
            for (int i = 0; i < D; ++i) {
                const Point& xi = points[static_cast<std::size_t>(nodes[static_cast<std::size_t>(i + 1)])];
                for (int j = 0; j < D; ++j) a[i][j] = xi[j] - x0[j];
                r[i] = psi[static_cast<std::size_t>(nodes[static_cast<std::size_t>(i + 1)])] -
                       psi[static_cast<std::size_t>(nodes[0])];
            }
            Point grad{};
            if (!solve(D, a, r, grad)) continue;
            out.value[e] = {-grad[0], -grad[1], -grad[2]};
            out.defined[e] = 1;
        }
        return out;
    }

    const auto edges = m.edges();
    const bool electrons = field != VectorField::hole_current;
    const bool holes = field != VectorField::electron_current;
    if ((electrons && point.edge_current_n.size() != edges.size()) ||
        (holes && point.edge_current_p.size() != edges.size())) {
        return std::unexpected(invalid(
            "the result holds no edge currents, one per mesh edge",
            electrons ? point.edge_current_n.size() : point.edge_current_p.size()));
    }
    std::unordered_map<std::uint64_t, std::size_t> edge_of;
    edge_of.reserve(edges.size());
    for (std::size_t k = 0; k < edges.size(); ++k) edge_of.emplace(edge_key(edges[k].first, edges[k].second), k);
    std::vector<std::uint8_t> semiconductor(m.node_count());
    for (std::size_t i = 0; i < semiconductor.size(); ++i) {
        semiconductor[i] = device.is_insulator(static_cast<NodeId>(i)) ? 0 : 1;
    }

    for (std::size_t e = 0; e < E; ++e) {
        const auto nodes = g.element(e);
        // Normal equations of the fit: sum t t^T J = sum t (I / A) over the usable edges.
        std::array<std::array<double, 3>, 3> mt{};
        Point rhs{};
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            for (std::size_t j = i + 1; j < nodes.size(); ++j) {
                const auto it = edge_of.find(edge_key(nodes[i], nodes[j]));
                if (it == edge_of.end()) continue;
                const mesh::Edge& ed = edges[it->second];
                if (!(ed.coupling_area > 0.0) || !(ed.length > 0.0) ||
                    !semiconductor[static_cast<std::size_t>(ed.first)] ||
                    !semiconductor[static_cast<std::size_t>(ed.second)]) {
                    continue;
                }
                const Point& a = points[static_cast<std::size_t>(ed.first)];
                const Point& b = points[static_cast<std::size_t>(ed.second)];
                Point t{};
                for (int c = 0; c < D; ++c) t[c] = (b[c] - a[c]) / ed.length;
                double current = 0.0;
                if (electrons) current += point.edge_current_n[it->second];
                if (holes) current += point.edge_current_p[it->second];
                const double density = current / ed.coupling_area;
                for (int r = 0; r < D; ++r) {
                    for (int c = 0; c < D; ++c) mt[r][c] += t[r] * t[c];
                    rhs[r] += t[r] * density;
                }
            }
        }
        Point j{};
        if (!solve(D, mt, rhs, j)) continue;
        out.value[e] = j;
        out.defined[e] = 1;
    }
    return out;
}

double element_volume(const Geometry& g, std::size_t e) noexcept {
    const auto nodes = g.element(e);
    const auto points = g.mesh().points();
    const int D = g.dimension();
    const Point& x0 = points[static_cast<std::size_t>(nodes[0])];
    const auto d = [&](std::size_t k, int c) { return points[static_cast<std::size_t>(nodes[k])][c] - x0[c]; };
    if (D == 1) return std::abs(d(1, 0));
    if (D == 2) return 0.5 * std::abs(d(1, 0) * d(2, 1) - d(1, 1) * d(2, 0));
    const double det = d(1, 0) * (d(2, 1) * d(3, 2) - d(2, 2) * d(3, 1)) -
                       d(1, 1) * (d(2, 0) * d(3, 2) - d(2, 2) * d(3, 0)) +
                       d(1, 2) * (d(2, 0) * d(3, 1) - d(2, 1) * d(3, 0));
    return std::abs(det) / 6.0;
}

std::array<double, 4> control_volume_shares(const Geometry& g, std::size_t e) noexcept {
    const int D = g.dimension();
    std::array<double, 4> share{};
    const auto nodes = g.element(e);
    if (g.kind() != GeometryKind::tensor) {
        for (std::size_t k = 0; k < nodes.size(); ++k) share[k] = 1.0 / static_cast<double>(nodes.size());
        return share;
    }
    // The vertex m steps along the Kuhn path is the corner with m upper bits: binomial(D, m) / 2^D.
    constexpr std::array<std::array<double, 4>, 4> binomial{
        {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 1, 0}, {1, 3, 3, 1}}};
    const auto corners = g.cell_corners(g.element_cell(e));
    for (std::size_t k = 0; k < nodes.size(); ++k) {
        const auto at = std::ranges::find(corners, nodes[k]);
        const auto bits = static_cast<unsigned>(at - corners.begin());
        share[k] = binomial[static_cast<std::size_t>(D)][static_cast<std::size_t>(std::popcount(bits))] /
                   static_cast<double>(1 << D);
    }
    return share;
}

std::vector<Point> nodal_vectors(const Geometry& g, const ElementVectors& v) {
    const std::size_t n = g.mesh().node_count();
    std::vector<Point> sum(n, Point{});
    std::vector<double> weight(n, 0.0);
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        if (!v.defined[e]) continue;
        const auto nodes = g.element(e);
        const double volume = element_volume(g, e);
        const auto share = control_volume_shares(g, e);
        for (std::size_t k = 0; k < nodes.size(); ++k) {
            const double w = volume * share[k];
            const auto i = static_cast<std::size_t>(nodes[k]);
            for (std::size_t c = 0; c < 3; ++c) sum[i][c] += w * v.value[e][c];
            weight[i] += w;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (weight[i] > 0.0) {
            for (std::size_t c = 0; c < 3; ++c) sum[i][c] /= weight[i];
        } else {
            sum[i] = {nan, nan, nan};
        }
    }
    return sum;
}

std::vector<double> magnitudes(std::span<const Point> v) {
    std::vector<double> out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        out[i] = std::sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
    }
    return out;
}

}  // namespace NiTCAD::render
