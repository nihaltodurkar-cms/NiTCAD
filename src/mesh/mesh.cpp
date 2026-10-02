#include "NiTCAD/mesh/mesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <utility>

namespace NiTCAD::mesh {

namespace {

base::Error error(base::ErrorCode code, std::string message, std::size_t index,
                  std::optional<double> value = std::nullopt) {
    return {code, std::move(message), base::ErrorContext{.index = index, .value = value}};
}

bool positive_finite(double v) noexcept { return std::isfinite(v) && v > 0.0; }

// Relative tolerance between a producer's edge length and its endpoints' distance: a few
// roundings of the coordinate differences.
constexpr double length_tolerance = 1e-12;

}  // namespace

std::expected<Mesh, base::Error> Mesh::from_parts(int dimension, std::vector<Point> points,
                                                  std::vector<double> volumes,
                                                  std::vector<Edge> edges,
                                                  std::vector<BoundaryPatch> boundary) {
    if (dimension < 1 || dimension > 3) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           "mesh dimension must be 1, 2 or 3", std::nullopt});
    }
    if (points.empty() || volumes.size() != points.size()) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           "mesh needs one volume per node and at least one node",
                                           std::nullopt});
    }
    if (points.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                           "too many mesh nodes for NodeId", std::nullopt});
    }
    const auto dim = static_cast<std::size_t>(dimension);
    const auto n = static_cast<NodeId>(points.size());

    for (std::size_t i = 0; i < points.size(); ++i) {
        for (std::size_t c = 0; c < 3; ++c) {
            const double v = points[i][c];
            if (!std::isfinite(v) || (c >= dim && v != 0.0)) {
                return std::unexpected(error(base::ErrorCode::invalid_input,
                                             "node coordinate is not finite, or nonzero beyond "
                                             "the mesh dimension",
                                             i, v));
            }
        }
        if (!positive_finite(volumes[i])) {
            return std::unexpected(error(base::ErrorCode::degenerate_mesh,
                                         "control volume is not positive", i, volumes[i]));
        }
    }

    std::vector<std::pair<NodeId, NodeId>> pairs;
    pairs.reserve(edges.size());
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const Edge& edge = edges[e];
        if (edge.first < 0 || edge.second >= n || edge.first >= edge.second) {
            return std::unexpected(error(base::ErrorCode::invalid_input,
                                         "edge endpoints out of range or not first < second", e));
        }
        if (!positive_finite(edge.length)) {
            return std::unexpected(error(base::ErrorCode::degenerate_mesh,
                                         "edge length is not positive", e, edge.length));
        }
        if (!positive_finite(edge.coupling_area)) {
            // Also catches the legacy hazard of negative dual areas (6.11, item 4).
            return std::unexpected(error(base::ErrorCode::degenerate_mesh,
                                         "edge coupling area is not positive", e,
                                         edge.coupling_area));
        }
        const Point& a = points[static_cast<std::size_t>(edge.first)];
        const Point& b = points[static_cast<std::size_t>(edge.second)];
        const double distance = std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
        if (std::abs(distance - edge.length) > length_tolerance * edge.length) {
            return std::unexpected(error(base::ErrorCode::degenerate_mesh,
                                         "edge length differs from its endpoints' distance", e,
                                         edge.length));
        }
        pairs.emplace_back(edge.first, edge.second);
    }
    {
        std::vector<std::size_t> order(pairs.size());
        for (std::size_t e = 0; e < order.size(); ++e) order[e] = e;
        std::ranges::sort(order, [&](std::size_t l, std::size_t r) {
            return pairs[l] != pairs[r] ? pairs[l] < pairs[r] : l < r;
        });
        for (std::size_t k = 1; k < order.size(); ++k) {
            if (pairs[order[k]] == pairs[order[k - 1]]) {
                return std::unexpected(
                    error(base::ErrorCode::invalid_input, "edge listed twice", order[k]));
            }
        }
    }

    std::set<std::string, std::less<>> names;
    for (std::size_t p = 0; p < boundary.size(); ++p) {
        const BoundaryPatch& patch = boundary[p];
        if (patch.name.empty() || !names.insert(patch.name).second) {
            return std::unexpected(error(base::ErrorCode::invalid_input,
                                         "boundary patch name is empty or repeated", p));
        }
        if (patch.areas.size() != patch.nodes.size()) {
            return std::unexpected(error(base::ErrorCode::invalid_input,
                                         "boundary patch needs one area per node", p));
        }
        for (std::size_t k = 0; k < patch.nodes.size(); ++k) {
            const NodeId node = patch.nodes[k];
            if (node < 0 || node >= n || (k > 0 && node <= patch.nodes[k - 1])) {
                return std::unexpected(error(base::ErrorCode::invalid_input,
                                             "boundary patch '" + patch.name +
                                                 "': nodes out of range or not increasing",
                                             k));
            }
            if (!positive_finite(patch.areas[k])) {
                return std::unexpected(error(base::ErrorCode::degenerate_mesh,
                                             "boundary patch '" + patch.name +
                                                 "': face area is not positive",
                                             k, patch.areas[k]));
            }
        }
    }

    Mesh mesh;
    mesh.dimension_ = dimension;
    mesh.points_ = std::move(points);
    mesh.volumes_ = std::move(volumes);
    mesh.edges_ = std::move(edges);
    mesh.boundary_ = std::move(boundary);
    return mesh;
}

const BoundaryPatch* Mesh::find_boundary(std::string_view name) const noexcept {
    const auto it = std::ranges::find(boundary_, name, &BoundaryPatch::name);
    return it == boundary_.end() ? nullptr : &*it;
}

}  // namespace NiTCAD::mesh
