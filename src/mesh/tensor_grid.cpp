#include "NiTCAD/mesh/tensor_grid.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace NiTCAD::mesh {

namespace {

constexpr std::array<const char*, 3> axis_names{"x", "y", "z"};

// Box-integration control-volume widths along one axis: half of each adjacent spacing (legacy
// mesh2d.control_volume_widths).
std::vector<double> control_volume_widths(std::span<const double> nodes) {
    std::vector<double> widths(nodes.size(), 0.0);
    for (std::size_t i = 0; i + 1 < nodes.size(); ++i) {
        const double half = 0.5 * (nodes[i + 1] - nodes[i]);
        widths[i] += half;
        widths[i + 1] += half;
    }
    return widths;
}

std::expected<Mesh, base::Error> build(std::span<const std::span<const double>> axes) {
    const std::size_t dim = axes.size();
    std::size_t count = 1;
    for (std::size_t a = 0; a < dim; ++a) {
        const std::span<const double> axis = axes[a];
        const std::string name = axis_names[a];
        if (axis.size() < 2) {
            return std::unexpected(base::Error{base::ErrorCode::degenerate_mesh,
                                               "axis " + name + " needs at least two nodes",
                                               base::ErrorContext{.index = axis.size(),
                                                                  .value = std::nullopt}});
        }
        for (std::size_t i = 0; i < axis.size(); ++i) {
            if (!std::isfinite(axis[i]) || (i > 0 && !(axis[i] > axis[i - 1]))) {
                return std::unexpected(base::Error{
                    base::ErrorCode::degenerate_mesh,
                    "axis " + name + " is not finite and strictly increasing",
                    base::ErrorContext{.index = i, .value = axis[i]}});
            }
        }
        if (axis.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max()) / count) {
            return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                               "tensor grid has too many nodes for NodeId",
                                               std::nullopt});
        }
        count *= axis.size();
    }

    std::array<std::vector<double>, 3> widths;
    std::array<std::size_t, 3> size{1, 1, 1};
    std::array<std::size_t, 3> stride{0, 0, 0};
    for (std::size_t a = 0, s = 1; a < dim; ++a) {
        widths[a] = control_volume_widths(axes[a]);
        size[a] = axes[a].size();
        stride[a] = s;
        s *= size[a];
    }

    // Nodes in x-fastest order. The order is this producer's choice; callers see only the graph.
    std::vector<Point> points(count, Point{0.0, 0.0, 0.0});
    std::vector<double> volumes(count, 1.0);
    std::vector<Edge> edges;
    edges.reserve(count * dim);
    for (std::size_t node = 0; node < count; ++node) {
        std::array<std::size_t, 3> at{0, 0, 0};
        for (std::size_t a = 0; a < dim; ++a) {
            at[a] = (node / stride[a]) % size[a];
            points[node][a] = axes[a][at[a]];
            volumes[node] *= widths[a][at[a]];
        }
        for (std::size_t a = 0; a < dim; ++a) {
            if (at[a] + 1 == size[a]) continue;
            double area = 1.0;
            for (std::size_t b = 0; b < dim; ++b) {
                if (b != a) area *= widths[b][at[b]];
            }
            edges.push_back(Edge{.first = static_cast<NodeId>(node),
                                 .second = static_cast<NodeId>(node + stride[a]),
                                 .length = axes[a][at[a] + 1] - axes[a][at[a]],
                                 .coupling_area = area});
        }
    }

    std::vector<BoundaryPatch> boundary;
    for (std::size_t a = 0; a < dim; ++a) {
        for (const bool high : {false, true}) {
            BoundaryPatch patch{.name = std::string(axis_names[a]) + (high ? "_max" : "_min"),
                                .nodes = {},
                                .areas = {}};
            const std::size_t face = high ? size[a] - 1 : 0;
            for (std::size_t node = 0; node < count; ++node) {
                if ((node / stride[a]) % size[a] != face) continue;
                double area = 1.0;
                for (std::size_t b = 0; b < dim; ++b) {
                    if (b != a) area *= widths[b][(node / stride[b]) % size[b]];
                }
                patch.nodes.push_back(static_cast<NodeId>(node));
                patch.areas.push_back(area);
            }
            boundary.push_back(std::move(patch));
        }
    }

    return Mesh::from_parts(static_cast<int>(dim), std::move(points), std::move(volumes),
                            std::move(edges), std::move(boundary));
}

}  // namespace

std::expected<Mesh, base::Error> make_tensor_grid(std::span<const double> x) {
    const std::array<std::span<const double>, 1> axes{x};
    return build(axes);
}

std::expected<Mesh, base::Error> make_tensor_grid(std::span<const double> x,
                                                  std::span<const double> y) {
    const std::array<std::span<const double>, 2> axes{x, y};
    return build(axes);
}

std::expected<Mesh, base::Error> make_tensor_grid(std::span<const double> x,
                                                  std::span<const double> y,
                                                  std::span<const double> z) {
    const std::array<std::span<const double>, 3> axes{x, y, z};
    return build(axes);
}

std::expected<std::vector<double>, base::Error> straddle_interface(std::vector<double> axis,
                                                                  double position,
                                                                  double spacing) {
    const auto invalid = [](const char* message, double value) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input, message,
                                           base::ErrorContext{.index = std::nullopt,
                                                              .value = value}});
    };
    for (std::size_t i = 0; i < axis.size(); ++i) {
        if (!std::isfinite(axis[i]) || (i > 0 && !(axis[i] > axis[i - 1]))) {
            return invalid("axis is not finite and strictly increasing", axis[i]);
        }
    }
    if (!(std::isfinite(spacing) && spacing > 0.0)) {
        return invalid("interface spacing must be finite and positive", spacing);
    }
    if (axis.size() < 2 ||
        !(position - spacing > axis.front() && position + spacing < axis.back())) {
        return invalid("interface position must lie inside the axis by more than the spacing",
                       position);
    }
    std::vector<double> out;
    out.reserve(axis.size() + 2);
    bool inserted = false;
    for (const double x : axis) {
        if (!inserted && x > position) {
            out.push_back(position - 0.5 * spacing);
            out.push_back(position + 0.5 * spacing);
            inserted = true;
        }
        if (std::abs(x - position) < spacing) continue;
        if (!inserted || x > position) out.push_back(x);
    }
    return out;
}

}  // namespace NiTCAD::mesh
