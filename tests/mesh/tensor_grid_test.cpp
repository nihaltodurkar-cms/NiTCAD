// Geometry gates for the tensor-grid producer (ARCHITECTURE.md section 11, Unit 4). Every check
// goes through the public graph (points, volumes, edges, boundary patches); none uses grid indices.
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

// Deliberately non-uniform axes (graded, with a 1e3 spacing ratio on x), in cm.
const std::vector<double> ax{0.0, 1e-7, 3e-7, 7e-7, 1.5e-6, 1.0e-5, 4.0e-5, 1.0e-4};
const std::vector<double> ay{0.0, 2.0e-6, 2.5e-6, 9.0e-6, 3.0e-5};
const std::vector<double> az{0.0, 5.0e-6, 6.0e-6, 2.0e-5};

double extent(const std::vector<double>& axis) { return axis.back() - axis.front(); }

double total(std::span<const double> v) {
    double s = 0.0;
    for (const double x : v) s += x;
    return s;
}

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// Outward unit normal of a tensor-grid boundary patch, from its name.
std::array<double, 3> outward_normal(const std::string& name) {
    std::array<double, 3> n{0.0, 0.0, 0.0};
    const std::size_t axis = static_cast<std::size_t>(name[0] - 'x');
    n[axis] = name.ends_with("_max") ? 1.0 : -1.0;
    return n;
}

// Discrete Gauss identity of the box method. For linear u(p) = g . p, the outward flux of grad u
// through each control volume's surface is zero: the internal faces give, over its edges,
// area * (u_other - u_self) / length, and its boundary faces give, over its patches, (g . n) * area. It ties volumes, edge geometry and
// boundary areas together, in any dimension.
double worst_flux_imbalance(const Mesh& mesh, const std::array<double, 3>& g) {
    const auto u = [&](NodeId k) {
        const Point& p = mesh.points()[static_cast<std::size_t>(k)];
        return g[0] * p[0] + g[1] * p[1] + g[2] * p[2];
    };
    std::vector<double> balance(mesh.node_count(), 0.0);
    std::vector<double> scale(mesh.node_count(), 0.0);
    for (const Edge& e : mesh.edges()) {
        const double flux = e.coupling_area * (u(e.second) - u(e.first)) / e.length;
        balance[static_cast<std::size_t>(e.first)] += flux;
        balance[static_cast<std::size_t>(e.second)] -= flux;
        scale[static_cast<std::size_t>(e.first)] += std::abs(flux);
        scale[static_cast<std::size_t>(e.second)] += std::abs(flux);
    }
    for (const BoundaryPatch& patch : mesh.boundary()) {
        const std::array<double, 3> n = outward_normal(patch.name);
        const double gn = g[0] * n[0] + g[1] * n[1] + g[2] * n[2];
        for (std::size_t k = 0; k < patch.nodes.size(); ++k) {
            const auto i = static_cast<std::size_t>(patch.nodes[k]);
            balance[i] += gn * patch.areas[k];
            scale[i] += std::abs(gn * patch.areas[k]);
        }
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < balance.size(); ++i) {
        if (scale[i] > 0.0) worst = std::max(worst, std::abs(balance[i]) / scale[i]);
    }
    return worst;
}

}  // namespace

TEST_CASE("tensor grid: 1D matches the legacy box method") {
    const Mesh m = make_tensor_grid(ax).value();
    REQUIRE(m.dimension() == 1);
    REQUIRE(m.node_count() == ax.size());
    REQUIRE(m.edges().size() == ax.size() - 1);
    // Control volume: half of each adjacent spacing (legacy control_volume_widths).
    REQUIRE(m.volumes()[0] == 0.5 * (ax[1] - ax[0]));
    REQUIRE(m.volumes()[3] == 0.5 * (ax[3] - ax[2]) + 0.5 * (ax[4] - ax[3]));
    REQUIRE(m.volumes().back() == 0.5 * (ax[7] - ax[6]));
    for (const Edge& e : m.edges()) {
        REQUIRE(e.coupling_area == 1.0);  // per unit cross-section
        REQUIRE(e.length == m.points()[static_cast<std::size_t>(e.second)][0] -
                                m.points()[static_cast<std::size_t>(e.first)][0]);
    }
    REQUIRE(m.boundary().size() == 2);
    const BoundaryPatch* left = m.find_boundary("x_min");
    REQUIRE(left != nullptr);
    REQUIRE(left->nodes.size() == 1);
    REQUIRE(m.points()[static_cast<std::size_t>(left->nodes[0])][0] == 0.0);
    REQUIRE(left->areas[0] == 1.0);
    REQUIRE(m.find_boundary("y_min") == nullptr);
}

TEST_CASE("tensor grid: node, edge and patch counts in 2D and 3D") {
    const Mesh m2 = make_tensor_grid(ax, ay).value();
    REQUIRE(m2.node_count() == ax.size() * ay.size());
    REQUIRE(m2.edges().size() == (ax.size() - 1) * ay.size() + ax.size() * (ay.size() - 1));
    REQUIRE(m2.boundary().size() == 4);
    REQUIRE(m2.find_boundary("y_max")->nodes.size() == ax.size());

    const Mesh m3 = make_tensor_grid(ax, ay, az).value();
    const std::size_t nx = ax.size(), ny = ay.size(), nz = az.size();
    REQUIRE(m3.node_count() == nx * ny * nz);
    REQUIRE(m3.edges().size() ==
            (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1));
    REQUIRE(m3.boundary().size() == 6);
    REQUIRE(m3.find_boundary("z_min")->nodes.size() == nx * ny);
    for (const Point& p : m2.points()) REQUIRE(p[2] == 0.0);
}

TEST_CASE("tensor grid: total volume equals the domain (section 10 gate)") {
    // Legacy gates: 1e-14 in 2D, 1e-10 in 3D (absolute, on O(1) domains). Here, relative, on
    // micrometre-scale graded axes.
    const Mesh m1 = make_tensor_grid(ax).value();
    const Mesh m2 = make_tensor_grid(ax, ay).value();
    const Mesh m3 = make_tensor_grid(ax, ay, az).value();
    REQUIRE(close(total(m1.volumes()), extent(ax), 1e-14));
    REQUIRE(close(total(m2.volumes()), extent(ax) * extent(ay), 1e-14));
    REQUIRE(close(total(m3.volumes()), extent(ax) * extent(ay) * extent(az), 1e-13));
}

TEST_CASE("tensor grid: volumes, lengths and areas are positive") {
    for (const Mesh& m : {make_tensor_grid(ax).value(), make_tensor_grid(ax, ay).value(),
                          make_tensor_grid(ax, ay, az).value()}) {
        for (const double v : m.volumes()) REQUIRE(v > 0.0);
        for (const Edge& e : m.edges()) {
            REQUIRE(e.length > 0.0);
            REQUIRE(e.coupling_area > 0.0);
            REQUIRE(e.first < e.second);
        }
        for (const BoundaryPatch& p : m.boundary()) {
            for (const double a : p.areas) REQUIRE(a > 0.0);
        }
    }
}

TEST_CASE("tensor grid: each boundary face's area equals the face's measure") {
    const Mesh m3 = make_tensor_grid(ax, ay, az).value();
    REQUIRE(close(total(m3.find_boundary("x_min")->areas), extent(ay) * extent(az), 1e-14));
    REQUIRE(close(total(m3.find_boundary("y_max")->areas), extent(ax) * extent(az), 1e-14));
    REQUIRE(close(total(m3.find_boundary("z_max")->areas), extent(ax) * extent(ay), 1e-14));
    const Mesh m2 = make_tensor_grid(ax, ay).value();
    REQUIRE(close(total(m2.find_boundary("x_max")->areas), extent(ay), 1e-14));
    REQUIRE(close(total(m2.find_boundary("y_min")->areas), extent(ax), 1e-14));
}

TEST_CASE("tensor grid: 3D and 2D reduce to 1D (consistent edge geometry across D)") {
    // Grouping by x coordinate only (no grid indices): the nodes of a 2D/3D grid at one x carry
    // the 1D control volume times the transverse area, and the x-edges between two x values carry
    // the 1D coupling (1) times the same area, with the same length.
    const Mesh m1 = make_tensor_grid(ax).value();
    for (const auto& [m, transverse] :
         {std::pair{make_tensor_grid(ax, ay).value(), extent(ay)},
          std::pair{make_tensor_grid(ax, ay, az).value(), extent(ay) * extent(az)}}) {
        std::map<double, double> volume_at_x;
        for (std::size_t k = 0; k < m.node_count(); ++k) {
            volume_at_x[m.points()[k][0]] += m.volumes()[k];
        }
        std::map<std::pair<double, double>, std::pair<double, double>> x_edges;  // area, length
        for (const Edge& e : m.edges()) {
            const Point& a = m.points()[static_cast<std::size_t>(e.first)];
            const Point& b = m.points()[static_cast<std::size_t>(e.second)];
            if (a[1] != b[1] || a[2] != b[2]) continue;  // not along x
            auto& [area, length] = x_edges[{a[0], b[0]}];
            area += e.coupling_area;
            length = e.length;
        }
        REQUIRE(volume_at_x.size() == m1.node_count());
        REQUIRE(x_edges.size() == m1.edges().size());
        for (std::size_t k = 0; k < m1.node_count(); ++k) {
            REQUIRE(close(volume_at_x.at(m1.points()[k][0]), m1.volumes()[k] * transverse, 1e-14));
        }
        for (const Edge& e : m1.edges()) {
            const double xa = m1.points()[static_cast<std::size_t>(e.first)][0];
            const double xb = m1.points()[static_cast<std::size_t>(e.second)][0];
            const auto& [area, length] = x_edges.at({xa, xb});
            REQUIRE(close(area, e.coupling_area * transverse, 1e-14));
            REQUIRE(length == e.length);
        }
    }
}

TEST_CASE("tensor grid: the box method is exact for linear fields (discrete Gauss identity)") {
    const std::array<double, 3> g{0.7, -1.3, 0.4};
    REQUIRE(worst_flux_imbalance(make_tensor_grid(ax).value(), g) <= 1e-12);
    REQUIRE(worst_flux_imbalance(make_tensor_grid(ax, ay).value(), g) <= 1e-12);
    REQUIRE(worst_flux_imbalance(make_tensor_grid(ax, ay, az).value(), g) <= 1e-12);
}

TEST_CASE("tensor grid: the Gauss identity detects a 1% error in one coupling area") {
    const Mesh good = make_tensor_grid(ax, ay, az).value();
    std::vector<Edge> edges(good.edges().begin(), good.edges().end());
    edges[edges.size() / 2].coupling_area *= 1.01;
    const Mesh bad = Mesh::from_parts(3, {good.points().begin(), good.points().end()},
                                      {good.volumes().begin(), good.volumes().end()},
                                      std::move(edges),
                                      {good.boundary().begin(), good.boundary().end()})
                         .value();
    REQUIRE(worst_flux_imbalance(bad, {0.7, -1.3, 0.4}) > 1e-4);
}

TEST_CASE("tensor grid: invalid axes are degenerate_mesh with their position") {
    const std::vector<double> one{0.0};
    const std::vector<double> flat{0.0, 1.0, 1.0, 2.0};
    const std::vector<double> back{0.0, 2.0, 1.0};
    const std::vector<double> nan{0.0, std::numeric_limits<double>::quiet_NaN(), 1.0};

    REQUIRE(make_tensor_grid(one).error().code == ErrorCode::degenerate_mesh);
    const auto f = make_tensor_grid(ax, flat);
    REQUIRE(f.error().code == ErrorCode::degenerate_mesh);
    REQUIRE(f.error().context->index == 2);
    REQUIRE(f.error().message.starts_with("axis y"));
    REQUIRE(make_tensor_grid(ax, ay, back).error().context->index == 2);
    REQUIRE(make_tensor_grid(nan).error().context->index == 1);
}
