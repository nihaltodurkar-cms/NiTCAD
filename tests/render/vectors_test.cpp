// Vectors reconstructed from results (Unit 25b, DECISIONS.md DS1-DS3, BG1): E of a linear potential
// and J of a uniform current exact on every element of tensor and generated meshes in 1D, 2D and 3D;
// the control-volume shares that average them to the nodes; a solved resistor's |J| = I / area at
// every node; insulators without current; results without edge currents refused.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/render/vectors.hpp"
#include "NiTCAD/solve/bias.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;

namespace {

// Nonuniform spacing, so no element is a scaled copy of another.
std::vector<double> graded(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) {
        const double s = static_cast<double>(k) / (n - 1);
        v.push_back(k == n - 1 ? b : a + (b - a) * (s + 0.15 * s * (1.0 - s)));
    }
    return v;
}

// A device on a mesh: silicon (or oxide where `oxide` says so), n-type 1e16, ohmic contacts on the
// two named boundary patches.
struct Built {
    device::Device device;
    std::optional<mesh::SimplexMesh> simplices;
};

device::Device make_device(mesh::Mesh m, std::optional<mesh::TensorCells> cells, const std::string& left,
                           const std::string& right, bool (*oxide)(const Point&) = nullptr) {
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, oxide ? 1 : 0);
    std::vector<double> donors(n, 1e16), acceptors(n, 0.0);
    if (oxide) {
        for (std::size_t i = 0; i < n; ++i) {
            if (oxide(m.points()[i])) {
                region[i] = 0;
                donors[i] = 0.0;
            }
        }
    }
    const auto contact_nodes = [&](const std::string& patch) {
        std::vector<mesh::NodeId> nodes;
        for (const mesh::NodeId v : m.find_boundary(patch)->nodes) {
            if (region[static_cast<std::size_t>(v)] == (oxide ? 1 : 0)) nodes.push_back(v);
        }
        return nodes;
    };
    auto a = contact_nodes(left);
    auto b = contact_nodes(right);
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = oxide ? std::vector<device::Region>{{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}}
                         : std::vector<device::Region>{{"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(a)},
                      {"right", device::ContactKind::ohmic, std::move(b)}},
         .cells = std::move(cells)});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

Built tensor(int D) {
    const auto x = graded(0.0, 2e-4, 7);
    const auto y = graded(0.0, 1e-4, 5);
    const auto z = graded(0.0, 0.7e-4, 4);
    if (D == 1) {
        mesh::Mesh m = *mesh::make_tensor_grid(x);
        auto cells = *mesh::TensorCells::create(m, x);
        return {make_device(std::move(m), std::move(cells), "x_min", "x_max"), std::nullopt};
    }
    if (D == 2) {
        mesh::Mesh m = *mesh::make_tensor_grid(x, y);
        auto cells = *mesh::TensorCells::create(m, x, y);
        return {make_device(std::move(m), std::move(cells), "x_min", "x_max"), std::nullopt};
    }
    mesh::Mesh m = *mesh::make_tensor_grid(x, y, z);
    auto cells = *mesh::TensorCells::create(m, x, y, z);
    return {make_device(std::move(m), std::move(cells), "x_min", "x_max"), std::nullopt};
}

Built generated(int D) {
    if (D == 2) {
        mesh::Geometry2D geo;
        geo.polygons = {{{{0, 0, 0}, {2e-4, 0, 0}, {2e-4, 1e-4, 0}, {0, 1e-4, 0}}, 0}};
        geo.patches = {{"left", {0, 0, 0}, {0, 1e-4, 0}}, {"right", {2e-4, 0, 0}, {2e-4, 1e-4, 0}}};
        mesh::MeshOptions2D o;
        o.max_size = 0.3e-4;
        auto gm = mesh::generate_mesh(geo, o);
        if (!gm) FAIL(gm.error().message);
        auto s = gm->simplices;
        return {make_device(std::move(gm->mesh), std::nullopt, "left", "right"), std::move(s)};
    }
    mesh::Geometry3D geo;
    const double L = 1e-4;
    geo.domain.vertices = {{0, 0, 0}, {L, 0, 0}, {L, L, 0}, {0, L, 0},
                           {0, 0, L}, {L, 0, L}, {L, L, L}, {0, L, L}};
    geo.domain.faces = {{{0, 3, 2, 1}, "z0", {}}, {{4, 5, 6, 7}, "z1", {}},
                        {{0, 1, 5, 4}, "y0", {}}, {{2, 3, 7, 6}, "y1", {}},
                        {{0, 4, 7, 3}, "x0", {}}, {{1, 2, 6, 5}, "x1", {}}};
    mesh::MeshOptions3D o;
    o.max_size = 0.4e-4;
    auto gm = mesh::generate_mesh(geo, o);
    if (!gm) FAIL(gm.error().message);
    auto s = gm->simplices;
    return {make_device(std::move(gm->mesh), std::nullopt, "x0", "x1"), std::move(s)};
}

Geometry geometry_of(const Built& b) {
    if (b.simplices) return *Geometry::simplex(b.device.mesh(), *b.simplices);
    return *Geometry::tensor(b.device.mesh(), *b.device.cells());
}

double dot(const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm(const Point& a) { return std::sqrt(dot(a, a)); }

bool near(const Point& a, const Point& b, double tol) {
    return norm({a[0] - b[0], a[1] - b[1], a[2] - b[2]}) <= tol * std::max(norm(b), 1e-300);
}

// The meshes BG1 runs on: 1D, 2D and 3D tensor; 2D and 3D generated.
std::vector<std::pair<std::string, Built>> all_meshes() {
    std::vector<std::pair<std::string, Built>> v;
    v.emplace_back("1D tensor", tensor(1));
    v.emplace_back("2D tensor", tensor(2));
    v.emplace_back("3D tensor", tensor(3));
    v.emplace_back("2D generated", generated(2));
    v.emplace_back("3D generated", generated(3));
    return v;
}

Point in_dimension(Point v, int D) {
    for (int c = D; c < 3; ++c) v[static_cast<std::size_t>(c)] = 0.0;
    return v;
}

}  // namespace

TEST_CASE("render vectors: E of a linear potential is exact on every element (BG1)") {
    for (auto& [name, b] : all_meshes()) {
        CAPTURE(name);
        const Geometry g = geometry_of(b);
        const mesh::Mesh& m = b.device.mesh();
        const int D = m.dimension();
        const Point grad = in_dimension({3.0e4, -1.25e4, 0.5e4}, D);  // V/cm
        results::BiasPoint p;
        for (const Point& x : m.points()) p.fields.potential_V.push_back(0.3 + dot(grad, x));
        const auto e = element_vectors(b.device, g, p, VectorField::electric);
        REQUIRE(e.has_value());
        const Point expected{-grad[0], -grad[1], -grad[2]};
        for (std::size_t k = 0; k < g.element_count(); ++k) {
            CAPTURE(k);
            REQUIRE(e->defined[k]);
            CHECK(near(e->value[k], expected, 1e-11));
        }
        const auto nodal = nodal_vectors(g, *e);
        for (const Point& v : nodal) CHECK(near(v, expected, 1e-11));
        const auto mag = *field_values(b.device, g, p, Field::electric_field);
        for (const double v : mag) CHECK(std::abs(v - norm(grad)) <= 1e-11 * norm(grad));
    }
}

TEST_CASE("render vectors: J of a uniform current is exact on every element (BG1)") {
    for (auto& [name, b] : all_meshes()) {
        CAPTURE(name);
        const Geometry g = geometry_of(b);
        const mesh::Mesh& m = b.device.mesh();
        const int D = m.dimension();
        const Point jn = in_dimension({12.0, -7.5, 3.0}, D);  // A/cm^2
        const Point jp = in_dimension({-2.0, 0.5, 9.0}, D);
        results::BiasPoint p;
        for (const mesh::Edge& ed : m.edges()) {
            const Point& a = m.points()[static_cast<std::size_t>(ed.first)];
            const Point& c = m.points()[static_cast<std::size_t>(ed.second)];
            const Point t{(c[0] - a[0]) / ed.length, (c[1] - a[1]) / ed.length, (c[2] - a[2]) / ed.length};
            p.edge_current_n.push_back(dot(jn, t) * ed.coupling_area);
            p.edge_current_p.push_back(dot(jp, t) * ed.coupling_area);
        }
        const Point j{jn[0] + jp[0], jn[1] + jp[1], jn[2] + jp[2]};
        for (const auto& [field, expected] : {std::pair{VectorField::electron_current, jn},
                                             std::pair{VectorField::hole_current, jp},
                                             std::pair{VectorField::current, j}}) {
            const auto e = element_vectors(b.device, g, p, field);
            REQUIRE(e.has_value());
            for (std::size_t k = 0; k < g.element_count(); ++k) {
                CAPTURE(k);
                REQUIRE(e->defined[k]);
                CHECK(near(e->value[k], expected, 1e-11));
            }
            for (const Point& v : nodal_vectors(g, *e)) CHECK(near(v, expected, 1e-11));
        }
        const auto mag = *field_values(b.device, g, p, Field::current_density);
        for (const double v : mag) CHECK(std::abs(v - norm(j)) <= 1e-11 * norm(j));
    }
}

TEST_CASE("render vectors: tensor shares are the box method's control volumes (DS3)") {
    for (const int D : {1, 2, 3}) {
        CAPTURE(D);
        const Built b = tensor(D);
        const Geometry g = geometry_of(b);
        const mesh::Mesh& m = b.device.mesh();
        std::vector<double> volume(m.node_count(), 0.0);
        for (std::size_t e = 0; e < g.element_count(); ++e) {
            const auto share = control_volume_shares(g, e);
            double sum = 0.0;
            for (std::size_t k = 0; k < g.element(e).size(); ++k) {
                sum += share[k];
                volume[static_cast<std::size_t>(g.element(e)[k])] += share[k] * element_volume(g, e);
            }
            CHECK(std::abs(sum - 1.0) <= 1e-15);
        }
        for (std::size_t i = 0; i < m.node_count(); ++i) {
            CAPTURE(i);
            CHECK(std::abs(volume[i] - m.volumes()[i]) <= 1e-12 * m.volumes()[i]);
        }
    }
    // A simplex's shares are equal and its volumes cover the domain.
    const Built b = generated(3);
    const Geometry g = geometry_of(b);
    double total = 0.0;
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        for (const double s : control_volume_shares(g, e)) CHECK(s == 0.25);
        total += element_volume(g, e);
    }
    CHECK(std::abs(total - 1e-12) <= 1e-12 * 1e-12);
}

TEST_CASE("render vectors: a nonuniform E averages within each node's control volume (DS3)") {
    // psi = x^2 on a 1D tensor mesh: each segment's E is -(x_a + x_b); a node's average weighs its
    // two half segments by their lengths.
    const Built b = tensor(1);
    const Geometry g = geometry_of(b);
    const auto x = b.device.mesh().points();
    results::BiasPoint p;
    for (const Point& q : x) p.fields.potential_V.push_back(q[0] * q[0] * 1e8);
    const auto e = *element_vectors(b.device, g, p, VectorField::electric);
    const auto nodal = nodal_vectors(g, e);
    for (std::size_t i = 1; i + 1 < x.size(); ++i) {
        const double hl = x[i][0] - x[i - 1][0], hr = x[i + 1][0] - x[i][0];
        const double el = -(x[i - 1][0] + x[i][0]) * 1e8, er = -(x[i][0] + x[i + 1][0]) * 1e8;
        CHECK(std::abs(nodal[i][0] - (hl * el + hr * er) / (hl + hr)) <= 1e-12 * std::abs(er));
    }
    CHECK(std::abs(nodal[0][0] + (x[0][0] + x[1][0]) * 1e8) <= 1e-6);
}

TEST_CASE("render vectors: insulators carry no current; refusals (DS2)") {
    // Oxide below and above a 2D tensor device (so an oxide edge has it at either end): no current
    // vector on elements with an oxide node, NaN |J| at oxide nodes, E everywhere.
    const auto x = graded(0.0, 2e-4, 7);
    const auto y = graded(0.0, 1e-4, 6);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    auto cells = *mesh::TensorCells::create(m, x, y);
    const device::Device d = make_device(std::move(m), std::move(cells), "x_min", "x_max",
                                         [](const Point& q) { return q[1] > 0.85e-4 || q[1] < 0.15e-4; });
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    results::BiasPoint p;
    p.fields.potential_V.assign(d.mesh().node_count(), 0.0);
    for (const mesh::Edge& ed : d.mesh().edges()) {
        const bool semi = !d.is_insulator(ed.first) && !d.is_insulator(ed.second);
        p.edge_current_n.push_back(semi ? 1e-3 : 0.0);
        p.edge_current_p.push_back(0.0);
    }
    const auto e = *element_vectors(d, g, p, VectorField::current);
    for (std::size_t k = 0; k < g.element_count(); ++k) {
        bool oxide = false;
        for (const NodeId v : g.element(k)) oxide = oxide || d.is_insulator(v);
        CHECK(static_cast<bool>(e.defined[k]) == !oxide);
    }
    const auto mag = *field_values(d, g, p, Field::current_density);
    for (std::size_t i = 0; i < mag.size(); ++i) {
        CHECK(std::isnan(mag[i]) == d.is_insulator(static_cast<NodeId>(i)));
    }
    CHECK(field_values(d, g, p, Field::electric_field).has_value());

    // No edge currents (a quasi-static result), a graph geometry, another mesh: refused.
    results::BiasPoint q = p;
    q.edge_current_n.clear();
    const auto none = element_vectors(d, g, q, VectorField::electron_current);
    REQUIRE_FALSE(none.has_value());
    CHECK(none.error().message.find("no edge currents") != std::string::npos);
    CHECK(element_vectors(d, g, q, VectorField::hole_current).has_value());
    const Geometry graph = Geometry::graph(d.mesh());
    CHECK(element_vectors(d, graph, p, VectorField::electric).error().code == base::ErrorCode::invalid_input);
    const Built other = tensor(2);
    CHECK(element_vectors(other.device, g, p, VectorField::electric).error().code ==
          base::ErrorCode::invalid_input);
}

TEST_CASE("render vectors: a solved resistor's |J| is I / area at every node (BG1)") {
    // A uniform n-type bar, 2 um by 1 um, 0.1 V across: the current is uniform, so every node's |J|
    // is the terminal current over the bar's height (per unit depth), and J points from + to -.
    const auto x = graded(0.0, 2e-4, 9);
    const auto y = graded(0.0, 1e-4, 5);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    auto cells = *mesh::TensorCells::create(m, x, y);
    const device::Device d = make_device(std::move(m), std::move(cells), "x_min", "x_max");
    const std::vector<double> bias{0.1, 0.0};
    const auto r = solve::solve_bias(d, bias);
    REQUIRE(r.has_value());
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const double i_left = r->terminal_current[0];  // A/cm, into the device at the left
    const double expected = std::abs(i_left) / 1e-4;
    const double resolution = r->terminal_current_resolution[0] / 1e-4;
    CAPTURE(i_left, resolution);
    const auto e = *element_vectors(d, g, *r, VectorField::current);
    for (const Point& v : nodal_vectors(g, e)) {
        CHECK(std::abs(norm(v) - expected) <= 1e-6 * expected + resolution);
        CHECK(v[0] > 0.0);  // from the biased left contact to the right
    }
    const auto jn = *element_vectors(d, g, *r, VectorField::electron_current);
    const auto jp = *element_vectors(d, g, *r, VectorField::hole_current);
    for (std::size_t k = 0; k < g.element_count(); ++k) {
        for (std::size_t c = 0; c < 3; ++c) {
            CHECK(std::abs(jn.value[k][c] + jp.value[k][c] - e.value[k][c]) <= 1e-12 * expected);
        }
    }
}
