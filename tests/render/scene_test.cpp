// The CPU scenes (Unit 25a, VR1-VR4, VR8): region maps whose areas are the control volumes, fields
// at their exact colormap positions, contours on their levels, 3D solids with translucent
// insulators and the semiconductor seen through them, contact slabs, slices, mesh edges, field
// values from results and plots with fitted axes and flagged points.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/scene.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;

namespace {

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

// A MOS-like device in cm: silicon on 0 <= y <= 0.9 um under an oxide to 1 um, x over 2 um (and z
// over 0.6 um in 3D); n+ pockets (1e19) at the top corners, a p channel (1e17) between them, a
// lightly doped substrate (1e15) below (net 1.99e19 in the pockets). An electrode on the oxide's top over the channel, a source
// on the left pocket's x_min face, the substrate contact on y_min.
struct Mos {
    device::Device device;
    std::vector<double> x, y, z;
};

Mos mos(bool three_d, bool mirrored = false) {
    const auto x = uniform(0.0, 2e-4, 9);
    const auto y = uniform(0.0, 1e-4, 11);
    const auto z = three_d ? uniform(0.0, 0.6e-4, 4) : std::vector<double>{};
    mesh::Mesh m = three_d ? *mesh::make_tensor_grid(x, y, z) : *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        const double h = mirrored ? 1e-4 - p[1] : p[1];  // height above the substrate contact
        region[i] = h > 0.95e-4 ? 0 : 1;
        if (region[i] == 0) continue;
        const bool top = h > 0.55e-4;
        acceptors[i] = top ? 1e17 : 1e15;
        if (top && (p[0] < 0.45e-4 || p[0] > 1.55e-4)) donors[i] = 2e19;
    }
    std::vector<mesh::NodeId> gate, source;
    for (const mesh::NodeId v : m.find_boundary(mirrored ? "y_min" : "y_max")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        if (px > 0.6e-4 && px < 1.4e-4) gate.push_back(v);
    }
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        const double py = m.points()[static_cast<std::size_t>(v)][1];
        const double hy = mirrored ? 1e-4 - py : py;
        if (hy > 0.55e-4 && hy < 0.95e-4) source.push_back(v);
    }
    auto substrate = m.find_boundary(mirrored ? "y_max" : "y_min")->nodes;
    device::Contact g{"gate", device::ContactKind::electrode, std::move(gate)};
    g.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto cells = three_d ? *mesh::TensorCells::create(m, x, y, z) : *mesh::TensorCells::create(m, x, y);
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(g),
                      {"source", device::ContactKind::ohmic, std::move(source)},
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .cells = std::move(cells)});
    if (!d) FAIL(d.error().message);
    return {std::move(*d), x, y, z};
}

double triangles_area(const TriangleBatch& b) {
    double a = 0.0;
    for (std::size_t i = 0; i + 2 < b.vertices.size(); i += 3) {
        std::vector<Point> p;
        for (std::size_t k = 0; k < 3; ++k) {
            const Float3& f = b.vertices[i + k].position;
            p.push_back({f[0], f[1], f[2]});
        }
        a += polygon_area(p);
    }
    return a;
}

// Every triangle is counterclockwise seen from the side its normal points to.
void check_winding(const FieldScene& s) {
    for (const TriangleBatch& b : s.triangles) {
        for (std::size_t i = 0; i + 2 < b.vertices.size(); i += 3) {
            const Float3 a = b.vertices[i].position, c = b.vertices[i + 1].position,
                         d = b.vertices[i + 2].position;
            const std::array<double, 3> u{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            const std::array<double, 3> w{d[0] - a[0], d[1] - a[1], d[2] - a[2]};
            const std::array<double, 3> n{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                                          u[0] * w[1] - u[1] * w[0]};
            const Float3 m = b.vertices[i].normal;
            CAPTURE(b.name, i);
            CHECK(n[0] * m[0] + n[1] * m[1] + n[2] * m[2] > 0.0);
        }
    }
}

// The horizontal triangles at height y (um), not translucent: their normals' y components.
std::vector<float> facing_at(const FieldScene& s, float y) {
    std::vector<float> out;
    for (const TriangleBatch& b : s.triangles) {
        if (b.transparent || b.name == "electrode") continue;
        for (std::size_t i = 0; i + 2 < b.vertices.size(); i += 3) {
            bool flat = true;
            for (std::size_t k = 0; k < 3; ++k) flat = flat && std::abs(b.vertices[i + k].position[1] - y) < 1e-5f;
            if (flat) out.push_back(b.vertices[i].normal[1]);
        }
    }
    return out;
}

const TriangleBatch* find(const FieldScene& s, std::string_view name) {
    for (const TriangleBatch& b : s.triangles) {
        if (b.name == name) return &b;
    }
    return nullptr;
}

}  // namespace

TEST_CASE("render scene: a 2D region map is the control volumes by material") {
    const Mos m = mos(false);
    const auto g = Geometry::tensor(m.device.mesh(), *m.device.cells());
    REQUIRE(g.has_value());
    const auto s = solid_view(m.device, *g);
    REQUIRE(s.has_value());
    CHECK(s->dimension == 2);
    check_winding(*s);
    // Expected area per material (um^2): the control volumes of its nodes.
    std::array<double, 7> expected{};
    const auto volumes = m.device.mesh().volumes();
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        const Material c = material_class(m.device, static_cast<NodeId>(i));
        expected[static_cast<std::size_t>(c)] += volumes[i] * 1e8;
    }
    double total = 0.0;
    for (int c = 0; c < 6; ++c) {
        const Material mat = static_cast<Material>(c);
        const TriangleBatch* b = find(*s, material_name(mat));
        const double area = b ? triangles_area(*b) : 0.0;
        CAPTURE(material_name(mat), area, expected[static_cast<std::size_t>(c)]);
        CHECK(std::abs(area - expected[static_cast<std::size_t>(c)]) <= 1e-5);
        if (b) {
            CHECK_FALSE(b->lit);
            for (const SceneVertex& v : b->vertices) CHECK(v.color.a == 1.0f);
        }
        total += area;
    }
    CHECK(std::abs(total - 2.0) <= 1e-5);
    // The outline: the boundary and the oxide-silicon interface halfway between node rows.
    REQUIRE(!s->lines.empty());
    CHECK(s->lines.front().on_top);
    bool interface = false;
    for (std::size_t i = 0; i + 1 < s->lines.front().points.size(); i += 2) {
        const Float3 a = s->lines.front().points[i], b = s->lines.front().points[i + 1];
        interface = interface || (std::abs(a[1] - 0.95f) < 1e-5f && std::abs(b[1] - 0.95f) < 1e-5f);
    }
    CHECK(interface);
    // Contacts as slabs outside the domain; the camera frames everything.
    const TriangleBatch* metal = find(*s, "electrode");
    REQUIRE(metal != nullptr);
    float top = -1.0f;
    for (const SceneVertex& v : metal->vertices) top = std::max(top, v.position[1]);
    CHECK(top > 1.0f);
    CHECK(s->hi[1] >= top);
    CHECK(s->lo[0] < 0.0f);  // the source pad on x_min
    CHECK(s->camera.elevation_deg == 0.0f);  // 2D: looking down -z
}

TEST_CASE("render scene: a 3D solid with a translucent oxide shows the silicon under it") {
    const Mos m = mos(true);
    const auto g = Geometry::tensor(m.device.mesh(), *m.device.cells());
    const auto s = solid_view(m.device, *g);
    REQUIRE(s.has_value());
    CHECK(s->dimension == 3);
    CHECK(s->axis_triad);
    const TriangleBatch* ox = find(*s, "insulator");
    REQUIRE(ox != nullptr);
    CHECK(ox->transparent);
    CHECK(ox->lit);
    for (const SceneVertex& v : ox->vertices) CHECK(v.color.a < 1.0f);
    // The silicon's surface under the oxide (y = 0.95) faces the oxide: up here, and down in the
    // mirrored device (oxide below, at y < 0.05).
    const auto up = facing_at(*s, 0.95f);
    REQUIRE(!up.empty());
    for (const float ny : up) CHECK(ny > 0.99f);
    const Mos mm = mos(true, true);
    const auto gm = Geometry::tensor(mm.device.mesh(), *mm.device.cells());
    const auto sm = solid_view(mm.device, *gm);
    REQUIRE(sm.has_value());
    const auto down = facing_at(*sm, 0.05f);
    REQUIRE(!down.empty());
    for (const float ny : down) CHECK(ny < -0.99f);
    check_winding(*s);
    check_winding(*sm);
    // Boundary surface pieces of the materials cover the box's surface (um^2).
    double surface = 0.0;
    for (const TriangleBatch& b : s->triangles) {
        if (b.name == "electrode" || b.name == "gate oxide" || b.name == "contact") continue;
        surface += triangles_area(b);
    }
    const double box = 2.0 * (2.0 * 1.0 + 2.0 * 0.6 + 1.0 * 0.6);
    const double interface = 2.0 * 0.6;  // the silicon under the oxide, drawn once more
    CHECK(std::abs(surface - box - interface) <= 1e-4);
    // Every lit normal is unit length.
    for (const TriangleBatch& b : s->triangles) {
        for (const SceneVertex& v : b.vertices) {
            const float n = std::sqrt(v.normal[0] * v.normal[0] + v.normal[1] * v.normal[1] +
                                      v.normal[2] * v.normal[2]);
            CHECK(std::abs(n - 1.0f) < 1e-5f);
        }
    }
}

TEST_CASE("render scene: a 2D field map carries exact colormap positions and contours") {
    const Mos m = mos(false);
    const auto g = Geometry::tensor(m.device.mesh(), *m.device.cells());
    // A potential linear in x: 0 V to 1 V over 2 um.
    std::vector<double> psi;
    for (const Point& p : m.device.mesh().points()) psi.push_back(p[0] / 2e-4);
    FieldOptions o;
    o.contours = true;
    o.title = "psi (V)";
    o.scale = *ColorScale::create(ScaleKind::linear, -1.0, 3.0);
    const auto s = field_view(m.device, *g, psi, o);
    REQUIRE(s.has_value());
    REQUIRE(s->colorbar.has_value());
    CHECK(s->colorbar->title == "psi (V)");
    const TriangleBatch* f = find(*s, "field");
    REQUIRE(f != nullptr);
    CHECK(f->shading == Shading::field);
    CHECK(f->vertices.size() == 3 * g->element_count());
    for (const SceneVertex& v : f->vertices) {
        // position of psi = x / 2 um on [-1, 3]
        CHECK(std::abs(v.value - (v.position[0] / 2.0f + 1.0f) / 4.0f) <= 1e-6f);
    }
    check_winding(*s);
    // Contours at the interior ticks -0.5, 0, ..., 2.5 that the field reaches: psi = 0.5 at
    // x = 1 um, and psi = 1 along the domain's edge x = 2 um (a level at a vertex counts as
    // above it); psi = 0 at x = 0 draws nothing, every triangle there being at or above it.
    const LineBatch* c = nullptr;
    for (const LineBatch& l : s->lines) {
        if (l.name == "contours") c = &l;
    }
    REQUIRE(c != nullptr);
    std::set<int> xs;
    for (const Float3& p : c->points) {
        const float x = p[0];
        CHECK(std::abs(x - std::round(x)) <= 1e-5f);
        xs.insert(static_cast<int>(std::lround(x)));
    }
    CHECK(xs == std::set<int>{1, 2});
    // A fitted scale, and the errors.
    const auto fitted = field_view(m.device, *g, psi);
    REQUIRE(fitted.has_value());
    CHECK(fitted->colorbar->scale.lo() == 0.0);
    CHECK(fitted->colorbar->scale.hi() == 1.0);
    CHECK(field_view(m.device, *g, std::vector<double>(3, 0.0)).error().code ==
          base::ErrorCode::invalid_input);
    const std::vector<double> nothing(m.device.mesh().node_count(), std::nan(""));
    CHECK(field_view(m.device, *g, nothing).error().code == base::ErrorCode::invalid_input);
}

TEST_CASE("render scene: a 3D field shows on the silicon; the oxide stays translucent") {
    const Mos m = mos(true);
    const auto g = Geometry::tensor(m.device.mesh(), *m.device.cells());
    std::vector<double> v;
    for (const Point& p : m.device.mesh().points()) v.push_back(p[1] / 1e-4);
    const auto s = field_view(m.device, *g, v);
    REQUIRE(s.has_value());
    const TriangleBatch* f = find(*s, "field");
    REQUIRE(f != nullptr);
    CHECK(f->lit);
    check_winding(*s);
    for (const SceneVertex& x : f->vertices) {
        CHECK(x.position[1] <= 0.95f + 1e-5f);  // only the silicon's pieces and its interface
        CHECK(std::abs(x.value - x.position[1]) <= 2e-6f);  // y in um equals the value
    }
    REQUIRE(find(*s, "insulator") != nullptr);
    CHECK(find(*s, "insulator")->transparent);
}

TEST_CASE("render scene: slices and the mesh view") {
    const Mos m = mos(true);
    const auto g = Geometry::tensor(m.device.mesh(), *m.device.cells());
    std::vector<double> v;
    for (const Point& p : m.device.mesh().points()) v.push_back(1.0 + p[0] / 2e-4);
    const std::vector<Plane> planes{{{0.0, 0.0, 1.0}, 0.3e-4}, {{1.0, 0.0, 0.0}, 1.0e-4}};
    const auto s = slice_view(m.device, *g, v, planes);
    REQUIRE(s.has_value());
    const TriangleBatch* slice = find(*s, "slice");
    REQUIRE(slice != nullptr);
    // The fitted scale is [1, 2], so the position of 1 + x / 2 um is x / 2 um.
    for (const SceneVertex& x : slice->vertices) CHECK(std::abs(x.value - x.position[0] / 2.0f) <= 1e-6f);
    check_winding(*s);
    // The two sections: 2 x 1 um and 1 x 0.6 um.
    CHECK(std::abs(triangles_area(*slice) - (2.0 + 0.6)) <= 1e-4);
    const TriangleBatch* shell = find(*s, "outline");
    REQUIRE(shell != nullptr);
    CHECK(shell->transparent);
    const std::vector<Plane> miss{{{0.0, 0.0, 1.0}, 5e-4}};
    CHECK(slice_view(m.device, *g, v, miss).error().code == base::ErrorCode::invalid_input);
    const Mos flat = mos(false);
    const auto g2 = Geometry::tensor(flat.device.mesh(), *flat.device.cells());
    std::vector<double> v2(flat.device.mesh().node_count(), 1.0);
    CHECK(slice_view(flat.device, *g2, v2, planes).error().code == base::ErrorCode::invalid_input);

    // 2D mesh view: every element edge once (8 x 10 cells: 9*10 + 8*11 grid edges + 80 diagonals).
    const auto mv = mesh_view(flat.device, *g2);
    REQUIRE(mv.has_value());
    const LineBatch* edges = nullptr;
    for (const LineBatch& l : mv->lines) {
        if (l.name == "mesh") edges = &l;
    }
    REQUIRE(edges != nullptr);
    CHECK(edges->points.size() == 2 * (9 * 10 + 8 * 11 + 80));
}

TEST_CASE("render scene: views refuse another mesh and 1D geometries") {
    const Mos a = mos(false);
    const Mos b = mos(true);
    const auto g = Geometry::tensor(b.device.mesh(), *b.device.cells());
    CHECK(solid_view(a.device, *g).error().code == base::ErrorCode::invalid_input);
    const auto line = *mesh::make_tensor_grid(uniform(0.0, 1e-4, 5));
    const Geometry g1 = Geometry::graph(line);
    CHECK(solid_view(a.device, g1).error().code == base::ErrorCode::invalid_input);
}

TEST_CASE("render scene: field values from a bias point, and material classes") {
    const Mos m = mos(false);
    const std::size_t n = m.device.mesh().node_count();
    results::BiasPoint p;
    p.fields.potential_V.assign(n, 0.25);
    p.fields.n_cm3.assign(n, 1e10);
    p.fields.p_cm3.assign(n, 1e17);
    p.bands.conduction_eV.assign(n, 0.5);
    p.bands.valence_eV.assign(n, -0.6);
    p.bands.electron_fermi_eV.assign(n, 0.0);
    p.bands.hole_fermi_eV.assign(n, 0.0);
    CHECK(field_values(m.device, p, Field::potential)->at(3) == 0.25);
    CHECK(field_values(m.device, p, Field::valence_band)->at(0) == -0.6);
    CHECK(field_values(m.device, p, Field::temperature)->at(5) == 300.0);
    // Unit 52: the carrier temperatures, an error that says so without energy balance.
    CHECK(field_values(m.device, p, Field::electron_temperature).error().message.find("energy_balance") !=
          std::string::npos);
    CHECK_FALSE(field_values(m.device, p, Field::hole_temperature).has_value());
    p.fields.electron_temperature_K.assign(n, 300.0);
    p.fields.electron_temperature_K[4] = 812.5;
    p.fields.hole_temperature_K.assign(n, 301.0);
    CHECK(field_values(m.device, p, Field::electron_temperature)->at(4) == 812.5);
    CHECK(field_values(m.device, p, Field::hole_temperature)->at(2) == 301.0);
    CHECK(find_field("Tn") == Field::electron_temperature);
    CHECK(find_field("Tp") == Field::hole_temperature);
    const auto net = *field_values(m.device, p, Field::net_doping);
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(net[i] == m.device.donors()[i] - m.device.acceptors()[i]);
    }
    p.fields.n_cm3.pop_back();
    CHECK(field_values(m.device, p, Field::electrons).error().code == base::ErrorCode::invalid_input);
    for (int f = 0; f <= static_cast<int>(Field::last); ++f) {
        CHECK(find_field(field_info(static_cast<Field>(f)).key) == static_cast<Field>(f));
    }
    CHECK_FALSE(find_field("Jz").has_value());
    // The vector magnitudes need the geometry (Unit 25b).
    CHECK(field_values(m.device, p, Field::electric_field).error().message.find("needs the geometry") !=
          std::string::npos);

    // Unit 45: the rates, generations as positive; a missing rate is an error that says so.
    const auto rf = [&]() {
        p.fields.n_cm3.assign(n, 1e10);
        return p;
    }();
    results::BiasPoint q = rf;
    q.rates.net_n.assign(n, 2e15);
    q.rates.srh.assign(n, 1.5e15);
    q.rates.impact.assign(n, -3e12);
    CHECK(field_values(m.device, q, Field::recombination)->at(4) == 2e15);
    CHECK(field_values(m.device, q, Field::srh)->at(4) == 1.5e15);
    CHECK(field_values(m.device, q, Field::impact_generation)->at(4) == 3e12);
    const auto missing = field_values(m.device, q, Field::btbt_generation);
    REQUIRE_FALSE(missing.has_value());
    CHECK(missing.error().message.find("holds no Gbtbt rates") != std::string::npos);
    CHECK(field_info(Field::recombination).default_scale == ScaleKind::signed_log10);
    CHECK(field_info(Field::impact_generation).default_scale == ScaleKind::log10);

    std::set<Material> seen;
    for (std::size_t i = 0; i < n; ++i) seen.insert(material_class(m.device, static_cast<NodeId>(i)));
    CHECK(seen == std::set<Material>{Material::insulator, Material::n_heavy, Material::p_channel,
                                     Material::substrate});
    const Palette strict{.heavy_cm3 = 1e20, .channel_cm3 = 1e16, .insulator_alpha = 0.35f};
    for (std::size_t i = 0; i < n; ++i) {
        // With a 1e20 heavy threshold the pockets (net 1.99e19) are n, the channel (1e17) still p.
        const Material c = material_class(m.device, static_cast<NodeId>(i), strict);
        const Material d = material_class(m.device, static_cast<NodeId>(i));
        if (d == Material::n_heavy) CHECK(c == Material::n_light);
        if (d == Material::p_channel) CHECK(c == Material::p_channel);
    }
}

TEST_CASE("render scene: plots fit their axes and flag points below resolution") {
    analysis::Curve iv{{0.0, 0.1, 0.2, 0.3}, {-1e-15, -2e-12, -3e-9, -4e-6}, {1e-14, 1e-14, 1e-14, 1e-14}};
    const std::vector<SeriesInput> in{{"I_D", iv}};
    const auto p = curve_plot(in, "V (V)", "|I| (A)", AxisScale::linear, AxisScale::log10);
    REQUIRE(p.has_value());
    CHECK(p->y.scale == AxisScale::log10);
    CHECK(p->y.lo == 1e-15);
    CHECK(p->y.hi == 1e-5);
    CHECK(p->series[0].y[3] == 4e-6);  // all negative: plotted as |y|
    CHECK(p->series[0].flagged == std::vector<std::uint8_t>{1, 0, 0, 0});
    CHECK_FALSE(p->legend);
    for (const Tick& t : p->x.ticks) {
        CHECK(t.value >= p->x.lo);
        CHECK(t.value <= p->x.hi);
    }
    CHECK(p->x.lo < 0.0);  // padded
    const std::vector<SeriesInput> none;
    CHECK(curve_plot(none, "x", "y").error().code == base::ErrorCode::invalid_input);
    const std::vector<SeriesInput> zero{{"z", {{0.0, 1.0}, {0.0, 0.0}, {}}}};
    CHECK(curve_plot(zero, "x", "y", AxisScale::linear, AxisScale::log10).error().code ==
          base::ErrorCode::invalid_input);

    // A 1D band diagram: four series, the insulator's NaN kept as a gap.
    const auto line = *mesh::make_tensor_grid(uniform(0.0, 1e-4, 5));
    results::BandDiagram b{{0.5, 0.5, std::nan(""), 0.4, 0.3},
                           {-0.6, -0.6, std::nan(""), -0.7, -0.8},
                           {0.0, 0.0, std::nan(""), 0.0, 0.0},
                           {0.0, 0.0, std::nan(""), -0.1, -0.1}};
    const auto bp = band_plot(line, b);
    REQUIRE(bp.has_value());
    REQUIRE(bp->series.size() == 4);
    CHECK(std::isnan(bp->series[1].y[2]));
    CHECK(bp->series[0].x[4] == 1.0);  // um
    CHECK(bp->y.lo < -0.8);
    CHECK(bp->y.hi > 0.5);
    b.valence_eV.pop_back();
    CHECK(band_plot(line, b).error().code == base::ErrorCode::invalid_input);
    const Mos m = mos(false);
    CHECK(profile_plot(m.device.mesh(), std::vector<double>(m.device.mesh().node_count(), 1.0), "y")
              .error()
              .code == base::ErrorCode::invalid_input);
    const auto prof = profile_plot(line, std::vector<double>{1e15, 1e16, 1e17, 1e18, 1e19}, "n",
                                   AxisScale::log10);
    REQUIRE(prof.has_value());
    CHECK(prof->y.lo == 1e15);
    CHECK(prof->y.hi == 1e19);
}

TEST_CASE("render scene: isosurfaces on the field's scale, contours at levels in 2D (Unit 25b, BQ3)") {
    // A field exponential in y, 1e10 at y = 0 to 1e20 at y = 1 um: linear on the log scale, so the
    // 10^14.3 isosurface is exactly the plane y = 0.43 um there (linear interpolation in the value would
    // put it near the upper node of each cell).
    const Mos m3 = mos(true);
    const Geometry g3 = *Geometry::tensor(m3.device.mesh(), *m3.device.cells());
    std::vector<double> f;
    for (const Point& x : m3.device.mesh().points()) f.push_back(std::pow(10.0, 10.0 + 10.0 * x[1] / 1e-4));
    FieldOptions o;
    o.scale = *ColorScale::create(ScaleKind::log10, 1e10, 1e20);
    o.view.contacts = false;
    const std::vector<double> one{std::pow(10.0, 14.3)};  // y = 0.43 um, between grid planes
    const auto s = iso_view(m3.device, g3, f, one, o);
    REQUIRE(s.has_value());
    const auto iso = std::ranges::find_if(s->triangles, [](const TriangleBatch& b) { return b.name == "isosurface"; });
    REQUIRE(iso != s->triangles.end());
    CHECK_FALSE(iso->transparent);
    for (const SceneVertex& v : iso->vertices) {
        CHECK(std::abs(v.position[1] - 0.43f) <= 1e-5f);  // um
        CHECK(std::abs(v.value - 0.43f) <= 1e-6f);         // the level's colormap position
        CHECK(v.normal[1] > 0.999f);                      // toward higher values
    }
    // Two levels: two translucent surfaces.
    const std::vector<double> two{1e12, 1e18};
    const auto t = iso_view(m3.device, g3, f, two, o);
    REQUIRE(t.has_value());
    CHECK(std::ranges::count_if(t->triangles, [](const TriangleBatch& b) {
              return b.name == "isosurface" && b.transparent;
          }) == 2);
    // Refused: no level, a level off the log scale.
    CHECK(iso_view(m3.device, g3, f, std::vector<double>{}, o).error().code == base::ErrorCode::invalid_input);
    CHECK(iso_view(m3.device, g3, f, std::vector<double>{-1.0}, o).error().code == base::ErrorCode::invalid_input);

    // 2D: the field map with contour lines at the levels only, on the log scale (y = 0.43 um).
    const Mos m2 = mos(false);
    const Geometry g2 = *Geometry::tensor(m2.device.mesh(), *m2.device.cells());
    std::vector<double> f2;
    for (const Point& x : m2.device.mesh().points()) f2.push_back(std::pow(10.0, 10.0 + 10.0 * x[1] / 1e-4));
    const auto c = iso_view(m2.device, g2, f2, one, o);
    REQUIRE(c.has_value());
    const auto lines = std::ranges::find_if(c->lines, [](const LineBatch& b) { return b.name == "contours"; });
    REQUIRE(lines != c->lines.end());
    CHECK_FALSE(lines->points.empty());
    for (const Float3& p : lines->points) CHECK(std::abs(p[1] - 0.43f) <= 1e-5f);
}
