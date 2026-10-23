// Region visibility (Unit 26c, YD4, YG2): with no region hidden every scene is as before; a hidden
// region's pieces, elements and contacts are absent and the shown regions' unchanged; the faces a hidden
// region shared with a shown one close the shown regions; scales are fitted over the shown nodes.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/quality.hpp"
#include "NiTCAD/render/scene.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;

namespace {

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

// Silicon (region 1) up to y = 0.8 um, oxide (region 0) above, so the interface is the plane y = 0.85 um
// between the node rows; a metal gate on the oxide's top over the middle, ohmic source and drain on the
// silicon's x faces, the body on y_min. 2D, or 3D 0.6 um deep.
device::Device mos(bool three) {
    const auto x = uniform(0.0, 2e-4, 9), y = uniform(0.0, 1e-4, 11);
    const auto z = three ? uniform(0.0, 0.6e-4, 4) : std::vector<double>{};
    mesh::Mesh m = three ? *mesh::make_tensor_grid(x, y, z) : *mesh::make_tensor_grid(x, y);
    auto cells = three ? *mesh::TensorCells::create(m, x, y, z) : *mesh::TensorCells::create(m, x, y);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 1);
    std::vector<double> donors(n, 0.0), acceptors(n, 1e17);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if (p[1] > 0.85e-4) region[i] = 0, acceptors[i] = 0.0;
        else if (p[1] > 0.5e-4 && (p[0] < 0.6e-4 || p[0] > 1.4e-4)) donors[i] = 2e19;
    }
    const auto on = [&](const char* patch, double x_lo, double x_hi, double y_lo, double y_hi) {
        std::vector<mesh::NodeId> v;
        for (const mesh::NodeId k : m.find_boundary(patch)->nodes) {
            const auto& p = m.points()[static_cast<std::size_t>(k)];
            if (p[0] >= x_lo && p[0] <= x_hi && p[1] >= y_lo && p[1] <= y_hi) v.push_back(k);
        }
        return v;
    };
    device::Contact gate{"gate", device::ContactKind::electrode, on("y_max", 0.6e-4, 1.4e-4, 0.0, 1.0)};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.1};
    auto source = on("x_min", -1.0, 1.0, 0.55e-4, 0.8e-4);
    auto drain = on("x_max", -1.0, 1.0, 0.55e-4, 0.8e-4);
    auto body = m.find_boundary("y_min")->nodes;
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
                                     .node_region = std::move(region),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {std::move(gate),
                                                  {"source", device::ContactKind::ohmic, std::move(source)},
                                                  {"drain", device::ContactKind::ohmic, std::move(drain)},
                                                  {"body", device::ContactKind::ohmic, std::move(body)}},
                                     .cells = std::move(cells)});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

// A field with every value in the oxide far above the silicon's: x in the silicon (cm), 1 in the oxide.
std::vector<double> field_of(const device::Device& d) {
    std::vector<double> v;
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        v.push_back(d.node_region()[i] == 0 ? 1.0 : d.mesh().points()[i][0]);
    }
    return v;
}

using P = std::array<double, 3>;
P at(const SceneVertex& v) { return {v.position[0], v.position[1], v.position[2]}; }
P sub(const P& a, const P& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
P cross(const P& a, const P& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
double dot(const P& a, const P& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

bool is_contact(const TriangleBatch& b) { return b.name == "electrode" || b.name == "gate oxide"; }

// The triangles' area (2D, um^2), or (3D) the volume they enclose and their vector area.
double area(const FieldScene& s) {
    double a = 0.0;
    for (const TriangleBatch& b : s.triangles) {
        if (is_contact(b)) continue;
        for (std::size_t i = 0; i + 2 < b.vertices.size(); i += 3) {
            const P n = cross(sub(at(b.vertices[i + 1]), at(b.vertices[i])), sub(at(b.vertices[i + 2]), at(b.vertices[i])));
            a += 0.5 * std::sqrt(dot(n, n));
        }
    }
    return a;
}
struct Closure {
    double volume = 0.0, area = 0.0;
    P vector_area{};
};
Closure closure(const FieldScene& s) {
    Closure c;
    for (const TriangleBatch& b : s.triangles) {
        if (is_contact(b)) continue;
        for (std::size_t i = 0; i + 2 < b.vertices.size(); i += 3) {
            const P p0 = at(b.vertices[i]), p1 = at(b.vertices[i + 1]), p2 = at(b.vertices[i + 2]);
            c.volume += dot(p0, cross(p1, p2)) / 6.0;
            const P n = cross(sub(p1, p0), sub(p2, p0));
            for (std::size_t k = 0; k < 3; ++k) c.vector_area[k] += 0.5 * n[k];
            c.area += 0.5 * std::sqrt(dot(n, n));
        }
    }
    return c;
}

// The control volumes of a region's nodes, um^2 (2D) or um^3 (3D).
double region_measure(const device::Device& d, std::int32_t r) {
    double v = 0.0;
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        if (d.node_region()[i] == r) v += d.mesh().volumes()[i];
    }
    return v * (d.mesh().dimension() == 2 ? 1e8 : 1e12);
}

const TriangleBatch* batch(const FieldScene& s, const std::string& name) {
    for (const TriangleBatch& b : s.triangles) {
        if (b.name == name) return &b;
    }
    return nullptr;
}

bool same_vertices(const TriangleBatch& a, const TriangleBatch& b) {
    if (a.vertices.size() != b.vertices.size()) return false;
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        const SceneVertex& u = a.vertices[i];
        const SceneVertex& v = b.vertices[i];
        if (u.position != v.position || u.normal != v.normal || u.color.r != v.color.r || u.color.a != v.color.a) return false;
        if (!(u.value == v.value || (std::isnan(u.value) && std::isnan(v.value)))) return false;
    }
    return true;
}

bool same_scene(const FieldScene& a, const FieldScene& b) {
    if (a.triangles.size() != b.triangles.size() || a.lines.size() != b.lines.size()) return false;
    for (std::size_t k = 0; k < a.triangles.size(); ++k) {
        if (a.triangles[k].name != b.triangles[k].name || !same_vertices(a.triangles[k], b.triangles[k])) return false;
    }
    for (std::size_t k = 0; k < a.lines.size(); ++k) {
        if (a.lines[k].points != b.lines[k].points) return false;
    }
    return a.lo == b.lo && a.hi == b.hi;
}

// The largest y (um) of every triangle and line point but the contacts'.
float top(const FieldScene& s, bool lines_too = true) {
    float y = -std::numeric_limits<float>::infinity();
    for (const TriangleBatch& b : s.triangles) {
        if (is_contact(b)) continue;
        for (const SceneVertex& v : b.vertices) y = std::max(y, v.position[1]);
    }
    if (lines_too) {
        for (const LineBatch& l : s.lines) for (const Float3& p : l.points) y = std::max(y, p[1]);
    }
    return y;
}

}  // namespace

TEST_CASE("render hidden regions: nothing hidden draws as before; a region without nodes hides nothing (YG2)") {
    for (const bool three : {false, true}) {
        CAPTURE(three);
        const device::Device d = mos(three);
        const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
        const std::vector<double> f = field_of(d);
        FieldOptions o;
        o.contours = true;
        o.view.mesh_edges = true;
        // Hiding nothing is the default; an empty list and a list of no region with nodes draw the same.
        const auto plain = field_view(d, g, f, o);
        REQUIRE(plain.has_value());
        o.view.hidden_regions = {};
        CHECK(same_scene(*plain, *field_view(d, g, f, o)));
        CHECK(same_scene(*solid_view(d, g), *solid_view(d, g, o.view)));
        CHECK(same_scene(*mesh_view(d, g), *mesh_view(d, g, o.view)));
        // An id that is not a region's is refused.
        o.view.hidden_regions = {2};
        CHECK(field_view(d, g, f, o).error().code == base::ErrorCode::invalid_input);
        CHECK(solid_view(d, g, o.view).error().code == base::ErrorCode::invalid_input);
        o.view.hidden_regions = {-1};
        CHECK_FALSE(mesh_view(d, g, o.view).has_value());
    }
}

TEST_CASE("render hidden regions: 2D views leave out the hidden region's pieces and contacts (YG2)") {
    const device::Device d = mos(false);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    ViewOptions all, no_oxide, no_silicon;
    no_oxide.hidden_regions = {0};
    no_silicon.hidden_regions = {1};
    const auto full = solid_view(d, g, all);
    const auto silicon = solid_view(d, g, no_oxide);
    const auto oxide = solid_view(d, g, no_silicon);
    REQUIRE((full && silicon && oxide));
    // The pieces drawn are the shown nodes' control volumes, exactly.
    CHECK(std::abs(area(*full) - 2.0) <= 1e-5);
    CHECK(std::abs(area(*silicon) - region_measure(d, 1)) <= 1e-5);
    CHECK(std::abs(area(*silicon) - 1.7) <= 1e-5);
    CHECK(std::abs(area(*oxide) - region_measure(d, 0)) <= 1e-5);
    // The silicon's batches as drawn with every region; no insulator.
    for (const char* name : {"n+", "p"}) {
        CAPTURE(name);
        REQUIRE(batch(*full, name) != nullptr);
        REQUIRE(batch(*silicon, name) != nullptr);
        CHECK(same_vertices(*batch(*full, name), *batch(*silicon, name)));
        CHECK(batch(*oxide, name) == nullptr);
    }
    CHECK(batch(*silicon, "insulator") == nullptr);
    CHECK(same_vertices(*batch(*full, "insulator"), *batch(*oxide, "insulator")));
    // The gate (on oxide nodes only) is left out with the oxide; source, drain and body stay.
    const auto slab_top = [](const FieldScene& s) {
        float y = -1.0f;
        for (const SceneVertex& v : batch(s, "electrode")->vertices) y = std::max(y, v.position[1]);
        return y;
    };
    CHECK(slab_top(*full) > 1.0f);
    CHECK(slab_top(*silicon) < 0.85f);
    CHECK(batch(*oxide, "electrode") != nullptr);  // the gate
    // Outlines: the shown halves of the boundary and the interface, nothing above it; the interface (y = 0.85 um)
    // drawn across the device, closing the silicon.
    CHECK(top(*silicon) <= 0.85f + 1e-5f);
    CHECK(top(*full) >= 1.0f - 1e-5f);
    const auto across = [](const FieldScene& s) {
        double length = 0.0;
        const auto& p = s.lines.front().points;
        for (std::size_t i = 0; i + 1 < p.size(); i += 2) {
            if (std::abs(p[i][1] - 0.85f) < 1e-5f && std::abs(p[i + 1][1] - 0.85f) < 1e-5f) length += std::abs(p[i + 1][0] - p[i][0]);
        }
        return length;
    };
    CHECK(std::abs(across(*silicon) - 2.0) <= 1e-5);
    CHECK(std::abs(across(*oxide) - 2.0) <= 1e-5);
    // The mesh's edges: the shown nodes' halves.
    const auto mesh = mesh_view(d, g, no_oxide);
    REQUIRE(mesh.has_value());
    CHECK(top(*mesh) <= 0.85f + 1e-5f);
    CHECK(mesh->lines.back().points.size() < mesh_view(d, g)->lines.back().points.size());
    // The field map: the silicon's pieces, its scale fitted over the silicon.
    const std::vector<double> f = field_of(d);
    FieldOptions o;
    o.view = no_oxide;
    o.contours = true;
    o.view.mesh_edges = true;
    const auto map = field_view(d, g, f, o);
    REQUIRE(map.has_value());
    CHECK(std::abs(area(*map) - 1.7) <= 1e-5);
    CHECK(top(*map) <= 0.85f + 1e-5f);
    REQUIRE(map->colorbar.has_value());
    CHECK(map->colorbar->scale.hi() < 1e-3);  // the oxide's 1 left out of the fit
    o.view.hidden_regions.clear();
    CHECK(field_view(d, g, f, o)->colorbar->scale.hi() >= 1.0);
    // A given scale is kept.
    o.view = no_oxide;
    o.scale = *ColorScale::create(ScaleKind::linear, 0.0, 5.0);
    CHECK(field_view(d, g, f, o)->colorbar->scale.hi() == 5.0);
    // A filmstrip: one scale over every frame's shown nodes.
    std::vector<FilmFrame> frames{{"a", f}, {"b", f}};
    for (double& v : frames[1].values) v *= 0.5;
    FieldOptions film;
    film.view = no_oxide;
    const auto strip = filmstrip_scenes(d, g, frames, ScaleKind::linear, film);
    REQUIRE(strip.has_value());
    for (const FieldScene& s : *strip) {
        CHECK(s.colorbar->scale.hi() < 1e-3);
        CHECK(std::abs(area(s) - 1.7) <= 1e-5);
    }
    // Mesh quality: the elements' parts in the silicon.
    const auto angles = *element_quality(g);
    FieldOptions q;
    q.view = no_oxide;
    const auto quality = quality_view(d, g, angles, {}, q);
    REQUIRE(quality.has_value());
    CHECK(std::abs(area(*quality) - 1.7) <= 1e-5);
    CHECK(top(*quality) <= 0.85f + 1e-5f);
}

TEST_CASE("render hidden regions: 3D surfaces close the shown regions with the faces they shared (YG2)") {
    const device::Device d = mos(true);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const std::vector<double> f = field_of(d);
    for (const std::int32_t hidden : {0, 1}) {
        CAPTURE(hidden);
        const std::int32_t kept = 1 - hidden;
        const double expected = region_measure(d, kept);
        ViewOptions v;
        v.hidden_regions = {hidden};
        // The solid view (the oxide translucent, the silicon opaque) and the field on the surface.
        FieldOptions o;
        o.view = v;
        for (const auto& s : {solid_view(d, g, v), field_view(d, g, f, o)}) {
            REQUIRE(s.has_value());
            const Closure c = closure(*s);
            CHECK(std::abs(c.volume - expected) <= 1e-5 * expected);  // closed and outward
            for (std::size_t k = 0; k < 3; ++k) CHECK(std::abs(c.vector_area[k]) <= 1e-5 * c.area);
        }
        // Opaque insulators: the same surface.
        v.transparent_insulators = false;
        const Closure opaque = closure(*solid_view(d, g, v));
        CHECK(std::abs(opaque.volume - expected) <= 1e-5 * expected);
    }
    // The oxide hidden: nothing above the interface, the gate gone, the silicon's boundary pieces as before.
    ViewOptions v;
    v.hidden_regions = {0};
    const auto silicon = solid_view(d, g, v);
    const auto full = solid_view(d, g);
    REQUIRE((silicon && full));
    CHECK(top(*silicon) <= 0.85f + 1e-5f);
    CHECK(batch(*silicon, "insulator") == nullptr);
    float slab = -1.0f;
    for (const SceneVertex& x : batch(*silicon, "electrode")->vertices) slab = std::max(slab, x.position[1]);
    CHECK(slab < 0.85f);
    // Clipped: what is left of the silicon, closed.
    v.clip = Box{{1.0e-4, 0.3e-4, 0.3e-4}, {3e-4, 3e-4, 3e-4}};
    const Closure clipped = closure(*solid_view(d, g, v));
    const double box = 1.0 * (0.85 - 0.3) * 0.3;
    CHECK(std::abs(clipped.volume - (region_measure(d, 1) - box)) <= 1e-5 * region_measure(d, 1));
    v.clip.reset();
    // Mesh edges: the silicon's halves.
    CHECK(top(*mesh_view(d, g, v)) <= 0.85f + 1e-5f);
    // Slices and isosurfaces: the silicon's parts only, in the silicon's outline.
    FieldOptions o;
    o.view = v;
    const std::vector<Plane> planes{{{0.0, 0.0, 1.0}, 0.3e-4}, {{1.0, 0.0, 0.0}, 1.0e-4}};
    const auto slices = slice_view(d, g, f, planes, o);
    REQUIRE(slices.has_value());
    CHECK(top(*slices) <= 0.85f + 1e-5f);
    FieldScene sections;
    sections.triangles.push_back(*batch(*slices, "slice"));
    CHECK(std::abs(area(sections) - (2.0 * 0.85 + 0.85 * 0.6)) <= 1e-5);
    // The outline around them: the silicon's shell, closed by the faces the oxide opened.
    FieldScene shell;
    shell.triangles.push_back(*batch(*slices, "outline"));
    const Closure sc = closure(shell);
    CHECK(std::abs(std::abs(sc.volume) - region_measure(d, 1)) <= 1e-5 * region_measure(d, 1));
    // A curved isosurface (y + 0.3 x^2 = 0.9 um) that crosses the interface: only its part in the silicon.
    std::vector<double> y;
    for (const Point& p : d.mesh().points()) y.push_back(p[1] * 1e4 + 0.3 * (p[0] * 1e4) * (p[0] * 1e4));
    o.scale = *ColorScale::create(ScaleKind::linear, 0.0, 1.5);
    const std::vector<double> levels{0.9};
    const auto iso = iso_view(d, g, y, levels, o);
    o.view.hidden_regions.clear();
    const auto iso_all = iso_view(d, g, y, levels, o);
    REQUIRE((iso && iso_all));
    REQUIRE(batch(*iso, "isosurface") != nullptr);
    float iso_top = -1.0f, all_top = -1.0f;
    for (const SceneVertex& x : batch(*iso, "isosurface")->vertices) iso_top = std::max(iso_top, x.position[1]);
    for (const SceneVertex& x : batch(*iso_all, "isosurface")->vertices) all_top = std::max(all_top, x.position[1]);
    CHECK(iso_top <= 0.85f + 1e-5f);
    CHECK(all_top > 0.88f);
    // The pieces' normals interpolated from the triangle's: where a piece's corner is one of the surface's, its
    // normal is that corner's (the surface's vertices are not shared between elements, so a position can hold
    // a corner of each, with its own normal).
    std::size_t shared = 0, matched = 0;
    for (const SceneVertex& x : batch(*iso, "isosurface")->vertices) {
        bool at = false, same = false;
        for (const SceneVertex& a : batch(*iso_all, "isosurface")->vertices) {
            if (std::abs(x.position[0] - a.position[0]) + std::abs(x.position[1] - a.position[1]) +
                    std::abs(x.position[2] - a.position[2]) > 1e-6f) {
                continue;
            }
            at = true;
            same = same || std::abs(x.normal[0] - a.normal[0]) + std::abs(x.normal[1] - a.normal[1]) +
                                   std::abs(x.normal[2] - a.normal[2]) <= 1e-4f;
        }
        shared += at ? 1 : 0;
        matched += same ? 1 : 0;
    }
    CHECK(shared > 0);
    CHECK(matched == shared);
    // Quality on the surface: closed over the silicon.
    FieldOptions q;
    q.view = v;
    const Closure qc = closure(*quality_view(d, g, *element_quality(g), {}, q));
    CHECK(std::abs(qc.volume - region_measure(d, 1)) <= 1e-5 * region_measure(d, 1));
    // The exploded view: the silicon where it was, the oxide and the gate gone.
    const auto exploded = exploded_view(d, g, 1, 0.2e-4, {}, q);
    const auto exploded_all = exploded_view(d, g, 1, 0.2e-4);
    REQUIRE((exploded && exploded_all));
    CHECK(batch(*exploded, "insulator") == nullptr);
    for (const char* name : {"n+", "p"}) CHECK(same_vertices(*batch(*exploded, name), *batch(*exploded_all, name)));
    CHECK(batch(*exploded, "electrode")->vertices.size() < batch(*exploded_all, "electrode")->vertices.size());
}

TEST_CASE("render hidden regions: shown_values puts NaN on the hidden nodes (YG2)") {
    const device::Device d = mos(false);
    const std::vector<double> f = field_of(d);
    const std::vector<std::int32_t> oxide{0};
    const auto v = shown_values(d, f, oxide);
    REQUIRE(v.size() == f.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (d.node_region()[i] == 0) CHECK(std::isnan(v[i]));
        else CHECK(v[i] == f[i]);
    }
    CHECK(shown_values(d, f, {}) == f);
}
