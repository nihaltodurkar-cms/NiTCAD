// The band surface (Unit 26o; DECISIONS.md BSD1-BSD4, BSC2, BSC3, BSC6, BSC8): Ec, Ev, EFn and EFp lifted as
// energy-height surfaces over a 2D device or a plane of a 3D one, steps kept at material boundaries.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/band_surface.hpp"
#include "NiTCAD/render/cut_line.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/results/solution.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;

namespace {

constexpr double um = 1e-4;  // cm

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

physics::Semiconductor silicon_with(double d_chi) {
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.electron_affinity_eV += d_chi;
    return *physics::Semiconductor::create(p);
}

// A tensor device 2 um (x) by 1 um (y) (and 0.6 um in z for 3D), 9 x 11 (x 4) nodes, regions chosen per node; an
// ohmic contact on y_min.
struct Dev {
    device::Device device;
};

Dev make_dev(std::vector<device::Region> regions, const std::function<int(const mesh::Point&)>& region_of,
             bool three_d = false) {
    const auto x = uniform(0.0, 2 * um, 9);
    const auto y = uniform(0.0, 1 * um, 11);
    const auto z = three_d ? uniform(0.0, 0.6 * um, 4) : std::vector<double>{};
    mesh::Mesh m = three_d ? *mesh::make_tensor_grid(x, y, z) : *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        region[i] = static_cast<device::RegionId>(region_of(m.points()[i]));
        if (!device::is_insulator(regions[static_cast<std::size_t>(region[i])])) acceptors[i] = 1e16;
    }
    auto contact = m.find_boundary("y_min")->nodes;  // the bottom: semiconductor in every device here
    auto cells = three_d ? *mesh::TensorCells::create(m, x, y, z) : *mesh::TensorCells::create(m, x, y);
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = std::move(regions),
                                     .node_region = std::move(region),
                                     .donors = std::vector<double>(n, 0.0),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"left", device::ContactKind::ohmic, std::move(contact)}},
                                     .cells = std::move(cells)});
    if (!d) FAIL(d.error().message);
    return {std::move(*d)};
}

// Silicon left of x = 1 um, silicon with a 0.2 eV larger affinity from x = 1 um (a heterojunction).
Dev hetero(bool three_d = false) {
    return make_dev({{"left", physics::silicon()}, {"right", silicon_with(0.2)}},
                    [](const mesh::Point& p) { return p[0] > 0.99 * um ? 1 : 0; }, three_d);
}

// Oxide above y = 0.95 um on silicon.
Dev mos_like(bool three_d = false) {
    return make_dev({{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
                    [](const mesh::Point& p) { return p[1] > 0.95 * um ? 0 : 1; }, three_d);
}

Geometry geometry_of(const device::Device& d) { return *Geometry::tensor(d.mesh(), *d.cells()); }

}  // namespace

TEST_CASE("band surface: material keys follow the solver's parameter rule (Unit 26o, BSD2)") {
    const Dev a = mos_like();
    const auto ka = band_material_keys(a.device);
    REQUIRE(ka.size() == a.device.mesh().node_count());
    for (std::size_t i = 0; i < ka.size(); ++i) CHECK(ka[i] == (a.device.node_region()[i] == 0 ? -1 : 0));

    const Dev b = make_dev({{"si1", physics::silicon()}, {"si2", physics::silicon()}, {"oxide", physics::silicon_dioxide()}},
                           [](const mesh::Point& p) { return p[1] > 0.95 * um ? 2 : (p[0] > 0.99 * um ? 1 : 0); });
    const auto kb = band_material_keys(b.device);
    for (std::size_t i = 0; i < kb.size(); ++i) CHECK(kb[i] == (b.device.node_region()[i] == 2 ? -1 : 0));

    const Dev c = hetero();
    const auto kc = band_material_keys(c.device);
    for (std::size_t i = 0; i < kc.size(); ++i) CHECK(kc[i] == c.device.node_region()[i]);
}

TEST_CASE("band surface: one material reads as the field views; a step reads each side's own nodes (Unit 26o, BSD2)") {
    const Dev h = hetero();
    const Geometry g = geometry_of(h.device);
    const auto keys = band_material_keys(h.device);
    const auto pts = h.device.mesh().points();
    std::vector<double> level(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) level[i] = pts[i][0] / um + (keys[i] == 1 ? 10.0 : 0.0);
    std::size_t single = 0, mixed = 0;
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        const Polygon p = g.element_polygon(e);
        std::vector<std::int32_t> ks;
        for (const NodeId v : g.element(e)) ks.push_back(keys[static_cast<std::size_t>(v)]);
        const bool one = std::ranges::all_of(ks, [&](std::int32_t k) { return k == ks[0]; });
        for (const OwnedPolygon& piece : g.owner_pieces(p)) {
            mesh::Point c{};
            for (const mesh::Point& q : piece.points) for (int a = 0; a < 3; ++a) c[a] += q[a] / piece.points.size();
            const double v = band_value(g, e, c, piece.owner, keys, level);
            if (one) {
                ++single;
                CHECK(v == interpolate(g.interpolant(e, c), level));
            } else {
                ++mixed;
                if (keys[static_cast<std::size_t>(piece.owner)] == 1) {
                    CHECK(v >= 10.0);
                } else {
                    CHECK(v < 10.0);
                }
            }
        }
    }
    CHECK(single > 0);
    CHECK(mixed > 0);
}

TEST_CASE("band surface: the value is finite in every corner piece of a tensor cell (Unit 26o, Review Focus 1)") {
    const Dev m = mos_like();
    const Geometry g = geometry_of(m.device);
    const auto keys = band_material_keys(m.device);
    const auto pts = m.device.mesh().points();
    std::vector<double> level(pts.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (keys[i] >= 0) level[i] = 0.3 + pts[i][0] / um + 2.0 * pts[i][1] / um;
    }
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        for (const OwnedPolygon& piece : g.owner_pieces(g.element_polygon(e))) {
            if (keys[static_cast<std::size_t>(piece.owner)] < 0) continue;
            mesh::Point c{};
            for (const mesh::Point& q : piece.points) for (int a = 0; a < 3; ++a) c[a] += q[a] / piece.points.size();
            const double v = band_value(g, e, c, piece.owner, keys, level);
            CHECK(std::isfinite(v));
        }
    }
}

namespace {

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();

// Bands from a function of the node's position (um) and material key; NaN on insulator nodes.
results::BandDiagram bands_of(const device::Device& d,
                              const std::function<double(double, double, std::int32_t)>& ec) {
    const auto keys = band_material_keys(d);
    results::BandDiagram b;
    for (std::size_t i = 0; i < keys.size(); ++i) {
        const auto& p = d.mesh().points()[i];
        const double c = keys[i] < 0 ? not_a_number : ec(p[0] / um, p[1] / um, keys[i]);
        b.conduction_eV.push_back(c);
        b.valence_eV.push_back(c - 1.12);
        b.electron_fermi_eV.push_back(keys[i] < 0 ? not_a_number : -0.1 * p[0] / um);
        b.hole_fermi_eV.push_back(keys[i] < 0 ? not_a_number : 0.0);
    }
    return b;
}

double linear_ec(double x, double y, std::int32_t) { return 0.5 + 0.2 * x + 0.1 * y; }

const TriangleBatch* batch(const FieldScene& s, std::string_view name) {
    for (const TriangleBatch& b : s.triangles) {
        if (b.name == name) return &b;
    }
    return nullptr;
}

bool is_wall(const SceneVertex& v) { return std::abs(v.normal[1]) < 1e-6f; }

}  // namespace

TEST_CASE("band surface: 2D surfaces agree with the cut-line band diagram (Unit 26o, BSC2)") {
    const Dev m = mos_like();
    const Geometry g = geometry_of(m.device);
    const auto keys = band_material_keys(m.device);
    // A smooth, non-linear Ec: the double-precision value against the cut line.
    const auto wavy =
        bands_of(m.device, [](double x, double y, std::int32_t) { return 0.5 + 0.3 * std::sin(x) + 0.1 * y; });
    for (const double y : {0.2, 0.4, 0.6, 0.8}) {
        const auto cut = cut_line(g, {0.0, y * um, 0.0}, {2 * um, y * um, 0.0});
        REQUIRE(cut.has_value());
        const auto plot = cut_band_plot(m.device.mesh(), *cut, wavy);
        REQUIRE(plot.has_value());
        REQUIRE(plot->series.size() == 4);
        REQUIRE(plot->series[0].y.size() == cut->samples.size());
        for (std::size_t i = 0; i < cut->samples.size(); ++i) {
            const CutSample& s = cut->samples[i];
            mesh::Point x{};
            for (std::size_t k = 0; k < 4; ++k) {
                if (s.at.weights[k] == 0.0) continue;
                const auto& p = m.device.mesh().points()[static_cast<std::size_t>(s.at.nodes[k])];
                for (int a = 0; a < 3; ++a) x[a] += s.at.weights[k] * p[a];
            }
            bool one = true;
            for (const NodeId v : g.element(s.element)) one = one && keys[static_cast<std::size_t>(v)] == 1;
            if (!one) continue;
            for (std::size_t lv = 0; lv < 4; ++lv) {
                const auto e = band_energy_at(m.device, g, wavy, static_cast<BandLevel>(lv), s.element, x);
                REQUIRE(e.has_value());
                CHECK(std::abs(*e - plot->series[lv].y[i]) < 1e-9);
            }
        }
    }
    // A linear Ec: every drawn surface vertex sits at the field's own value (float rounding).
    const auto flat = bands_of(m.device, linear_ec);
    BandSurfaceOptions o;
    o.levels = {true, true, true, true};
    const auto s = band_surface_view(m.device, g, flat, o);
    REQUIRE(s.has_value());
    const double span = s->range_eV.second - s->range_eV.first;
    for (const auto& [name, shift] : std::vector<std::pair<std::string_view, double>>{{"Ec", 0.0}, {"Ev", -1.12}}) {
        const TriangleBatch* b = batch(s->scene, name);
        REQUIRE(b != nullptr);
        REQUIRE_FALSE(b->vertices.empty());
        for (const SceneVertex& v : b->vertices) {
            REQUIRE(std::isfinite(v.position[0]));
            REQUIRE(std::isfinite(v.position[1]));
            REQUIRE(std::isfinite(v.position[2]));
            // Single-material elements only: the silicon's top row is y = 0.9 um; above it the elements reach the
            // oxide and read the silicon nodes alone (BSD2).
            if (v.position[2] > 0.9f + 1e-4f) continue;
            const double e = v.position[1] / s->um_per_eV;
            CHECK(std::abs(e - (linear_ec(v.position[0], v.position[2], 0) + shift)) < 1e-6 * span + 1e-6);
        }
    }
    CHECK(batch(s->scene, "EFn") != nullptr);
    CHECK(batch(s->scene, "EFp") != nullptr);
    REQUIRE(s->scene.axes.has_value());
    CHECK((*s->scene.axes)[0].title == "x (\xC2\xB5m)");
    CHECK((*s->scene.axes)[1].title == "E (eV)");
    CHECK((*s->scene.axes)[1].per_unit == s->um_per_eV);
    CHECK((*s->scene.axes)[2].title == "y (\xC2\xB5m)");
    CHECK(s->scene.dimension == 3);
    CHECK_FALSE(s->scene.ground_shadow);
}

TEST_CASE("band surface: a heterojunction is a wall at the region boundary (Unit 26o, BSC3)") {
    const Dev h = hetero();
    const Geometry g = geometry_of(h.device);
    // Ec linear in x on the left, a 0.2 eV step to the right.
    const auto b = bands_of(h.device, [](double x, double, std::int32_t k) { return 0.4 + 0.05 * x + (k == 1 ? 0.2 : 0.0); });
    BandSurfaceOptions o;
    o.levels = {true, false, false, false};
    const auto s = band_surface_view(h.device, g, b, o);
    REQUIRE(s.has_value());
    const TriangleBatch* ec = batch(s->scene, "Ec");
    REQUIRE(ec != nullptr);
    // The boundary is midway between the left nodes at x = 0.75 um and the right ones at 1 um; each side reads only
    // its own nodes there (renormalised), so the left at x = 0.75 um and the right at x = 1 um.
    const double bx = 0.875, low = 0.4 + 0.05 * 0.75, high = 0.4 + 0.05 * 1.0 + 0.2;
    std::size_t walls = 0;
    bool saw_low = false, saw_high = false;
    const auto& t = ec->vertices;
    for (std::size_t i = 0; i + 2 < t.size(); i += 3) {
        if (is_wall(t[i])) {
            ++walls;
            for (std::size_t k = i; k < i + 3; ++k) {
                CHECK(std::abs(t[k].position[0] - bx) < 1e-4);
                const double e = t[k].position[1] / s->um_per_eV;
                const bool at_low = std::abs(e - low) < 1e-5, at_high = std::abs(e - high) < 1e-5;
                CHECK((at_low || at_high));
                saw_low = saw_low || at_low;
                saw_high = saw_high || at_high;
            }
        } else {
            // A surface triangle lies on one side of the boundary.
            bool left = true, right = true;
            for (std::size_t k = i; k < i + 3; ++k) {
                left = left && t[k].position[0] <= bx + 1e-4;
                right = right && t[k].position[0] >= bx - 1e-4;
            }
            CHECK((left || right));
        }
    }
    CHECK(walls > 0);
    CHECK(saw_low);
    CHECK(saw_high);
}

TEST_CASE("band surface: the oxide is a gap (Unit 26o, BSC3, BSQ3)") {
    const Dev m = mos_like();
    const auto s = band_surface_view(m.device, geometry_of(m.device), bands_of(m.device, linear_ec));
    REQUIRE(s.has_value());
    for (const TriangleBatch& b : s->scene.triangles) {
        for (const SceneVertex& v : b.vertices) {
            CHECK(v.position[2] <= 0.95f + 1e-4f);
            CHECK_FALSE(is_wall(v));
        }
    }
}

TEST_CASE("band surface: regions of one material make no wall (Unit 26o, Review Focus 2)") {
    const Dev d = make_dev({{"si1", physics::silicon()}, {"si2", physics::silicon()}},
                           [](const mesh::Point& p) { return p[0] > 0.99 * um ? 1 : 0; });
    const auto s = band_surface_view(d.device, geometry_of(d.device), bands_of(d.device, linear_ec));
    REQUIRE(s.has_value());
    for (const TriangleBatch& b : s->scene.triangles) {
        for (const SceneVertex& v : b.vertices) CHECK_FALSE(is_wall(v));
    }
}

TEST_CASE("band surface: refusals and dropped pieces (Unit 26o, BSC6)") {
    const Dev m = mos_like();
    const Geometry g = geometry_of(m.device);
    const auto good = bands_of(m.device, linear_ec);
    const auto refused = [&](const results::BandDiagram& b, const BandSurfaceOptions& o) {
        const auto s = band_surface_view(m.device, g, b, o);
        return !s.has_value() && s.error().code == base::ErrorCode::invalid_input;
    };
    const BandSurfaceOptions o;
    BandSurfaceOptions plane = o;
    plane.plane = Plane{{0.0, 0.0, 1.0}, 0.0};
    CHECK(refused(good, plane));  // a plane on a 2D device
    auto short_bands = good;
    short_bands.valence_eV.pop_back();
    CHECK(refused(short_bands, o));
    BandSurfaceOptions none = o;
    none.levels = {false, false, false, false};
    CHECK(refused(good, none));
    BandSurfaceOptions ex = o;
    ex.exaggeration = 0.05;
    CHECK(refused(good, ex));
    ex.exaggeration = 11.0;
    CHECK(refused(good, ex));
    BandSurfaceOptions op = o;
    op.opacity = 0.1f;
    CHECK(refused(good, op));
    BandSurfaceOptions range = o;
    range.range_eV = std::pair{1.0, 1.0};
    CHECK(refused(good, range));
    auto all_nan = good;
    for (auto* v : {&all_nan.conduction_eV, &all_nan.valence_eV, &all_nan.electron_fermi_eV, &all_nan.hole_fermi_eV}) {
        std::ranges::fill(*v, not_a_number);
    }
    CHECK(refused(all_nan, o));
    // A NaN and an infinity at two semiconductor nodes: their pieces dropped and counted, nothing non-finite drawn.
    auto bad = good;
    const auto keys = band_material_keys(m.device);
    std::size_t marked = 0;
    for (std::size_t i = 0; i < keys.size() && marked < 2; ++i) {
        if (keys[i] < 0) continue;
        bad.conduction_eV[i] = marked == 0 ? not_a_number : std::numeric_limits<double>::infinity();
        ++marked;
    }
    const auto s = band_surface_view(m.device, g, bad, o);
    REQUIRE(s.has_value());
    CHECK(s->dropped > 0);
    for (const TriangleBatch& b : s->scene.triangles) {
        for (const SceneVertex& v : b.vertices) {
            CHECK(std::isfinite(v.position[0]));
            CHECK(std::isfinite(v.position[1]));
            CHECK(std::isfinite(v.position[2]));
            if (b.shading == Shading::field) CHECK(std::isfinite(v.value));
        }
    }
    // Flat bands (Review Focus 5): a finite height scale over a widened range.
    const auto level =
        band_surface_view(m.device, g, bands_of(m.device, [](double, double, std::int32_t) { return 0.3; }), o);
    REQUIRE(level.has_value());
    CHECK(std::isfinite(level->um_per_eV));
    CHECK(level->um_per_eV > 0.0);
    CHECK(level->range_eV.first < level->range_eV.second);
}

TEST_CASE("band surface: colour by energy or by level (Unit 26o, BSD4)") {
    const Dev m = mos_like();
    const Geometry g = geometry_of(m.device);
    const auto b = bands_of(m.device, linear_ec);
    BandSurfaceOptions o;
    o.opacity = 0.5f;
    const auto by_energy = band_surface_view(m.device, g, b, o);
    REQUIRE(by_energy.has_value());
    REQUIRE(by_energy->scene.colorbar.has_value());
    CHECK(by_energy->scene.colorbar->title == "E (eV)");
    CHECK(by_energy->scene.colorbar->map == Colormap::viridis);
    for (const TriangleBatch& t : by_energy->scene.triangles) {
        CHECK(t.shading == Shading::field);
        CHECK(t.transparent);
        for (const SceneVertex& v : t.vertices) CHECK(v.color.a == 0.5f);
    }
    o.by_level = true;
    o.opacity = 1.0f;
    const auto by_level = band_surface_view(m.device, g, b, o);
    REQUIRE(by_level.has_value());
    CHECK_FALSE(by_level->scene.colorbar.has_value());
    const TriangleBatch* ec = batch(by_level->scene, "Ec");
    const TriangleBatch* ev = batch(by_level->scene, "Ev");
    REQUIRE(ec != nullptr);
    REQUIRE(ev != nullptr);
    CHECK(ec->shading == Shading::flat);
    CHECK_FALSE(ec->transparent);
    CHECK(ec->vertices.front().color.r == 0.12f);
    CHECK(ec->vertices.front().color.b == 0.71f);
    CHECK(ev->vertices.front().color.r == 0.84f);
}

TEST_CASE("band surface: exaggeration scales heights only (Unit 26o, BSD3)") {
    const Dev m = mos_like();
    const Geometry g = geometry_of(m.device);
    const auto b = bands_of(m.device, linear_ec);
    BandSurfaceOptions o;
    const auto one = band_surface_view(m.device, g, b, o);
    o.exaggeration = 2.0;
    const auto two = band_surface_view(m.device, g, b, o);
    REQUIRE(one.has_value());
    REQUIRE(two.has_value());
    CHECK(two->um_per_eV == 2.0 * one->um_per_eV);
    REQUIRE(one->scene.triangles.size() == two->scene.triangles.size());
    for (std::size_t t = 0; t < one->scene.triangles.size(); ++t) {
        const auto& a = one->scene.triangles[t].vertices;
        const auto& c = two->scene.triangles[t].vertices;
        REQUIRE(a.size() == c.size());
        for (std::size_t i = 0; i < a.size(); ++i) {
            CHECK(a[i].position[0] == c[i].position[0]);
            CHECK(a[i].position[2] == c[i].position[2]);
            CHECK(std::abs(c[i].position[1] - 2.0f * a[i].position[1]) <= 1e-5f * std::abs(c[i].position[1]) + 1e-6f);
        }
    }
}
