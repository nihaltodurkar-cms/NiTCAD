// The band surface (Unit 26o; DECISIONS.md BSD1-BSD4, BSC2, BSC3, BSC6, BSC8): Ec, Ev, EFn and EFp lifted as
// energy-height surfaces over a 2D device or a plane of a 3D one, steps kept at material boundaries.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
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
#include "NiTCAD/render/axes3d.hpp"
#include "NiTCAD/render/band_surface.hpp"
#include "NiTCAD/render/cut_line.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/render/edges.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/results/solution.hpp"
#include "image_check.hpp"

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
    // Flat bands (Review Focus 5): Ec alone, one value everywhere; a finite height scale over a widened range.
    BandSurfaceOptions ec_only = o;
    ec_only.levels = {true, false, false, false};
    const auto level =
        band_surface_view(m.device, g, bands_of(m.device, [](double, double, std::int32_t) { return 0.3; }), ec_only);
    REQUIRE(level.has_value());
    CHECK(std::isfinite(level->um_per_eV));
    CHECK(level->um_per_eV > 0.0);
    CHECK(level->range_eV.second - level->range_eV.first >= 1.0 - 1e-9);  // widened by 0.5 eV each side
    CHECK(level->um_per_eV < 10.0);
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

namespace {

// For each cut line (y in um, along x at z in um), each sample in a silicon-only element: band_energy_at equals the
// band diagram's value there to 1e-9 eV, for all four levels.
void check_against_cuts(const device::Device& d, const Geometry& g, const results::BandDiagram& b, double z) {
    const auto keys = band_material_keys(d);
    std::size_t checked = 0;
    for (const double y : {0.2, 0.4, 0.6, 0.8}) {
        const auto cut = cut_line(g, {0.0, y * um, z * um}, {2 * um, y * um, z * um});
        REQUIRE(cut.has_value());
        const auto plot = cut_band_plot(d.mesh(), *cut, b);
        REQUIRE(plot.has_value());
        REQUIRE(plot->series[0].y.size() == cut->samples.size());
        for (std::size_t i = 0; i < cut->samples.size(); ++i) {
            const CutSample& s = cut->samples[i];
            mesh::Point x{};
            for (std::size_t k = 0; k < 4; ++k) {
                if (s.at.weights[k] == 0.0) continue;
                const auto& p = d.mesh().points()[static_cast<std::size_t>(s.at.nodes[k])];
                for (int a = 0; a < 3; ++a) x[a] += s.at.weights[k] * p[a];
            }
            bool one = true;
            for (const NodeId v : g.element(s.element)) one = one && keys[static_cast<std::size_t>(v)] == keys[static_cast<std::size_t>(g.element(s.element)[0])] && keys[static_cast<std::size_t>(v)] >= 0;
            if (!one) continue;
            for (std::size_t lv = 0; lv < 4; ++lv) {
                const auto e = band_energy_at(d, g, b, static_cast<BandLevel>(lv), s.element, x);
                REQUIRE(e.has_value());
                CHECK(std::abs(*e - plot->series[lv].y[i]) < 1e-9);
                ++checked;
            }
        }
    }
    CHECK(checked > 0);
}

}  // namespace

TEST_CASE("band surface: a 3D plane agrees with the cut line (Unit 26o, BSC2)") {
    const Dev m = mos_like(true);
    const Geometry g = geometry_of(m.device);
    const auto wavy =
        bands_of(m.device, [](double x, double y, std::int32_t) { return 0.5 + 0.3 * std::sin(x) + 0.1 * y; });
    check_against_cuts(m.device, g, wavy, 0.3);
    // The plane z = 0.3 um: base axes x and y; a linear Ec is drawn at its own value.
    BandSurfaceOptions o;
    o.plane = Plane{{0.0, 0.0, 1.0}, 0.3 * um};
    const auto s = band_surface_view(m.device, g, bands_of(m.device, linear_ec), o);
    REQUIRE(s.has_value());
    CHECK((*s->scene.axes)[0].title == "x (\xC2\xB5m)");
    CHECK((*s->scene.axes)[2].title == "y (\xC2\xB5m)");
    const double span = s->range_eV.second - s->range_eV.first;
    const TriangleBatch* ec = batch(s->scene, "Ec");
    REQUIRE(ec != nullptr);
    REQUIRE_FALSE(ec->vertices.empty());
    for (const SceneVertex& v : ec->vertices) {
        if (v.position[2] > 0.9f + 1e-4f) continue;
        CHECK(std::abs(v.position[1] / s->um_per_eV - linear_ec(v.position[0], v.position[2], 0)) < 1e-6 * span + 1e-6);
    }
    // A plane normal to x: base axes y and z.
    o.plane = Plane{{1.0, 0.0, 0.0}, 1.1 * um};
    const auto across = band_surface_view(m.device, g, bands_of(m.device, linear_ec), o);
    REQUIRE(across.has_value());
    CHECK((*across->scene.axes)[0].title == "y (\xC2\xB5m)");
    CHECK((*across->scene.axes)[2].title == "z (\xC2\xB5m)");
    // Off the device, and only in the oxide.
    o.plane = Plane{{0.0, 0.0, 1.0}, 5.0 * um};
    const auto off = band_surface_view(m.device, g, bands_of(m.device, linear_ec), o);
    REQUIRE_FALSE(off.has_value());
    CHECK(off.error().message.find("misses the device") != std::string::npos);
    o.plane = Plane{{0.0, 1.0, 0.0}, 0.98 * um};
    const auto oxide = band_surface_view(m.device, g, bands_of(m.device, linear_ec), o);
    REQUIRE_FALSE(oxide.has_value());
    CHECK(oxide.error().message.find("no band values on the plane") != std::string::npos);
}

TEST_CASE("band surface: a plane on a mesh plane between materials is nudged (Unit 26o, Review Focus 4)") {
    const Dev h = hetero(true);
    const Geometry g = geometry_of(h.device);
    const auto b = bands_of(h.device, [](double x, double, std::int32_t k) { return 0.4 + 0.05 * x + (k == 1 ? 0.2 : 0.0); });
    BandSurfaceOptions o;
    o.levels = {true, false, false, false};
    o.plane = Plane{{1.0, 0.0, 0.0}, 1.0 * um};  // the right material's first node plane
    const auto s = band_surface_view(h.device, g, b, o);
    REQUIRE(s.has_value());
    const TriangleBatch* ec = batch(s->scene, "Ec");
    REQUIRE(ec != nullptr);
    REQUIRE_FALSE(ec->vertices.empty());
    // Nudged a hair toward lower x: still the right nodes' control volumes, which read x = 1 um's 0.65 eV.
    for (const SceneVertex& v : ec->vertices) CHECK(std::abs(v.position[1] / s->um_per_eV - 0.65) < 1e-5);
}

TEST_CASE("band surface: the energy range clips pieces (Unit 26o, BSD3)") {
    const Dev m = mos_like();
    const Geometry g = geometry_of(m.device);
    const auto b = bands_of(m.device, linear_ec);
    BandSurfaceOptions o;
    o.levels = {true, false, false, false};
    const auto whole = band_surface_view(m.device, g, b, o);
    o.range_eV = std::pair{0.6, 0.8};
    const auto clipped = band_surface_view(m.device, g, b, o);
    REQUIRE(whole.has_value());
    REQUIRE(clipped.has_value());
    CHECK(clipped->range_eV == std::pair{0.6, 0.8});
    CHECK(clipped->dropped == whole->dropped);
    const TriangleBatch* ec = batch(clipped->scene, "Ec");
    REQUIRE(ec != nullptr);
    REQUIRE_FALSE(ec->vertices.empty());
    for (const SceneVertex& v : ec->vertices) {
        const double e = v.position[1] / clipped->um_per_eV;
        CHECK(e >= 0.6 - 1e-5);
        CHECK(e <= 0.8 + 1e-5);
    }
    // Less of the base is covered than unclipped (projected on the base plane).
    const auto base_area = [](const TriangleBatch& t) {
        double a = 0.0;
        for (std::size_t i = 0; i + 2 < t.vertices.size(); i += 3) {
            const auto& p = t.vertices;
            const double ux = p[i + 1].position[0] - p[i].position[0], uz = p[i + 1].position[2] - p[i].position[2];
            const double vx = p[i + 2].position[0] - p[i].position[0], vz = p[i + 2].position[2] - p[i].position[2];
            a += 0.5 * std::abs(ux * vz - uz * vx);
        }
        return a;
    };
    CHECK(base_area(*ec) < base_area(*batch(whole->scene, "Ec")));
}

TEST_CASE("band surface: reference images (Unit 26o, BSC8)") {
    constexpr int W = 960, H = 720;
    Renderer& r = render_test::warp();
    const render_test::Mosfet& m2 = render_test::mosfet(false);
    const render_test::Mosfet& m3 = render_test::mosfet(true);
    {
        const Geometry g = geometry_of(m2.device);
        auto s = band_surface_view(m2.device, g, m2.point.bands);
        REQUIRE(s.has_value());
        set_look(s->scene, Look::enhanced);
        render_test::check_image(*r.draw(s->scene, W, H), "bandsurface_mosfet2d");
    }
    {
        const Geometry g = geometry_of(m3.device);
        BandSurfaceOptions o;
        o.plane = Plane{{0.0, 0.0, 1.0}, 0.4 * um};  // the gate's centre, source-drain by depth
        auto s = band_surface_view(m3.device, g, m3.point.bands, o);
        REQUIRE(s.has_value());
        set_look(s->scene, Look::enhanced);
        render_test::check_image(*r.draw(s->scene, W, H), "bandsurface_mosfet3d");
    }
    {
        const Dev h = hetero();
        BandSurfaceOptions o;
        o.by_level = true;
        auto s = band_surface_view(h.device, geometry_of(h.device),
                                   bands_of(h.device, [](double x, double, std::int32_t k) {
                                       return 0.4 + 0.05 * x + (k == 1 ? 0.2 : 0.0);
                                   }),
                                   o);
        REQUIRE(s.has_value());
        set_look(s->scene, Look::enhanced);
        render_test::check_image(*r.draw(s->scene, W, H), "bandsurface_hetero");
    }
}

#ifdef NDEBUG
TEST_CASE("band surface: building the scene is quick (Unit 26o, BSC9; Release only)") {
    const auto median_ms = [](const std::function<void()>& build) {
        std::vector<double> t;
        for (int k = 0; k < 20; ++k) {
            const auto a = std::chrono::steady_clock::now();
            build();
            t.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count());
        }
        std::ranges::sort(t);
        return 0.5 * (t[9] + t[10]);
    };
    const render_test::Mosfet& m2 = render_test::mosfet(false);
    const render_test::Mosfet& m3 = render_test::mosfet(true);
    const Geometry g2 = geometry_of(m2.device);
    const Geometry g3 = geometry_of(m3.device);
    BandSurfaceOptions o3;
    o3.plane = Plane{{0.0, 0.0, 1.0}, 0.4 * um};
    const double two = median_ms([&] { REQUIRE(band_surface_view(m2.device, g2, m2.point.bands).has_value()); });
    const double three = median_ms([&] { REQUIRE(band_surface_view(m3.device, g3, m3.point.bands, o3).has_value()); });
    WARN("band surface build median: 2D " << two << " ms, 3D plane " << three << " ms");
    CHECK(two <= 150.0);
    CHECK(three <= 400.0);
}
#endif

TEST_CASE("band surface: the plane's outline on a device view (Unit 26o, BSD6)") {
    const Dev m = mos_like(true);
    const Geometry g = geometry_of(m.device);
    FieldScene s;
    add_plane_outline(s, g, Plane{{0.0, 0.0, 1.0}, 0.3 * um});
    REQUIRE(s.lines.size() == 1);
    const LineBatch& l = s.lines[0];
    CHECK(l.name == "band plane");
    CHECK(l.on_top);
    REQUIRE(l.points.size() % 2 == 0);
    REQUIRE_FALSE(l.points.empty());
    double length = 0.0;
    for (std::size_t i = 0; i < l.points.size(); i += 2) {
        CHECK(std::abs(l.points[i][2] - 0.3f) < 1e-5f);
        CHECK(std::abs(l.points[i + 1][2] - 0.3f) < 1e-5f);
        length += std::hypot(l.points[i + 1][0] - l.points[i][0], l.points[i + 1][1] - l.points[i][1]);
    }
    CHECK(std::abs(length - 6.0) < 1e-4);  // the section's perimeter, 2 um x 1 um, the oxide included
    // A plane off the device draws nothing.
    FieldScene none;
    add_plane_outline(none, g, Plane{{0.0, 0.0, 1.0}, 5.0 * um});
    CHECK(none.lines.empty());
}

TEST_CASE("band surface: the plane inset (Unit 26o, BSD6)") {
    const Dev m = mos_like(true);
    const Geometry g = geometry_of(m.device);
    BandSurfaceOptions o;
    o.plane = Plane{{0.0, 0.0, 1.0}, 0.3 * um};
    auto s = band_surface_view(m.device, g, bands_of(m.device, linear_ec), o);
    REQUIRE(s.has_value());
    REQUIRE(s->scene.plane_inset.has_value());
    const PlaneInset& p = *s->scene.plane_inset;
    CHECK(p.axis == 2);
    CHECK(std::abs(p.at_um - 0.3) < 1e-12);
    CHECK(std::abs(p.hi_um[0] - 2.0) < 1e-9);
    CHECK(std::abs(p.hi_um[1] - 1.0) < 1e-9);
    CHECK(std::abs(p.hi_um[2] - 0.6) < 1e-9);
    CHECK(p.lo_um == std::array<double, 3>{0.0, 0.0, 0.0});
    set_look(s->scene, Look::enhanced);
    constexpr int W = 960, H = 720;
    const InsetFigure f = plane_inset_figure(s->scene, W, H);
    REQUIRE(f.edges.size() == 12);
    // Top right of the view, left of the colour bar, clear of the view cube (bottom left).
    CHECK(f.frame[2] <= W - 130.0f);
    CHECK(f.frame[0] > W / 2.0f);
    CHECK(f.frame[3] < H / 2.0f);
    const auto inside = [&](const std::array<float, 2>& q) {
        return q[0] >= f.frame[0] - 1e-3f && q[0] <= f.frame[2] + 1e-3f && q[1] >= f.frame[1] - 1e-3f && q[1] <= f.frame[3] + 1e-3f;
    };
    for (const auto& e : f.edges) {
        CHECK(inside(e[0]));
        CHECK(inside(e[1]));
    }
    double area = 0.0;
    for (std::size_t k = 0; k < 4; ++k) {
        CHECK(inside(f.plane[k]));
        const auto& a = f.plane[k];
        const auto& b = f.plane[(k + 1) % 4];
        area += 0.5 * (double(a[0]) * b[1] - double(b[0]) * a[1]);
    }
    CHECK(std::abs(area) > 100.0);
    // It turns with the camera.
    FieldScene turned = s->scene;
    turned.camera.azimuth_deg += 40.0f;
    CHECK(plane_inset_figure(turned, W, H).plane != f.plane);
    // 2D: none.
    auto flat = band_surface_view(mos_like().device, geometry_of(mos_like().device), bands_of(mos_like().device, linear_ec));
    REQUIRE(flat.has_value());
    CHECK_FALSE(flat->scene.plane_inset.has_value());
    CHECK(plane_inset_figure(flat->scene, W, H).edges.empty());
}

TEST_CASE("band surface: the inset's reference image (Unit 26o, BSC8)") {
    Renderer& r = render_test::warp();
    const render_test::Mosfet& m3 = render_test::mosfet(true);
    const Geometry g = geometry_of(m3.device);
    BandSurfaceOptions o;
    o.plane = Plane{{1.0, 0.0, 0.0}, 1.0 * um};  // across the device under the gate
    o.view.azimuth_deg = 25.0f;
    auto s = band_surface_view(m3.device, g, m3.point.bands, o);
    REQUIRE(s.has_value());
    set_look(s->scene, Look::enhanced);
    render_test::check_image(*r.draw(s->scene, 960, 720), "bandsurface_inset");
}

TEST_CASE("band surface: one energy range over a run's points (Unit 26o, BSD8)") {
    const Dev m = mos_like();
    const auto make = [&](double shift) {
        results::BiasPoint p;
        p.bands = bands_of(m.device, [shift](double x, double, std::int32_t) { return 0.2 + 0.1 * x + shift; });
        return p;
    };
    const std::vector<results::BiasPoint> run{make(0.0), make(0.3), make(-0.2)};
    const auto r = band_energy_range(m.device, run, {true, true, false, false});
    REQUIRE(r);
    // Ec spans 0.2 + [0, 0.2] + shifts; Ev is Ec - 1.12.
    CHECK(std::abs(r->first - (0.2 - 0.2 - 1.12)) < 1e-12);
    CHECK(std::abs(r->second - (0.2 + 0.2 + 0.3)) < 1e-12);
    // EFn alone (-0.1 x, the same at every point): its own span.
    const auto efn = band_energy_range(m.device, run, {false, false, true, false});
    REQUIRE(efn);
    CHECK(std::abs(efn->first + 0.2) < 1e-12);
    CHECK(std::abs(efn->second - 0.0) < 1e-12);
    // EFp (0 everywhere): widened as a flat surface is.
    const auto efp = band_energy_range(m.device, run, {false, false, false, true});
    REQUIRE(efp);
    CHECK(*efp == std::pair{-0.5, 0.5});
    // Nothing finite: none.
    results::BiasPoint empty;
    empty.bands = bands_of(m.device, [](double, double, std::int32_t) { return not_a_number; });
    CHECK_FALSE(band_energy_range(m.device, std::vector<results::BiasPoint>{empty}, {true, false, false, false}));
}

TEST_CASE("band surface: a cut line in the plane drawn at the band heights (Unit 26o, BSD9)") {
    const auto check_cut = [](const Dev& m, const std::optional<Plane>& plane, double z_um) {
        const Geometry g = geometry_of(m.device);
        const auto b = bands_of(m.device, linear_ec);
        auto cut = cut_line(g, {0.1 * um, 0.4 * um, z_um * um}, {1.9 * um, 0.4 * um, z_um * um});
        REQUIRE(cut.has_value());
        BandSurfaceOptions o;
        o.plane = plane;
        o.cut = *cut;
        const auto s = band_surface_view(m.device, g, b, o);
        REQUIRE(s.has_value());
        const LineBatch* line = nullptr;
        for (const LineBatch& l : s->scene.lines) {
            if (l.name == "band cut") line = &l;
        }
        REQUIRE(line != nullptr);
        CHECK(line->on_top);
        REQUIRE(line->points.size() % 2 == 0);
        REQUIRE(line->points.size() >= 4);
        // Ec and Ev both drawn, each point at a level's height for its (x, y).
        std::size_t on_ec = 0, on_ev = 0;
        for (const Float3& p : line->points) {
            CHECK(std::abs(p[2] - 0.4f) < 1e-5f);  // y = 0.4 um is scene z
            const double e = p[1] / s->um_per_eV;
            const double ec = linear_ec(p[0], 0.4, 0);
            if (std::abs(e - ec) < 1e-5) ++on_ec;
            else if (std::abs(e - (ec - 1.12)) < 1e-5) ++on_ev;
            else FAIL_CHECK("a point off both levels: E = " << e);
        }
        CHECK(on_ec > 0);
        CHECK(on_ev > 0);
        // Without the cut, no such batch.
        o.cut.reset();
        const auto none = band_surface_view(m.device, g, b, o);
        REQUIRE(none.has_value());
        CHECK(std::ranges::none_of(none->scene.lines, [](const LineBatch& l) { return l.name == "band cut"; }));
    };
    check_cut(mos_like(), std::nullopt, 0.0);
    check_cut(mos_like(true), Plane{{0.0, 0.0, 1.0}, 0.3 * um}, 0.3);
    // A 3D cut line off the plane is not drawn.
    const Dev m = mos_like(true);
    const Geometry g = geometry_of(m.device);
    auto off = cut_line(g, {0.1 * um, 0.4 * um, 0.1 * um}, {1.9 * um, 0.4 * um, 0.1 * um});
    REQUIRE(off.has_value());
    BandSurfaceOptions o;
    o.plane = Plane{{0.0, 0.0, 1.0}, 0.3 * um};
    o.cut = *off;
    const auto s = band_surface_view(m.device, g, bands_of(m.device, linear_ec), o);
    REQUIRE(s.has_value());
    CHECK(std::ranges::none_of(s->scene.lines, [](const LineBatch& l) { return l.name == "band cut"; }));
}

TEST_CASE("band surface: the cut line steps with the bands at a heterojunction (Unit 26o, BSD9)") {
    const Dev h = hetero();
    const Geometry g = geometry_of(h.device);
    const auto b = bands_of(h.device, [](double x, double, std::int32_t k) { return 0.4 + 0.05 * x + (k == 1 ? 0.2 : 0.0); });
    auto cut = cut_line(g, {0.1 * um, 0.5 * um, 0.0}, {1.9 * um, 0.5 * um, 0.0});
    REQUIRE(cut.has_value());
    BandSurfaceOptions o;
    o.levels = {true, false, false, false};
    o.cut = *cut;
    const auto s = band_surface_view(h.device, g, b, o);
    REQUIRE(s.has_value());
    const LineBatch* line = nullptr;
    for (const LineBatch& l : s->scene.lines) {
        if (l.name == "band cut") line = &l;
    }
    REQUIRE(line != nullptr);
    // A vertical segment of 0.2 eV near the wall (x = 0.875 um, within 1/16 of an element).
    bool step = false;
    for (std::size_t i = 0; i + 1 < line->points.size(); i += 2) {
        const Float3& a = line->points[i];
        const Float3& c = line->points[i + 1];
        if (a[0] == c[0] && std::abs(std::abs(c[1] - a[1]) / s->um_per_eV - 0.2) < 0.02 && std::abs(a[0] - 0.875f) < 0.02f) step = true;
    }
    CHECK(step);
}

namespace {

const LineBatch* line_batch(const FieldScene& s, std::string_view name) {
    for (const LineBatch& l : s.lines) {
        if (l.name == name) return &l;
    }
    return nullptr;
}

// The outline's length in the base plane (scene x and z, um), and whether every point lies on one of the given
// lines x = const or z = const (um).
struct Outline {
    double base_length = 0.0;
    bool on_lines = true;
};

Outline measure(const LineBatch& l, const std::vector<double>& xs, const std::vector<double>& zs) {
    Outline o;
    for (std::size_t i = 0; i + 1 < l.points.size(); i += 2) {
        const Float3& a = l.points[i];
        const Float3& b = l.points[i + 1];
        o.base_length += std::hypot(double(b[0]) - a[0], double(b[2]) - a[2]);
        for (const Float3& p : {a, b}) {
            const bool on = std::ranges::any_of(xs, [&](double x) { return std::abs(p[0] - x) < 1e-4; }) ||
                            std::ranges::any_of(zs, [&](double z) { return std::abs(p[2] - z) < 1e-4; });
            o.on_lines = o.on_lines && on;
        }
    }
    return o;
}

}  // namespace

TEST_CASE("band surface: its own outline, not the look's creases (Unit 26o, the stray-lines fix)") {
    const Dev h = hetero();
    const Geometry g = geometry_of(h.device);
    // A curved Ec, so the enhanced look's 30-degree creases would mark it.
    const auto b = bands_of(h.device, [](double x, double y, std::int32_t k) {
        return 0.4 + 0.8 * std::sin(3.0 * x) * y + (k == 1 ? 0.2 : 0.0);
    });
    auto s = band_surface_view(h.device, g, b);
    REQUIRE(s.has_value());
    for (const TriangleBatch& t : s->scene.triangles) CHECK_FALSE(t.feature_edges);
    set_look(s->scene, Look::enhanced);
    CHECK(line_batch(s->scene, "feature edges") == nullptr);
    // The outline: each material's rim on Ec and Ev, the step at the interface (x = 0.875 um) drawn on both sides.
    // Base length per level: left 2 x (0.875 + 1) plus right 2 x (1.125 + 1) = 8 um; two levels.
    const LineBatch* outline = line_batch(s->scene, "band outline");
    REQUIRE(outline != nullptr);
    const Outline o = measure(*outline, {0.0, 0.875, 2.0}, {0.0, 1.0});
    CHECK(o.on_lines);
    CHECK(std::abs(o.base_length - 16.0) < 1e-3);
    // The classic look keeps the outline (it is the view's, not the look's).
    set_look(s->scene, Look::classic);
    CHECK(line_batch(s->scene, "band outline") != nullptr);
    // Other views still take feature edges (BSC1).
    TriangleBatch other{"field", Shading::field, true, false, {}};
    CHECK(takes_feature_edges(other));
}

TEST_CASE("band surface: hidden regions are left out, a gap, and out of the run's range (Unit 26o, YD4)") {
    const Dev h = hetero();
    const Geometry g = geometry_of(h.device);
    const auto b = bands_of(h.device, [](double x, double, std::int32_t k) { return 0.4 + 0.05 * x + (k == 1 ? 0.5 : 0.0); });
    BandSurfaceOptions o;
    o.levels = {true, false, false, false};
    o.view.hidden_regions = {1};  // "right"
    const auto s = band_surface_view(h.device, g, b, o);
    REQUIRE(s.has_value());
    const TriangleBatch* ec = batch(s->scene, "Ec");
    REQUIRE(ec != nullptr);
    for (const SceneVertex& v : ec->vertices) CHECK(v.position[0] <= 0.875f + 1e-4f);
    // No wall: every vertex sits at its own (left) value, no step up to the right's. In the element the interface
    // cuts, the left piece reads the left nodes only (BSD2): flat at the x = 0.75 um node's value.
    for (const SceneVertex& v : ec->vertices) {
        CHECK(std::abs(v.position[1] / s->um_per_eV - (0.4 + 0.05 * std::min(double(v.position[0]), 0.75))) < 1e-5);
    }
    const LineBatch* outline = line_batch(s->scene, "band outline");
    REQUIRE(outline != nullptr);
    CHECK(std::abs(measure(*outline, {0.0, 0.875}, {0.0, 1.0}).base_length - 3.75) < 1e-3);
    // Hiding the left too leaves nothing to draw.
    o.view.hidden_regions = {0, 1};
    const auto none = band_surface_view(h.device, g, b, o);
    REQUIRE_FALSE(none.has_value());
    CHECK(none.error().message.find("no band values") != std::string::npos);
    // The run's range over the shown region only.
    results::BiasPoint p;
    p.bands = b;
    const std::vector<results::BiasPoint> run{p};
    const auto all = band_energy_range(h.device, run, {true, false, false, false});
    const std::vector<std::int32_t> right{1};
    const auto left = band_energy_range(h.device, run, {true, false, false, false}, right);
    REQUIRE(all);
    REQUIRE(left);
    CHECK(all->second > 0.9);
    CHECK(std::abs(left->first - 0.4) < 1e-12);
    CHECK(left->second < 0.45);
}
