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
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/band_surface.hpp"
#include "NiTCAD/render/geometry.hpp"

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
