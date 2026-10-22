// Run comparison (Unit 26b, XD6, XG5): B - A of a run with itself is zero, of two runs the nodal
// difference in the field's units (mV, meV, decades); the matching bias point; the symmetric scale.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/compare.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

// A silicon bar on a tensor mesh with a donor profile, at a fixed address with its geometry.
struct Bar {
    std::unique_ptr<device::Device> device;
    std::optional<Geometry> geometry;
};

Bar bar(double donors) {
    const auto x = uniform(0.0, 2e-4, 5), y = uniform(0.0, 1e-4, 3);
    auto m = mesh::make_tensor_grid(x, y);
    REQUIRE(m.has_value());
    auto c = mesh::TensorCells::create(*m, x, y);
    REQUIRE(c.has_value());
    const std::size_t n = m->node_count();
    std::vector<double> nd(n);
    for (std::size_t i = 0; i < n; ++i) nd[i] = donors * (1.0 + m->points()[i][0] * 1e4);
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(nd),
                                     .acceptors = std::vector<double>(n, 1e15),
                                     .contacts = {{"left", device::ContactKind::ohmic, {mesh::NodeId{0}}}},
                                     .cells = std::move(*c)});
    if (!d) FAIL(d.error().message);
    Bar out{std::make_unique<device::Device>(std::move(*d)), std::nullopt};
    auto g = Geometry::tensor(out.device->mesh(), *out.device->cells());
    REQUIRE(g.has_value());
    out.geometry = std::move(*g);
    return out;
}

// A bias point with fields varying over the nodes; `shift` moves the potential and bands (V, eV) and
// scales n by 10^shift.
results::BiasPoint point(const mesh::Mesh& m, std::vector<double> bias, double shift) {
    results::BiasPoint p;
    p.bias_V = std::move(bias);
    for (const Point& x : m.points()) {
        const double s = x[0] * 1e4 + 2.0 * x[1] * 1e4;  // 0 .. 4
        p.fields.potential_V.push_back(0.1 * s + shift);
        p.fields.n_cm3.push_back(1e10 * std::pow(10.0, s + shift));
        p.fields.p_cm3.push_back(1e5 * std::pow(10.0, -s));
        p.bands.conduction_eV.push_back(0.56 - 0.1 * s - shift);
        p.bands.valence_eV.push_back(-0.56 - 0.1 * s - shift);
        p.bands.electron_fermi_eV.push_back(0.2 * s);
        p.bands.hole_fermi_eV.push_back(-0.2 * s + 2.0 * shift);
    }
    return p;
}

}  // namespace

TEST_CASE("render compare: the bias point both runs have", "[render]") {
    const Bar b = bar(1e16);
    const auto& m = b.device->mesh();
    const std::vector<results::BiasPoint> runs{point(m, {0.0, 0.1}, 0.0), point(m, {0.0, 0.2}, 0.0),
                                               point(m, {0.0, 0.2}, 0.0)};
    CHECK(matching_point(point(m, {0.0, 0.2}, 0.0), runs) == 1);  // the first of equal ones
    CHECK(matching_point(point(m, {0.0, 0.2 + 5e-10}, 0.0), runs) == 1);
    CHECK_FALSE(matching_point(point(m, {0.0, 0.2 + 2e-9}, 0.0), runs));
    CHECK_FALSE(matching_point(point(m, {0.1}, 0.0), runs));  // another number of contacts
    CHECK_FALSE(matching_point(point(m, {0.0}, 0.0), runs));  // equal as far as A goes
    CHECK(matching_point(point(m, {}, 0.0), std::vector<results::BiasPoint>{point(m, {}, 1.0)}) == 0);
    CHECK_FALSE(matching_point(point(m, {0.0, 0.1}, 0.0), {}));
}

TEST_CASE("render compare: B - A of a run with itself is zero (XG5)", "[render]") {
    const Bar b = bar(1e16);
    const results::BiasPoint a = point(b.device->mesh(), {0.0, 0.5}, 0.0);
    std::size_t compared = 0;
    for (int k = 0; k <= static_cast<int>(Field::last); ++k) {
        const Field f = static_cast<Field>(k);
        if (!field_values(*b.device, *b.geometry, a, f)) continue;  // not in this point (rates)
        auto d = difference_field(*b.device, *b.device, *b.geometry, a, a, f);
        REQUIRE(d.has_value());
        ++compared;
        for (const double v : d->values) CHECK((v == 0.0 || std::isnan(v)));
        const ColorScale s = difference_scale(d->values);
        CHECK(s.lo() == -1.0);
        CHECK(s.hi() == 1.0);
        CHECK(s.position(0.0) == 0.5);
    }
    CHECK(compared >= 11);  // potential, n, p, the bands, T, the doping, |E| ...
}

TEST_CASE("render compare: B - A of two runs is the nodal difference in mV, meV and decades (XG5)", "[render]") {
    const Bar ba = bar(1e16), bb = bar(4e16);
    const auto& m = ba.device->mesh();
    const results::BiasPoint a = point(m, {0.0, 0.5}, 0.0), b = point(m, {0.0, 0.5}, 0.25);
    const std::size_t n = m.node_count();

    auto psi = difference_field(*ba.device, *bb.device, *ba.geometry, a, b, Field::potential);
    REQUIRE(psi.has_value());
    CHECK(psi->title == "\xCE\x94\xCF\x88 (mV)");
    REQUIRE(psi->values.size() == n);
    for (std::size_t i = 0; i < n; ++i) {
        CHECK_THAT(psi->values[i], WithinAbs((b.fields.potential_V[i] - a.fields.potential_V[i]) * 1e3, 1e-9));
        CHECK_THAT(psi->values[i], WithinAbs(250.0, 1e-9));
    }
    auto ec = difference_field(*ba.device, *bb.device, *ba.geometry, a, b, Field::conduction_band);
    REQUIRE(ec.has_value());
    CHECK(ec->title == "\xCE\x94" "Ec (meV)");
    for (const double v : ec->values) CHECK_THAT(v, WithinAbs(-250.0, 1e-9));
    auto efp = difference_field(*ba.device, *bb.device, *ba.geometry, a, b, Field::hole_fermi);
    REQUIRE(efp.has_value());
    for (const double v : efp->values) CHECK_THAT(v, WithinAbs(500.0, 1e-9));
    // Densities in decades: n is 10^0.25 times run A's everywhere, p unchanged.
    auto dn = difference_field(*ba.device, *bb.device, *ba.geometry, a, b, Field::electrons);
    REQUIRE(dn.has_value());
    CHECK(dn->title == "\xCE\x94log\xE2\x82\x81\xE2\x82\x80 n (decades)");
    for (const double v : dn->values) CHECK_THAT(v, WithinAbs(0.25, 1e-12));
    auto dp = difference_field(*ba.device, *bb.device, *ba.geometry, a, b, Field::holes);
    REQUIRE(dp.has_value());
    for (const double v : dp->values) CHECK(v == 0.0);
    // The doping of each run's own device: N_D four times run A's, N_A equal; the net doping linear.
    auto nd = difference_field(*ba.device, *bb.device, *ba.geometry, a, b, Field::donors);
    REQUIRE(nd.has_value());
    for (const double v : nd->values) CHECK_THAT(v, WithinAbs(std::log10(4.0), 1e-12));
    auto net = difference_field(*ba.device, *bb.device, *ba.geometry, a, b, Field::net_doping);
    REQUIRE(net.has_value());
    CHECK(net->title == "\xCE\x94(ND \xE2\x88\x92 NA) (cm\xE2\x81\xBB\xC2\xB3)");
    for (std::size_t i = 0; i < n; ++i) {
        CHECK_THAT(net->values[i], WithinRel(3.0 * ba.device->donors()[i], 1e-12));
    }
    // Not positive in either run: no density difference there.
    results::BiasPoint zero = b;
    zero.fields.n_cm3[3] = 0.0;
    auto dz = difference_field(*ba.device, *bb.device, *ba.geometry, a, zero, Field::electrons);
    REQUIRE(dz.has_value());
    CHECK(std::isnan(dz->values[3]));
    CHECK_THAT(dz->values[4], WithinAbs(0.25, 1e-12));
    // Symmetric about 0 over the largest |difference|, the NaN left out.
    const ColorScale s = difference_scale(std::vector<double>{-3.0, 1.5, std::nan(""), 2.0});
    CHECK(s.lo() == -3.0);
    CHECK(s.hi() == 3.0);
    CHECK(difference_scale(dn->values).hi() == 0.25);
    // Fields not one per node: field_values' error.
    results::BiasPoint short_ = b;
    short_.fields.potential_V.pop_back();
    CHECK_FALSE(difference_field(*ba.device, *bb.device, *ba.geometry, a, short_, Field::potential).has_value());
}
