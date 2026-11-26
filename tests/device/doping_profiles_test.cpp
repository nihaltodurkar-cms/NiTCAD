// Analytic doping profiles (Unit 16; DECISIONS.md M15): closed-form values, windows, sums per
// species, insulator nodes, and the refusals.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::device;
using NiTCAD::base::ErrorCode;
using NiTCAD::mesh::Point;

namespace {
bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }
}  // namespace

TEST_CASE("doping profiles: depth shapes and windows", "[device][doping]") {
    DopingProfile c{.dopant = Dopant::acceptor, .peak_cm3 = 1e16};
    CHECK(evaluate(c, {3e-4, -1, 0}) == 1e16);
    DopingProfile g{.dopant = Dopant::donor, .peak_cm3 = 1e19, .shape = DepthShape::gaussian, .axis = 1,
                    .center = 1e-5, .width = 2e-6};
    CHECK(evaluate(g, {0, 1e-5, 0}) == 1e19);
    CHECK(close(evaluate(g, {0, 1.2e-5, 0}), 1e19 * std::exp(-0.5), 1e-14));  // one sigma away
    DopingProfile e{.dopant = Dopant::donor, .peak_cm3 = 1e20, .shape = DepthShape::erfc, .axis = 0, .center = 2e-5,
                    .width = 1e-5};
    CHECK(close(evaluate(e, {2e-5, 0, 0}), 0.5e20, 1e-14));  // the 50 % point
    CHECK(close(evaluate(e, {0, 0, 0}), 0.5e20 * std::erfc(-2.0), 1e-14));
    // A sharp window in y: 1 inside, 1/2 on its faces, 0 outside.
    DopingProfile box = c;
    box.lo = {0, 1e-5, 0};
    box.hi = {0, 3e-5, 0};
    box.lo[0] = -INFINITY;
    box.hi[0] = INFINITY;
    box.lo[2] = -INFINITY;
    box.hi[2] = INFINITY;
    CHECK(evaluate(box, {0, 2e-5, 0}) == 1e16);
    CHECK(evaluate(box, {0, 1e-5, 0}) == 0.5e16);
    CHECK(evaluate(box, {0, 4e-5, 0}) == 0.0);
    // With an erfc fall-off: 1/2 at the window's edge, nearly 1 deep inside.
    box.lateral = 1e-7;
    CHECK(close(evaluate(box, {0, 1e-5, 0}), 0.5e16, 1e-12));
    CHECK(close(evaluate(box, {0, 2e-5, 0}), 1e16, 1e-12));
}

TEST_CASE("doping profiles: sampled per species, zero in insulators", "[device][doping]") {
    const std::vector<double> x{0.0, 1e-5, 2e-5, 3e-5};
    const auto m = *NiTCAD::mesh::make_tensor_grid(x);
    const std::vector<Region> regions{{"oxide", NiTCAD::physics::silicon_dioxide()}, {"silicon", NiTCAD::physics::silicon()}};
    const std::vector<RegionId> node_region{0, 1, 1, 1};
    const std::vector<DopingProfile> profiles{
        {.dopant = Dopant::acceptor, .peak_cm3 = 1e16},
        {.dopant = Dopant::donor, .peak_cm3 = 1e18, .shape = DepthShape::erfc, .center = 2e-5, .width = 1e-6},
        {.dopant = Dopant::acceptor, .peak_cm3 = 1e15}};
    const auto d = sample(profiles, m, node_region, regions);
    REQUIRE(d.has_value());
    CHECK(d->acceptors[0] == 0.0);
    CHECK(d->donors[0] == 0.0);
    CHECK(d->acceptors[1] == 1.1e16);
    CHECK(close(d->donors[2], 0.5e18, 1e-14));
    CHECK(d->donors[3] < 1e-10 * 1e18);
}

TEST_CASE("doping profiles: refusals", "[device][doping]") {
    const std::vector<double> x{0.0, 1e-5};
    const auto m = *NiTCAD::mesh::make_tensor_grid(x);
    const std::vector<Region> regions{{"silicon", NiTCAD::physics::silicon()}};
    const std::vector<RegionId> nodes{0, 0};
    const auto code = [&](const DopingProfile& p) {
        const std::vector<DopingProfile> v{p};
        const auto r = sample(v, m, nodes, regions);
        return r.has_value() ? ErrorCode::cancelled : r.error().code;
    };
    CHECK(code({.peak_cm3 = -1}) == ErrorCode::invalid_input);
    CHECK(code({.peak_cm3 = NAN}) == ErrorCode::invalid_input);
    CHECK(code({.peak_cm3 = 1, .axis = 1}) == ErrorCode::invalid_input);  // a 1D mesh has only x
    CHECK(code({.peak_cm3 = 1, .shape = DepthShape::gaussian, .width = 0}) == ErrorCode::invalid_input);
    CHECK(code({.peak_cm3 = 1, .lateral = -1}) == ErrorCode::invalid_input);
// Unit 67 TAC8: a Gaussian or erfc window on its own depth axis is refused, no longer dropped.    DopingProfile gw{.peak_cm3 = 1, .shape = DepthShape::gaussian, .width = 1e-6};    gw.lo[0] = 0.0;    CHECK(code(gw) == ErrorCode::invalid_input);    gw.shape = DepthShape::constant;    CHECK(code(gw) == ErrorCode::cancelled);  // a constant's window on any axis is accepted
    DopingProfile w{.peak_cm3 = 1};
    w.lo[0] = 2;
    w.hi[0] = 1;
    CHECK(code(w) == ErrorCode::invalid_input);
    const std::vector<DopingProfile> ok{{.peak_cm3 = 1}};
    CHECK_FALSE(sample(ok, m, std::vector<RegionId>{0}, regions).has_value());
    CHECK_FALSE(sample(ok, m, std::vector<RegionId>{0, 5}, regions).has_value());
}

TEST_CASE("doping profiles: a constant is bounded on its depth axis too, a 1D step (Unit 67 TAC8)",
          "[device][doping]") {
    // Two constants windowed on x, the default depth axis: the 1D abrupt junction. Before Unit 67
    // the depth axis's window was dropped and both covered the whole line (a compensated resistor).
    const std::vector<double> x{0.0, 0.5e-4, 1e-4, 1.5e-4, 2e-4};
    const auto m = *NiTCAD::mesh::make_tensor_grid(x);
    const std::vector<Region> regions{{"silicon", NiTCAD::physics::silicon()}};
    const std::vector<RegionId> nodes(x.size(), 0);
    DopingProfile a{.dopant = Dopant::acceptor, .peak_cm3 = 1e17};
    a.lo[0] = -1.0;
    a.hi[0] = 1e-4;
    DopingProfile d{.dopant = Dopant::donor, .peak_cm3 = 1e17};
    d.lo[0] = 1e-4;
    d.hi[0] = 3e-4;
    const std::vector<DopingProfile> v{a, d};
    const auto s = sample(v, m, nodes, regions);
    REQUIRE(s.has_value());
    CHECK(s->acceptors == std::vector<double>{1e17, 1e17, 0.5e17, 0.0, 0.0});
    CHECK(s->donors == std::vector<double>{0.0, 0.0, 0.5e17, 1e17, 1e17});
    // A Gaussian keeps its window on the other axes only (unchanged).
    DopingProfile g{.peak_cm3 = 1e18, .shape = DepthShape::gaussian, .axis = 1, .center = 0.0, .width = 1e-5};
    g.hi[0] = 1e-4;
    CHECK(evaluate(g, {2e-4, 0, 0}) == 0.0);
    CHECK(evaluate(g, {0.5e-4, 0, 0}) == 1e18);
}
