// Unit 26n slice 4 (DT8): the models' domain guardrails (the 3D modeller's and the 2D editor's).
#include <catch2/catch_test_macros.hpp>

#include <string>

#include "NiTCAD/app/geometry.hpp"
#include "NiTCAD/app/solid.hpp"

using namespace NiTCAD::app;
using solid::Box3;
using solid::P3;

namespace {

constexpr solid::Coord um = solid::per_um;

const std::string layered3d = R"([mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
max_size = 0.15um

[region si]
material = silicon

[region ox]
material = sio2
box = 0 1um 0 0.6um 0.4um 0.6um
)";

}  // namespace

TEST_CASE("app guardrails: a solid past the domain is named with its face and distance", "[app][guardrails]") {
    auto m = solid::read_model(layered3d);
    REQUIRE(m);
    CHECK(solid::outside_domain(*m).empty());
    solid::Model t = *m;
    REQUIRE(solid::move_solid(t, {1, 0}, {0, 0, um / 10}));
    CHECK(solid::outside_domain(t) == "region 'ox' reaches 0.1um past the domain's z_max face");
    t = *m;
    REQUIRE(solid::move_solid(t, {1, 0}, {-um / 20, 0, 0}));
    CHECK(solid::outside_domain(t) == "region 'ox' reaches 0.05um past the domain's x_min face");
    t = *m;
    REQUIRE(solid::add_box(t, {{0, 0, 0}, {um / 2, um / 2, um / 2}}, "silicon") >= 0);
    CHECK(solid::outside_domain(t).empty());  // on the faces is inside
}

TEST_CASE("app guardrails: clamps hold points, offsets and scalars at the domain", "[app][guardrails]") {
    const Box3 d{{0, 0, 0}, {um, 6000, 6000}};
    CHECK(solid::clamp_to(d, {-5, 7000, 3000}) == P3{0, 6000, 3000});
    const solid::Solid box{solid::SolidKind::box, {{1000, 1000, 1000}, {3000, 2000, 2000}}, 2, 0, 0, {}};
    CHECK(solid::clamp_offset(d, {box}, {9000, -5000, 100}) == P3{um - 3000, -1000, 100});
    CHECK(solid::clamp_offset(d, {box}, {-1000, 4000, 4000}) == P3{-1000, 4000, 4000});
    CHECK(solid::clamp_offset(d, {}, {9, 9, 9}) == P3{9, 9, 9});
    // A solid already past a face is not pulled back, nor sent farther out.
    const solid::Solid out{solid::SolidKind::box, {{-500, 1000, 1000}, {3000, 2000, 2000}}, 2, 0, 0, {}};
    CHECK(solid::clamp_offset(d, {out}, {-10, 0, 0}).x == 0);
    CHECK(solid::clamp_offset(d, {out}, {300, 0, 0}).x == 300);
    // The bisection: the last value inside, on either side of 0.
    CHECK(solid::clamp_scalar(5000, [](solid::Coord v) { return v <= 1234; }) == 1234);
    CHECK(solid::clamp_scalar(-5000, [](solid::Coord v) { return v >= -77; }) == -77);
    CHECK(solid::clamp_scalar(40, [](solid::Coord) { return true; }) == 40);
    CHECK(solid::clamp_scalar(40, [](solid::Coord v) { return v == 0; }) == 0);
    CHECK(solid::clamp_scalar(0, [](solid::Coord) { return false; }) == 0);
}

TEST_CASE("app guardrails: the 2D editor's shapes held only by a domain that is given", "[app][guardrails]") {
    const std::string layered2d = "[mesh]\ntype = generate2d\nlayered = on\ndomain = 0 2um 0 1um\n\n[region p]\nmaterial = silicon\n\n"
                                  "[region ox]\nmaterial = sio2\nbox = 0 2um 0.8um 1um\n";
    auto g = geo::read_geometry(layered2d);
    REQUIRE(g);
    CHECK(geo::domain_holds(*g));
    CHECK(geo::outside_domain(*g).empty());
    g->regions[1].shapes[0].points[1].y = 12000;
    CHECK(geo::outside_domain(*g) == "region 'ox' reaches 0.2um past the domain's y_max side");
    // Without a domain line the domain is the shapes' bounds: nothing to hold them.
    auto free = geo::read_geometry("[mesh]\ntype = generate2d\n\n[region p]\nmaterial = silicon\nbox = 0 1um 0 1um\n");
    REQUIRE(free);
    CHECK_FALSE(geo::domain_holds(*free));
    CHECK(geo::outside_domain(*free).empty());
}
