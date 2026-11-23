// Unlock shapes (Unit 26n slice 2, DT5, DC8): a 3D tensor or extrude3d deck becomes a layered generate3d deck the
// modeller edits, its regions where they were.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "NiTCAD/app/geometry.hpp"
#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/app/unlock.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/mesh/locate.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
using app::unlock::Kind;

namespace {

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

const std::string tensor_deck = R"([parameters]
L = 1um  # Length
H = 0.6um
W = 0.4um
t = 0.2um  # Oxide
y_ox = $H - $t

[mesh]
type = tensor
x = linspace(0, $L, 11)
y = linspace(0, $y_ox, 5) + linspace($y_ox, $H, 3)
z = 0, 0.2um, $W
straddle y = $y_ox 1nm

[region si]
material = silicon

[region ox]
material = sio2
box = 0 $L  $y_ox $H  0 $W   # the oxide

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = y_min

[contact gate]
kind = electrode
patch = y_max
box = 0.2um 0.8um  $H $H  0 $W
electrode = metal
work_function = 4.1

[equilibrium]
)";

// The region a point lies in by the deck's own boxes: the last box holding it in deck order, else the background
// (the tensor deck's rule, away from box faces); -1 within `margin` cm of a face (left out of the comparison).
int painted_region(const deck::Deck& substituted, const mesh::Point& p, double margin) {
    int region = -1, background = -1, index = 0;
    bool near = false;
    for (const auto& s : substituted.sections) {
        if (s.kind != "region") continue;
        bool any = false;
        for (const auto& e : s.entries) {
            if (e.key != "box") continue;
            any = true;
            const auto w = deck::value_words(e.value);
            double b[6];
            for (int k = 0; k < 6; ++k) b[k] = static_cast<double>(*app::geo::parse_length(w[static_cast<std::size_t>(k)])) * app::geo::quantum_cm;
            bool in = true;
            for (int a = 0; a < 3; ++a) {
                in = in && p[static_cast<std::size_t>(a)] >= b[2 * a] && p[static_cast<std::size_t>(a)] <= b[2 * a + 1];
            }
            for (int a = 0; a < 3; ++a) {
                for (int j = 0; j < 2; ++j) near = near || std::abs(p[static_cast<std::size_t>(a)] - b[2 * a + j]) < margin;
            }
            if (in) region = index;
        }
        if (!any) background = index;
        ++index;
    }
    if (near) return -1;
    return region >= 0 ? region : background;
}

}  // namespace

TEST_CASE("app unlock: the deck kinds it takes", "[app][unlock]") {
    CHECK(app::unlock::kind_of(tensor_deck) == Kind::tensor);
    CHECK(app::unlock::kind_of("[mesh]\ntype = tensor\nx = 0 1um\ny = 0 1um\n") == Kind::none);
    CHECK(app::unlock::kind_of("[mesh]\ntype = extrude3d\n") == Kind::extruded);
    CHECK(app::unlock::kind_of("[mesh]\ntype = generate3d\nlayered = on\n") == Kind::none);
    CHECK(app::unlock::kind_of("[mesh") == Kind::none);
    CHECK(!app::unlock::unlock("[mesh]\ntype = tensor\nx = 0 1um\ny = 0 1um\n"));
}

TEST_CASE("app unlock: a tensor deck becomes a layered deck, references kept, the rest byte for byte", "[app][unlock]") {
    REQUIRE(deck::build(*deck::parse(tensor_deck)));  // the deck itself builds
    const auto u = app::unlock::unlock(tensor_deck);
    REQUIRE(u);
    const std::string& t = u->text;
    INFO(t);
    CHECK(app::unlock::kind_of(t) == Kind::none);
    CHECK(t.find("type = generate3d\nlayered = on\ndomain = 0 $L 0 $H 0 $W\n") != std::string::npos);
    CHECK(t.find("patch body = 0 0 0  $L 0 $W\n") != std::string::npos);
    CHECK(t.find("patch gate = 0.2um $H 0  0.8um $H $W\n") != std::string::npos);
    CHECK(t.find("max_size = 0.2um\n") != std::string::npos);  // the widest spacing: z's 200 nm
    CHECK(t.find("interface_resolution = ") != std::string::npos);
    CHECK(t.find("straddle") == std::string::npos);
    CHECK(t.find("linspace") == std::string::npos);
    CHECK(t.find("box = 0 $L  $y_ox $H  0 $W   # the oxide\n") != std::string::npos);  // the region box kept as written
    CHECK(t.find("box = 0.2um") == std::string::npos);  // the contact's box is its patch now
    CHECK(t.find("[contact gate]\nkind = electrode\npatch = gate\nelectrode = metal\n") != std::string::npos);
    CHECK(t.find("[parameters]\nL = 1um  # Length\n") == 0);
    CHECK(t.find("[doping]\nprofile = acceptor constant 1e16\n") != std::string::npos);
    CHECK(u->nodes > 0.0);
    // The modeller reads it, its parameters kept in the solids.
    const auto m = app::solid::read_model(t);
    REQUIRE(m);
    REQUIRE(m->mesh == app::solid::MeshKind::layered);
    REQUIRE(m->regions.size() == 2);
    CHECK(m->regions[0].solids.empty());
    REQUIRE(m->regions[1].solids.size() == 1);
    CHECK(m->regions[1].solids[0].refs.size() == 6);
    CHECK(m->regions[1].solids[0].refs[3] == "$H");
    REQUIRE(m->contacts.size() == 2);
    CHECK(m->contacts[1].editable);
    CHECK(m->contacts[1].box.lo.x == 2000);
    // And it builds.
    const auto d = deck::parse(t);
    REQUIRE(d);
    const auto plan = deck::build(*d);
    INFO((plan ? std::string() : plan.error().message));
    CHECK(plan);
}

TEST_CASE("app unlock: line endings, boxes past the domain, refusals", "[app][unlock]") {
    std::string crlf;
    for (const char c : tensor_deck) {
        if (c == '\n') crlf += '\r';
        crlf += c;
    }
    const auto u = app::unlock::unlock(crlf);
    REQUIRE(u);
    CHECK(u->text.find("layered = on\r\n") != std::string::npos);
    CHECK(u->text.find("\n\n") == std::string::npos);

    std::string past = tensor_deck;
    past.replace(past.find("box = 0 $L  $y_ox"), 10, "box = -1um 2um");
    const auto p = app::unlock::unlock(past);
    REQUIRE(p);
    CHECK(p->text.find("box = 0 $L $y_ox $H 0 $W\n") != std::string::npos);  // cut to the domain (its comment goes)
    CHECK(p->note.find("cut to the domain") != std::string::npos);

    std::string twice = tensor_deck;
    twice.replace(twice.find("electrode = metal"), 0, "box = 0 0.1um  $H $H  0 $W\n");
    const auto two = app::unlock::unlock(twice);
    REQUIRE(!two);
    CHECK(two.error().find("2 boxes") != std::string::npos);

    const auto adapt = app::unlock::unlock(tensor_deck + "\n[adapt]\nmax_passes = 2\n");
    REQUIRE(!adapt);
    CHECK(adapt.error().find("[adapt]") != std::string::npos);
}

TEST_CASE("app unlock: an extrude3d deck goes through the modeller's conversion", "[app][unlock]") {
    const std::string text = read_file(fs::path(NITCAD_EXAMPLES_DIR) / ".." / "decks" / "finfet_3d_reduced.deck");
    REQUIRE(!text.empty());
    REQUIRE(app::unlock::kind_of(text) == Kind::extruded);
    auto m = app::solid::read_model(text);
    REQUIRE(m);
    auto want = *m;
    const auto note = app::solid::convert_to_layered(want);
    const auto u = app::unlock::unlock(text);
    REQUIRE(bool(u) == bool(note));
    if (u) CHECK(u->text == *app::solid::write_model(text, want));
}

// DC8: each 3D library deck converts and builds with its regions where they were at 200 seeded points, or is
// refused with the reason (recorded per deck).
TEST_CASE("app unlock: the 3D library decks (DC8)", "[app][unlock][.unlock_library]") {
    for (const char* name : {"mosfet_3d.deck", "finfet_bulk_3d.deck", "finfet_soi_3d.deck"}) {
        const std::string text = read_file(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / name);
        REQUIRE(!text.empty());
        const auto t0 = std::chrono::steady_clock::now();
        const auto u = app::unlock::unlock(text);
        if (!u) {
            std::printf("DC8 %s: refused: %s\n", name, u.error().c_str());
            continue;
        }
        std::printf("DC8 %s: %s\n", name, u->note.c_str());
        std::fflush(stdout);
        const auto d = deck::parse(u->text);
        REQUIRE(d);
        const auto plan = deck::build(*d);
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (!plan) {
            std::printf("DC8 %s: refused by the build (%.1f s): %s\n", name, secs, plan.error().message.c_str());
            continue;
        }
        REQUIRE(plan->generated);
        const auto& g = *plan->generated;
        std::printf("DC8 %s: built, %zu nodes, %.1f s\n", name, g.node_region.size(), secs);
        const auto loc = mesh::PointLocator::create(g.simplices);
        REQUIRE(loc);
        const auto original = deck::substituted(*deck::parse(text));
        REQUIRE(original);
        const auto& lo = g.mesh.points().front();
        mesh::Point plo = lo, phi = lo;
        for (const auto& q : g.mesh.points()) {
            for (std::size_t a = 0; a < 3; ++a) {
                plo[a] = std::min(plo[a], q[a]);
                phi[a] = std::max(phi[a], q[a]);
            }
        }
        std::mt19937 rng(2026u);
        int compared = 0, tries = 0;
        while (compared < 200 && tries < 100000) {
            ++tries;
            mesh::Point p{};
            for (std::size_t a = 0; a < 3; ++a) p[a] = std::uniform_real_distribution<double>(plo[a], phi[a])(rng);
            const int want = painted_region(*original, p, 1e-9);
            if (want < 0) continue;
            const auto at = loc->locate(p);
            REQUIRE(at.inside);
            INFO(name << " at " << p[0] << " " << p[1] << " " << p[2]);
            CHECK(g.node_region[static_cast<std::size_t>(at.nodes[0])] == want);
            ++compared;
        }
        CHECK(compared == 200);
    }
}
