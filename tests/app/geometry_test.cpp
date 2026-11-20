// The geometry editor's model (Unit 26k slice 1; DECISIONS.md GE3-GE6, GE10): the deck read and written by
// minimal edits, edits undone byte for byte, the write's time, snapping and hit tests at every dpi and zoom.
#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <format>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "NiTCAD/app/document.hpp"
#include "NiTCAD/app/geometry.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"

using namespace NiTCAD::app::geo;

namespace {

constexpr Coord um = per_um;

// The layered deck of the deck build test (a junction under an oxide), with comments and odd spacing to keep.
const std::string layered_deck = R"(# a layered device
[mesh]
type = generate2d
layered = on
domain = 0 2um 0 1um
patch source = 0 0, 0 0.6um
patch body = 0.2um 0, 0.8um 0   # the body contact
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon

[region n]
material = silicon
# the well
polygon = 1um 0, 2um 0, 2um 0.5um, \
          2um 0.8um, 1um 0.8um

[region oxide]
material = sio2
box = 0   2um  0.8um 1um

[doping]
profile = acceptor constant 1e16

[contact source]
kind = ohmic
patch = source

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";

// A tiled generated deck (its polygons tile the domain).
const std::string tiled_deck = R"([mesh]
type = generate2d
patch source = 0 0, 0 0.3um
max_size = 0.25um
interface_resolution = 20nm

[region low]
material = silicon
polygon = 0 0, 2um 0, 2um 0.5um, 0 0.5um

[region high]
material = sio2
polygon = 0 0.5um, 2um 0.5um, 2um 1um, 0 1um

[doping]
profile = acceptor constant 1e16

[contact source]
kind = ohmic
patch = source

[equilibrium]
)";

std::string read_file(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    std::stringstream s;
    s << in.rdbuf();
    return s.str();
}

std::string example(const char* name) { return read_file(std::filesystem::path(NITCAD_EXAMPLES_DIR) / name); }

const Rect mosfet_bounds{{0, 0}, {2 * um, 1 * um}};

Geometry must_read(const std::string& text, std::optional<Rect> bounds = {}) {
    auto g = read_geometry(text, bounds);
    if (!g) FAIL(g.error());
    return std::move(*g);
}

std::string must_write(const std::string& text, const Geometry& g) {
    auto t = write_geometry(text, g);
    if (!t) FAIL(t.error());
    return std::move(*t);
}

// The model without where its parts came from (a write and a read renumber them).
Geometry strip(Geometry g) {
    for (Region& r : g.regions) r.origin = -1;
    for (Contact& c : g.contacts) c.origin = -1;
    return g;
}

std::vector<std::string> lines_of(const std::string& t) {
    std::vector<std::string> out;
    std::istringstream in(t);
    for (std::string l; std::getline(in, l);) out.push_back(l);
    return out;
}

// Lines in one text and not the other (as multisets, both ways).
std::vector<std::string> changed_lines(const std::string& a, const std::string& b) {
    auto x = lines_of(a), y = lines_of(b);
    std::ranges::sort(x);
    std::ranges::sort(y);
    std::vector<std::string> out;
    std::ranges::set_symmetric_difference(x, y, std::back_inserter(out));
    return out;
}

bool geometry_line(std::string_view l) {
    while (!l.empty() && (l.back() == '\r' || l.back() == ' ')) l.remove_suffix(1);
    if (l.empty()) return true;
    for (const char* k : {"box =", "polygon =", "material =", "[region ", "[contact ", "patch", "layered =", "domain =", "kind =", "          "}) {
        if (l.starts_with(k)) return true;
    }
    return false;
}

// Whether every line `after` adds is a geometry line, and every line it drops one or a comment (a removed section's).
std::string only_geometry(const std::string& before, const std::string& after, bool removes = false) {
    auto x = lines_of(before), y = lines_of(after);
    std::ranges::sort(x);
    std::ranges::sort(y);
    std::vector<std::string> added, removed;
    std::ranges::set_difference(y, x, std::back_inserter(added));
    std::ranges::set_difference(x, y, std::back_inserter(removed));
    for (const auto& l : added) {
        if (!geometry_line(l)) return "added: " + l;
    }
    for (const auto& l : removed) {
        if (!removes && !geometry_line(l) && !l.starts_with("#")) return "removed: " + l;
    }
    return {};
}

}  // namespace

TEST_CASE("app geometry: lengths in 0.1 nm quanta, written as the deck writes them", "[app][geometry]") {
    CHECK(parse_length("0.95um") == 9500);
    CHECK(parse_length("20nm") == 200);
    CHECK(parse_length("1e-4") == 10000);  // cm
    CHECK(parse_length("2cm") == 200000000);
    CHECK(parse_length("-0.5um") == -5000);
    CHECK(parse_length("0") == 0);
    CHECK(!parse_length("abc"));
    CHECK(!parse_length("1umx"));
    CHECK(format_length(0) == "0");
    CHECK(format_length(9500) == "0.95um");
    CHECK(format_length(20000) == "2um");
    CHECK(format_length(-1) == "-0.0001um");
    CHECK(format_length(12345678) == "1234.5678um");
    std::mt19937_64 rng(7);
    std::uniform_int_distribution<Coord> q(-Coord{1} << 40, Coord{1} << 40);
    for (int k = 0; k < 10000; ++k) {
        const Coord c = q(rng);
        REQUIRE(parse_length(format_length(c)) == c);
    }
}

TEST_CASE("app geometry: the decks read; without edits the text comes back byte for byte (GE3)", "[app][geometry]") {
    const std::string mosfet = example("mosfet_2d.deck");
    SECTION("a tensor deck") {
        const Geometry g = must_read(mosfet, mosfet_bounds);
        REQUIRE(g.mesh == MeshKind::tensor);
        REQUIRE(g.regions.size() == 2);
        CHECK(g.regions[0].name == "silicon");
        CHECK(g.regions[0].shapes.empty());
        CHECK(g.background() == 0);
        REQUIRE(g.regions[1].shapes.size() == 1);
        CHECK(g.regions[1].shapes[0] == Shape{ShapeKind::box, {{0, 9500}, {2 * um, 1 * um}}});
        REQUIRE(g.contacts.size() == 4);
        CHECK(g.contacts[0].name == "source");
        CHECK(g.contacts[0].a == Pt{0, 6000});
        CHECK(g.contacts[0].b == Pt{0, 9400});
        CHECK(g.contacts[2].a == Pt{5000, um});
        CHECK(g.contacts[2].b == Pt{15000, um});
        CHECK(g.contacts[3].a == Pt{0, 0});  // the whole face
        CHECK(g.contacts[3].b == Pt{2 * um, 0});
        CHECK(must_write(mosfet, g) == mosfet);
    }
    SECTION("a layered deck") {
        const Geometry g = must_read(layered_deck);
        REQUIRE(g.mesh == MeshKind::layered);
        CHECK(g.domain_given);
        CHECK(g.domain == Rect{{0, 0}, {2 * um, um}});
        REQUIRE(g.regions.size() == 3);
        CHECK(g.regions[1].shapes[0].points.size() == 5);
        CHECK(g.contacts[1].a == Pt{2000, 0});
        CHECK(g.contacts[1].b == Pt{8000, 0});
        CHECK(must_write(layered_deck, g) == layered_deck);
    }
    SECTION("a tiled deck, CRLF and a last line without its end") {
        std::string crlf;
        for (const char c : tiled_deck) crlf += c == '\n' ? std::string("\r\n") : std::string(1, c);
        crlf += "# end";
        const Geometry g = must_read(crlf);
        CHECK(g.mesh == MeshKind::tiled);
        CHECK(g.domain == Rect{{0, 0}, {2 * um, um}});
        CHECK(must_write(crlf, g) == crlf);
    }
    SECTION("what the editor does not draw") {
        const Geometry one = must_read(example("pn_diode_1d.deck"));
        CHECK(one.mesh == MeshKind::unsupported);
        CHECK(!one.why.empty());
        CHECK(!write_geometry(example("pn_diode_1d.deck"), one));
        Geometry copy = one;
        CHECK(add_box(copy, {{0, 0}, {1, 1}}, "silicon") == -1);
        CHECK(must_read(example("mosfet_3d.deck")).mesh == MeshKind::unsupported);
        CHECK(!read_geometry("[mesh\n"));
    }
}

TEST_CASE("app geometry: an edit rewrites only geometry lines, and the deck still builds (GE3)", "[app][geometry]") {
    SECTION("layered: resize the oxide, reshape the well, add and move a box, move a contact") {
        Geometry g = must_read(layered_deck);
        REQUIRE(set_box(g, {2, 0}, {{0, 7500}, {2 * um, um}}));
        REQUIRE(set_vertex(g, {1, 0}, 4, {12000, 7500}));  // its neighbours follow
        CHECK(g.regions[1].shapes[0].points[0] == Pt{12000, 0});
        CHECK(g.regions[1].shapes[0].points[3] == Pt{2 * um, 7500});
        REQUIRE(add_box(g, {{3000, 5500}, {6000, 7500}}, "silicon") == 3);  // under the oxide
        CHECK(g.regions[3].name == "silicon_1");
        REQUIRE(move_shape(g, {3, 0}, 1000, 0));
        REQUIRE(set_contact(g, 1, {2000, 0}, {9000, 0}));
        const std::string out = must_write(layered_deck, g);
        for (const auto& l : changed_lines(layered_deck, out)) {
            INFO(l);
            CHECK(geometry_line(l));
        }
        CHECK(out.find("# the well\n") != std::string::npos);
        CHECK(out.find("# a layered device\n") == 0);
        CHECK(out.find("box = 0 2um 0.75um 1um") != std::string::npos);
        CHECK(out.find("box = 0.4um 0.7um 0.55um 0.75um") != std::string::npos);
        CHECK(out.find("patch body = 0.2um 0, 0.9um 0") != std::string::npos);
        CHECK(strip(must_read(out)) == strip(g));
        const auto d = NiTCAD::deck::parse(out, "t.deck");
        REQUIRE(d);
        const auto plan = NiTCAD::deck::build(*d);
        INFO((plan ? std::string() : plan.error().message));
        CHECK(plan);
    }
    SECTION("tiled: the first edit makes it layered; the deck builds") {
        Geometry g = must_read(tiled_deck);
        REQUIRE(set_vertex(g, {1, 0}, 0, {0, 4000}));  // the oxide down over the silicon: still covered
        CHECK(g.mesh == MeshKind::layered);
        const std::string out = must_write(tiled_deck, g);
        CHECK(out.find("type = generate2d\nlayered = on\n") != std::string::npos);
        for (const auto& l : changed_lines(tiled_deck, out)) CHECK(geometry_line(l));
        const auto d = NiTCAD::deck::parse(out, "t.deck");
        REQUIRE(d);
        CHECK(NiTCAD::deck::build(*d));
        // A slanted tiling is edited too since slice 4 (GD10): it becomes layered and builds.
        const std::string slanted_deck = R"([mesh]
type = generate2d
patch source = 0 0, 0 0.3um
max_size = 0.25um
interface_resolution = 20nm

[region a]
material = silicon
polygon = 0 0, 2um 0, 2um 0.3um, 0 0.7um

[region b]
material = sio2
polygon = 0 0.7um, 2um 0.3um, 2um 1um, 0 1um

[doping]
profile = acceptor constant 1e16

[contact source]
kind = ohmic
patch = source

[equilibrium]
)";
        Geometry slanted = must_read(slanted_deck);
        REQUIRE(insert_vertex(slanted, {1, 0}, 0, {um, 5000}));  // a vertex on the oxide's slanted edge
        CHECK(slanted.mesh == MeshKind::layered);
        const auto ds = NiTCAD::deck::parse(must_write(slanted_deck, slanted), "s.deck");
        REQUIRE(ds);
        const auto bs = NiTCAD::deck::build(*ds);
        if (!bs) FAIL(bs.error().message);
    }
    SECTION("tensor: boxes and contact boxes") {
        const std::string mosfet = example("mosfet_2d.deck");
        Geometry g = must_read(mosfet, mosfet_bounds);
        CHECK(add_polygon(g, {{0, 0}, {1, 0}, {1, 1}, {0, 1}}, "silicon") == -1);
        REQUIRE(set_box(g, {1, 0}, {{0, 9000}, {2 * um, um}}));
        REQUIRE(set_contact(g, 2, {6000, um}, {14000, um}));
        REQUIRE(set_contact(g, 3, {0, 0}, {um, 0}));       // a box added
        REQUIRE(set_contact(g, 0, {0, 0}, {0, um}));       // the whole face: its box kept, still read so
        REQUIRE(add_contact(g, "sub", "ohmic", {um, 0}, {2 * um, 0}) >= 0);
        const std::string out = must_write(mosfet, g);
        for (const auto& l : changed_lines(mosfet, out)) {
            INFO(l);
            CHECK(geometry_line(l));
        }
        CHECK(strip(must_read(out, mosfet_bounds)) == strip(g));
        std::string lf = out;
        std::erase(lf, '\r');
        CHECK(lf.find("[contact sub]\nkind = ohmic\npatch = y_min\nbox = 1um 2um 0 0\n") != std::string::npos);
        CHECK(out.find("\r\n[contact sub]\r\n") != std::string::npos);  // the deck's line ends kept
    }
    SECTION("regions reordered, renamed, deleted; contacts deleted and renamed") {
        Geometry g = must_read(layered_deck);
        REQUIRE(lower_region(g, 2));  // the oxide under the well
        REQUIRE(rename_region(g, 1, "ox"));
        CHECK(!rename_region(g, 0, "ox"));
        CHECK(!rename_region(g, 0, "9bad"));
        REQUIRE(set_material(g, 1, "hfo2"));
        REQUIRE(delete_contact(g, 0));
        REQUIRE(rename_contact(g, 0, "substrate"));
        std::string out = must_write(layered_deck, g);
        CHECK(strip(must_read(out)) == strip(g));
        CHECK(out.find("patch source") == std::string::npos);
        CHECK(out.find("[contact source]") == std::string::npos);
        CHECK(out.find("[region ox]\nmaterial = hfo2\nbox =") != std::string::npos);
        CHECK(out.find("\n\n\n") == std::string::npos);
        // Deleting a region's last shape deletes the region.
        Geometry h = must_read(out);
        REQUIRE(delete_shape(h, {2, 0}));
        CHECK(h.regions.size() == 2);
        out = must_write(out, h);
        CHECK(out.find("[region n]") == std::string::npos);
        CHECK(out.find("\n\n\n") == std::string::npos);
        CHECK(strip(must_read(out)) == strip(h));
    }
    SECTION("contacts stay on one side of the domain") {
        Geometry g = must_read(layered_deck);
        CHECK(add_contact(g, "gate", "ohmic", {0, um}, {um, um}) >= 0);
        CHECK(add_contact(g, "gate", "ohmic", {0, um}, {um, um}) == -1);
        CHECK(add_contact(g, "x", "ohmic", {0, 5000}, {um, 5000}) == -1);  // inside
        CHECK(add_contact(g, "y", "ohmic", {0, 0}, {0, 0}) == -1);
        CHECK(!set_contact(g, 0, {0, 0}, {um, um}));
        const std::string out = must_write(layered_deck, g);
        CHECK(out.find("patch gate = 0 1um, 1um 1um") != std::string::npos);
        CHECK(out.find("[contact gate]\nkind = ohmic\npatch = gate\n") != std::string::npos);
    }
}

TEST_CASE("app geometry: 200 random edits, undone one by one, give the deck back byte for byte (GE4, GE5)", "[app][geometry]") {
    for (const bool tensor : {false, true}) {
        const std::string start = tensor ? example("mosfet_2d.deck") : layered_deck;
        const std::optional<Rect> bounds = tensor ? std::optional(mosfet_bounds) : std::nullopt;
        std::mt19937 rng(tensor ? 11u : 5u);
        const auto pick = [&](int n) { return std::uniform_int_distribution<int>(0, n - 1)(rng); };
        const auto len = [&] { return static_cast<Coord>(pick(40) - 20) * 500; };
        std::string text = start;
        std::vector<std::string> texts{text};
        EditHistory h;
        int done = 0, tries = 0;
        double write_ms = 0.0;
        while (done < 200 && tries < 20000) {
            ++tries;
            Geometry g = must_read(text, bounds);
            const auto shapes = [&] {
                std::vector<ShapeRef> out;
                for (int r = 0; r < static_cast<int>(g.regions.size()); ++r) {
                    for (int s = 0; s < static_cast<int>(g.regions[static_cast<std::size_t>(r)].shapes.size()); ++s) out.push_back({r, s});
                }
                return out;
            }();
            const auto any_shape = [&]() -> ShapeRef { return shapes.empty() ? ShapeRef{} : shapes[static_cast<std::size_t>(pick(static_cast<int>(shapes.size())))]; };
            const int nr = static_cast<int>(g.regions.size()), nc = static_cast<int>(g.contacts.size());
            const auto side = [&](Pt& a, Pt& b) {
                const Rect d = g.domain;
                const Coord t0 = pick(10) * 1000, t1 = t0 + (1 + pick(5)) * 1000;
                switch (pick(4)) {
                case 0: a = {d.lo.x + t0 * (d.hi.x - d.lo.x) / 20000, d.lo.y}; b = {d.lo.x + t1 * (d.hi.x - d.lo.x) / 20000, d.lo.y}; break;
                case 1: a = {d.lo.x + t0 * (d.hi.x - d.lo.x) / 20000, d.hi.y}; b = {d.lo.x + t1 * (d.hi.x - d.lo.x) / 20000, d.hi.y}; break;
                case 2: a = {d.lo.x, d.lo.y + t0 * (d.hi.y - d.lo.y) / 20000}; b = {d.lo.x, d.lo.y + t1 * (d.hi.y - d.lo.y) / 20000}; break;
                default: a = {d.hi.x, d.lo.y + t0 * (d.hi.y - d.lo.y) / 20000}; b = {d.hi.x, d.lo.y + t1 * (d.hi.y - d.lo.y) / 20000}; break;
                }
            };
            bool ok = false, removes = false;  // removes: a whole section may go, its other lines with it
            switch (pick(16)) {
            case 0: ok = move_shape(g, any_shape(), len(), len()); break;
            case 1: {
                const Coord x = pick(20) * 1000, y = pick(10) * 1000;
                ok = set_box(g, any_shape(), {{x, y}, {x + 1000 + pick(10) * 500, y + 500 + pick(10) * 500}});
                break;
            }
            case 2: {
                const Coord x = pick(20) * 1000, y = pick(10) * 1000;
                ok = add_box(g, {{x, y}, {x + 2000, y + 1000}}, pick(2) ? "silicon" : "sio2") >= 0;
                break;
            }
            case 3: ok = removes = shapes.size() > 1 && delete_shape(g, any_shape()); break;
            case 4: ok = raise_region(g, pick(std::max(nr, 1))); break;
            case 5: ok = lower_region(g, pick(std::max(nr, 1))); break;
            case 6: ok = rename_region(g, pick(std::max(nr, 1)), "r" + std::to_string(pick(1000))); break;
            case 7: ok = set_material(g, pick(std::max(nr, 1)), pick(2) ? "silicon" : "hfo2"); break;
            case 8: {
                Pt a, b;
                side(a, b);
                ok = add_contact(g, "c" + std::to_string(pick(1000)), "ohmic", a, b) >= 0;
                break;
            }
            case 9: {
                Pt a, b;
                side(a, b);
                ok = nc > 0 && set_contact(g, pick(nc), a, b);
                break;
            }
            case 10: ok = removes = nc > 1 && delete_contact(g, pick(nc)); break;
            case 11: ok = nc > 0 && rename_contact(g, pick(nc), "k" + std::to_string(pick(1000))); break;
            case 12: ok = duplicate_shape(g, any_shape(), len(), len()) >= 0; break;
            case 13: ok = mirror_shape(g, any_shape(), pick(2) == 0); break;
            case 14: {
                const ShapeRef s = any_shape();
                ok = s.region >= 0 && g.regions[static_cast<std::size_t>(s.region)].shapes[static_cast<std::size_t>(s.shape)].kind == ShapeKind::polygon &&
                     set_vertex(g, s, pick(4), {pick(20) * 1000, pick(10) * 1000});
                break;
            }
            default: {
                const Coord x = pick(15) * 1000, y = pick(8) * 1000;
                ok = !tensor && add_polygon(g, {{x, y}, {x + 3000, y}, {x + 3000, y + 1000}, {x + 1000, y + 1000}, {x + 1000, y + 2000}, {x, y + 2000}},
                                            "silicon") >= 0;
                break;
            }
            }
            if (!ok) continue;
            const auto t0 = std::chrono::steady_clock::now();
            std::string next = must_write(text, g);
            write_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            {
                INFO("edit " << done << "\n" << next);
                REQUIRE(strip(must_read(next, bounds)) == strip(g));
                REQUIRE(only_geometry(text, next, removes) == "");
            }
            h.record(text);
            text = std::move(next);
            texts.push_back(text);
            ++done;
        }
        REQUIRE(done == 200);
        WARN("GE5: " << (tensor ? "tensor" : "layered") << " deck, mean write " << write_ms / done << " ms over " << done << " edits");
        CHECK(write_ms / done <= 2.0);  // GE5
        for (int k = done; k > 0; --k) {
            auto t = h.undo(text);
            REQUIRE(t);
            text = std::move(*t);
            REQUIRE(text == texts[static_cast<std::size_t>(k - 1)]);
        }
        CHECK(text == start);
        CHECK(!h.can_undo());
        for (int k = 1; k <= done; ++k) {
            auto t = h.redo(text);
            REQUIRE(t);
            text = std::move(*t);
        }
        CHECK(text == texts.back());
        CHECK(!h.can_redo());
    }
}

TEST_CASE("app geometry: a deck of 300 regions is written within 2 ms per edit (GE5)", "[app][geometry]") {
    std::string deck = "[mesh]\ntype = generate2d\nlayered = on\ndomain = 0 100um 0 100um\npatch a = 0 0, 0 1um\n\n[region bg]\nmaterial = silicon\n";
    for (int k = 0; k < 300; ++k) {
        deck += std::format("\n# block {}\n[region r{}]\nmaterial = sio2\nbox = {}um {}um {}um {}um\n", k, k, k % 20 * 5, k % 20 * 5 + 2, k / 20 * 5,
                            k / 20 * 5 + 2);
    }
    deck += "\n[contact a]\nkind = ohmic\npatch = a\n";
    Geometry g = must_read(deck);
    REQUIRE(g.regions.size() == 301);
    // GE5 bounds the write; the read that follows each edit is timed apart and only recorded.
    double write_ms = 0, read_ms = 0;
    std::string text = deck;
    for (int k = 0; k < 50; ++k) {
        REQUIRE(move_shape(g, {1 + k * 5, 0}, 1000, 0));
        const auto t0 = std::chrono::steady_clock::now();
        text = must_write(text, g);
        const auto t1 = std::chrono::steady_clock::now();
        g = must_read(text);
        write_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
        read_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    }
    const double ms = write_ms / 50;
    WARN("GE5: 300 regions, write " << ms << " ms, read " << read_ms / 50 << " ms per edit");
    // Release only, as KG6 (owner, 2026-11-19): the unoptimised Debug build is many times slower; there it is recorded.
#ifdef NDEBUG
    CHECK(ms <= 2.0);
#endif
}

TEST_CASE("app geometry: the grid, snapping and the Shift constraint", "[app][geometry]") {
    CHECK(grid_step(View{100.0, 0, 0, 96}) == 2000);   // 12 px = 0.12 µm -> 0.2 µm
    CHECK(grid_step(View{100.0, 0, 0, 192}) == 5000);  // 24 px = 0.24 µm -> 0.5 µm
    CHECK(grid_step(View{1e4, 0, 0, 96}) == 20);
    CHECK(grid_step(View{1e6, 0, 0, 96}) == 1);        // never below a quantum
    CHECK(grid_step(View{1e-2, 0, 0, 96}) == 2000LL * per_um);  // 12 px = 1200 µm
    for (const int dpi : {96, 144, 192}) {
        for (const double z : {1e-2, 1.0, 100.0, 1e4}) {
            const Coord s = grid_step(View{z, 0, 0, dpi});
            const double px = static_cast<double>(s) / per_um * z;
            INFO("dpi " << dpi << " zoom " << z);
            CHECK(px >= 12.0 * dpi / 96 - 1e-9);
            CHECK((s == 1 || px < 12.0 * dpi / 96 * 2.5 + 1e-9));
        }
    }
    Geometry g = must_read(layered_deck);
    for (const int dpi : {96, 144, 192}) {
        const View v{100.0, 10.0, 500.0, dpi};
        const double tol = 8.0 * dpi / 96;  // px
        const auto off = [&](Pt p, double px) { return v.model(v.sx(p.x) + px, v.sy(p.y)); };
        // A vertex (the oxide's corner) within 8 px, beyond it not.
        const Pt corner{0, 8000};
        CHECK(snap(g, v, off(corner, tol - 0.6)).kind == SnapKind::vertex);
        CHECK(snap(g, v, off(corner, tol - 0.6)).p == corner);
        CHECK(snap(g, v, off(corner, tol + 0.6)).kind != SnapKind::vertex);
        // An edge's midpoint, then a point on the edge, exact on the edge's line.
        CHECK(snap(g, v, {10000, um}).kind != SnapKind::none);
        const Snap e = snap(g, v, v.model(v.sx(3333), v.sy(8000) + tol / 2));
        CHECK(e.kind == SnapKind::edge);
        CHECK(e.p.y == 8000);
        const Snap m = snap(g, v, v.model(v.sx(15000) + 2, v.sy(8000) - 1));  // the well's top edge (1..2um) midpoint
        CHECK(m.kind == SnapKind::midpoint);
        CHECK(m.p == Pt{15000, 8000});
        // Off every edge: the grid, or the raw point without it.
        const Pt raw{3123, 3456};
        const Snap gr = snap(g, v, raw);
        CHECK(gr.kind == SnapKind::grid);
        CHECK(gr.p.x % grid_step(v) == 0);
        CHECK(gr.p.y % grid_step(v) == 0);
        CHECK(snap(g, v, raw, {}, false).p == raw);
        // A dragged shape does not snap to itself.
        CHECK(snap(g, v, off(corner, 1.0), ShapeRef{2, 0}).kind == SnapKind::edge);  // the domain's side, no longer the corner
    }
    CHECK(constrain({0, 0}, {10, 3}) == Pt{10, 0});
    CHECK(constrain({0, 0}, {-3, 10}) == Pt{0, 10});
    CHECK(drag_handle({{0, 0}, {10, 10}}, 4, {20, 30}) == Rect{{0, 0}, {20, 30}});
    CHECK(drag_handle({{0, 0}, {10, 10}}, 1, {99, 4}) == Rect{{0, 4}, {10, 10}});
    CHECK(drag_handle({{0, 0}, {10, 10}}, 7, {15, 99}) == Rect{{10, 0}, {15, 10}});  // past the other side: normalized
    CHECK(drag_handle({{0, 0}, {10, 10}}, 2, {20, -5}) == Rect{{0, -5}, {20, 10}});
}

TEST_CASE("app geometry: hit tests at 96, 144 and 192 dpi and zooms 1e-2 to 1e4 px/um (GE6)", "[app][geometry]") {
    for (const int dpi : {96, 144, 192}) {
        for (const double z : {1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4}) {
            INFO("dpi " << dpi << " zoom " << z);
            const View v{z, 40.0, 900.0, dpi};
            const Coord s = std::llround(300.0 / z * per_um);  // 300 px
            Geometry g;
            g.mesh = MeshKind::layered;
            g.domain = {{-s, -s}, {2 * s, 2 * s}};
            g.regions.push_back({"bg", "silicon", {}, -1});
            REQUIRE(add_box(g, {{0, 0}, {s, s}}, "sio2") == 1);
            REQUIRE(add_polygon(g, {{s / 2, s / 2}, {2 * s, s / 2}, {2 * s, 2 * s}, {s / 2, 2 * s}}, "silicon") == 2);
            REQUIRE(add_contact(g, "c", "ohmic", {-s, 0}, {-s, s}) == 0);
            const ShapeRef box{1, 0}, poly{2, 0};
            const auto at = [&](Pt p, double dx = 0, double dy = 0, std::optional<ShapeRef> sel = {}) {
                return hit_test(g, v, v.sx(p.x) + dx, v.sy(p.y) + dy, sel);
            };
            // The model point under a pixel comes back.
            for (const Pt p : {Pt{0, 0}, Pt{s, s}, Pt{-s, 2 * s}}) CHECK(v.model(v.sx(p.x), v.sy(p.y)) == p);
            // Handles of the selected box, within 4.5 px, the lower-left one 0, counter-clockwise.
            const double in = 4.0 * dpi / 96, out = 5.0 * dpi / 96;
            CHECK(at({0, 0}, -in, in, box) == Hit{HitKind::handle, box, 0, -1});
            CHECK(at({s / 2, 0}, 0, in, box) == Hit{HitKind::handle, box, 1, -1});
            CHECK(at({0, s}, -in, -in, box) == Hit{HitKind::handle, box, 6, -1});
            CHECK(at({0, 0}, -out, out, box).kind != HitKind::handle);
            // Unselected: no handles; the topmost shape (the polygon over the box's corner).
            CHECK(at({s * 3 / 4, s * 3 / 4}) == Hit{HitKind::shape, poly, -1, -1});
            CHECK(at({s / 4, s / 4}) == Hit{HitKind::shape, box, -1, -1});
            CHECK(at({-s / 2, -s / 2}) == Hit{HitKind::shape, {0, -1}, -1, -1});  // the background
            CHECK(at({3 * s, 3 * s}).kind == HitKind::none);                       // outside the domain
            // The selected polygon's vertices and edges.
            CHECK(at({2 * s, 2 * s}, in, in, poly) == Hit{HitKind::vertex, poly, 2, -1});
            CHECK(at({s, s / 2}, 0, 3.5 * dpi / 96, poly) == Hit{HitKind::edge, poly, 0, -1});
            CHECK(at({s, s / 2}, 0, out, poly).kind != HitKind::edge);
            // A contact's end before the contact; within 5 px.
            CHECK(at({-s, s}, 4.0 * dpi / 96, 0) == Hit{HitKind::contact_end, {}, 0, 1});
            CHECK(at({-s, s / 2}, 4.0 * dpi / 96, 0) == Hit{HitKind::contact, {}, 0, -1});
            CHECK(at({-s, s / 2}, 6.0 * dpi / 96, 0).kind == HitKind::shape);
        }
    }
}

TEST_CASE("app geometry: a contact's kind and a layered deck's domain (slice 2)", "[app][geometry]") {
    Geometry g = must_read(layered_deck);
    REQUIRE(set_contact_kind(g, 1, "schottky"));
    CHECK(!set_contact_kind(g, 1, "9x"));
    // Wider and taller: the contacts stay on their sides.
    REQUIRE(set_domain(g, {{0, 0}, {3 * um, 15000}}));
    std::string out = must_write(layered_deck, g);
    CHECK(out.find("[contact body]\nkind = schottky\npatch = body\n") != std::string::npos);
    CHECK(out.find("domain = 0 3um 0 1.5um\n") != std::string::npos);
    CHECK(out.find("patch body = 0.2um 0, 0.8um 0") != std::string::npos);
    CHECK(strip(must_read(out)) == strip(g));
    // Narrower: the body contact clipped to it; then too narrow for it: refused, nothing changed.
    REQUIRE(set_domain(g, {{0, 0}, {5000, um}}));
    CHECK(g.contacts[1].b == Pt{5000, 0});
    const Geometry before = g;
    CHECK(!set_domain(g, {{0, 0}, {1000, um}}));
    CHECK(g == before);
    // Only a layered deck's domain is drawn by its line.
    Geometry t = must_read(example("mosfet_2d.deck"), mosfet_bounds);
    CHECK(!set_domain(t, {{0, 0}, {um, um}}));
}

// --- slice 3: live validation (GD9, GE7) ------------------------------------------------------------------------

namespace {

// The deck as the live check builds it, and its first problem placed on the model.
struct Checked {
    std::string text;
    NiTCAD::app::CheckResult result;
    Problem problem;
};
Checked check_layout(const Geometry& g, const std::string& base = layered_deck) {
    Checked c;
    c.text = must_write(base, g);
    c.result = NiTCAD::app::check_deck(c.text, {});
    if (!c.result.diagnostics.empty()) {
        const auto& d = c.result.diagnostics.front();
        c.problem = locate_problem(must_read(c.text), c.text, d.line, d.message);
    }
    return c;
}

int line_of(const std::string& text, std::string_view what) {
    const auto at = text.find(what);
    REQUIRE(at != std::string::npos);
    return 1 + static_cast<int>(std::count(text.begin(), text.begin() + static_cast<std::ptrdiff_t>(at), '\n'));
}

bool has_marker(const Problem& p, Marker m) { return std::ranges::find(p.markers, m) != p.markers.end(); }

}  // namespace

TEST_CASE("app geometry: an invalid layout's problem is the build's message, placed where it lies (GE7)", "[app][geometry]") {
    const Geometry base = must_read(layered_deck);
    REQUIRE(check_layout(base).result.diagnostics.empty());
    using K = Marker::Kind;

    SECTION("two contacts overlapping along the bottom: both, the error at the later patch's line") {
        Geometry g = base;
        REQUIRE(add_contact(g, "c2", "ohmic", {5000, 0}, {um, 0}) == 2);
        const Checked c = check_layout(g);
        REQUIRE(c.result.diagnostics.size() == 1);
        const auto& d = c.result.diagnostics.front();
        INFO(d.message);
        CHECK(d.message.find("patches overlap") != std::string::npos);
        CHECK(d.line == line_of(c.text, "patch c2 ="));
        CHECK(c.problem.message == d.message);
        CHECK(has_marker(c.problem, {K::contact, {}, 2, {}}));
        CHECK(has_marker(c.problem, {K::contact, {}, 1, {}}));
        CHECK(c.problem.markers.size() == 2);
    }
    SECTION("a contact's end too near the junction under it: that contact") {
        Geometry g = base;
        REQUIRE(add_contact(g, "near", "ohmic", {10200, 0}, {15000, 0}) == 2);
        const Checked c = check_layout(g);
        REQUIRE(c.result.diagnostics.size() == 1);
        const auto& d = c.result.diagnostics.front();
        INFO(d.message);
        CHECK(d.message.find("a patch end within 3 s of an interface junction") != std::string::npos);
        CHECK(d.line == line_of(c.text, "patch near ="));
        CHECK(c.problem.message == d.message);
        CHECK(c.problem.markers == std::vector<Marker>{{K::contact, {}, 2, {}}});
    }
    SECTION("a box too small for the interface resolution: the points the mesher names") {
        Geometry g = base;
        REQUIRE(add_box(g, {{15000, 5000}, {16000, 6000}}, "sio2") >= 0);
        const Checked c = check_layout(g);
        REQUIRE(c.result.diagnostics.size() == 1);
        const auto& d = c.result.diagnostics.front();
        INFO(d.message);
        CHECK(d.line == line_of(c.text, "[mesh]"));
        CHECK(c.problem.message == d.message);
        REQUIRE_FALSE(c.problem.markers.empty());
        for (const Marker& m : c.problem.markers) {
            CHECK(m.kind == K::point);
            // On the small box's outline (the message's six decimals in cm: within 0.01 µm).
            CHECK(m.at.x >= 15000 - 100);
            CHECK(m.at.x <= 16000 + 100);
            CHECK(m.at.y >= 5000 - 100);
            CHECK(m.at.y <= 6000 + 100);
        }
    }
    SECTION("an unknown material: its region") {
        Geometry g = base;
        REQUIRE(set_material(g, 2, "unobtainium"));
        const Checked c = check_layout(g);
        REQUIRE(c.result.diagnostics.size() == 1);
        const auto& d = c.result.diagnostics.front();
        INFO(d.message);
        CHECK(c.problem.message == d.message);
        CHECK(c.problem.markers == std::vector<Marker>{{K::shape, {2, -1}, -1, {}}});
    }
    SECTION("the background's material: the domain") {
        Geometry g = base;
        REQUIRE(set_material(g, 0, "unobtainium"));
        const Checked c = check_layout(g);
        REQUIRE(c.result.diagnostics.size() == 1);
        CHECK(c.problem.markers == std::vector<Marker>{{K::domain, {0, -1}, -1, {}}});
    }
    SECTION("a contact's unknown kind: that contact") {
        Geometry g = base;
        REQUIRE(set_contact_kind(g, 0, "omic"));
        const Checked c = check_layout(g);
        REQUIRE(c.result.diagnostics.size() == 1);
        const auto& d = c.result.diagnostics.front();
        INFO(d.message);
        CHECK(d.line == line_of(c.text, "kind = omic"));
        CHECK(c.problem.markers == std::vector<Marker>{{K::contact, {}, 0, {}}});
    }
    SECTION("a shape the build refuses (no height): that shape, at its line") {
        std::string text = layered_deck;
        text.replace(text.find("box = 0   2um  0.8um 1um"), 24, "box = 0   2um  0.8um 0.8um");
        const auto r = NiTCAD::app::check_deck(text, {});
        REQUIRE(r.diagnostics.size() == 1);
        INFO(r.diagnostics[0].message);
        CHECK(r.diagnostics[0].line == line_of(text, "box = 0   2um  0.8um 0.8um"));
        const Problem p = locate_problem(must_read(text), text, r.diagnostics[0].line, r.diagnostics[0].message);
        CHECK(p.markers == std::vector<Marker>{{K::shape, {2, 0}, -1, {}}});
    }
    SECTION("a problem outside the geometry (a doping profile): no marker, the message kept") {
        std::string text = layered_deck;
        text.replace(text.find("constant 1e16"), 13, "constant -1");
        const auto r = NiTCAD::app::check_deck(text, {});
        REQUIRE(r.diagnostics.size() == 1);
        const Problem p = locate_problem(must_read(text), text, r.diagnostics[0].line, r.diagnostics[0].message);
        CHECK(p.message == r.diagnostics[0].message);
        CHECK(p.markers.empty());
    }
}

TEST_CASE("app geometry: valid layouts build (GE7)", "[app][geometry]") {
    Geometry g = must_read(layered_deck);
    // A moved well, a new oxide layer over it, a contact on the top's silicon-free edge away from the rest.
    REQUIRE(move_shape(g, {1, 0}, -5000, 0));
    REQUIRE(add_box(g, {{0, 9000}, {2 * um, um}}, "sio2") >= 0);
    REQUIRE(add_contact(g, "gate", "electrode", {5000, um}, {15000, um}) == 2);
    REQUIRE(set_contact(g, 1, {500, 0}, {3000, 0}));
    for (int step = 0; step < 4; ++step) {
        if (step == 1) REQUIRE(raise_region(g, 1));
        if (step == 2) REQUIRE(rename_region(g, 3, "cap"));
        if (step == 3) REQUIRE(set_contact(g, 2, {4000, um}, {16000, um}));  // the gate wider
        std::string text = must_write(layered_deck, g);
        // The gate needs its electrode's keys.
        if (const auto at = text.find("kind = electrode\n"); at != std::string::npos) {
            text.insert(at + 17, "electrode = metal\nwork_function = 4.1\n");
        }
        const auto r = NiTCAD::app::check_deck(text, {});
        INFO("step " << step << "\n" << text);
        for (const auto& d : r.diagnostics) UNSCOPED_INFO(d.message);
        CHECK(r.diagnostics.empty());
        CHECK(r.plan.has_value());
    }
}

TEST_CASE("app geometry: slanted polygons drawn and edited exactly, and the deck builds (slice 4, GD10)", "[app][geometry]") {
    Geometry g = must_read(layered_deck);
    // A trapezoid (a gate with slanted sidewalls) and a triangle are taken; a crossing, a flat or a doubled outline is not.
    const int gate = add_polygon(g, {{6000, 8000}, {14000, 8000}, {13000, um}, {7000, um}}, "sio2");
    REQUIRE(gate == 3);
    CHECK(add_polygon(g, {{0, 0}, {um, um}, {um, 0}, {0, um}}, "silicon") == -1);         // a bow tie
    CHECK(add_polygon(g, {{0, 0}, {um, 0}, {2 * um, 0}}, "silicon") == -1);               // no area
    CHECK(add_polygon(g, {{0, 0}, {um, 0}, {um, 0}, {0, um}}, "silicon") == -1);          // an empty edge
    CHECK(add_polygon(g, {{0, 0}, {um, 0}, {um / 2, 0}, {0, um}}, "silicon") == -1);      // an edge turning back
    CHECK(add_polygon(g, {{0, 0}, {um, 0}}, "silicon") == -1);
    CHECK(add_polygon(g, {{0, 0}, {2'000'000'000, 0}, {0, um}}, "silicon") == -1);       // beyond 10 cm
    // A vertex on the slanted sidewall, exactly on it (the cross product exact); one off it refused.
    CHECK(insert_vertex(g, {gate, 0}, 1, {13500, 9000}));
    CHECK_FALSE(insert_vertex(g, {gate, 0}, 2, {12000, 9600}));
    CHECK_FALSE(insert_vertex(g, {gate, 0}, 2, {13000, um}));                            // an end of the edge
    CHECK(g.regions[3].shapes[0].points.size() == 5);
    CHECK(remove_vertex(g, {gate, 0}, 2));
    REQUIRE(g.regions[3].shapes[0].points.size() == 4);
    // A triangle is taken and kept (built apart: its corners are below the mesher's 60 degrees); two vertices are not.
    Geometry h = must_read(layered_deck);
    REQUIRE(add_polygon(h, {{2000, 2000}, {6000, 2000}, {2000, 5000}, {2000, 4000}}, "silicon") == 3);
    CHECK(remove_vertex(h, {3, 0}, 3));
    CHECK_FALSE(remove_vertex(h, {3, 0}, 0));
    CHECK(h.regions[3].shapes[0].points.size() == 3);
    // Alone, a vertex moves by itself; following, a rectilinear one keeps its horizontal and vertical edges.
    REQUIRE(set_vertex(g, {1, 0}, 4, {11000, 8000}, false));                                 // the well's corner: a slant
    CHECK(g.regions[1].shapes[0].points[3] == Pt{2 * um, 8000});
    CHECK(g.regions[1].shapes[0].points[0] == Pt{um, 0});                                    // alone: not followed
    REQUIRE(set_vertex(g, {1, 0}, 0, {12000, 0}, true));
    CHECK(g.regions[1].shapes[0].points[1] == Pt{2 * um, 0});                               // along y = 0: not followed (x differs)
    CHECK_FALSE(set_vertex(g, {1, 0}, 4, {2 * um, 8000}, false));                            // onto its neighbour
    const std::string out = must_write(layered_deck, g);
    for (const auto& l : changed_lines(layered_deck, out)) CHECK(geometry_line(l));
    const auto d = NiTCAD::deck::parse(out, "t.deck");
    REQUIRE(d);
    const auto b = NiTCAD::deck::build(*d);
    if (!b) FAIL(b.error().message);
}
