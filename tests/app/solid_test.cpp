// The 3D solid modeller's model (Unit 26l slice 1; DECISIONS.md SD3, SE3-SE5, SE10): the deck read and written by
// minimal edits, extrude3d decks read as prisms and converted, edits undone byte for byte, the write's time, snapping.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cstdio>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <numbers>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "NiTCAD/app/modeller.hpp"
#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "fixtures.hpp"

using namespace NiTCAD::app::solid;

namespace {

constexpr Coord um = per_um;

// A layered deck with a prism along x written over two lines, comments and odd spacing to keep.
const std::string well_deck = R"(# a layered 3D stack
[mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
patch gate = 0.3um 0 0.6um  0.7um 0.6um 0.6um   # the gate on top
patch body = 0 0 0  1um 0.6um 0
max_size = 0.15um
interface_resolution = 0.02um

[region si]
material = silicon

[region ox]
material = sio2
# the oxide
box = 0   1um 0 0.6um 0.4um 0.6um

[region well]
material = silicon
prism = x 0.2um 0.8um, 0 0, \
        0.6um 0, 0.6um 0.2um, 0 0.2um

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[contact gate]
kind = electrode
patch = gate
electrode = metal
work_function = 4.4

[equilibrium]
)";

// The deck build test's layered MOS stack (Unit 26l slice 0), which meshes in seconds.
const std::string mos_deck = R"([mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
patch gate = 0.3um 0 0.6um  0.7um 0.6um 0.6um
patch body = 0 0 0  1um 0.6um 0
max_size = 0.15um
interface_resolution = 0.02um

[region si]
material = silicon

[region ox]
material = sio2
box = 0 1um 0 0.6um 0.4um 0.6um

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[contact gate]
kind = electrode
patch = gate
electrode = metal
work_function = 4.4

[equilibrium]
)";

// The same stack as an extrusion along z through two cells.
const std::string extruded_mos_deck = R"([mesh]
type = extrude3d
axis = z
stations = 0, 0.4um, 0.6um
cell left = 0 0, 0.5um 0, 0.5um 0.6um, 0 0.6um
cell right = 0.5um 0, 1um 0, 1um 0.6um, 0.5um 0.6um
patch gate = 0.5um 0 0.6um  1um 0.6um 0.6um
patch body = 0 0 0  1um 0.6um 0
max_size = 0.15um
interface_resolution = 0.02um

[region si]
material = silicon
extrude = left right 0:0.4um

[region ox]
material = sio2
extrude = left right 0.4um:0.6um

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[contact gate]
kind = electrode
patch = gate
electrode = metal
work_function = 4.4

[equilibrium]
)";

std::string read_file(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    std::stringstream s;
    s << in.rdbuf();
    return s.str();
}

std::string deck_file(const char* name) { return read_file(std::filesystem::path(NITCAD_DECKS_DIR) / name); }

// The text with CRLF line ends made LF (the decks may be checked out either way).
std::string lf(std::string t) {
    std::erase(t, '\r');
    return t;
}

Model must_read(const std::string& text) {
    auto m = read_model(text);
    if (!m) FAIL(m.error());
    return std::move(*m);
}

std::string must_write(const std::string& text, const Model& m) {
    auto t = write_model(text, m);
    if (!t) FAIL(t.error());
    return std::move(*t);
}

// The model without where its parts came from (a write and a read renumber them).
Model strip(Model m) {
    for (Region& r : m.regions) r.origin = -1;
    for (Contact& c : m.contacts) c.origin = -1;
    return m;
}

std::vector<std::string> lines_of(const std::string& t) {
    std::vector<std::string> out;
    std::istringstream in(t);
    for (std::string l; std::getline(in, l);) out.push_back(l);
    return out;
}

bool geometry_line(std::string_view l) {
    while (!l.empty() && (l.back() == '\r' || l.back() == ' ')) l.remove_suffix(1);
    if (l.empty()) return true;
    for (const char* k : {"box =", "prism =", "material =", "overlap =", "[region ", "[contact ", "patch", "layered =", "domain =", "kind =", "type =",
                          "axis =", "stations =", "cell ", "extrude =", "sweep =", "lines =", "straddle", "layer =", "min_angle", "        "}) {
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

std::string build_error(const std::string& text) {
    const auto d = NiTCAD::deck::parse(text, "t.deck");
    if (!d) return d.error().message;
    const auto p = NiTCAD::deck::build(*d);
    return p ? std::string() : p.error().message;
}

// The exact volume per region (quanta³), layered: the topmost solid at each cell of the solids' grid. Prisms of the
// model are rectilinear, so the grid of every coordinate resolves them; a cell's centre decides.
std::map<int, double> exact_volumes(const Model& m) {
    std::vector<Coord> g[3];
    for (int a = 0; a < 3; ++a) g[a] = {m.domain.lo[a], m.domain.hi[a]};
    for (const Region& r : m.regions) {
        for (const Solid& s : r.solids) {
            for (const P3& v : s.vertices()) {
                for (int a = 0; a < 3; ++a) g[a].push_back(std::clamp(v[a], m.domain.lo[a], m.domain.hi[a]));
            }
        }
    }
    for (auto& v : g) {
        std::ranges::sort(v);
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    const auto inside = [](const std::vector<Pt>& poly, double u, double v) {
        bool in = false;
        for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
            const double ax = static_cast<double>(poly[i].x), ay = static_cast<double>(poly[i].y), bx = static_cast<double>(poly[j].x),
                         by = static_cast<double>(poly[j].y);
            if ((ay > v) != (by > v) && u < (bx - ax) * (v - ay) / (by - ay) + ax) in = !in;
        }
        return in;
    };
    std::map<int, double> out;
    for (std::size_t i = 0; i + 1 < g[0].size(); ++i) {
        for (std::size_t j = 0; j + 1 < g[1].size(); ++j) {
            for (std::size_t k = 0; k + 1 < g[2].size(); ++k) {
                const double c[3] = {(g[0][i] + g[0][i + 1]) / 2.0, (g[1][j] + g[1][j + 1]) / 2.0, (g[2][k] + g[2][k + 1]) / 2.0};
                int region = m.background();
                bool held = false;  // Unit 26m MT4: a fill region takes only what no earlier solid holds
                for (std::size_t r = 0; r < m.regions.size(); ++r) {
                    bool in_region = false;
                    for (const Solid& s : m.regions[r].solids) {
                        bool in = false;
                        if (s.kind == SolidKind::box) {
                            in = c[0] > s.box.lo.x && c[0] < s.box.hi.x && c[1] > s.box.lo.y && c[1] < s.box.hi.y && c[2] > s.box.lo.z && c[2] < s.box.hi.z;
                        } else {
                            in = c[s.axis] > s.from && c[s.axis] < s.to && inside(s.polygon, c[(s.axis + 1) % 3], c[(s.axis + 2) % 3]);
                        }
                        in_region = in_region || in;
                    }
                    if (in_region && (!m.regions[r].fill || !held)) region = static_cast<int>(r);
                    held = held || in_region;
                }
                out[region] += static_cast<double>(g[0][i + 1] - g[0][i]) * static_cast<double>(g[1][j + 1] - g[1][j]) *
                               static_cast<double>(g[2][k + 1] - g[2][k]);
            }
        }
    }
    return out;
}

// The partition kernel's complex of the model passes check_complex, and its volume per region is the exact one.
void check_partition(const Model& m) {
    const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(m));
    if (!complex) FAIL(complex.error().message);
    const auto report = NiTCAD::mesh::check_complex(*complex);
    if (!report) FAIL(report.error().message);
    std::map<int, double> got;
    for (std::size_t v = 0; v < complex->volumes.size(); ++v) got[complex->volumes[v].region] += report->volume[v];
    const auto want = exact_volumes(m);
    REQUIRE(got.size() == want.size());
    for (const auto& [r, q3] : want) {
        const double cm3 = q3 * 1e-24;
        INFO("region " << r);
        CHECK(std::abs(got[r] - cm3) <= 1e-12 * cm3);
    }
}

}  // namespace

TEST_CASE("app solid: decks read; without edits the text comes back byte for byte (SE3)", "[app][solid]") {
    std::string crlf;
    for (const char c : well_deck) {
        if (c == '\n') crlf += '\r';
        crlf += c;
    }
    for (const std::string& text : {well_deck, mos_deck, extruded_mos_deck, crlf, deck_file("finfet_3d_reduced.deck"),
                                    deck_file("finfet_3d_shaped.deck"), deck_file("finfet_3d_swept.deck"), deck_file("finfet_3d_shaped_swept.deck")}) {
        const Model m = must_read(text);
        INFO(m.why);
        REQUIRE(m.mesh != MeshKind::unsupported);
        CHECK(must_write(text, m) == text);
        CHECK(must_write(must_write(text, m), m) == text);
    }
    const Model w = must_read(well_deck);
    CHECK(w.mesh == MeshKind::layered);
    CHECK(w.domain == Box3{{0, 0, 0}, {um, 6000, 6000}});
    REQUIRE(w.regions.size() == 3);
    CHECK(w.background() == 0);
    CHECK(w.regions[1].solids[0] == Solid{SolidKind::box, {{0, 0, 4000}, {um, 6000, 6000}}, 2, 0, 0, {}});
    const Solid& well = w.regions[2].solids[0];
    CHECK(well.kind == SolidKind::prism);
    CHECK(well.axis == 0);
    CHECK((well.from == 2000 && well.to == 8000));
    CHECK(well.polygon == std::vector<Pt>{{0, 0}, {6000, 0}, {6000, 2000}, {0, 2000}});
    CHECK(well.bounds() == Box3{{2000, 0, 0}, {8000, 6000, 2000}});
    REQUIRE(w.contacts.size() == 2);
    CHECK(w.contacts[0].box == Box3{{0, 0, 0}, {um, 6000, 0}});
    CHECK(w.contacts[1].box == Box3{{3000, 0, 6000}, {7000, 6000, 6000}});
    CHECK(w.contacts[1].kind == "electrode");

    // Not drawn here, each with why.
    for (const auto& [text, why] : std::vector<std::pair<std::string, std::string>>{
             {deck_file("finfet_3d.deck"), "3D tensor"},
             {read_file(std::filesystem::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck"), "2D editor"},
             {"[mesh]\ntype = generate3d\ncuboid = 0 1um 0 1um 0 1um\n[region a]\nmaterial = silicon\n", "without layers"}}) {
        const Model m = must_read(text);
        CHECK(m.mesh == MeshKind::unsupported);
        CHECK(m.why.find(why) != std::string::npos);
        CHECK(!write_model(text, m));
    }
}

TEST_CASE("app solid: the partition of the models read is exact (SE2's model side)", "[app][solid]") {
    check_partition(must_read(well_deck));
    check_partition(must_read(mos_deck));
}

TEST_CASE("app solid: an edit rewrites only geometry lines, and the deck still builds (SE3)", "[app][solid]") {
    Model m = must_read(mos_deck);
    REQUIRE(move_face(m, {1, 0}, 4, 500));  // the oxide's bottom up to 0.45um
    REQUIRE(move_contact(m, 1, {500, 0, 0}));
    REQUIRE(rename_region(m, 1, "oxide"));
    const std::string out = must_write(mos_deck, m);
    CHECK(only_geometry(mos_deck, out) == "");
    CHECK(out.find("[region oxide]\nmaterial = sio2\nbox = 0 1um 0 0.6um 0.45um 0.6um\n") != std::string::npos);
    CHECK(out.find("patch gate = 0.35um 0 0.6um  0.75um 0.6um 0.6um\n") != std::string::npos);
    CHECK(strip(must_read(out)) == strip(m));
    CHECK(build_error(out) == "");
    // A patch box no contact names (a named boundary) is the deck's, kept through edits.
    std::string spare = mos_deck;
    spare.replace(spare.find("max_size"), 0, "patch side = 0 0 0  0 0.6um 0.6um\n");
    Model s = must_read(spare);
    REQUIRE(delete_contact(s, 0));
    const std::string spare_out = must_write(spare, s);
    CHECK(spare_out.find("patch side = 0 0 0  0 0.6um 0.6um\n") != std::string::npos);
    CHECK(spare_out.find("patch body") == std::string::npos);

    // The well deck: comments and the continued prism kept when other solids change; the prism rewritten on one line.
    Model w = must_read(well_deck);
    REQUIRE(add_box(w, {{1000, 1000, 2000}, {3000, 3000, 4000}}, "silicon") == 3);
    CHECK(w.regions[3].name == "silicon_1");
    std::string t = must_write(well_deck, w);
    CHECK(t.find("# the oxide\n") != std::string::npos);
    CHECK(t.find("prism = x 0.2um 0.8um, 0 0, \\\n        0.6um 0, 0.6um 0.2um, 0 0.2um\n") != std::string::npos);
    CHECK(t.find("\n[region silicon_1]\nmaterial = silicon\nbox = 0.1um 0.3um 0.1um 0.3um 0.2um 0.4um\n\n[doping]") != std::string::npos);
    REQUIRE(move_face(w, {2, 0}, 1, 1000));  // the well's far cap to 0.9um
    t = must_write(well_deck, w);
    CHECK(t.find("prism = x 0.2um 0.9um, 0 0, 0.6um 0, 0.6um 0.2um, 0 0.2um\n\n[region silicon_1]") != std::string::npos);
    CHECK(only_geometry(well_deck, t) == "");
    CHECK(strip(must_read(t)) == strip(w));
    check_partition(w);

    // A new contact: its patch box in [mesh] after the last, its section at the end; deleting one drops its patch.
    REQUIRE(add_contact(w, "source", "ohmic", {{0, 1000, 1000}, {0, 5000, 3000}}) == 2);
    REQUIRE(delete_contact(w, 0));
    t = must_write(well_deck, w);
    CHECK(t.find("patch body") == std::string::npos);
    CHECK(t.find("patch gate = 0.3um 0 0.6um  0.7um 0.6um 0.6um   # the gate on top\npatch source = 0 0.1um 0.1um  0 0.5um 0.3um\n") !=
          std::string::npos);
    CHECK(t.ends_with("[equilibrium]\n\n[contact source]\nkind = ohmic\npatch = source\n"));
    CHECK(strip(must_read(t)) == strip(w));
}

TEST_CASE("app solid: an extrude3d deck reads as prisms, exact, and converts to a layered deck (SQ5, SE3)", "[app][solid]") {
    SECTION("the reduced FinFET") {
        const std::string text = deck_file("finfet_3d_reduced.deck");
        Model m = must_read(text);
        REQUIRE(m.mesh == MeshKind::extruded);
        CHECK(m.convert_why == "");
        CHECK(m.domain == Box3{{0, 0, 0}, {850, 520, 140}});
        std::vector<std::size_t> counts;
        for (const Region& r : m.regions) counts.push_back(r.solids.size());
        CHECK(counts == std::vector<std::size_t>{1, 4, 1, 1, 1});
        const Solid& fin = m.regions[4].solids[0];
        CHECK((fin.kind == SolidKind::prism && fin.axis == 0 && fin.from == 0 && fin.to == 850));
        CHECK(fin.polygon == std::vector<Pt>{{200, 20}, {500, 20}, {500, 120}, {200, 120}});
        CHECK(m.contacts[2].box == Box3{{300, 200, 0}, {550, 520, 140}});
        // Not edited before the conversion.
        CHECK(!move_solid(m, {4, 0}, {10, 0, 0}));
        CHECK(add_box(m, {{0, 0, 0}, {10, 10, 10}}, "silicon") < 0);
        CHECK(must_write(text, m) == text);

        const auto note = convert_to_layered(m);
        REQUIRE(note);
        CHECK(note->find("layered generate3d") != std::string::npos);
        const std::string out = must_write(text, m);
        CHECK(only_geometry(text, out) == "");
        CHECK(lf(out).find("type = generate3d\nlayered = on\ndomain = 0 0.085um 0 0.052um 0 0.014um\n") != std::string::npos);
        const std::string mesh_part = lf(out).substr(0, lf(out).find("[region "));
        for (const char* key : {"\nstations", "\ncell ", "\naxis ="}) CHECK(mesh_part.find(key) == std::string::npos);
        CHECK(out.find("extrude =") == std::string::npos);
        CHECK(lf(out).find("[region silicon]\nmaterial = silicon\nprism = x 0 0.085um, 0.02um 0.002um, 0.05um 0.002um, 0.05um 0.012um, 0.02um 0.012um\n") !=
              std::string::npos);
        CHECK(lf(out).find("patch gate = 30nm 20nm 0  55nm 52nm 14nm\n") != std::string::npos);  // kept as written
        CHECK(out.find("# Region 4: the fin, the whole length") != std::string::npos);  // comments kept
        const Model again = must_read(out);
        CHECK(again.mesh == MeshKind::layered);
        CHECK(strip(again) == strip(m));
        check_partition(again);  // the layered complex has the extrusion's volumes, exactly (SE2)
    }
    SECTION("a coarse stack converts, builds, and keeps its regions and contacts") {
        Model m = must_read(extruded_mos_deck);
        REQUIRE(convert_to_layered(m));
        const std::string out = must_write(extruded_mos_deck, m);
        CHECK(out.find("[region si]\nmaterial = silicon\nprism = z 0 0.4um, 0 0, 0.5um 0, 0.5um 0.6um, 0 0.6um\n"
                       "prism = z 0 0.4um, 0.5um 0, 1um 0, 1um 0.6um, 0.5um 0.6um\n") != std::string::npos);
        CHECK(build_error(extruded_mos_deck) == "");
        CHECK(build_error(out) == "");
        check_partition(must_read(out));
    }
    SECTION("a swept extrusion's keys are dropped with the sweep") {
        const std::string text = deck_file("finfet_3d_swept.deck");
        Model m = must_read(text);
        CHECK(m.swept);
        const auto note = convert_to_layered(m);
        REQUIRE(note);
        CHECK(note->find("swept") != std::string::npos);
        const std::string out = must_write(text, m);
        CHECK(only_geometry(text, out) == "");
        const std::string mesh_part = lf(out).substr(0, lf(out).find("[region "));
        for (const char* key : {"\nsweep =", "\nlines =", "\nstraddle", "\nmin_angle", "\naxis =", "\nlayer ="}) CHECK(mesh_part.find(key) == std::string::npos);
        CHECK(strip(must_read(out)) == strip(m));
    }
    SECTION("the shaped FinFET (slanted fin cells) converts and resolves through the slanted kernel (slice 4, SD6)") {
        const std::string text = deck_file("finfet_3d_shaped.deck");
        Model m = must_read(text);
        const auto note = convert_to_layered(m);
        REQUIRE(note);
        // The source and drain boxes cut across the spacer's slanted cells at the ends: said.
        CHECK(note->find("contact 'source' cuts across the end of region 'spacer'") != std::string::npos);
        CHECK(note->find("contact 'drain' cuts across the end of region 'spacer'") != std::string::npos);
        CHECK(note->find("of region 'silicon'") == std::string::npos);
        const Model again = must_read(must_write(text, m));
        CHECK(strip(again) == strip(m));
        const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(again));
        if (!complex) FAIL(complex.error().message);
        const auto report = NiTCAD::mesh::check_complex(*complex);
        if (!report) FAIL(report.error().message);
        // The extrusion's volume per region: its prisms' areas times lengths (cells tile, sections are disjoint).
        std::map<int, double> want, got;
        for (std::size_t r = 0; r < again.regions.size(); ++r) {
            for (const Solid& s : again.regions[r].solids) {
                std::int64_t a2 = 0;
                for (std::size_t k = 0; k < s.polygon.size(); ++k) {
                    const Pt p = s.polygon[k], q = s.polygon[(k + 1) % s.polygon.size()];
                    a2 += p.x * q.y - q.x * p.y;
                }
                want[static_cast<int>(r)] += 0.5 * static_cast<double>(a2 < 0 ? -a2 : a2) * static_cast<double>(s.to - s.from) * 1e-24;
            }
        }
        for (std::size_t v = 0; v < complex->volumes.size(); ++v) got[complex->volumes[v].region] += report->volume[v];
        REQUIRE(got.size() == want.size());
        for (const auto& [r, cm3] : want) {
            INFO("region " << r);
            CHECK(std::abs(got[r] - cm3) <= 1e-12 * cm3);
        }
    }
    SECTION("a void in the extrusion is refused with why") {
        std::string text = lf(deck_file("finfet_3d_reduced.deck"));
        text.replace(text.find("extrude = fin 0:85nm"), 20, "extrude = fin 0:55nm");
        Model m = must_read(text);
        CHECK(m.convert_why == "cell 'fin' is void in section 2 (a layered deck fills its domain)");
        CHECK(convert_to_layered(m).error() == m.convert_why);
        CHECK(m.mesh == MeshKind::extruded);
        text = lf(deck_file("finfet_3d_reduced.deck"));
        text.replace(text.find("extrude = box 0:85nm\n"), 21, "");
        CHECK(must_read(text).convert_why.starts_with("region 'buried_oxide' holds no cell"));
    }
}

TEST_CASE("app solid: the edits, refused where they would leave a solid invalid", "[app][solid]") {
    Model m = must_read(well_deck);
    const SolidRef ox{1, 0}, well{2, 0};
    CHECK(move_face(m, ox, 4, 500));
    CHECK(m.regions[1].solids[0].box.lo.z == 4500);
    CHECK(!move_face(m, ox, 4, 1500));  // empty
    CHECK(!move_face(m, ox, 6, 10));
    CHECK(move_face(m, well, 4, 500));  // the side over edge 2 (v = 0.2um) up
    CHECK(m.regions[2].solids[0].polygon == std::vector<Pt>{{0, 0}, {6000, 0}, {6000, 2500}, {0, 2500}});
    CHECK(move_face(m, well, 1, 1000));
    CHECK(m.regions[2].solids[0].to == 9000);
    CHECK(!set_extent(m, well, 5000, 5000));
    CHECK(set_vertex(m, well, 2, {5000, 2000}));  // its neighbours follow
    CHECK(m.regions[2].solids[0].polygon == std::vector<Pt>{{0, 0}, {5000, 0}, {5000, 2000}, {0, 2000}});
    // Slanted since slice 4 (SD6): a vertex moved alone, a vertex inserted exactly on the slanted edge (one off it
    // refused), the slanted side not moved as a face; then put back.
    REQUIRE(set_vertex(m, well, 2, {4000, 3000}, false));
    CHECK(m.regions[2].solids[0].polygon == std::vector<Pt>{{0, 0}, {5000, 0}, {4000, 3000}, {0, 2000}});
    CHECK(!move_face(m, well, 3, 100));  // the side over edge 1, from (5000, 0) to (4000, 3000)
    CHECK(!insert_vertex(m, well, 1, {4500, 1501}));
    REQUIRE(insert_vertex(m, well, 1, {4500, 1500}));
    CHECK(remove_vertex(m, well, 2));
    REQUIRE(set_vertex(m, well, 2, {5000, 2000}, false));
    CHECK(m.regions[2].solids[0].polygon == std::vector<Pt>{{0, 0}, {5000, 0}, {5000, 2000}, {0, 2000}});
    CHECK(insert_vertex(m, well, 0, {3000, 0}));
    CHECK(!insert_vertex(m, well, 0, {3000, 10}));
    CHECK(m.regions[2].solids[0].polygon.size() == 5);
    CHECK(remove_vertex(m, well, 1));
    CHECK(!remove_vertex(m, well, 9));
    CHECK(move_solid(m, well, {100, 200, 300}));
    CHECK(m.regions[2].solids[0].from == 2100);
    CHECK(m.regions[2].solids[0].polygon[0] == Pt{200, 300});
    CHECK(!set_box(m, well, {{0, 0, 0}, {1, 1, 1}}));

    const int l = add_prism(m, 2, 0, 1000, {{0, 0}, {2000, 0}, {2000, 1000}, {1000, 1000}, {1000, 2000}, {0, 2000}}, "sio2");
    REQUIRE(l == 3);
    CHECK(m.regions[3].name == "sio2_1");
    CHECK(mirror_solid(m, {l, 0}, 0));
    auto poly = m.regions[3].solids[0].polygon;
    std::ranges::sort(poly);
    CHECK(poly == std::vector<Pt>{{0, 0}, {0, 1000}, {1000, 1000}, {1000, 2000}, {2000, 0}, {2000, 2000}});
    CHECK(m.regions[3].solids[0].bounds() == Box3{{0, 0, 0}, {2000, 2000, 1000}});
    CHECK(add_prism(m, 2, 0, 1000, {{0, 0}, {2000, 2000}, {2000, 0}, {0, 2000}}, "sio2") < 0);  // crossing itself
    CHECK(add_prism(m, 2, 1000, 0, {{0, 0}, {2000, 0}, {2000, 1000}, {0, 1000}}, "sio2") < 0);
    CHECK(add_box(m, {{0, 0, 0}, {0, 10, 10}}, "sio2") < 0);
    const int d = duplicate_solid(m, ox, {0, 0, 100});
    CHECK(d == 4);
    CHECK(m.regions[4].solids[0].box.lo.z == 4600);
    CHECK(m.regions[4].material == "sio2");
    CHECK(!rename_region(m, 4, "ox"));
    CHECK(!rename_region(m, 4, "9a"));
    CHECK(set_material(m, 4, "hfo2"));
    CHECK(raise_region(m, 0));
    CHECK(m.regions[1].name == "si");
    CHECK(lower_region(m, 1));
    CHECK(!lower_region(m, 0));
    CHECK(!raise_region(m, 4));
    CHECK(delete_solid(m, {4, 0}));
    CHECK(m.regions.size() == 4);

    // Contacts: on the boundary with area; the domain's faces carry their contacts along.
    CHECK(add_contact(m, "inner", "ohmic", {{1000, 1000, 1000}, {2000, 2000, 2000}}) < 0);
    CHECK(add_contact(m, "edge", "ohmic", {{0, 0, 0}, {0, 0, 1000}}) < 0);  // no area on a face
    CHECK(add_contact(m, "gate", "ohmic", {{0, 0, 0}, {0, 1000, 1000}}) < 0);  // the name taken
    CHECK(add_contact(m, "side", "ohmic", {{0, 0, 0}, {0, 1000, 1000}}) == 2);
    CHECK(!move_contact(m, 2, {10, 0, 0}));  // off the face
    CHECK(set_contact_kind(m, 2, "schottky"));
    CHECK(set_domain(m, {{0, 0, 0}, {um, 6000, 7000}}));
    CHECK(m.contacts[1].box == Box3{{3000, 0, 7000}, {7000, 6000, 7000}});
    CHECK(m.contacts[0].box == Box3{{0, 0, 0}, {um, 6000, 0}});
    CHECK(!set_domain(m, {{0, 0, 0}, {um, 6000, 0}}));
    CHECK(set_domain(m, {{0, 0, -1000}, {um, 6000, 7000}}));  // the low face carries its contacts too
    CHECK(m.contacts[0].box == Box3{{0, 0, -1000}, {um, 6000, -1000}});
    CHECK(m.contacts[2].box == Box3{{0, 0, -1000}, {0, 1000, 1000}});
}

TEST_CASE("app solid: 200 random edits, undone one by one, give the deck back byte for byte (SE4, SE5)", "[app][solid]") {
    for (const bool extruded : {false, true}) {
        const std::string start = extruded ? deck_file("finfet_3d_reduced.deck") : well_deck;
        std::mt19937 rng(extruded ? 26u : 12u);
        const auto pick = [&](int n) { return std::uniform_int_distribution<int>(0, n - 1)(rng); };
        std::string text = start;
        std::vector<std::string> texts{text};
        EditHistory h;
        int done = 0, tries = 0;
        double write_ms = 0.0;
        if (extruded) {  // the first step converts it
            Model m = must_read(text);
            REQUIRE(convert_to_layered(m));
            h.record(text);
            text = must_write(text, m);
            texts.push_back(text);
        }
        while (done < 200 && tries < 20000) {
            ++tries;
            Model m = must_read(text);
            const Box3 dom = m.domain;
            const auto at = [&](int a) { return dom.lo[a] + (dom.hi[a] - dom.lo[a]) * pick(11) / 10; };
            const auto len = [&](int a) { return static_cast<Coord>(pick(9) - 4) * std::max<Coord>(1, (dom.hi[a] - dom.lo[a]) / 20); };
            std::vector<SolidRef> solids;
            for (int r = 0; r < static_cast<int>(m.regions.size()); ++r) {
                for (int s = 0; s < static_cast<int>(m.regions[static_cast<std::size_t>(r)].solids.size()); ++s) solids.push_back({r, s});
            }
            const auto any = [&]() -> SolidRef { return solids.empty() ? SolidRef{} : solids[static_cast<std::size_t>(pick(static_cast<int>(solids.size())))]; };
            const int nr = static_cast<int>(m.regions.size()), nc = static_cast<int>(m.contacts.size());
            const auto face_box = [&] {  // a patch box on a face of the domain
                const int a = pick(3), u = (a + 1) % 3, v = (a + 2) % 3;
                Box3 b;
                b.lo[a] = b.hi[a] = pick(2) ? dom.lo[a] : dom.hi[a];
                b.lo[u] = at(u);
                b.hi[u] = std::min(dom.hi[u], b.lo[u] + (dom.hi[u] - dom.lo[u]) / 5 + 1);
                b.lo[v] = at(v);
                b.hi[v] = std::min(dom.hi[v], b.lo[v] + (dom.hi[v] - dom.lo[v]) / 5 + 1);
                return b;
            };
            bool ok = false, removes = false;
            switch (pick(16)) {
            case 0: ok = move_solid(m, any(), {len(0), len(1), len(2)}); break;
            case 1: {
                const SolidRef s = any();
                ok = s.region >= 0 && move_face(m, s, pick(8), len(pick(3)));
                break;
            }
            case 2: {
                const Box3 b = normalized({at(0), at(1), at(2)}, {at(0), at(1), at(2)});
                ok = add_box(m, b, pick(2) ? "silicon" : "sio2") >= 0;
                break;
            }
            case 3: {
                const int axis = pick(3), u = (axis + 1) % 3, v = (axis + 2) % 3;
                const Coord u0 = at(u), v0 = at(v), du = std::max<Coord>(2, (dom.hi[u] - dom.lo[u]) / 4), dv = std::max<Coord>(2, (dom.hi[v] - dom.lo[v]) / 4);
                const Coord f = at(axis);
                ok = add_prism(m, axis, f, f + std::max<Coord>(1, (dom.hi[axis] - dom.lo[axis]) / 3),
                               {{u0, v0}, {u0 + du, v0}, {u0 + du, v0 + dv / 2}, {u0 + du / 2, v0 + dv / 2}, {u0 + du / 2, v0 + dv}, {u0, v0 + dv}},
                               "silicon") >= 0;
                break;
            }
            case 4: ok = removes = solids.size() > 1 && delete_solid(m, any()); break;
            case 5: ok = raise_region(m, pick(std::max(nr, 1))); break;
            case 6: ok = lower_region(m, pick(std::max(nr, 1))); break;
            case 7: ok = rename_region(m, pick(std::max(nr, 1)), "r" + std::to_string(pick(1000))); break;
            case 8: ok = set_material(m, pick(std::max(nr, 1)), pick(2) ? "silicon" : "hfo2"); break;
            case 9: ok = add_contact(m, "c" + std::to_string(pick(1000)), "ohmic", face_box()) >= 0; break;
            case 10: ok = nc > 0 && set_contact(m, pick(nc), face_box()); break;
            case 11: ok = removes = nc > 1 && delete_contact(m, pick(nc)); break;
            case 12: ok = nc > 0 && rename_contact(m, pick(nc), "k" + std::to_string(pick(1000))); break;
            case 13: ok = duplicate_solid(m, any(), {len(0), len(1), len(2)}) >= 0; break;
            case 14: ok = mirror_solid(m, any(), pick(3)); break;
            default: {
                const SolidRef s = any();
                if (s.region < 0) break;
                const Solid& sol = m.regions[static_cast<std::size_t>(s.region)].solids[static_cast<std::size_t>(s.solid)];
                if (sol.kind != SolidKind::prism) break;
                const int k = pick(static_cast<int>(sol.polygon.size()));
                const Pt p = sol.polygon[static_cast<std::size_t>(k)];
                ok = set_vertex(m, s, k, {p.x + len((sol.axis + 1) % 3), p.y + len((sol.axis + 2) % 3)});
                break;
            }
            }
            if (!ok) continue;
            const auto t0 = std::chrono::steady_clock::now();
            std::string next = must_write(text, m);
            write_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            {
                INFO("edit " << done << "\n" << next);
                REQUIRE(strip(must_read(next)) == strip(m));
                REQUIRE(only_geometry(text, next, removes) == "");
            }
            h.record(text);
            text = std::move(next);
            texts.push_back(text);
            ++done;
        }
        REQUIRE(done == 200);
        WARN("SE5: " << (extruded ? "converted FinFET" : "well") << " deck, mean write " << write_ms / done << " ms over " << done << " edits");
        // Release only, as KG6 (owner, 2026-11-19): CI's Debug runner wrote 3.57 ms (fd6a64a); there it is recorded.
#ifdef NDEBUG
        CHECK(write_ms / done <= 2.0);  // SE5
#endif
        for (std::size_t k = texts.size() - 1; k > 0; --k) {
            auto t = h.undo(text);
            REQUIRE(t);
            text = std::move(*t);
            REQUIRE(text == texts[k - 1]);
        }
        CHECK(text == start);
        CHECK(!h.can_undo());
        for (std::size_t k = 1; k < texts.size(); ++k) {
            auto t = h.redo(text);
            REQUIRE(t);
            text = std::move(*t);
            REQUIRE(text == texts[k]);
        }
        CHECK(!h.can_redo());
        if (!extruded) check_partition(must_read(text));  // the FinFET has no background: a move may uncover (slice 3 says so)
    }
}

TEST_CASE("app solid: a deck of 100 solids is written within 2 ms per edit (SE5)", "[app][solid]") {
    std::string deck = "[mesh]\ntype = generate3d\nlayered = on\ndomain = 0 100um 0 100um 0 10um\npatch a = 0 0 0  0 1um 1um\n\n[region bg]\nmaterial = silicon\n";
    for (int k = 0; k < 100; ++k) {
        deck += std::format("\n# block {}\n[region r{}]\nmaterial = sio2\n{}\n", k, k,
                            k % 2 ? std::format("box = {}um {}um {}um {}um 1um 3um", k % 10 * 10, k % 10 * 10 + 4, k / 10 * 10, k / 10 * 10 + 4)
                                  : std::format("prism = z 2um 6um, {}um {}um, {}um {}um, {}um {}um, {}um {}um", k % 10 * 10, k / 10 * 10, k % 10 * 10 + 4,
                                                k / 10 * 10, k % 10 * 10 + 4, k / 10 * 10 + 4, k % 10 * 10, k / 10 * 10 + 4));
    }
    deck += "\n[contact a]\nkind = ohmic\npatch = a\n";
    Model m = must_read(deck);
    REQUIRE(m.regions.size() == 101);
    const auto t0 = std::chrono::steady_clock::now();
    std::string text = deck;
    for (int k = 0; k < 50; ++k) {
        REQUIRE(move_solid(m, {1 + k * 2, 0}, {1000, 0, 0}));
        text = must_write(text, m);
        m = must_read(text);
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() / 50;
    WARN("SE5: 100 solids, write and read " << ms << " ms per edit");
    // Release only (owner, 2026-11-19, as KG6): the unoptimised Debug build records the time.
#ifdef NDEBUG
    CHECK(ms <= 2.0);
#endif
}

TEST_CASE("app solid: snapping to vertices, midpoints, edges, faces and the work plane's grid", "[app][solid]") {
    const Model m = must_read(well_deck);
    const WorkPlane floor{2, 0};
    CHECK(snap(m, floor, {2010, 10, 20}, 50) == Snap{{2000, 0, 0}, SnapKind::vertex});  // the well's corner
    CHECK(snap(m, floor, {2010, 10, 20}, 50, SolidRef{2, 0}) == Snap{{2010, 0, 0}, SnapKind::edge});  // the domain's edge
    CHECK(snap(m, floor, {5010, 20, 6000}, 50) == Snap{{5000, 0, 6000}, SnapKind::midpoint});
    CHECK(snap(m, floor, {1000, 15, 6010}, 50) == Snap{{1000, 0, 6000}, SnapKind::edge});
    CHECK(snap(m, floor, {1000, 3000, 6030}, 50) == Snap{{1000, 3000, 6000}, SnapKind::face});
    CHECK(snap(m, floor, {5000, 1000, 2030}, 50) == Snap{{5000, 1000, 2000}, SnapKind::face});  // the well's top side
    // On the work plane's own plane its faces do not hold the point: the grid does.
    CHECK(snap(m, {2, 6000}, {1100, 3100, 6000}, 50, {}, 1000) == Snap{{1000, 3000, 6000}, SnapKind::grid});
    CHECK(snap(m, {2, 6000}, {1100, 3100, 6000}, 50) == Snap{{1100, 3100, 6000}, SnapKind::none});
    CHECK(snap(m, {0, 3000}, {3000, -1600, 2400}, 50, {}, 1000) == Snap{{3000, -2000, 2000}, SnapKind::grid});
    CHECK(snap(m, floor, {1500, 3000, 3000}, 50) == Snap{{1500, 3000, 3000}, SnapKind::none});
    CHECK(grid_step(100.0) == NiTCAD::app::geo::grid_step(NiTCAD::app::geo::View{100.0, 0, 0, 96}));
    CHECK(grid_step(100.0, 192) >= grid_step(100.0));
}

TEST_CASE("app solid: the model names no Win32 type (SE10)", "[app][solid]") {
    for (const char* f : {"include/NiTCAD/app/solid.hpp", "src/app/solid.cpp", "src/app/deck_text.hpp", "include/NiTCAD/app/solid_view.hpp", "src/app/solid_view.cpp"}) {
        const std::string t = read_file(std::filesystem::path(NITCAD_SOURCE_DIR) / f);
        REQUIRE(!t.empty());
        for (const char* bad : {"windows.h", "Windows.h", "win32.hpp", "HWND", "d3d12"}) {
            INFO(f << ": " << bad);
            CHECK(t.find(bad) == std::string::npos);
        }
    }
}

TEST_CASE("app solid: the build's errors placed on the model; the region at a point (SD5)", "[app][solid][check]") {
    const Model m = must_read(well_deck);
    const std::string t = lf(well_deck);
    const auto line_of = [&](const std::string& what) {
        const auto at = t.find(what);
        REQUIRE(at != std::string::npos);
        return 1 + static_cast<int>(std::count(t.begin(), t.begin() + static_cast<std::ptrdiff_t>(at), '\n'));
    };
    using K = Marker::Kind;
    // A region's solid line: that solid; another line of the region: all its solids; the background: the domain.
    Problem p = locate_problem(m, t, line_of("box = 0   1um"), "a box");
    REQUIRE(p.markers.size() == 1);
    CHECK(p.markers[0] == Marker{K::solid, {1, 0}, -1, {}});
    p = locate_problem(m, t, line_of("material = sio2"), "a material");
    REQUIRE(p.markers.size() == 1);
    CHECK(p.markers[0] == Marker{K::solid, {1, -1}, -1, {}});
    p = locate_problem(m, t, line_of("[region si]") + 1, "the background");
    REQUIRE(p.markers.size() == 1);
    CHECK(p.markers[0].kind == K::domain);
    // A contact's section; a [mesh] patch line: the contacts using it; another [mesh] line: the domain.
    p = locate_problem(m, t, line_of("[contact gate]") + 1, "a contact");
    REQUIRE(p.markers.size() == 1);
    CHECK(p.markers[0] == Marker{K::contact, {}, 1, {}});
    p = locate_problem(m, t, line_of("patch body"), "a patch");
    REQUIRE(p.markers.size() == 1);
    CHECK(p.markers[0] == Marker{K::contact, {}, 0, {}});
    p = locate_problem(m, t, line_of("max_size"), "a size");
    REQUIRE(p.markers.size() == 1);
    CHECK(p.markers[0].kind == K::domain);
    // A point the message names (cm): marked, with the solids whose bounds hold it; the message unchanged.
    const std::string named = "an interface piece has holes (unsupported) near (5e-05, 1e-05, 4e-05)";
    p = locate_problem(m, t, 2, named);
    CHECK(p.message == named);
    CHECK(std::ranges::count(p.markers, Marker{K::point, {}, -1, {5000, 1000, 4000}}) == 1);
    CHECK(std::ranges::count(p.markers, Marker{K::solid, {1, 0}, -1, {}}) == 1);  // the oxide's lower face
    CHECK(std::ranges::count(p.markers, K::solid, &Marker::kind) == 1);  // the well (z up to 0.2 um) does not hold it
    // Outside the geometry: none.
    CHECK(locate_problem(m, t, line_of("profile = acceptor"), "a profile").markers.empty());
    CHECK(locate_problem(m, t, 0, "no line").markers.empty());

    // The region at a point: the last solid holding it strictly, else the background.
    CHECK(region_at(m, {5000, 3000, 5000}) == 1);  // in the oxide
    CHECK(region_at(m, {5000, 1000, 1000}) == 2);  // in the well (over the substrate)
    CHECK(region_at(m, {5000, 1000, 4000}) == 0);  // on the oxide's lower face: not strictly inside it
    CHECK(region_at(m, {9000, 5000, 1000}) == 0);  // the substrate
    // Overlapping solids: the later one.
    Model over = m;
    REQUIRE(add_box(over, {{4000, 0, 3000}, {6000, 6000, 5000}}, "silicon") == 3);
    CHECK(region_at(over, {5000, 3000, 4500}) == 3);  // in the oxide too
    CHECK(region_at(over, {5000, 3000, 3500}) == 3);
    CHECK(region_at(over, {7000, 3000, 4500}) == 1);
}

// --- Unit 26m slice 1 (MT4, MT5; MC5's model side, MC6) ----------------------------------------------------------

namespace {

// A silicon block under an oxide cap, one sio2 fin buried in the silicon (features >= 6 s, gaps >= 4 s at s = 0.02um).
const std::string fin_deck = R"([mesh]
type = generate3d
layered = on
domain = 0 1.2um 0 0.6um 0 0.6um
patch gate = 0.3um 0 0.6um  0.9um 0.6um 0.6um
patch body = 0 0 0  1.2um 0.6um 0
max_size = 0.15um
interface_resolution = 0.02um

[region si]
material = silicon

[region ox]
material = sio2
box = 0 1.2um 0 0.6um 0.45um 0.6um

[region fin]
material = sio2
box = 0.15um 0.3um 0 0.6um 0.1um 0.3um

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[contact gate]
kind = electrode
patch = gate
electrode = metal
work_function = 4.4

[equilibrium]
)";

}  // namespace

TEST_CASE("app solid: overlap = fill read, written by a minimal edit, resolved and built (Unit 26m MT4)", "[app][solid][overlap]") {
    // On a new region: its line after the material; the solid takes only the background (the oxide keeps its part).
    Model m = must_read(mos_deck);
    const int r = add_box(m, {{2000, 0, 2500}, {8000, 6000, 6000}}, "silicon");
    REQUIRE(r == 2);
    CHECK_FALSE(set_overlap(m, r, false));  // no change
    REQUIRE(set_overlap(m, r, true));
    std::string out = must_write(mos_deck, m);
    CHECK(out.find("\n[region silicon_1]\nmaterial = silicon\noverlap = fill\nbox = 0.2um 0.8um 0 0.6um 0.25um 0.6um\n") != std::string::npos);
    CHECK(only_geometry(mos_deck, out) == "");
    const Model back = must_read(out);
    CHECK(back.regions[2].fill);
    CHECK(strip(back) == strip(m));
    check_partition(m);
    const auto v = exact_volumes(m);
    CHECK(v.at(2) == 6000.0 * 6000.0 * 1500.0);  // x 0.2-0.8, all y, z 0.25-0.4: under the oxide only
    CHECK(region_at(m, {5000, 3000, 3500}) == 2);
    CHECK(region_at(m, {5000, 3000, 5000}) == 1);  // the oxide keeps what it held
    CHECK(paint_order(m) == std::vector<int>{2, 1});
    CHECK(build_error(out) == "");

    // On a kept region: the line inserted after its material, then removed again, giving the deck back byte for byte.
    Model k = must_read(mos_deck);
    CHECK_FALSE(set_overlap(k, 0, true));  // the background has no solid
    REQUIRE(set_overlap(k, 1, true));
    out = must_write(mos_deck, k);
    CHECK(out.find("[region ox]\nmaterial = sio2\noverlap = fill\nbox = 0 1um 0 0.6um 0.4um 0.6um\n") != std::string::npos);
    REQUIRE(set_overlap(k, 1, false));
    CHECK(must_write(mos_deck, k) == mos_deck);
    // A deck's `overlap = replace` kept as it is; switched to fill, that line is the one rewritten.
    std::string said = mos_deck;
    said.replace(said.find("material = sio2\n"), 16, "material = sio2\noverlap = replace\n");
    Model s = must_read(said);
    CHECK_FALSE(s.regions[1].fill);
    CHECK(must_write(said, s) == said);
    REQUIRE(set_overlap(s, 1, true));
    out = must_write(said, s);
    CHECK(out.find("material = sio2\noverlap = fill\nbox") != std::string::npos);
    CHECK(out.find("overlap = replace") == std::string::npos);
}

TEST_CASE("app solid: Array, Mirror, Align and the quarter turn, each one edit that builds (Unit 26m MT5, MC6)", "[app][solid][transform]") {
    const Model m0 = must_read(fin_deck);
    REQUIRE(build_error(fin_deck) == "");
    const SolidRef fin{2, 0};
    // Array of 3 fins, 0.35um apart: two new regions on top, as sio2 fins.
    {
        Model m = m0;
        CHECK(array_solid(m, fin, {0, 0, 0}, 3) == -1);
        CHECK(array_solid(m, fin, {3500, 0, 0}, 1) == -1);
        REQUIRE(array_solid(m, fin, {3500, 0, 0}, 3) == 3);
        REQUIRE(m.regions.size() == 5);
        CHECK(m.regions[4].solids[0].box == Box3{{8500, 0, 1000}, {10000, 6000, 3000}});
        CHECK(m.regions[4].material == "sio2");
        const std::string out = must_write(fin_deck, m);
        CHECK(only_geometry(fin_deck, out) == "");
        CHECK(strip(must_read(out)) == strip(m));
        check_partition(m);
        CHECK(build_error(out) == "");
        // A fill region's copies fill too.
        Model f = m0;
        REQUIRE(set_overlap(f, 2, true));
        REQUIRE(array_solid(f, fin, {3500, 0, 0}, 2) == 3);
        CHECK(f.regions[3].fill);
    }
    // Mirror across the domain's centre plane x = 0.6um, and about its own centre (no change for a box).
    {
        Model m = m0;
        REQUIRE(mirror_solid_at(m, fin, 0, 12000));
        CHECK(m.regions[2].solids[0].box == Box3{{9000, 0, 1000}, {10500, 6000, 3000}});
        const std::string out = must_write(fin_deck, m);
        check_partition(m);
        CHECK(build_error(out) == "");
        Model c = m0;
        REQUIRE(mirror_solid_at(c, fin, 2, 1000 + 3000));
        CHECK(c == m0);
    }
    // Align: a second fin, lower, its top brought to the first's.
    {
        Model m = m0;
        const int r = add_box(m, {{6000, 0, 1000}, {7500, 6000, 2500}}, "sio2");
        REQUIRE(r == 3);
        CHECK_FALSE(align_solids(m, {{r, 0}}, 5, 2500));  // there already
        REQUIRE(align_solids(m, {{r, 0}, {r, 0}}, 5, 3000));
        CHECK(m.regions[3].solids[0].box == Box3{{6000, 0, 1500}, {7500, 6000, 3000}});
        check_partition(m);
        CHECK(build_error(must_write(fin_deck, m)) == "");
        CHECK_FALSE(align_solids(m, {{9, 0}}, 5, 0));
        CHECK_FALSE(align_solids(m, {{r, 0}}, 6, 0));
    }
    // The quarter turn: a box about its centre; four turns give it back; an L prism along x turned about z runs along y.
    {
        Model m = m0;
        REQUIRE(rotate_solid(m, fin, 2));  // 0.15 x 0.6 in (x, y) about (0.225, 0.3): 0.6 x 0.15
        CHECK(m.regions[2].solids[0].box == Box3{{-750, 2250, 1000}, {5250, 3750, 3000}});
        for (int k = 0; k < 3; ++k) REQUIRE(rotate_solid(m, fin, 2));
        CHECK(m == m0);
        Model w = must_read(well_deck);
        const SolidRef well{2, 0};
        const double before = exact_volumes(w).at(2);
        REQUIRE(rotate_solid(w, well, 2));
        const Solid& t = w.regions[2].solids[0];
        CHECK(t.kind == SolidKind::prism);
        CHECK(t.axis == 1);
        CHECK(t.from == 0);
        CHECK(t.to == 6000);
        CHECK(prism_ok(t.polygon));
        check_partition(w);
        CHECK(exact_volumes(w).at(2) == before);
        for (int k = 0; k < 3; ++k) REQUIRE(rotate_solid(w, well, 2));
        CHECK(w.regions[2].solids[0].axis == 0);
        CHECK(w.regions[2].solids[0].bounds() == must_read(well_deck).regions[2].solids[0].bounds());
    }
}

// --- Unit 26m slice 2 (MT6, MT7, MT12; MC6, MC13) -------------------------------------------------------------------

namespace {

// A silicon block under an oxide cap with a wide sio2 fin across it along y; s = 0.01um (6 s = 600 quanta, 4 s = 400).
const std::string tool_deck = R"([mesh]
type = generate3d
layered = on
domain = 0 2um 0 1um 0 1um
patch body = 0 0 0  2um 1um 0
max_size = 0.25um
interface_resolution = 0.01um

[region si]
material = silicon

[region ox]
material = sio2
box = 0 2um 0 1um 0.7um 1um

[region fin]
material = sio2
box = 0.4um 1um 0 1um 0.1um 0.5um

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";

// The box edge between corners i and i | bit.
int box_edge(const Solid& s, int i, int bit) {
    const auto e = s.edges();
    for (std::size_t k = 0; k < e.size(); ++k) {
        if (e[k] == std::pair{i, i | bit}) return static_cast<int>(k);
    }
    return -1;
}

double area_q2(const std::vector<Pt>& p) {
    double a = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const Pt q = p[i], r = p[(i + 1) % p.size()];
        a += static_cast<double>(q.x) * static_cast<double>(r.y) - static_cast<double>(r.x) * static_cast<double>(q.y);
    }
    return std::abs(a) / 2;
}

// The partition's volume of a region, cm³.
double region_volume(const Model& m, int region) {
    const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(m));
    if (!complex) FAIL(complex.error().message);
    const auto report = NiTCAD::mesh::check_complex(*complex);
    if (!report) FAIL(report.error().message);
    double v = 0;
    for (std::size_t k = 0; k < complex->volumes.size(); ++k) {
        if (complex->volumes[k].region == region) v += report->volume[k];
    }
    return v;
}

double prism_cm3(const Solid& s) { return area_q2(s.polygon) * static_cast<double>(s.to - s.from) * 1e-24; }

// An edit's deck: only geometry lines changed, read back to the same model, its undo the deck byte for byte.
std::string round_trip(const std::string& deck, const Model& m) {
    const std::string out = must_write(deck, m);
    CHECK(only_geometry(deck, out) == "");
    CHECK(strip(must_read(out)) == strip(m));
    CHECK(must_write(out, must_read(deck)) == deck);
    return out;
}

}  // namespace

TEST_CASE("app solid: arcs and corners on the quantum (Unit 26m MT7, MQ4)", "[app][solid][edges]") {
    using namespace NiTCAD::app::geo;
    const double right = std::numbers::pi / 2;
    CHECK(arc_segments(right, 1e6, 0) == 6);  // ⌈90° / 15°⌉
    CHECK(arc_segments(std::numbers::pi / 3, 1e6, 0) == 4);
    CHECK(arc_segments(0.1, 1e6, 0) == 2);
    CHECK(arc_segments(right, 2000, 600) == 5);  // 2 r sin(9°) = 626 >= 600, 2 r sin(7.5°) = 522 < 600
    CHECK(arc_segments(right, 300, 600) == 0);   // even one chord (424) is short
    CHECK(polygon_sides(2000, 600) == 16);
    CHECK(polygon_sides(1000, 600) == 10);  // 2 r sin(18°) = 618
    CHECK(polygon_sides(100, 600) == 0);
    const std::vector<Pt> sq{{0, 0}, {4000, 0}, {4000, 4000}, {0, 4000}};
    CHECK(std::abs(corner_turn(sq, 2) - right) < 1e-12);
    CHECK(fillet_tangent(sq, 2, 1000) == Catch::Approx(1000));
    CHECK(*chamfer_corner(sq, 2, 1000) == std::vector<Pt>{{0, 0}, {4000, 0}, {4000, 3000}, {3000, 4000}, {0, 4000}});
    CHECK_FALSE(chamfer_corner(sq, 2, 4001));
    CHECK(*chamfer_corner(sq, 2, 4000) == std::vector<Pt>{{0, 0}, {4000, 0}, {0, 4000}});  // merged with both neighbours
    const auto f = fillet_corner(sq, 2, 1000, 6);
    REQUIRE(f);
    REQUIRE(f->size() == 10);
    CHECK((*f)[2] == Pt{4000, 3000});
    CHECK((*f)[8] == Pt{3000, 4000});
    for (std::size_t k = 2; k <= 8; ++k) {  // on the circle about (3000, 3000), to the rounding
        CHECK(std::abs(std::hypot((*f)[k].x - 3000.0, (*f)[k].y - 3000.0) - 1000) <= 0.75);
    }
    CHECK_FALSE(fillet_corner(sq, 2, 4001, 6));
    CHECK_FALSE(fillet_corner({{0, 0}, {2000, 0}, {4000, 0}, {0, 4000}}, 1, 100, 2));  // a straight corner
    // A reflex corner of an L takes a fillet too (the arc bulging into the notch).
    const std::vector<Pt> ell{{0, 0}, {4000, 0}, {4000, 2000}, {2000, 2000}, {2000, 4000}, {0, 4000}};
    const auto r = fillet_corner(ell, 3, 500, 3);
    REQUIRE(r);
    CHECK(polygon_ok(*r));
    CHECK(area_q2(*r) > area_q2(ell));
    const auto c = regular_polygon({0, 0}, 1000, 4);
    CHECK(c == std::vector<Pt>{{1000, 0}, {0, 1000}, {-1000, 0}, {0, -1000}});
    // 2D shapes: a box's corner chamfered becomes a polygon; a tensor deck refuses.
    Geometry g;
    g.mesh = NiTCAD::app::geo::MeshKind::layered;
    g.domain = {{0, 0}, {10000, 10000}};
    g.regions.push_back({"bg", "silicon", {}, -1});
    REQUIRE(add_box(g, {{1000, 1000}, {5000, 5000}}, "sio2") == 1);
    REQUIRE(chamfer_shape(g, {1, 0}, 2, 1000));
    CHECK(g.regions[1].shapes[0].kind == ShapeKind::polygon);
    CHECK(g.regions[1].shapes[0].points.size() == 5);
    REQUIRE(fillet_shape(g, {1, 0}, 0, 1000, 6));
    CHECK(g.regions[1].shapes[0].points.size() == 11);
    g.mesh = NiTCAD::app::geo::MeshKind::tensor;
    CHECK_FALSE(chamfer_shape(g, {1, 0}, 1, 500));
}

TEST_CASE("app solid: Chamfer, Fillet, Round vertex, Taper and Cylinder, each one edit that builds (Unit 26m MT7, MC6)",
          "[app][solid][edges]") {
    const Model m0 = must_read(tool_deck);
    REQUIRE(build_error(tool_deck) == "");
    const SolidRef fin{2, 0};
    const Solid& box = m0.regions[2].solids[0];
    const int top_right = box_edge(box, 1 | 4, 2);  // x hi, z hi, along y
    REQUIRE(top_right >= 0);
    // The edge as a prism's corner: along y, (u, v) = (z, x).
    const auto ec = edge_corner(box, top_right);
    REQUIRE(ec);
    CHECK(ec->prism.axis == 1);
    CHECK(ec->prism.polygon == std::vector<Pt>{{1000, 4000}, {5000, 4000}, {5000, 10000}, {1000, 10000}});
    CHECK(ec->corner == 2);
    {
        Model m = m0;
        CHECK_FALSE(chamfer_edge(m, fin, top_right, 4500));  // longer than the 4000 wall
        REQUIRE(chamfer_edge(m, fin, top_right, 1000));
        const Solid& t = m.regions[2].solids[0];
        CHECK(t.polygon == std::vector<Pt>{{1000, 4000}, {5000, 4000}, {5000, 9000}, {4000, 10000}, {1000, 10000}});
        CHECK(region_volume(m, 2) == Catch::Approx((6000.0 * 4000 - 500'000) * 10000 * 1e-24).epsilon(1e-12));
        const std::string out = round_trip(tool_deck, m);
        CHECK(build_error(out) == "");
        // A cap's edge is not along the axis.
        CHECK_FALSE(chamfer_edge(m, fin, 0, 500));
    }
    {
        Model m = m0;
        const int n = NiTCAD::app::geo::arc_segments(std::numbers::pi / 2, 2000, 6 * 100);
        REQUIRE(n == 5);
        REQUIRE(fillet_edge(m, fin, top_right, 2000, n));
        const Solid& t = m.regions[2].solids[0];
        CHECK(t.polygon.size() == 9);
        CHECK(region_volume(m, 2) == Catch::Approx(prism_cm3(t)).epsilon(1e-12));
        const std::string out = round_trip(tool_deck, m);
        CHECK(build_error(out) == "");
        // Round vertex: the prism's other top corner (vertex 3 of the original four, now the last).
        Model v = m;
        REQUIRE(fillet_vertex(v, fin, static_cast<int>(t.polygon.size()) - 1, 2000, n));
        CHECK(v.regions[2].solids[0].polygon.size() == 13);  // its tangent point lands on the first arc's end: merged
        CHECK(build_error(must_write(tool_deck, v)) == "");
        CHECK_FALSE(fillet_vertex(v, fin, 99, 2000, n));
    }
    {
        Model m = m0;  // Taper: the top-right edge 1000 in along x; the right wall leans in by atan(1000 / 4000)
        REQUIRE(taper_edge(m, fin, top_right, {0, -1000}));
        CHECK(m.regions[2].solids[0].polygon[2] == Pt{5000, 9000});
        CHECK(region_volume(m, 2) == Catch::Approx((6000.0 * 4000 - 2'000'000) * 10000 * 1e-24).epsilon(1e-12));
        CHECK(build_error(round_trip(tool_deck, m)) == "");
        CHECK_FALSE(taper_edge(m, fin, top_right, {0, 0}));
        CHECK_FALSE(taper_edge(m, fin, top_right, {0, -20000}));  // the polygon folds
    }
    {
        Model m = m0;  // a nanowire of 16 sides through the silicon along y, beside the fin
        const int sides = NiTCAD::app::geo::polygon_sides(2000, 600);
        const int r = add_cylinder(m, 1, 0, 10000, {3500, 15000}, 2000, sides, "sio2");
        REQUIRE(r == 3);
        const Solid& t = m.regions[3].solids[0];
        CHECK(t.polygon.size() == 16);
        CHECK(region_volume(m, 3) == Catch::Approx(prism_cm3(t)).epsilon(1e-12));
        CHECK(build_error(round_trip(tool_deck, m)) == "");
        CHECK(add_cylinder(m, 1, 0, 10000, {3500, 15000}, 0, 16, "sio2") == -1);
        CHECK(add_cylinder(m, 1, 0, 10000, {3500, 15000}, 2000, 2, "sio2") == -1);
    }
}

TEST_CASE("app solid: Inset and Extrude, new solids on a face (Unit 26m MT6, MC6)", "[app][solid][faces]") {
    const Model m0 = must_read(tool_deck);
    const SolidRef fin{2, 0};
    // Face normals: a box's, a prism's caps and sides (counter-clockwise or not), none for a slanted side.
    CHECK(face_normal(m0.regions[2].solids[0], 5)->axis == 2);
    CHECK(face_normal(m0.regions[2].solids[0], 0)->sign == -1);
    {
        const Solid p{SolidKind::prism, {}, 2, 0, 1000, {{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}};
        CHECK((face_normal(p, 2)->axis == 1 && face_normal(p, 2)->sign == -1));  // the edge along +x: outward is -y
        CHECK((face_normal(p, 3)->axis == 0 && face_normal(p, 3)->sign == 1));
        Solid cw = p;
        std::ranges::reverse(cw.polygon);  // (0,1000) -> (1000,1000) -> ...: the first edge's outward is +y
        CHECK((face_normal(cw, 2)->axis == 1 && face_normal(cw, 2)->sign == 1));
        const Solid tri{SolidKind::prism, {}, 2, 0, 1000, {{0, 0}, {1000, 0}, {0, 1000}}};
        CHECK_FALSE(face_normal(tri, 3));
        // A cap: its polygon kept without a margin; a margin only on a rectangle.
        const auto cap = face_solid(tri, 1, 0, 500);
        REQUIRE(cap);
        CHECK((cap->kind == SolidKind::prism && cap->from == 1000 && cap->to == 1500 && cap->polygon == tri.polygon));
        CHECK_FALSE(face_solid(tri, 1, 100, 500));
        const auto inset = face_solid(p, 1, 100, -300);
        REQUIRE(inset);
        CHECK(inset->box == Box3{{100, 100, 700}, {900, 900, 1000}});
    }
    {
        Model m = m0;  // Extrude the fin's top by 0.1um as a new region: a full-face piece, so it builds
        CHECK(extrude_face(m, fin, 5, -1000, "silicon") == -1);
        const int r = extrude_face(m, fin, 5, 1000, "silicon");
        REQUIRE(r == 3);
        CHECK(m.regions[3].solids[0].box == Box3{{4000, 0, 5000}, {10000, 10000, 6000}});
        check_partition(m);
        const std::string out = round_trip(tool_deck, m);
        CHECK(build_error(out) == "");
    }
    {
        Model m = m0;  // Inset of the fin's own material joins its region (a bump, no interface between them)
        const int r = inset_face(m, fin, 5, 1000, 1000, "sio2");
        REQUIRE(r == 2);
        REQUIRE(m.regions[2].solids.size() == 2);
        CHECK(m.regions[2].solids[1].box == Box3{{5000, 1000, 5000}, {9000, 9000, 6000}});
        check_partition(m);
        const std::string out = round_trip(tool_deck, m);
        // The fin's top is then a ring around the bump: a holed interface piece, which meshes since Unit 64.
        CHECK(build_error(out) == "");
        // Pushed into the oxide from the domain's top face, the inset's pieces are its walls and its floor: it builds.
        Model ox = m0;
        REQUIRE(inset_face(ox, {1, 0}, 5, 1000, -1500, "insulator") == 3);
        CHECK(ox.regions[3].solids[0].box == Box3{{1000, 1000, 8500}, {19000, 9000, 10000}});
        std::string nitride = round_trip(tool_deck, ox);
        nitride.replace(nitride.find("material = insulator\n"), 21, "material = insulator\neps_r = 7.5\n");
        CHECK(build_error(nitride) == "");
        // Of another material on the fin it is a new region meeting the fin inside its face: refused alike.
        Model o = m0;
        REQUIRE(inset_face(o, fin, 5, 1000, -1000, "silicon") == 3);
        // Of another material on the fin it is a new region meeting the fin inside its face: it builds alike.
        CHECK(build_error(must_write(tool_deck, o)) == "");
        CHECK(inset_face(o, fin, 5, 3000, 1000, "sio2") == -1);  // the margin empties it
        CHECK(inset_face(o, fin, 9, 0, 1000, "sio2") == -1);
    }
}

TEST_CASE("app solid: keep_clear and the automatic s (Unit 26m MT12)", "[app][solid][meshable]") {
    const Model m0 = must_read(tool_deck);
    const Solid fin = m0.regions[2].solids[0];  // x 0.4-1.0um
    // A copy moved along x so its x_min falls 100 quanta past the fin's x_max: margin 400 brings it back to touch.
    const std::vector<SolidRef> none;
    CHECK(keep_clear(m0, {fin}, none, 0, 6100, 400) == 6000);
    CHECK(keep_clear(m0, {fin}, none, 0, 6300, 400) == 6400);  // 300 past: apart by the margin is nearer
    CHECK(keep_clear(m0, {fin}, none, 0, 7000, 400) == 7000);  // already clear
    CHECK(keep_clear(m0, {fin}, none, 0, 6100, 0) == 6100);
    // The domain's faces count: the fin moved to 50 quanta from x = 0 snaps onto it.
    CHECK(keep_clear(m0, {fin}, {{2, 0}}, 0, -3950, 400) == -4000);
    // A solid far off on the other axes does not count.
    Model far = m0;
    REQUIRE(add_box(far, {{12000, 0, 0}, {13000, 10000, 500}}, "sio2") >= 0);
    CHECK(keep_clear(far, {fin}, {{2, 0}}, 2, 3800, 400) == 3800);  // z 0.48-0.88um: nothing within 4 s on x

    // The automatic s: the largest the rules take; 1 % more is refused by the build's own check.
    const auto max = mesh_length(tool_deck, "max_size");
    REQUIRE(max);
    CHECK(*max == Catch::Approx(0.25e-4));
    CHECK(mesh_length(tool_deck, "interface_resolution") == Catch::Approx(1e-6));
    CHECK_FALSE(mesh_length(tool_deck, "refine"));
    const auto s = automatic_resolution(m0, *max);
    REQUIRE(s);
    INFO("s = " << *s);
    // The tightest gap is the fin's 0.1um from the oxide and from the bottom: 4 s <= 0.1um.
    CHECK(*s <= 0.25e-5 * 1.0000001);
    CHECK(*s >= 0.25e-5 * 0.998);
    const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(m0));
    REQUIRE(complex);
    NiTCAD::mesh::MeshOptions3D o;
    o.max_size = *max;
    o.interface_resolution = *s;
    CHECK(NiTCAD::mesh::generate_mesh(*complex, o).has_value());
    o.interface_resolution = *s * 1.01;
    CHECK_FALSE(NiTCAD::mesh::generate_mesh(*complex, o).has_value());
    // Without interfaces any s goes.
    Model plain = m0;
    plain.regions.resize(1);
    CHECK(std::isinf(*automatic_resolution(plain, *max)));
}

TEST_CASE("app solid: mesh_length and set_mesh_length take one [mesh] and keep the deck's last line end (Unit 26m MT12)", "[app][solid][meshable][meshlen]") {
    // Two [mesh] sections (the build refuses such a deck; the editor may hold one): both read and write the last, as
    // read_model and write_model do.
    const std::string two = "[mesh]\nmax_size = 0.1um\n\n[mesh]\nmax_size = 0.2um\n";
    CHECK(mesh_length(two, "max_size") == Catch::Approx(0.2e-4));
    const auto set = set_mesh_length(two, "max_size", 3000);
    REQUIRE(set);
    CHECK(*set == "[mesh]\nmax_size = 0.1um\n\n[mesh]\nmax_size = 0.3um\n");
    CHECK(mesh_length(*set, "max_size") == Catch::Approx(0.3e-4));
    // An entry continued over lines, last in a deck without a final newline: replaced, still without one.
    const std::string tail = "[mesh]\ntype = generate3d\ninterface_resolution = \\\n    0.01um";
    const auto one = set_mesh_length(tail, "interface_resolution", 50);
    REQUIRE(one);
    CHECK(*one == "[mesh]\ntype = generate3d\ninterface_resolution = 0.005um");
    // With CRLF ends the new line keeps them.
    const auto crlf = set_mesh_length("[mesh]\r\nmax_size = 0.1um\r\n", "interface_resolution", 100);
    REQUIRE(crlf);
    CHECK(*crlf == "[mesh]\r\nmax_size = 0.1um\r\ninterface_resolution = 0.01um\r\n");
}

namespace {

// The MT12 table: tool_deck's layouts drawn with the tools, and the max_size each is meshed at.
struct Layout {
    const char* name;
    Model m;
    double h;  // cm
};
std::vector<Layout> mt12_table() {
    const Model m0 = must_read(tool_deck);
    const SolidRef fin{2, 0};
    const int top_right = box_edge(m0.regions[2].solids[0], 1 | 4, 2);
    std::vector<Layout> out;
    out.push_back({"fin", m0, 0.25e-4});
    out.push_back({"fin, h 0.15um", m0, 0.15e-4});
    Model c = m0;
    REQUIRE(chamfer_edge(c, fin, top_right, 1000));
    out.push_back({"chamfer", c, 0.25e-4});
    Model f = m0;
    REQUIRE(fillet_edge(f, fin, top_right, 2000, 5));
    out.push_back({"fillet", f, 0.25e-4});
    Model w = m0;
    REQUIRE(add_cylinder(w, 1, 0, 10000, {3500, 15000}, 2000, 16, "sio2") >= 0);
    out.push_back({"nanowire", w, 0.25e-4});
    Model a = m0;
    REQUIRE(set_box(a, fin, {{4000, 0, 1000}, {7000, 10000, 5000}}));
    REQUIRE(array_solid(a, fin, {5000, 0, 0}, 3) >= 0);
    out.push_back({"three fins", a, 0.25e-4});
    Model e = m0;
    REQUIRE(extrude_face(e, fin, 5, 1000, "silicon") >= 0);
    out.push_back({"extruded gate", e, 0.25e-4});
    return out;
}

}  // namespace

TEST_CASE("app solid: the automatic s on a table of drawn layouts; 1 % more is refused (Unit 26m MT12, MC13)", "[app][solid][meshable]") {
    for (const Layout& l : mt12_table()) {
        INFO(l.name);
        const auto s = automatic_resolution(l.m, l.h);
        REQUIRE(s);
        REQUIRE(std::isfinite(*s));
        const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(l.m));
        REQUIRE(complex);
        NiTCAD::mesh::MeshOptions3D o;
        o.max_size = l.h;
        // The automatic s meshes; 1 % more the generator refuses.
        o.interface_resolution = *s;
        CHECK(NiTCAD::mesh::generate_mesh(*complex, o).has_value());
        o.interface_resolution = *s * 1.01;
        CHECK_FALSE(NiTCAD::mesh::generate_mesh(*complex, o).has_value());
        // The rules alone: their largest s passes them, 1 % more does not; the automatic s is no larger.
        const auto rules = NiTCAD::mesh::largest_resolution(*complex, o);
        REQUIRE(rules);
        CHECK(*s <= *rules);
        o.interface_resolution = *rules;
        CHECK(NiTCAD::mesh::check_rules(*complex, o).has_value());
        o.interface_resolution = *rules * 1.01;
        CHECK_FALSE(NiTCAD::mesh::check_rules(*complex, o).has_value());
    }
}

// The estimate against built meshes at the deck's s (recorded, not a gate: MC13 says ±30 %); also the fit's data.
TEST_CASE("app solid: the node estimate against built meshes (Unit 26m MT12, MC13; recorded)", "[.][app][solid][estimate]") {
    for (const Layout& l : mt12_table()) {
        const double s = 1e-6;
        const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(l.m));
        REQUIRE(complex);
        NiTCAD::mesh::MeshOptions3D o;
        o.max_size = l.h;
        o.interface_resolution = s;
        const auto built = NiTCAD::mesh::generate_mesh(*complex, o);
        INFO(l.name);
        REQUIRE(built);
        const auto e = extent(l.m);
        REQUIRE(e);
        const double hp = std::min(l.h, 4 * s);
        const double n = static_cast<double>(built->mesh.node_count()), est = node_estimate(*e, l.h, s);
        WARN(std::format("estimate {}: V/h3 {:.1f}, A/hp2 {:.1f}, L/hp {:.1f}, nodes {}, estimate {:.0f} ({:+.1f} %)", l.name, e->volume / (l.h * l.h * l.h),
                         e->interface_area / (hp * hp), e->feature_length / hp, n, est, 100.0 * (est - n) / n));
        CHECK(std::abs(est - n) <= 0.3 * n);  // MC13: ±30 % (in-sample: the constants were fitted on this table)
    }
}

// MC13's layouts: seed 26231, drawn in order with the tools on a 2 x 1 x 1 µm silicon block (one to four solids each:
// a box, chamfered, filleted, extruded or arrayed, or a cylinder along y).
std::vector<Model> mc13_layouts(int count) {
    std::mt19937 rng(26231);
    const auto pick = [&](int lo, int hi) { return std::uniform_int_distribution<int>(lo, hi)(rng); };
    Model base;
    base.mesh = MeshKind::layered;
    base.domain = {{0, 0, 0}, {20000, 10000, 10000}};
    base.regions.push_back({"si", "silicon", {}, -1});
    std::vector<Model> out;
    for (int n = 0; n < count; ++n) {
        Model m = base;
        const int solids = pick(1, 4);
        for (int k = 0; k < solids; ++k) {
            const int tool = pick(0, 5);
            const P3 lo{500 * pick(0, 36), 500 * pick(0, 16), 500 * pick(0, 16)};
            const P3 hi{std::min<Coord>(20000, lo.x + 500 * pick(2, 12)), std::min<Coord>(10000, lo.y + 500 * pick(2, 12)),
                        std::min<Coord>(10000, lo.z + 500 * pick(2, 8))};
            if (tool == 5) {  // a cylinder along y through the device
                const Coord r = 250 * pick(2, 6);
                const Pt c{std::clamp<Coord>(lo.z, r + 500, 10000 - r - 500), std::clamp<Coord>(lo.x, r + 500, 20000 - r - 500)};
                add_cylinder(m, 1, 0, 10000, c, r, 16, "sio2");
                continue;
            }
            const int r = add_box(m, {lo, hi}, "sio2");
            if (r < 0) continue;
            const SolidRef s{r, 0};
            const int edge = pick(0, 11);
            if (tool == 1) chamfer_edge(m, s, edge, 250 * pick(1, 4));
            if (tool == 2) fillet_edge(m, s, edge, 250 * pick(2, 6), pick(1, 6));
            if (tool == 3) extrude_face(m, s, pick(0, 5), 500 * pick(1, 4), "silicon");
            if (tool == 4) array_solid(m, s, {static_cast<Coord>(500 * pick(4, 10)), 0, 0}, pick(2, 3));
        }
        out.push_back(std::move(m));
    }
    return out;
}

// MC13: 500 seeded random layouts drawn with the tools. Accepted (the automatic s found: it meshes there): the generator
// meshes it, and refuses 1 % more. Refused (no s meshes): the generator refuses it too. Counts are printed for the
// record; hidden, as it meshes every layout (sharded by NITCAD_SHARD=k/n over processes; NITCAD_LAYOUT=n prints one).
TEST_CASE("app solid: 500 random drawn layouts, accepted ones mesh, refused ones refused alike (Unit 26m MT12, MC13)", "[.][app][solid][mc13]") {
    const std::vector<Model> layouts = mc13_layouts(500);
    int accepted = 0, refused = 0, false_accepts = 0, false_refusals = 0, mismatched = 0, above_passes = 0;
    std::map<std::string, int> reasons;
    const auto env = [](const char* name) {
        char* value = nullptr;
        std::size_t length = 0;
        std::string text;
        if (_dupenv_s(&value, &length, name) == 0 && value != nullptr) text = value;
        std::free(value);
        return text;
    };
    int shard = 0, shards = 1, only = -1;
    if (const std::string t = env("NITCAD_SHARD"); t.find('/') != std::string::npos) {
        shard = std::stoi(t.substr(0, t.find('/')));
        shards = std::max(1, std::stoi(t.substr(t.find('/') + 1)));
    }
    if (const std::string t = env("NITCAD_LAYOUT"); !t.empty()) only = std::stoi(t);
    const auto stage = [&](const char* what) {
        if (only < 0) return;
        std::printf("MC13STAGE %s\n", what);
        std::fflush(stdout);
    };
    for (int n = 0; n < 500; ++n) {
        const Model& m = layouts[static_cast<std::size_t>(n)];
        if (n % shards != shard || (only >= 0 && n != only)) continue;
        if (only >= 0) {  // the layout's solids, to rebuild it elsewhere
            for (std::size_t r = 1; r < m.regions.size(); ++r) {
                for (const Solid& s : m.regions[r].solids) {
                    std::printf("MC13SOLID region %zu %s: ", r, m.regions[r].material.c_str());
                    if (s.kind == SolidKind::box) {
                        std::printf("box %lld %lld %lld  %lld %lld %lld\n", s.box.lo.x, s.box.lo.y, s.box.lo.z, s.box.hi.x, s.box.hi.y, s.box.hi.z);
                    } else {
                        std::printf("prism axis %d from %lld to %lld:", s.axis, s.from, s.to);
                        for (const Pt& p : s.polygon) std::printf(" (%lld,%lld)", p.x, p.y);
                        std::printf("\n");
                    }
                }
            }
        }
        if (shards > 1 || only >= 0) {  // the log names the layout it is on (a crash ends the process)
            std::printf("MC13AT %d\n", n);
            std::fflush(stdout);
        }
        const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(m));
        if (!complex) {  // the layers' own refusal (a solid leaving the domain, ...): not the generator's
            ++reasons["partition: " + complex.error().message.substr(0, 60)];
            continue;
        }
        NiTCAD::mesh::MeshOptions3D o;
        o.max_size = 0.25e-4;  // max_nodes as automatic_resolution meshes (the default)
        stage("automatic s");
        const auto s = automatic_resolution(m, o.max_size);
        if (s && std::isinf(*s)) continue;  // no interfaces
        if (!s) {
            ++refused;
            ++reasons[s.error().substr(0, 70)];
            o.interface_resolution = 0.5e-6;
            stage("build at 0.005um");
            const auto built = NiTCAD::mesh::generate_mesh(*complex, o);
            const auto rules = NiTCAD::mesh::check_rules(*complex, o);
            if (built) {
                ++false_refusals;
                std::printf("MC13NOTE layout %d refused (%s) but meshes at 0.005um\n", n, s.error().c_str());
            } else if (!rules && rules.error().message != s.error()) {
                // A refusal by the rules: a rule that does not scale with s says the same at any s. (A refusal by
                // builds is not compared: the generator's reason at another s may differ.)
                ++mismatched;
                std::printf("MC13NOTE layout %d refused: %s | rules at 0.005um: %s\n", n, s.error().c_str(), rules.error().message.c_str());
            }
            continue;
        }
        ++accepted;
        o.interface_resolution = *s;
        stage("build at s");
        const auto built = NiTCAD::mesh::generate_mesh(*complex, o);
        if (!built) {
            ++false_accepts;
            std::printf("MC13NOTE layout %d accepted at s = %g but: %s\n", n, *s, built.error().message.c_str());
        }
        // 1 % above: the generator refuses (by its rules, or a check only the mesh makes).
        o.interface_resolution = *s * 1.01;
        stage("build at 1.01 s");
        if (NiTCAD::mesh::generate_mesh(*complex, o)) {
            ++above_passes;
            std::printf("MC13NOTE layout %d meshes at 1.01 s = %g\n", n, *s * 1.01);
        }
    }
    std::string why;
    for (const auto& [r, k] : reasons) why += std::format("\n  {} x {}", k, r);
    WARN(std::format("MC13: accepted {}, refused {}, false accepts {}, false refusals {}, messages differing {}, 1.01 s passing {}{}",
                     accepted, refused, false_accepts, false_refusals, mismatched, above_passes, why));
    std::printf("MC13SUM %d %d %d %d %d %d\n", accepted, refused, false_accepts, false_refusals, mismatched, above_passes);
    for (const auto& [r, k] : reasons) std::printf("MC13WHY %d|%s\n", k, r.c_str());
    CHECK(false_accepts == 0);
    // mismatched: informational, not MC13 (a refusal by one rule at the automatic s and by another at 0.005um are
    // both the rules'; which rule fires first depends on s).
    CHECK(above_passes == 0);
}

// MC13 layout 284: a 16-sided cylinder along y whose axis lies on a box's face. At s = 0.5 nm the rules pass and a
// pair set's 2D mesh met two copies of one vertex an ulp apart: splitting the segment between them inserted a point
// twice, ending in the 2D triangulation's contract. Now refused as a degenerate mesh; the automatic s search ends too.
TEST_CASE("app solid: MC13 layout 284, a cylinder halved by a box face, is refused, not a crash (Unit 26m MC13)", "[app][solid][meshable][l284]") {
    const Model m = mc13_layouts(285).back();
    REQUIRE(m.regions.size() == 3);
    const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(m));
    REQUIRE(complex);
    NiTCAD::mesh::MeshOptions3D o;
    o.max_size = 0.25e-4;
    o.interface_resolution = 5e-8;
    REQUIRE(NiTCAD::mesh::check_rules(*complex, o));
    const auto built = NiTCAD::mesh::generate_mesh(*complex, o);
    REQUIRE_FALSE(built);
    // The pair set's refusal, as the 3D generator reports it (invalid_input: a finer s does not cure it).
    CHECK(built.error().code == NiTCAD::base::ErrorCode::invalid_input);
    CHECK(built.error().message.find("the pair set of an interface could not be meshed") != std::string::npos);
    CHECK(built.error().message.find("too short to split") != std::string::npos);
    const auto s = automatic_resolution(m, o.max_size);  // ends, with an s that meshes or a refusal
    if (s) {
        o.interface_resolution = *s;
        CHECK(NiTCAD::mesh::generate_mesh(*complex, o));
    }
}

// The automatic s's cost per table layout (NITCAD_TABLE=i picks one; recorded for MT12's button, hidden).
TEST_CASE("app solid: the automatic s's time on each table layout (Unit 26m MT12; recorded)", "[.][app][solid][autotime]") {
    char* value = nullptr;
    std::size_t length = 0;
    int only = -1;
    if (_dupenv_s(&value, &length, "NITCAD_TABLE") == 0 && value != nullptr) only = std::stoi(value);
    std::free(value);
    const auto table = mt12_table();
    for (std::size_t i = 0; i < table.size(); ++i) {
        if (only >= 0 && static_cast<int>(i) != only) continue;
        const auto complex = NiTCAD::mesh::layered_partition3d(to_layered(table[i].m));
        REQUIRE(complex);
        NiTCAD::mesh::MeshOptions3D o;
        o.max_size = table[i].h;
        auto t0 = std::chrono::steady_clock::now();
        const auto rules = NiTCAD::mesh::largest_resolution(*complex, o);
        const double rules_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        t0 = std::chrono::steady_clock::now();
        const auto s = automatic_resolution(table[i].m, table[i].h);
        const double auto_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("AUTOTIME %zu %s: rules s %.4g cm in %.2f s; automatic s %s in %.2f s\n", i, table[i].name, rules ? *rules : -1.0, rules_s,
                    s ? std::format("{:.4g} cm", *s).c_str() : s.error().c_str(), auto_s);
        std::fflush(stdout);
    }
}

// --- Unit 26m slice 3 (MT8, MC7): parameters in the modeller ------------------------------------------------------

namespace {

const std::string tied_deck = R"([parameters]
x0 = 30nm
Lg = 20nm

[mesh]
type = generate3d
layered = on
domain = 0 100nm 0 60nm 0 14nm
patch body = 0 0 0  100nm 60nm 0
max_size = 0.05um
interface_resolution = 2nm

[region si]
material = silicon

[region gate]
material = sio2
box = $x0 $x0 + $Lg 40nm 50nm 0 14nm

[region cap]
material = sio2
box = $x0 $x0 + $Lg 50nm 60nm 0 14nm

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";

constexpr Coord nm = 10;  // quanta per nm

bool has(const std::string& s, std::string_view part) { return s.find(part) != std::string::npos; }

std::size_t changed_lines(const std::string& a, const std::string& b) {
    std::vector<std::string> x, y;
    for (const auto& [s, v] : {std::pair{&a, &x}, std::pair{&b, &y}}) {
        std::size_t at = 0;
        while (at <= s->size()) {
            const auto end = s->find('\n', at);
            v->push_back(s->substr(at, end == std::string::npos ? std::string::npos : end - at));
            if (end == std::string::npos) break;
            at = end + 1;
        }
    }
    std::size_t n = x.size() == y.size() ? 0 : 1000;
    for (std::size_t k = 0; k < std::min(x.size(), y.size()); ++k) n += x[k] != y[k];
    return n;
}

}  // namespace

TEST_CASE("parameters: the modeller reads ties, writes them back, and Lg moves every tied face (MC7)", "[app][solid][parameters]") {
    auto m = read_model(tied_deck);
    REQUIRE(m);
    REQUIRE(m->editable());
    REQUIRE(m->parameters.size() == 2);
    CHECK(m->parameters[1].name == "Lg");
    CHECK(m->parameters[1].expr == "20nm");
    CHECK(m->parameters[1].value == 20 * nm);
    const Solid& gate = m->regions[1].solids[0];
    CHECK(gate.box.lo.x == 30 * nm);
    CHECK(gate.box.hi.x == 50 * nm);
    REQUIRE(gate.refs.size() == 6);
    CHECK(gate.refs[0] == "$x0");
    CHECK(gate.refs[1] == "$x0 + $Lg");
    CHECK(gate.refs[2].empty());
    // Unchanged: the text byte for byte.
    CHECK(*write_model(tied_deck, *m) == tied_deck);

    SECTION("changing Lg moves both tied faces; only its line changes") {
        REQUIRE(set_parameter(*m, "Lg", "25nm"));
        CHECK(m->regions[1].solids[0].box.hi.x == 55 * nm);
        CHECK(m->regions[2].solids[0].box.hi.x == 55 * nm);
        CHECK(m->regions[1].solids[0].box.lo.x == 30 * nm);
        const auto out = write_model(tied_deck, *m);
        REQUIRE(out);
        CHECK(changed_lines(tied_deck, *out) == 1);
        CHECK(has(*out, "Lg = 25nm\n"));
        CHECK(has(*out, "box = $x0 $x0 + $Lg 40nm 50nm 0 14nm"));
        // The written deck reads back to the same model, and builds as the deck with values written out would.
        CHECK(*read_model(*out) == *m);
    }
    SECTION("a parameter's errors leave the model as it was") {
        const Model before = *m;
        CHECK_FALSE(set_parameter(*m, "Lg", "$Lg + 1nm"));          // a cycle
        CHECK_FALSE(set_parameter(*m, "Lg", "$nope"));              // undefined
        CHECK_FALSE(set_parameter(*m, "Lg", "silicon"));            // neither
        CHECK_FALSE(set_parameter(*m, "Lg", "-40nm"));              // the gate's far face past its near one
        CHECK_FALSE(set_parameter(*m, "2x", "1nm"));                // a bad name
        CHECK(*m == before);
        REQUIRE(set_parameter(*m, "W", "5nm"));                     // a new one, at the end
        CHECK(m->parameters.back().name == "W");
        CHECK(has(*write_model(tied_deck, *m), "Lg = 20nm\nW = 5nm\n"));
    }
    SECTION("dragging a tied face edits its parameter: the far face sets Lg, every tied face follows") {
        const Model before = *m;
        REQUIRE(move_face(*m, {1, 0}, 1, 5 * nm));
        tie_edits(before, *m);
        CHECK(m->parameters[1].value == 25 * nm);
        CHECK(m->parameters[0].value == 30 * nm);
        CHECK(m->regions[2].solids[0].box.hi.x == 55 * nm);
        CHECK(m->regions[1].solids[0].refs[1] == "$x0 + $Lg");
    }
    SECTION("the near face sets x0; the far face stays (Lg takes the rest)") {
        const Model before = *m;
        REQUIRE(move_face(*m, {1, 0}, 0, -10 * nm));
        tie_edits(before, *m);
        CHECK(m->parameters[0].value == 20 * nm);
        CHECK(m->parameters[1].value == 30 * nm);
        CHECK(m->regions[1].solids[0].box.lo.x == 20 * nm);
        CHECK(m->regions[1].solids[0].box.hi.x == 50 * nm);
        CHECK(m->regions[2].solids[0].box.lo.x == 20 * nm);  // the cap's tied faces too
        CHECK(m->regions[2].solids[0].box.hi.x == 50 * nm);
    }
    SECTION("moving the solid moves x0 only") {
        const Model before = *m;
        REQUIRE(move_solid(*m, {1, 0}, {10 * nm, 0, 0}));
        tie_edits(before, *m);
        CHECK(m->parameters[0].value == 40 * nm);
        CHECK(m->parameters[1].value == 20 * nm);
        CHECK(m->regions[2].solids[0].box.lo.x == 40 * nm);
    }
    SECTION("a tie whose parameter is an expression is broken, the edit kept") {
        REQUIRE(set_parameter(*m, "x0", "10nm + 20nm"));
        const Model before = *m;
        REQUIRE(move_face(*m, {1, 0}, 0, 5 * nm));
        tie_edits(before, *m);
        CHECK(m->regions[1].solids[0].refs[0].empty());
        CHECK(m->regions[1].solids[0].box.lo.x == 35 * nm);
        CHECK(m->regions[2].solids[0].box.lo.x == 30 * nm);  // the cap keeps its tie
        CHECK(has(*write_model(tied_deck, *m), "box = " + NiTCAD::app::geo::format_length(35 * nm) + " $x0 + $Lg " + NiTCAD::app::geo::format_length(40 * nm) + " " + NiTCAD::app::geo::format_length(50 * nm) + " 0 " + NiTCAD::app::geo::format_length(14 * nm)));
    }
    SECTION("N names an extent: a new parameter, the far face written near + $name") {
        REQUIRE(name_extent(*m, {1, 0}, 4, 5, "H"));  // the gate's z extent
        CHECK(m->parameters.back().name == "H");
        CHECK(m->parameters.back().value == 14 * nm);
        CHECK(m->regions[1].solids[0].refs[5] == "0 + $H");
        const auto out = write_model(tied_deck, *m);
        CHECK(has(*out, "H = " + NiTCAD::app::geo::format_length(14 * nm) + "\n"));
        CHECK(has(*out, "box = $x0 $x0 + $Lg " + NiTCAD::app::geo::format_length(40 * nm) + " " + NiTCAD::app::geo::format_length(50 * nm) + " 0 0 + $H"));
        // Naming the cap's z extent with the same H ties it too; changing H moves both tops.
        REQUIRE(name_extent(*m, {2, 0}, 4, 5, "H"));
        REQUIRE(set_parameter(*m, "H", "12nm"));
        CHECK(m->regions[1].solids[0].box.hi.z == 12 * nm);
        CHECK(m->regions[2].solids[0].box.hi.z == 12 * nm);
        CHECK_FALSE(name_extent(*m, {1, 0}, 4, 5, "bad name"));
    }
    SECTION("a low face named: written near - $name, and dragging it sets the parameter by the minus") {
        REQUIRE(name_extent(*m, {1, 0}, 5, 4, "B"));  // the gate's bottom from its top
        CHECK(m->regions[1].solids[0].refs[4] == NiTCAD::app::geo::format_length(14 * nm) + " - $B");
        CHECK(m->parameters.back().value == 14 * nm);
        const Model before = *m;
        REQUIRE(move_face(*m, {1, 0}, 4, 2 * nm));  // the bottom up 2nm: B is 12nm
        tie_edits(before, *m);
        CHECK(m->regions[1].solids[0].refs[4] == NiTCAD::app::geo::format_length(14 * nm) + " - $B");
        CHECK(m->parameters.back().value == 12 * nm);
        CHECK(m->regions[1].solids[0].box.lo.z == 2 * nm);
    }
}

TEST_CASE("parameters: a tied drag frame with 100 tied solids adds at most 1 ms (MC9, MT8)", "[app][solid][parameters][perf]") {
    // 100 boxes in a row, each one's far x face `$x<k> + $W` (one W for all): a frame of the drag of box 0's far face
    // is what SolidPage::after_drag_step does: the edit, tie_edits, write_model and the deck line.
    std::string deck = "[parameters]\nW = 4nm\n";
    for (int k = 0; k < 100; ++k) deck += std::format("x{} = {}nm\n", k, 10 * k);
    deck += "\n[mesh]\ntype = generate3d\nlayered = on\ndomain = 0 1010nm 0 60nm 0 14nm\nmax_size = 0.05um\ninterface_resolution = 1nm\n\n"
            "[region si]\nmaterial = silicon\n";
    for (int k = 0; k < 100; ++k) deck += std::format("\n[region b{}]\nmaterial = sio2\nbox = $x{} $x{} + $W 0 10nm 0 10nm\n", k, k, k);
    deck += "\n[equilibrium]\n";
    const auto m0 = read_model(deck);
    REQUIRE(m0);
    REQUIRE(m0->editable());
    std::vector<double> ms;
    std::string last;
    for (int f = 0; f < 60; ++f) {
        const auto t0 = std::chrono::steady_clock::now();
        Model m = *m0;
        REQUIRE(move_face(m, {1, 0}, 1, (1 + f % 4) * 10));
        tie_edits(*m0, m);
        const auto out = write_model(deck, m);
        REQUIRE(out);
        last = NiTCAD::app::modeller::deck_line(deck, *out);
        ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    }
    std::ranges::sort(ms);
    std::printf("MC9 slice 3, 100 tied solids: a tied drag frame median %.3f ms (max %.3f)\n", ms[ms.size() / 2], ms.back());
    CHECK(has(last, "W = "));  // W moved: the first line written is its
    // Every box's far face followed W.
    Model m = *m0;
    REQUIRE(move_face(m, {1, 0}, 1, 40));
    tie_edits(*m0, m);
    for (int k = 0; k < 100; ++k) CHECK(m.regions[static_cast<std::size_t>(k) + 1].solids[0].box.hi.x - m.regions[static_cast<std::size_t>(k) + 1].solids[0].box.lo.x == 80);
    // The CI runner records it (owner, 2026-11-22, as SE5 and KG6 (b)): there it took 1.055 ms on PR #116, four tests
    // sharing its four cores; this workstation is the gate.
    if (!app_test::on_ci()) CHECK(ms[ms.size() / 2] <= 1.0);
}

// --- Unit 26m slice 4 (MT9, MT10, MC8): doping and mesh tools -------------------------------------------------------

TEST_CASE("doping and mesh tools: each writes exactly its line, and the built doping is the profiles' formula (MC8)",
          "[app][solid][doping]") {
    const auto m = read_model(well_deck);
    REQUIRE(m);
    REQUIRE(m->editable());
    using NiTCAD::device::Dopant;
    using NiTCAD::device::DepthShape;
    std::string deck = well_deck;
    const auto add = [&](const NiTCAD::device::DopingProfile& p, const std::string& line) {
        const auto out = add_profile(deck, p);
        REQUIRE(out);
        INFO(*out);
        CHECK(has(*out, "\n" + line + "\n"));
        // Minimal edit: one line more, the rest byte for byte.
        CHECK(out->size() == deck.size() + line.size() + 1);
        deck = *out;
    };
    // Region profile on the well (a prism along x over y 0..0.6um, z 0..0.2um): it spans the domain along y, so y is
    // the depth axis and x, z the windows.
    const auto region = region_profile(*m, {2, 0}, Dopant::donor, 1e18);
    CHECK(region.shape == DepthShape::constant);
    CHECK(region.axis == 1);
    add(region, "profile = donor constant 1e+18 axis=y x=2e-05cm:8e-05cm z=0cm:2e-05cm");
    // Implant through the top (z = 0.6um, looking up z): its peak 50nm down, σ 25nm, over x 0.3..0.7um, y 0.1..0.5um.
    const Box3 window{{3000, 1000, 6000}, {7000, 5000, 6000}};
    const auto implant = implant_profile({2, 1}, 6000, window, 50 * nm, 10 * nm, Dopant::acceptor, 1e19);
    CHECK(implant.center == 5500 / 1e8);
    CHECK(implant.width == 500 / 1e8 / 2.0);
    add(implant, "profile = acceptor gaussian 1e+19 axis=z center=" + std::format("{}", implant.center) + "cm width=" +
                     std::format("{}", implant.width) + "cm x=3e-05cm:7e-05cm y=1e-05cm:5e-05cm lateral=1e-06cm");
    // Junction from the bottom (z = 0, looking down z) meeting the deck's 1e16 acceptors 100nm up.
    const Box3 bottom{{0, 0, 0}, {10000, 6000, 0}};
    const auto existing = profiles(deck);
    REQUIRE(existing.size() == 3);  // the deck's own, the two above
    const auto junction = junction_profile({existing.front()}, {2, -1}, 0, bottom, 100 * nm, Dopant::donor, 1e19);
    REQUIRE(junction);
    CHECK(junction->shape == DepthShape::erfc);
    const NiTCAD::mesh::Point xj{0.5e-4, 0.3e-4, 100 * nm / 1e8};
    CHECK(std::abs(NiTCAD::device::evaluate(*junction, xj) - 1e16) <= 1e-9 * 1e16);
    CHECK_FALSE(junction_profile({existing.front()}, {2, 1}, 6000, bottom, 100 * nm, Dopant::donor, 1e19));   // looking up
    CHECK_FALSE(junction_profile({existing.front()}, {2, -1}, 0, bottom, 100 * nm, Dopant::acceptor, 1e19));  // no opposite
    CHECK_FALSE(junction_profile({existing.front()}, {2, -1}, 0, bottom, 100 * nm, Dopant::donor, 1e15));     // below it
    add(*junction, "profile = " + NiTCAD::deck::profile_text(*junction));
    CHECK(has(NiTCAD::deck::profile_text(*junction), "donor erfc 1e+19 axis=z center="));
    // The profiles read back as written.
    const auto back = profiles(deck);
    REQUIRE(back.size() == 4);
    CHECK(NiTCAD::deck::profile_text(back[1]) == NiTCAD::deck::profile_text(region));
    CHECK(NiTCAD::deck::profile_text(back[2]) == NiTCAD::deck::profile_text(implant));
    CHECK(NiTCAD::deck::profile_text(back[3]) == NiTCAD::deck::profile_text(*junction));
    // Max size and a refine box: the [mesh] lines.
    auto meshed = set_mesh_length(deck, "max_size", 1200);
    REQUIRE(meshed);
    CHECK(has(*meshed, "max_size = 0.12um\n"));
    const SizeBox3 refine{{{2000, 0, 4000}, {8000, 6000, 6000}}, 400};
    meshed = add_size_box(*meshed, refine);
    REQUIRE(meshed);
    CHECK(has(*meshed, "box = 0.2um 0 0.4um  0.8um 0.6um 0.6um  0.04um\n"));
    REQUIRE(size_boxes(*meshed).size() == 1);
    CHECK(size_boxes(*meshed)[0] == refine);
    CHECK_FALSE(add_size_box(*meshed, {{{0, 0, 0}, {0, 10, 10}}, 10}));
    // MC8: the built device's net doping at every node is the profiles' formula, to 1e-12 relative.
    const auto parsed = NiTCAD::deck::parse(*meshed);
    REQUIRE(parsed);
    const auto plan = NiTCAD::deck::build(*parsed);
    REQUIRE(plan);
    const auto& dev = plan->device;
    const auto pts = dev.mesh().points();
    const auto all = profiles(*meshed);
    std::size_t checked = 0;
    for (std::size_t k = 0; k < pts.size(); ++k) {
        if (NiTCAD::device::is_insulator(dev.regions()[dev.node_region()[k]])) continue;
        double nd = 0.0, na = 0.0;
        for (const auto& p : all) (p.dopant == Dopant::donor ? nd : na) += NiTCAD::device::evaluate(p, pts[k]);
        CHECK(std::abs(dev.donors()[k] - nd) <= 1e-12 * std::max(1.0, nd));
        CHECK(std::abs(dev.acceptors()[k] - na) <= 1e-12 * std::max(1.0, na));
        ++checked;
    }
    CHECK(checked > 100);
    std::printf("MC8: %zu semiconductor nodes checked against the formulas\n", checked);
}
