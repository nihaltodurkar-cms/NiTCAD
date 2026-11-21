// The 3D solid modeller's model (Unit 26l slice 1; DECISIONS.md SD3, SE3-SE5, SE10): the deck read and written by
// minimal edits, extrude3d decks read as prisms and converted, edits undone byte for byte, the write's time, snapping.
#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/mesh/complex.hpp"

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
    std::erase(t, '');
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
    for (const char* k : {"box =", "prism =", "material =", "[region ", "[contact ", "patch", "layered =", "domain =", "kind =", "type =",
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
                for (std::size_t r = 0; r < m.regions.size(); ++r) {
                    for (const Solid& s : m.regions[r].solids) {
                        bool in = false;
                        if (s.kind == SolidKind::box) {
                            in = c[0] > s.box.lo.x && c[0] < s.box.hi.x && c[1] > s.box.lo.y && c[1] < s.box.hi.y && c[2] > s.box.lo.z && c[2] < s.box.hi.z;
                        } else {
                            in = c[s.axis] > s.from && c[s.axis] < s.to && inside(s.polygon, c[(s.axis + 1) % 3], c[(s.axis + 2) % 3]);
                        }
                        if (in) region = static_cast<int>(r);
                    }
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
        CHECK(write_ms / done <= 2.0);
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
    CHECK(ms <= 2.0);
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
