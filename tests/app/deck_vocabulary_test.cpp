// The deck's vocabulary (Unit 26b, XD2, XD3, XG2): the same keys as deck/build.hpp's reference and
// deck::build, misspellings suggested back from real check errors, completion at every caret context.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "NiTCAD/app/deck_sections.hpp"
#include "NiTCAD/app/deck_vocabulary.hpp"
#include "NiTCAD/app/document.hpp"
#include "NiTCAD/io/files.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;

namespace {

std::string read_file(const std::filesystem::path& p) {
    std::ifstream f(p, std::ios::binary);
    std::string s{std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
    std::erase(s, '\r');
    return s;
}

// A deck that builds, for keys of a section kind in a context, and the section's header to add to.
struct Base {
    std::string text;
    std::string header;  // "contact anode", "mesh", ...
};

std::string diode_with(const std::string& more) { return app_test::diode_sweep() + more; }

Base base_for(std::string_view kind, std::string_view context) {
    const std::string diode = app_test::diode_sweep();
    if (kind == "mesh") {
        if (context == "generate2d") return {app_test::generate2d, "mesh"};
        if (context == "generate3d") return {app_test::generate3d, "mesh"};
        if (context == "extrude3d") return {app_test::extrude3d, "mesh"};
        if (context == "file") {
            // The generated 2D resistor's mesh, saved once.
            static const std::filesystem::path file = [] {
                const auto p = app_test::scratch("vocabulary_mesh") / "resistor.nitmesh";
                auto r = app::check_deck(app_test::generate2d, {});
                REQUIRE(r.plan);
                REQUIRE(r.plan->generated);
                REQUIRE(io::save_generated_mesh(p, *r.plan->generated));
                return p;
            }();
            std::string t = app_test::generate2d;
            const auto from = t.find("type = generate2d");
            const auto to = t.find("\n\n", from);
            t.replace(from, to - from, "type = file\npath = \"" + file.generic_string() + "\"");
            // A mesh file carries its regions: no polygon.
            const auto polygon = t.find("polygon = ");
            t.erase(polygon, t.find('\n', polygon) + 1 - polygon);
            return {t, "mesh"};
        }
        return {diode, "mesh"};
    }
    if (kind == "contact") {
        if (context == "gate") return {read_file(std::filesystem::path(NITCAD_EXAMPLES_DIR) / "moscap_cv.deck"), "contact gate"};
        if (context == "electrode") return {app_test::mosfet2d, "contact gate"};
        if (context == "schottky") {
            // The diode with a Schottky anode (Unit 48).
            std::string t = diode;
            const std::string ohmic = "[contact anode]\nkind = ohmic\n";
            t.replace(t.find(ohmic), ohmic.size(), "[contact anode]\nkind = schottky\nbarrier = 0.7\n");
            return {t, "contact anode"};
        }
        return {diode, "contact anode"};
    }
    if (kind == "region") return {diode, "region silicon"};
    if (kind == "interface") return {app_test::mosfet2d + "\n[interface silicon oxide]\n", "interface silicon oxide"};
    if (kind == "thermal_contact") {
        return {diode_with("\n[models]\nelectrothermal = on\n\n[thermal_contact sink]\npatch = x_min\nkind = isothermal\n"
                           "temperature = 300\n"),
                "thermal_contact sink"};
    }
    if (kind == "transient") return {app_test::transient_diode(), "transient"};
    if (kind == "small_signal") {
        return {app_test::diode_device + "\n[small_signal]\ncontact anode = 0\nfrequencies = 1e6\n", "small_signal"};
    }
    if (kind == "trace") return {app_test::diode_device + "\n[trace]\ncontact = anode\nend = 0.5\n", "trace"};
    if (kind == "equilibrium") return {app_test::diode_device + "\n[equilibrium]\n", "equilibrium"};
    return {diode, std::string(kind)};
}

// The text with `entry` added as the last entry of the section `header`, and its line.
struct Edited {
    std::string text;
    int line = 0;
    int header_line = 0;
};

Edited add_entry(const std::string& text, const std::string& header, const std::string& entry) {
    const auto found = app::find_sections(text);
    for (const auto& s : found) {
        const std::string h = s.name.empty() ? s.kind : s.kind + " " + s.name;
        if (h != header) continue;
        std::vector<std::string> lines;
        std::size_t start = 0;
        while (start <= text.size()) {
            const auto end = text.find('\n', start);
            lines.push_back(text.substr(start, end == std::string::npos ? std::string::npos : end - start));
            if (end == std::string::npos) break;
            start = end + 1;
        }
        lines.insert(lines.begin() + s.last_line, entry);
        std::string out;
        for (std::size_t k = 0; k < lines.size(); ++k) out += lines[k] + (k + 1 < lines.size() ? "\n" : "");
        return {out, s.last_line + 1, s.line};
    }
    // The section added at the end.
    std::string out = text;
    if (!out.empty() && out.back() != '\n') out += '\n';
    out += "\n[" + header + "]\n";
    const int header_line = static_cast<int>(std::ranges::count(out, '\n'));
    out += entry + "\n";
    return {out, header_line + 1, header_line};
}

std::map<std::string, std::string> samples() {
    return {{"x", "0, 1um"},
            {"y", "0, 1um"},
            {"z", "0, 1um"},
            {"straddle", "1um 10nm"},
            {"max_size", "0.3um"},
            {"interface_resolution", "0.1um"},
            {"vertex", "0 0 0"},
            {"face", "0 1 2"},
            {"plane", "0 0 0.5um, 1um 0 0.5um, 0 1um 0.5um"},
            {"path", "\"missing.nitmesh\""},
            {"polygon", "0 0, 1um 0, 1um 1um"},
            {"seed", "0 0 0"},
            {"axis", "x"},
            {"stations", "0, 1um"},
            {"cell", "0 0, 1um 0, 1um 1um"},
            {"extrude", "c 0:1um"},
            {"profile", "donor constant 1e15"},
            {"oxide_thickness", "10nm"},
            {"work_function", "4.1"},
            {"trap_level", "acceptor 1e10 0.1"},
            {"trap_band", "acceptor 1e10 -0.2 0.2"},
            {"thermal_velocity_n", "1e7"},
            {"thermal_velocity_p", "1e7"},
            {"resistance", "1"},
            {"mkl", "\"missing_mkl_rt.dll\""},
            {"Ns", "1e17"},
            {"restart", "\"missing.nitchk\""},
            {"t_end", "1e-9"},
            {"dt_initial", "1e-12"},
            {"dt_min", "1e-15"},
            {"dt_max", "1e-10"},
            {"output_times", "1e-9"},
            {"frequencies", "1e6"},
            {"field_frequencies", "1e6"},
            {"end", "0.5"},
            {"step", "0.01"},
            {"max_step", "0.1"},
            {"min_step", "1e-6"},
            {"current_limit", "1"},
            {"max_points", "10"},
            {"density_floor", "1e10"},
            {"sweep", "x.nitsweep"},
            {"device", "x.nitdev"},
            {"mesh", "x.nitmesh"},
            {"checkpoint", "x.nitchk"},
            {"table", "x.csv"},
            {"transient", "x.nittran"},
            {"small_signal", "x.nitac"}};
}

std::string sample_value(const app::KeyInfo& k, std::string_view kind) {
    const auto s = samples();
    if (kind == "transient" && k.key == "contact") return "constant(0)";
    if (kind == "trace" && k.key == "contact") return "anode";
    if (k.key == "box") return kind == "mesh" ? "0 0 1um 1um 0.1um" : kind == "region" ? "-1um 3um" : "0 1um";
    if (k.key == "patch") return kind == "mesh" ? "0 0, 1um 0" : "x_min";
    if (const auto it = s.find(std::string(k.key)); it != s.end()) return it->second;
    if (k.key == "contact" || k.key == "thermal" || k.key == "start") return "0";
    return std::string(k.value);
}

std::string name_for(const app::KeyInfo& k) {
    if (k.names == 2) return " x_min a";
    if (k.names == 0) return {};
    if (k.key == "contact" || k.key == "start") return " cathode";
    if (k.key == "thermal") return " sink";
    if (k.key == "straddle") return " x";
    return " extra";
}

// The contexts of a section: its keys' contexts, or none.
std::vector<std::string> contexts_of(const app::SectionInfo& s) {
    std::set<std::string> out;
    for (const auto& k : s.keys) {
        std::string c(k.contexts);
        std::size_t i = 0;
        while (i < c.size()) {
            const auto j = c.find(' ', i);
            out.insert(c.substr(i, j == std::string::npos ? std::string::npos : j - i));
            if (j == std::string::npos) break;
            i = j + 1;
        }
    }
    if (s.kind == "mesh") out.insert("tensor");
    if (s.kind == "contact") out.insert("ohmic");
    if (out.empty()) out.insert("");
    return {out.begin(), out.end()};
}

// Whether deck::build takes the entry in the section: the deck builds, or its error is about the
// entry's line or section and is not "unknown key".
bool accepted(const Base& base, const std::string& entry, std::string* why = nullptr) {
    const Edited e = add_entry(base.text, base.header, entry);
    const auto r = app::check_deck(e.text, {});
    if (r.plan) return true;
    const auto& d = r.diagnostics.front();
    if (why) *why = d.message;
    if (d.message.find("unknown key") != std::string::npos) return false;
    return d.line == e.line || d.line == e.header_line;
}

}  // namespace

TEST_CASE("app vocabulary: every key is one deck::build takes (XG2)", "[app]") {
    int probed = 0;
    for (const app::SectionInfo& s : app::deck_vocabulary()) {
        for (const std::string& context : contexts_of(s)) {
            const Base base = base_for(s.kind, context);
            {
                const auto r = app::check_deck(base.text, {});
                INFO(s.kind << " " << context << ": " << (r.plan ? "" : r.diagnostics.front().message));
                REQUIRE(r.plan);
            }
            for (const app::KeyInfo* k : app::section_keys(s, context)) {
                ++probed;
                // The key that selects the context (type, kind) is in the base, which builds: a second
                // one would only switch the context.
                if (k->key == s.context_key) continue;
                const std::string entry = std::string(k->key) + name_for(*k) + " = " + sample_value(*k, s.kind);
                std::string why;
                INFO("[" << s.kind << "] (" << context << ") " << entry);
                const bool ok = accepted(base, entry, &why);
                INFO(why);
                CHECK(ok);
            }
        }
    }
    CHECK(probed > 130);
}

TEST_CASE("app vocabulary: every key deck/build.hpp names that deck::build takes is listed (XG2)", "[app]") {
    const std::string header = read_file(std::filesystem::path(NITCAD_SOURCE_DIR) / "include/NiTCAD/deck/build.hpp");
    const auto from = header.find("// Sections");
    const auto to = header.find("// Exactly one analysis section");
    REQUIRE(from != std::string::npos);
    REQUIRE(to != std::string::npos);
    // The words of each section's part of the reference.
    std::map<std::string, std::set<std::string>> words;
    std::string kind;
    std::size_t i = from;
    while (i < to) {
        const auto end = header.find('\n', i);
        const std::string line = header.substr(i, end - i);
        i = end + 1;
        if (line.rfind("//   [", 0) == 0) kind = line.substr(6, line.find_first_of(" ]", 6) - 6);
        if (kind.empty()) continue;
        std::size_t j = 2;
        while (j < line.size()) {
            if (!std::isalpha(static_cast<unsigned char>(line[j]))) {
                ++j;
                continue;
            }
            std::size_t e = j;
            while (e < line.size() && (std::isalnum(static_cast<unsigned char>(line[e])) || line[e] == '_' || line[e] == '.')) {
                ++e;
            }
            std::string w = line.substr(j, e - j);
            while (!w.empty() && w.back() == '.') w.pop_back();
            words[kind].insert(w);
            j = e;
        }
    }
    REQUIRE(words.size() == app::deck_vocabulary().size());
    // A region's material parameters by member path are open-ended and not completed (XD3).
    const std::set<std::string> parameters{"eps_r", "electron_mobility.mu_max", "thermal.conductivity_W_cmK",
                                           "passler_gap.Theta_K", "energy_relaxation.tau_n_s"};
    int taken = 0;
    for (const auto& [section, list] : words) {
        const app::SectionInfo* s = app::find_section_info(section);
        REQUIRE(s != nullptr);
        for (const std::string& w : list) {
            if (parameters.contains(w)) continue;
            bool listed = false;
            for (const auto& k : s->keys) listed = listed || k.key == w;
            if (listed) {
                ++taken;
                continue;
            }
            for (const std::string& context : contexts_of(*s)) {
                const Base base = base_for(section, context);
                for (const std::string& entry : {w + " = 1", w + " cathode = 1"}) {
                    INFO("[" << section << "] (" << context << ") " << entry << " is taken but not listed");
                    CHECK_FALSE(accepted(base, entry));
                }
            }
        }
    }
    CHECK(taken > 100);
}

TEST_CASE("app vocabulary: the menus' defaults are the vocabulary's", "[app]") {
    for (const app::DeckMenu m : {app::DeckMenu::device, app::DeckMenu::mesh, app::DeckMenu::physics,
                                  app::DeckMenu::simulation}) {
        for (const app::SectionSpec& s : app::menu_sections(m)) {
            const app::SectionInfo* info = app::find_section_info(s.kind);
            REQUIRE(info != nullptr);
            for (const app::KeySpec& k : s.keys) {
                const app::KeyInfo* v = app::find_key_info(*info, k.key, {});
                REQUIRE(v != nullptr);
                CHECK(v->value == k.value);
            }
        }
    }
}

TEST_CASE("app vocabulary: edit distance and the nearest name", "[app]") {
    CHECK(app::edit_distance("", "abc") == 3);
    CHECK(app::edit_distance("kitten", "sitting") == 3);
    CHECK(app::edit_distance("srh", "srh") == 0);
    CHECK(app::edit_distance("max_iteration", "max_iterations") == 1);
    const std::vector<std::string_view> names{"tol_update", "max_update", "max_iterations"};
    CHECK(app::nearest("max_updte", names) == "max_update");
    CHECK(app::nearest("max_update", names) == std::nullopt);       // exact
    CHECK(app::nearest("min_upd", names) == std::nullopt);          // 3 edits
    CHECK(app::nearest("max_iteratn", names) == std::nullopt);      // 3 edits: never more than 2
    CHECK(app::nearest("max_itrations", names) == "max_iterations");
    // At most a third of the word: one edit for words of up to five letters.
    const std::vector<std::string_view> short_names{"srh", "bgn"};
    CHECK(app::nearest("sr", short_names) == "srh");
    CHECK(app::nearest("s", short_names) == std::nullopt);
}

namespace {

// The word with its middle letter dropped: one edit.
std::string misspelt(std::string_view w) {
    std::string s(w);
    s.erase(s.size() / 2, 1);
    return s;
}

}  // namespace

TEST_CASE("app vocabulary: misspelt sections, keys and choices are suggested back (XG2)", "[app]") {
    int suggested = 0;
    // Sections: an unknown header.
    for (const app::SectionInfo& s : app::deck_vocabulary()) {
        if (s.kind.size() < 4) continue;
        const std::string text = app_test::diode_sweep() + "\n[" + misspelt(s.kind) + (s.names ? " a b" : "") + "]\n";
        const auto r = app::check_deck(text, {});
        REQUIRE_FALSE(r.plan);
        INFO(s.kind << ": " << r.diagnostics.front().message);
        CHECK(app::suggestion(text, r.diagnostics.front()) == "did you mean '" + std::string(s.kind) + "'?");
        ++suggested;
    }
    // Keys and their choices, each in a deck of its section and context.
    for (const app::SectionInfo& s : app::deck_vocabulary()) {
        for (const std::string& context : contexts_of(s)) {
            const Base base = base_for(s.kind, context);
            for (const app::KeyInfo* k : app::section_keys(s, context)) {
                if (k->key.size() >= 4) {
                    const std::string entry = misspelt(k->key) + name_for(*k) + " = " + sample_value(*k, s.kind);
                    const Edited e = add_entry(base.text, base.header, entry);
                    const auto r = app::check_deck(e.text, {});
                    INFO("[" << s.kind << "] (" << context << ") " << entry);
                    REQUIRE_FALSE(r.plan);
                    INFO(r.diagnostics.front().message);
                    CHECK(r.diagnostics.front().line == e.line);
                    CHECK(app::suggestion(e.text, r.diagnostics.front()) == "did you mean '" + std::string(k->key) + "'?");
                    ++suggested;
                }
                for (int position = 0; position < 2; ++position) {
                    const auto& choices = position == 0 ? k->choices : k->choices2;
                    for (const std::string_view choice : choices) {
                        if (choice.size() < 4) continue;
                        std::string value = sample_value(*k, s.kind);
                        auto words = std::vector<std::string>{};
                        std::size_t at = 0;
                        while (at < value.size()) {
                            const auto sp = value.find(' ', at);
                            words.push_back(value.substr(at, sp == std::string::npos ? std::string::npos : sp - at));
                            if (sp == std::string::npos) break;
                            at = sp + 1;
                        }
                        if (static_cast<std::size_t>(position) >= words.size()) continue;
                        std::string bad = misspelt(choice);
                        if (choice == "algaas") bad += " 0.3";
                        const auto paren = words[static_cast<std::size_t>(position)].find('(');
                        words[static_cast<std::size_t>(position)] =
                            bad + (paren == std::string::npos ? "" : words[static_cast<std::size_t>(position)].substr(paren));
                        value.clear();
                        for (std::size_t w = 0; w < words.size(); ++w) value += (w ? " " : "") + words[w];
                        const std::string entry = std::string(k->key) + name_for(*k) + " = " + value;
                        const Edited e = add_entry(base.text, base.header, entry);
                        const auto r = app::check_deck(e.text, {});
                        INFO("[" << s.kind << "] (" << context << ") " << entry);
                        REQUIRE_FALSE(r.plan);
                        INFO(r.diagnostics.front().message);
                        const app::Diagnostic& d = r.diagnostics.front();
                        CHECK(d.line == e.line);
                        CHECK(app::suggestion(e.text, d) == "did you mean '" + std::string(choice) + "'?");
                        ++suggested;
                    }
                }
            }
        }
    }
    CHECK(suggested > 150);
    // Nothing for an error that is not a misspelling.
    const std::string text = app_test::diode_sweep() + "\n[models]\nsrh = maybe\n";
    const auto r = app::check_deck(text, {});
    REQUIRE_FALSE(r.plan);
    CHECK(app::suggestion(text, r.diagnostics.front()) == std::nullopt);
}

namespace {

std::vector<std::string> labels(const std::optional<app::Completion>& c) {
    std::vector<std::string> out;
    if (c) {
        for (const auto& i : c->items) out.push_back(i.label);
    }
    return out;
}

}  // namespace

TEST_CASE("app vocabulary: completion at the caret (XD3, XG2)", "[app]") {
    using V = std::vector<std::string>;
    // Section kinds after '['.
    {
        const std::string text = "[mesh]\ntype = tensor\nx = 0, 1um\n[s";
        const auto c = app::complete(text, 4, 2);
        REQUIRE(c);
        CHECK(labels(c) == V{"solver", "sweep", "small_signal"});
        CHECK(c->from == 1);
        CHECK(c->to == 2);
        CHECK(c->items[0].insert == "solver]");
        const auto all = app::complete("[", 1, 1);
        CHECK(all->items.size() == app::deck_vocabulary().size());
        CHECK(app::complete("[con", 1, 4)->items[0].insert == "contact ");      // named: the name follows
        CHECK(app::complete("[mo]", 1, 3)->items[0].insert == "models");        // already closed
        CHECK(labels(app::complete("[region si", 1, 10)).empty());             // a name: nothing
    }
    // A section's keys at an entry's start: those of its context, present ones not again.
    {
        const std::string text = "[mesh]\ntype = generate2d\nmax_size = 1um\nm\n";
        const auto c = app::complete(text, 4, 1);
        REQUIRE(c);
        CHECK(labels(c) == V{"min_angle", "max_nodes"});
        CHECK(c->items[0].insert == "min_angle = ");
        CHECK(labels(app::complete("[mesh]\ntype = tensor\n\n", 3, 0)) == V{"x", "y", "z", "straddle"});
        // Repeatable keys stay offered; a named key inserts its word and a blank.
        const std::string sweep = "[sweep]\ncontact anode = 0\nc";
        const auto s = app::complete(sweep, 3, 1);
        CHECK(labels(s) == V{"contact"});
        CHECK(s->items[0].insert == "contact ");
        // Contact keys by kind.
        CHECK(labels(app::complete("[contact a]\nkind = ohmic\n\n", 3, 0)) == V{"patch", "box"});
        CHECK(labels(app::complete("[contact a]\nkind = gate\no\n", 3, 1)) ==
              V{"oxide_thickness", "oxide_permittivity", "oxide_affinity", "oxide_band_gap",
                "oxide_mass_n", "oxide_mass_p", "oxide_mass_evb"});  // Unit 51 (GD10) added five
        // The word at the caret is replaced whole.
        const auto mid = app::complete("[solver]\nmax_iter = 5\n", 2, 4);
        REQUIRE(mid);
        CHECK(mid->from == 0);
        CHECK(mid->to == 8);
        CHECK(labels(mid) == V{"max_iterations", "max_update", "max_cuts", "max_backward_error",
                               "max_refinement_steps", "max_linear_iterations"});
    }
    // Value choices after '='.
    {
        CHECK(labels(app::complete("[models]\nsrh = \n", 2, 6)) == V{"on", "off"});
        CHECK(labels(app::complete("[models]\nhole_mobility = ma\n", 2, 18)) == V{"masetti"});
        CHECK(labels(app::complete("[doping]\nprofile = donor g\n", 2, 17)) == V{"gaussian"});
        CHECK(labels(app::complete("[doping]\nprofile = \n", 2, 10)) == V{"donor", "acceptor"});
        const auto w = app::complete("[transient]\ncontact a = st\n", 2, 14);
        REQUIRE(w);
        CHECK(labels(w) == V{"step"});
        CHECK(w->items[0].insert == "step(");
        CHECK(labels(app::complete("[solver]\nmax_iterations = 1\n", 2, 18)).empty());  // no choices
    }
    // Nothing in a comment, a continued value, a header's line or before the first section.
    CHECK_FALSE(app::complete("[models]\n# srh\n", 2, 5));
    CHECK_FALSE(app::complete("[models]\nsrh = #o\n", 2, 8));  // a value's choices would fit "o"
    CHECK_FALSE(app::complete("[region a]\npolygon = 0 0, \\\n  m", 3, 3));
    CHECK_FALSE(app::complete("ty\n[mesh]\n", 1, 2));
    CHECK_FALSE(app::complete("[mesh]  \n", 1, 8));
    CHECK_FALSE(app::complete("[mesh]\n", 5, 0));
}

TEST_CASE("app vocabulary: the Schottky kind in the deck, the app and the CLI (Unit 48, SC15)", "[app]") {
    const app::SectionInfo* contact = app::find_section_info("contact");
    REQUIRE(contact != nullptr);
    const app::KeyInfo* kind = app::find_key_info(*contact, "kind", "");
    REQUIRE(kind != nullptr);
    CHECK(std::ranges::find(kind->choices, std::string_view{"schottky"}) != kind->choices.end());
    for (const std::string_view key : {"barrier", "richardson_n", "richardson_p", "work_function"}) {
        CAPTURE(key);
        CHECK(app::find_key_info(*contact, key, "schottky") != nullptr);
        if (key != "work_function") CHECK(app::find_key_info(*contact, key, "ohmic") == nullptr);
    }
    std::string text = app_test::diode_sweep();
    const std::string ohmic = "[contact anode]\nkind = ohmic\n";
    text.replace(text.find(ohmic), ohmic.size(), "[contact anode]\nkind = schottky\nbarrier = 0.7\n");
    const auto r = app::check_deck(text, {});
    REQUIRE(r.plan);
    CHECK(r.plan->device.contacts()[0].kind == device::ContactKind::schottky);
    const auto deck = app_test::write(app_test::scratch("vocabulary_schottky"), "sbd.deck", text);
    const auto cli = app_test::nitcad("check " + app_test::quoted(deck));
    CAPTURE(cli.output);
    REQUIRE(cli.status == 0);
    CHECK(cli.output.find("contact anode: schottky, 1 nodes") != std::string::npos);
}

TEST_CASE("app vocabulary: a mesh patch's form follows the mesh type (Unit 53, FD13)", "[app]") {
    const app::SectionInfo* mesh = app::find_section_info("mesh");
    REQUIRE(mesh != nullptr);
    const app::KeyInfo* extruded = app::find_key_info(*mesh, "patch far", "extrude3d");
    REQUIRE(extruded != nullptr);
    CHECK(extruded->value == "<x0 y0 z0 x1 y1 z1>");
    const app::KeyInfo* planar = app::find_key_info(*mesh, "patch far", "generate2d");
    REQUIRE(planar != nullptr);
    CHECK(planar->value == "<x y>, <x y>");
    const app::KeyInfo* faced = app::find_key_info(*mesh, "patch x_min far", "generate3d");
    REQUIRE(faced != nullptr);
    CHECK(faced->value == "<x y z>, ...");
}
