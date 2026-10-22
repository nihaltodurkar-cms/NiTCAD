// The menus' deck sections (Unit 26, VQ13, WG7): sections found in the text, the keys offered, and
// every offered key inserted with its default giving a deck that builds with the same meaning.
#include <catch2/catch_test_macros.hpp>

#include <format>
#include <string>
#include <vector>

#include "NiTCAD/app/deck_sections.hpp"
#include "NiTCAD/app/document.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;

namespace {

// What a plan means for the keys the menus insert: the device's temperature and mesh, and every
// option those keys set. Two plans with equal fingerprints solve the same problem the same way.
std::string fingerprint(const deck::Plan& p) {
    const auto& m = p.options.models;
    const auto& n = p.options.newton;
    const auto& s = p.options.step_control;
    const auto& l = p.options.linear;
    const auto& t = p.transient;
    return std::format(
        "T {} mesh {:x} | models {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} | newton {} {} {} {} {} {} {} | "
        "steps {} {} {} | linear {} {} {} {} {} {} {} {} {} {} {} | eq {} w {} r {} | tr {} {} {} {} {}",
        p.device.temperature_K(), mesh::graph_digest(p.device.mesh()), m.doping_mobility,
        static_cast<int>(m.hole_mobility), static_cast<int>(m.electron_mobility), m.srh, m.auger, m.bgn,
        m.field_mobility, m.fermi_dirac, m.radiative, m.incomplete_ionization, static_cast<int>(m.ionization_levels),
        m.impact_ionization, m.impact_current_resolution, m.btbt_local, static_cast<int>(m.btbt_nonlocal),
        m.electrothermal, n.max_iterations, n.tol_update, n.max_update, n.line_search.enabled,
        n.line_search.max_backtracks, n.line_search.sufficient_decrease, static_cast<int>(n.line_search.test),
        s.max_cuts, s.growth, s.max_iterations, static_cast<int>(l.backend), l.threads, l.equilibrate,
        l.max_backward_error, l.max_refinement_steps, l.min_pivot_ratio, l.restart, l.max_iterations,
        l.ilut_drop_tolerance, l.ilut_fill_factor, static_cast<int>(l.fallback), static_cast<int>(p.options.equations),
        p.options.weighted_current, p.options.store_rates, static_cast<int>(t.integrator), t.adaptive, t.rtol,
        t.density_ref_cm3, t.fields_every_step);
}

deck::Plan built(const std::string& text) {
    auto r = app::check_deck(text, {});
    if (!r.plan) FAIL(r.diagnostics.front().message);
    return std::move(*r.plan);
}

using app_test::generate2d;
using app_test::generate3d;
using app_test::transient_diode;

}  // namespace

TEST_CASE("app deck sections: headers, names, entries and their lines", "[app]") {
    const std::string text = "# a comment\n[mesh]\ntype = tensor   # trailing\nx = 0, 1um, \\\n    2um\n\n"
                             "[region  silicon]\nmaterial = silicon\n[interface silicon   oxide]\n"
                             "[contact anode]\nkind = ohmic\npatch = \"x_min # not a comment\"\n";
    const auto s = app::find_sections(text);
    REQUIRE(s.size() == 4);
    CHECK(s[0].kind == "mesh");
    CHECK(s[0].line == 2);
    CHECK(s[0].last_line == 5);  // the continued value's second line
    REQUIRE(s[0].entries.size() == 2);
    CHECK(s[0].entries[0].key == "type");
    CHECK(s[0].entries[0].value == "tensor");
    CHECK(s[0].entries[1].line == 4);
    CHECK(s[1].kind == "region");
    CHECK(s[1].name == "silicon");
    CHECK(s[2].kind == "interface");
    CHECK(s[2].name == "silicon oxide");
    CHECK(s[2].line == s[2].last_line);
    CHECK(s[3].entries[1].value == "\"x_min # not a comment\"");
    CHECK(app::find_section(text, "contact")->name == "anode");
    CHECK_FALSE(app::find_section(text, "solver"));
}

TEST_CASE("app deck sections: the menus and the keys they offer", "[app]") {
    CHECK(app::to_string(app::DeckMenu::physics) == "Physics");
    CHECK(app::menu_sections(app::DeckMenu::device).front().kind == "device");
    CHECK(app::menu_sections(app::DeckMenu::physics).size() == 1);
    CHECK(app::find_section_spec("solver")->keys.size() == 23);
    CHECK(app::find_section_spec("unknown") == nullptr);
    const auto keys_of = [](const std::string& text, const char* kind) {
        std::vector<std::string> out;
        for (const auto& k : app::insertable_keys(*app::find_section(text, kind))) out.emplace_back(k.key);
        return out;
    };
    using K = std::vector<std::string>;
    CHECK(keys_of(app_test::diode_device, "mesh").empty());  // a tensor mesh: none
    CHECK(keys_of(generate2d, "mesh") == K{"min_angle", "max_nodes"});
    CHECK(keys_of(generate3d, "mesh") == K{"radius_edge", "max_nodes"});
    // A key the section has is not offered again.
    const std::string with = app_test::diode_device + "[solver]\nmax_iterations = 50\n";
    const auto solver = keys_of(with, "solver");
    CHECK(solver.size() == 22);
    CHECK(std::ranges::find(solver, "max_iterations") == solver.end());
}

TEST_CASE("app deck sections: inserting a key", "[app]") {
    const app::KeySpec key{"max_iterations", "100", ""};
    // Into an existing section, after its last entry, before the next section.
    const std::string text = "[solver]\ntol_update = 1e-8\n\n[sweep]\ncontact a = 0\n";
    const auto in = app::insert_key(text, "solver", key);
    REQUIRE(in);
    CHECK(in->text == "[solver]\ntol_update = 1e-8\nmax_iterations = 100\n\n[sweep]\ncontact a = 0\n");
    CHECK(in->line == 3);
    // A section without entries; a text without a final newline.
    const auto bare = app::insert_key("[solver]", "solver", key);
    REQUIRE(bare);
    CHECK(bare->text == "[solver]\nmax_iterations = 100\n");
    CHECK(bare->line == 2);
    // A missing section is added at the end.
    const auto added = app::insert_key("[sweep]\ncontact a = 0", "solver", key);
    REQUIRE(added);
    CHECK(added->text == "[sweep]\ncontact a = 0\n\n[solver]\nmax_iterations = 100\n");
    CHECK(added->line == 5);
    // Not a second analysis.
    CHECK_FALSE(app::insert_key("[sweep]\ncontact a = 0\n", "transient", {"rtol", "1e-3", ""}));
}

TEST_CASE("app deck sections: every offered key keeps the deck's meaning (WG7)", "[app]") {
    for (const std::string& text : {app_test::diode_sweep(), app_test::mosfet2d, generate2d, generate3d,
                                    transient_diode()}) {
        CAPTURE(text.substr(0, text.find('\n')));
        const std::string reference = fingerprint(built(text));
        int inserted = 0;
        for (const app::DeckMenu menu : {app::DeckMenu::device, app::DeckMenu::mesh, app::DeckMenu::physics,
                                         app::DeckMenu::simulation}) {
            for (const app::SectionSpec& spec : app::menu_sections(menu)) {
                // The section as it is, or as the insertion will create it.
                const auto here = app::find_section(text, spec.kind);
                app::SectionLocation location;
                location.kind = std::string(spec.kind);
                const auto keys = app::insertable_keys(here ? *here : location);
                for (const app::KeySpec& key : keys) {
                    CAPTURE(spec.kind, key.key);
                    const auto in = app::insert_key(text, spec.kind, key);
                    if (!in) {
                        CHECK_FALSE(here);  // only an absent analysis section is refused
                        continue;
                    }
                    const auto r = app::check_deck(in->text, {});
                    INFO((r.diagnostics.empty() ? std::string() : r.diagnostics.front().message));
                    REQUIRE(r.plan);
                    CHECK(fingerprint(*r.plan) == reference);
                    // The inserted line holds the key.
                    const auto sections = app::find_sections(in->text);
                    bool found = false;
                    for (const auto& s : sections) {
                        for (const auto& e : s.entries) found = found || (e.line == in->line && e.key == key.key);
                    }
                    CHECK(found);
                    ++inserted;
                }
            }
        }
        CHECK(inserted >= 40);
    }
}
