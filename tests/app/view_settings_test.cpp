// View settings (Unit 26, VA4, VA8, WG3): the views and fields offered, the errors, and the scenes,
// drawn on WARP, equal to `nitcad plot`'s images for the same options (WG4's choices).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/render/image.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
using app::ViewKind;
using app::ViewSettings;

namespace {

struct Ran {
    fs::path deck;
    app::Workspace workspace;
};

// The deck run once per test process, in its own directory.
const Ran& ran(const std::string& name, const std::string& text) {
    static std::vector<std::unique_ptr<Ran>> cache;
    for (const auto& r : cache) {
        if (r->deck.parent_path().filename() == name) return *r;
    }
    const fs::path dir = app_test::scratch(name);
    const fs::path deck = app_test::write(dir, "d.deck", text);
    app::RunController rc;
    REQUIRE(rc.start_run(text, deck));
    rc.wait();
    rc.process();
    REQUIRE(rc.state() == app::RunState::done);
    auto taken = rc.take_results();
    auto w = app::workspace_from_run(std::move(taken->plan), std::move(taken->report));
    REQUIRE(w);
    cache.push_back(std::make_unique<Ran>(Ran{deck, std::move(*w)}));
    return *cache.back();
}

const Ran& diode() { return ran("views_diode", app_test::diode_sweep()); }
const Ran& mosfet() { return ran("views_mosfet", app_test::mosfet2d); }

render::Renderer& warp() {
    static auto r = render::Renderer::create({.adapter = render::Adapter::warp});
    REQUIRE(r);
    return *r;
}

render::Image draw(const app::Scene& s, int w, int h) {
    auto image = std::visit([&](const auto& scene) { return warp().draw(scene, w, h); }, s);
    REQUIRE(image);
    return std::move(*image);
}

bool has(const std::vector<std::string>& v, const std::string& x) { return std::ranges::find(v, x) != v.end(); }

}  // namespace

TEST_CASE("app views: names, and the views and fields a workspace offers", "[app]") {
    for (const ViewKind v : {ViewKind::solid, ViewKind::mesh, ViewKind::field, ViewKind::slice, ViewKind::iso,
                             ViewKind::bands, ViewKind::profile, ViewKind::iv, ViewKind::cv, ViewKind::current}) {
        CHECK(app::find_view(app::to_string(v)) == v);
    }
    CHECK_FALSE(app::find_view("exploded"));
    CHECK(app::is_plot(ViewKind::iv));
    CHECK_FALSE(app::is_plot(ViewKind::field));

    using V = std::vector<ViewKind>;
    CHECK(app::available_views(diode().workspace) == V{ViewKind::profile, ViewKind::bands, ViewKind::iv});
    CHECK(app::available_views(mosfet().workspace) ==
          V{ViewKind::solid, ViewKind::mesh, ViewKind::field, ViewKind::iso, ViewKind::iv});
    const auto fields = app::available_fields(mosfet().workspace);
    for (const char* f : {"potential", "n", "p", "Ec", "doping", "R", "srh", "E", "J"}) CHECK(has(fields, f));
    CHECK_FALSE(has(fields, "Gii"));  // impact ionization is off: no rate stored for it
    // Before a run: the doping only.
    const fs::path dir = app_test::scratch("views_unrun");
    auto unrun = app::workspace_from_deck(app_test::write(dir, "d.deck", app_test::mosfet2d));
    REQUIRE(unrun);
    CHECK(app::available_fields(*unrun) == std::vector<std::string>{"doping", "ND", "NA"});
    CHECK(app::available_views(*unrun) == V{ViewKind::solid, ViewKind::mesh, ViewKind::field});
}

TEST_CASE("app views: errors as nitcad plot reports them", "[app]") {
    const app::Workspace& w = mosfet().workspace;
    ViewSettings s;
    s.kind = ViewKind::field;
    s.field = "nonsense";
    CHECK(app::build_scene(w, s).error().message == "unknown field nonsense");
    s.field = "n";
    s.scale = "cubic";
    CHECK(app::build_scene(w, s).error().message == "unknown scale cubic");
    s.scale.clear();
    s.colormap = "rainbow";
    CHECK(app::build_scene(w, s).error().message == "unknown colormap rainbow");
    s.colormap = "turbo";
    s.point = 99;
    CHECK(app::build_scene(w, s).error().message.find("out of range") != std::string::npos);
    s.point.reset();
    s.kind = ViewKind::slice;
    CHECK_FALSE(app::build_scene(w, s));  // no slice plane given
    s.kind = ViewKind::field;
    s.clip = ViewSettings::Clip{{0, 0, 0}, {1, 1, 1}};
    CHECK(app::build_scene(w, s).error().message == "a clipping box needs a 3D device");
    s.clip.reset();
    s.kind = ViewKind::cv;
    CHECK(app::build_scene(w, s).error().message == "a C-V plot needs a small-signal run's results");
    s.kind = ViewKind::field;
    CHECK(app::build_scene(w, s).has_value());
    CHECK(std::holds_alternative<render::FieldScene>(*app::build_scene(w, s)));
    s.kind = ViewKind::iv;
    CHECK(std::holds_alternative<render::PlotScene>(*app::build_scene(w, s)));
}

TEST_CASE("app views: the options of nitcad plot for a setting", "[app]") {
    ViewSettings s;
    using A = std::vector<std::string>;
    CHECK(app::plot_arguments(s) == A{"--view", "solid"});
    s.kind = ViewKind::field;
    s.field = "n";
    s.point = 2;
    s.scale = "log";
    s.colormap = "viridis";
    s.range = std::pair{1e10, 1e20};
    s.contours = s.edges = s.opaque = s.perspective = true;
    s.azimuth = 30.5f;
    s.elevation = -10.0f;
    CHECK(app::plot_arguments(s) == A{"--view", "field", "--field", "n", "--point", "2", "--scale", "log",
                                      "--colormap", "viridis", "--range", "1e+10,1e+20",
                                      "--contours", "--edges", "--opaque", "--perspective", "--azimuth", "30.5",
                                      "--elevation", "-10"});
    ViewSettings t;
    t.kind = ViewKind::slice;
    t.slices = {{0, 0.25}, {2, 1.5}};
    t.levels = {1e17, 2.5e18};
    t.clip = ViewSettings::Clip{{0, 0, 0}, {1, 0.5, 2}};
    t.contact = "gate";
    t.sense = "drain";
    t.frequency = 1;
    t.log = true;
    CHECK(app::plot_arguments(t) == A{"--view", "slice", "--field", "potential", "--slice", "x=0.25", "--slice",
                                      "z=1.5", "--clip", "0,0,0:1,0.5,2", "--levels", "1e+17,2.5e+18",
                                      "--contact", "gate", "--sense", "drain", "--frequency", "1", "--log"});
}

TEST_CASE("app views: the same images as nitcad plot (WG4's choices)", "[app]") {
    const fs::path out = app_test::scratch("views_images");
    struct Case {
        const Ran* ran;
        ViewSettings s;
        const char* name;
    };
    const auto make = [](ViewKind k, std::string field = "potential") {
        ViewSettings s;
        s.kind = k;
        s.field = std::move(field);
        return s;
    };
    std::vector<Case> cases;
    cases.push_back({&mosfet(), make(ViewKind::solid), "solid"});
    cases.push_back({&mosfet(), make(ViewKind::mesh), "mesh"});
    {
        ViewSettings s = make(ViewKind::field);
        s.contours = true;
        cases.push_back({&mosfet(), s, "potential"});
    }
    {
        ViewSettings s = make(ViewKind::field, "n");
        s.point = 1;
        s.edges = true;
        s.colormap = "viridis";
        cases.push_back({&mosfet(), s, "n"});
    }
    cases.push_back({&mosfet(), make(ViewKind::field, "J"), "J"});
    {
        ViewSettings s = make(ViewKind::iso, "n");
        s.levels = {1e10, 1e16};
        cases.push_back({&mosfet(), s, "iso"});
    }
    {
        ViewSettings s = make(ViewKind::iv);
        s.log = true;
        s.sense = "drain";
        cases.push_back({&mosfet(), s, "iv"});
    }
    {
        ViewSettings s = make(ViewKind::profile, "n");
        s.log = true;
        cases.push_back({&diode(), s, "profile"});
    }
    cases.push_back({&diode(), make(ViewKind::bands), "bands"});
    cases.push_back({&diode(), make(ViewKind::iv), "diode_iv"});
    for (const Case& c : cases) {
        CAPTURE(c.name);
        auto scene = app::build_scene(c.ran->workspace, c.s);
        REQUIRE(scene);
        const render::Image mine = draw(*scene, 480, 360);
        const fs::path png = out / (std::string(c.name) + ".png");
        std::string args = "plot " + app_test::quoted(c.ran->deck);
        for (const std::string& a : app::plot_arguments(c.s)) args += " \"" + a + "\"";
        args += " --size 480x360 --adapter warp --out " + app_test::quoted(png);
        const auto cli = app_test::nitcad(args);
        CAPTURE(args, cli.output);
        REQUIRE(cli.status == 0);
        auto theirs = render::load_png(png);
        REQUIRE(theirs);
        const auto d = render::compare(mine, *theirs, 0);
        REQUIRE(d);
        CHECK(d->max_channel == 0);
    }
}
