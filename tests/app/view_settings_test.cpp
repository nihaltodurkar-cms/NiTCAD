// View settings (Unit 26, VA4, VA8, WG3): the views and fields offered, the errors, and the scenes,
// drawn on WARP, equal to `nitcad plot`'s images for the same options (WG4's choices).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <string>
#include <utility>
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
    w->source = deck;  // as the window sets it (run B is named by it)
    cache.push_back(std::make_unique<Ran>(Ran{deck, std::move(*w)}));
    return *cache.back();
}

const Ran& diode() { return ran("views_diode", app_test::diode_sweep()); }
const Ran& mosfet() { return ran("views_mosfet", app_test::mosfet2d); }
// Run B (Unit 26b, XD6): the transistor with a higher gate work function, on the same mesh (or, with
// another x, on a finer one).
std::string variant(const std::string& from, const std::string& to) {
    std::string t = app_test::mosfet2d;
    t.replace(t.find(from), from.size(), to);
    return t;
}
const Ran& mosfet_b() { return ran("views_mosfet_b", variant("work_function = 4.1", "work_function = 4.3")); }

// Unit 26c: the small 3D transistor, the diode's C-V at two frequencies and its transient.
const Ran& mosfet3d() { return ran("views_mosfet3d", app_test::mosfet3d); }
const Ran& diode_cv() { return ran("views_diode_cv", app_test::diode_cv()); }
const Ran& diode_transient() {
    return ran("views_diode_transient", app_test::transient_diode() + "\n[output]\ntransient = t.nittran\ndevice = t.nitdev\n");
}

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
                             ViewKind::bands, ViewKind::profile, ViewKind::iv, ViewKind::cv, ViewKind::current,
                             ViewKind::quality, ViewKind::quality_histogram, ViewKind::difference, ViewKind::exploded}) {
        CHECK(app::find_view(app::to_string(v)) == v);
    }
    CHECK_FALSE(app::find_view("filmstrip"));  // nitcad plot's alone (Unit 26c: playback instead)
    CHECK(app::is_plot(ViewKind::iv));
    CHECK_FALSE(app::is_plot(ViewKind::field));

    using V = std::vector<ViewKind>;
    CHECK(app::available_views(diode().workspace) == V{ViewKind::profile, ViewKind::bands, ViewKind::iv});
    CHECK(app::available_views(mosfet().workspace) ==
          V{ViewKind::solid, ViewKind::mesh, ViewKind::field, ViewKind::iso, ViewKind::profile, ViewKind::bands,
            ViewKind::iv, ViewKind::quality,
            ViewKind::quality_histogram});  // profile and bands along a cut line, mesh quality (Unit 26b)
    const auto fields = app::available_fields(mosfet().workspace);
    for (const char* f : {"potential", "n", "p", "Ec", "doping", "R", "srh", "E", "J"}) CHECK(has(fields, f));
    CHECK_FALSE(has(fields, "Gii"));  // impact ionization is off: no rate stored for it
    // Before a run: the doping only.
    const fs::path dir = app_test::scratch("views_unrun");
    auto unrun = app::workspace_from_deck(app_test::write(dir, "d.deck", app_test::mosfet2d));
    REQUIRE(unrun);
    CHECK(app::available_fields(*unrun) == std::vector<std::string>{"doping", "ND", "NA"});
    CHECK(app::available_views(*unrun) == V{ViewKind::solid, ViewKind::mesh, ViewKind::field, ViewKind::profile,
                                            ViewKind::quality, ViewKind::quality_histogram});
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

TEST_CASE("app views 26b: cut lines, as nitcad plot takes them (XD5, XG4)", "[app]") {
    const app::Workspace& w = mosfet().workspace;
    const ViewSettings::Cut channel{{0.0, 0.87}, {2.0, 0.87}};  // source to drain under the oxide, um
    ViewSettings s;
    s.kind = ViewKind::bands;
    CHECK(app::build_scene(w, s).error().message == "--view bands in 2D and 3D needs --cut");
    s.kind = ViewKind::profile;
    CHECK(app::build_scene(w, s).error().message == "--view profile in 2D and 3D needs --cut");
    s.kind = ViewKind::iv;
    s.cut = channel;
    CHECK(app::build_scene(w, s).error().message ==
          "a cut line applies to the profile and bands views and to 2D and 3D device views");
    s.kind = ViewKind::bands;
    s.cut = ViewSettings::Cut{{0.0, 0.87, 0.0}, {2.0, 0.87, 0.0}};
    CHECK(app::build_scene(w, s).error().message == "a cut line needs its end points in the device's dimension");
    s.cut = ViewSettings::Cut{{0.0, 0.87}, {2.0, 0.87, 0.0}};  // the end alone
    CHECK(app::build_scene(w, s).error().message == "a cut line needs its end points in the device's dimension");
    s.cut = ViewSettings::Cut{{1.0, 0.5}, {1.0, 0.5}};
    CHECK(app::build_scene(w, s).error().message == "a cut line needs two different end points");
    s.cut = ViewSettings::Cut{{3.0, 0.5}, {4.0, 0.5}};
    CHECK(app::build_scene(w, s).error().message == "the cut line misses the device");
    ViewSettings d;
    d.kind = ViewKind::profile;
    d.cut = channel;
    CHECK(app::build_scene(diode().workspace, d).error().message ==
          "a cut line needs its end points in the device's dimension");

    // The band diagram along the line is render's, at the last point.
    s.cut = channel;
    const auto bands = app::build_scene(w, s);
    REQUIRE(bands);
    const auto cut = render::cut_line(*w.geometry, {0.0, 0.87e-4, 0.0}, {2.0e-4, 0.87e-4, 0.0});
    REQUIRE(cut);
    const auto expected = render::cut_band_plot(w.device->mesh(), *cut, w.points.back().bands);
    REQUIRE(expected);
    const auto& plot = std::get<render::PlotScene>(*bands);
    REQUIRE(plot.series.size() == 4);
    CHECK(plot.series[0].x == expected->series[0].x);
    CHECK(plot.series[0].y == expected->series[0].y);
    CHECK(plot.series[0].x.size() == cut->samples.size());
    CHECK(plot.x.label == "distance (\xC2\xB5m)");
    // A device view draws the line over everything: two batches, black under white.
    s.kind = ViewKind::field;
    const auto map = app::build_scene(w, s);
    REQUIRE(map);
    const auto& lines = std::get<render::FieldScene>(*map).lines;
    REQUIRE(lines.size() >= 2);
    CHECK(lines[lines.size() - 2].name == "cut line");
    CHECK(lines.back().name == "cut line");
    CHECK(lines.back().on_top);
    CHECK(lines.back().points.size() == 2 * (cut->samples.size() - 1));
    // The doping along it needs no results.
    const fs::path dir = app_test::scratch("views_cut_unrun");
    auto unrun = app::workspace_from_deck(app_test::write(dir, "d.deck", app_test::mosfet2d));
    REQUIRE(unrun);
    ViewSettings doping;
    doping.kind = ViewKind::profile;
    doping.field = "doping";
    doping.cut = channel;
    CHECK(app::build_scene(*unrun, doping).has_value());
    // nitcad plot's option.
    using A = std::vector<std::string>;
    CHECK(app::plot_arguments(s) == A{"--view", "field", "--field", "potential", "--cut", "0,0.87:2,0.87"});
}

TEST_CASE("app views 26b: cut lines drawn as nitcad plot draws them (XG4, WG4's choices)", "[app]") {
    const fs::path out = app_test::scratch("views_cut_images");
    const ViewSettings::Cut channel{{0.0, 0.87}, {2.0, 0.87}};
    std::vector<std::pair<ViewSettings, const char*>> cases;
    {
        ViewSettings s;
        s.kind = ViewKind::bands;
        s.cut = channel;
        cases.emplace_back(s, "bands");
    }
    {
        ViewSettings s;
        s.kind = ViewKind::profile;
        s.field = "n";
        s.log = true;
        s.point = 2;
        s.cut = ViewSettings::Cut{{1.0, 0.0}, {1.0, 1.0}};  // across the oxide: a gap in n
        cases.emplace_back(s, "profile");
    }
    {
        ViewSettings s;
        s.kind = ViewKind::field;
        s.cut = channel;
        cases.emplace_back(s, "map");
    }
    for (const auto& [s, name] : cases) {
        CAPTURE(name);
        auto scene = app::build_scene(mosfet().workspace, s);
        REQUIRE(scene);
        const render::Image mine = draw(*scene, 480, 360);
        const fs::path png = out / (std::string(name) + ".png");
        std::string args = "plot " + app_test::quoted(mosfet().deck);
        for (const std::string& a : app::plot_arguments(s)) args += " \"" + a + "\"";
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

TEST_CASE("app views 26b: run B and mesh quality, as nitcad plot takes them (XD6, XD7, XG5)", "[app]") {
    const app::Workspace& a = mosfet().workspace;
    const app::Workspace& b = mosfet_b().workspace;
    // The difference is offered with a run B on the same mesh only.
    CHECK(app::available_views(a, &b).back() == ViewKind::difference);
    const fs::path dir = app_test::scratch("views_compare_unrun");
    auto finer = app::workspace_from_deck(
        app_test::write(dir, "finer.deck", variant("x = linspace(0, 2um, 21)", "x = linspace(0, 2um, 31)")));
    REQUIRE(finer);
    CHECK(app::available_views(a, &*finer) == app::available_views(a));
    auto unrun = app::workspace_from_deck(app_test::write(dir, "unrun.deck", app_test::mosfet2d));
    REQUIRE(unrun);
    CHECK(app::available_views(a, &*unrun).back() == ViewKind::difference);

    // The refusals.
    ViewSettings s;
    s.kind = ViewKind::difference;
    CHECK(app::build_scene(a, s).error().message == "a difference map needs run B (--compare)");
    CHECK(app::build_scene(a, s, &*finer).error().message ==
          "run B was computed on another mesh: a difference needs one mesh");
    CHECK(app::build_scene(a, s, &*unrun).error().message == "run B has no bias point at run A's biases");
    s.field = "doping";  // the doping needs no bias point
    CHECK(app::build_scene(a, s, &*unrun).has_value());
    s.field = "potential";
    const auto map = app::build_scene(a, s, &b);
    REQUIRE(map);
    const auto& field = std::get<render::FieldScene>(*map);
    REQUIRE(field.colorbar);
    CHECK(field.colorbar->scale.lo() == -field.colorbar->scale.hi());  // symmetric about 0
    CHECK(field.colorbar->map == render::Colormap::coolwarm);
    CHECK(field.colorbar->title == "Δψ (mV)");

    // The curves of both runs, named by run.
    ViewSettings iv;
    iv.kind = ViewKind::iv;
    iv.sense = "drain";
    const auto curves = app::build_scene(a, iv, &b);
    REQUIRE(curves);
    const auto& plot = std::get<render::PlotScene>(*curves);
    REQUIRE(plot.series.size() == 2);
    CHECK(plot.series[0].name == "I(drain) (A)");
    CHECK(plot.series[1].name == "I(drain) (B)");
    CHECK(plot.series[1].y != plot.series[0].y);
    CHECK(app::build_scene(a, iv, &*unrun).error().message == "run B holds no sweep to compare");

    // Mesh quality: the map and the histogram; a 1D device has neither.
    ViewSettings q;
    q.kind = ViewKind::quality;
    REQUIRE(app::build_scene(a, q));
    CHECK(std::holds_alternative<render::FieldScene>(*app::build_scene(a, q)));
    q.kind = ViewKind::quality_histogram;
    CHECK(std::holds_alternative<render::PlotScene>(*app::build_scene(a, q)));
    CHECK(app::is_plot(ViewKind::quality_histogram));
    CHECK_FALSE(app::is_plot(ViewKind::quality));
    CHECK(app::build_scene(diode().workspace, q).error().message == "mesh quality needs a 2D or 3D mesh");

    // nitcad plot's options: --compare names run B, for the views that use it.
    using A = std::vector<std::string>;
    const std::string bs = b.source.string();
    CHECK(app::plot_arguments(s, &b) == A{"--view", "difference", "--compare", bs, "--field", "potential"});
    CHECK(app::plot_arguments(iv, &b) == A{"--view", "iv", "--compare", bs, "--sense", "drain"});
    CHECK(app::plot_arguments(q, &b) == A{"--view", "quality-histogram"});
}

TEST_CASE("app views 26b: run B and mesh quality drawn as nitcad plot draws them (XG5, WG4's choices)", "[app]") {
    const fs::path out = app_test::scratch("views_compare_images");
    const app::Workspace& b = mosfet_b().workspace;
    std::vector<std::pair<ViewSettings, const char*>> cases;
    {
        ViewSettings s;
        s.kind = ViewKind::quality;
        s.edges = true;
        cases.emplace_back(s, "quality");
    }
    {
        ViewSettings s;
        s.kind = ViewKind::quality_histogram;
        cases.emplace_back(s, "quality_histogram");
    }
    {
        ViewSettings s;
        s.kind = ViewKind::difference;
        s.point = 3;
        s.contours = true;
        cases.emplace_back(s, "difference_potential");
    }
    {
        ViewSettings s;
        s.kind = ViewKind::difference;
        s.field = "n";
        cases.emplace_back(s, "difference_n");
    }
    {
        ViewSettings s;
        s.kind = ViewKind::iv;
        s.sense = "drain";
        s.log = true;
        cases.emplace_back(s, "iv_overlay");
    }
    for (const auto& [s, name] : cases) {
        CAPTURE(name);
        auto scene = app::build_scene(mosfet().workspace, s, &b);
        REQUIRE(scene);
        const render::Image mine = draw(*scene, 480, 360);
        const fs::path png = out / (std::string(name) + ".png");
        std::string args = "plot " + app_test::quoted(mosfet().deck);
        for (const std::string& a : app::plot_arguments(s, &b)) args += " \"" + a + "\"";
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

TEST_CASE("app views 26c: arrows, streamlines, the exploded view and the curve pickers, as nitcad plot takes them (YD1-YD3)", "[app]") {
    using V = std::vector<ViewKind>;
    using A = std::vector<std::string>;
    // The exploded view in 3D only, after the slices.
    const auto views3d = app::available_views(mosfet3d().workspace);
    REQUIRE(views3d.size() >= 5);
    CHECK(V(views3d.begin(), views3d.begin() + 5) ==
          V{ViewKind::solid, ViewKind::mesh, ViewKind::field, ViewKind::slice, ViewKind::exploded});
    const auto views2d = app::available_views(mosfet().workspace);
    CHECK(std::ranges::find(views2d, ViewKind::exploded) == views2d.end());

    // nitcad plot's messages.
    const app::Workspace& w = mosfet().workspace;
    ViewSettings s;
    s.kind = ViewKind::field;
    const auto message = [&](const ViewSettings& t) {
        auto r = app::build_scene(w, t);
        REQUIRE_FALSE(r);
        return r.error().message;
    };
    s.vectors = "B";
    CHECK(message(s) == "--vectors is E, J, Jn or Jp");
    s.vectors.clear();
    s.streamlines = "x";
    CHECK(message(s) == "--streamlines is E, J, Jn or Jp");
    s.streamlines = "E";
    s.lines = 20000;
    CHECK(message(s) == "--lines needs a count from 1 to 10000");
    s.lines = 0;
    s.seeds = "x=1";
    CHECK(message(s) == "a seed plane needs a 3D device");
    s.seeds = "1,0.5";
    CHECK(message(s) == "--seeds needs <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>] in the device's dimension");
    s.seeds = "x=y";
    CHECK(message(s) == "--seeds needs a segment or x=<um>, y=<um> or z=<um>");
    s.seeds.clear();
    s.from = "nowhere";
    CHECK(message(s) == "there is no contact nowhere");
    s = {};
    s.kind = ViewKind::field;
    s.vectors = "E";
    s.arrows = 200000;
    CHECK(message(s) == "--arrows needs a count from 1 to 100000");
    s.arrows = 0;
    // A 1D device's field view is its profile: no device to draw them on.
    ViewSettings one;
    one.kind = ViewKind::field;
    one.vectors = "E";
    CHECK(app::build_scene(diode().workspace, one).error().message ==
          "arrows need a 2D or 3D device and a field, slice, iso or solid view");
    one.vectors.clear();
    one.streamlines = "E";
    CHECK(app::build_scene(diode().workspace, one).error().message ==
          "streamlines need a 2D or 3D device and a field, slice or iso view");
    // Off the field views they are left out (takes_vectors), as plot_arguments leaves them out.
    s.kind = ViewKind::iv;
    s.streamlines = "J";
    CHECK(app::build_scene(w, s).has_value());
    s = {};
    s.kind = ViewKind::exploded;
    CHECK(message(s) == "the exploded view needs a 3D device");
    const app::Workspace& w3 = mosfet3d().workspace;
    s.axis = 3;
    CHECK(app::build_scene(w3, s).error().message == "--axis is x, y or z");
    s.axis = 1;
    s.gap_um = -0.1;
    CHECK(app::build_scene(w3, s).error().message == "--gap needs a length in um, at least 0");
    s.gap_um.reset();
    s.vectors = "E";
    const auto with = app::build_scene(w3, s);
    s.vectors.clear();
    const auto without = app::build_scene(w3, s);
    REQUIRE(with);
    REQUIRE(without);
    CHECK(std::get<render::FieldScene>(*with).triangles.size() == std::get<render::FieldScene>(*without).triangles.size());

    // The options of nitcad plot.
    ViewSettings v;
    v.kind = ViewKind::field;
    v.field = "J";
    v.vectors = "J";
    v.arrows = 120;
    v.streamlines = "E";
    v.lines = 10;
    v.from = "gate";
    CHECK(app::plot_arguments(v) == A{"--view", "field", "--field", "J", "--streamlines", "E", "--lines", "10", "--from",
                                      "gate", "--vectors", "J", "--arrows", "120"});
    v = {};
    v.kind = ViewKind::slice;
    v.streamlines = "J";
    v.seeds = "x=1";
    CHECK(app::plot_arguments(v) == A{"--view", "slice", "--field", "potential", "--streamlines", "J", "--seeds", "x=1"});
    ViewSettings e;
    e.kind = ViewKind::exploded;
    CHECK(app::plot_arguments(e) == A{"--view", "exploded"});
    e.axis = 2;
    e.gap_um = 0.3;
    e.exploded_field = true;
    e.field = "n";
    CHECK(app::plot_arguments(e) == A{"--view", "exploded", "--field", "n", "--axis", "z", "--gap", "0.3"});
    // Arrows and streamlines kept in a setting but left out off the field views.
    e.vectors = "E";
    e.streamlines = "J";
    e.lines = 5;
    CHECK(app::plot_arguments(e) == A{"--view", "exploded", "--field", "n", "--axis", "z", "--gap", "0.3"});
    CHECK(app::takes_vectors(ViewKind::iso));
    CHECK_FALSE(app::takes_vectors(ViewKind::solid));
    // The curve pickers (YD3) are the existing options.
    ViewSettings c;
    c.kind = ViewKind::cv;
    c.contact = "anode";
    c.sense = "cathode";
    c.frequency = 1;
    c.log = true;
    CHECK(app::plot_arguments(c) == A{"--view", "cv", "--contact", "anode", "--sense", "cathode", "--frequency", "1", "--log"});
}

TEST_CASE("app views 26c: arrows, streamlines, the exploded view and the curve pickers drawn as nitcad plot draws them (YG1)", "[app]") {
    const fs::path out = app_test::scratch("views_26c_images");
    struct Case {
        const Ran* ran;
        ViewSettings s;
        const char* name;
    };
    std::vector<Case> cases;
    const auto make = [](ViewKind k, std::string field = "potential") {
        ViewSettings s;
        s.kind = k;
        s.field = std::move(field);
        return s;
    };
    // 2D: arrows on their own magnitude and over contour levels; streamlines from a contact and from seeds.
    {
        ViewSettings s = make(ViewKind::field, "J");
        s.vectors = "J";
        cases.push_back({&mosfet(), s, "2d_arrows_J"});
    }
    {
        ViewSettings s = make(ViewKind::iso, "n");
        s.levels = {1e16};
        s.vectors = "E";
        s.arrows = 50;
        cases.push_back({&mosfet(), s, "2d_iso_arrows_E"});
    }
    {
        ViewSettings s = make(ViewKind::field);
        s.streamlines = "E";
        cases.push_back({&mosfet(), s, "2d_lines_E"});
    }
    {
        ViewSettings s = make(ViewKind::field);
        s.streamlines = "J";
        s.lines = 6;
        cases.push_back({&mosfet(), s, "2d_lines_J"});  // from the contact carrying the most current
    }
    {
        ViewSettings s = make(ViewKind::field, "n");
        s.streamlines = "Jn";
        s.lines = 8;
        s.from = "source";
        s.vectors = "Jn";
        s.arrows = 80;
        cases.push_back({&mosfet(), s, "2d_lines_and_arrows_Jn"});
    }
    {
        ViewSettings s = make(ViewKind::field);
        s.streamlines = "J";
        s.seeds = "0.2,0.8:1.8,0.8";
        s.lines = 6;
        cases.push_back({&mosfet(), s, "2d_lines_seeded"});
    }
    // 3D: colored arrows on the solid, colored lines in the outline, arrows on a slice, lines from a plane.
    {
        ViewSettings s = make(ViewKind::field, "J");
        s.vectors = "J";
        cases.push_back({&mosfet3d(), s, "3d_arrows_J"});
    }
    {
        ViewSettings s = make(ViewKind::field, "E");
        s.streamlines = "E";
        s.lines = 12;
        cases.push_back({&mosfet3d(), s, "3d_lines_E"});
    }
    {
        // Arrows and lines of one field over its magnitude: the lines colored in the outline, the arrows dark.
        ViewSettings s = make(ViewKind::field, "J");
        s.vectors = "J";
        s.arrows = 50;
        s.streamlines = "J";
        s.lines = 6;
        cases.push_back({&mosfet3d(), s, "3d_lines_and_arrows_J"});
    }
    {
        ViewSettings s = make(ViewKind::slice);
        s.slices = {{2, 0.4}};
        s.vectors = "E";
        s.arrows = 60;
        cases.push_back({&mosfet3d(), s, "3d_slice_arrows_E"});
    }
    {
        ViewSettings s = make(ViewKind::field);
        s.streamlines = "J";
        s.seeds = "x=1";
        s.lines = 9;
        cases.push_back({&mosfet3d(), s, "3d_lines_plane"});
    }
    // The exploded view: the regions, and a field along x with a gap.
    cases.push_back({&mosfet3d(), make(ViewKind::exploded), "3d_exploded"});
    {
        ViewSettings s = make(ViewKind::exploded, "n");
        s.exploded_field = true;
        s.axis = 0;
        s.gap_um = 0.2;
        cases.push_back({&mosfet3d(), s, "3d_exploded_n"});
    }
    // The curves with their contacts and frequency chosen.
    {
        ViewSettings s = make(ViewKind::iv);
        s.sense = "source";
        cases.push_back({&mosfet(), s, "iv_source"});
    }
    {
        ViewSettings s = make(ViewKind::cv);
        s.sense = "cathode";
        s.frequency = 1;
        cases.push_back({&diode_cv(), s, "cv_1GHz"});
    }
    {
        ViewSettings s = make(ViewKind::current);
        s.contact = "cathode";
        cases.push_back({&diode_transient(), s, "current_cathode"});
    }
    for (const Case& c : cases) {
        CAPTURE(c.name);
        auto scene = app::build_scene(c.ran->workspace, c.s);
        if (!scene) FAIL(scene.error().message);
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
