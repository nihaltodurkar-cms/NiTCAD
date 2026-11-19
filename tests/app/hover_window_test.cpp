// The hover readout in the window (Unit 26e, DECISIONS.md HD1-HD3, HG4, HG7): a mouse move over a field view
// shows the node Probe picks there with its value, over a curve the nearest data point; it hides off the device,
// on leaving the view, while dragging, on a tab change and when switched off; the switch is saved.

#include <catch2/catch_test_macros.hpp>

#include <windows.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include "NiTCAD/app/app_state.hpp"
#include "NiTCAD/app/probe.hpp"
#include "NiTCAD/render/edges.hpp"
#include "NiTCAD/render/image.hpp"
#include "fixtures.hpp"
#include "main_window.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
using app::Command;
using app::MainWindow;

namespace {

bool pump_until(const std::function<bool()>& done, double seconds = 600.0) {
    const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    MSG m;
    while (!done()) {
        if (std::chrono::steady_clock::now() > end) return false;
        if (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
            TranslateMessage(&m);
            DispatchMessageW(&m);
        } else {
            MsgWaitForMultipleObjects(0, nullptr, FALSE, 10, QS_ALLINPUT);
        }
    }
    return true;
}

std::unique_ptr<MainWindow> window() {
    app::WindowOptions o;
    o.adapter = render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    return w;
}

// A move to the view's client pixel, as the mouse sends it (no button down).
void move(MainWindow& w, float x, float y) {
    SendMessageW(w.viewport().hwnd(), WM_MOUSEMOVE, 0,
                 MAKELPARAM(static_cast<int>(std::lround(x)), static_cast<int>(std::lround(y))));
}

std::unique_ptr<MainWindow> solved(const std::string& name, const std::string& deck_text) {
    const fs::path dir = app_test::scratch(name);
    const fs::path deck = app_test::write(dir, "d.deck", deck_text);
    auto w = window();
    REQUIRE(w->renderer());
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    return w;
}

}  // namespace

TEST_CASE("app window 26e: hover readout on a 2D field view (HD1, HD3, HG4)", "[app]") {
    auto w = solved("hover_2d", app_test::generate2d);
    w->select_tab(MainWindow::tab_view);
    REQUIRE(w->settings().field == "potential");
    REQUIRE(w->scene());
    const auto& f = std::get<render::FieldScene>(*w->scene());
    const int W = w->viewport().width(), H = w->viewport().height();
    CHECK(w->hover_readout());
    CHECK_FALSE(w->readout_shown());

    // Onto a node in the middle of the device: Probe's node and its potential.
    const auto& m = w->workspace()->device->mesh();
    std::size_t node = 0;
    double best = 1e300;
    for (std::size_t i = 0; i < m.node_count(); ++i) {
        const double d = std::hypot(m.points()[i][0] - 1e-4, m.points()[i][1] - 0.5e-4);
        if (d < best) best = d, node = i;
    }
    const auto& p = m.points()[node];
    const render::ScreenPoint at =
        render::screen_position(f, W, H, {static_cast<float>(p[0] * 1e4), static_cast<float>(p[1] * 1e4), 0.0f});
    move(*w, at.x, at.y);
    REQUIRE(w->readout_shown());
    w->probe_at(std::round(at.x), std::round(at.y));
    REQUIRE(w->probe());
    CHECK(static_cast<std::size_t>(w->probe()->node) == node);
    CHECK(std::isfinite(w->probe()->value));
    CHECK(w->readout_text() == app::hover_text("ψ (V)", *w->probe(), 2));
    CHECK(w->readout_text().starts_with("ψ = "));
    CHECK(w->status(3).starts_with("ψ = "));

    // Over empty space: hidden.
    move(*w, 1.0f, 1.0f);
    CHECK_FALSE(w->readout_shown());
    // Leaving the view hides it.
    move(*w, at.x, at.y);
    REQUIRE(w->readout_shown());
    SendMessageW(w->viewport().hwnd(), WM_MOUSELEAVE, 0, 0);
    CHECK_FALSE(w->readout_shown());
    // A drag hides it and shows nothing while it lasts.
    move(*w, at.x, at.y);
    REQUIRE(w->readout_shown());
    const LPARAM xy = MAKELPARAM(static_cast<int>(std::lround(at.x)), static_cast<int>(std::lround(at.y)));
    SendMessageW(w->viewport().hwnd(), WM_LBUTTONDOWN, MK_LBUTTON, xy);
    CHECK_FALSE(w->readout_shown());
    SendMessageW(w->viewport().hwnd(), WM_MOUSEMOVE, MK_LBUTTON, xy);
    CHECK_FALSE(w->readout_shown());
    SendMessageW(w->viewport().hwnd(), WM_LBUTTONUP, 0, xy);
    // A tab change hides it.
    move(*w, at.x, at.y);
    REQUIRE(w->readout_shown());
    w->select_tab(MainWindow::tab_curves);
    CHECK_FALSE(w->readout_shown());
    w->select_tab(MainWindow::tab_view);

    // Switched off: nothing shows; the switch is in the saved state.
    w->execute(Command::hover_readout);
    CHECK_FALSE(w->hover_readout());
    REQUIRE(w->scene());
    move(*w, at.x, at.y);
    CHECK_FALSE(w->readout_shown());
    CHECK_FALSE(app::parse_state(app::serialize(w->app_state())).hover_readout);
    HMENU bar = GetMenu(w->hwnd());
    CHECK((GetMenuState(bar, static_cast<UINT>(Command::hover_readout), MF_BYCOMMAND) & MF_CHECKED) == 0);
    w->execute(Command::hover_readout);
    CHECK(w->hover_readout());
    CHECK((GetMenuState(bar, static_cast<UINT>(Command::hover_readout), MF_BYCOMMAND) & MF_CHECKED) != 0);
    CHECK(app::parse_state(app::serialize(w->app_state())).hover_readout);
    move(*w, at.x, at.y);
    CHECK(w->readout_shown());
}

TEST_CASE("app window 26e: hover readout on curves, the 1D device view and the I-V (HD2, HG4)", "[app]") {
    auto w = solved("hover_1d", app_test::diode_sweep());
    for (const int tab : {static_cast<int>(MainWindow::tab_view), static_cast<int>(MainWindow::tab_curves)}) {
        w->select_tab(tab);
        REQUIRE(w->scene());
        REQUIRE(std::holds_alternative<render::PlotScene>(*w->scene()));
        const auto& s = std::get<render::PlotScene>(*w->scene());
        REQUIRE_FALSE(s.series.empty());
        const auto& series = s.series.back();
        REQUIRE_FALSE(series.x.empty());
        const std::size_t i = series.x.size() / 2;
        const int W = w->viewport().width(), H = w->viewport().height();
        const render::ScreenPoint at = render::plot_position(s, W, H, series.x[i], series.y[i]);
        REQUIRE(std::isfinite(at.x));
        move(*w, at.x + 1.0f, at.y);
        REQUIRE(w->readout_shown());
        const auto hit = app::pick_curve(s, W, H, std::round(at.x + 1.0f), std::round(at.y));
        REQUIRE(hit);
        CHECK(w->readout_text() == app::curve_text(s, *hit));
        CHECK(w->readout_text().find(s.x.label + " = ") != std::string::npos);
        move(*w, 1.0f, 1.0f);
        CHECK_FALSE(w->readout_shown());
    }
}

TEST_CASE("app window 26e: a hover costs less than a frame on the largest 3D example (HG7)", "[app]") {
    const fs::path deck = fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_3d.deck";
    REQUIRE(fs::exists(deck));
    auto w = window();
    REQUIRE(w->renderer());
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_view);
    w->settings().kind = app::ViewKind::field;
    w->settings().field = "doping";
    w->refresh();
    REQUIRE(w->scene());
    REQUIRE(std::holds_alternative<render::FieldScene>(*w->scene()));
    const int W = w->viewport().width(), H = w->viewport().height();
    std::vector<double> ms;
    int shown = 0;
    for (int k = 0; k < 200; ++k) {
        const float x = W * (0.2f + 0.6f * static_cast<float>(k % 20) / 19.0f);
        const float y = H * (0.2f + 0.6f * static_cast<float>(k / 20) / 9.0f);
        const auto t0 = std::chrono::steady_clock::now();
        move(*w, x, y);
        ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        shown += w->readout_shown() ? 1 : 0;
    }
    std::ranges::sort(ms);
    const double median = 0.5 * (ms[99] + ms[100]);
    WARN("HG7: median hover " << median << " ms over 200 positions (" << shown << " shown), max " << ms.back()
                               << " ms, " << w->workspace()->device->mesh().node_count() << " nodes");
    CHECK(shown > 0);
    CHECK(median < 16.0);
}

// Unit 26f (KG6), on the largest bundled 3D example's doping surface. Times recorded in DECISIONS.md.
namespace {

double median_of(std::vector<double> ms) {
    std::ranges::sort(ms);
    return 0.5 * (ms[(ms.size() - 1) / 2] + ms[ms.size() / 2]);
}

std::unique_ptr<MainWindow> example_3d_doping() {
    const fs::path deck = fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_3d.deck";
    REQUIRE(fs::exists(deck));
    auto w = window();
    REQUIRE(w->renderer());
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_view);
    w->settings().kind = app::ViewKind::field;
    w->settings().field = "doping";
    w->refresh();
    REQUIRE(w->scene());
    REQUIRE(std::holds_alternative<render::FieldScene>(*w->scene()));
    return w;
}

}  // namespace

TEST_CASE("app window 26f: the enhanced look's edges are built with the scene, under 50 ms (KG6)", "[app]") {
    auto w = example_3d_doping();
    std::vector<double> build;
    std::size_t points = 0;
    for (int k = 0; k < 21; ++k) {
        render::FieldScene f = std::get<render::FieldScene>(*w->scene());
        render::set_look(f, render::Look::classic);
        const auto t0 = std::chrono::steady_clock::now();
        render::set_look(f, render::Look::enhanced);
        build.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        points = f.lines.back().points.size();
    }
    const double edges = median_of(build);
    WARN("KG6: edges built in " << edges << " ms (median of 21, " << points / 2 << " segments, "
                                << w->workspace()->device->mesh().node_count() << " nodes)");
    // Release only (owner, 2026-11-19): the unoptimised Debug build takes about 20x longer; there it is recorded.
#ifdef NDEBUG
    CHECK(edges < 50.0);
#endif
}

// The redraw clause, measured on the hardware adapter by the owner's ruling (KG6 (b), 2026-11-18): the
// scene drawn off screen (no vsync), its camera turned 1 degree before each draw and nothing rebuilt.
// Skipped where there is no hardware adapter (CI). The test window's WARP drag is reported only:
// there the per-pixel lighting runs on the CPU (about 1.25x classic).
TEST_CASE("app window 26f: a camera-move redraw within 10% of classic on the hardware adapter (KG6)", "[app]") {
    auto hw = render::Renderer::create({.adapter = render::Adapter::hardware});
    if (!hw) SKIP("no Direct3D 12 hardware adapter: " + hw.error().message);
    auto w = example_3d_doping();
    const int W = w->viewport().width(), H = w->viewport().height();
    render::FieldScene f = std::get<render::FieldScene>(*w->scene());
    double t[2]{};
    for (int mode = 0; mode < 2; ++mode) {
        render::set_look(f, mode ? render::Look::enhanced : render::Look::classic);
        std::vector<double> ms;
        for (int k = 0; k < 55; ++k) {
            f.camera.azimuth_deg += 1.0f;
            const auto t0 = std::chrono::steady_clock::now();
            REQUIRE(hw->draw(f, W, H));
            if (k >= 5) ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        }
        t[mode] = median_of(ms);
    }
    // For the record: the window's own drag on WARP, each look.
    const HWND view = w->viewport().hwnd();
    const int cy = w->viewport().height() / 2;
    const auto orbits = [&] {
        std::vector<double> ms;
        int x = 100;
        SendMessageW(view, WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(x, cy));
        for (int k = 0; k < 55; ++k) {
            x += 3;
            const auto t0 = std::chrono::steady_clock::now();
            SendMessageW(view, WM_MOUSEMOVE, MK_LBUTTON, MAKELPARAM(x, cy));
            const double d = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            if (k >= 5) ms.push_back(d);
        }
        SendMessageW(view, WM_LBUTTONUP, 0, MAKELPARAM(x, cy));
        return median_of(ms);
    };
    REQUIRE(w->settings().enhanced);
    const double warp_enhanced = orbits();
    w->settings().enhanced = false;
    w->refresh();
    const double warp_classic = orbits();
    WARN("KG6: redraw on " << hw->adapter_name() << " " << t[1] << " ms enhanced, " << t[0]
                           << " ms classic (medians of 50), ratio " << t[1] / t[0] << "; WARP drag " << warp_enhanced
                           << " ms enhanced, " << warp_classic << " ms classic, ratio " << warp_enhanced / warp_classic);
    // Release only (owner, 2026-11-19): in the unoptimised Debug build the ratio wanders past 1.10 (1.11 and 1.19
    // against a pass); there it is recorded.
#ifdef NDEBUG
    CHECK(t[1] <= 1.10 * t[0]);
#endif
}

// Unit 26g (MG7), by KG6's method on the hardware adapter: the 26g enhanced draw within 1.25x of 26f's.
// 26f's draw is the same scene with its batches matte and no ground shadow: exactly 26f's triangles and
// lines (set_look puts every edge in the ink), leaving only 26g's shader branches, which matte skips.
// WARP is recorded only.
TEST_CASE("app window 26g: the enhanced draw within 1.25x of 26f's on the hardware adapter (MG7)", "[app]") {
    auto w = example_3d_doping();
    const int W = w->viewport().width(), H = w->viewport().height();
    render::FieldScene now = std::get<render::FieldScene>(*w->scene());
    render::set_look(now, render::Look::enhanced);
    render::FieldScene before = now;
    for (render::TriangleBatch& b : before.triangles) b.finish = render::Finish::matte;
    before.ground_shadow = false;
    render::set_look(before, render::Look::enhanced);
    REQUIRE(std::ranges::any_of(now.triangles, [](const render::TriangleBatch& b) { return b.finish != render::Finish::matte; }));
    const auto time = [&](render::Renderer& r) {
        double t[2]{};
        for (int mode = 0; mode < 2; ++mode) {
            render::FieldScene& f = mode ? now : before;
            std::vector<double> ms;
            for (int k = 0; k < 55; ++k) {
                f.camera.azimuth_deg += 1.0f;
                const auto t0 = std::chrono::steady_clock::now();
                REQUIRE(r.draw(f, W, H));
                if (k >= 5) ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
            }
            t[mode] = median_of(ms);
        }
        return std::array<double, 2>{t[0], t[1]};
    };
    auto warp = render::Renderer::create({.adapter = render::Adapter::warp});
    REQUIRE(warp);
    const auto tw = time(*warp);
    WARN("MG7: WARP " << tw[1] << " ms 26g, " << tw[0] << " ms 26f, ratio " << tw[1] / tw[0] << " (recorded only)");
    auto hw = render::Renderer::create({.adapter = render::Adapter::hardware});
    if (!hw) SKIP("no Direct3D 12 hardware adapter: " + hw.error().message);
    const auto th = time(*hw);
    WARN("MG7: " << hw->adapter_name() << " " << th[1] << " ms 26g, " << th[0] << " ms 26f (medians of 50), ratio "
                 << th[1] / th[0]);
    CHECK(th[1] <= 1.25 * th[0]);
}
