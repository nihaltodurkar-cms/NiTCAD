// The main window (Unit 26, WG5, WG7): every command in the menus and the accelerators, the panes
// re-laid out at a new DPI, the deck menus opening sections and inserting keys, a deck opened, run,
// shown and probed, the Direct3D 12 debug layer clean and no graphics object left at exit.
#include <catch2/catch_test_macros.hpp>

// The application's headers before the COM ones (objbase.h defines `interface` as a macro).
#include "main_window.hpp"
#include "win32.hpp"
#include "NiTCAD/app/deck_sections.hpp"
#include "NiTCAD/app/run_monitor.hpp"
#include "fixtures.hpp"

#include <windows.h>
#include <commctrl.h>
#include <dxgi1_3.h>
#include <dxgidebug.h>
#include <wrl/client.h>

#include <chrono>
#include <filesystem>
#include <functional>
#include <set>
#include <string>
#include <variant>
#include <vector>

using namespace NiTCAD;
namespace fs = std::filesystem;
using app::Command;
using app::MainWindow;

namespace {

// Pumps the thread's messages until the condition holds (or the time runs out).
bool pump_until(const std::function<bool()>& done, double seconds = 120.0) {
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

std::unique_ptr<MainWindow> window(bool debug_layer = false) {
    app::WindowOptions o;
    o.adapter = render::Adapter::warp;
    o.debug_layer = debug_layer;
    o.show = false;
    o.quit_on_close = false;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    return w;
}

// Every command id in a menu and its submenus, after letting the window fill the dynamic ones.
void collect(HMENU menu, std::set<UINT>& ids, std::vector<std::wstring>* texts = nullptr) {
    for (int i = 0; i < GetMenuItemCount(menu); ++i) {
        if (HMENU sub = GetSubMenu(menu, i)) {
            collect(sub, ids, texts);
            continue;
        }
        const UINT id = GetMenuItemID(menu, i);
        if (id != 0 && id != static_cast<UINT>(-1)) ids.insert(id);
        if (texts) {
            wchar_t buffer[256]{};
            GetMenuStringW(menu, static_cast<UINT>(i), buffer, 256, MF_BYPOSITION);
            texts->push_back(buffer);
        }
    }
}

void open_popup(MainWindow& w, int index) {
    HMENU popup = GetSubMenu(GetMenu(w.hwnd()), index);
    SendMessageW(w.hwnd(), WM_INITMENUPOPUP, reinterpret_cast<WPARAM>(popup), MAKELPARAM(index, FALSE));
}

// The id of the first item of a popup whose text starts with `prefix` (in a submenu too).
UINT find_item(HMENU menu, const std::wstring& prefix) {
    for (int i = 0; i < GetMenuItemCount(menu); ++i) {
        wchar_t buffer[256]{};
        GetMenuStringW(menu, static_cast<UINT>(i), buffer, 256, MF_BYPOSITION);
        if (HMENU sub = GetSubMenu(menu, i)) {
            if (const UINT id = find_item(sub, prefix)) return id;
            continue;
        }
        if (std::wstring(buffer).starts_with(prefix)) return GetMenuItemID(menu, i);
    }
    return 0;
}

RECT client_rect_of(HWND child, HWND parent) {
    RECT r{};
    GetWindowRect(child, &r);
    MapWindowPoints(nullptr, parent, reinterpret_cast<POINT*>(&r), 2);
    return r;
}

}  // namespace

TEST_CASE("app window: every command in the menus and the accelerators (WG5)", "[app]") {
    auto w = window();
    HMENU bar = GetMenu(w->hwnd());
    REQUIRE(GetMenuItemCount(bar) == 7);
    for (int i = 0; i < 7; ++i) open_popup(*w, i);
    std::set<UINT> ids;
    collect(bar, ids);
    for (const app::CommandInfo& c : app::commands()) {
        CAPTURE(c.text);
        CHECK(ids.contains(static_cast<UINT>(c.id)));
    }
    // The accelerator table holds every command's keys (Alt+F4 is the system's).
    const int n = CopyAcceleratorTableW(w->accelerators(), nullptr, 0);
    std::vector<ACCEL> table(static_cast<std::size_t>(n));
    CopyAcceleratorTableW(w->accelerators(), table.data(), n);
    std::size_t expected = 0;
    for (const app::CommandInfo& c : app::commands()) {
        if (c.accelerator.key == 0 || c.id == Command::exit) continue;
        ++expected;
        const bool found = std::ranges::any_of(table, [&](const ACCEL& a) {
            return a.cmd == static_cast<WORD>(c.id) && a.key == c.accelerator.key &&
                   ((a.fVirt & FCONTROL) != 0) == c.accelerator.ctrl && ((a.fVirt & FSHIFT) != 0) == c.accelerator.shift;
        });
        CAPTURE(c.text);
        CHECK(found);
    }
    CHECK(table.size() == expected);
    // The deck menus open their sections: each section of the reference's menus is listed.
    std::vector<std::wstring> texts;
    std::set<UINT> device_ids;
    collect(GetSubMenu(bar, 2), device_ids, &texts);
    CHECK(std::ranges::find(texts, L"Go to [device] (not in the deck)") != texts.end());
    CHECK(std::ranges::find(texts, L"Go to [contact] (not in the deck)") != texts.end());
}

TEST_CASE("app window: per-monitor DPI changes re-lay the window out (WG5)", "[app]") {
    auto w = window();
    const int before = w->dpi();
    RECT r{};
    GetWindowRect(w->hwnd(), &r);
    // As Windows sends it when the window moves to a 200% monitor: the new DPI and a suggested rectangle.
    RECT suggested{r.left, r.top, r.left + 2 * (r.right - r.left), r.top + 2 * (r.bottom - r.top)};
    SendMessageW(w->hwnd(), WM_DPICHANGED, MAKEWPARAM(192, 192), reinterpret_cast<LPARAM>(&suggested));
    CHECK(w->dpi() == 192);
    RECT client{};
    GetClientRect(w->hwnd(), &client);
    // The panes are the layout's at 192 DPI, and the children sit in them.
    const app::Panes& p = w->panes();
    CHECK(p.tree.width == app::scaled(app::tree_width_96, 192));
    CHECK(p.inspector.width == app::scaled(app::inspector_width_96, 192));
    CHECK(p.inspector.right() == client.right);
    const RECT tree = client_rect_of(w->child(MainWindow::id_tree), w->hwnd());
    CHECK(tree.left == p.tree.x);
    CHECK(tree.right == p.tree.right());
    const RECT inspector = client_rect_of(w->child(MainWindow::id_inspector), w->hwnd());
    CHECK(inspector.left == p.inspector.x);
    CHECK(inspector.right - inspector.left == p.inspector.width);
    // Fonts at the new DPI: 9 points.
    LOGFONTW lf{};
    GetObjectW(reinterpret_cast<HFONT>(SendMessageW(w->child(MainWindow::id_tree), WM_GETFONT, 0, 0)), sizeof(lf), &lf);
    CHECK(lf.lfHeight == -MulDiv(9, 192, 72));
    // And back.
    SendMessageW(w->hwnd(), WM_DPICHANGED, MAKEWPARAM(before, before), reinterpret_cast<LPARAM>(&r));
    CHECK(w->panes().tree.width == app::scaled(app::tree_width_96, before));
}

TEST_CASE("app window: a deck opened, run, shown, probed; the deck menus (WG7)", "[app]") {
    const fs::path dir = app_test::scratch("window_run");
    const fs::path deck = app_test::write(dir, "m.deck", app_test::mosfet2d);
    auto w = window();
    REQUIRE(w->renderer());
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr; }));
    CHECK_FALSE(w->workspace()->has_results());
    REQUIRE(w->scene());
    CHECK(std::holds_alternative<render::FieldScene>(*w->scene()));  // the regions, before a run
    CHECK_FALSE(w->viewport().error());
    CHECK(w->document().text() == app_test::mosfet2d);

    // Run (F5): the worker solves, the window shows the potential and the I-V curve.
    w->execute(Command::run);
    CHECK(w->runner().busy());
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    CHECK(w->status(0).starts_with("Done in"));
    CHECK(w->tab() == MainWindow::tab_view);
    CHECK(w->settings().kind == app::ViewKind::field);
    CHECK(w->settings().field == "potential");
    REQUIRE(w->scene());
    CHECK_FALSE(w->viewport().error());
    w->select_tab(MainWindow::tab_curves);
    CHECK(w->settings().kind == app::ViewKind::iv);
    REQUIRE(w->scene());
    CHECK(std::holds_alternative<render::PlotScene>(*w->scene()));
    w->select_tab(MainWindow::tab_view);
    // A probe at the middle of the view: a node and its potential.
    w->probe_at(static_cast<float>(w->viewport().width()) / 2, static_cast<float>(w->viewport().height()) / 2);
    REQUIRE(w->probe());
    CHECK(std::isfinite(w->probe()->value));
    // Camera moves redraw without errors; Fit restores the view.
    w->execute(Command::zoom_in);
    CHECK(w->settings().zoom == 1.25f);
    w->execute(Command::pan_left);
    CHECK(w->settings().pan_x == -0.1f);
    w->execute(Command::fit);
    CHECK(w->settings().zoom == 1.0f);
    CHECK_FALSE(w->viewport().error());

    // The Physics menu inserts a [models] key with its default (the deck has no [models]: it is added).
    open_popup(*w, 4);
    HMENU physics = GetSubMenu(GetMenu(w->hwnd()), 4);
    const UINT srh = find_item(physics, L"srh = on");
    REQUIRE(srh != 0);
    w->execute(static_cast<Command>(srh));
    CHECK(w->tab() == MainWindow::tab_deck);
    const std::string text = w->deck().text();
    const auto models = app::find_section(text, "models");
    REQUIRE(models);
    REQUIRE(models->entries.size() == 1);
    CHECK(models->entries[0].key == "srh");
    CHECK(w->deck().caret_line() == models->entries[0].line);
    CHECK(w->document().modified());
    CHECK(app::check_deck(text, deck).diagnostics.empty());
    // The Simulation menu goes to [sweep].
    open_popup(*w, 5);
    const UINT sweep = find_item(GetSubMenu(GetMenu(w->hwnd()), 5), L"Go to [sweep]");
    REQUIRE(sweep != 0);
    w->execute(static_cast<Command>(sweep));
    CHECK(w->deck().caret_line() == app::find_section(text, "sweep")->line);
    // A deck error goes to its line (F7 then F8).
    std::string bad = text;
    bad.insert(bad.find("[contact source]"), "[bogus]\n");
    w->deck().set_text(bad);
    w->execute(Command::check);
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    REQUIRE(w->runner().diagnostics().size() == 1);
    CHECK(w->tab() == MainWindow::tab_deck);
    CHECK(w->deck().caret_line() == w->runner().diagnostics()[0].line);
    CHECK(w->runner().diagnostics()[0].line > 0);
}

TEST_CASE("app window: the start page, the run monitor, Stop and the results (WG1, WG8, VA3, VA5)", "[app]") {
    const fs::path copies = app_test::scratch("window_examples");
    app::WindowOptions o;
    o.adapter = render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    o.examples_copy_dir = copies;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    // The start page lists exactly the bundled examples (WG8).
    CHECK(w->tab() == MainWindow::tab_start);
    const auto expected = app::list_examples(NITCAD_EXAMPLES_DIR);
    REQUIRE(expected.size() == 4);
    REQUIRE(w->start_page().example_count() == 4);
    std::size_t mosfet = 4;
    for (std::size_t k = 0; k < 4; ++k) {
        CHECK(w->start_page().example_title(static_cast<int>(k)) == app::win32::widen(expected[k].title));
        if (expected[k].path.filename() == "mosfet_2d.deck") mosfet = k;
    }
    REQUIRE(mosfet < 4);
    // Opening an example opens a copy (the installed deck stays as shipped).
    REQUIRE(w->open_example(mosfet));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr; }));
    CHECK(w->document().path() == copies / "mosfet_2d.deck");
    CHECK(fs::exists(copies / "mosfet_2d.deck"));

    // Run: the run monitor follows the points; when done, the potential and the figures.
    w->execute(Command::run);
    CHECK(w->tab() == MainWindow::tab_run);
    int rows_while_running = 0;
    REQUIRE(pump_until([&] {
        if (w->runner().busy()) rows_while_running = std::max(rows_while_running, w->run_page().row_count());
        return !w->runner().busy();
    }));
    CHECK(rows_while_running > 1);  // the monitor followed the run, not only its end
    CHECK(w->tab() == MainWindow::tab_view);
    CHECK(w->run_page().row_count() == static_cast<int>(w->runner().points().size()));
    CHECK(w->run_page().cell(0, 1) == L"equilibrium");
    CHECK(w->run_page().cell(1, 2) == L"1 of 16");
    // Every cell as the final rows format it: the cells kept between paints are not stale (WG1).
    const auto& final_rows = w->runner().points();
    for (std::size_t r = 0; r < final_rows.size(); ++r) {
        const auto cells = app::monitor_cells(final_rows[r], r);
        for (std::size_t c = 0; c < cells.size(); ++c) {
            CAPTURE(r, c);
            CHECK(w->run_page().cell(static_cast<int>(r), static_cast<int>(c)) == app::win32::widen(cells[c]));
        }
    }
    CHECK_FALSE(w->run_page().plot().error());
    w->select_tab(MainWindow::tab_results);
    CHECK(w->results_page().figure_count() >= 4);
    CHECK(w->results_page().figure_cell(0, 0) == L"Threshold (max gm)");
    CHECK(w->results_page().figure_cell(0, 3) == L"I(drain) vs V(gate)");
    CHECK(w->results_page().record_count() > 10);

    // Stop after three points: the run ends early, its completed part kept and shown (WG1).
    w->execute(Command::run);
    bool asked = false;
    REQUIRE(pump_until([&] {
        const auto converged = std::ranges::count_if(w->runner().points(), [](const app::PointRow& r) {
            return r.phase == solve::Phase::bias && r.converged;
        });
        if (!asked && converged >= 3) {
            w->execute(Command::stop);
            asked = true;
        }
        return asked && !w->runner().busy();
    }));
    CHECK(w->status(0).starts_with("Stopped early"));
    REQUIRE(w->workspace());
    REQUIRE(w->workspace()->sweep);
    CHECK(w->workspace()->sweep->points.size() >= 3);
    CHECK(w->workspace()->sweep->points.size() < 16);
    CHECK(w->tab() == MainWindow::tab_view);
    CHECK(w->scene());
    CHECK_FALSE(w->viewport().error());
    // Every tab can be shown, from the keyboard too.
    for (int k = 0; k < MainWindow::tab_count; ++k) w->execute(Command::next_tab);
    CHECK(w->tab() == MainWindow::tab_view);
}

TEST_CASE("app window: the start page lists the examples of a given directory (WG8)", "[app]") {
    const fs::path dir = app_test::scratch("window_examples_dir");
    app_test::write(dir, "b.deck", "# Second\n# the second deck\n");
    app_test::write(dir, "a.deck", "# First\n# the first deck\n");
    app::WindowOptions o;
    o.adapter = render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    o.examples_dir = dir;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    REQUIRE(w->start_page().example_count() == 2);
    CHECK(w->start_page().example_title(0) == L"First");
    CHECK(w->start_page().example_title(1) == L"Second");
}

TEST_CASE("app window: the debug layer stays clean and no graphics object outlives the window (WG5)", "[app]") {
    const fs::path dir = app_test::scratch("window_debug");
    const fs::path deck = app_test::write(dir, "m.deck", app_test::mosfet2d);
    auto w = window(true);
    if (!w->renderer()) SKIP("the Direct3D 12 debug layer is not installed (Graphics Tools): " + w->renderer_error().value_or(""));
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr; }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace()->has_results(); }));
    for (const Command c : {Command::zoom_in, Command::pan_up, Command::fit, Command::tab_curves, Command::tab_view}) {
        w->execute(c);
    }
    CHECK_FALSE(w->viewport().error());
    const auto messages = w->renderer()->debug_messages();
    for (const auto& m : messages) UNSCOPED_INFO(m);
    CHECK(messages.empty());
    const int live_with_window = app::win32::live_graphics_objects();
    if (live_with_window < 0) SKIP("the DXGI debug interface is not available");
    CHECK(live_with_window > 0);  // the counter sees the renderer's objects
    w.reset();
    CHECK(app::win32::live_graphics_objects() == 0);
}

TEST_CASE("app window: the self-test's window takes no input (WG1)", "[app]") {
    app::WindowOptions o;
    o.adapter = render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    auto usual = MainWindow::create(o);
    REQUIRE(usual);
    CHECK(IsWindowEnabled(usual->hwnd()));
    // Shown without being activated and disabled: the desktop's mouse and keys do not reach it, so the
    // message loop's gaps measure the window's own work; it still paints and runs its timers.
    o.input = false;
    auto quiet = MainWindow::create(o);
    REQUIRE(quiet);
    CHECK_FALSE(IsWindowEnabled(quiet->hwnd()));
    quiet->refresh();
}
