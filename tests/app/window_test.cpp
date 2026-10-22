// The main window (Unit 26, WG5, WG7): every command in the menus and the accelerators, the panes
// re-laid out at a new DPI, the deck menus opening sections and inserting keys, a deck opened, run,
// shown and probed, the Direct3D 12 debug layer clean and no graphics object left at exit.
#include <catch2/catch_test_macros.hpp>

// The application's headers before the COM ones (objbase.h defines `interface` as a macro).
#include "main_window.hpp"
#include "NiTCAD/app/deck_sections.hpp"
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

// Live Direct3D and DXGI objects the debug layer reports, or -1 without the DXGI debug interface.
int live_objects() {
    Microsoft::WRL::ComPtr<IDXGIDebug1> debug;
    Microsoft::WRL::ComPtr<IDXGIInfoQueue> queue;
    if (FAILED(DXGIGetDebugInterface1(0, IID_PPV_ARGS(&debug))) || FAILED(DXGIGetDebugInterface1(0, IID_PPV_ARGS(&queue)))) {
        return -1;
    }
    queue->ClearStoredMessages(DXGI_DEBUG_ALL);
    debug->ReportLiveObjects(DXGI_DEBUG_ALL, static_cast<DXGI_DEBUG_RLO_FLAGS>(DXGI_DEBUG_RLO_DETAIL |
                                                                              DXGI_DEBUG_RLO_IGNORE_INTERNAL));
    int live = 0;
    const UINT64 n = queue->GetNumStoredMessages(DXGI_DEBUG_ALL);
    for (UINT64 i = 0; i < n; ++i) {
        SIZE_T size = 0;
        queue->GetMessage(DXGI_DEBUG_ALL, i, nullptr, &size);
        std::vector<char> buffer(size);
        auto* m = reinterpret_cast<DXGI_INFO_QUEUE_MESSAGE*>(buffer.data());
        if (SUCCEEDED(queue->GetMessage(DXGI_DEBUG_ALL, i, m, &size)) &&
            std::string_view(m->pDescription).starts_with("Live ")) {
            ++live;
        }
    }
    return live;
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
    const int live_with_window = live_objects();
    if (live_with_window < 0) SKIP("the DXGI debug interface is not available");
    CHECK(live_with_window > 0);  // the counter sees the renderer's objects
    w.reset();
    CHECK(live_objects() == 0);
}
