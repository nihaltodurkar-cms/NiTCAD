// The analysis workspace in the window (Unit 26b, XG1, XG3, XG7): the deck checked while it is typed
// (also during a run) with suggestions and the structure preview; completion by the keyboard.
#include <catch2/catch_test_macros.hpp>

// The application's headers before the COM ones (objbase.h defines `interface` as a macro).
#include "main_window.hpp"
#include "fixtures.hpp"

#include <windows.h>
#include <richedit.h>

#include <chrono>
#include <filesystem>
#include <functional>
#include <string>

using namespace NiTCAD;
namespace fs = std::filesystem;
using app::Command;
using app::MainWindow;

namespace {

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

void pump_for(double seconds) {
    const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    pump_until([&] { return std::chrono::steady_clock::now() > end; }, seconds + 1);
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

// Replaces the first `from` in the deck editor's text by `to`, as typing would (the host is told).
void replace(MainWindow& w, const std::string& from, const std::string& to) {
    const std::string text = w.deck().text();
    const auto at = text.find(from);
    REQUIRE(at != std::string::npos);
    // ASCII text: a character per byte, a line end ('\r' in Rich Edit) one character.
    SendMessageW(w.deck().editor(), EM_SETSEL, at, at + from.size());
    const std::wstring wide(to.begin(), to.end());
    SendMessageW(w.deck().editor(), EM_REPLACESEL, TRUE, reinterpret_cast<LPARAM>(wide.c_str()));
}

void type(MainWindow& w, const std::wstring& keys) {
    for (const wchar_t c : keys) SendMessageW(w.deck().editor(), WM_CHAR, c, 0);
}

void key(MainWindow& w, UINT vk) {
    SendMessageW(w.deck().editor(), WM_KEYDOWN, vk, 0);
    if (vk == VK_RETURN || vk == VK_TAB || vk == VK_ESCAPE) {
        SendMessageW(w.deck().editor(), WM_CHAR, vk == VK_RETURN ? L'\r' : vk == VK_TAB ? L'\t' : 0x1b, 0);
    }
}

}  // namespace

TEST_CASE("app window 26b: the deck checked while typed, suggestions and the structure preview (XG1, XG3)", "[app]") {
    const fs::path dir = app_test::scratch("window_live");
    const fs::path deck = app_test::write(dir, "mosfet.deck", app_test::mosfet2d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    const auto checks_before = w->live_checker().submitted();
    w->select_tab(MainWindow::tab_deck);

    // A misspelt material: its line marked, the suggestion given, without Check (the run controller idle).
    replace(*w, "material = sio2", "material = sio3");
    replace(*w, "material = sio3", "material = sioo2");
    REQUIRE(pump_until([&] { return !w->deck().diagnostics().empty(); }));
    CHECK(w->runner().state() == app::RunState::idle);
    CHECK(w->live_checker().submitted() == checks_before + 1);  // one check for both edits (the pause)
    CHECK(w->deck().diagnostics().front().line == 11);
    CHECK(w->deck().diagnostics().front().message.ends_with("(did you mean 'sio2'?)"));
    CHECK(w->status(0) == w->deck().diagnostics().front().message);

    // Fixed: no error, and the deck's structure shown, not solved.
    replace(*w, "material = sioo2", "material = sio2");
    REQUIRE(pump_until([&] { return w->deck().diagnostics().empty() && w->preview(); }));
    CHECK_FALSE(w->workspace()->has_results());
    CHECK(w->workspace()->device->mesh().node_count() == 21 * 21);
    // A check that ends after a newer edit is not shown: its error is no longer in the text.
    {
        const auto submitted = w->live_checker().submitted();
        replace(*w, "material = sio2", "material = sio9");
        REQUIRE(pump_until([&] { return w->live_checker().submitted() == submitted + 1; }));
        replace(*w, "material = sio9", "material = sio2");  // before the first check's result is taken
        REQUIRE(pump_until([&] { return w->live_checker().applied() >= submitted + 1; }));
        CHECK(w->deck().diagnostics().empty());
        REQUIRE(pump_until([&] { return w->live_checker().applied() == submitted + 2; }));
        CHECK(w->deck().diagnostics().empty());
    }
    // A finer mesh: the preview follows, in the view chosen.
    w->select_tab(MainWindow::tab_view);
    w->settings().kind = app::ViewKind::mesh;
    w->select_tab(MainWindow::tab_deck);
    replace(*w, "linspace(0, 2um, 21)", "linspace(0, 2um, 41)");
    REQUIRE(pump_until([&] { return w->workspace()->device->mesh().node_count() == 41 * 21; }));
    CHECK(w->preview());
    w->select_tab(MainWindow::tab_view);
    CHECK(w->settings().kind == app::ViewKind::mesh);
    REQUIRE(w->scene());
    CHECK(std::holds_alternative<render::FieldScene>(*w->scene()));

    // Run it; then an edit keeping the mesh keeps the results, one changing it previews.
    w->execute(Command::save);
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace()->has_results(); }));
    CHECK_FALSE(w->preview());
    w->select_tab(MainWindow::tab_deck);
    const auto before = w->live_checker().applied();
    replace(*w, "contact drain = 0.5", "contact drain = 0.4");
    REQUIRE(pump_until([&] { return w->live_checker().applied() > before; }));
    pump_for(0.1);
    CHECK(w->workspace()->has_results());
    CHECK_FALSE(w->preview());
    replace(*w, "linspace(0, 1um, 21)", "linspace(0, 1um, 41)");
    REQUIRE(pump_until([&] { return w->preview(); }));
    CHECK(w->workspace()->device->mesh().node_count() == 41 * 41);
}

TEST_CASE("app window 26b: the deck is checked while a run goes on (XG1)", "[app]") {
    const fs::path dir = app_test::scratch("window_live_run");
    const fs::path deck = app_test::write(dir, "diode.deck", app_test::diode_sweep("0 : 0.8 : 0.0005"));
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return w->runner().state() == app::RunState::running; }));
    w->select_tab(MainWindow::tab_deck);
    replace(*w, "kind = ohmic", "kind = omic");
    REQUIRE(pump_until([&] { return !w->deck().diagnostics().empty(); }));
    CHECK(w->runner().state() == app::RunState::running);
    CHECK(w->deck().diagnostics().front().line == 14);
    CHECK(w->deck().diagnostics().front().message.ends_with("(did you mean 'ohmic'?)"));
    w->execute(Command::stop);
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
}

TEST_CASE("app window 26b: completion by the keyboard (XD3, XG7)", "[app]") {
    auto w = window();
    w->execute(Command::new_deck);
    REQUIRE(w->tab() == MainWindow::tab_deck);
    // Ctrl+Space (the Complete command) after "[mo": the section kind.
    type(*w, L"[mo");
    CHECK(w->deck().completion_visible());  // typing opens it too
    w->deck().hide_completion();
    w->execute(Command::complete);
    REQUIRE(w->deck().completion_visible());
    REQUIRE(w->deck().completion()->items.size() == 1);
    key(*w, VK_RETURN);
    CHECK_FALSE(w->deck().completion_visible());
    CHECK(w->deck().text() == "[models]");  // Enter took the item, not a new line
    // A key typed at a line's start, then its value.
    SendMessageW(w->deck().editor(), EM_REPLACESEL, TRUE, reinterpret_cast<LPARAM>(L"\r"));
    type(*w, L"sr");
    REQUIRE(w->deck().completion_visible());
    CHECK(w->deck().completion()->items.front().label == "srh");
    key(*w, VK_TAB);
    CHECK(w->deck().text() == "[models]\nsrh = ");
    // The value's choices follow at once.
    REQUIRE(w->deck().completion_visible());
    CHECK(w->deck().completion()->items.size() == 2);
    key(*w, VK_DOWN);
    CHECK(w->deck().completion_selection() == 1);
    key(*w, VK_DOWN);
    CHECK(w->deck().completion_selection() == 1);  // the last stays selected
    key(*w, VK_UP);
    key(*w, VK_DOWN);
    key(*w, VK_RETURN);
    CHECK(w->deck().text() == "[models]\nsrh = off");
    // Ctrl+Space on an empty line: every key not yet there.
    SendMessageW(w->deck().editor(), EM_REPLACESEL, TRUE, reinterpret_cast<LPARAM>(L"\r"));
    w->execute(Command::complete);
    REQUIRE(w->deck().completion_visible());
    CHECK(w->deck().completion()->items.size() == 15);
    key(*w, VK_DOWN);
    key(*w, VK_DOWN);
    CHECK(w->deck().completion_selection() == 2);  // one item per key
    // Escape closes the list and types nothing.
    key(*w, VK_ESCAPE);
    type(*w, L"au");
    REQUIRE(w->deck().completion_visible());
    key(*w, VK_ESCAPE);
    CHECK_FALSE(w->deck().completion_visible());
    CHECK(w->deck().text() == "[models]\nsrh = off\nau");
    // Nothing fits: the status says so.
    type(*w, L"#");
    w->execute(Command::complete);
    CHECK_FALSE(w->deck().completion_visible());
    CHECK(w->status(0) == "Nothing to complete here.");
    // A double click on an item takes it.
    SendMessageW(w->deck().editor(), EM_REPLACESEL, TRUE, reinterpret_cast<LPARAM>(L"\r"));
    type(*w, L"[");
    REQUIRE(w->deck().completion_visible());
    SendMessageW(w->deck().completion_list(), LB_SETCURSEL, 1, 0);
    SendMessageW(w->deck().hwnd(), WM_COMMAND, MAKEWPARAM(3, LBN_DBLCLK),
                 reinterpret_cast<LPARAM>(w->deck().completion_list()));
    CHECK(w->deck().text() == "[models]\nsrh = off\nau#\n[mesh]");
}
