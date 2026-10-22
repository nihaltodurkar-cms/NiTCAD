// The analysis workspace in the window (Unit 26b, XG1, XG3, XG7): the deck checked while it is typed
// (also during a run) with suggestions and the structure preview; completion by the keyboard; cut lines (XG4).
#include <catch2/catch_test_macros.hpp>

// The application's headers before the COM ones (objbase.h defines `interface` as a macro).
#include "main_window.hpp"
#include "fixtures.hpp"
#include "NiTCAD/render/image.hpp"

#include <windows.h>
#include <richedit.h>
#include <commctrl.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <string>
#include <variant>
#include <vector>

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

namespace {

// A click on the device view at a point of the device (um), through the viewport's own messages.
void click_at(MainWindow& w, const render::Float3& um) {
    REQUIRE(w.scene());
    const auto& f = std::get<render::FieldScene>(*w.scene());
    const auto s = render::screen_position(f, w.viewport().width(), w.viewport().height(), um);
    const LPARAM at = MAKELPARAM(static_cast<int>(std::lround(s.x)), static_cast<int>(std::lround(s.y)));
    SendMessageW(w.viewport().hwnd(), WM_LBUTTONDOWN, MK_LBUTTON, at);
    SendMessageW(w.viewport().hwnd(), WM_LBUTTONUP, 0, at);
}

void apply_cut(MainWindow& w, const std::string& text) {
    SetWindowTextW(w.child(MainWindow::id_cut), std::wstring(text.begin(), text.end()).c_str());
    SendMessageW(GetParent(w.child(MainWindow::id_apply)), WM_COMMAND, MAKEWPARAM(MainWindow::id_apply, BN_CLICKED),
                 reinterpret_cast<LPARAM>(w.child(MainWindow::id_apply)));
}

}  // namespace

TEST_CASE("app window 26b: a cut line drawn by two clicks and typed (XD5, XG4, XG7)", "[app]") {
    const fs::path dir = app_test::scratch("window_cut");
    const fs::path deck = app_test::write(dir, "mosfet.deck", app_test::mosfet2d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_view);

    // Before a run: the doping along the line, typed; the Curves tab offers the profile in 2D.
    apply_cut(*w, "0,0.87:2,0.87");
    REQUIRE(w->cut());
    CHECK(w->tab() == MainWindow::tab_curves);
    CHECK(w->settings().kind == app::ViewKind::profile);
    CHECK(w->settings().field == "doping");
    CHECK(SendMessageW(w->child(MainWindow::id_view), CB_FINDSTRINGEXACT, static_cast<WPARAM>(-1),
                       reinterpret_cast<LPARAM>(L"profile")) != CB_ERR);
    // The run keeps it on the device it was drawn on.
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace()->has_results(); }));
    CHECK(w->cut());
    apply_cut(*w, "");
    w->select_tab(MainWindow::tab_view);

    // The tool: two clicks on the map, from the source to the drain under the oxide.
    w->execute(Command::cut_line_mode);
    CHECK(w->status(0) == "Draw cut line: click the line's start, then its end.");
    click_at(*w, {0.0f, 0.87f, 0.0f});
    CHECK_FALSE(w->cut());
    CHECK(w->tab() == MainWindow::tab_view);
    click_at(*w, {2.0f, 0.87f, 0.0f});
    REQUIRE(w->cut());
    // A pixel is about 4 nm here: the ends within two pixels, to the nanometre.
    for (std::size_t k = 0; k < 2; ++k) {
        CHECK(std::abs(w->cut()->from_um[k] - (k == 0 ? 0.0 : 0.87)) <= 0.01);
        CHECK(std::abs(w->cut()->to_um[k] - (k == 0 ? 2.0 : 0.87)) <= 0.01);
        CHECK(w->cut()->from_um[k] == std::round(w->cut()->from_um[k] * 1e3) / 1e3);
    }
    // Its band diagram in the Curves tab.
    CHECK(w->tab() == MainWindow::tab_curves);
    CHECK(w->settings().kind == app::ViewKind::bands);
    REQUIRE(w->scene());
    const auto& plot = std::get<render::PlotScene>(*w->scene());
    CHECK(plot.title == "Band diagram");
    CHECK(plot.series.size() == 4);
    // And on the device view, over the map.
    w->select_tab(MainWindow::tab_view);
    REQUIRE(w->scene());
    CHECK(std::get<render::FieldScene>(*w->scene()).lines.back().name == "cut line");

    // Typed in the inspector: the line moves; empty removes it.
    apply_cut(*w, "0.5,0.2:1.5,0.2");
    REQUIRE(w->cut());
    CHECK(w->cut()->from_um == std::vector<double>{0.5, 0.2});
    CHECK(w->cut()->to_um == std::vector<double>{1.5, 0.2});
    CHECK(w->tab() == MainWindow::tab_curves);
    w->select_tab(MainWindow::tab_view);
    apply_cut(*w, "0.5,0.2");
    CHECK(w->status(0).starts_with("A cut line is"));
    REQUIRE(w->cut());
    apply_cut(*w, "");
    CHECK_FALSE(w->cut());
    REQUIRE(w->scene());
    const auto& lines = std::get<render::FieldScene>(*w->scene()).lines;
    CHECK(std::ranges::none_of(lines, [](const render::LineBatch& b) { return b.name == "cut line"; }));
    // Probe mode leaves the tool.
    w->execute(Command::probe_mode);
    click_at(*w, {1.0f, 0.5f, 0.0f});
    CHECK_FALSE(w->cut());
    CHECK(w->probe());
}

TEST_CASE("app window 26b: a cut line on a 3D preview, between surface nodes (XD5)", "[app]") {
    const fs::path dir = app_test::scratch("window_cut3d");
    const fs::path deck = app_test::write(dir, "resistor.deck", app_test::generate3d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_view);
    w->execute(Command::cut_line_mode);
    // Two surface nodes in front, apart on the screen: the front-most of all, then of those 60 px from it.
    REQUIRE(w->scene());
    const auto& f = std::get<render::FieldScene>(*w->scene());
    const auto& m = w->workspace()->device->mesh();
    std::vector<char> surface(m.node_count(), 0);
    for (const render::Polygon& p : w->workspace()->geometry->boundary_faces()) {
        for (const render::SurfaceVertex& v : p.vertices) surface[static_cast<std::size_t>(v.at.nodes[0])] = 1;
    }
    const auto um = [&](std::size_t i) {
        const auto& p = m.points()[i];
        return render::Float3{static_cast<float>(p[0] * 1e4), static_cast<float>(p[1] * 1e4), static_cast<float>(p[2] * 1e4)};
    };
    const auto front = [&](const std::function<bool(const render::ScreenPoint&)>& allowed) {
        std::size_t best = 0;
        float depth = 2.0f;
        for (std::size_t i = 0; i < m.node_count(); ++i) {
            if (!surface[i]) continue;
            const auto s = render::screen_position(f, w->viewport().width(), w->viewport().height(), um(i));
            if (allowed(s) && s.depth < depth) {
                depth = s.depth;
                best = i;
            }
        }
        return best;
    };
    const std::size_t a = front([](const render::ScreenPoint&) { return true; });
    const auto sa = render::screen_position(f, w->viewport().width(), w->viewport().height(), um(a));
    const std::size_t b = front([&](const render::ScreenPoint& s) { return std::hypot(s.x - sa.x, s.y - sa.y) > 60.0f; });
    REQUIRE(a != b);
    // What each click picks: the front-most surface node within 12 px of the pixel clicked.
    const auto picked = [&](std::size_t i) {
        const auto s = render::screen_position(f, w->viewport().width(), w->viewport().height(), um(i));
        const float x = std::round(s.x), y = std::round(s.y);
        return front([&](const render::ScreenPoint& t) { return std::hypot(t.x - x, t.y - y) <= 12.0f; });
    };
    const std::size_t pa = picked(a), pb = picked(b);
    click_at(*w, um(a));
    click_at(*w, um(b));
    REQUIRE(w->cut());
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(w->cut()->from_um[k] == std::round(m.points()[pa][k] * 1e4 * 1e3) / 1e3);
        CHECK(w->cut()->to_um[k] == std::round(m.points()[pb][k] * 1e4 * 1e3) / 1e3);
    }
    // Not solved: the doping along it.
    CHECK(w->tab() == MainWindow::tab_curves);
    CHECK(w->settings().kind == app::ViewKind::profile);
    CHECK(w->settings().field == "doping");
    REQUIRE(w->scene());
    CHECK(std::holds_alternative<render::PlotScene>(*w->scene()));
}

namespace {

std::wstring metrics(MainWindow& w) {
    wchar_t text[512]{};
    GetWindowTextW(w.child(MainWindow::id_metrics), text, 512);
    return text;
}

std::string narrow(const std::wstring& s) {
    std::string out(static_cast<std::size_t>(WideCharToMultiByte(CP_UTF8, 0, s.c_str(), -1, nullptr, 0, nullptr, nullptr)),
                    '\0');
    WideCharToMultiByte(CP_UTF8, 0, s.c_str(), -1, out.data(), static_cast<int>(out.size()), nullptr, nullptr);
    out.pop_back();
    return out;
}

}  // namespace

TEST_CASE("app window 26b: run B from a result file and from the previous run; mesh quality (XD6, XD7, XG7)", "[app]") {
    const fs::path dir = app_test::scratch("window_compare");
    const fs::path deck = app_test::write(dir, "mosfet.deck", app_test::mosfet2d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    CHECK_FALSE(w->has_previous_run());
    CHECK_FALSE(w->run_b());
    // The first run's files, kept as a result to compare with.
    fs::copy_file(dir / "mosfet.nitsweep", dir / "first.nitsweep");
    fs::copy_file(dir / "mosfet.nitdev", dir / "first.nitdev");

    // Mesh quality: the map, its figures in the metrics strip (0.1 by 0.05 um cells: every triangle atan 1/2,
    // 26.6°).
    w->execute(Command::mesh_quality);
    CHECK(w->tab() == MainWindow::tab_view);
    CHECK(w->settings().kind == app::ViewKind::quality);
    REQUIRE(w->scene());
    CHECK(std::get<render::FieldScene>(*w->scene()).title == "Mesh quality");
    CAPTURE(narrow(metrics(*w)));
    CHECK(metrics(*w).find(L"Min angle 26.6°, median 26.6°, 0 below 20°") != std::wstring::npos);

    // Another gate work function, run: the first run becomes the previous one.
    replace(*w, "work_function = 4.1", "work_function = 4.3");
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->has_previous_run(); }));
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    w->execute(Command::compare_previous);
    REQUIRE(w->run_b());
    CHECK_FALSE(w->has_previous_run());
    CHECK(w->status(0) == "Comparing with the previous run");
    // The difference map, the metrics naming run B with ΔV_th, the results side by side.
    w->select_tab(MainWindow::tab_view);
    w->settings().kind = app::ViewKind::difference;
    w->settings().field = "potential";
    w->refresh();
    REQUIRE(w->scene());
    CHECK(std::get<render::FieldScene>(*w->scene()).colorbar->title == "\xCE\x94\xCF\x88 (mV)");
    CHECK(metrics(*w).find(L"Run B mosfet.deck, ΔVth -") != std::wstring::npos);  // B's threshold is lower
    w->select_tab(MainWindow::tab_results);
    CHECK(w->results_page().figure_cell(0, 0) == L"Threshold (max gm)");
    CHECK_FALSE(w->results_page().figure_cell(0, 5).empty());
    CHECK(w->results_page().figure_cell(0, 6).starts_with(L"-"));
    bool marked = false;
    for (int r = 0; r < w->results_page().record_count(); ++r) {
        if (w->results_page().record_cell(r, 0) == L"≠ Input identity") marked = true;
    }
    CHECK(marked);

    // Stop comparing: the difference map goes, the previous run stays to compare with again.
    w->execute(Command::stop_comparing);
    CHECK_FALSE(w->run_b());
    CHECK(w->has_previous_run());
    w->select_tab(MainWindow::tab_view);
    CHECK(w->settings().kind == app::ViewKind::field);
    CHECK(metrics(*w).find(L"Run B") == std::wstring::npos);

    // A result file as run B: the first run, the same mesh.
    REQUIRE(w->compare_result(dir / "first.nitsweep", dir / "first.nitdev"));
    REQUIRE(w->run_b());
    CHECK(w->status(0) == "Comparing with first.nitsweep");
    const auto& views = app::available_views(*w->workspace(), w->run_b());
    CHECK(std::ranges::find(views, app::ViewKind::difference) != views.end());
    // Opening another deck forgets both.
    REQUIRE(w->open_deck(deck));
    CHECK_FALSE(w->run_b());
    CHECK_FALSE(w->has_previous_run());
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
}

namespace {

// Unit 26c: an inspector list's entry chosen, or a box's text applied, as the user would.
void choose(MainWindow& w, int id, const wchar_t* item) {
    HWND list = w.child(id);
    const LRESULT k = SendMessageW(list, CB_FINDSTRINGEXACT, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(item));
    REQUIRE(k != CB_ERR);
    SendMessageW(list, CB_SETCURSEL, static_cast<WPARAM>(k), 0);
    SendMessageW(GetParent(list), WM_COMMAND, MAKEWPARAM(id, CBN_SELCHANGE), reinterpret_cast<LPARAM>(list));
}

void click(MainWindow& w, int id) {
    HWND b = w.child(id);
    if ((GetWindowLongPtrW(b, GWL_STYLE) & BS_TYPEMASK) == BS_AUTOCHECKBOX) {
        SendMessageW(b, BM_SETCHECK, SendMessageW(b, BM_GETCHECK, 0, 0) == BST_CHECKED ? BST_UNCHECKED : BST_CHECKED, 0);
    }
    SendMessageW(GetParent(b), WM_COMMAND, MAKEWPARAM(id, BN_CLICKED), reinterpret_cast<LPARAM>(b));
}

void apply_text(MainWindow& w, int id, const std::string& text) {
    SetWindowTextW(w.child(id), std::wstring(text.begin(), text.end()).c_str());
    click(w, MainWindow::id_apply);
}

bool shown(MainWindow& w, int id) { return (GetWindowLongPtrW(w.child(id), GWL_STYLE) & WS_VISIBLE) != 0; }

// The window's view equals nitcad plot's image for its settings, drawn on the same adapter (WARP).
void same_as_plot(MainWindow& w, const fs::path& deck, const fs::path& png) {
    REQUIRE(w.scene());
    auto mine = std::visit([&](const auto& s) { return w.renderer()->draw(s, 480, 360); }, *w.scene());
    REQUIRE(mine);
    std::string args = "plot " + app_test::quoted(deck);
    for (const std::string& a : app::plot_arguments(w.settings())) args += " \"" + a + "\"";
    args += " --size 480x360 --adapter warp --out " + app_test::quoted(png);
    const auto cli = app_test::nitcad(args);
    CAPTURE(args, cli.output);
    REQUIRE(cli.status == 0);
    auto theirs = render::load_png(png);
    REQUIRE(theirs);
    const auto d = render::compare(*mine, *theirs, 0);
    REQUIRE(d);
    CHECK(d->max_channel == 0);
}

}  // namespace

TEST_CASE("app window 26c: arrows, streamlines and the exploded view in the inspector (YD1, YD2, YG1, YG5)", "[app]") {
    const fs::path dir = app_test::scratch("window_vectors");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_view);
    // Before a run: no arrows (they need the results); the exploded view's controls with that view only.
    CHECK_FALSE(shown(*w, MainWindow::id_arrows));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    REQUIRE(w->settings().kind == app::ViewKind::field);
    for (const int id : {MainWindow::id_arrows, MainWindow::id_arrow_count, MainWindow::id_lines, MainWindow::id_line_count,
                         MainWindow::id_from, MainWindow::id_seeds}) {
        CAPTURE(id);
        CHECK(shown(*w, id));
    }
    for (const int id : {MainWindow::id_axis, MainWindow::id_gap, MainWindow::id_exploded_field, MainWindow::id_swept,
                         MainWindow::id_sensed, MainWindow::id_frequency, MainWindow::id_log}) {
        CAPTURE(id);
        CHECK_FALSE(shown(*w, id));
    }

    // Arrows of J: the map shows |J|, the arrows colored by it on the solid device (YD1).
    choose(*w, MainWindow::id_arrows, L"J");
    CHECK(w->settings().vectors == "J");
    CHECK(w->settings().field == "J");
    apply_text(*w, MainWindow::id_arrow_count, "40");
    CHECK(w->settings().arrows == 40);
    apply_text(*w, MainWindow::id_arrow_count, "0");
    CHECK(w->status(0) == "Arrows are a count from 1 to 100000.");
    CHECK(w->settings().arrows == 40);
    same_as_plot(*w, deck, dir / "arrows.png");
    // Streamlines of E from the gate, then from a seed plane, beside the arrows.
    choose(*w, MainWindow::id_lines, L"E");
    CHECK(w->settings().streamlines == "E");
    CHECK(w->settings().field == "E");
    choose(*w, MainWindow::id_from, L"gate");
    CHECK(w->settings().from == "gate");
    apply_text(*w, MainWindow::id_line_count, "7");
    CHECK(w->settings().lines == 7);
    same_as_plot(*w, deck, dir / "lines_gate.png");
    apply_text(*w, MainWindow::id_seeds, " x=1 ");
    CHECK(w->settings().seeds == "x=1");
    same_as_plot(*w, deck, dir / "lines_plane.png");
    choose(*w, MainWindow::id_arrows, L"none");
    CHECK(w->settings().vectors.empty());
    CHECK(w->settings().field == "E");

    // The exploded view: the arrows' controls hidden (none drawn there), its own shown.
    choose(*w, MainWindow::id_view, L"exploded");
    CHECK(w->settings().kind == app::ViewKind::exploded);
    CHECK_FALSE(shown(*w, MainWindow::id_lines));
    CHECK(shown(*w, MainWindow::id_axis));
    CHECK(shown(*w, MainWindow::id_gap));
    REQUIRE(w->scene());
    same_as_plot(*w, deck, dir / "exploded.png");
    choose(*w, MainWindow::id_axis, L"x");
    CHECK(w->settings().axis == 0);
    click(*w, MainWindow::id_exploded_field);
    CHECK(w->settings().exploded_field);
    apply_text(*w, MainWindow::id_gap, "0.2");
    CHECK(w->settings().gap_um == 0.2);
    apply_text(*w, MainWindow::id_gap, "-1");
    CHECK(w->status(0) == "The gap is a length in µm, at least 0.");
    CHECK(w->settings().gap_um == 0.2);
    same_as_plot(*w, deck, dir / "exploded_field.png");
    // No probe on the moved regions.
    // (A click where the top face's centre node projects would pick it on any other view.)
    click_at(*w, {1.0f, 1.0f, 0.4f});
    CHECK_FALSE(w->probe());
    choose(*w, MainWindow::id_view, L"solid");
    click_at(*w, {1.0f, 1.0f, 0.4f});
    CHECK(w->probe());
    choose(*w, MainWindow::id_view, L"exploded");
    // Back on the field view the streamlines are drawn again.
    choose(*w, MainWindow::id_view, L"field");
    CHECK(w->settings().streamlines == "E");
    CHECK(shown(*w, MainWindow::id_lines));
    same_as_plot(*w, deck, dir / "lines_again.png");

    // Keyboard: every new control is a tab stop (YG5); at 200% the pairs stay side by side in the inspector.
    for (int id = MainWindow::id_arrows; id <= MainWindow::id_log; ++id) {
        CAPTURE(id);
        CHECK((GetWindowLongPtrW(w->child(id), GWL_STYLE) & WS_TABSTOP) != 0);
    }
    RECT r{};
    GetWindowRect(w->hwnd(), &r);
    RECT suggested{r.left, r.top, r.left + 2 * (r.right - r.left), r.top + 2 * (r.bottom - r.top)};
    SendMessageW(w->hwnd(), WM_DPICHANGED, MAKEWPARAM(192, 192), reinterpret_cast<LPARAM>(&suggested));
    RECT list{}, box{}, inspector{};
    GetWindowRect(w->child(MainWindow::id_lines), &list);
    GetWindowRect(w->child(MainWindow::id_line_count), &box);
    GetWindowRect(w->child(MainWindow::id_inspector), &inspector);
    CHECK(list.right <= box.left);
    CHECK(list.top == box.top);
    CHECK(box.right <= inspector.right);
    CHECK(box.bottom - box.top >= 40);  // 23 px rows at 96 dpi, doubled
}

TEST_CASE("app window 26c: the curve pickers (YD3, YG1)", "[app]") {
    const fs::path dir = app_test::scratch("window_pickers");
    const fs::path deck = app_test::write(dir, "cv.deck", app_test::diode_cv());
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_curves);
    REQUIRE(w->settings().kind == app::ViewKind::cv);
    for (const int id : {MainWindow::id_swept, MainWindow::id_sensed, MainWindow::id_frequency, MainWindow::id_log}) {
        CAPTURE(id);
        CHECK(shown(*w, id));
    }
    CHECK_FALSE(shown(*w, MainWindow::id_arrows));
    // The two frequencies, the contacts by name; the defaults first.
    CHECK(SendMessageW(w->child(MainWindow::id_frequency), CB_GETCOUNT, 0, 0) == 2);
    CHECK(SendMessageW(w->child(MainWindow::id_sensed), CB_GETCOUNT, 0, 0) == 3);
    choose(*w, MainWindow::id_frequency, L"1e+09 Hz");
    CHECK(w->settings().frequency == 1);
    choose(*w, MainWindow::id_sensed, L"cathode");
    CHECK(w->settings().sense == "cathode");
    choose(*w, MainWindow::id_swept, L"anode");
    CHECK(w->settings().contact == "anode");
    click(*w, MainWindow::id_log);
    CHECK(w->settings().log);
    click(*w, MainWindow::id_log);
    CHECK_FALSE(w->settings().log);
    same_as_plot(*w, deck, dir / "cv.png");
    // The default again: nitcad plot's choice.
    choose(*w, MainWindow::id_sensed, L"(automatic)");
    CHECK(w->settings().sense.empty());
    same_as_plot(*w, deck, dir / "cv_default.png");
}

namespace {

// Unit 26c (YD4): the scene tree's items under the device, by the start of their text.
HTREEITEM tree_item(MainWindow& w, const std::wstring& prefix) {
    HWND tree = w.child(MainWindow::id_tree);
    for (HTREEITEM root = TreeView_GetRoot(tree); root; root = TreeView_GetNextSibling(tree, root)) {
        for (HTREEITEM i = TreeView_GetChild(tree, root); i; i = TreeView_GetNextSibling(tree, i)) {
            wchar_t text[128]{};
            TVITEMW t{};
            t.mask = TVIF_TEXT;
            t.hItem = i;
            t.pszText = text;
            t.cchTextMax = 128;
            TreeView_GetItem(tree, &t);
            if (std::wstring(text).starts_with(prefix)) return i;
        }
    }
    return nullptr;
}

int check_state(MainWindow& w, HTREEITEM item) {
    return static_cast<int>(TreeView_GetCheckState(w.child(MainWindow::id_tree), item)) + 1;  // 0 none, 1 off, 2 on
}

// The item's box toggled by the space bar, as the user would.
void press_space(MainWindow& w, HTREEITEM item) {
    HWND tree = w.child(MainWindow::id_tree);
    TreeView_SelectItem(tree, item);
    SendMessageW(tree, WM_KEYDOWN, VK_SPACE, 0);
    SendMessageW(tree, WM_CHAR, L' ', 0);
    SendMessageW(tree, WM_KEYUP, VK_SPACE, 0);
    pump_for(0.2);  // the window reads the boxes when its posted message comes
}

}  // namespace

TEST_CASE("app window 26c: regions hidden from the scene tree (YD4, YG1, YG5)", "[app]") {
    const fs::path dir = app_test::scratch("window_regions");
    const fs::path deck = app_test::write(dir, "mosfet.deck", app_test::mosfet2d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    REQUIRE(w->settings().kind == app::ViewKind::field);
    // A box on each region, checked; none on the contacts or the other items.
    HTREEITEM silicon = tree_item(*w, L"silicon");
    HTREEITEM oxide = tree_item(*w, L"oxide");
    HTREEITEM gate = tree_item(*w, L"contact gate");
    REQUIRE((silicon && oxide && gate));
    CHECK(check_state(*w, silicon) == 2);
    CHECK(check_state(*w, oxide) == 2);
    CHECK(check_state(*w, gate) == 0);
    HMENU bar = GetMenu(w->hwnd());
    CHECK((GetMenuState(bar, static_cast<UINT>(Command::show_regions), MF_BYCOMMAND) & MF_GRAYED) != 0);
    // A probe on the oxide's top before it is hidden.
    click_at(*w, {1.0f, 1.0f, 0.0f});
    REQUIRE(w->probe());
    const auto before = w->probe()->node;

    // The oxide hidden by its box: the view as nitcad plot --hide oxide draws it, no probe where it was.
    press_space(*w, oxide);
    CHECK(check_state(*w, oxide) == 1);
    CHECK(w->settings().hidden == std::vector<std::string>{"oxide"});
    CHECK(w->status(0) == "Hiding the region oxide.");
    CHECK((GetMenuState(bar, static_cast<UINT>(Command::show_regions), MF_BYCOMMAND) & MF_GRAYED) == 0);
    same_as_plot(*w, deck, dir / "field_no_oxide.png");
    click_at(*w, {1.0f, 1.0f, 0.0f});
    CHECK(w->status(0) == "No node within reach.");
    CHECK(w->probe()->node == before);  // the last probe, kept
    click_at(*w, {1.0f, 0.5f, 0.0f});
    REQUIRE(w->probe());
    CHECK(w->probe()->position_um[1] < 0.93);
    // The setting stays on the other device views.
    choose(*w, MainWindow::id_view, L"mesh");
    CHECK(w->settings().hidden == std::vector<std::string>{"oxide"});
    same_as_plot(*w, deck, dir / "mesh_no_oxide.png");
    // The space bar on a contact gives it no box and hides nothing.
    press_space(*w, gate);
    CHECK(check_state(*w, gate) == 0);
    CHECK(w->settings().hidden == std::vector<std::string>{"oxide"});
    // Both hidden, in the device's order; the silicon shown again by its box.
    press_space(*w, silicon);
    CHECK(w->settings().hidden == std::vector<std::string>{"silicon", "oxide"});
    press_space(*w, silicon);
    CHECK(w->settings().hidden == std::vector<std::string>{"oxide"});
    CHECK(check_state(*w, silicon) == 2);
    // Show every region clears it and checks the boxes.
    w->execute(Command::show_regions);
    CHECK(w->settings().hidden.empty());
    CHECK(check_state(*w, tree_item(*w, L"oxide")) == 2);
    CHECK(w->status(0) == "Showing every region.");
    same_as_plot(*w, deck, dir / "mesh_all.png");
    // A new run resets the view settings, the hidden regions with them, and the boxes show it.
    press_space(*w, tree_item(*w, L"oxide"));
    CHECK(w->settings().hidden == std::vector<std::string>{"oxide"});
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    CHECK(w->settings().hidden.empty());
    CHECK(check_state(*w, tree_item(*w, L"oxide")) == 2);
    // In the View menu with its keys (YG5).
    const auto* info = app::find_command(Command::show_regions);
    REQUIRE(info);
    CHECK(info->menu == "&View");
    CHECK(info->accelerator.key == 'H');
}
