// The analysis workspace in the window (Unit 26b, XG1, XG3, XG7): the deck checked while it is typed
// (also during a run) with suggestions and the structure preview; completion by the keyboard; cut lines (XG4).
#include <catch2/catch_test_macros.hpp>

// The application's headers before the COM ones (objbase.h defines `interface` as a macro).
#include "main_window.hpp"
#include "viewport.hpp"
#include "fixtures.hpp"
#include "NiTCAD/app/probe.hpp"
#include "NiTCAD/render/band_surface.hpp"
#include "NiTCAD/render/image.hpp"

#include <windows.h>
#include <richedit.h>
#include <shlobj.h>  // DROPFILES (Unit 26d, ZG4)
#include <commctrl.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <filesystem>
#include <format>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <variant>
#include <vector>

using namespace NiTCAD;
namespace fs = std::filesystem;
using app::Command;
using app::MainWindow;

namespace {

// A guard against a hang, not a gate: the hosted CI runner's Debug build ran the playback test's transient
// MOSFET for more than 120 s.
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


// Unit 26d (ZD1): text typed into an inspector box (its change flag set, as typing sets it), then Enter through
// the window's message filter (as the message loop passes it) or focus leaving the box.
void type_box(MainWindow& w, int id, const std::string& text) {
    HWND box = w.child(id);
    SetWindowTextW(box, std::wstring(text.begin(), text.end()).c_str());
    SendMessageW(box, EM_SETMODIFY, TRUE, 0);
}

void enter_box(MainWindow& w, int id, const std::string& text) {
    type_box(w, id, text);
    MSG m{};
    m.hwnd = w.child(id);
    m.message = WM_KEYDOWN;
    m.wParam = VK_RETURN;
    REQUIRE(w.pre_translate(m));
}

void leave_box(MainWindow& w, int id, const std::string& text) {
    type_box(w, id, text);
    HWND box = w.child(id);
    SendMessageW(GetParent(box), WM_COMMAND, MAKEWPARAM(id, EN_KILLFOCUS), reinterpret_cast<LPARAM>(box));
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
    CHECK(w->deck().completion()->items.size() == 31);  // the 32 [models] keys less srh (Unit 62: + 4 dg_ keys; Unit 67: + 3 TAT keys; Unit 68: + 4 polarization keys)
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

void apply_cut(MainWindow& w, const std::string& text) { enter_box(w, MainWindow::id_cut, text); }

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

void apply_text(MainWindow& w, int id, const std::string& text) { enter_box(w, id, text); }

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

namespace {

// A slider moved by the user to a step (its notification to the inspector, as a drag ends).
void slide(MainWindow& w, int id, int step) {
    HWND bar = w.child(id);
    SendMessageW(bar, TBM_SETPOS, TRUE, step);
    SendMessageW(GetParent(bar), WM_HSCROLL, MAKEWPARAM(TB_THUMBPOSITION, step), reinterpret_cast<LPARAM>(bar));
}

int position(MainWindow& w, int id) { return static_cast<int>(SendMessageW(w.child(id), TBM_GETPOS, 0, 0)); }

std::string text_of(MainWindow& w, int id) {
    wchar_t buffer[256]{};
    GetWindowTextW(w.child(id), buffer, 256);
    std::string out;
    for (const wchar_t* c = buffer; *c; ++c) out += static_cast<char>(*c < 128 ? *c : '?');
    return out;
}

bool grayed(MainWindow& w, Command c) {
    return (GetMenuState(GetMenu(w.hwnd()), static_cast<UINT>(c), MF_BYCOMMAND) & MF_GRAYED) != 0;
}

// The window's view equals the image at `png` (drawn on WARP at 480 x 360).
void same_as_png(MainWindow& w, const fs::path& png) {
    REQUIRE(w.scene());
    auto mine = std::visit([&](const auto& s) { return w.renderer()->draw(s, 480, 360); }, *w.scene());
    REQUIRE(mine);
    auto theirs = render::load_png(png);
    REQUIRE(theirs);
    const auto d = render::compare(*mine, *theirs, 0);
    REQUIRE(d);
    CHECK(d->max_channel == 0);
}

}  // namespace

TEST_CASE("app window 26c: playing a transient run's snapshots (YD5, YG3, YG5)", "[app]") {
    const fs::path dir = app_test::scratch("window_playback");
    const fs::path deck = app_test::write(dir, "mosfet.deck", app_test::mosfet2d_transient());
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    REQUIRE(w->settings().kind == app::ViewKind::field);
    REQUIRE(w->settings().field == "potential");
    const std::size_t n = w->workspace()->points.size();
    REQUIRE(n >= 3);
    // The frames nitcad plot draws for the run, each snapshot's.
    const fs::path frames = dir / "frames";
    const auto cli = app_test::nitcad("plot " + app_test::quoted(deck) +
                                      " --view filmstrip --field potential --every 1 --colormap viridis --size 480x360 --adapter warp --frames " +
                                      app_test::quoted(frames) + " --out " + app_test::quoted(dir / "strip.png"));
    CAPTURE(cli.output);
    REQUIRE(cli.status == 0);
    const auto frame = [&](std::size_t k) {
        char name[32];
        std::snprintf(name, sizeof name, "frame_%03zu.png", k);
        return frames / name;
    };
    // The controls and the commands on a transient run's field view; the last snapshot shown, as its frame.
    CHECK(shown(*w, MainWindow::id_time));
    CHECK(shown(*w, MainWindow::id_play));
    CHECK_FALSE(grayed(*w, Command::play_pause));
    CHECK(w->settings().film.has_value());
    CHECK(position(*w, MainWindow::id_time) == static_cast<int>(n - 1));
    same_as_png(*w, frame(n - 1));
    // The step keys land on the snapshot asked for.
    w->execute(Command::snapshot_first);
    CHECK(w->settings().point == 0u);
    CHECK(text_of(*w, MainWindow::id_time_text).starts_with("t = 0 s (1 of "));
    CHECK(position(*w, MainWindow::id_time) == 0);
    same_as_png(*w, frame(0));
    w->execute(Command::snapshot_previous);
    CHECK(w->settings().point == 0u);
    w->execute(Command::snapshot_next);
    w->execute(Command::snapshot_next);
    CHECK(w->settings().point == 2u);
    CHECK(SendMessageW(w->child(MainWindow::id_point), CB_GETCURSEL, 0, 0) == 2);  // the point list follows
    same_as_png(*w, frame(2));
    w->execute(Command::snapshot_previous);
    CHECK(w->settings().point == 1u);
    w->execute(Command::snapshot_last);
    CHECK(w->settings().point == n - 1);
    w->execute(Command::snapshot_next);
    CHECK(w->settings().point == n - 1);
    // The slider: its snapshot shown.
    slide(*w, MainWindow::id_time, 1);
    CHECK(w->settings().point == 1u);
    same_as_png(*w, frame(1));

    // Play from the first: five snapshots a second to the last, every frame shown its snapshot's, and the
    // message loop never held 100 ms (WG1).
    w->execute(Command::snapshot_first);
    const auto start = std::chrono::steady_clock::now();
    w->execute(Command::play_pause);
    CHECK(w->playing());
    CHECK(text_of(*w, MainWindow::id_play) == "Pa&use");
    std::vector<std::size_t> seen{0};
    double longest_ms = 0.0;
    MSG m;
    while (w->playing() && std::chrono::steady_clock::now() - start < std::chrono::seconds(20)) {
        if (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
            const auto t0 = std::chrono::steady_clock::now();
            TranslateMessage(&m);
            DispatchMessageW(&m);
            longest_ms = std::max(longest_ms,
                                  std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        } else {
            MsgWaitForMultipleObjects(0, nullptr, FALSE, 10, QS_ALLINPUT);
        }
        const std::size_t k = *w->settings().point;
        if (k != seen.back()) {
            // No snapshot before its time (five a second).
            const double at = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            CHECK(at >= static_cast<double>(k) / app::playback_fps - 0.05);
            seen.push_back(k);
            same_as_png(*w, frame(k));
        }
    }
    const double played_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    CAPTURE(seen, longest_ms, played_s);
    CHECK_FALSE(w->playing());
    CHECK(seen.back() == n - 1);
    CHECK(std::ranges::is_sorted(seen));
    CHECK(played_s >= static_cast<double>(n - 1) / app::playback_fps - 0.05);
    CHECK(longest_ms < app_test::wg1_limit_ms());
    CHECK(text_of(*w, MainWindow::id_play) == "Pla&y");
    // Play from the last starts again at the first; Pause holds the snapshot; the button and the slider pause too.
    w->execute(Command::play_pause);
    CHECK(w->playing());
    CHECK(w->settings().point == 0u);
    pump_for(0.5);
    w->execute(Command::play_pause);
    CHECK_FALSE(w->playing());
    const std::size_t held = *w->settings().point;
    CHECK(held >= 1u);
    pump_for(0.5);
    CHECK(w->settings().point == held);
    same_as_png(*w, frame(held));
    click(*w, MainWindow::id_play);
    CHECK(w->playing());
    slide(*w, MainWindow::id_time, 0);
    CHECK_FALSE(w->playing());
    CHECK(w->settings().point == 0u);
    // Another tab pauses it.
    w->execute(Command::play_pause);
    CHECK(w->playing());
    w->select_tab(MainWindow::tab_deck);
    pump_for(0.2);
    CHECK_FALSE(w->playing());
    w->select_tab(MainWindow::tab_view);

    // Only on the field view: not on the mesh view (the commands greyed).
    choose(*w, MainWindow::id_view, L"mesh");
    CHECK_FALSE(shown(*w, MainWindow::id_time));
    CHECK(grayed(*w, Command::play_pause));
    w->execute(Command::play_pause);
    CHECK_FALSE(w->playing());
    choose(*w, MainWindow::id_view, L"field");
    CHECK(shown(*w, MainWindow::id_time));
    // In the Simulation menu with their keys; the new controls are tab stops (YG5).
    for (const Command c : {Command::play_pause, Command::snapshot_first, Command::snapshot_previous,
                            Command::snapshot_next, Command::snapshot_last}) {
        const auto* info = app::find_command(c);
        REQUIRE(info);
        CHECK(info->menu == "&Simulation");
        CHECK(info->accelerator.key != 0);
    }
    for (const int id : {MainWindow::id_time, MainWindow::id_play}) {
        CHECK((GetWindowLongPtrW(w->child(id), GWL_STYLE) & WS_TABSTOP) != 0);
    }
    // The oxide hidden: the scale fitted again over the silicon, each frame as nitcad plot --hide oxide draws it.
    const fs::path hidden_frames = dir / "frames_no_oxide";
    const auto hidden_cli = app_test::nitcad("plot " + app_test::quoted(deck) +
                                             " --view filmstrip --field potential --every 1 --hide oxide --colormap viridis --size 480x360"
                                             " --adapter warp --frames " + app_test::quoted(hidden_frames) + " --out " +
                                             app_test::quoted(dir / "strip_no_oxide.png"));
    CAPTURE(hidden_cli.output);
    REQUIRE(hidden_cli.status == 0);
    press_space(*w, tree_item(*w, L"oxide"));
    REQUIRE(w->settings().hidden == std::vector<std::string>{"oxide"});
    w->execute(Command::snapshot_last);
    char name[32];
    std::snprintf(name, sizeof name, "frame_%03zu.png", n - 1);
    same_as_png(*w, hidden_frames / name);
}

TEST_CASE("app window 26c: the clipping box's sliders (YD6, YG4, YG5)", "[app]") {
    const fs::path dir = app_test::scratch("window_clip_sliders");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    choose(*w, MainWindow::id_view, L"solid");
    // Six sliders on the solid view of a 3D device, at the whole extent; no playback (a sweep).
    for (int id = MainWindow::id_clip_x0; id <= MainWindow::id_clip_z1; ++id) {
        CAPTURE(id);
        CHECK(shown(*w, id));
        CHECK(position(*w, id) == ((id - MainWindow::id_clip_x0) % 2 == 0 ? 0 : app::clip_steps));
        CHECK((GetWindowLongPtrW(w->child(id), GWL_STYLE) & WS_TABSTOP) != 0);
    }
    CHECK_FALSE(shown(*w, MainWindow::id_time));
    CHECK(grayed(*w, Command::play_pause));
    // The low x slider at a quarter: a box over the device from x = 0.5 um, checked and typed.
    slide(*w, MainWindow::id_clip_x0, 250);
    REQUIRE(w->settings().clip);
    CHECK(w->settings().clip->lo_um == std::array{0.5, 0.0, 0.0});
    CHECK(w->settings().clip->hi_um == std::array{2.0, 1.0, 0.8});
    CHECK(SendMessageW(w->child(MainWindow::id_clip), BM_GETCHECK, 0, 0) == BST_CHECKED);
    CHECK(text_of(*w, MainWindow::id_clip_box) == "0.5,0,0:2,1,0.8");
    same_as_plot(*w, deck, dir / "clip_x.png");
    // The high z slider below the low one: the low follows.
    slide(*w, MainWindow::id_clip_z0, 500);
    slide(*w, MainWindow::id_clip_z1, 300);
    CHECK(w->settings().clip->lo_um[2] == 0.24);
    CHECK(w->settings().clip->hi_um[2] == 0.24);
    CHECK(position(*w, MainWindow::id_clip_z0) == 300);
    // A typed box moves the sliders to its steps; the view as nitcad plot --clip draws it.
    apply_text(*w, MainWindow::id_clip_box, "0.2,0.1,0.1:1.5,0.75,0.6");
    REQUIRE(w->settings().clip);
    CHECK(position(*w, MainWindow::id_clip_x0) == 100);
    CHECK(position(*w, MainWindow::id_clip_x1) == 750);
    CHECK(position(*w, MainWindow::id_clip_y0) == 100);
    CHECK(position(*w, MainWindow::id_clip_y1) == 750);
    CHECK(position(*w, MainWindow::id_clip_z0) == 125);
    CHECK(position(*w, MainWindow::id_clip_z1) == 750);
    same_as_plot(*w, deck, dir / "clip_typed.png");
    // The field view keeps them; the mesh view has none.
    choose(*w, MainWindow::id_view, L"field");
    CHECK(shown(*w, MainWindow::id_clip_y1));
    slide(*w, MainWindow::id_clip_y1, 500);
    CHECK(w->settings().clip->hi_um[1] == 0.5);
    same_as_plot(*w, deck, dir / "clip_field.png");
    choose(*w, MainWindow::id_view, L"mesh");
    CHECK_FALSE(shown(*w, MainWindow::id_clip_x0));
    // At 200% the low and high sliders stay side by side in the inspector (YG5).
    choose(*w, MainWindow::id_view, L"solid");
    RECT r{};
    GetWindowRect(w->hwnd(), &r);
    RECT suggested{r.left, r.top, r.left + 2 * (r.right - r.left), r.top + 2 * (r.bottom - r.top)};
    SendMessageW(w->hwnd(), WM_DPICHANGED, MAKEWPARAM(192, 192), reinterpret_cast<LPARAM>(&suggested));
    RECT low{}, high{}, inspector{};
    GetWindowRect(w->child(MainWindow::id_clip_x0), &low);
    GetWindowRect(w->child(MainWindow::id_clip_x1), &high);
    GetWindowRect(w->child(MainWindow::id_inspector), &inspector);
    CHECK(low.right <= high.left);
    CHECK(low.top == high.top);
    CHECK(high.right <= inspector.right);
    CHECK(high.bottom - high.top >= 40);
}

TEST_CASE("app window: an adaptive run's passes listed, a pass's mesh shown and the results back (Unit 17, AD8)", "[app]") {
    const fs::path dir = app_test::scratch("window_adapt");
    const fs::path deck = app_test::write(dir, "a.deck", app_test::adapt_diode);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_results);
    CHECK(w->results_page().pass_count() == 0);  // not run: no passes
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace()->has_results(); }));
    REQUIRE(w->workspace()->adapt);
    const std::size_t passes = w->workspace()->adapt->passes.size();
    const std::size_t final_nodes = w->workspace()->device->mesh().node_count();
    REQUIRE(passes >= 2);
    CHECK(final_nodes > 42);

    // The list: a row per pass, the last the shown and final one.
    w->select_tab(MainWindow::tab_results);
    REQUIRE(w->results_page().pass_count() == static_cast<int>(passes));
    CHECK(w->results_page().pass_cell(0, 1) == L"start");
    CHECK(w->results_page().pass_cell(0, 2) == L"42");
    const int last = static_cast<int>(passes) - 1;
    CHECK(w->results_page().pass_cell(last, 7).starts_with(L"shown, final"));
    bool adaptation_line = false;
    for (int r = 0; r < w->results_page().record_count(); ++r) {
        if (w->results_page().record_cell(r, 0) == L"Adaptation") adaptation_line = true;
    }
    CHECK(adaptation_line);

    // Pass 0's mesh in the device view; the list stays, marking it.
    REQUIRE(w->results_show_pass(0));
    CHECK(w->tab() == MainWindow::tab_view);
    CHECK(w->settings().kind == app::ViewKind::mesh);
    REQUIRE(w->workspace()->pass == std::size_t{0});
    CHECK(w->workspace()->device->mesh().node_count() == 42);
    CHECK_FALSE(w->workspace()->has_results());
    REQUIRE(w->scene());
    CHECK(w->results_page().pass_count() == static_cast<int>(passes));
    CHECK(w->results_page().pass_cell(0, 7).starts_with(L"shown"));
    // Another pass, then the last: the results back on the final mesh.
    REQUIRE(w->results_show_pass(1));
    CHECK(w->workspace()->pass == std::size_t{1});
    REQUIRE(w->results_show_pass(passes - 1));
    CHECK_FALSE(w->workspace()->pass);
    CHECK(w->workspace()->has_results());
    CHECK(w->workspace()->device->mesh().node_count() == final_nodes);
    CHECK_FALSE(w->results_show_pass(passes));  // out of range
}

namespace {

// Unit 26d: whether the inspector has a button reading Apply (ZD1 removed it).
bool has_apply(MainWindow& w) {
    for (HWND c = GetWindow(GetParent(w.child(MainWindow::id_view)), GW_CHILD); c; c = GetWindow(c, GW_HWNDNEXT)) {
        wchar_t text[64]{};
        GetWindowTextW(c, text, 64);
        if (std::wstring_view(text) == L"&Apply" || std::wstring_view(text) == L"Apply") return true;
    }
    return false;
}

// ZG1 for one box: a valid text taken on Enter, another on focus loss; focus leaving it unchanged takes nothing;
// an invalid text, on either, leaves the setting, gives the message and shows the setting again. `setting` is
// the setting as the box shows it.
void check_box(MainWindow& w, int id, const std::string& on_enter, const std::string& on_leave, const std::string& invalid,
               const std::string& message, const std::function<std::string()>& setting) {
    CAPTURE(id);
    enter_box(w, id, on_enter);
    CHECK(setting() == on_enter);
    CHECK(text_of(w, id) == on_enter);
    leave_box(w, id, on_leave);
    CHECK(setting() == on_leave);
    CHECK(text_of(w, id) == on_leave);
    // Text put in by a program, not typed: focus leaving does not take it.
    HWND box = w.child(id);
    SetWindowTextW(box, std::wstring(on_enter.begin(), on_enter.end()).c_str());
    SendMessageW(box, EM_SETMODIFY, FALSE, 0);
    SendMessageW(GetParent(box), WM_COMMAND, MAKEWPARAM(id, EN_KILLFOCUS), reinterpret_cast<LPARAM>(box));
    CHECK(setting() == on_leave);
    for (const bool enter : {true, false}) {
        CAPTURE(enter);
        (enter ? enter_box : leave_box)(w, id, invalid);
        CHECK(w.status(0) == message);
        CHECK(setting() == on_leave);
        CHECK(text_of(w, id) == on_leave);
    }
}

std::string list(const std::vector<double>& v) {
    std::string t;
    for (const double x : v) t += (t.empty() ? "" : ",") + std::format("{}", x);
    return t;
}

const render::ColorScale& drawn_scale(MainWindow& w) {
    REQUIRE(w.scene());
    const auto& f = std::get<render::FieldScene>(*w.scene());
    REQUIRE(f.colorbar);
    return f.colorbar->scale;
}

bool checked(MainWindow& w, int id) { return SendMessageW(w.child(id), BM_GETCHECK, 0, 0) == BST_CHECKED; }

}  // namespace

TEST_CASE("app window 26d: the inspector's boxes on Enter and focus loss; the colour range (ZD1, ZD3, ZG1, ZG3)", "[app]") {
    const fs::path dir = app_test::scratch("window_boxes");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    auto w = window();
    CHECK_FALSE(has_apply(*w));
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    REQUIRE(w->settings().kind == app::ViewKind::field);
    const app::ViewSettings& s = w->settings();

    // ZG10: the device view and the curves view start with viridis, first in the list (ZD8).
    const auto colormap_shown = [&] {
        HWND list = w->child(MainWindow::id_colormap);
        wchar_t text[32]{};
        SendMessageW(list, CB_GETLBTEXT, static_cast<WPARAM>(SendMessageW(list, CB_GETCURSEL, 0, 0)), reinterpret_cast<LPARAM>(text));
        return std::wstring(text);
    };
    CHECK(s.colormap == "viridis");
    CHECK(colormap_shown() == L"viridis");
    CHECK(SendMessageW(w->child(MainWindow::id_colormap), CB_GETCURSEL, 0, 0) == 0);
    // ZG7: the metrics strip's maximum by the field's symbol, with its unit.
    {
        wchar_t metrics[512]{};
        GetWindowTextW(w->child(MainWindow::id_metrics), metrics, 512);
        const std::wstring_view m = metrics;
        const auto at = m.find(L"Max \x03C8 ");
        REQUIRE(at != std::wstring_view::npos);
        const auto unit = m.find(L" V", at + 6);
        CHECK(unit != std::wstring_view::npos);
        CHECK(m.find(L"Max potential") == std::wstring_view::npos);
    }
    w->select_tab(MainWindow::tab_curves);
    CHECK(w->settings().colormap == "viridis");
    CHECK(colormap_shown() == L"viridis");
    w->select_tab(MainWindow::tab_view);

    // ZG1: every box, each taking effect by itself.
    choose(*w, MainWindow::id_arrows, L"E");
    check_box(*w, MainWindow::id_arrow_count, "40", "60", "0", "Arrows are a count from 1 to 100000.",
              [&] { return std::to_string(s.arrows); });
    same_as_plot(*w, deck, dir / "arrows.png");
    choose(*w, MainWindow::id_lines, L"E");
    check_box(*w, MainWindow::id_line_count, "7", "9", "2.5", "Streamlines are a count from 1 to 10000.",
              [&] { return std::to_string(s.lines); });
    enter_box(*w, MainWindow::id_seeds, " x=1 ");
    CHECK(s.seeds == "x=1");
    leave_box(*w, MainWindow::id_seeds, "x=0.5");
    CHECK(s.seeds == "x=0.5");
    CHECK(s.lines == 9);  // a box takes effect by itself
    same_as_plot(*w, deck, dir / "lines.png");
    choose(*w, MainWindow::id_arrows, L"none");
    choose(*w, MainWindow::id_lines, L"none");
    // The clipping box, while checked; unchecked, none, its text kept.
    SendMessageW(w->child(MainWindow::id_clip), BM_SETCHECK, BST_CHECKED, 0);
    check_box(*w, MainWindow::id_clip_box, "0.2,0.1,0.1:1.5,0.75,0.6", "0,0,0:1,1,0.8", "1,2",
              "A clipping box is <x0>,<y0>,<z0>:<x1>,<y1>,<z1> in µm.", [&] {
                  if (!s.clip) return std::string("none");
                  return list({s.clip->lo_um.begin(), s.clip->lo_um.end()}) + ":" +
                         list({s.clip->hi_um.begin(), s.clip->hi_um.end()});
              });
    CHECK(position(*w, MainWindow::id_clip_x1) == 500);  // YD6: the box on the sliders
    same_as_plot(*w, deck, dir / "clip.png");
    click(*w, MainWindow::id_clip);
    CHECK_FALSE(s.clip);
    CHECK(text_of(*w, MainWindow::id_clip_box) == "0,0,0:1,1,0.8");
    click(*w, MainWindow::id_clip);
    REQUIRE(s.clip);
    CHECK(s.clip->hi_um[0] == 1.0);
    click(*w, MainWindow::id_clip);
    CHECK_FALSE(s.clip);
    // The slice view's plane, the iso view's levels, the exploded view's gap.
    choose(*w, MainWindow::id_view, L"slice");
    check_box(*w, MainWindow::id_slice_at, "0.3", "0.5", "0.3,0.4", "A slice plane is at one position in µm.", [&] {
        return s.slices.size() == 1 && s.slices[0].axis == 2 ? std::format("{}", s.slices[0].at_um) : std::string("?");
    });
    same_as_plot(*w, deck, dir / "slice.png");
    choose(*w, MainWindow::id_slice_axis, L"x");
    REQUIRE(s.slices.size() == 1);
    CHECK(s.slices[0].axis == 0);
    CHECK(s.slices[0].at_um == 0.5);
    choose(*w, MainWindow::id_view, L"iso");
    check_box(*w, MainWindow::id_levels, "0.2,0.4", "0.3", "0.2,x", "Iso levels are numbers separated by commas.",
              [&] { return list(s.levels); });
    choose(*w, MainWindow::id_view, L"exploded");
    check_box(*w, MainWindow::id_gap, "0.2", "0.1", "-1", "The gap is a length in µm, at least 0.",
              [&] { return s.gap_um ? std::format("{}", *s.gap_um) : std::string(); });
    // The cut line: taken on focus loss as on Enter (it is drawn in the Curves tab).
    choose(*w, MainWindow::id_view, L"field");
    leave_box(*w, MainWindow::id_cut, "0.5,0.2,0.4:1.5,0.2,0.4");
    REQUIRE(w->cut());
    CHECK(w->cut()->to_um == std::vector<double>{1.5, 0.2, 0.4});
    CHECK(w->tab() == MainWindow::tab_curves);
    leave_box(*w, MainWindow::id_cut, "0.5");
    CHECK(w->status(0).starts_with("A cut line is"));
    CHECK(text_of(*w, MainWindow::id_cut) == "0.5,0.2,0.4:1.5,0.2,0.4");
    leave_box(*w, MainWindow::id_cut, "");
    CHECK_FALSE(w->cut());
    w->select_tab(MainWindow::tab_view);

    // ZG3: the range on the field view; Auto shows the scale drawn.
    const app::ViewSettings& d = w->settings();
    REQUIRE(d.kind == app::ViewKind::field);
    choose(*w, MainWindow::id_field, L"potential");  // the arrows made it E's magnitude
    CHECK(shown(*w, MainWindow::id_range_min));
    CHECK(checked(*w, MainWindow::id_range_auto));
    CHECK_FALSE(d.range);
    CHECK(text_of(*w, MainWindow::id_range_min) == std::format("{:.6g}", drawn_scale(*w).lo()));
    CHECK(text_of(*w, MainWindow::id_range_max) == std::format("{:.6g}", drawn_scale(*w).hi()));
    // Typed: the colour bar's scale is the range, as nitcad plot --range draws it.
    type_box(*w, MainWindow::id_range_min, "-0.25");
    enter_box(*w, MainWindow::id_range_max, "1.25");
    REQUIRE(d.range);
    CHECK(*d.range == std::pair{-0.25, 1.25});
    CHECK(drawn_scale(*w).lo() == -0.25);
    CHECK(drawn_scale(*w).hi() == 1.25);
    CHECK_FALSE(checked(*w, MainWindow::id_range_auto));
    same_as_plot(*w, deck, dir / "range.png");
    // Focus leaving the minimum takes it, with the maximum as shown.
    leave_box(*w, MainWindow::id_range_min, "0");
    REQUIRE(d.range);
    CHECK(*d.range == std::pair{0.0, 1.25});
    // Refused: a minimum not below the maximum, a value not a number or not finite.
    for (const auto& [lo, hi] : {std::pair{"2", "1.25"}, std::pair{"1", "1"}, std::pair{"a", "1"}, std::pair{"0", "inf"}}) {
        CAPTURE(lo, hi);
        type_box(*w, MainWindow::id_range_min, lo);
        enter_box(*w, MainWindow::id_range_max, hi);
        CHECK(w->status(0) == "The range is a minimum and a larger maximum.");
        REQUIRE(d.range);
        CHECK(*d.range == std::pair{0.0, 1.25});
        CHECK(text_of(*w, MainWindow::id_range_min) == "0");
        CHECK(text_of(*w, MainWindow::id_range_max) == "1.25");
    }
    // Auto clears it; unchecked, the scale drawn is kept as the range.
    click(*w, MainWindow::id_range_auto);
    CHECK_FALSE(d.range);
    CHECK(checked(*w, MainWindow::id_range_auto));
    const double auto_lo = drawn_scale(*w).lo(), auto_hi = drawn_scale(*w).hi();
    CHECK(text_of(*w, MainWindow::id_range_max) == std::format("{:.6g}", auto_hi));
    click(*w, MainWindow::id_range_auto);
    REQUIRE(d.range);
    CHECK(*d.range == std::pair{auto_lo, auto_hi});
    // Both boxes emptied: Auto.
    type_box(*w, MainWindow::id_range_min, "");
    enter_box(*w, MainWindow::id_range_max, "");
    CHECK_FALSE(d.range);
    // A new field resets the range to Auto; on a log scale (the electron density's own) both are above 0.
    type_box(*w, MainWindow::id_range_min, "1");
    enter_box(*w, MainWindow::id_range_max, "2");
    REQUIRE(d.range);
    choose(*w, MainWindow::id_field, L"n");
    CHECK_FALSE(d.range);
    CHECK(checked(*w, MainWindow::id_range_auto));
    type_box(*w, MainWindow::id_range_min, "0");
    enter_box(*w, MainWindow::id_range_max, "1e18");
    CHECK(w->status(0) == "The range is a minimum and a larger maximum, both above 0 on a log scale.");
    CHECK_FALSE(d.range);
    type_box(*w, MainWindow::id_range_min, "1e10");
    enter_box(*w, MainWindow::id_range_max, "1e18");
    REQUIRE(d.range);
    CHECK(drawn_scale(*w).kind() == render::ScaleKind::log10);
    CHECK(drawn_scale(*w).lo() == 1e10);
    same_as_plot(*w, deck, dir / "range_log.png");
    // Shown on the views drawn on a colour scale of the field, not on the others nor in the Curves tab.
    choose(*w, MainWindow::id_view, L"solid");
    CHECK_FALSE(shown(*w, MainWindow::id_range_min));
    CHECK_FALSE(shown(*w, MainWindow::id_range_auto));
    for (const wchar_t* view : {L"slice", L"iso", L"exploded", L"field"}) {
        CAPTURE(view == std::wstring_view(L"slice") ? 0 : view == std::wstring_view(L"iso") ? 1 : 2);
        choose(*w, MainWindow::id_view, view);
        CHECK(shown(*w, MainWindow::id_range_min));
    }
    w->select_tab(MainWindow::tab_curves);
    CHECK_FALSE(shown(*w, MainWindow::id_range_min));
}

namespace {

// Unit 26d (ZG4): the files of a drop as the shell hands them over (DROPFILES with wide names); the window frees
// it (DragFinish).
HDROP make_drop(const std::vector<fs::path>& files) {
    std::wstring list;
    for (const fs::path& f : files) {
        list += f.wstring();
        list.push_back(L'\0');
    }
    list.push_back(L'\0');
    const SIZE_T bytes = sizeof(DROPFILES) + list.size() * sizeof(wchar_t);
    HGLOBAL g = GlobalAlloc(GHND, bytes);
    REQUIRE(g);
    auto* d = static_cast<DROPFILES*>(GlobalLock(g));
    d->pFiles = sizeof(DROPFILES);
    d->fWide = TRUE;
    std::memcpy(reinterpret_cast<char*>(d) + sizeof(DROPFILES), list.data(), list.size() * sizeof(wchar_t));
    GlobalUnlock(g);
    return static_cast<HDROP>(g);
}

void drop(MainWindow& w, const std::vector<fs::path>& files) {
    SendMessageW(w.hwnd(), WM_DROPFILES, reinterpret_cast<WPARAM>(make_drop(files)), 0);
}

// Answers the window's next message box with Cancel, from another thread; true once it has.
struct CancelNextBox {
    std::atomic<bool> answered{false};
    std::thread thread;
    CancelNextBox() {
        thread = std::thread([this] {
            const auto end = std::chrono::steady_clock::now() + std::chrono::seconds(20);
            while (std::chrono::steady_clock::now() < end) {
                if (HWND box = FindWindowW(L"#32770", L"NiTCAD")) {
                    answered = true;
                    PostMessageW(box, WM_COMMAND, IDCANCEL, 0);
                    return;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
            }
        });
    }
    ~CancelNextBox() { thread.join(); }
};

}  // namespace

TEST_CASE("app window 26d: files dropped on the window (ZD4, ZG4)", "[app]") {
    const fs::path dir = app_test::scratch("window_drop");
    const fs::path deck = app_test::write(dir, "diode.deck", app_test::diode_sweep());
    const fs::path other = app_test::write(dir, "other.deck", app_test::diode_sweep());
    auto w = window();
    // A deck: opened as Open deck opens it.
    drop(*w, {deck});
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    CHECK(w->document().path() == deck);
    // Refused: two files, a folder; and while a run goes on.
    drop(*w, {deck, other});
    CHECK(w->status(0) == "Drop one deck or result file at a time.");
    drop(*w, {dir});
    CHECK(w->status(0) == "A folder does not open: drop a deck or a result file.");
    CHECK(w->document().path() == deck);
    w->execute(Command::run);
    REQUIRE(w->runner().busy());
    drop(*w, {other});
    CHECK(w->status(0) == "A run or check is in progress: drop the file when it ends.");
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    CHECK(w->document().path() == deck);
    // A result file: with the device file of its name beside it.
    const fs::path result = dir / "diode.nitsweep";
    REQUIRE(fs::exists(result));
    REQUIRE(fs::exists(dir / "diode.nitdev"));
    drop(*w, {result});
    REQUIRE(w->workspace());
    CHECK(w->workspace()->source == result);
    CHECK(w->workspace()->has_results());
    // Unsaved edits: asked first; Cancel keeps the deck being edited.
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_deck);
    replace(*w, "[sweep]", "[sweep]\n# edited");
    REQUIRE(w->document().modified());
    {
        CancelNextBox cancel;
        drop(*w, {other});
        CHECK(cancel.answered);
    }
    CHECK(w->document().path() == deck);
    CHECK(w->document().modified());
}

TEST_CASE("app window 26d: the taskbar's progress and the title bar's theme (ZD5, ZD6, ZG5, ZG6)", "[app]") {
    const fs::path dir = app_test::scratch("window_shell");
    const fs::path deck = app_test::write(dir, "diode.deck", app_test::diode_sweep());
    auto w = window();
    // ZG6: the theme given to the title bar once at creation, as the user's app theme says; again per theme change.
    DWORD value = 1, size = sizeof(value);
    const bool read = RegGetValueW(HKEY_CURRENT_USER, L"Software\\Microsoft\\Windows\\CurrentVersion\\Themes\\Personalize",
                                   L"AppsUseLightTheme", RRF_RT_REG_DWORD, nullptr, &value, &size) == ERROR_SUCCESS;
    CHECK(w->dark_title() == app::dark_title_bar(read ? std::optional<unsigned long>(value) : std::nullopt));
    CHECK(w->title_theme_updates() == 1);
    SendMessageW(w->hwnd(), WM_SETTINGCHANGE, 0, reinterpret_cast<LPARAM>(L"ImmersiveColorSet"));
    CHECK(w->title_theme_updates() == 2);
    SendMessageW(w->hwnd(), WM_SETTINGCHANGE, 0, reinterpret_cast<LPARAM>(L"Environment"));
    SendMessageW(w->hwnd(), WM_SETTINGCHANGE, 0, 0);
    CHECK(w->title_theme_updates() == 2);
    // ZG5: a hidden window has no taskbar button; the progress is kept all the same and nothing fails.
    CHECK(w->taskbar_progress() == app::TaskbarProgress{});
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    CHECK((w->taskbar_progress().flag == app::TaskbarFlag::normal || w->taskbar_progress().flag == app::TaskbarFlag::indeterminate));
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    CHECK(w->taskbar_progress() == app::TaskbarProgress{});
}

TEST_CASE("app window 26f: the Enhanced 3D look check box (KD5, KG5)", "[app]") {
    const fs::path dir = app_test::scratch("window_enhanced_look");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    choose(*w, MainWindow::id_view, L"solid");
    const auto look = [&] {
        REQUIRE(w->scene());
        return std::get<render::FieldScene>(*w->scene()).look;
    };
    // On by default, on a 3D view; the command names no look and draws the app's image.
    CHECK(w->app_state().enhanced_look);
    CHECK(shown(*w, MainWindow::id_enhanced));
    CHECK(SendMessageW(w->child(MainWindow::id_enhanced), BM_GETCHECK, 0, 0) == BST_CHECKED);
    CHECK(w->settings().enhanced);
    CHECK(look() == render::Look::enhanced);
    same_as_plot(*w, deck, dir / "enhanced.png");
    // Off: classic, the command says so and still draws the app's image; remembered in the saved state.
    click(*w, MainWindow::id_enhanced);
    CHECK_FALSE(w->settings().enhanced);
    CHECK(look() == render::Look::classic);
    const auto args = app::plot_arguments(w->settings());
    CHECK(std::ranges::search(args, std::vector<std::string>{"--look", "classic"}).begin() != args.end());
    same_as_plot(*w, deck, dir / "classic.png");
    CHECK_FALSE(w->app_state().enhanced_look);
    CHECK_FALSE(app::parse_state(app::serialize(w->app_state())).enhanced_look);
    // Kept on the other device views; hidden on the curves.
    choose(*w, MainWindow::id_view, L"field");
    CHECK(SendMessageW(w->child(MainWindow::id_enhanced), BM_GETCHECK, 0, 0) == BST_UNCHECKED);
    CHECK(look() == render::Look::classic);
    w->select_tab(MainWindow::tab_curves);
    CHECK_FALSE(shown(*w, MainWindow::id_enhanced));
    // A 2D device: hidden, and its new views start with the remembered choice.
    const fs::path deck2 = app_test::write(dir, "generate2d.deck", app_test::generate2d);
    REQUIRE(w->open_deck(deck2));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_view);
    CHECK_FALSE(shown(*w, MainWindow::id_enhanced));
    CHECK_FALSE(w->settings().enhanced);
}

TEST_CASE("app window 26h: a click on the view cube snaps; a drag from it orbits; Probe still probes (JD3, JG6)",
          "[app]") {
    const fs::path dir = app_test::scratch("window_view_cube");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    choose(*w, MainWindow::id_view, L"solid");
    const int W = w->viewport().width(), H = w->viewport().height();
    const auto face = [&](std::string_view name) {
        REQUIRE(w->scene());
        const auto cube = render::view_cube(std::get<render::FieldScene>(*w->scene()), W, H);
        const auto it = std::ranges::find_if(cube.faces, [&](const render::CubeFace& f) { return f.label == name; });
        REQUIRE(it != cube.faces.end());
        float x = 0.0f, y = 0.0f;
        for (const auto& p : it->polygon) {
            x += 0.25f * p[0];
            y += 0.25f * p[1];
        }
        return MAKELPARAM(static_cast<int>(std::lround(x)), static_cast<int>(std::lround(y)));
    };
    const HWND view = w->viewport().hwnd();
    // Zoomed and panned first: the snap resets them, as Fit does.
    w->settings().zoom = 2.0f;
    w->settings().pan_x = 0.1f;
    // A click on Top: straight down, zoom 1, no pan; the scene's camera says so.
    const LPARAM top = face("Top");
    SendMessageW(view, WM_LBUTTONDOWN, MK_LBUTTON, top);
    SendMessageW(view, WM_LBUTTONUP, 0, top);
    CHECK(w->settings().azimuth == 0.0f);
    CHECK(w->settings().elevation == 90.0f);
    CHECK(w->settings().zoom == 1.0f);
    CHECK(w->settings().pan_x == 0.0f);
    CHECK(std::get<render::FieldScene>(*w->scene()).camera.elevation_deg == 90.0f);
    // From Top, the Top face is the one shown: a click on Front's place now is nothing of the cube's.
    // Back to the default camera, then a drag that starts on the cube orbits and does not snap.
    w->settings().azimuth.reset();
    w->settings().elevation.reset();
    w->refresh();
    const LPARAM front = face("Front");
    const float before = std::get<render::FieldScene>(*w->scene()).camera.azimuth_deg;
    SendMessageW(view, WM_LBUTTONDOWN, MK_LBUTTON, front);
    SendMessageW(view, WM_MOUSEMOVE, MK_LBUTTON, MAKELPARAM(static_cast<short>(LOWORD(front)) + 30, static_cast<short>(HIWORD(front))));
    SendMessageW(view, WM_LBUTTONUP, 0, MAKELPARAM(static_cast<short>(LOWORD(front)) + 30, static_cast<short>(HIWORD(front))));
    REQUIRE(w->settings().azimuth.has_value());
    CHECK(*w->settings().azimuth != before);
    CHECK(*w->settings().azimuth != 0.0f);
    CHECK(w->settings().elevation != 0.0f);
    // Probe mode (on by default): a click on the device probes; a click on the cube snaps, probe kept.
    w->settings().azimuth.reset();
    w->settings().elevation.reset();
    w->refresh();
    // The front-most node at the device's middle: the top surface's node nearest x and z mid.
    const auto& mesh = w->workspace()->device->mesh();
    std::size_t pick = 0;
    double best = 1e300;
    for (std::size_t n = 0; n < mesh.points().size(); ++n) {
        const auto& p = mesh.points()[n];
        const double d = std::abs(p[0] * 1e4 - 1.0) + std::abs(p[2] * 1e4 - 0.6) - p[1] * 1e4;
        if (d < best) {
            best = d;
            pick = n;
        }
    }
    const auto& node = mesh.points()[pick];
    click_at(*w, {float(node[0] * 1e4), float(node[1] * 1e4), float(node[2] * 1e4)});
    CHECK(w->probe());
    const auto probed = w->probe()->node;
    const LPARAM right = face("Right");
    SendMessageW(view, WM_LBUTTONDOWN, MK_LBUTTON, right);
    SendMessageW(view, WM_LBUTTONUP, 0, right);
    CHECK(w->settings().azimuth == 90.0f);
    CHECK(w->settings().elevation == 0.0f);
    REQUIRE(w->probe());
    CHECK(w->probe()->node == probed);
}

TEST_CASE("app window 26h: the named views from the menu and the numpad, exact, in the command (JD4, JG4)",
          "[app]") {
    const fs::path dir = app_test::scratch("window_named_views");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    choose(*w, MainWindow::id_view, L"solid");
    const auto named = [](std::string_view name) {
        return *std::ranges::find_if(render::named_views(), [&](const render::CameraPreset& p) { return p.name == name; });
    };
    const auto lands_on = [&](std::string_view name) {
        CAPTURE(name);
        const auto p = named(name);
        CHECK_FALSE(w->camera_moving());
        CHECK(w->settings().azimuth == p.azimuth_deg);
        CHECK(w->settings().elevation == p.elevation_deg);
        CHECK(w->settings().zoom == 1.0f);
        CHECK(w->settings().pan_x == 0.0f);
        CHECK(w->settings().pan_y == 0.0f);
        const auto& camera = std::get<render::FieldScene>(*w->scene()).camera;
        CHECK(camera.azimuth_deg == p.azimuth_deg);
        CHECK(camera.elevation_deg == p.elevation_deg);
        // The equivalent command names them.
        const auto args = app::plot_arguments(w->settings());
        CHECK(std::ranges::search(args, std::vector<std::string>{"--azimuth", std::format("{}", p.azimuth_deg)}).begin() !=
              args.end());
        CHECK(std::ranges::search(args, std::vector<std::string>{"--elevation", std::format("{}", p.elevation_deg)})
                  .begin() != args.end());
    };
    // View ▸ Look from: the seven, in order, with their numpad keys shown.
    HMENU bar = GetMenu(w->hwnd());
    REQUIRE(bar);
    HMENU view_menu = GetSubMenu(bar, GetMenuItemCount(bar) - 1);
    REQUIRE(view_menu);
    HMENU look = nullptr;
    for (int k = 0; k < GetMenuItemCount(view_menu); ++k) {
        wchar_t text[64]{};
        GetMenuStringW(view_menu, static_cast<UINT>(k), text, 64, MF_BYPOSITION);
        if (std::wstring_view(text) == L"Loo&k from") look = GetSubMenu(view_menu, k);
    }
    REQUIRE(look);
    REQUIRE(GetMenuItemCount(look) == 7);
    CHECK(GetMenuItemID(look, 0) == static_cast<UINT>(Command::view_front));
    CHECK(GetMenuItemID(look, 6) == static_cast<UINT>(Command::view_isometric));
    wchar_t first[64]{};
    GetMenuStringW(look, 0, first, 64, MF_BYPOSITION);
    CHECK(std::wstring_view(first) == L"&Front\tNum 1");
    // The View menu's items, through WM_COMMAND as the menu sends them; zoomed and panned first each time.
    const std::vector<std::pair<Command, std::string_view>> items{
        {Command::view_front, "Front"}, {Command::view_back, "Back"},     {Command::view_left, "Left"},
        {Command::view_right, "Right"}, {Command::view_top, "Top"},       {Command::view_bottom, "Bottom"},
        {Command::view_isometric, "Isometric"}};
    for (const auto& [command, name] : items) {
        w->settings().zoom = 3.0f;
        w->settings().pan_y = -0.2f;
        SendMessageW(w->hwnd(), WM_COMMAND, MAKEWPARAM(static_cast<int>(command), 0), 0);
        lands_on(name);
    }
    // The numpad keys, in the view.
    const HWND view = w->viewport().hwnd();
    const std::vector<std::pair<WPARAM, std::string_view>> keys{
        {VK_NUMPAD1, "Front"}, {VK_NUMPAD3, "Right"}, {VK_NUMPAD7, "Top"}, {VK_NUMPAD0, "Isometric"}};
    for (const auto& [key, name] : keys) {
        SendMessageW(view, WM_KEYDOWN, key, 0);
        lands_on(name);
    }
    // Top and Bottom draw (make_view's world-up switch): the app's image is the command's.
    w->execute(Command::view_top);
    same_as_plot(*w, deck, dir / "top.png");
    w->execute(Command::view_bottom);
    same_as_plot(*w, deck, dir / "bottom.png");
    // Both looks: classic too.
    w->settings().enhanced = false;
    w->refresh();
    w->execute(Command::view_right);
    lands_on("Right");
}

TEST_CASE("app window 26h: smooth camera moves animate, land exactly, and a drag stops them (JD5, JG5, JG8)",
          "[app]") {
    const fs::path dir = app_test::scratch("window_smooth_moves");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    choose(*w, MainWindow::id_view, L"solid");
    // Off (the tests' default): instant.
    w->execute(Command::view_top);
    CHECK_FALSE(w->camera_moving());
    CHECK(w->settings().elevation == 90.0f);

    w->set_animate(true);
    std::vector<double> frame_ms;
    // Follows a move to its end: the poses on the way, each frame's drawing time.
    const auto follow = [&] {
        std::vector<app::CameraPose> seen;
        REQUIRE(pump_until(
            [&] {
                seen.push_back(app::camera_pose(w->settings()));
                frame_ms.push_back(w->viewport().last_present_ms());
                return !w->camera_moving();
            },
            5.0));
        return seen;
    };
    // A named view: frames between the ends, then the preset exactly.
    const app::CameraPose from = app::camera_pose(w->settings());
    w->execute(Command::view_right);
    CHECK(w->camera_moving());
    const auto seen = follow();
    const app::CameraPose to{90.0f, 0.0f, 1.0f, 0.0f, 0.0f};
    CHECK(app::camera_pose(w->settings()) == to);
    CHECK(std::get<render::FieldScene>(*w->scene()).camera.azimuth_deg == 90.0f);
    CHECK(std::ranges::any_of(seen, [&](const app::CameraPose& p) { return !(p == from) && !(p == to); }));
    // Fit: the zoom and pan move, the angles stay; it lands exactly.
    w->settings().zoom = 4.0f;
    w->settings().pan_x = 0.3f;
    w->refresh();
    w->execute(Command::fit);
    CHECK(w->camera_moving());
    const auto fitted = follow();
    CHECK(app::camera_pose(w->settings()) == to);
    CHECK(std::ranges::any_of(fitted, [](const app::CameraPose& p) { return p.zoom > 1.0f && p.zoom < 4.0f; }));
    // A drag during a move ends it where it is: not at the target, and no more frames.
    w->execute(Command::view_front);
    REQUIRE(w->camera_moving());
    REQUIRE(pump_until([&] { return !(app::camera_pose(w->settings()) == to); }, 5.0));
    const HWND view = w->viewport().hwnd();
    const int cx = w->viewport().width() / 2, cy = w->viewport().height() / 2;
    SendMessageW(view, WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(cx, cy));
    SendMessageW(view, WM_MOUSEMOVE, MK_LBUTTON, MAKELPARAM(cx + 20, cy));
    SendMessageW(view, WM_LBUTTONUP, 0, MAKELPARAM(cx + 20, cy));
    CHECK_FALSE(w->camera_moving());
    const app::CameraPose stopped = app::camera_pose(w->settings());
    pump_until([] { return false; }, 0.4);
    CHECK(app::camera_pose(w->settings()) == stopped);
    CHECK_FALSE(stopped.azimuth_deg == 0.0f);
    // JG8: each frame drawn under 16 ms (WARP in the tests; the owner's hardware run is the gate).
    const double worst = frame_ms.empty() ? 0.0 : *std::ranges::max_element(frame_ms);
    WARN("JG8: a smooth move's slowest frame drew in " << worst << " ms (" << frame_ms.size() << " samples)");
}

TEST_CASE("app window 26h: the Axis box check box (JD1, JQ3, JG6)", "[app]") {
    const fs::path dir = app_test::scratch("window_axis_box");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    choose(*w, MainWindow::id_view, L"solid");
    const auto boxed = [&] {
        REQUIRE(w->scene());
        return std::get<render::FieldScene>(*w->scene()).show_axis_box;
    };
    // On by default on a 3D view; the command names no --axes and draws the app's image.
    CHECK(w->app_state().axis_box);
    CHECK(shown(*w, MainWindow::id_axis_box));
    CHECK(SendMessageW(w->child(MainWindow::id_axis_box), BM_GETCHECK, 0, 0) == BST_CHECKED);
    CHECK(boxed());
    same_as_plot(*w, deck, dir / "boxed.png");
    // Off: no box, --axes off in the command, the same image; remembered in the saved state.
    click(*w, MainWindow::id_axis_box);
    CHECK_FALSE(w->settings().axis_box);
    CHECK_FALSE(boxed());
    const auto args = app::plot_arguments(w->settings());
    CHECK(std::ranges::search(args, std::vector<std::string>{"--axes", "off"}).begin() != args.end());
    same_as_plot(*w, deck, dir / "unboxed.png");
    CHECK_FALSE(w->app_state().axis_box);
    CHECK_FALSE(app::parse_state(app::serialize(w->app_state())).axis_box);
    // Hidden on the curves and on a 2D device.
    w->select_tab(MainWindow::tab_curves);
    CHECK_FALSE(shown(*w, MainWindow::id_axis_box));
    const fs::path deck2 = app_test::write(dir, "generate2d.deck", app_test::generate2d);
    REQUIRE(w->open_deck(deck2));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_view);
    CHECK_FALSE(shown(*w, MainWindow::id_axis_box));
}


namespace {

// Unit 26i (NG3): a Viewport's host that turns its scene's camera as MainWindow does (the camera, then a
// camera-only present).
struct OrbitHost final : app::Viewport::Host {
    render::FieldScene* scene = nullptr;
    app::Viewport* view = nullptr;
    void viewport_orbit(float azimuth_deg, float elevation_deg) override {
        scene->camera.azimuth_deg += azimuth_deg;
        scene->camera.elevation_deg = std::clamp(scene->camera.elevation_deg + elevation_deg, -89.0f, 89.0f);
        view->present_camera();
    }
    void viewport_pan(float, float) override {}
    void viewport_zoom(float factor) override {
        scene->camera.zoom *= factor;
        view->present_camera();
    }
    void viewport_fit() override {}
    void viewport_click(float, float) override {}
    void viewport_hover(float, float) override {}
    void viewport_leave() override {}
    bool viewport_is_3d() const override { return true; }
    bool viewport_over_cube(float, float) const override { return false; }
    void viewport_named_view(int) override {}
};

}  // namespace

// Unit 26i (ND4, NG3): in a shown window, a burst of drag moves draws at most 3 frames, at the camera the
// same moves give one by one; a wheel burst likewise; each frame is the full draw's, bit for bit.
TEST_CASE("app viewport 26i: a burst of camera moves draws once, at the last camera (ND4, NG3)", "[app]") {
    // The scene: the 3D MOSFET's enhanced solid, built by a (hidden) window.
    const fs::path dir = app_test::scratch("viewport_coalesced_moves");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    std::optional<render::FieldScene> built;
    {
        auto w = window();
        REQUIRE(w->open_deck(deck));
        REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
        w->select_tab(MainWindow::tab_view);
        choose(*w, MainWindow::id_view, L"solid");
        REQUIRE(w->scene());
        REQUIRE(std::holds_alternative<render::FieldScene>(*w->scene()));
        built = std::get<render::FieldScene>(*w->scene());
    }
    auto renderer = render::Renderer::create({.adapter = render::Adapter::warp});
    REQUIRE(renderer);
    using Scene = std::variant<render::FieldScene, render::PlotScene>;
    Scene scene = *built;
    auto& field = std::get<render::FieldScene>(scene);
    HWND top = CreateWindowExW(0, L"STATIC", L"", WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN, 0, 0, 520, 420, nullptr, nullptr,
                               nullptr, nullptr);
    REQUIRE(top != nullptr);
    {
        OrbitHost host;
        app::Viewport view(top, host, 1);
        host.scene = &field;
        host.view = &view;
        MoveWindow(view.hwnd(), 0, 0, 480, 360, TRUE);
        ShowWindow(top, SW_SHOWNOACTIVATE);
        EnableWindow(top, FALSE);  // shown, but the desktop mouse does not reach it: only the posted moves
        view.show(&*renderer, &scene);
        pump_for(0.2);
        if (view.error()) SKIP("no swap chain here: " + *view.error());
        HWND v = view.hwnd();
        const render::Camera start = field.camera;
        constexpr int moves = 200;
        const auto drag = [&](bool burst) {
            SendMessageW(v, WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(100, 100));
            for (int k = 1; k <= moves; ++k) {
                const LPARAM at = MAKELPARAM(100 + k, 100 + k / 4);
                if (burst) {
                    PostMessageW(v, WM_MOUSEMOVE, MK_LBUTTON, at);
                } else {
                    SendMessageW(v, WM_MOUSEMOVE, MK_LBUTTON, at);
                    view.flush();
                }
            }
            pump_for(0.3);
            SendMessageW(v, WM_LBUTTONUP, 0, MAKELPARAM(100 + moves, 100 + moves / 4));
            pump_for(0.1);
        };
        std::size_t before = view.presents();
        drag(true);
        const std::size_t burst_frames = view.presents() - before;
        CAPTURE(burst_frames);
        CHECK(burst_frames >= 1);
        CHECK(burst_frames <= 3);
        const render::Camera after_burst = field.camera;
        CHECK(after_burst.azimuth_deg != start.azimuth_deg);
        CHECK(renderer->last_present_reused());
        // The last frame is the full draw's at that camera.
        const auto shown = renderer->last_image();
        REQUIRE(shown.has_value());
        const auto drawn = renderer->draw(field, view.width(), view.height());
        REQUIRE(drawn.has_value());
        CHECK(shown->rgba == drawn->rgba);
        // The same moves one at a time, each drawn: the same camera, to the bit.
        field.camera = start;
        before = view.presents();
        drag(false);
        CHECK(view.presents() - before >= static_cast<std::size_t>(moves) - 3);
        CHECK(field.camera.azimuth_deg == after_burst.azimuth_deg);
        CHECK(field.camera.elevation_deg == after_burst.elevation_deg);
        // A burst of wheel notches: few frames, every notch taken.
        const float zoom = field.camera.zoom;
        before = view.presents();
        for (int k = 0; k < 20; ++k) PostMessageW(v, WM_MOUSEWHEEL, MAKEWPARAM(0, WHEEL_DELTA), MAKELPARAM(200, 200));
        pump_for(0.3);
        CHECK(view.presents() - before <= 3);
        float expect = zoom;
        for (int k = 0; k < 20; ++k) expect *= std::pow(1.1f, 1.0f);
        CHECK(field.camera.zoom == expect);
        CHECK_FALSE(view.error());
        view.release();
    }
    DestroyWindow(top);
}

// Unit 26i (ND5, NG6): the metrics strip, computed once per workspace and field, reads as when it was
// computed each frame: the same after camera moves, after switching kinds and back, and anew for another
// workspace and another field.
TEST_CASE("app window 26i: the metrics strip is cached and unchanged (ND5, NG6)", "[app]") {
    const fs::path dir = app_test::scratch("window_cached_metrics");
    const fs::path deck3 = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    const fs::path deck2 = app_test::write(dir, "mosfet.deck", app_test::mosfet2d);
    auto w = window();
    // Without the "Drawn in" time, which changes every frame.
    const auto strip = [&] {
        std::wstring s = metrics(*w);
        if (const auto at = s.find(L"     Drawn in"); at != std::wstring::npos) s.erase(at);
        return s;
    };
    REQUIRE(w->open_deck(deck3));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    auto& host = static_cast<app::Viewport::Host&>(*w);
    // The window's frame is the scene's, drawn in full (a missed full present would show the old field).
    const auto frame_is_scene = [&] {
        const auto shown = w->renderer()->last_image();
        REQUIRE(shown.has_value());
        const auto drawn = w->renderer()->draw(std::get<render::FieldScene>(*w->scene()), shown->width, shown->height);
        REQUIRE(drawn.has_value());
        return shown->rgba == drawn->rgba;
    };
    const auto moved = [&] {
        for (int k = 0; k < 3; ++k) host.viewport_orbit(2.0f, 1.0f);
        return strip();
    };
    // The potential, then the electrons: each field's maximum, kept over camera moves.
    choose(*w, MainWindow::id_view, L"field");
    const std::wstring potential = strip();
    CAPTURE(narrow(potential));
    CHECK(potential.find(L"tetrahedra") != std::wstring::npos);
    CHECK(moved() == potential);
    w->settings().field = "n";
    w->refresh();
    const std::wstring electrons = strip();
    CAPTURE(narrow(electrons));
    CHECK(electrons != potential);  // the field's maximum
    CHECK(frame_is_scene());
    CHECK(moved() == electrons);
    CHECK(frame_is_scene());
    // Mesh quality, then back: its figures, the field's again.
    w->execute(Command::mesh_quality);
    const std::wstring quality = strip();
    CAPTURE(narrow(quality));
    CHECK(quality.find(L"Min angle") != std::wstring::npos);
    CHECK(moved() == quality);
    w->settings() = app::ViewSettings{};
    w->settings().kind = app::ViewKind::field;
    w->settings().field = "n";
    w->refresh();
    CHECK(strip() == electrons);
    // Another workspace: a 2D device, its own size and triangles.
    REQUIRE(w->open_deck(deck2));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->select_tab(MainWindow::tab_view);
    const std::wstring two = strip();
    CAPTURE(narrow(two));
    CHECK(two.find(L"triangles") != std::wstring::npos);
    CHECK(two.find(L"tetrahedra") == std::wstring::npos);
    CHECK(two.substr(0, two.find(L"Mesh")) != potential.substr(0, potential.find(L"Mesh")));
    w->execute(Command::mesh_quality);
    CHECK(strip().find(L"Min angle 26.6°, median 26.6°, 0 below 20°") != std::wstring::npos);
}

// A shown window (it gets a taskbar button, so the taskbar list) with a deck built and its 3D view drawn closes
// cleanly in a host that did not initialize COM itself, as this test executable (Unit 26i: it crashed releasing
// the taskbar list after the run thread's implicit apartment had gone).
TEST_CASE("app window: a shown window with a deck closes cleanly without the host's COM", "[app]") {
    const fs::path dir = app_test::scratch("window_shown_close");
    const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
    for (const bool input : {false, true}) {
        CAPTURE(input);
        app::WindowOptions o;
        o.adapter = render::Adapter::warp;
        o.show = true;
        o.input = input;
        o.quit_on_close = false;
        auto w = MainWindow::create(o);
        REQUIRE(w);
        REQUIRE(w->open_deck(deck));
        REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
        w->select_tab(MainWindow::tab_view);
        choose(*w, MainWindow::id_view, L"solid");
        pump_for(0.2);
        CHECK(w->scene());
        w.reset();
        pump_for(0.1);
    }
    SUCCEED("closed");
}

namespace {

// A move to the view's client pixel, as the mouse sends it (no button down).
void hover(MainWindow& w, float x, float y) {
    SendMessageW(w.viewport().hwnd(), WM_MOUSEMOVE, 0,
                 MAKELPARAM(static_cast<int>(std::lround(x)), static_cast<int>(std::lround(y))));
}

}  // namespace

TEST_CASE("app window 26o: the band surface in the inspector, stepped, hovered; its plane shared with Slice (BSD5-BSD7, BSC10)", "[app]") {
    using app::ViewKind;
    const fs::path dir = app_test::scratch("window_band_surface");
    {
        const fs::path deck = app_test::write(dir, "mosfet.deck", app_test::mosfet2d);
        auto w = window();
        REQUIRE(w->open_deck(deck));
        REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
        w->execute(Command::run);
        REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
        w->select_tab(MainWindow::tab_view);
        const std::array<int, 8> band_controls{MainWindow::id_band_ec, MainWindow::id_band_ev, MainWindow::id_band_efn,
                                               MainWindow::id_band_efp, MainWindow::id_by_level, MainWindow::id_band_fit,
                                               MainWindow::id_exaggeration, MainWindow::id_band_opacity};
        for (const int id : band_controls) CHECK_FALSE(shown(*w, id));
        choose(*w, MainWindow::id_view, L"band-surface");
        REQUIRE(w->settings().kind == ViewKind::band_surface);
        for (const int id : band_controls) {
            CAPTURE(id);
            CHECK(shown(*w, id));
        }
        CHECK(shown(*w, MainWindow::id_range_min));
        // Ev off, the heights doubled; a bad value is refused and the setting kept.
        click(*w, MainWindow::id_band_ev);
        CHECK(w->settings().band_levels == std::array<bool, 4>{true, false, false, false});
        apply_text(*w, MainWindow::id_exaggeration, "2");
        CHECK(w->settings().exaggeration == 2.0);
        apply_text(*w, MainWindow::id_exaggeration, "20");
        CHECK(w->settings().exaggeration == 2.0);
        apply_text(*w, MainWindow::id_band_opacity, "0.5");
        CHECK(w->settings().band_opacity == 0.5f);
        apply_text(*w, MainWindow::id_band_opacity, "1");
        // The first bias point: the scene is that point's.
        HWND list = w->child(MainWindow::id_point);
        SendMessageW(list, CB_SETCURSEL, 0, 0);
        SendMessageW(GetParent(list), WM_COMMAND, MAKEWPARAM(MainWindow::id_point, CBN_SELCHANGE), reinterpret_cast<LPARAM>(list));
        REQUIRE(w->settings().point == std::optional<std::size_t>{0});
        REQUIRE(w->scene() != nullptr);
        const auto& f = std::get<render::FieldScene>(*w->scene());
        auto expected = app::build_scene(*w->workspace(), w->settings());
        REQUIRE(expected);
        const auto& e = std::get<render::FieldScene>(*expected);
        REQUIRE(f.triangles.size() == e.triangles.size());
        for (std::size_t k = 0; k < f.triangles.size(); ++k) {
            CHECK(f.triangles[k].name == e.triangles[k].name);
            CHECK(f.triangles[k].vertices.size() == e.triangles[k].vertices.size());
        }
        // Hover over a substrate node's Ec: the readout's Ec is that node's.
        const auto& m = w->workspace()->device->mesh();
        std::size_t node = 0;
        double best = 1e300;
        for (std::size_t i = 0; i < m.node_count(); ++i) {
            const double d = std::hypot(m.points()[i][0] - 1.0e-4, m.points()[i][1] - 0.3e-4);
            if (d < best) best = d, node = i;
        }
        const double ec = w->workspace()->points[0].bands.conduction_eV[node];
        REQUIRE(std::isfinite(ec));
        const int W = w->viewport().width(), H = w->viewport().height();
        const double k = (*f.axes)[1].per_unit;
        const auto& p = m.points()[node];
        const render::ScreenPoint at = render::screen_position(
            f, W, H, {static_cast<float>(p[0] * 1e4), static_cast<float>(k * ec), static_cast<float>(p[1] * 1e4)});
        hover(*w, at.x, at.y);
        REQUIRE(w->readout_shown());
        const auto hit = app::band_pick(*w->workspace(), w->settings(), f, W, H, std::round(at.x), std::round(at.y));
        REQUIRE(hit);
        REQUIRE(hit->energy_eV[0]);
        CHECK(std::abs(*hit->energy_eV[0] - ec) < 2e-3);
        CHECK(w->readout_text() == app::band_hover_text(*hit, 2));
        CHECK(w->readout_text().starts_with("Ec = "));
        CHECK(w->status(3).starts_with("Ec = "));
        hover(*w, 1.0f, 1.0f);
        CHECK_FALSE(w->readout_shown());
        // Slice 2 (BSD8, BSD9): one range over the run until Fit this point; a cut line drawn on the surfaces.
        CHECK_FALSE(w->settings().band_fit_point);
        const auto lo_of = [&] { return std::get<render::FieldScene>(*w->scene()).colorbar->scale.lo(); };
        const auto run = render::band_energy_range(*w->workspace()->device, w->workspace()->points,
                                                   std::array<bool, 4>{true, false, false, false});
        REQUIRE(run);
        CHECK(lo_of() == run->first);
        click(*w, MainWindow::id_band_fit);
        CHECK(w->settings().band_fit_point);
        auto own = app::build_scene(*w->workspace(), w->settings());
        REQUIRE(own);
        CHECK(lo_of() == std::get<render::FieldScene>(*own).colorbar->scale.lo());
        CHECK(lo_of() != run->first);
        click(*w, MainWindow::id_band_fit);
        CHECK_FALSE(w->settings().band_fit_point);
        apply_cut(*w, "0.1,0.5:1.9,0.5");
        CHECK(w->tab() == MainWindow::tab_curves);  // the line's band diagram, as for every view
        w->select_tab(MainWindow::tab_view);
        REQUIRE(w->settings().kind == ViewKind::band_surface);
        const auto& lines = std::get<render::FieldScene>(*w->scene()).lines;
        CHECK(std::ranges::any_of(lines, [](const render::LineBatch& l) { return l.name == "band cut"; }));
        CHECK(std::ranges::none_of(lines, [](const render::LineBatch& l) { return l.name == "cut line"; }));
        apply_cut(*w, "");
    }
    {
        const fs::path deck = app_test::write(dir, "mosfet3d.deck", app_test::mosfet3d);
        auto w = window();
        REQUIRE(w->open_deck(deck));
        REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
        w->execute(Command::run);
        REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
        w->select_tab(MainWindow::tab_view);
        choose(*w, MainWindow::id_view, L"band-surface");
        REQUIRE(w->settings().kind == ViewKind::band_surface);
        // The plane through the device's centre normal to z, by default; the slice controls set it.
        REQUIRE(w->settings().band_plane);
        CHECK(w->settings().band_plane->axis == 2);
        CHECK(std::abs(w->settings().band_plane->at_um - 0.4) < 1e-9);
        CHECK(shown(*w, MainWindow::id_slice_axis));
        // On this view the slice controls are labelled for the plane they set.
        CHECK(text_of(*w, MainWindow::id_first_label + 5) == "Plane");
        CHECK(text_of(*w, MainWindow::id_first_label + 6) == "Plane at (?m)");
        choose(*w, MainWindow::id_slice_axis, L"x");
        apply_text(*w, MainWindow::id_slice_at, "1");
        REQUIRE(w->settings().band_plane);
        CHECK(w->settings().band_plane->axis == 0);
        CHECK(w->settings().band_plane->at_um == 1.0);
        REQUIRE(w->scene());
        REQUIRE(std::get<render::FieldScene>(*w->scene()).plane_inset);
        CHECK(std::get<render::FieldScene>(*w->scene()).plane_inset->axis == 0);
        // The Slice view shows the same plane; the Solid view outlines it.
        choose(*w, MainWindow::id_view, L"slice");
        CHECK(text_of(*w, MainWindow::id_first_label + 5) == "Slice plane");
        CHECK(text_of(*w, MainWindow::id_first_label + 6) == "Slice at (?m)");
        REQUIRE(w->settings().slices.size() == 1);
        CHECK(w->settings().slices.front().axis == 0);
        CHECK(w->settings().slices.front().at_um == 1.0);
        choose(*w, MainWindow::id_view, L"solid");
        REQUIRE(w->scene());
        const auto& solid = std::get<render::FieldScene>(*w->scene());
        CHECK(std::ranges::any_of(solid.lines, [](const render::LineBatch& l) { return l.name == "band plane"; }));
        // A plane picked in the Slice view is the band surface's when it is chosen.
        choose(*w, MainWindow::id_view, L"slice");
        choose(*w, MainWindow::id_slice_axis, L"y");
        apply_text(*w, MainWindow::id_slice_at, "0.5");
        choose(*w, MainWindow::id_view, L"band-surface");
        REQUIRE(w->settings().band_plane);
        CHECK(w->settings().band_plane->axis == 1);
        CHECK(w->settings().band_plane->at_um == 0.5);
    }
}

TEST_CASE("app window 26o: the band surface plays a sweep's bias points (BSD11, BSC13)", "[app]") {
    using app::ViewKind;
    const fs::path dir = app_test::scratch("window_band_playback");
    const fs::path deck = app_test::write(dir, "mosfet.deck", app_test::mosfet2d);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    REQUIRE(pump_until([&] { return w->workspace() != nullptr && !w->runner().busy(); }));
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }));
    w->select_tab(MainWindow::tab_view);
    const std::size_t n = w->workspace()->points.size();
    REQUIRE(n >= 3);
    // A sweep's solid view does not play; its band surface does.
    choose(*w, MainWindow::id_view, L"solid");
    CHECK_FALSE(shown(*w, MainWindow::id_play));
    CHECK(grayed(*w, Command::play_pause));
    choose(*w, MainWindow::id_view, L"band-surface");
    REQUIRE(w->settings().kind == ViewKind::band_surface);
    CHECK(shown(*w, MainWindow::id_time));
    CHECK(shown(*w, MainWindow::id_play));
    CHECK_FALSE(grayed(*w, Command::play_pause));
    // Each shown point's scene is build_scene's for that point.
    const auto same_scene = [&] {
        REQUIRE(w->scene() != nullptr);
        const auto& f = std::get<render::FieldScene>(*w->scene());
        auto expected = app::build_scene(*w->workspace(), w->settings());
        REQUIRE(expected);
        const auto& e = std::get<render::FieldScene>(*expected);
        REQUIRE(f.triangles.size() == e.triangles.size());
        for (std::size_t t = 0; t < f.triangles.size(); ++t) {
            REQUIRE(f.triangles[t].vertices.size() == e.triangles[t].vertices.size());
            for (std::size_t i = 0; i < f.triangles[t].vertices.size(); ++i) {
                if (f.triangles[t].vertices[i].position != e.triangles[t].vertices[i].position) {
                    FAIL("a vertex differs from build_scene's");
                }
            }
        }
    };
    // The step keys land on the point asked for; the label reads "(k of n)".
    w->execute(Command::snapshot_first);
    CHECK(w->settings().point == 0u);
    CHECK(position(*w, MainWindow::id_time) == 0);
    CHECK(text_of(*w, MainWindow::id_time_text).ends_with(std::format("(1 of {})", n)));
    CHECK(text_of(*w, MainWindow::id_time_text) == app::playback_label(*w->workspace(), 0));
    same_scene();
    w->execute(Command::snapshot_next);
    w->execute(Command::snapshot_next);
    CHECK(w->settings().point == 2u);
    CHECK(SendMessageW(w->child(MainWindow::id_point), CB_GETCURSEL, 0, 0) == 2);
    same_scene();
    w->execute(Command::snapshot_last);
    CHECK(w->settings().point == n - 1);

    // Play from the first: every point in order, none before its time, each its own scene; the loop never held long.
    w->execute(Command::snapshot_first);
    const auto start = std::chrono::steady_clock::now();
    w->execute(Command::play_pause);
    CHECK(w->playing());
    CHECK(text_of(*w, MainWindow::id_play) == "Pa&use");
    std::vector<std::size_t> seen{0};
    double longest_ms = 0.0;
    MSG m;
    while (w->playing() && std::chrono::steady_clock::now() - start < std::chrono::seconds(30)) {
        if (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
            const auto t0 = std::chrono::steady_clock::now();
            TranslateMessage(&m);
            DispatchMessageW(&m);
            longest_ms = std::max(longest_ms,
                                  std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        } else {
            MsgWaitForMultipleObjects(0, nullptr, FALSE, 10, QS_ALLINPUT);
        }
        const std::size_t k = *w->settings().point;
        if (k != seen.back()) {
            const double at = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            CHECK(at >= static_cast<double>(k) / app::playback_fps - 0.05);
            seen.push_back(k);
            same_scene();
        }
    }
    CAPTURE(seen, longest_ms);
    CHECK_FALSE(w->playing());
    CHECK(seen.back() == n - 1);
    CHECK(std::ranges::is_sorted(seen));
    // Release only (as KG6): the unoptimised Debug build records the time (284 ms on CI, 9a85819).
#ifdef NDEBUG
    CHECK(longest_ms < app_test::wg1_limit_ms());
#endif
    CHECK(text_of(*w, MainWindow::id_play) == "Pla&y");
    // Pause holds the point shown.
    w->execute(Command::snapshot_first);
    w->execute(Command::play_pause);
    CHECK(w->playing());
    w->execute(Command::play_pause);
    CHECK_FALSE(w->playing());
    const auto held = w->settings().point;
    (void)pump_until([] { return false; }, 0.5);
    CHECK(w->settings().point == held);
}
