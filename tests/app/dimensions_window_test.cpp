// Unit 26n slice 0 (DT1, DT2): the dimensions panel and the Geometry tab's dimension lines, by window messages only.
#include <catch2/catch_test_macros.hpp>

#include "dimensions_panel.hpp"

#include <windows.h>

#include <algorithm>
#include <string>
#include <vector>

using namespace NiTCAD::app;

namespace {

const std::string mos2d = R"([parameters]
# group: Gate
Lg = 20nm  # Gate length
Lsd = 50nm
x_d = $Lsd + $Lg
Nsd = 1e20  # S/D doping

[mesh]
type = tensor
x = linspace(0, $Lsd, 6) + linspace($Lsd, $x_d, 5) + linspace($x_d, $x_d + $Lsd, 6)
y = linspace(0, 40nm, 9) + linspace(40nm, 42nm, 3)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0 $x_d + $Lsd  40nm 42nm

[doping]
profile = acceptor constant 1e16
profile = donor constant $Nsd x=0:$Lsd
profile = donor constant $Nsd x=$x_d:$x_d + $Lsd

[contact gate]
kind = electrode
patch = y_max
box = $Lsd $x_d  42nm 42nm
electrode = metal
work_function = 4.1

[equilibrium]
)";

struct FakeHost : DimensionsPanel::Host {
    std::string deck = mos2d;
    std::vector<std::string> edits;
    std::vector<std::string> selected;
    std::string dimensions_deck() override { return deck; }
    void dimensions_edited(const std::string& text) override {
        edits.push_back(text);
        deck = text;
    }
    void dimensions_selected(const std::string& name) override { selected.push_back(name); }
};

struct Parent {
    HWND hwnd = CreateWindowExW(0, L"STATIC", L"", WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN, 0, 0, 400, 600, nullptr, nullptr,
                                GetModuleHandleW(nullptr), nullptr);
    ~Parent() { DestroyWindow(hwnd); }
};

std::size_t changed_lines(const std::string& a, const std::string& b) {
    std::size_t n = 0, i = 0, j = 0;
    while (i <= a.size() && j <= b.size()) {
        const auto ea = std::min(a.find('\n', i), a.size()), eb = std::min(b.find('\n', j), b.size());
        if (a.compare(i, ea - i, b, j, eb - j) != 0) ++n;
        i = ea + 1;
        j = eb + 1;
    }
    return n;
}

}  // namespace

TEST_CASE("dimensions panel: rows from the deck; Enter writes one line; a refusal is shown, not written", "[app][window][dimensions]") {
    Parent p;
    FakeHost host;
    DimensionsPanel panel(p.hwnd, host, 4000);
    panel.set_dpi(96, static_cast<HFONT>(GetStockObject(DEFAULT_GUI_FONT)), static_cast<HFONT>(GetStockObject(DEFAULT_GUI_FONT)));
    panel.reload();
    REQUIRE(panel.rows() == 4);
    CHECK(panel.preferred_height(300) > 4 * 20);

    SECTION("a typed value, by the Enter key in the row's edit") {
        HWND edit = panel.row_edit(0);
        REQUIRE(edit);
        SetWindowTextW(edit, L"30nm");
        SendMessageW(edit, WM_KEYDOWN, VK_RETURN, 0);
        REQUIRE(host.edits.size() == 1);
        CHECK(changed_lines(mos2d, host.edits[0]) == 1);
        CHECK(host.edits[0].find("Lg = 30nm  # Gate length\n") != std::string::npos);
        CHECK(panel.row_error(0).empty());
    }
    SECTION("the edit asks the dialog manager for Enter") {
        CHECK((SendMessageW(panel.row_edit(0), WM_GETDLGCODE, VK_RETURN, 0) & DLGC_WANTALLKEYS) != 0);
    }
    SECTION("a refused value") {
        CHECK_FALSE(panel.commit(0, "abc"));
        CHECK(host.edits.empty());
        CHECK(panel.row_error(0).find("Lg") != std::string::npos);
        // The next good value clears it.
        CHECK(panel.commit(0, "25nm"));
        CHECK(panel.row_error(0).empty());
    }
    SECTION("the wheel over a number row") {
        CHECK(panel.wheel(3, 1));
        REQUIRE(host.edits.size() == 1);
        CHECK(host.edits[0].find("Nsd = 1.2589254117941") != std::string::npos);
        // A derived row is not stepped.
        CHECK_FALSE(panel.wheel(2, 1));
        CHECK(host.edits.size() == 1);
    }
    SECTION("focusing a row selects its dimension") {
        SendMessageW(panel.hwnd(), WM_COMMAND, MAKEWPARAM(GetDlgCtrlID(panel.row_edit(1)), EN_SETFOCUS),
                     reinterpret_cast<LPARAM>(panel.row_edit(1)));
        REQUIRE_FALSE(host.selected.empty());
        CHECK(host.selected.back() == "Lsd");
        panel.select("x_d");
        CHECK(panel.selected() == "x_d");
    }
    SECTION("a deck that stops parsing keeps the rows, disabled, and writes nothing") {
        host.deck = mos2d + "\n[broken\n";
        panel.reload();
        CHECK(panel.rows() == 4);
        CHECK_FALSE(IsWindowEnabled(panel.row_edit(0)));
        CHECK_FALSE(panel.message().empty());
        CHECK_FALSE(panel.commit(0, "30nm"));
        CHECK(host.edits.empty());
        host.deck = mos2d;
        panel.reload();
        CHECK(IsWindowEnabled(panel.row_edit(0)));
        CHECK(panel.message().empty());
    }
    SECTION("a deck without [parameters]: no rows, the hint") {
        host.deck = "[mesh]\ntype = tensor\nx = linspace(0, 1um, 11)\n\n[region si]\nmaterial = silicon\n";
        panel.reload();
        CHECK(panel.rows() == 0);
        CHECK(panel.message().find("No named dimensions") != std::string::npos);
    }
}

// --- Task 5: the editors host the panel and draw the dimension lines --------------------------------------------------

#include "main_window.hpp"
#include "fixtures.hpp"

#include "NiTCAD/app/commands.hpp"
#include "NiTCAD/app/geometry.hpp"
#include "NiTCAD/app/solid_view.hpp"

#include <cmath>
#include <fstream>

namespace {

void pump() {
    MSG m;
    while (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
        TranslateMessage(&m);
        DispatchMessageW(&m);
    }
}

std::unique_ptr<MainWindow> window() {
    WindowOptions o;
    o.adapter = NiTCAD::render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    return w;
}

std::filesystem::path write_deck(const std::string& name, const std::string& text) {
    const auto dir = app_test::scratch("dimensions_window_" + name);
    const auto p = dir / (name + ".deck");
    std::ofstream(p, std::ios::binary) << text;
    return p;
}

const std::string tied3d = R"([parameters]
x0 = 30nm
Lg = 20nm  # Gate length

[mesh]
type = generate3d
layered = on
domain = 0 100nm 0 60nm 0 14nm
patch body = 0 0 0  100nm 60nm 0
max_size = 0.05um
interface_resolution = 2nm

[region si]
material = silicon

[region gate]
material = sio2
box = $x0 $x0 + $Lg 40nm 50nm 0 14nm

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";

}  // namespace

TEST_CASE("dimensions in the 2D editor: the panel, the lines over their faces, a drag of one line, one undo", "[app][window][dimensions]") {
    auto w = window();
    REQUIRE(w->open_deck(write_deck("mos2d", mos2d)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE_FALSE(w->solid_shown());
    GeometryPage& page = w->geometry_page();
    INFO(page.geometry().why);
    REQUIRE(page.geometry().editable());
    CHECK(page.view_only());  // its shapes hold $references: drawn from the substituted text
    CHECK(page.dimensions_panel().rows() == 4);
    REQUIRE(page.render() >= 0.0);
    const auto& marks = page.dimension_marks();
    REQUIRE(marks.size() == 3);
    const auto lines = page.dimension_lines();
    REQUIRE(lines.size() == marks.size());
    // Lg: along the gate contact's x extent, 12 px above its y.
    const auto& v = page.view();
    const double lift = v.px(12.0);
    CHECK(std::abs(lines[0][0] - v.sx(500)) <= 1.0);
    CHECK(std::abs(lines[0][2] - v.sx(700)) <= 1.0);
    CHECK(std::abs(lines[0][1] - (v.sy(420) - lift)) <= 1.0);
    CHECK(std::abs(lines[0][3] - lines[0][1]) <= 0.01);

    // A click on the line selects its row; a drag along x by some pixels writes Lg + those quanta.
    const int mx = static_cast<int>(std::lround(0.5 * (lines[0][0] + lines[0][2]))), my = static_cast<int>(std::lround(lines[0][1]));
    const std::string before = w->deck().text();
    SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(mx, my));
    CHECK(page.dimensions_panel().selected() == "Lg");
    SendMessageW(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, MAKELPARAM(mx + 3, my));
    SendMessageW(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, MAKELPARAM(mx + 6, my));
    SendMessageW(page.canvas(), WM_LBUTTONUP, 0, MAKELPARAM(mx + 6, my));
    const geo::Coord dq = v.model(mx + 6, my).x - v.model(mx, my).x;
    REQUIRE(dq > 0);
    const std::string expected = "Lg = " + std::to_string((200 + dq) / 10) +
                                 ((200 + dq) % 10 ? "." + std::to_string((200 + dq) % 10) : std::string()) + "nm  # Gate length\n";
    CAPTURE(dq, expected);
    CHECK(w->deck().text().find(expected) != std::string::npos);
    CHECK(changed_lines(before, w->deck().text()) == 1);
    // One undo restores the deck byte for byte.
    w->execute(Command::undo);
    pump();
    CHECK(w->deck().text() == before);
}

TEST_CASE("dimensions in the 3D modeller: the panel replaces the parameters box; the lines and labels are in the scene", "[app][window][dimensions]") {
    auto w = window();
    REQUIRE(w->open_deck(write_deck("tied3d", tied3d)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    CHECK(FindWindowExW(page.inspector(), nullptr, L"BUTTON", L"Apply parameters") == nullptr);
    CHECK(page.dimensions_panel().rows() == 2);
    REQUIRE(page.render() >= 0.0);
    const auto& s = page.scene();
    const auto batch = std::ranges::find(s.lines, std::string("dimensions"), &NiTCAD::render::LineBatch::name);
    REQUIRE(batch != s.lines.end());
    CHECK(batch->on_top);
    REQUIRE(batch->points.size() >= 2);
    CHECK(std::abs(batch->points[0][0] - 0.03f) < 1e-6f);
    CHECK(std::abs(batch->points[1][0] - 0.05f) < 1e-6f);
    CHECK(std::abs(batch->points[0][1] - 0.04f) < 1e-6f);
    REQUIRE(s.labels.size() == 1);  // x0 has no span (its gate box ends both move)
    CHECK(s.labels[0].text == "Gate length = 20 nm");
    // A click on the line selects its row.
    const solid::V3 mid{0.04, 0.04, 0.0};
    const auto p = solid::screen_of(page.model(), page.view(), mid);
    SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(static_cast<int>(std::lround(p[0])), static_cast<int>(std::lround(p[1]))));
    SendMessageW(page.canvas(), WM_LBUTTONUP, 0, MAKELPARAM(static_cast<int>(std::lround(p[0])), static_cast<int>(std::lround(p[1]))));
    CHECK(page.dimensions_panel().selected() == "Lg");
}

// --- Task 6: the Geometry tab for decks neither editor draws ---------------------------------------------------------

#include <chrono>
#include <cstdio>
#include <functional>

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
            MsgWaitForMultipleObjects(0, nullptr, FALSE, 5, QS_ALLINPUT);
        }
    }
    return true;
}

bool shown(HWND h) { return (GetWindowLongPtrW(h, GWL_STYLE) & WS_VISIBLE) != 0; }

const std::string pn1d_deck = R"([parameters]
Ln = 0.8um  # n side
Lp = 1.2um  # p side

[mesh]
type = tensor
x = linspace(0, $Ln, 81) + linspace($Ln, $Ln + $Lp, 121)

[region silicon]
material = silicon

[doping]
profile = donor constant 1e18 x=0:$Ln
profile = acceptor constant 1e16

[contact cathode]
kind = ohmic
patch = x_min

[contact anode]
kind = ohmic
patch = x_max

[equilibrium]
)";

const std::string mos3d_deck = R"([parameters]
tox = 2nm  # Oxide
L = 100nm

[mesh]
type = tensor
x = linspace(0, $L, 11)
y = linspace(0, 40nm, 9) + linspace(40nm, 40nm + $tox, 3)
z = linspace(0, 30nm, 4)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0 $L  40nm 40nm + $tox  0 30nm

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = y_min

[equilibrium]
)";

const NiTCAD::render::FieldScene* field_scene(MainWindow& w) {
    const auto* s = w.scene();
    return s ? std::get_if<NiTCAD::render::FieldScene>(s) : nullptr;
}

bool has_label(MainWindow& w, const std::string& text) {
    const auto* f = field_scene(w);
    return f && std::ranges::any_of(f->labels, [&](const auto& l) { return l.text == text; });
}

double median(std::vector<double> v) {
    std::ranges::sort(v);
    return v[v.size() / 2];
}

}  // namespace

TEST_CASE("dimensions on the Geometry tab for a 1D deck: the preview in the canvas's place, the panel beside it", "[app][window][dimensions]") {
    auto w = window();
    REQUIRE(w->open_deck(write_deck("pn1d", pn1d_deck)));
    w->execute(Command::tab_geometry);
    REQUIRE(pump_until([&] { return w->preview() && w->scene() != nullptr; }));
    CHECK(shown(w->viewport().hwnd()));
    CHECK_FALSE(shown(w->geometry_page().hwnd()));
    CHECK(shown(w->geometry_page().inspector()));
    CHECK(w->geometry_page().dimensions_panel().rows() == 2);
}

TEST_CASE("dimensions on the Geometry tab for a 3D tensor deck: the lines and labels in the preview; an edit redraws it (DC9)", "[app][window][dimensions]") {
    auto w = window();
    REQUIRE(w->open_deck(write_deck("mos3d", mos3d_deck)));
    w->execute(Command::tab_geometry);
    REQUIRE(pump_until([&] { return w->preview() && has_label(*w, "Oxide = 2 nm"); }));
    CHECK(shown(w->viewport().hwnd()));
    CHECK_FALSE(shown(w->geometry_page().hwnd()));
    const auto* f = field_scene(*w);
    REQUIRE(f);
    CHECK(std::ranges::any_of(f->lines, [](const auto& b) { return b.name == "dimensions" && b.on_top; }));
    DimensionsPanel& panel = w->geometry_page().dimensions_panel();
    REQUIRE(panel.rows() == 2);
    std::vector<double> ms;
    for (int k = 0; k < 6; ++k) {
        const std::string tox = std::to_string(3 + k % 2) + "nm";
        const std::string label = "Oxide = " + std::to_string(3 + k % 2) + " nm";
        const auto t0 = std::chrono::steady_clock::now();
        REQUIRE(panel.commit(0, tox));
        REQUIRE(pump_until([&] { return has_label(*w, label); }, 30.0));
        ms.push_back(1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    std::printf("DC9: 3D tensor panel edit to refreshed preview: median %.0f ms\n", median(ms));
    CHECK(median(ms) <= 1000.0);
}

TEST_CASE("dimensions in the 2D editor: a panel edit to the redrawn canvas (DC9)", "[app][window][dimensions]") {
    auto w = window();
    REQUIRE(w->open_deck(write_deck("mos2d_speed", mos2d)));
    pump();
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    REQUIRE(page.view_only());
    std::vector<double> ms;
    for (int k = 0; k < 10; ++k) {
        const auto t0 = std::chrono::steady_clock::now();
        REQUIRE(page.dimensions_panel().commit(0, std::to_string(21 + k % 2) + "nm"));
        REQUIRE(page.render() >= 0.0);
        ms.push_back(1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    std::printf("DC9: 2D panel edit to redrawn canvas: median %.0f ms\n", median(ms));
    CHECK(median(ms) <= 300.0);
}

TEST_CASE("dimensions: a deck that does not parse keeps the editor's message on the Geometry tab, not a stale preview", "[app][window][dimensions]") {
    auto w = window();
    REQUIRE(w->open_deck(write_deck("mos2d_broken", mos2d)));
    pump();
    const std::string text = w->deck().text();
    const auto at = text.find("[mesh]");
    REQUIRE(at != std::string::npos);
    SendMessageW(w->deck().editor(), EM_SETSEL, at, at + 6);
    SendMessageW(w->deck().editor(), EM_REPLACESEL, TRUE, reinterpret_cast<LPARAM>(L"[mesh"));
    pump();
    w->execute(Command::tab_geometry);
    pump();
    CHECK(w->geometry_page().message().starts_with("The deck does not parse"));
    CHECK(shown(w->geometry_page().hwnd()));
    CHECK_FALSE(shown(w->viewport().hwnd()));
}
