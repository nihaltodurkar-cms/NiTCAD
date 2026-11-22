// Unit 26n slice 0 (DT1, DT2): the dimensions panel and the Geometry tab's dimension lines, by window messages only.
#include <catch2/catch_test_macros.hpp>

#include "dimensions_panel.hpp"

#include <windows.h>

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
