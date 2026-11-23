// Unit 26n slice 2 (DT5): Unlock shapes on the dimensions panel, by window messages: a 3D tensor deck shown as the
// preview becomes a layered deck in the 3D modeller, and the modeller's Undo takes it back; a refusal stays on the
// panel with the deck unchanged; an extrude3d deck converts from the modeller's panel.
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "fixtures.hpp"
#include "NiTCAD/app/commands.hpp"
#include "NiTCAD/app/unlock.hpp"

#include <windows.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>

using namespace NiTCAD::app;
namespace fs = std::filesystem;

namespace {

bool pump_until(const std::function<bool()>& done, double seconds = 60.0) {
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

std::unique_ptr<MainWindow> window() {
    WindowOptions o;
    o.adapter = NiTCAD::render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    return w;
}

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

void click(const DimensionsPanel& p) {
    HWND b = p.unlock_button();
    SendMessageW(p.hwnd(), WM_COMMAND, MAKEWPARAM(GetDlgCtrlID(b), BN_CLICKED), reinterpret_cast<LPARAM>(b));
}

const std::string tensor_deck = R"([parameters]
L = 1um  # Length
H = 0.6um
W = 0.4um
t = 0.2um  # Oxide
y_ox = $H - $t

[mesh]
type = tensor
x = linspace(0, $L, 11)
y = linspace(0, $y_ox, 5) + linspace($y_ox, $H, 3)
z = 0, 0.2um, $W

[region si]
material = silicon

[region ox]
material = sio2
box = 0 $L  $y_ox $H  0 $W

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = y_min

[contact gate]
kind = electrode
patch = y_max
box = 0.2um 0.8um  $H $H  0 $W
electrode = metal
work_function = 4.1

[equilibrium]
)";

}  // namespace

TEST_CASE("unlock window: Unlock shapes takes a 3D tensor deck to the modeller; Undo brings it back (DT5)", "[app][window][unlock]") {
    const fs::path dir = app_test::scratch("unlock_window");
    const fs::path deck = app_test::write(dir, "box.deck", tensor_deck);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    w->execute(Command::tab_geometry);
    pump_until([&] { return w->preview() && w->scene() != nullptr; });
    INFO("preview " << w->preview() << " scene " << (w->scene() != nullptr) << " solid " << w->solid_shown() << " msg " << w->geometry_page().message());
    REQUIRE((w->preview() && w->scene() != nullptr));
    CHECK_FALSE(w->solid_shown());
    DimensionsPanel& panel = w->geometry_page().dimensions_panel();
    REQUIRE(panel.unlock_button() != nullptr);
    CHECK(panel.unlocks());
    click(panel);
    CHECK(panel.unlock_error().empty());
    const auto want = unlock::unlock(tensor_deck);
    REQUIRE(want);
    CHECK(w->deck().text() == want->text);
    CHECK(w->solid_shown());
    CHECK(w->solid_page().model().mesh == solid::MeshKind::layered);
    CHECK(w->solid_page().model().regions.size() == 2);
    CHECK(read_file(deck) == tensor_deck);  // edited, not saved
    // The modeller's panel offers no Unlock for a layered deck.
    CHECK_FALSE(w->solid_page().dimensions_panel().unlocks());
    // Undo: the tensor deck again, shown as the preview.
    REQUIRE(w->solid_page().undo());
    CHECK(w->deck().text() == tensor_deck);
    CHECK_FALSE(w->solid_shown());
    REQUIRE(pump_until([&] { return w->preview() && w->scene() != nullptr; }));
    CHECK(panel.unlocks());
}

TEST_CASE("unlock window: a refusal stays on the panel and the deck is unchanged (DT5)", "[app][window][unlock]") {
    const fs::path dir = app_test::scratch("unlock_window_refused");
    std::string text = tensor_deck;  // a contact with two boxes: a layered deck takes one patch box a contact
    text.replace(text.find("electrode = metal"), 0, "box = 0 0.1um  $H $H  0 $W\n");
    const fs::path deck = app_test::write(dir, "two.deck", text);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    w->execute(Command::tab_geometry);
    REQUIRE(pump_until([&] { return w->preview() && w->scene() != nullptr; }));
    DimensionsPanel& panel = w->geometry_page().dimensions_panel();
    click(panel);
    CHECK(panel.unlock_error().find("2 boxes") != std::string::npos);
    CHECK(w->deck().text() == text);
    CHECK_FALSE(w->solid_shown());
}

TEST_CASE("unlock window: an extrude3d deck converts from the modeller's panel (DT5)", "[app][window][unlock]") {
    const std::string text = read_file(fs::path(NITCAD_EXAMPLES_DIR) / ".." / "decks" / "finfet_3d_reduced.deck");
    REQUIRE(!text.empty());
    const fs::path dir = app_test::scratch("unlock_window_extruded");
    const fs::path deck = app_test::write(dir, "fin.deck", text);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    DimensionsPanel& panel = w->solid_page().dimensions_panel();
    const std::string shown = w->deck().text();  // the pane's own line endings
    REQUIRE(panel.unlocks());
    click(panel);
    const auto want = unlock::unlock(shown);
    REQUIRE(bool(want) == panel.unlock_error().empty());
    if (want) {
        CHECK(w->deck().text() == want->text);
        CHECK(w->solid_page().model().mesh == solid::MeshKind::layered);
        CHECK_FALSE(panel.unlocks());
        REQUIRE(w->solid_page().undo());
        CHECK(w->deck().text() == shown);
        CHECK(w->solid_page().model().mesh == solid::MeshKind::extruded);
    }
}
