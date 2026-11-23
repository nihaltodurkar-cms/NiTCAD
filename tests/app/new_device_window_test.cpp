// Unit 26n slice 1 (DT4): File > New device..., by window messages and the dialog's own actions.
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "new_device_dialog.hpp"
#include "fixtures.hpp"

#include "NiTCAD/app/commands.hpp"

#include <windows.h>

#include <filesystem>
#include <fstream>
#include <string>

using namespace NiTCAD::app;
namespace fs = std::filesystem;

namespace {

void pump() {
    MSG m;
    while (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
        TranslateMessage(&m);
        DispatchMessageW(&m);
    }
}

std::unique_ptr<MainWindow> window(const fs::path& examples = {}, const fs::path& devices_copy = {}) {
    WindowOptions o;
    o.adapter = NiTCAD::render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    o.examples_dir = examples;
    o.devices_copy_dir = devices_copy;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    return w;
}

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

std::size_t family_index(NewDeviceDialog& d, const std::string& title) {
    for (std::size_t k = 0; k < d.families().size(); ++k) {
        if (d.families()[k].example.title == title) return k;
    }
    FAIL_CHECK("no family " << title);  // FAIL_CHECK: a FAIL before the return is C4702 in Debug
    return 0;
}

int changed_lines(const std::string& a, const std::string& b) {
    int n = 0;
    std::size_t i = 0, j = 0;
    while (i < a.size() || j < b.size()) {
        const auto ea = a.find('\n', i), eb = b.find('\n', j);
        n += a.substr(i, ea - i) != b.substr(j, eb - j);
        if (ea == std::string::npos || eb == std::string::npos) return n + (ea != eb);
        i = ea + 1;
        j = eb + 1;
    }
    return n;
}

}  // namespace

TEST_CASE("new device window: the command opens the dialog; Create writes the deck and opens it on the Geometry tab", "[app][window][new_device]") {
    const fs::path out = app_test::scratch("new_device_window_create");
    auto w = window();
    w->execute(Command::new_device);
    pump();
    NewDeviceDialog* d = w->new_device_dialog();
    REQUIRE(d);
    CHECK_FALSE(IsWindowEnabled(w->hwnd()));
    d->select(family_index(*d, "FinFET (bulk, 3D)"));
    CHECK(d->can_create());
    CHECK(d->message().empty());
    REQUIRE(d->panel().rows() == 5);
    const std::string before = d->text();
    REQUIRE(d->panel().commit(0, "30nm"));
    CHECK(changed_lines(before, d->text()) == 1);
    CHECK(d->text().find("\nLg = 30nm  # Gate length\n") != std::string::npos);
    d->set_folder(out);
    d->set_name("fet");
    REQUIRE(d->create());
    pump();
    CHECK(w->new_device_dialog() == nullptr);
    CHECK(IsWindowEnabled(w->hwnd()));
    REQUIRE(fs::exists(out / "fet.deck"));
    CHECK(read_file(out / "fet.deck").find("\nLg = 30nm  # Gate length\n") != std::string::npos);
    CHECK(w->document().path() == out / "fet.deck");
    CHECK(w->tab() == MainWindow::tab_geometry);
}

TEST_CASE("new device window: Cancel writes nothing; a second Create of the same name is refused", "[app][window][new_device]") {
    const fs::path out = app_test::scratch("new_device_window_cancel");
    auto w = window({}, out);
    w->execute(Command::new_device);
    pump();
    NewDeviceDialog* d = w->new_device_dialog();
    REQUIRE(d);
    d->select(family_index(*d, "PN diode (1D)"));
    CHECK(d->name() == "pn_diode_1d");  // the family's own name by default
    CHECK(d->text() == read_file(d->families()[d->selected()].example.path));  // its own deck, not the first family's
    CHECK(d->panel().rows() == 4);  // Ln, Lp, Nd, Na
    CHECK(d->folder() == out);          // the options' folder
    d->cancel();
    pump();
    CHECK(w->new_device_dialog() == nullptr);
    CHECK(IsWindowEnabled(w->hwnd()));
    CHECK_FALSE(fs::exists(out / "pn_diode_1d.deck"));
    // A file of that name is there: refused, the dialog stays open with the reason.
    std::ofstream(out / "pn_diode_1d.deck", std::ios::binary) << "# mine\n";
    w->execute(Command::new_device);
    pump();
    d = w->new_device_dialog();
    REQUIRE(d);
    d->select(family_index(*d, "PN diode (1D)"));
    CHECK_FALSE(d->create());
    pump();
    CHECK(w->new_device_dialog() == d);
    CHECK(d->message().find("already in") != std::string::npos);
    CHECK(read_file(out / "pn_diode_1d.deck") == "# mine\n");
    d->cancel();
    pump();
}

TEST_CASE("new device window: a deck that does not build is listed with its message and cannot be created", "[app][window][new_device]") {
    const fs::path examples = app_test::scratch("new_device_window_broken");
    fs::create_directories(examples / "devices");
    std::ofstream(examples / "devices" / "broken.deck", std::ios::binary) << "# Broken\n# It does not build.\n[mesh]\ntype = tensor\n";
    auto w = window(examples);
    w->execute(Command::new_device);
    pump();
    NewDeviceDialog* d = w->new_device_dialog();
    REQUIRE(d);
    REQUIRE(d->families().size() == 1);
    d->select(0);
    CHECK_FALSE(d->can_create());
    CHECK_FALSE(d->message().empty());
    CHECK_FALSE(d->create());
    d->cancel();
    pump();
}

TEST_CASE("new device window: a main dimension the build refuses shows on its row and is not written", "[app][window][new_device]") {
    auto w = window();
    w->execute(Command::new_device);
    pump();
    NewDeviceDialog* d = w->new_device_dialog();
    REQUIRE(d);
    d->select(family_index(*d, "Planar MOSFET (2D)"));
    const std::string before = d->text();
    const auto tox = std::ranges::find_if(d->panel().dimensions(), [](const auto& r) { return r.name == "tox"; }) - d->panel().dimensions().begin();
    REQUIRE(static_cast<std::size_t>(tox) < d->panel().rows());
    CHECK_FALSE(d->panel().commit(static_cast<std::size_t>(tox), "-1nm"));
    CHECK_FALSE(d->panel().row_error(static_cast<std::size_t>(tox)).empty());
    CHECK(d->text() == before);
    CHECK(d->can_create());
    d->cancel();
    pump();
}

TEST_CASE("new device window: the Start page's first button and the File menu run New device", "[app][window][new_device]") {
    auto w = window();
    HWND first = GetWindow(w->start_page().hwnd(), GW_CHILD);
    HWND button = nullptr;
    for (HWND h = first; h; h = GetWindow(h, GW_HWNDNEXT)) {
        wchar_t cls[16]{};
        GetClassNameW(h, cls, 16);
        if (std::wstring(cls) == L"Button") {
            button = h;
            break;
        }
    }
    REQUIRE(button);
    CHECK(GetDlgCtrlID(button) == static_cast<int>(Command::new_device));
    SendMessageW(button, BM_CLICK, 0, 0);
    pump();
    REQUIRE(w->new_device_dialog());
    w->new_device_dialog()->cancel();
    pump();
    const CommandInfo* c = find_command(Command::new_device);
    REQUIRE(c);
    CHECK(c->menu == "&File");
    CHECK(c->accelerator_text == "Ctrl+Shift+N");
}
