// Unit 26n slice 3 (DT6, DT7, DC10): Sweep... on a dimension row, the sweep dialog, the run, and the Curves tab's
// sweep views, by window messages and the dialog's own actions.
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "new_device_dialog.hpp"
#include "sweep_dialog.hpp"
#include "fixtures.hpp"

#include "NiTCAD/app/commands.hpp"

#include <windows.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <variant>

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

bool pump_until(const std::function<bool()>& done, double seconds = 300.0) {
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

std::size_t row_of(const DimensionsPanel& p, const std::string& name) {
    for (std::size_t k = 0; k < p.rows(); ++k) {
        if (p.dimensions()[k].name == name) return k;
    }
    FAIL_CHECK("no row " << name);  // FAIL_CHECK: a FAIL before the return is C4702 in Debug
    return 0;
}

const NiTCAD::render::PlotScene* plot_of(const MainWindow& w) {
    const Scene* s = w.scene();
    return s ? std::get_if<NiTCAD::render::PlotScene>(s) : nullptr;
}

}  // namespace

TEST_CASE("sweep window: Sweep... runs the sweep; the Curves tab overlays its runs; figures.csv (DT6, DT7)", "[app][window][sweep]") {
    const fs::path dir = app_test::scratch("sweep_window_run");
    const std::string text = read_file(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "pn_diode_1d.deck");
    const fs::path deck = app_test::write(dir, "pn.deck", text);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    w->execute(Command::tab_geometry);
    REQUIRE(pump_until([&] { return w->preview() && w->scene() != nullptr; }));
    DimensionsPanel& panel = w->geometry_page().dimensions_panel();
    const std::size_t lp = row_of(panel, "Lp");
    REQUIRE(panel.sweep_button(lp) != nullptr);
    SendMessageW(panel.hwnd(), WM_COMMAND, MAKEWPARAM(GetDlgCtrlID(panel.sweep_button(lp)), BN_CLICKED),
                 reinterpret_cast<LPARAM>(panel.sweep_button(lp)));
    pump();
    SweepDialog* d = w->sweep_dialog();
    REQUIRE(d);
    CHECK(d->name() == "Lp");
    CHECK_FALSE(IsWindowEnabled(w->hwnd()));
    d->set_values("1um, abc");
    CHECK_FALSE(d->can_run());
    CHECK(d->message().find("abc") != std::string::npos);
    CHECK_FALSE(d->run());
    d->set_values("1, 1.5");  // a bare number takes the row's unit (um)
    CHECK(d->can_run());
    CHECK(d->values() == std::vector<std::string>{"1um", "1.5um"});
    CHECK(d->message().find("2 runs") != std::string::npos);
    REQUIRE(d->run());
    pump();
    CHECK(w->sweep_dialog() == nullptr);
    CHECK(IsWindowEnabled(w->hwnd()));
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->sweep_results() != nullptr; }));
    CHECK(w->runner().state() == RunState::done);
    CHECK(read_file(deck) == text);  // the deck is not edited
    const SweepResults& r = *w->sweep_results();
    REQUIRE(r.runs.size() == 2);
    CHECK(r.runs[0].workspace);
    CHECK(r.runs[1].workspace);
    CHECK(w->tab() == MainWindow::tab_curves);
    REQUIRE(pump_until([&] { return plot_of(*w) != nullptr; }));
    const auto* overlay = plot_of(*w);
    REQUIRE(overlay->series.size() == 2);
    CHECK(overlay->series[0].name == "Lp = 1 \xC2\xB5m");
    CHECK(overlay->series[1].name == "Lp = 1.5 \xC2\xB5m");
    // figures.csv: a header and a row per run.
    const std::string csv = read_file(dir / "pn_sweep_Lp" / "figures.csv");
    CHECK(csv.starts_with("Lp (nm),"));
    CHECK(std::ranges::count(csv, '\n') == 3);
    // The view list: the sweep's curves, then each figure against Lp.
    REQUIRE(w->sweep_views().size() == 2);
    CHECK(w->sweep_views()[0] == "sweep: curves");
    CHECK(w->sweep_views()[1] == "sweep: Current at the last bias vs Lp");
    w->select_sweep_view(1);
    REQUIRE(pump_until([&] { return plot_of(*w) && plot_of(*w)->series.size() == 1; }));
    CHECK(plot_of(*w)->series[0].x == std::vector<double>{1000.0, 1500.0});
    CHECK(plot_of(*w)->series[0].markers);
}

TEST_CASE("sweep window: a derived row has no Sweep...; Cancel runs nothing", "[app][window][sweep]") {
    const fs::path dir = app_test::scratch("sweep_window_cancel");
    std::string deck_text = read_file(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "pn_diode_1d.deck");
    deck_text.replace(deck_text.find("Lp = 1.2um"), 10, "Lp = $Ln + 0.4um");
    const fs::path deck = app_test::write(dir, "pn.deck", deck_text);
    auto w = window();
    REQUIRE(w->open_deck(deck));
    w->execute(Command::tab_geometry);
    REQUIRE(pump_until([&] { return w->preview() && w->scene() != nullptr; }));
    DimensionsPanel& panel = w->geometry_page().dimensions_panel();
    CHECK(panel.sweep_button(row_of(panel, "Lp")) == nullptr);
    const std::size_t ln = row_of(panel, "Ln");
    REQUIRE(panel.sweep(ln));
    pump();
    SweepDialog* d = w->sweep_dialog();
    REQUIRE(d);
    CHECK(d->name() == "Ln");
    d->cancel();
    pump();
    CHECK(w->sweep_dialog() == nullptr);
    CHECK(IsWindowEnabled(w->hwnd()));
    CHECK_FALSE(w->runner().busy());
    CHECK_FALSE(fs::exists(dir / "pn_sweep_Ln"));
}

// DC10 as the criterion words it: New device > FinFET (3D), Lg = 30 nm and Wfin = 8 nm in the panel, Sweep... Lg =
// 20, 30 nm, run. Hidden: two 3D runs (about a minute in Release).
TEST_CASE("sweep window DC10: a FinFET from New device, two dimensions set, a sweep of Lg", "[app][window][sweep][.sweep_dc10]") {
    const fs::path out = app_test::scratch("sweep_window_dc10");
    auto w = window();
    w->execute(Command::new_device);
    pump();
    NewDeviceDialog* nd = w->new_device_dialog();
    REQUIRE(nd);
    std::size_t family = nd->families().size();
    for (std::size_t k = 0; k < nd->families().size(); ++k) {
        if (nd->families()[k].example.title == "FinFET (bulk, 3D)") family = k;
    }
    REQUIRE(family < nd->families().size());
    nd->select(family);
    nd->set_folder(out);
    nd->set_name("fet");
    REQUIRE(nd->create());
    pump();
    REQUIRE(w->tab() == MainWindow::tab_geometry);
    REQUIRE(pump_until([&] { return w->preview() && w->scene() != nullptr; }));
    DimensionsPanel& panel = w->geometry_page().dimensions_panel();
    REQUIRE(panel.commit(row_of(panel, "Lg"), "30nm"));
    REQUIRE(panel.commit(row_of(panel, "Wfin"), "8nm"));
    REQUIRE(panel.sweep(row_of(panel, "Lg")));
    pump();
    SweepDialog* d = w->sweep_dialog();
    REQUIRE(d);
    d->set_values("20, 30");
    REQUIRE(d->run());
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->sweep_results() != nullptr; }, 1800.0));
    CHECK(w->runner().state() == RunState::done);
    REQUIRE(pump_until([&] { return plot_of(*w) != nullptr; }));
    CHECK(plot_of(*w)->series.size() == 2);
    const std::string csv = read_file(out / "fet_sweep_Lg" / "figures.csv");
    CHECK(std::ranges::count(csv, '\n') == 3);
}
