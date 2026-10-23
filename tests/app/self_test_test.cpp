// nitcad-app --self-test (Unit 26, WG1, WG4): the window runs a deck with its message loop running, and
// writes the image `nitcad plot` writes for the same view, on WARP and on the GPU. Where no window or
// swap chain can be made (a session without a desktop), the self-test exits 77 and this test is skipped.
#include <catch2/catch_test_macros.hpp>

#include "NiTCAD/app/examples.hpp"
#include "NiTCAD/render/image.hpp"
#include "fixtures.hpp"

#include <windows.h>

#include <algorithm>
#include <filesystem>
#include <format>
#include <string>
#include <vector>

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

struct Ran {
    DWORD status = 0;
    std::string output;
};

// Runs a command line, its standard output and error collected, and waits for it.
Ran run(const std::wstring& command) {
    SECURITY_ATTRIBUTES sa{sizeof(sa), nullptr, TRUE};
    HANDLE read = nullptr, write = nullptr;
    REQUIRE(CreatePipe(&read, &write, &sa, 0));
    SetHandleInformation(read, HANDLE_FLAG_INHERIT, 0);
    STARTUPINFOW si{sizeof(si)};
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdOutput = si.hStdError = write;
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    PROCESS_INFORMATION pi{};
    std::wstring line = command;
    REQUIRE(CreateProcessW(nullptr, line.data(), nullptr, nullptr, TRUE, 0, nullptr, nullptr, &si, &pi));
    CloseHandle(write);
    Ran r;
    char buffer[4096];
    DWORD got = 0;
    while (ReadFile(read, buffer, sizeof buffer, &got, nullptr) && got > 0) r.output.append(buffer, got);
    WaitForSingleObject(pi.hProcess, INFINITE);
    GetExitCodeProcess(pi.hProcess, &r.status);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    CloseHandle(read);
    return r;
}

std::wstring q(const fs::path& p) { return L"\"" + p.wstring() + L"\""; }
std::wstring w(const std::string& s) { return std::wstring(s.begin(), s.end()); }
// WG1's bound for the self-test (fixtures.hpp: 250 ms on the hosted CI runner).
std::wstring gap_limit() { return L" --gap-limit-ms " + w(std::format("{:g}", app_test::wg1_limit_ms())); }

// The self-test's image and nitcad plot's for one view of a deck; false when skipped (77).
bool same_view(const fs::path& deck, const std::string& view, const std::string& adapter, const fs::path& out) {
    const fs::path app_png = out / "app.png", plot_png = out / "plot.png";
    const Ran a = run(q(NITCAD_APP_PATH) + L" --self-test " + q(deck) + L" --out " + q(app_png) + L" --size 640x480 --adapter " +
                      w(adapter) + gap_limit() + L" " + w(view));
    UNSCOPED_INFO(a.output);
    if (a.status == 77) return false;
    REQUIRE(a.status == 0);
    CHECK(a.output.find("self-test: passed") != std::string::npos);
    CHECK(a.output.find("longest gap in the window's message loop") != std::string::npos);
    const Ran p = run(q(NITCAD_CLI_PATH) + L" plot " + q(deck) + L" " + w(view) + L" --size 640x480 --adapter " + w(adapter) +
                      L" --out " + q(plot_png));
    UNSCOPED_INFO(p.output);
    REQUIRE(p.status == 0);
    const auto x = render::load_png(app_png);
    const auto y = render::load_png(plot_png);
    REQUIRE(x);
    REQUIRE(y);
    const auto d = render::compare(*x, *y, 0);
    REQUIRE(d);
    CHECK(d->max_channel == 0);
    return true;
}

fs::path example(const std::string& name, const fs::path& dir) {
    const auto examples = app::list_examples(NITCAD_EXAMPLES_DIR);
    const auto it = std::ranges::find_if(examples, [&](const app::Example& e) { return e.path.filename() == name; });
    REQUIRE(it != examples.end());
    auto copy = app::copy_example(*it, dir);
    REQUIRE(copy);
    return *copy;
}

}  // namespace

TEST_CASE("app self-test: the window's image equals nitcad plot's, on WARP (WG1, WG4)", "[app]") {
    const fs::path dir = app_test::scratch("self_test_warp");
    const fs::path m2 = example("mosfet_2d.deck", dir);
    for (const std::string view : {"--view field --field potential", "--view iv --sense drain --log",
                                   "--view field --field n --colormap viridis --contours"}) {
        CAPTURE(view);
        if (!same_view(m2, view, "warp", dir)) SKIP("no window or swap chain in this session (exit 77)");
    }
    // Usage errors.
    CHECK(run(q(NITCAD_APP_PATH) + L" --self-test").status == 64);
    CHECK(run(q(NITCAD_APP_PATH) + L" --self-test " + q(m2) + L" --out x.png --view bogus").status == 64);
    CHECK(run(q(NITCAD_APP_PATH) + L" --self-test " + q(m2) + L" --out x.png --gap-limit-ms 0").status == 64);
    CHECK(run(q(NITCAD_APP_PATH) + L" --self-test " + q(m2) + L" --out x.png --gap-limit-ms").status == 64);
    // WG1's bound is the one given: no message loop keeps every gap under a microsecond (a 1D diode, quick;
    // 701 points, so the run spans many 10 ms timer ticks in Release too).
    const fs::path diode = app_test::write(dir, "diode.deck", app_test::diode_sweep("0 : 0.7 : 0.001"));
    const Ran tight = run(q(NITCAD_APP_PATH) + L" --self-test " + q(diode) + L" --out " + q(dir / "diode.png") +
                          L" --adapter warp --view iv --gap-limit-ms 0.001");
    UNSCOPED_INFO(tight.output);
    CHECK(tight.status == 1);
    CHECK(tight.output.find("did not run for 0.001 ms or more during the run (WG1)") != std::string::npos);
}

TEST_CASE("app self-test: on the GPU too, and with the debug layer (WG4, WG5)", "[app]") {
    const fs::path dir = app_test::scratch("self_test_gpu");
    const fs::path m2 = example("mosfet_2d.deck", dir);
    // A hardware adapter, when this machine has one.
    const Ran probe = run(q(NITCAD_CLI_PATH) + L" plot " + q(m2) + L" --view solid --adapter gpu --size 64x64 --out " +
                          q(dir / "probe.png"));
    if (probe.status != 0) SKIP("no Direct3D 12 hardware adapter: " + probe.output);
    if (!same_view(m2, "--view field --field potential", "gpu", dir)) SKIP("no window or swap chain in this session");
    // The debug layer (Graphics Tools): clean, and nothing left alive after the window.
    const Ran d = run(q(NITCAD_APP_PATH) + L" --self-test " + q(m2) + L" --out " + q(dir / "debug.png") +
                      L" --adapter warp --debug-layer --view solid" + gap_limit());
    UNSCOPED_INFO(d.output);
    if (d.output.find("D3D12GetDebugInterface") != std::string::npos) SKIP("the D3D12 debug layer is not installed");
    CHECK(d.status == 0);
    CHECK(d.output.find("no graphics object outlived the window") != std::string::npos);
}

// The 3D MOSFET runs in 4 s in Release and about 5 minutes in Debug: its self-tests run in the Release
// suite (ctest app_window_3d, hidden tag [.selftest3d]).
TEST_CASE("app self-test: the 3D MOSFET, on WARP and on the GPU (WG4)", "[.selftest3d]") {
    const fs::path dir = app_test::scratch("self_test_3d");
    const fs::path m3 = example("mosfet_3d.deck", dir);
    if (!same_view(m3, "--view field --field n", "warp", dir)) SKIP("no window or swap chain in this session (exit 77)");
    const Ran probe = run(q(NITCAD_CLI_PATH) + L" plot " + q(m3) + L" --view solid --adapter gpu --size 64x64 --out " +
                          q(dir / "probe.png"));
    if (probe.status != 0) SKIP("no Direct3D 12 hardware adapter: " + probe.output);
    CHECK(same_view(m3, "--view field --field potential", "gpu", dir));
}
