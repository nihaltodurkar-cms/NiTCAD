// nitcad-app: the NiTCAD Windows application (Unit 26, DECISIONS.md VA1-VA9, AQ4). A separate executable
// from the `nitcad` batch runner, which keeps no window code.
//
//   nitcad-app [<deck>] [--adapter auto|gpu|warp] [--debug-layer]
//
// Per-monitor DPI aware (v2); the common controls of version 6.
// The application's headers come first: objbase.h defines `interface` as a macro, which a member of that
// name in the headers they include would not survive.
#include "main_window.hpp"
#include "win32.hpp"

#include <windows.h>
#include <objbase.h>
#include <shellapi.h>

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#pragma comment(linker, "/manifestdependency:\"type='win32' name='Microsoft.Windows.Common-Controls' " \
                        "version='6.0.0.0' processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

int WINAPI wWinMain(HINSTANCE, HINSTANCE, PWSTR, int) {
    using namespace NiTCAD;
    SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
    if (FAILED(CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED))) return 1;

    int argc = 0;
    LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
    std::vector<std::wstring> args(argv + 1, argv + argc);
    LocalFree(argv);

    app::WindowOptions options;
    options.state_file = app::win32::state_file();
    std::filesystem::path deck;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::wstring& a = args[i];
        if (a == L"--debug-layer") {
            options.debug_layer = true;
        } else if (a == L"--adapter" && i + 1 < args.size()) {
            const std::wstring& v = args[++i];
            options.adapter = v == L"warp" ? render::Adapter::warp : v == L"gpu" ? render::Adapter::hardware
                                                                                  : render::Adapter::automatic;
        } else {
            deck = a;
        }
    }

    auto window = app::MainWindow::create(options);
    if (!window) {
        MessageBoxW(nullptr, L"The main window could not be created.", L"NiTCAD", MB_ICONERROR);
        CoUninitialize();
        return 1;
    }
    if (!deck.empty()) window->open_deck(deck);
    MSG m;
    while (GetMessageW(&m, nullptr, 0, 0) > 0) {
        if (window->pre_translate(m)) continue;
        TranslateMessage(&m);
        DispatchMessageW(&m);
    }
    window.reset();
    CoUninitialize();
    return 0;
}
