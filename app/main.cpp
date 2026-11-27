// nitcad-app: the NiTCAD Windows application (Unit 26, DECISIONS.md VA1-VA9, AQ4). A separate executable
// from the `nitcad` batch runner, which keeps no window code.
//
//   nitcad-app [<deck>] [--adapter auto|gpu|warp] [--debug-layer]
//   nitcad-app --self-test <deck> --out <image.png> [view options]   (self_test.hpp; exit 0, 1, 64 or 77)
//   nitcad-app --capture <image.png> <deck> [--size <w>x<h>] [--aurora]   (the window drawn into a PNG, off-screen; 1600x968)
//
// Per-monitor DPI aware (v2); the common controls of version 6.
// The application's headers come first: objbase.h defines `interface` as a macro, which a member of that
// name in the headers they include would not survive.
#include "main_window.hpp"
#include "self_test.hpp"
#include "win32.hpp"

#include <windows.h>
#include <objbase.h>
#include <shellapi.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <format>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#pragma comment(linker, "/manifestdependency:\"type='win32' name='Microsoft.Windows.Common-Controls' " \
                        "version='6.0.0.0' processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")


namespace {

// stdout on the caller's console when it gave none (a GUI-subsystem process): --self-test's report and --capture's
// errors. A missing handle is NULL or INVALID_HANDLE_VALUE depending on how the process was started.
void attach_parent_console() {
    const HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
    if ((out == nullptr || out == INVALID_HANDLE_VALUE) && AttachConsole(ATTACH_PARENT_PROCESS)) {
        FILE* f = nullptr;
        freopen_s(&f, "CONOUT$", "w", stdout);
    }
}

int fail(const std::string& message) {
    std::fprintf(stdout, "nitcad-app: %s\n", message.c_str());
    std::fflush(stdout);
    return 1;
}

// --capture <image.png> <deck> [--size <w>x<h>] (docs/images/app.png): the window, off every monitor, draws itself into
// the PNG. Its client area is copied by PrintWindow from its own window (not the screen: no other application is in
// the image), and the device view, a Direct3D swap chain PrintWindow does not see, is the renderer's read-back image
// put in its place. The deck is opened as the application opens it (its results when its run wrote them).
int capture(const std::filesystem::path& png, const std::filesystem::path& deck, int width, int height, bool aurora) {
    using namespace NiTCAD;
    app::WindowOptions options;
    options.show = false;
    options.quit_on_close = false;
    auto window = app::MainWindow::create(options);
    if (!window) return fail("the main window could not be created");
    const HWND h = window->hwnd();
    // Shown off every monitor, without a taskbar button or the focus: a window never shown has nothing composed for
    // PrintWindow to copy, and this one is never on the screen. The frame is measured after the tool-window style
    // (a smaller caption) has taken effect, so the client area comes out at the size asked for.
    SetWindowLongPtrW(h, GWL_EXSTYLE, GetWindowLongPtrW(h, GWL_EXSTYLE) | WS_EX_TOOLWINDOW);
    SetWindowPos(h, HWND_BOTTOM, -32000, -32000, 0, 0, SWP_NOACTIVATE | SWP_NOSIZE | SWP_FRAMECHANGED);
    RECT outer{}, inner{};
    GetWindowRect(h, &outer);
    GetClientRect(h, &inner);
    const int frame_w = (outer.right - outer.left) - inner.right, frame_h = (outer.bottom - outer.top) - inner.bottom;
    SetWindowPos(h, HWND_BOTTOM, -32000, -32000, width + frame_w, height + frame_h, SWP_NOACTIVATE | SWP_SHOWWINDOW);
    // Windows caps a window at the virtual screen's size (WM_GETMINMAXINFO): a larger request is refused, not shrunk.
    GetClientRect(h, &inner);
    if (inner.right != width || inner.bottom != height) {
        return fail(std::format("the client area came out {}x{}, not the {}x{} asked for (larger than the screen?)",
                                inner.right, inner.bottom, width, height));
    }
    const auto pump_for = [&](double seconds, const auto& done) {
        const ULONGLONG end = GetTickCount64() + static_cast<ULONGLONG>(seconds * 1000);
        MSG m;
        while (!done() && GetTickCount64() < end) {
            if (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
                TranslateMessage(&m);
                DispatchMessageW(&m);
            } else {
                MsgWaitForMultipleObjects(0, nullptr, FALSE, 10, QS_ALLINPUT);
            }
        }
        return done();
    };
    if (aurora) window->execute(app::Command::aurora_theme);  // no state file here: nothing is saved
    if (!window->open_deck(deck)) return fail("cannot open " + deck.string() + ": " + window->status(0));
    // A deck that does not build ends the wait at once (as --self-test's), not after the time limit.
    const auto failed = [&] { return window->runner().state() == app::RunState::failed; };
    if (!pump_for(300.0, [&] { return failed() || (!window->runner().busy() && window->scene() != nullptr); }) ||
        failed()) {
        return fail(window->status(0));  // the state part says why ("The deck does not build: ...")
    }
    window->execute(app::Command::tab_view);
    // The device tree from its top (opening scrolls it to the item selected).
    EnumChildWindows(
        h,
        [](HWND c, LPARAM) -> BOOL {
            wchar_t name[32]{};
            GetClassNameW(c, name, 32);
            if (std::wstring_view(name) == L"SysTreeView32") SendMessageW(c, WM_VSCROLL, SB_TOP, 0);
            return TRUE;
        },
        0);
    window->refresh();
    pump_for(1.0, [] { return false; });
    window->refresh();

    // The controls: PrintWindow into a top-down 32-bit DIB.
    GetClientRect(h, &inner);
    const int w = inner.right, ht = inner.bottom;
    BITMAPINFO bi{};
    bi.bmiHeader.biSize = sizeof(bi.bmiHeader);
    bi.bmiHeader.biWidth = w;
    bi.bmiHeader.biHeight = -ht;
    bi.bmiHeader.biPlanes = 1;
    bi.bmiHeader.biBitCount = 32;
    bi.bmiHeader.biCompression = BI_RGB;
    void* bits = nullptr;
    const HDC screen = GetDC(nullptr);
    const HDC dc = CreateCompatibleDC(screen);
    ReleaseDC(nullptr, screen);
    const HBITMAP bmp = CreateDIBSection(dc, &bi, DIB_RGB_COLORS, &bits, nullptr, 0);
    if (!bmp || !bits) {
        if (bmp) DeleteObject(bmp);
        DeleteDC(dc);
        return fail(std::format("no {}x{} bitmap for the capture", w, ht));
    }
    const HGDIOBJ old = SelectObject(dc, bmp);
    PrintWindow(h, dc, PW_CLIENTONLY | PW_RENDERFULLCONTENT);
    GdiFlush();
    render::Image out;
    out.width = w;
    out.height = ht;
    out.rgba.resize(static_cast<std::size_t>(w) * static_cast<std::size_t>(ht) * 4);
    const auto* src = static_cast<const std::uint8_t*>(bits);
    for (std::size_t i = 0; i < out.rgba.size(); i += 4) {
        out.rgba[i] = src[i + 2];
        out.rgba[i + 1] = src[i + 1];
        out.rgba[i + 2] = src[i];
        out.rgba[i + 3] = 255;
    }
    SelectObject(dc, old);
    DeleteObject(bmp);
    DeleteDC(dc);

    // The device view: the renderer's last image, at the view's place in the client area.
    if (render::Renderer* r = window->renderer()) {
        if (auto view = r->last_image()) {
            RECT vr{};
            GetClientRect(window->viewport().hwnd(), &vr);
            // An image of another size is a frame from before the last resize: pasting it would spill over the
            // panels or leave part of the view black.
            if (view->width != vr.right || view->height != vr.bottom) {
                return fail(std::format("the device view's last image is {}x{}, the view {}x{}", view->width,
                                        view->height, vr.right, vr.bottom));
            }
            MapWindowPoints(window->viewport().hwnd(), h, reinterpret_cast<POINT*>(&vr), 2);
            const int bottom = std::min<int>(vr.bottom, ht), right = std::min<int>(vr.right, w);
            for (int y = 0; y < view->height && vr.top + y < bottom; ++y) {
                for (int x = 0; x < view->width && vr.left + x < right; ++x) {
                    if (vr.top + y < 0 || vr.left + x < 0) continue;
                    const std::size_t s = (static_cast<std::size_t>(y) * static_cast<std::size_t>(view->width) + static_cast<std::size_t>(x)) * 4;
                    const std::size_t d = (static_cast<std::size_t>(vr.top + y) * static_cast<std::size_t>(w) + static_cast<std::size_t>(vr.left + x)) * 4;
                    for (int c = 0; c < 3; ++c) out.rgba[d + static_cast<std::size_t>(c)] = view->rgba[s + static_cast<std::size_t>(c)];
                }
            }
        }
    }
    // The status bar left out: it shows the deck's full path (the user's folders).
    HWND status = FindWindowExW(h, nullptr, L"msctls_statusbar32", nullptr);
    if (status && IsWindowVisible(status)) {
        RECT sr{};
        GetWindowRect(status, &sr);
        MapWindowPoints(nullptr, h, reinterpret_cast<POINT*>(&sr), 2);
        if (sr.top > 0 && sr.top < out.height) {
            out.height = sr.top;
            out.rgba.resize(static_cast<std::size_t>(out.width) * static_cast<std::size_t>(out.height) * 4);
        }
    }
    return render::save_png(out, png) ? 0 : fail("cannot write " + png.string());
}

// "<w>x<h>", both positive, nothing after.
bool parse_size(const std::wstring& s, int& width, int& height) {
    int consumed = 0;
    return swscanf_s(s.c_str(), L"%dx%d%n", &width, &height, &consumed) == 2 &&
           consumed == static_cast<int>(s.size()) && width > 0 && height > 0;
}

int usage(const std::wstring& message) {
    MessageBoxW(nullptr,
                (message + L"\n\nnitcad-app [<deck>] [--adapter auto|gpu|warp] [--debug-layer]\n"
                           L"nitcad-app --self-test <deck> --out <image.png> [view options]\n"
                           L"nitcad-app --capture <image.png> <deck> [--size <w>x<h>] [--aurora]")
                    .c_str(),
                L"NiTCAD", MB_ICONERROR);
    return 64;
}

}  // namespace

int WINAPI wWinMain(HINSTANCE, HINSTANCE, PWSTR, int) {
    using namespace NiTCAD;
    SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
    if (FAILED(CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED))) return 1;

    int argc = 0;
    LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
    if (!argv) {
        CoUninitialize();
        return 1;
    }
    std::vector<std::wstring> args(argv + std::min(argc, 1), argv + argc);
    LocalFree(argv);

    // --self-test (26-3, WG4): its report on the caller's standard output (a pipe or the console).
    if (!args.empty() && args[0] == L"--self-test") {
        attach_parent_console();
        auto options = app::parse_self_test(std::span<const std::wstring>(args).subspan(1));
        int status = 64;
        if (options) status = app::run_self_test(*options, stdout);
        else std::fprintf(stdout, "%s\n", options.error().c_str());
        std::fflush(stdout);
        CoUninitialize();
        return status;
    }

    // --capture (docs/images/app.png): the window draws itself into a PNG, off-screen.
    if (!args.empty() && args[0] == L"--capture") {
        attach_parent_console();
        int width = 1600, height = 968;
        const bool aurora = !args.empty() && args.back() == L"--aurora";  // the Aurora theme (View > Aurora theme)
        if (aurora) args.pop_back();
        const bool sized = args.size() == 5 && args[3] == L"--size";
        if (!(args.size() == 3 || sized) || (sized && !parse_size(args[4], width, height))) {
            std::fprintf(stdout, "usage: nitcad-app --capture <image.png> <deck> [--size <w>x<h>] [--aurora] (w, h > 0)\n");
            std::fflush(stdout);
            CoUninitialize();
            return 64;
        }
        const int status = capture(args[1], args[2], width, height, aurora);
        CoUninitialize();
        return status;
    }

    app::WindowOptions options;
    options.animate = true;  // Unit 26h (JD5): smooth camera moves, if Windows shows animations
    options.state_file = app::win32::state_file();
    std::filesystem::path deck;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::wstring& a = args[i];
        if (a == L"--debug-layer") {
            options.debug_layer = true;
        } else if (a == L"--adapter") {
            const std::wstring v = i + 1 < args.size() ? args[++i] : std::wstring();
            if (v == L"warp") options.adapter = render::Adapter::warp;
            else if (v == L"gpu") options.adapter = render::Adapter::hardware;
            else if (v == L"auto") options.adapter = render::Adapter::automatic;
            else {
                CoUninitialize();
                return usage(L"--adapter takes auto, gpu or warp.");
            }
        } else if (a.starts_with(L"--") || !deck.empty()) {
            CoUninitialize();
            return usage(a.starts_with(L"--") ? L"Unknown option " + a + L"." : L"Only one deck can be opened.");
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
    MSG m{};
    BOOL got = 0;
    while ((got = GetMessageW(&m, nullptr, 0, 0)) > 0) {
        if (window->pre_translate(m)) continue;
        TranslateMessage(&m);
        DispatchMessageW(&m);
    }
    window.reset();
    CoUninitialize();
    return got < 0 ? 1 : static_cast<int>(m.wParam);  // -1: the loop failed; else WM_QUIT's exit code
}
