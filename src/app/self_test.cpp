// The self-test (self_test.hpp).
#include "self_test.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <format>
#include <variant>

#include "main_window.hpp"
#include "win32.hpp"

namespace NiTCAD::app {

namespace {

constexpr UINT_PTR gap_timer = 0x5E1F;

std::optional<double> number(const std::string& s) {
    double v = 0.0;
    const auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size()) return std::nullopt;
    return v;
}

void say(std::FILE* log, const std::string& line) {
    if (!log) return;
    std::fprintf(log, "self-test: %s\n", line.c_str());
    std::fflush(log);
}

}  // namespace

std::expected<SelfTestOptions, std::string> parse_self_test(std::span<const std::wstring> args) {
    SelfTestOptions o;
    o.view.kind = ViewKind::field;
    bool view_given = false;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string a = win32::narrow(args[i]);
        const auto value = [&]() -> std::optional<std::string> {
            if (i + 1 >= args.size()) return std::nullopt;
            return win32::narrow(args[++i]);
        };
        if (a == "--out") {
            const auto v = value();
            if (!v) return std::unexpected(std::string("--out needs a path"));
            o.out = std::filesystem::path(args[i]);
        } else if (a == "--view") {
            const auto v = value();
            const auto k = v ? find_view(*v) : std::nullopt;
            if (!k) return std::unexpected(std::string("--view is solid, mesh, field, slice, iso, bands, profile, iv, cv or current"));
            o.view.kind = *k;
            view_given = true;
        } else if (a == "--field" || a == "--scale" || a == "--colormap" || a == "--sense" || a == "--contact") {
            const auto v = value();
            if (!v) return std::unexpected(a + " needs a value");
            (a == "--field" ? o.view.field : a == "--scale" ? o.view.scale : a == "--colormap" ? o.view.colormap
                                                                  : a == "--sense"  ? o.view.sense
                                                                                    : o.view.contact) = *v;
        } else if (a == "--point") {
            const auto v = value();
            const auto n = v ? number(*v) : std::nullopt;
            if (!n || *n < 0) return std::unexpected(std::string("--point needs an index"));
            o.view.point = static_cast<std::size_t>(*n);
        } else if (a == "--size") {
            const auto v = value();
            const auto x = v ? v->find('x') : std::string::npos;
            const auto w = x != std::string::npos ? number(v->substr(0, x)) : std::nullopt;
            const auto h = x != std::string::npos ? number(v->substr(x + 1)) : std::nullopt;
            if (!w || !h) return std::unexpected(std::string("--size needs <w>x<h>"));
            o.width = static_cast<int>(*w);
            o.height = static_cast<int>(*h);
        } else if (a == "--adapter") {
            const auto v = value();
            if (!v || (*v != "auto" && *v != "gpu" && *v != "warp")) return std::unexpected(std::string("--adapter is auto, gpu or warp"));
            o.adapter = *v == "gpu" ? render::Adapter::hardware : *v == "warp" ? render::Adapter::warp : render::Adapter::automatic;
        } else if (a == "--debug-layer") o.debug_layer = true;
        else if (a == "--contours") o.view.contours = true;
        else if (a == "--edges") o.view.edges = true;
        else if (a == "--opaque") o.view.opaque = true;
        else if (a == "--perspective") o.view.perspective = true;
        else if (a == "--log") o.view.log = true;
        else if (!a.starts_with("--") && o.deck.empty()) o.deck = std::filesystem::path(args[i]);
        else return std::unexpected("unknown option " + a);
    }
    if (o.deck.empty() || o.out.empty()) return std::unexpected(std::string("usage: nitcad-app --self-test <deck> --out <image.png> [view options]"));
    (void)view_given;
    return o;
}

int run_self_test(const SelfTestOptions& o, std::FILE* log) {
    WindowOptions wo;
    wo.adapter = o.adapter;
    wo.debug_layer = o.debug_layer;
    wo.show = true;
    wo.quit_on_close = false;
    auto window = MainWindow::create(wo);
    if (!window) {
        say(log, "no window can be made here (skipped)");
        return self_test_skipped;
    }
    if (!window->renderer()) {
        say(log, "no Direct3D 12 device: " + window->renderer_error().value_or("unknown"));
        return 1;
    }
    say(log, "window open, " + window->renderer()->adapter_name());
    // The message loop, with a 10 ms timer whose longest gap is how long the window could not run (WG1).
    double longest_gap = 0.0, after_run = 0.0;
    auto last = std::chrono::steady_clock::now();
    bool timing = false;
    const auto pump_until = [&](const auto& done, double seconds) {
        const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
        MSG m;
        while (!done()) {
            if (std::chrono::steady_clock::now() > end) return false;
            if (!GetMessageW(&m, nullptr, 0, 0)) return false;
            if (m.message == WM_TIMER && m.wParam == gap_timer) {
                const auto now = std::chrono::steady_clock::now();
                if (timing) longest_gap = std::max(longest_gap, std::chrono::duration<double>(now - last).count());
                last = now;
                continue;
            }
            // The message that ends the run loads its results and draws the first view on the UI thread: that
            // work follows the solve and is measured on its own, not as a gap during the run.
            const bool running = window->runner().busy();
            const auto started = std::chrono::steady_clock::now();
            if (!window->pre_translate(m)) {
                TranslateMessage(&m);
                DispatchMessageW(&m);
            }
            if (timing && running && !window->runner().busy()) {
                timing = false;
                after_run = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            }
        }
        return true;
    };
    SetTimer(window->hwnd(), gap_timer, 10, nullptr);
    if (!window->open_deck(o.deck)) {
        say(log, "cannot open " + o.deck.string() + ": " + window->status(0));
        return 1;
    }
    if (!pump_until([&] { return window->workspace() != nullptr || window->runner().state() == RunState::failed; }, 600.0) ||
        !window->workspace()) {
        say(log, "the deck does not build: " + window->status(0));
        return 1;
    }
    if (window->viewport().error()) {
        say(log, "no swap chain here (skipped): " + *window->viewport().error());
        return self_test_skipped;
    }
    // Run on the worker thread; the window keeps pumping.
    window->execute(Command::run);
    timing = true;
    last = std::chrono::steady_clock::now();
    if (!pump_until([&] { return !window->runner().busy(); }, 3600.0)) {
        say(log, "the run did not finish");
        return 1;
    }
    timing = false;
    if (!window->workspace() || !window->workspace()->has_results()) {
        say(log, "the run failed: " + window->status(0));
        return 1;
    }
    say(log, std::format("ran: {} points in {:.2f} s; the longest gap in the window's message loop {:.1f} ms (the run "
                         "monitor's plot last drawn in {:.1f} ms); the results loaded and shown {:.0f} ms after the run",
                         window->workspace()->points.size(), window->runner().elapsed_seconds(), 1e3 * longest_gap,
                         window->run_page().plot().last_present_ms(), 1e3 * after_run));
    if (longest_gap >= 0.1) {
        say(log, "the window did not run for 100 ms or more during the run (WG1)");
        return 1;
    }
    // The view, in the window and as an image.
    window->select_tab(is_plot(o.view.kind) && o.view.kind != ViewKind::profile ? MainWindow::tab_curves : MainWindow::tab_view);
    window->settings() = o.view;
    window->refresh();
    // Let the window paint what it has.
    for (MSG m; PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE);) {
        TranslateMessage(&m);
        DispatchMessageW(&m);
    }
    if (!window->scene()) {
        say(log, "the view cannot be built: " + window->status(0));
        return 1;
    }
    if (window->viewport().error()) {
        say(log, "the view cannot be shown: " + *window->viewport().error());
        return 1;
    }
    auto image = std::visit([&](const auto& s) { return window->renderer()->draw(s, o.width, o.height); }, *window->scene());
    if (!image) {
        say(log, "the view cannot be drawn: " + image.error().message);
        return 1;
    }
    if (auto saved = render::save_png(*image, o.out); !saved) {
        say(log, saved.error().message);
        return 1;
    }
    say(log, std::format("wrote {} ({}x{}), shown in the window at {}x{}", o.out.string(), image->width, image->height,
                         window->viewport().width(), window->viewport().height()));
    if (o.debug_layer) {
        const auto messages = window->renderer()->debug_messages();
        for (const auto& m : messages) say(log, "debug layer: " + m);
        if (!messages.empty()) return 1;
    }
    KillTimer(window->hwnd(), gap_timer);
    window.reset();
    if (o.debug_layer) {
        const int live = win32::live_graphics_objects();
        if (live > 0) {
            say(log, std::format("{} graphics objects outlived the window", live));
            return 1;
        }
        say(log, live < 0 ? "live objects not checked (no DXGI debug interface)" : "no graphics object outlived the window");
    }
    say(log, "passed");
    return 0;
}

}  // namespace NiTCAD::app
