// nitcad-app --self-test (Unit 26, DECISIONS.md VA6, WG1, WG4): the window opened, a deck loaded and run
// on the worker thread while the window's message loop keeps running (no gap of 100 ms while the solver
// runs; the loading of the results when it ends is reported apart), a view shown in the window's swap
// chain and written as a PNG equal to `nitcad plot`'s for the same view; then the window closed.
//
//   nitcad-app --self-test <deck> --out <image.png> [--view v] [--field f] [--point k] [--scale s]
//              [--colormap c] [--contours] [--edges] [--opaque] [--perspective] [--log] [--sense c]
//              [--contact c] [--size WxH] [--adapter auto|gpu|warp] [--debug-layer]
//
// Exit status: 0 passed; 1 failed (the reason on stdout); 64 usage; 77 no window or swap chain can be made
// here (a session without a desktop: the test that runs it counts that as skipped).
#pragma once

#include <cstdio>
#include <expected>
#include <filesystem>
#include <span>
#include <string>

#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/render/image.hpp"

namespace NiTCAD::app {

struct SelfTestOptions {
    std::filesystem::path deck;
    std::filesystem::path out;
    ViewSettings view;  // the default: the potential's field view
    int width = 640, height = 480;
    render::Adapter adapter = render::Adapter::automatic;
    bool debug_layer = false;
};

inline constexpr int self_test_skipped = 77;

// The options after "--self-test". Errors: the usage message.
[[nodiscard]] std::expected<SelfTestOptions, std::string> parse_self_test(std::span<const std::wstring> args);

// Runs it, writing a line per step to `log`; returns the exit status.
[[nodiscard]] int run_self_test(const SelfTestOptions& options, std::FILE* log);

}  // namespace NiTCAD::app
