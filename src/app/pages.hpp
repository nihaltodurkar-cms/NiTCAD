// The main window's start page, run monitor and results pages (Unit 26, DECISIONS.md VA3, VA5, VA9,
// page Figures 3, 5, 6).
//
// - Start: new, open and open-result buttons, the recent files and the bundled examples (VQ11).
// - Run monitor: a row per solver point (a virtual list), the progress, Stop, and the live convergence
//   plot drawn by the renderer in its own swap chain (at most four times a second while running, and at
//   most a twentieth of the time).
// - Results: Unit 24's figures and the run record with its provenance (Unit 42).
#pragma once

#include <windows.h>

#include <chrono>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/app/examples.hpp"
#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "viewport.hpp"

namespace NiTCAD::app {

class StartPage {
public:
    struct Host {
        virtual void start_command(int command) = 0;      // a button: a Command's id
        virtual void start_recent(std::size_t index) = 0;
        virtual void start_example(std::size_t index) = 0;

    protected:
        ~Host() = default;
    };

    StartPage(HWND parent, Host& host, int id);
    ~StartPage();
    StartPage(const StartPage&) = delete;
    StartPage& operator=(const StartPage&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return panel_; }
    void set(const std::vector<std::filesystem::path>& recent, const std::vector<Example>& examples);
    void set_dpi(int dpi, HFONT font, HFONT bold);
    // What the lists show (for the tests).
    [[nodiscard]] int example_count() const;
    [[nodiscard]] std::wstring example_title(int index) const;
    [[nodiscard]] int recent_count() const;

private:
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    void layout();

    Host& host_;
    HWND panel_ = nullptr, examples_ = nullptr, recent_ = nullptr;
    std::vector<HWND> headings_, buttons_;
    int dpi_ = 96;
};

class RunPage final : public Viewport::Host {
public:
    struct Host {
        virtual void run_page_stop() = 0;

    protected:
        ~Host() = default;
    };

    RunPage(HWND parent, Host& host, int id);
    ~RunPage();
    RunPage(const RunPage&) = delete;
    RunPage& operator=(const RunPage&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return panel_; }
    [[nodiscard]] Viewport& plot() noexcept { return *plot_; }
    // The controller's rows, progress and state; the plot when `force` or a quarter second has passed.
    void update(const RunController& runner, render::Renderer* renderer, bool force = false);
    void set_dpi(int dpi, HFONT font);
    void release();  // the plot's swap chain, before the renderer
    [[nodiscard]] int row_count() const;
    [[nodiscard]] std::wstring cell(int row, int column) const;

private:
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    void layout();
    // Viewport::Host: the plot does not move.
    void viewport_orbit(float, float) override {}
    void viewport_pan(float, float) override {}
    void viewport_zoom(float) override {}
    void viewport_fit() override {}
    void viewport_click(float, float) override {}
    void viewport_hover(float, float) override {}
    bool viewport_is_3d() const override { return false; }

    Host& host_;
    HWND panel_ = nullptr, list_ = nullptr, progress_ = nullptr, label_ = nullptr, stop_ = nullptr;
    std::unique_ptr<Viewport> plot_;
    std::optional<Scene> scene_;
    std::vector<PointRow> rows_;
    std::chrono::steady_clock::time_point plotted_{};
    int dpi_ = 96;
};

class ResultsPage {
public:
    ResultsPage(HWND parent, int id);
    ~ResultsPage();
    ResultsPage(const ResultsPage&) = delete;
    ResultsPage& operator=(const ResultsPage&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return panel_; }
    void show(const Workspace* w);
    void set_dpi(int dpi, HFONT font, HFONT bold);
    [[nodiscard]] int figure_count() const;
    [[nodiscard]] int record_count() const;
    [[nodiscard]] std::wstring figure_cell(int row, int column) const;

private:
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    void layout();

    HWND panel_ = nullptr, figures_ = nullptr, record_ = nullptr, figures_label_ = nullptr, record_label_ = nullptr;
    int dpi_ = 96;
};

}  // namespace NiTCAD::app
