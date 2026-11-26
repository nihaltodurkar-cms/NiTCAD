// Sweep... (Unit 26n slice 3, DT6, DECISIONS.md): an owned popup over the main window, which is disabled while it is
// open (no nested message loop), for one dimension row: the values (`20nm, 30nm, 45nm` or `20nm : 45nm : 5nm`, a bare
// number taking the row's unit), what they give (the runs and the folder, or why not), Run and Cancel. Run hands the
// values to the host, which runs the sweep (RunController::start_sweep). The model is app/sweep_results.hpp.
#pragma once

#include <windows.h>

#include <filesystem>
#include <string>
#include <vector>

namespace NiTCAD::app {

class SweepDialog {
public:
    struct Host {
        virtual void sweep_run(const std::string& name, const std::vector<std::string>& values) = 0;  // closed later (posted)
        virtual void sweep_closed() = 0;                                                               // Cancel or the close box
    protected:
        ~Host() = default;
    };

    // `deck` the text swept, saved at `path` (the folder named beside it); `name` and `label` the row's.
    SweepDialog(HWND owner, Host& host, std::string deck, std::filesystem::path path, std::string name, const std::string& label,
                const std::string& current, int dpi, HFONT font, HFONT bold);
    ~SweepDialog();
    SweepDialog(const SweepDialog&) = delete;
    SweepDialog& operator=(const SweepDialog&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return hwnd_; }
    [[nodiscard]] const std::string& name() const noexcept { return name_; }
    void set_values(const std::string& typed);  // as typing in the values box
    [[nodiscard]] const std::vector<std::string>& values() const noexcept { return values_; }
    [[nodiscard]] bool can_run() const noexcept { return !values_.empty(); }
    [[nodiscard]] const std::string& message() const noexcept { return message_; }
    bool run();     // as the Run button
    void cancel();  // as Cancel

private:
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    void read_values();
    void layout();

    HWND owner_ = nullptr;
    Host& host_;
    HWND hwnd_ = nullptr, title_ = nullptr, edit_ = nullptr, hint_ = nullptr, message_label_ = nullptr;
    HWND run_ = nullptr, cancel_ = nullptr;
    std::string deck_;
    std::filesystem::path path_;
    std::string name_, message_;
    std::vector<std::string> values_;
    int dpi_ = 96;
    HFONT font_ = nullptr, bold_ = nullptr;
};

}  // namespace NiTCAD::app
