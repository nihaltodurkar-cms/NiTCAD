// The deck pane (Unit 26, DECISIONS.md VA2): the deck's text in a Rich Edit control with line numbers
// beside it, and below it the check's errors, each taking the caret to its line.
#pragma once

#include <windows.h>

#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/document.hpp"

namespace NiTCAD::app {

class DeckPane {
public:
    struct Host {
        virtual void deck_changed() = 0;  // the user edited the text

    protected:
        ~Host() = default;
    };

    DeckPane(HWND parent, Host& host, int id);
    ~DeckPane();
    DeckPane(const DeckPane&) = delete;
    DeckPane& operator=(const DeckPane&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return panel_; }
    [[nodiscard]] HWND editor() const noexcept { return edit_; }
    [[nodiscard]] HWND errors() const noexcept { return list_; }

    // The text with '\n' line ends (UTF-8). Setting it does not notify the host.
    [[nodiscard]] std::string text() const;
    void set_text(std::string_view utf8);
    void set_diagnostics(const std::vector<Diagnostic>& d);
    // Selects the line (1-based) and scrolls it into view; the editor takes the focus.
    void go_to_line(int line);
    [[nodiscard]] int caret_line() const;
    void set_dpi(int dpi);

private:
    static LRESULT CALLBACK panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK gutter_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    void layout();
    void paint_gutter(HDC dc);
    int gutter_width() const;

    Host& host_;
    HWND panel_ = nullptr;
    HWND gutter_ = nullptr;
    HWND edit_ = nullptr;
    HWND list_ = nullptr;
    HFONT code_ = nullptr;
    HFONT ui_ = nullptr;
    int dpi_ = 96;
    bool setting_ = false;
    std::vector<Diagnostic> diagnostics_;
};

}  // namespace NiTCAD::app
