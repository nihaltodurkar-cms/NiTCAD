// The deck pane (Unit 26, DECISIONS.md VA2): the deck's text in a Rich Edit control with line numbers
// beside it, and below it the check's errors, each taking the caret to its line.
//
// Unit 26b (XD1-XD3): the lines with errors underlined in the text; completion from the deck's
// vocabulary in a list at the caret, opened by Ctrl+Space or by typing, moved by Up and Down,
// accepted by Enter, Tab or a double click, closed by Escape.
#pragma once

#include <windows.h>

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/deck_vocabulary.hpp"
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
    [[nodiscard]] HWND completion_list() const noexcept { return completion_; }

    // The text with '\n' line ends (UTF-8). Setting it does not notify the host.
    [[nodiscard]] std::string text() const;
    void set_text(std::string_view utf8);
    void set_diagnostics(const std::vector<Diagnostic>& d);
    [[nodiscard]] const std::vector<Diagnostic>& diagnostics() const noexcept { return diagnostics_; }
    // Selects the line (1-based) and scrolls it into view; the editor takes the focus.
    void go_to_line(int line);
    [[nodiscard]] int caret_line() const;
    // The caret's line (1-based) and its byte column in the line's UTF-8.
    [[nodiscard]] std::pair<int, std::size_t> caret() const;
    // Puts the caret at a line (1-based) and byte column.
    void set_caret(int line, std::size_t column);
    void set_dpi(int dpi);

    // Completion (XD3): opens the list for the caret (forced: even with nothing typed); false when
    // nothing fits. accept_completion replaces the word at the caret with the item's text.
    bool show_completion(bool forced);
    void hide_completion();
    [[nodiscard]] bool completion_visible() const noexcept { return completion_shown_.has_value(); }
    [[nodiscard]] const std::optional<Completion>& completion() const noexcept { return completion_shown_; }
    [[nodiscard]] int completion_selection() const;
    void accept_completion(int index);

private:
    static LRESULT CALLBACK panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK gutter_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK edit_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l, UINT_PTR id, DWORD_PTR data);
    void layout();
    void paint_gutter(HDC dc);
    void paint_underlines();
    int gutter_width() const;
    [[nodiscard]] std::wstring line_text(int line) const;  // 1-based, without its end
    [[nodiscard]] LONG char_of(int line, std::size_t column) const;
    void after_typing(wchar_t typed);

    Host& host_;
    HWND panel_ = nullptr;
    HWND gutter_ = nullptr;
    HWND edit_ = nullptr;
    HWND list_ = nullptr;
    HWND completion_ = nullptr;
    HFONT code_ = nullptr;
    HFONT ui_ = nullptr;
    int dpi_ = 96;
    bool setting_ = false;
    bool swallow_char_ = false;  // the WM_CHAR of a key the completion list took
    std::vector<Diagnostic> diagnostics_;
    std::optional<Completion> completion_shown_;
};

}  // namespace NiTCAD::app
