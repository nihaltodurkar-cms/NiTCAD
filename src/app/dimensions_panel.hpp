// The dimensions panel (Unit 26n, DT1, DECISIONS.md): a deck's [parameters] as labelled rows in the Geometry tab's
// inspector. A row's value typed and Enter (or the mouse wheel over it) writes its one [parameters] line through the
// host (the deck editor's undo); a refused value shows under its row and is not written. While the deck does not
// parse, the rows stay, disabled, with the message. The model is app/dimensions.hpp.
#pragma once

#include "NiTCAD/app/dimensions.hpp"

#include <windows.h>

#include <string>
#include <vector>

namespace NiTCAD::app {

class DimensionsPanel {
public:
    struct Host {
        [[nodiscard]] virtual std::string dimensions_deck() = 0;           // the deck now
        virtual void dimensions_edited(const std::string& text) = 0;      // through the deck editor (undo)
        virtual void dimensions_selected(const std::string& name) = 0;    // a row picked: its line lit
        // Unit 26n slice 1: why the deck with an edit is not taken (shown on the row, not written); empty: taken.
        [[nodiscard]] virtual std::string dimensions_refuse(const std::string& /*text*/) { return {}; }

    protected:
        ~Host() = default;
    };

    DimensionsPanel(HWND parent, Host& host, int id);
    ~DimensionsPanel();
    DimensionsPanel(const DimensionsPanel&) = delete;
    DimensionsPanel& operator=(const DimensionsPanel&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return hwnd_; }
    void reload();                          // re-read the rows from the deck
    void select(const std::string& name);   // light a row (its dimension line was clicked)
    void set_group(std::string group);       // only that group's rows (Unit 26n slice 1, the New device dialog); empty: all
    [[nodiscard]] const std::string& selected() const noexcept { return selected_; }
    [[nodiscard]] int preferred_height(int width) const;
    void set_dpi(int dpi, HFONT font, HFONT bold);
    void layout();

    // The rows' state, and the actions the messages take (the tests drive them too).
    [[nodiscard]] std::size_t rows() const noexcept { return rows_.size(); }
    [[nodiscard]] const std::vector<dims::Dimension>& dimensions() const noexcept { return rows_; }
    [[nodiscard]] std::string row_error(std::size_t row) const;
    [[nodiscard]] const std::string& message() const noexcept { return message_; }
    [[nodiscard]] HWND row_edit(std::size_t row) const;
    bool commit(std::size_t row, const std::string& typed);  // as Enter in the row's edit
    bool wheel(std::size_t row, int notches);

private:
    struct Row {
        HWND header = nullptr;  // a group's title above the row, or null
        HWND label = nullptr, edit = nullptr, error = nullptr;
        std::string error_text;
    };
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK edit_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l, UINT_PTR id, DWORD_PTR self);
    void rebuild();  // the child windows for rows_
    void show_rows();
    [[nodiscard]] std::size_t row_of(HWND edit) const;

    Host& host_;
    HWND hwnd_ = nullptr;
    HWND message_label_ = nullptr;
    int dpi_ = 96;
    HFONT font_ = nullptr, bold_ = nullptr;
    std::vector<dims::Dimension> rows_;
    std::vector<Row> controls_;
    std::string message_, selected_, group_;
    bool enabled_ = true;
};

}  // namespace NiTCAD::app
