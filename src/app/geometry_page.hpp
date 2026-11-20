// The Geometry tab (Unit 26k slice 2; DECISIONS.md GD6-GD8): a Direct2D canvas over the editor's model
// (NiTCAD/app/geometry.hpp), its tools, and the Geometry section shown in the inspector's place.
//
// - The canvas: rulers in µm, the grid, each region's shapes filled by material in layer order (a later layer
//   covers an earlier one, as the deck paints them), outlines, contacts as thick segments with their names, the
//   selection with its handles (a box's eight, a polygon's vertices; the background's: the domain's), rubber-band
//   selection, zoom at the cursor (wheel), pan (middle drag, Space+drag, the Pan tool), Fit.
// - Tools, each with a key while the canvas has the focus: Select (V), Rectangle (R), Polygon (P: click the
//   vertices, each kept horizontal or vertical from the last; double-click or Enter closes), Contact (C: drag along
//   a side of the domain), Pan (H), Measure (M). Arrows nudge by a grid step (Shift: ×10); Delete; Escape cancels.
// - Every edit goes to the deck at once, by minimal edits (VQ13: the deck stays the only source of truth); undo
//   and redo step through the deck texts the page wrote, and start again when the deck is edited elsewhere.
//
// - Live validation (slice 3, GD9): the window hands the page each live check of the deck (the build's own
//   messages); a problem of the text shown is drawn where it lies (a shape, a contact, a point, the domain), with
//   its message along the canvas's foot; Run stays off on this tab while the layout does not build. Mesh preview
//   overlays the mesh the last check built, while it is of the text shown.
//
// The tests drive the canvas by its messages (GE8; no screen automation).
#pragma once

#include <windows.h>

#include <chrono>
#include <optional>
#include <string>
#include <vector>

#include <span>
#include <utility>

#include "NiTCAD/app/document.hpp"
#include "NiTCAD/app/geometry.hpp"

struct ID2D1Factory;
struct ID2D1HwndRenderTarget;
struct IDWriteFactory;
struct IDWriteTextFormat;

namespace NiTCAD::app {

class GeometryPage {
public:
    struct Host {
        [[nodiscard]] virtual std::string geometry_deck() = 0;       // the deck's text now
        virtual void geometry_edited(const std::string& text) = 0;    // the deck after an edit, undo or redo
        [[nodiscard]] virtual std::optional<geo::Rect> geometry_bounds() = 0;  // the built device's extent
        virtual void geometry_status(const std::string& cursor) = 0;  // the cursor's position, for the status bar
        virtual void geometry_selected(const std::string& region) = 0;  // a region selected on the canvas

    protected:
        ~Host() = default;
    };

    enum class Tool : int { select, rectangle, polygon, contact, pan, measure };

    // `software`: draw with Direct2D's software rasterizer (the WARP adapter's counterpart; the tests).
    GeometryPage(HWND parent, Host& host, int id, int inspector_id, bool software);
    ~GeometryPage();
    GeometryPage(const GeometryPage&) = delete;
    GeometryPage& operator=(const GeometryPage&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return panel_; }
    [[nodiscard]] HWND canvas() const noexcept { return canvas_; }
    [[nodiscard]] HWND inspector() const noexcept { return inspector_; }  // shown in the inspector's place
    void set_dpi(int dpi, HFONT font, HFONT bold);
    void layout_inspector();

    // Reads the deck again (the tab shown, the deck edited elsewhere); `fit` frames the domain.
    void reload(bool fit = false, bool force = false);
    void fit();
    void set_tool(Tool t);
    void set_view(const geo::View& v);  // the tests: a known zoom and origin
    [[nodiscard]] Tool tool() const noexcept { return tool_; }
    bool undo();
    bool redo();
    bool copy();
    bool paste();
    bool duplicate();
    bool remove();
    // Selects a region's first shape (the tree's selection); false when the deck has none of that name.
    bool select_region(const std::string& name);
    // An inspector box's or button's value taken (Enter, a button, a choice); false when not one of its controls.
    bool commit(int control);

    // State, for the tests.
    [[nodiscard]] const geo::Geometry& geometry() const noexcept { return g_; }
    [[nodiscard]] const geo::View& view() const noexcept { return view_; }
    [[nodiscard]] const std::vector<geo::ShapeRef>& selection() const noexcept { return selection_; }
    [[nodiscard]] int selected_contact() const noexcept { return contact_; }
    [[nodiscard]] std::string material() const { return material_; }
    void set_material(std::string m) { material_ = std::move(m); }
    [[nodiscard]] bool can_undo() const noexcept { return history_.can_undo(); }
    [[nodiscard]] bool can_redo() const noexcept { return history_.can_redo(); }
    [[nodiscard]] const std::string& message() const noexcept { return message_; }  // why the deck is not drawn
    // Draws now (not waiting for WM_PAINT), its time in ms; -1 when nothing could be drawn.
    double render();
    [[nodiscard]] std::string ruler_text(bool x_axis) const;  // the labels the last frame's ruler drew

    // --- live validation and the mesh preview (slice 3, GD9) ---
    // A live check of `text` ended: its diagnostics (empty: the deck builds). Kept while `text` is the text shown.
    void checked(const std::string& text, std::span<const Diagnostic> diagnostics);
    enum class Check { pending, valid, invalid };
    // The state of the text shown: not checked yet, builds, or does not.
    [[nodiscard]] Check check() const noexcept;
    [[nodiscard]] const std::vector<geo::Problem>& problems() const noexcept { return problems_; }
    // Why Run is off on this tab ("" when it is not): the layout's first problem.
    [[nodiscard]] std::string run_blocked() const;
    // The mesh a live check of `text` built (points in quanta, edges as point indices).
    struct MeshOverlay {
        std::vector<geo::Pt> points;
        std::vector<std::pair<int, int>> edges;
    };
    void set_mesh(const std::string& text, MeshOverlay mesh);
    void set_mesh_preview(bool on);
    [[nodiscard]] bool mesh_preview() const noexcept { return mesh_preview_; }
    // The overlay drawn now: the preview on and the mesh of the text shown.
    [[nodiscard]] bool mesh_shown() const noexcept { return mesh_preview_ && !mesh_text_.empty() && mesh_text_ == text_seen_; }
    [[nodiscard]] const MeshOverlay& mesh() const noexcept { return mesh_; }
    [[nodiscard]] std::string banner() const;  // the canvas's foot line drawn last ("" none)

    // The inspector's controls' ids (offsets from inspector_id).
    enum Control : int {
        c_heading = 0, c_name_label, c_name, c_material_label, c_material, c_raise, c_lower,
        c_x_label, c_x, c_y_label, c_y, c_w_label, c_w, c_h_label, c_h, c_vertices_label, c_vertices, c_apply,
        c_kind_label, c_kind, c_draw_label, c_draw, c_hint, c_new_layout, c_count
    };
    [[nodiscard]] HWND control(int c) const;

private:
    static LRESULT CALLBACK panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK canvas_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK inspector_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    LRESULT on_canvas(UINT message, WPARAM w, LPARAM l);
    void layout();
    void on_down(double x, double y, bool shift, bool ctrl, bool double_click);
    void on_move(double x, double y, bool shift);
    void on_up(double x, double y, bool shift);
    bool on_key(UINT key, bool shift, bool ctrl);
    void cancel();
    // The model written to the deck (an edit's end): true when the deck changed.
    bool commit_model(geo::Geometry next);
    void apply_text(const std::string& text);
    void update_inspector();
    void update_tool_buttons();
    void invalidate();
    bool ensure_target();
    void release_target();
    [[nodiscard]] double ruler_px() const noexcept;
    [[nodiscard]] geo::Pt snapped(double x, double y, bool shift, std::optional<geo::ShapeRef> skip = {}) const;

    Host& host_;
    HWND panel_ = nullptr, strip_ = nullptr, canvas_ = nullptr, inspector_ = nullptr;
    std::vector<HWND> tool_buttons_;
    std::vector<HWND> controls_;
    int inspector_id_ = 0;
    int dpi_ = 96;
    HFONT font_ = nullptr, bold_ = nullptr;
    bool software_ = false;

    ID2D1Factory* d2d_ = nullptr;
    ID2D1HwndRenderTarget* target_ = nullptr;
    IDWriteFactory* dwrite_ = nullptr;
    IDWriteTextFormat* text_ = nullptr;

    geo::Geometry g_;
    std::string text_seen_;  // the deck text g_ was read from
    std::string message_;
    geo::EditHistory history_;
    geo::View view_;
    bool fitted_ = false;
    Tool tool_ = Tool::select;
    std::string material_ = "silicon";
    std::vector<geo::ShapeRef> selection_;
    int contact_ = -1;
    std::optional<geo::Shape> clipboard_;
    std::string clipboard_material_;

    // A drag in progress.
    enum class Drag { none, move, handle, vertex, domain, contact_end, rubber, rectangle, contact, pan, measure };
    Drag drag_ = Drag::none;
    geo::Geometry g_start_;  // the model when the drag began
    geo::Pt anchor_{}, at_{};
    double press_x_ = 0, press_y_ = 0, last_x_ = 0, last_y_ = 0;
    int drag_index_ = -1;
    bool space_ = false;
    std::vector<geo::Pt> polygon_;  // the polygon being drawn
    std::optional<geo::Pt> hover_;   // the polygon tool's next point
    mutable std::string ruler_x_, ruler_y_;

    // Slice 3: the last live check (its text, problems), the mesh it built, the preview toggle.
    std::string checked_text_;
    std::vector<Diagnostic> diagnostics_;
    std::vector<geo::Problem> problems_;  // diagnostics_ placed on g_, while checked_text_ is the text shown
    std::string mesh_text_;
    MeshOverlay mesh_;
    bool mesh_preview_ = false;
    HWND mesh_button_ = nullptr;
    void locate_problems();
    void prune_selection();
};

}  // namespace NiTCAD::app
