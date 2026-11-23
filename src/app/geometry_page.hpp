// The Geometry tab (Unit 26k slice 2; DECISIONS.md GD6-GD8): a Direct2D canvas over the editor's model
// (NiTCAD/app/geometry.hpp), its tools, and the Geometry section shown in the inspector's place.
//
// - The canvas: rulers in µm, the grid, each region's shapes filled by material in layer order (a later layer
//   covers an earlier one, as the deck paints them), outlines, contacts as thick segments with their names, the
//   selection with its handles (a box's eight, a polygon's vertices; the background's: the domain's), rubber-band
//   selection, zoom at the cursor (wheel), pan (middle drag, Space+drag, the Pan tool), Fit.
// - Tools, each with a key while the canvas has the focus: Select (V), Rectangle (R), Polygon (P: click the
//   vertices, Shift keeping one horizontal or vertical from the last; double-click or Enter closes; Ctrl+drag moves
//   a vertex alone, slanting its edges, slice 4), Contact (C: drag along
//   a side of the domain), Pan (H), Measure (M). Arrows nudge by a grid step (Shift: ×10); Delete; Escape cancels.
// - Every edit goes to the deck at once, by minimal edits (VQ13: the deck stays the only source of truth); undo
//   and redo step through the deck texts the page wrote, and start again when the deck is edited elsewhere.
//
// - Live validation (slice 3, GD9): the window hands the page each live check of the deck (the build's own
//   messages); a problem of the text shown is drawn where it lies (a shape, a contact, a point, the domain), with
//   its message along the canvas's foot; Run stays off on this tab while the layout does not build. Mesh preview
//   overlays the mesh the last check built, while it is of the text shown.
//
// - Unit 26m slice 0 (MT1-MT3, MQ6): the group switcher and the floating bar, Camera free / fixed (L: free, a left
//   drag on empty background pans, Shift+drag keeps the rubber band; fixed, no drag or wheel moves the view), the
//   deck line a drag will write in the status bar, its value in the chip at the cursor, typed digits setting it.
// - Unit 26m slice 1 (MT4, MT5): Fill / Replace and the inspector's overlap switch; Copy (or Alt+drag), Mirror
//   (left-right about the centre; Ctrl: up-down; a side handle: across that side; Shift: a quarter turn), Array
//   (drag: pitch; wheel: count) and Align (the side nearest the click). Each is one deck edit and one undo step.
//
// The tests drive the canvas by its messages (GE8; no screen automation).
#pragma once

#include <windows.h>

#include <chrono>
#include <expected>
#include <optional>
#include <string>
#include <vector>

#include <span>
#include <utility>

#include "NiTCAD/app/document.hpp"
#include "NiTCAD/app/geometry.hpp"
#include "dimensions_panel.hpp"
#include "modeller_bar.hpp"

#include <memory>

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
        virtual void geometry_sweep(const std::string& /*name*/) {}  // Unit 26n slice 3: a dimension row's Sweep...

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

    // --- Unit 26m slice 0 (MT1-MT3, MQ6) ---
    void set_camera_fixed(bool fixed);
    [[nodiscard]] bool camera_fixed() const noexcept { return camera_fixed_; }
    void set_group(modeller::Group g);
    [[nodiscard]] modeller::Group group() const noexcept { return bar_->group(); }
    [[nodiscard]] const ModellerBar& bar() const noexcept { return *bar_; }
    [[nodiscard]] modeller::Placement bar_placement() const;              // canvas pixels
    [[nodiscard]] std::optional<modeller::Box> selection_bounds() const;  // canvas pixels
    [[nodiscard]] const std::string& live_line() const noexcept { return live_line_; }  // "" when no drag
    [[nodiscard]] std::string chip() const { return bar_->chip_text(); }
    [[nodiscard]] const modeller::Entry& typed() const noexcept { return entry_; }
    // --- Unit 26m slice 1 (MT4, MT5) ---
    void set_mode(modeller::Mode m);  // a Transform mode of the Select tool (modeller.hpp)
    [[nodiscard]] modeller::Mode mode() const noexcept { return mode_; }
    bool toggle_overlap();  // the picked shapes' regions: fill <-> replace
    [[nodiscard]] int array_count() const noexcept { return array_count_; }

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
    // A deck with no mesh: a drag on the canvas draws a new device, no button. A drag makes a 2D one; Ctrl+drag a
    // 3D one, whose height follows the cursor after the release until a click ends it (Esc cancels).
    [[nodiscard]] bool can_start() const;
    // The size being drawn ("2 × 1 µm", "2 × 1 × 0.5 µm"; "" none), shown at the cursor.
    [[nodiscard]] std::string start_readout() const;

    // The inspector's controls' ids (offsets from inspector_id).
    enum Control : int {
        c_heading = 0, c_name_label, c_name, c_material_label, c_material, c_raise, c_lower,
        c_x_label, c_x, c_y_label, c_y, c_w_label, c_w, c_h_label, c_h, c_vertices_label, c_vertices, c_apply,
        c_kind_label, c_kind, c_draw_label, c_draw, c_hint, c_overlap, c_count
    };
    [[nodiscard]] HWND control(int c) const;

    // Unit 26n (DT1, DT2): the dimensions panel in the inspector, the dimension lines on the canvas.
    [[nodiscard]] DimensionsPanel& dimensions_panel() noexcept { return *dims_; }
    [[nodiscard]] const std::vector<dims::Mark>& dimension_marks() const noexcept { return marks_; }
    [[nodiscard]] std::vector<std::array<float, 4>> dimension_lines() const;  // screen segments, x0 y0 x1 y1, one per mark
    [[nodiscard]] bool view_only() const noexcept { return view_only_; }

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
    enum class Drag { none, move, handle, vertex, domain, contact_end, rubber, rectangle, contact, pan, measure, array, corner };
    // Unit 26m slice 2 (MQ6): a corner cut (Chamfer, Round vertex) dragged in from a shape's corner.
    geo::ShapeRef corner_shape_;
    geo::Coord corner_r_ = 0;
    int corner_chords_ = 0;
    Drag drag_ = Drag::none;
    geo::Geometry g_start_;  // the model when the drag began
    geo::Pt anchor_{}, at_{};
    double press_x_ = 0, press_y_ = 0, last_x_ = 0, last_y_ = 0;
    int drag_index_ = -1;
    bool alone_ = false;  // a vertex dragged with Ctrl (slice 4)
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
    void locate_problems();

    // A new device being drawn (can_start).
    enum class Start { none, footprint, height };
    Start start_ = Start::none;
    bool start3d_ = false;
    geo::Pt start_a_{}, start_b_{};
    double start_y_ = 0;  // the cursor's y when the footprint ended
    geo::Coord start_h_ = 0;
    [[nodiscard]] geo::Pt start_point(double x, double y) const;
    void start_finish();
    void prune_selection();

    // Unit 26m slice 0.
    std::unique_ptr<ModellerBar> bar_;
    // Unit 26n: the dimensions panel's host (the deck through this page's undo), the panel, the marks, a line's drag.
    struct DimsHost final : DimensionsPanel::Host {
        GeometryPage* page = nullptr;
        std::string dimensions_deck() override;
        void dimensions_edited(const std::string& text) override;
        void dimensions_selected(const std::string& name) override;
        bool dimensions_sweeps() override { return true; }
        void dimensions_sweep(const std::string& name) override;
        bool dimensions_unlocks() override;            // Unit 26n slice 2 (DT5)
        std::string dimensions_unlock() override;
    };
    DimsHost dims_host_;
    std::unique_ptr<DimensionsPanel> dims_;
    std::vector<dims::Mark> marks_;
    int dim_drag_ = -1;  // the mark being dragged, or -1
    double dim_x_ = 0.0, dim_y_ = 0.0;
    std::string dim_literal_;
    bool view_only_ = false;  // a parametric deck drawn from its substituted text: dimensions only, no shape edits
    [[nodiscard]] std::expected<geo::Geometry, std::string> read_page(const std::string& text);
    void reload_dimensions();
    bool dimension_down(double x, double y);
    bool dimension_move(double x, double y);
    bool dimension_up();
    bool camera_fixed_ = false;  // MQ2: free on open
    modeller::Entry entry_;
    int typed_axis_ = -1;  // the axis (0 x, 1 y) a typed value sets
    geo::Pt drag_d_{};     // the drag's offset now
    std::string live_line_;
    int canvas_top_ = 0;
    void after_drag_step(double x, double y);
    [[nodiscard]] std::optional<geo::Geometry> pending() const;
    [[nodiscard]] std::string chip_value() const;
    [[nodiscard]] std::string fresh_contact() const;
    [[nodiscard]] geo::Coord typed_or(geo::Coord dragged) const;
    bool on_char(wchar_t c);
    void place_bar();
    void end_drag_ui();

    // Unit 26m slice 1.
    modeller::Mode mode_ = modeller::Mode::none;
    int array_count_ = 3;
    std::optional<geo::Geometry> copy_from_;  // a copy drag: the geometry before the copies were made
    std::vector<geo::ShapeRef> copy_sel_;
    bool transform_click(const geo::Hit& hit, geo::Pt at, bool shift, bool ctrl);  // Mirror and Align
};

}  // namespace NiTCAD::app
