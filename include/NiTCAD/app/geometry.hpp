// The geometry editor's model (Unit 26k, slice 1: GD3-GD5). No Win32: the canvas (slice 2) draws it and sends it
// the user's edits.
//
// The deck stays the only source of truth (VQ13): read_geometry takes a 2D device's regions, shapes and contacts
// from the deck's text; write_geometry gives the text with only the lines that changed rewritten (the regions'
// box and polygon lines, their headers and materials, the [mesh] patch, layered and domain lines, the contacts'
// sections), every other line, comment and blank kept byte for byte. Layers (DECISIONS.md Unit 26k, GD1): a region
// is painted in deck order, each of its shapes in turn, each replacing what it covers; the one region without
// shapes is the background. A region's shapes are one section, so the layer order is the regions' order.
//
// Coordinates are integers in units of 0.1 nm (exact: a decimal deck value round-trips; snapping and hit tests
// compare integers). y grows upward, as in the deck.
#pragma once

#include <compare>
#include <cstdint>
#include <expected>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace NiTCAD::app::geo {

using Coord = std::int64_t;                  // 0.1 nm
inline constexpr Coord per_um = 10000;        // quanta per µm
inline constexpr double quantum_cm = 1e-8;

struct Pt {
    Coord x = 0, y = 0;
    friend auto operator<=>(const Pt&, const Pt&) = default;
};

struct Rect {
    Pt lo, hi;  // lo.x <= hi.x, lo.y <= hi.y
    friend bool operator==(const Rect&, const Rect&) = default;
    [[nodiscard]] bool contains(Pt p) const noexcept { return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y; }
};
[[nodiscard]] Rect normalized(Pt a, Pt b) noexcept;

enum class MeshKind : unsigned char {
    tensor,       // a tensor deck: boxes only, the domain the mesh's (given by the caller)
    tiled,        // generate2d with polygons that tile the domain (editing writes it layered)
    layered,      // generate2d with layered = on
    unsupported,  // 1D, 3D, a mesh file, ...: shown, not edited (Geometry::why)
};

enum class ShapeKind : unsigned char { box, polygon };

struct Shape {
    ShapeKind kind = ShapeKind::box;
    std::vector<Pt> points;  // box: lo and hi; polygon: its vertices
    friend bool operator==(const Shape&, const Shape&) = default;
    [[nodiscard]] Rect bounds() const noexcept;
    [[nodiscard]] std::vector<Pt> outline() const;  // a box's four corners counter-clockwise, or the polygon
};

struct Region {
    std::string name, material;
    std::vector<Shape> shapes;  // empty: the background
    int origin = -1;            // the deck's region section it came from (its index among them), -1 new
    friend bool operator==(const Region&, const Region&) = default;
};

// A contact on the domain's boundary: a segment (a generated deck's [mesh] patch), or a tensor deck's face with a
// box along it. `editable` false for what the editor cannot rewrite (a face patch without a box covering the
// whole face is shown as the whole face, still editable; an unknown form is not).
struct Contact {
    std::string name, kind;
    Pt a, b;
    bool editable = true;
    int origin = -1;  // the deck's contact section it came from (its index among them), -1 new
    friend bool operator==(const Contact&, const Contact&) = default;
};

struct Geometry {
    MeshKind mesh = MeshKind::unsupported;
    std::string why;           // unsupported: why
    Rect domain;
    bool domain_given = false;  // a layered deck's `domain` line
    std::vector<Region> regions;
    std::vector<Contact> contacts;
    friend bool operator==(const Geometry&, const Geometry&) = default;
    [[nodiscard]] bool editable() const noexcept { return mesh != MeshKind::unsupported; }
    [[nodiscard]] int background() const noexcept;  // the region without shapes, or -1
};

// The deck's geometry. `bounds` is the built device's extent (a tensor deck's domain, which its axes give; the
// shapes' bounds otherwise). Errors: the deck does not parse.
[[nodiscard]] std::expected<Geometry, std::string> read_geometry(std::string_view deck, std::optional<Rect> bounds = {});

// The deck's text with `g` written into it by minimal edits (see above). Errors: the deck does not parse, or `g`
// is not editable.
[[nodiscard]] std::expected<std::string, std::string> write_geometry(std::string_view deck, const Geometry& g);

// A length as the deck writes it: µm with up to four decimals ("0.95um", "2um", "0").
[[nodiscard]] std::string format_length(Coord q);
// A deck length (a number, optionally cm, um or nm) in quanta, rounded.
[[nodiscard]] std::optional<Coord> parse_length(std::string_view text);

// --- edits (GD4): each changes `g` and says whether it could ---------------------------------------------------

struct ShapeRef {
    int region = -1, shape = -1;
    friend bool operator==(const ShapeRef&, const ShapeRef&) = default;
};

// A new region on top with one shape, named "<material>_<n>" (unique), its index returned (-1 refused: a tensor
// deck takes no polygon; a polygon has four or more vertices and axis-parallel edges until slice 4). Any edit of a
// tiled deck first makes it layered, refused when a polygon has a slanted edge.
int add_box(Geometry& g, Rect r, std::string material);
int add_polygon(Geometry& g, std::vector<Pt> points, std::string material);
bool move_shape(Geometry& g, ShapeRef s, Coord dx, Coord dy);
bool set_box(Geometry& g, ShapeRef s, Rect r);
bool set_vertex(Geometry& g, ShapeRef s, int k, Pt p);  // its neighbours follow along their shared edges
bool insert_vertex(Geometry& g, ShapeRef s, int after, Pt p);  // a point on the edge after vertex `after`
bool remove_vertex(Geometry& g, ShapeRef s, int k);  // one between collinear edges; four kept
bool delete_shape(Geometry& g, ShapeRef s);          // a region left without shapes goes too
bool rename_region(Geometry& g, int region, std::string name);  // unique, a deck identifier
bool set_material(Geometry& g, int region, std::string material);
bool raise_region(Geometry& g, int region);  // one layer up (later in the deck)
bool lower_region(Geometry& g, int region);
int duplicate_shape(Geometry& g, ShapeRef s, Coord dx, Coord dy);  // a new region on top, its index
bool mirror_shape(Geometry& g, ShapeRef s, bool horizontal);       // about its bounds' centre
// Contacts: a segment along one side of the domain (refused otherwise).
int add_contact(Geometry& g, std::string name, std::string kind, Pt a, Pt b);
bool set_contact(Geometry& g, int contact, Pt a, Pt b);
bool delete_contact(Geometry& g, int contact);
bool rename_contact(Geometry& g, int contact, std::string name);
bool set_contact_kind(Geometry& g, int contact, std::string kind);
// A layered deck's domain (its `domain` line); each contact moves with its side, clipped (refused if one empties).
bool set_domain(Geometry& g, Rect domain);
[[nodiscard]] bool on_one_side(const Rect& domain, Pt a, Pt b) noexcept;

// --- undo (GD4): the deck's text before each edit -------------------------------------------------------------

class EditHistory {
public:
    static constexpr std::size_t limit = 1000;
    void record(std::string before);  // before an edit; clears redo
    [[nodiscard]] std::optional<std::string> undo(std::string current);
    [[nodiscard]] std::optional<std::string> redo(std::string current);
    [[nodiscard]] bool can_undo() const noexcept { return !undo_.empty(); }
    [[nodiscard]] bool can_redo() const noexcept { return !redo_.empty(); }
    void clear() noexcept;

private:
    std::vector<std::string> undo_, redo_;
};

// --- the view, snapping and hit tests (GD5) --------------------------------------------------------------------

// Model to screen: x right, y down on screen; `origin_x/y` the screen position of model (0, 0).
struct View {
    double px_per_um = 100.0;
    double origin_x = 0.0, origin_y = 0.0;
    int dpi = 96;
    [[nodiscard]] double sx(Coord x) const noexcept { return origin_x + px_per_um * static_cast<double>(x) / per_um; }
    [[nodiscard]] double sy(Coord y) const noexcept { return origin_y - px_per_um * static_cast<double>(y) / per_um; }
    [[nodiscard]] Pt model(double x, double y) const noexcept;
    [[nodiscard]] double px(double at_96) const noexcept { return at_96 * dpi / 96.0; }
};

// The grid step: the smallest 1, 2 or 5 × 10ⁿ µm (at least one quantum) at least 12 px (at 96 dpi) apart.
[[nodiscard]] Coord grid_step(const View& v) noexcept;

enum class SnapKind : unsigned char { none, grid, edge, midpoint, vertex };
struct Snap {
    Pt p;
    SnapKind kind = SnapKind::none;
};
// The point snapped within 8 px (at 96 dpi): to a vertex (the shapes', the contacts' ends, the domain's corners),
// else an edge's midpoint, else the nearest point of an edge or the domain's boundary, else the grid (when on).
// `skip` excludes a shape being dragged.
[[nodiscard]] Snap snap(const Geometry& g, const View& v, Pt raw, std::optional<ShapeRef> skip = {}, bool grid = true);
// The point kept horizontal or vertical from `anchor`, along its longer offset (Shift).
[[nodiscard]] Pt constrain(Pt anchor, Pt p) noexcept;

enum class HitKind : unsigned char { none, handle, vertex, edge, contact_end, contact, shape };
struct Hit {
    HitKind kind = HitKind::none;
    ShapeRef shape;
    int index = -1;  // the handle (0..7 counter-clockwise from the lower left), vertex, edge, contact or end (0, 1)
    int end = -1;
    friend bool operator==(const Hit&, const Hit&) = default;
};
// What lies under a screen point: the selected shape's handles (a box) or vertices (a polygon), within 4.5 px; a
// contact's end, then the contact, within 5 px; the selected polygon's edges, within 4 px; else the topmost shape
// whose outline contains the point, else the background inside the domain.
[[nodiscard]] Hit hit_test(const Geometry& g, const View& v, double x, double y, std::optional<ShapeRef> selected = {});

// The box with handle `h` dragged to `p` (the opposite side fixed), normalized.
[[nodiscard]] Rect drag_handle(const Rect& r, int h, Pt p) noexcept;

}  // namespace NiTCAD::app::geo
