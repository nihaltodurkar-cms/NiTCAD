// The 3D solid modeller's model (Unit 26l, slice 1: SD3). No Win32: the viewport (slice 2) draws it and sends it the
// user's edits.
//
// The deck stays the only source of truth: read_model takes a 3D device's regions, solids and contacts from the
// deck's text; write_model gives the text with only the lines that changed rewritten (the regions' box and prism
// lines, their headers and materials, the [mesh] type, layered, domain and patch lines, the contacts' sections), every
// other line, comment and blank kept byte for byte. Layers (SD1): a region's solids are painted in deck order, each
// replacing what it covers; the one region without solids is the background. A region's solids are one section, so
// the layer order is the regions' order.
//
// An extrude3d deck (SQ5) reads as one prism per cell and `extrude` line, exact, along the extrusion's axis. It is
// not edited in place: convert_to_layered makes it a layered generate3d model, and the next write rewrites its
// [mesh] and its regions' extrude lines.
//
// Coordinates are integers in units of 0.1 nm, as the 2D editor's (geo::Coord). A prism's polygon is in (u, v), the
// coordinates along axes (axis + 1) % 3 and (axis + 2) % 3 (axis x: (y, z); y: (z, x); z: (x, y)), as the deck's.
#pragma once

#include <compare>
#include <cstdint>
#include <expected>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "NiTCAD/app/geometry.hpp"
#include "NiTCAD/mesh/complex.hpp"

namespace NiTCAD::app::solid {

using geo::Coord;
using geo::EditHistory;
using geo::per_um;
using geo::Pt;

struct P3 {
    Coord x = 0, y = 0, z = 0;
    friend auto operator<=>(const P3&, const P3&) = default;
    [[nodiscard]] Coord operator[](int axis) const noexcept { return axis == 0 ? x : axis == 1 ? y : z; }
    [[nodiscard]] Coord& operator[](int axis) noexcept { return axis == 0 ? x : axis == 1 ? y : z; }
};

struct Box3 {
    P3 lo, hi;  // lo <= hi on each axis
    friend bool operator==(const Box3&, const Box3&) = default;
    [[nodiscard]] bool contains(P3 p) const noexcept {
        return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z;
    }
};
[[nodiscard]] Box3 normalized(P3 a, P3 b) noexcept;

// A point of the plane across `axis` at `t`, from its (u, v); and a point's (u, v).
[[nodiscard]] P3 from_uv(int axis, Pt uv, Coord t) noexcept;
[[nodiscard]] Pt uv_of(int axis, P3 p) noexcept;

enum class MeshKind : unsigned char {
    layered,      // generate3d with layered = on
    extruded,     // extrude3d: read as prisms, edited after convert_to_layered (SQ5)
    unsupported,  // 1D, 2D, a tensor or unlayered generate3d mesh, a mesh file, ...: Model::why
};

enum class SolidKind : unsigned char { box, prism };

struct Solid {
    SolidKind kind = SolidKind::box;
    Box3 box;                 // a box
    int axis = 2;             // a prism: its axis, extent and polygon in (u, v)
    Coord from = 0, to = 0;
    std::vector<Pt> polygon;
    // Unit 26m MT8: the deck's references by coordinate, in the line's order (a box's x0 x1 y0 y1 z0 z1; a prism's
    // from, to, then u v of each vertex): "" for a value, else the reference as written (`$Lg`, `$x0 + $Lg`). Empty
    // when none. A coordinate whose reference no longer gives it is written as its value (the tie broken).
    std::vector<std::string> refs;
    friend bool operator==(const Solid&, const Solid&) = default;
    [[nodiscard]] Box3 bounds() const noexcept;
    // Corners: a box's eight (bit k of the index: hi on axis k); a prism's polygon at `from`, then at `to`.
    [[nodiscard]] std::vector<P3> vertices() const;
    // Edges as pairs of vertices' indices.
    [[nodiscard]] std::vector<std::pair<int, int>> edges() const;
};

struct Region {
    std::string name, material;
    std::vector<Solid> solids;  // empty: the background
    int origin = -1;            // the deck's region section it came from (its index among them), -1 new
    bool fill = false;          // Unit 26m (MT4): `overlap = fill`, its shapes taking only what no earlier shape covers
    friend bool operator==(const Region&, const Region&) = default;
};

// A contact: the boundary facets inside its patch box (the deck's `patch <name> = x0 y0 z0 x1 y1 z1`). `editable`
// false for a contact whose patch the modeller cannot read (it is kept as it is).
struct Contact {
    std::string name, kind;
    Box3 box;
    bool editable = true;
    int origin = -1;  // the deck's contact section it came from (its index among them), -1 new
    friend bool operator==(const Contact&, const Contact&) = default;
};

// A named parameter (Unit 26m MT8, the deck's [parameters]): its expression as written and what it resolves to.
struct Param {
    std::string name, expr;
    std::string text;     // what it resolves to, written out (deck::Parameter::text)
    bool length = false;  // a length (value in quanta), else a number
    Coord value = 0;
    double number = 0;    // a number's value
    friend bool operator==(const Param&, const Param&) = default;
};

struct Model {
    MeshKind mesh = MeshKind::unsupported;
    std::string why;          // unsupported: why
    std::string convert_why;  // extruded: why it cannot become layered (empty: it can)
    bool swept = false;       // extruded: a swept extrusion (sweep = on), whose mesh keys the conversion drops
    Box3 domain;
    std::vector<Region> regions;
    std::vector<Contact> contacts;
    // Unit 26m MT8: the deck's [parameters] in its order, and the domain line's references (as Solid::refs, x0 x1 y0
    // y1 z0 z1).
    std::vector<Param> parameters;
    std::vector<std::string> domain_refs;
    friend bool operator==(const Model&, const Model&) = default;
    [[nodiscard]] bool editable() const noexcept { return mesh == MeshKind::layered; }
    [[nodiscard]] int background() const noexcept;  // the region without solids, or -1
};

// The deck's 3D geometry. Errors: the deck does not parse.
[[nodiscard]] std::expected<Model, std::string> read_model(std::string_view deck);

// The deck's text with `m` written into it by minimal edits (see above). An extruded model not converted gives the
// text back unchanged. Errors: the deck does not parse, or `m` is unsupported.
[[nodiscard]] std::expected<std::string, std::string> write_model(std::string_view deck, const Model& m);

// An extruded model made layered (SQ5): its prisms kept, the domain the cells' bounds over the stations. Returns what
// the conversion changes in the deck, for the button to say; refused (Model::convert_why) when the extrusion leaves
// part of its bounds void, which a layered deck cannot.
[[nodiscard]] std::expected<std::string, std::string> convert_to_layered(Model& m);

// The layers as the partition kernel takes them (lengths in cm; the contacts' boxes as patch boxes by contact name).
[[nodiscard]] mesh::LayeredGeometry3D to_layered(const Model& m);

// --- edits (SD3): each changes `m` and says whether it could; refused on a model not layered (convert first) -------

struct SolidRef {
    int region = -1, solid = -1;
    friend bool operator==(const SolidRef&, const SolidRef&) = default;
};

// Coordinates within ±10 cm. A prism's polygon is simple, of non-zero area, three vertices or more, its edges in any
// direction (slanted faces since slice 4, SD6: the deck build refuses slanted prisms along two axes); a box non-empty
// on each axis; a prism's from < to.
inline constexpr Coord coord_limit = 1'000'000'000;
[[nodiscard]] bool prism_ok(const std::vector<Pt>& polygon) noexcept;

// A new region on top with one solid, named "<material>_<n>" (unique), its index returned (-1 refused).
int add_box(Model& m, Box3 b, std::string material);
int add_prism(Model& m, int axis, Coord from, Coord to, std::vector<Pt> polygon, std::string material);
bool move_solid(Model& m, SolidRef s, P3 delta);
bool set_box(Model& m, SolidRef s, Box3 b);
bool set_extent(Model& m, SolidRef s, Coord from, Coord to);  // a prism's
// A face moved along its normal by `offset` (positive: along the axis): a box's faces 0..5 are x_min, x_max, y_min,
// y_max, z_min, z_max; a prism's 0 and 1 its caps at `from` and `to`, 2 + k the side over the polygon's edge k (from
// vertex k to k + 1), moved with its two vertices. Refused when the solid would be empty or invalid, and for a slanted
// side (its normal off the quanta's lattice: its vertices are moved instead).
bool move_face(Model& m, SolidRef s, int face, Coord offset);
// A prism's vertex; `follow`: its neighbours follow along the u- and v-parallel edges they share with it.
bool set_vertex(Model& m, SolidRef s, int k, Pt uv, bool follow = true);
bool insert_vertex(Model& m, SolidRef s, int after, Pt uv);  // a point exactly on the edge after vertex `after`
bool remove_vertex(Model& m, SolidRef s, int k);              // the polygon left valid
bool delete_solid(Model& m, SolidRef s);                     // a region left without solids goes too
bool rename_region(Model& m, int region, std::string name);  // unique, a deck identifier
bool set_material(Model& m, int region, std::string material);
bool raise_region(Model& m, int region);  // one layer up (later in the deck)
bool lower_region(Model& m, int region);
int duplicate_solid(Model& m, SolidRef s, P3 delta);  // a new region on top, its index
bool mirror_solid(Model& m, SolidRef s, int axis);    // about its bounds' centre across `axis`
// --- Unit 26m slice 1 (MT4, MT5) ---
// The region's overlap rule (refused for the background or no change).
bool set_overlap(Model& m, int region, bool fill);
// Mirrored across the plane across `axis` at twice / 2 quanta (a face's plane, or its centre's).
bool mirror_solid_at(Model& m, SolidRef s, int axis, Coord twice);
// A quarter turn about `axis` (counter-clockwise looking down it) through its bounds' centre, rounded down to the
// quantum when the centre is off it; a prism across another axis becomes one across the turned axis.
bool rotate_solid(Model& m, SolidRef s, int axis);
// `count` - 1 copies, the k-th moved k × pitch, each a new region on top (its overlap rule kept); the first one's
// index (-1 refused: count 2..100, a non-zero pitch, every copy valid).
int array_solid(Model& m, SolidRef s, P3 pitch, int count);
// Each solid moved along the face's axis so its bounds' face (0..5: x_min, x_max, ..., z_max) lies at `at` (refused if
// any is invalid or none would move).
bool align_solids(Model& m, const std::vector<SolidRef>& solids, int face, Coord at);

// --- Unit 26m slice 2 (MT6, MT7, MT12) ---
// An edge (Solid::edges) parallel to a prism's axis (a box's: any edge) as a corner of a prism: the solid written as
// the prism along the edge's axis (a box as its rectangle, counter-clockwise from (lo u, lo v)) and the polygon's
// vertex at the edge. None for a prism's cap edges.
struct EdgeCorner {
    Solid prism;
    int corner = -1;
};
[[nodiscard]] std::optional<EdgeCorner> edge_corner(const Solid& s, int edge);
// The edge's corner chamfered (d along both faces) or filleted (radius r in `segments` chords, geo::fillet_corner),
// a box first written as its prism; a prism's vertex k likewise (Round vertex).
bool chamfer_edge(Model& m, SolidRef s, int edge, Coord d);
bool fillet_edge(Model& m, SolidRef s, int edge, Coord r, int segments);
bool fillet_vertex(Model& m, SolidRef s, int k, Coord r, int segments);
// Taper: the edge's corner moved alone by `delta` in the prism's (u, v), so the walls meeting there slant.
bool taper_edge(Model& m, SolidRef s, int edge, Pt delta);
// A cylinder: a prism along `axis` from `from` to `to` whose polygon is geo::regular_polygon(centre, r, sides).
int add_cylinder(Model& m, int axis, Coord from, Coord to, Pt centre, Coord r, int sides, std::string material);
// A face's outward normal: its axis and sign (+1: along the axis); none for a slanted prism side.
struct Normal {
    int axis = 2, sign = 1;
};
[[nodiscard]] std::optional<Normal> face_normal(const Solid& s, int face);
// The new solid on a face (move_face's numbering) inset by `margin` from its edges and reaching `depth` along its
// outward normal (negative: into the solid). A box face, a u/v-parallel prism side, or a cap: a prism's cap keeps its
// polygon with no margin, and takes a margin only when the polygon is its bounds' rectangle. None when it would be
// empty or invalid.
[[nodiscard]] std::optional<Solid> face_solid(const Solid& s, int face, Coord margin, Coord depth);
// Inset and Extrude (MT6): that solid as a new region on top of `material` (Extrude: margin 0, depth > 0), or, of the
// face's own material, a new solid of its region (no interface between them); the region's index, -1 refused. A margin
// on every side of a face inside the device leaves that face a ring, a holed interface piece the generator refuses.
int inset_face(Model& m, SolidRef s, int face, Coord margin, Coord depth, std::string material);
int extrude_face(Model& m, SolidRef s, int face, Coord depth, std::string material);
// MT12 (2): an offset of `moving` along `axis` kept clear: the nearest offset to `delta` (within 2 × margin of it) at
// which no face of the moved solids' bounds across `axis` lies strictly within `margin` of a face of the domain or of
// another solid's bounds whose bounds come within `margin` of it on the other two axes (`skip`: solids that move
// with them). `delta` when no such offset is near.
// MT12 (1): the largest interface_resolution s (cm) the layout's rules allow, within 0.1 % (mesh::largest_resolution
// over the model's complex with `max_size` cm); +infinity without interfaces; the build's own message when no s
// passes (a rule that does not scale with s) or the layers do not resolve.
[[nodiscard]] std::expected<double, std::string> automatic_resolution(const Model& m, double max_size);
// MT12 (3): the expected node count of the mesh at max_size h and interface_resolution s: the bulk at h, and the
// interfaces' pair sets at hp = min(h, 4 s) on both sides, and the graded rows along their feature lines:
// N = a V / h³ + b A / hp² + c L / hp, a, b and c fitted to built meshes (DECISIONS.md, slice 2). 0 when the layers do not resolve.
struct Extent {
    double volume = 0, interface_area = 0;  // cm³, cm²
    double feature_length = 0;              // cm: edges where an interface bends, meets another, or meets the boundary
};
[[nodiscard]] std::optional<Extent> extent(const Model& m);
[[nodiscard]] double node_estimate(const Extent& e, double max_size, double s);
// --- parameters (Unit 26m MT8) -----------------------------------------------------------------------------------
// A coordinate slot of a solid, in its line's order (Solid::refs); none past the end.
[[nodiscard]] Coord* slot(Solid& s, int k);
[[nodiscard]] int slot_count(const Solid& s);
// `name` set to `expr` (added at the end when new), then every parameter and every tied coordinate resolved again:
// changing Lg moves every face tied to it. Errors (the model unchanged): a bad name, an expression that does not
// resolve (undefined, a cycle, neither a length nor a number), or a tied coordinate it would leave off a valid solid.
[[nodiscard]] std::expected<void, std::string> set_parameter(Model& m, std::string_view name, std::string expr);
// After an edit took `before` to `after` (the same solids, coordinates changed): a moved coordinate tied by a lone
// `$p` sets p (a length written as a value) to it; then one tied by a sum `... + $q` (or `- $q`) moves its last
// parameter q by the rest; every tied coordinate is then resolved again. A tie whose parameter is itself an
// expression, or that disagrees with another moved coordinate, is broken instead (written as its value). Solids added,
// removed or with another slot count keep their refs as they are.
void tie_edits(const Model& before, Model& after);
// MT8's N: the solid's extent between slots `near` and `far` (one axis) named `name`: a new length parameter of that
// extent, or an existing one (the far slot then moves to it), and the far slot written `<near> + $name` (`- $name` for a
// low face: a box's lo, a prism's from; near its
// reference if it has one, else its value). False: a bad name, a number parameter, a slot out of range.
[[nodiscard]] bool name_extent(Model& m, SolidRef s, int near, int far, std::string_view name);

// A [mesh] length key of the deck (max_size, interface_resolution) in cm; none when absent or not a length.
[[nodiscard]] std::optional<double> mesh_length(std::string_view deck, std::string_view key);
// The deck with its [mesh] `key` set to `q` quanta (the line replaced, or added after the section's last entry), every
// other line kept byte for byte. Errors: the deck does not parse or has no [mesh].
[[nodiscard]] std::expected<std::string, std::string> set_mesh_length(std::string_view deck, std::string_view key, Coord q);
[[nodiscard]] Coord keep_clear(const Model& m, const std::vector<Solid>& moving, const std::vector<SolidRef>& skip, int axis,
                               Coord delta, Coord margin);

// Contacts: a patch box inside the domain that reaches a face of it with some area (it selects boundary facets).
[[nodiscard]] bool selects_boundary(const Box3& domain, const Box3& b) noexcept;
int add_contact(Model& m, std::string name, std::string kind, Box3 b);
bool set_contact(Model& m, int contact, Box3 b);
bool move_contact(Model& m, int contact, P3 delta);
bool delete_contact(Model& m, int contact);
bool rename_contact(Model& m, int contact, std::string name);
bool set_contact_kind(Model& m, int contact, std::string kind);
// The domain; each contact clipped to it (refused if one would select nothing).
bool set_domain(Model& m, Box3 domain);

// --- snapping (SD3) ------------------------------------------------------------------------------------------------

// The plane across `normal` at `offset` where new solids are drawn; its (u, v) as a prism's along `normal`.
struct WorkPlane {
    int normal = 2;
    Coord offset = 0;
    friend bool operator==(const WorkPlane&, const WorkPlane&) = default;
};

// The grid step on the work plane: the 2D editor's (1, 2 or 5 × 10ⁿ µm, at least 12 px apart at 96 dpi).
[[nodiscard]] Coord grid_step(double px_per_um, int dpi = 96) noexcept;

enum class SnapKind : unsigned char { none, grid, face, edge, midpoint, vertex };
struct Snap {
    P3 p;
    SnapKind kind = SnapKind::none;
    friend bool operator==(const Snap&, const Snap&) = default;
};
// The point snapped within `tolerance` (quanta; the viewport gives 8 px at the point's depth): to a vertex (the
// solids', the contacts' boxes', the domain's corners), else an edge's midpoint, else the nearest point of an edge,
// else onto a face (a box's or the domain's face, a prism's cap or u/v-parallel side, the point inside it), else the
// work plane's grid (`grid` > 0: its step), else the point as given. `skip` excludes a solid being dragged.
[[nodiscard]] Snap snap(const Model& m, const WorkPlane& plane, P3 raw, Coord tolerance, std::optional<SolidRef> skip = {},
                        Coord grid = 0);

}  // namespace NiTCAD::app::solid

namespace NiTCAD::app::solid {

// --- validation (Unit 26l slice 3, SD5): the deck build's own error, placed on the model -------------------------

// Where a problem is drawn: a solid ({region, -1}: every solid of the region, or the domain for the background), a
// contact, a point the message names, or the whole domain.
struct Marker {
    enum class Kind : unsigned char { solid, contact, point, domain };
    Kind kind = Kind::domain;
    SolidRef solid;
    int contact = -1;
    P3 at{};
    friend bool operator==(const Marker&, const Marker&) = default;
};

struct Problem {
    int line = 0;         // the deck line the build named (0: none)
    std::string message;  // the build's message, unchanged
    std::vector<Marker> markers;
};

// The build's error at `line` of `deck` (read into `m`) placed on `m`, as the 2D editor's: the line's region solid
// (its box, prism or extrude line; the region's solids for another line of it), its contact section, or the contacts
// of its [mesh] patch line; the points "(x, y, z)" in cm the message names, each with the solids whose closed
// bounds hold it; the domain for the [mesh] section's other lines. No markers outside the geometry.
[[nodiscard]] Problem locate_problem(const Model& m, std::string_view deck, int line, std::string message);

// The region a point strictly inside one solid lies in: the last solid in paint_order holding it, else the background
// (-1 when there is none). Points on a solid's boundary count as outside it (the mesh preview asks at tet centroids).
[[nodiscard]] int region_at(const Model& m, P3 p);
// The regions with solids in the order that paints them so the last painted holds each point (Unit 26m MT4): the
// `fill` regions last first, then the others in deck order. A view drawing in this order shows the build's layers.
[[nodiscard]] std::vector<int> paint_order(const Model& m);

}  // namespace NiTCAD::app::solid
