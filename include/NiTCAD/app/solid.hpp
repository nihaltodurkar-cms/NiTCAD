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

struct Model {
    MeshKind mesh = MeshKind::unsupported;
    std::string why;          // unsupported: why
    std::string convert_why;  // extruded: why it cannot become layered (empty: it can)
    bool swept = false;       // extruded: a swept extrusion (sweep = on), whose mesh keys the conversion drops
    Box3 domain;
    std::vector<Region> regions;
    std::vector<Contact> contacts;
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

// Coordinates within ±10 cm. A prism's polygon is simple, of non-zero area, its edges parallel to u or v (slanted
// faces are slice 4, SD6); a box non-empty on each axis; a prism's from < to.
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
// vertex k to k + 1), moved with its two vertices. Refused when the solid would be empty or invalid.
bool move_face(Model& m, SolidRef s, int face, Coord offset);
// A prism's vertex; `follow`: its neighbours follow along the u- and v-parallel edges they share with it.
bool set_vertex(Model& m, SolidRef s, int k, Pt uv, bool follow = true);
bool insert_vertex(Model& m, SolidRef s, int after, Pt uv);  // a point on the edge after vertex `after`
bool remove_vertex(Model& m, SolidRef s, int k);              // the polygon left valid
bool delete_solid(Model& m, SolidRef s);                     // a region left without solids goes too
bool rename_region(Model& m, int region, std::string name);  // unique, a deck identifier
bool set_material(Model& m, int region, std::string material);
bool raise_region(Model& m, int region);  // one layer up (later in the deck)
bool lower_region(Model& m, int region);
int duplicate_solid(Model& m, SolidRef s, P3 delta);  // a new region on top, its index
bool mirror_solid(Model& m, SolidRef s, int axis);    // about its bounds' centre across `axis`

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
