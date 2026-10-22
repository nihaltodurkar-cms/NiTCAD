#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/render/vectors.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <map>
#include <numbers>
#include <string>
#include <utility>

namespace NiTCAD::render {

namespace {

constexpr double um_per_cm = 1e4;
constexpr float nan_f = std::numeric_limits<float>::quiet_NaN();

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

Float3 um(const Point& p) noexcept {
    return {static_cast<float>(p[0] * um_per_cm), static_cast<float>(p[1] * um_per_cm),
            static_cast<float>(p[2] * um_per_cm)};
}

Point sub(const Point& a, const Point& b) noexcept { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Point add(const Point& a, const Point& b) noexcept { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Point scale(const Point& a, double s) noexcept { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Point& a, const Point& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Point cross(const Point& a, const Point& b) noexcept {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Point unit(const Point& a) noexcept {
    const double n = std::sqrt(dot(a, a));
    return n > 0.0 ? scale(a, 1.0 / n) : Point{0.0, 0.0, 1.0};
}
Float3 to_f(const Point& p) noexcept {
    return {static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2])};
}

// The (unnormalized) normal of a polygon (Newell).
Point polygon_normal(std::span<const Point> pts) noexcept {
    Point n{};
    for (std::size_t i = 0; i < pts.size(); ++i) n = add(n, cross(pts[i], pts[(i + 1) % pts.size()]));
    return n;
}

// Appends a polygon as a triangle fan facing `normal` (reversed when it faces away).
void fan(TriangleBatch& batch, std::vector<Point> pts, const std::vector<float>& values,
         const Rgba& color, const Point& normal) {
    if (pts.size() < 3) return;
    std::vector<float> vals = values;
    // Clipping can leave vertices a rounding apart: merge those within 1e-9 of the polygon's size.
    {
        double size = 0.0;
        for (const Point& p : pts) size = std::max(size, std::sqrt(dot(sub(p, pts[0]), sub(p, pts[0]))));
        const double tol = 1e-9 * size;
        std::vector<Point> kept;
        std::vector<float> kept_vals;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            const Point& next = pts[(i + 1) % pts.size()];
            if (std::sqrt(dot(sub(next, pts[i]), sub(next, pts[i]))) <= tol && i + 1 < pts.size()) continue;
            if (!kept.empty() && i + 1 == pts.size() &&
                std::sqrt(dot(sub(kept.front(), pts[i]), sub(kept.front(), pts[i]))) <= tol) {
                continue;
            }
            kept.push_back(pts[i]);
            if (!vals.empty()) kept_vals.push_back(vals[i]);
        }
        if (kept.size() < 3) return;
        pts = std::move(kept);
        if (!vals.empty()) vals = std::move(kept_vals);
    }
    if (dot(polygon_normal(pts), normal) < 0.0) {
        std::ranges::reverse(pts);
        std::ranges::reverse(vals);
    }
    const Float3 n = to_f(unit(normal));
    for (std::size_t i = 1; i + 1 < pts.size(); ++i) {
        for (const std::size_t k : {std::size_t{0}, i, i + 1}) {
            batch.vertices.push_back({um(pts[k]), n, vals.empty() ? nan_f : vals[k], color});
        }
    }
}

struct Bounds {
    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
    Point hi{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity()};
    void add(const Float3& p) noexcept {
        for (std::size_t a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], static_cast<double>(p[a]));
            hi[a] = std::max(hi[a], static_cast<double>(p[a]));
        }
    }
};

// Bounds and the default camera of a scene from everything in it.
void frame(FieldScene& s, const ViewOptions& o) {
    Bounds b;
    for (const TriangleBatch& t : s.triangles) for (const SceneVertex& v : t.vertices) b.add(v.position);
    for (const LineBatch& l : s.lines) for (const Float3& p : l.points) b.add(p);
    if (!(b.lo[0] <= b.hi[0])) b = Bounds{{0, 0, 0}, {1, 1, 1}};
    s.lo = to_f(b.lo);
    s.hi = to_f(b.hi);
    const Point c = scale(add(b.lo, b.hi), 0.5);
    const Point d = sub(b.hi, b.lo);
    s.camera.target = to_f(c);
    s.camera.projection = o.projection;
    if (s.dimension == 2) {
        s.camera.azimuth_deg = 0.0f;
        s.camera.elevation_deg = 0.0f;
        s.camera.extent_um = static_cast<float>(0.55 * std::max(d[0], d[1]));
        s.axis_triad = false;
    } else {
        s.camera.azimuth_deg = o.azimuth_deg;
        s.camera.elevation_deg = o.elevation_deg;
        s.camera.extent_um = static_cast<float>(0.5 * std::sqrt(dot(d, d)));
    }
}

bool same_mesh(const device::Device& d, const Geometry& g) {
    if (&d.mesh() == &g.mesh()) return true;
    const auto a = d.mesh().points(), b = g.mesh().points();
    return a.size() == b.size() && std::ranges::equal(a, b) &&
           d.mesh().dimension() == g.mesh().dimension();
}

std::expected<void, base::Error> check_view(const device::Device& d, const Geometry& g) {
    if (!same_mesh(d, g)) return std::unexpected(invalid("the geometry's mesh is not the device's"));
    if (g.dimension() < 2 || g.element_count() == 0) {
        return std::unexpected(invalid("a 2D or 3D geometry with elements is needed for this view"));
    }
    return {};
}

std::vector<std::int32_t> regions_of(const device::Device& d) {
    return {d.node_region().begin(), d.node_region().end()};
}

// The batch of each material, in the order of the Material enum.
struct Batches {
    std::array<TriangleBatch, 7> by{};
    Batches(const Palette& p, bool lit, bool transparent_insulators) {
        for (std::size_t m = 0; m < by.size(); ++m) {
            by[m].name = std::string(material_name(static_cast<Material>(m)));
            by[m].lit = lit;
            by[m].transparent = transparent_insulators && static_cast<Material>(m) == Material::insulator;
        }
        (void)p;
    }
    TriangleBatch& operator[](Material m) { return by[static_cast<std::size_t>(m)]; }
    void move_into(std::vector<TriangleBatch>& out) {
        for (TriangleBatch& b : by) {
            if (!b.vertices.empty()) out.push_back(std::move(b));
        }
    }
};

// Boundary edges of a 2D geometry: element edges that appear once, oriented as in their element
// (counterclockwise, so the outward normal is (dy, -dx)).
std::vector<std::array<NodeId, 2>> boundary_edges_2d(const Geometry& g) {
    std::map<std::pair<NodeId, NodeId>, std::pair<int, std::array<NodeId, 2>>> count;
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        const auto v = g.element(e);
        for (std::size_t i = 0; i < 3; ++i) {
            const NodeId a = v[i], b = v[(i + 1) % 3];
            auto& c = count[{std::min(a, b), std::max(a, b)}];
            ++c.first;
            c.second = {a, b};
        }
    }
    std::vector<std::array<NodeId, 2>> out;
    for (const auto& [key, c] : count) {
        if (c.first == 1) out.push_back(c.second);
    }
    return out;
}

// Every element edge once (2D mesh view).
std::vector<std::array<NodeId, 2>> element_edges(const Geometry& g) {
    std::vector<std::array<NodeId, 2>> out;
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        const auto v = g.element(e);
        for (std::size_t i = 0; i < v.size(); ++i) {
            for (std::size_t j = i + 1; j < v.size(); ++j) {
                out.push_back({std::min(v[i], v[j]), std::max(v[i], v[j])});
            }
        }
    }
    std::ranges::sort(out);
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

double largest_extent_cm(const mesh::Mesh& m) {
    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
    Point hi = scale(lo, -1.0);
    for (const Point& p : m.points()) {
        for (std::size_t a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], p[a]);
            hi[a] = std::max(hi[a], p[a]);
        }
    }
    return std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
}

// Contacts as slabs of nominal thickness on their boundary (VR8): a lumped gate as its oxide (pale)
// under the electrode, an electrode as metal, an ohmic contact as a thin metal pad.
void contact_slabs(const device::Device& d, const Geometry& g, const Palette& palette,
                   std::vector<TriangleBatch>& out, const std::optional<Box>& clip = std::nullopt,
                   std::span<const Point> contact_offset = {}) {
    const double t = 0.03 * largest_extent_cm(d.mesh());
    const auto points = d.mesh().points();
    TriangleBatch metal{"electrode", Shading::flat, g.dimension() == 3, false, {}};
    TriangleBatch oxide{"gate oxide", Shading::flat, g.dimension() == 3, true, {}};
    const Rgba metal_color = material_color(Material::metal, palette);
    Rgba pad_color = metal_color;
    pad_color.r += 0.12f;
    pad_color.g += 0.12f;
    pad_color.b += 0.12f;
    const Rgba oxide_color = material_color(Material::insulator, palette);
    // Prism over a polygon (outward normal n) from offset a to offset b.
    const auto prism = [&](TriangleBatch& batch, const std::vector<Point>& base, const Point& n,
                           double a, double b, const Rgba& c,
                           const std::vector<std::array<std::size_t, 2>>& sides) {
        std::vector<Point> top, bottom;
        for (const Point& p : base) {
            bottom.push_back(add(p, scale(n, a)));
            top.push_back(add(p, scale(n, b)));
        }
        if (base.size() >= 3) {
            fan(batch, top, {}, c, n);
            if (a > 0.0) fan(batch, bottom, {}, c, scale(n, -1.0));
        }
        for (const auto& [i, j] : sides) {
            const Point side_n = cross(sub(base[j], base[i]), n);
            if (g.dimension() == 2) {
                // 2D: the slab is the quad bottom_i, bottom_j, top_j, top_i.
                fan(batch, {bottom[i], bottom[j], top[j], top[i]}, {}, c, {0.0, 0.0, 1.0});
            } else {
                fan(batch, {bottom[i], bottom[j], top[j], top[i]}, {}, c, side_n);
            }
        }
    };
    for (std::size_t ci = 0; ci < d.contacts().size(); ++ci) {
        const device::Contact& contact = d.contacts()[ci];
        // The exploded view moves each contact with its region.
        const Point off = contact_offset.empty() ? Point{} : contact_offset[ci];
        std::vector<bool> in(d.mesh().node_count(), false);
        for (const NodeId n : contact.nodes) in[static_cast<std::size_t>(n)] = true;
        const bool gate = contact.kind == device::ContactKind::gate;
        const bool pad = contact.kind == device::ContactKind::ohmic;
        const double lo = gate ? 0.5 * t : 0.0;
        const double hi = pad ? t / 3.0 : lo + t;
        TriangleBatch& target = metal;
        const Rgba color = pad ? pad_color : metal_color;
        if (g.dimension() == 2) {
            for (const auto& [a, b] : boundary_edges_2d(g)) {
                if (!in[static_cast<std::size_t>(a)] || !in[static_cast<std::size_t>(b)]) continue;
                const Point pa = add(points[static_cast<std::size_t>(a)], off), pb = add(points[static_cast<std::size_t>(b)], off);
                const Point n = unit({pb[1] - pa[1], -(pb[0] - pa[0]), 0.0});
                const std::vector<Point> base{pa, pb};
                if (gate) prism(oxide, base, n, 0.0, lo, oxide_color, {{0, 1}});
                prism(target, base, n, lo, hi, color, {{0, 1}});
            }
            continue;
        }
        // 3D: the boundary faces wholly on the contact; sides along the patch's outline.
        std::vector<Polygon> faces;
        for (Polygon& f : g.boundary_faces()) {
            bool all = true;
            for (const SurfaceVertex& v : f.vertices) all = all && in[static_cast<std::size_t>(v.at.nodes[0])];
            if (all && clip) {
                // Under a clipping box: a face whose centroid lies in the box is left out.
                Point c{};
                for (const SurfaceVertex& v : f.vertices) c = add(c, scale(v.position, 1.0 / 3.0));
                bool in_box = true;
                for (std::size_t a = 0; a < 3; ++a) in_box = in_box && c[a] >= clip->lo[a] && c[a] <= clip->hi[a];
                all = !in_box;
            }
            if (all) faces.push_back(std::move(f));
        }
        std::map<std::pair<NodeId, NodeId>, int> edge_count;
        for (const Polygon& f : faces) {
            for (std::size_t i = 0; i < 3; ++i) {
                const NodeId a = f.vertices[i].at.nodes[0], b = f.vertices[(i + 1) % 3].at.nodes[0];
                ++edge_count[{std::min(a, b), std::max(a, b)}];
            }
        }
        for (const Polygon& f : faces) {
            std::vector<Point> base;
            for (const SurfaceVertex& v : f.vertices) base.push_back(add(v.position, off));
            const Point n = unit(polygon_normal(base));
            std::vector<std::array<std::size_t, 2>> sides;
            for (std::size_t i = 0; i < 3; ++i) {
                const NodeId a = f.vertices[i].at.nodes[0], b = f.vertices[(i + 1) % 3].at.nodes[0];
                if (edge_count[{std::min(a, b), std::max(a, b)}] == 1) sides.push_back({i, (i + 1) % 3});
            }
            if (gate) prism(oxide, base, n, 0.0, lo, oxide_color, sides);
            prism(target, base, n, lo, hi, color, sides);
        }
    }
    if (!oxide.vertices.empty()) out.push_back(std::move(oxide));
    if (!metal.vertices.empty()) out.push_back(std::move(metal));
}

// The device's boundary, faint and translucent, for context in 3D.
TriangleBatch outline_shell(const Geometry& g) {
    TriangleBatch shell{"outline", Shading::flat, true, true, {}};
    for (const Polygon& f : g.boundary_faces()) {
        std::vector<Point> pts;
        for (const SurfaceVertex& v : f.vertices) pts.push_back(v.position);
        fan(shell, pts, {}, {0.82f, 0.85f, 0.88f, 0.12f}, polygon_normal(pts));
    }
    return shell;
}

// 2D outlines: the boundary and the interfaces between regions.
void outlines_2d(const device::Device& d, const Geometry& g, std::vector<LineBatch>& out) {
    LineBatch lines{"outline", {0.12f, 0.14f, 0.16f, 1.0f}, 1.5f, true, {}};
    const auto points = d.mesh().points();
    for (const auto& [a, b] : boundary_edges_2d(g)) {
        lines.points.push_back(um(points[static_cast<std::size_t>(a)]));
        lines.points.push_back(um(points[static_cast<std::size_t>(b)]));
    }
    for (const InterfacePiece& p : g.interfaces(regions_of(d))) {
        lines.points.push_back(um(p.points[0]));
        lines.points.push_back(um(p.points[1]));
    }
    out.push_back(std::move(lines));
}

// The value of a field (positions) at a point of a face polygon, by barycentric weights in the
// face's triangle fan.
float position_at(const Geometry& g, std::size_t element, const Point& x, std::span<const double> values,
                  const ColorScale& scale_) {
    return static_cast<float>(scale_.position(interpolate(g.interpolant(element, x), values)));
}

// Contour lines at the colorbar's major ticks, or at the given colormap positions.
void contours(const TriangleBatch& field, const Colorbar& bar, std::vector<LineBatch>& out,
              std::span<const double> positions = {}) {
    LineBatch lines{"contours", {0.1f, 0.1f, 0.1f, 0.85f}, 1.0f, true, {}};
    std::vector<double> levels(positions.begin(), positions.end());
    if (positions.empty()) {
        for (const Tick& t : bar.scale.ticks()) {
            if (t.major && t.position > 0.0 && t.position < 1.0) levels.push_back(t.position);
        }
    }
    for (std::size_t i = 0; i + 2 < field.vertices.size(); i += 3) {
        std::array<Point, 3> p{};
        std::array<double, 3> v{};
        for (std::size_t k = 0; k < 3; ++k) {
            const SceneVertex& s = field.vertices[i + k];
            p[k] = {s.position[0], s.position[1], s.position[2]};
            v[k] = s.value;
        }
        for (const double level : levels) {
            if (const auto seg = contour_segment(p, v, level)) {
                lines.points.push_back(to_f((*seg)[0]));
                lines.points.push_back(to_f((*seg)[1]));
            }
        }
    }
    if (!lines.points.empty()) out.push_back(std::move(lines));
}

// The solid layer of a 3D view: boundary pieces by material, and the semiconductor side of each
// interface with an insulator (seen through it). With `field`, semiconductor pieces carry the
// field's positions instead of their material color.
void solid_3d(const device::Device& d, const Geometry& g, const ViewOptions& o, FieldScene& s,
              std::span<const double> values, const ColorScale* field_scale) {
    Batches batches(o.palette, true, o.transparent_insulators);
    TriangleBatch field_batch{"field", Shading::field, true, false, {}};
    const auto is_insulator = [&](NodeId n) { return d.is_insulator(n); };
    // With a clipping box (Unit 25c), the surface of what is left: cut faces drawn like the boundary.
    for (const Polygon& face : o.clip ? g.clipped_surface(*o.clip) : g.boundary_faces()) {
        std::vector<Point> pts;
        for (const SurfaceVertex& v : face.vertices) pts.push_back(v.position);
        const Point n = polygon_normal(pts);
        for (const OwnedPolygon& piece : g.owner_pieces(face)) {
            const Material m = material_class(d, piece.owner, o.palette);
            if (field_scale && (!is_insulator(piece.owner) || !o.transparent_insulators)) {
                std::vector<float> vals;
                for (const Point& x : piece.points) vals.push_back(position_at(g, face.element, x, values, *field_scale));
                fan(field_batch, piece.points, vals, {}, n);
            } else {
                fan(batches[m], piece.points, {}, material_color(m, o.palette), n);
            }
        }
    }
    if (o.transparent_insulators) {
        for (const InterfacePiece& p : g.interfaces(regions_of(d))) {
            const bool ia = is_insulator(p.node_a), ib = is_insulator(p.node_b);
            if (ia == ib) continue;
            const NodeId semi = ia ? p.node_b : p.node_a;
            const NodeId ins = ia ? p.node_a : p.node_b;
            const Point toward = sub(d.mesh().points()[static_cast<std::size_t>(ins)],
                                     d.mesh().points()[static_cast<std::size_t>(semi)]);
            Point n = polygon_normal(p.points);
            if (dot(n, toward) < 0.0) n = scale(n, -1.0);
            std::vector<std::vector<Point>> parts;
            if (o.clip) parts = outside_box(p.points, *o.clip);
            else parts.push_back(p.points);
            for (const std::vector<Point>& part : parts) {
                if (field_scale) {
                    std::vector<float> vals;
                    for (const Point& x : part) vals.push_back(position_at(g, p.element, x, values, *field_scale));
                    fan(field_batch, part, vals, {}, n);
                } else {
                    const Material m = material_class(d, semi, o.palette);
                    fan(batches[m], part, {}, material_color(m, o.palette), n);
                }
            }
        }
    }
    if (!field_batch.vertices.empty()) s.triangles.push_back(std::move(field_batch));
    batches.move_into(s.triangles);
}

}  // namespace

// --- fields -------------------------------------------------------------------------------------

FieldInfo field_info(Field f) noexcept {
    switch (f) {
        case Field::potential: return {"potential", "\xCF\x88 (V)", ScaleKind::linear};
        case Field::electrons: return {"n", "n (cm\xE2\x81\xBB\xC2\xB3)", ScaleKind::log10};
        case Field::holes: return {"p", "p (cm\xE2\x81\xBB\xC2\xB3)", ScaleKind::log10};
        case Field::conduction_band: return {"Ec", "Ec (eV)", ScaleKind::linear};
        case Field::valence_band: return {"Ev", "Ev (eV)", ScaleKind::linear};
        case Field::electron_fermi: return {"Efn", "EFn (eV)", ScaleKind::linear};
        case Field::hole_fermi: return {"Efp", "EFp (eV)", ScaleKind::linear};
        case Field::temperature: return {"T", "T (K)", ScaleKind::linear};
        case Field::net_doping: return {"doping", "ND \xE2\x88\x92 NA (cm\xE2\x81\xBB\xC2\xB3)", ScaleKind::signed_log10};
        case Field::donors: return {"ND", "ND (cm\xE2\x81\xBB\xC2\xB3)", ScaleKind::log10};
        case Field::acceptors: return {"NA", "NA (cm\xE2\x81\xBB\xC2\xB3)", ScaleKind::log10};
        // Rates in cm^-3 s^-1: signed ones on the signed log scale, generations on log10.
        case Field::recombination: return {"R", "R \xE2\x88\x92 G (cm\xE2\x81\xBB\xC2\xB3 s\xE2\x81\xBB\xC2\xB9)", ScaleKind::signed_log10};
        case Field::srh: return {"srh", "SRH (cm\xE2\x81\xBB\xC2\xB3 s\xE2\x81\xBB\xC2\xB9)", ScaleKind::signed_log10};
        case Field::auger: return {"auger", "Auger (cm\xE2\x81\xBB\xC2\xB3 s\xE2\x81\xBB\xC2\xB9)", ScaleKind::signed_log10};
        case Field::radiative: return {"radiative", "radiative (cm\xE2\x81\xBB\xC2\xB3 s\xE2\x81\xBB\xC2\xB9)", ScaleKind::signed_log10};
        case Field::surface: return {"surface", "surface, per volume (cm\xE2\x81\xBB\xC2\xB3 s\xE2\x81\xBB\xC2\xB9)", ScaleKind::signed_log10};
        case Field::impact_generation: return {"Gii", "G impact (cm\xE2\x81\xBB\xC2\xB3 s\xE2\x81\xBB\xC2\xB9)", ScaleKind::log10};
        case Field::btbt_generation: return {"Gbtbt", "G band-to-band (cm\xE2\x81\xBB\xC2\xB3 s\xE2\x81\xBB\xC2\xB9)", ScaleKind::log10};
        // Vector magnitudes, log10 (Unit 25b).
        case Field::electric_field: return {"E", "|E| (V/cm)", ScaleKind::log10};
        case Field::current_density: return {"J", "|J| (A/cm\xC2\xB2)", ScaleKind::log10};
        case Field::electron_current: return {"Jn", "|Jn| (A/cm\xC2\xB2)", ScaleKind::log10};
        case Field::hole_current: return {"Jp", "|Jp| (A/cm\xC2\xB2)", ScaleKind::log10};
    }
    return {"?", "?", ScaleKind::linear};
}

std::optional<Field> find_field(std::string_view key) noexcept {
    for (int f = 0; f <= static_cast<int>(Field::last); ++f) {
        if (field_info(static_cast<Field>(f)).key == key) return static_cast<Field>(f);
    }
    return std::nullopt;
}

std::expected<std::vector<double>, base::Error> field_values(const device::Device& device,
                                                             const results::BiasPoint& point,
                                                             Field f) {
    const std::size_t n = device.mesh().node_count();
    const auto take = [&](const std::vector<double>& v) -> std::expected<std::vector<double>, base::Error> {
        if (v.size() != n) return std::unexpected(invalid("the result's field is not one value per mesh node", v.size()));
        return v;
    };
    switch (f) {
        case Field::potential: return take(point.fields.potential_V);
        case Field::electrons: return take(point.fields.n_cm3);
        case Field::holes: return take(point.fields.p_cm3);
        case Field::conduction_band: return take(point.bands.conduction_eV);
        case Field::valence_band: return take(point.bands.valence_eV);
        case Field::electron_fermi: return take(point.bands.electron_fermi_eV);
        case Field::hole_fermi: return take(point.bands.hole_fermi_eV);
        case Field::temperature:
            if (point.fields.temperature_K.empty()) return std::vector<double>(n, device.temperature_K());
            return take(point.fields.temperature_K);
        case Field::net_doping:
        case Field::donors:
        case Field::acceptors: {
            std::vector<double> v(n);
            for (std::size_t i = 0; i < n; ++i) {
                v[i] = f == Field::donors      ? device.donors()[i]
                       : f == Field::acceptors ? device.acceptors()[i]
                                               : device.donors()[i] - device.acceptors()[i];
            }
            return v;
        }
        case Field::recombination:
        case Field::srh:
        case Field::auger:
        case Field::radiative:
        case Field::surface:
        case Field::impact_generation:
        case Field::btbt_generation: {
            const results::RateFields& r = point.rates;
            const std::vector<double>& v = f == Field::recombination ? r.net_n
                                           : f == Field::srh         ? r.srh
                                           : f == Field::auger       ? r.auger
                                           : f == Field::radiative   ? r.radiative
                                           : f == Field::surface     ? r.interface_n
                                           : f == Field::impact_generation ? r.impact
                                                                           : r.band_to_band_n;
            if (v.size() != n) {
                return std::unexpected(invalid(
                    "the result holds no " + std::string(field_info(f).key) +
                        " rates (stored with [solver] rates = on, for a mechanism that is on)",
                    v.size()));
            }
            if (f != Field::impact_generation && f != Field::btbt_generation) return v;
            std::vector<double> g(n);
            for (std::size_t i = 0; i < n; ++i) g[i] = -v[i];  // a generation, positive
            return g;
        }
        case Field::electric_field:
        case Field::current_density:
        case Field::electron_current:
        case Field::hole_current:
            return std::unexpected(invalid(std::string(field_info(f).key) +
                                           " is reconstructed on the mesh's elements: it needs the geometry"));
    }
    return std::unexpected(invalid("unknown field"));
}

bool is_vector_magnitude(Field f) noexcept {
    return f == Field::electric_field || f == Field::current_density || f == Field::electron_current ||
           f == Field::hole_current;
}

std::expected<std::vector<double>, base::Error> field_values(const device::Device& device,
                                                             const Geometry& geometry,
                                                             const results::BiasPoint& point,
                                                             Field f) {
    if (!is_vector_magnitude(f)) return field_values(device, point, f);
    const VectorField v = f == Field::electric_field     ? VectorField::electric
                          : f == Field::current_density  ? VectorField::current
                          : f == Field::electron_current ? VectorField::electron_current
                                                         : VectorField::hole_current;
    auto e = element_vectors(device, geometry, point, v);
    if (!e) return std::unexpected(std::move(e.error()));
    return magnitudes(nodal_vectors(geometry, *e));
}

// --- region colors ------------------------------------------------------------------------------

Material material_class(const device::Device& device, NodeId node, const Palette& p) noexcept {
    if (device.is_insulator(node)) return Material::insulator;
    const double net = device.net_doping(node);
    if (net >= p.heavy_cm3) return Material::n_heavy;
    if (net > 0.0) return Material::n_light;
    if (-net >= p.heavy_cm3) return Material::p_heavy;
    if (-net >= p.channel_cm3) return Material::p_channel;
    return Material::substrate;
}

Rgba material_color(Material m, const Palette& p) noexcept {
    switch (m) {
        case Material::insulator: return {0.80f, 0.87f, 0.94f, p.insulator_alpha};
        case Material::n_heavy: return {0.84f, 0.38f, 0.29f, 1.0f};
        case Material::n_light: return {0.93f, 0.66f, 0.58f, 1.0f};
        case Material::p_channel: return {0.42f, 0.58f, 0.82f, 1.0f};
        case Material::p_heavy: return {0.22f, 0.36f, 0.66f, 1.0f};
        case Material::substrate: return {0.89f, 0.83f, 0.69f, 1.0f};
        case Material::metal: return {0.33f, 0.35f, 0.38f, 1.0f};
    }
    return {0.5f, 0.5f, 0.5f, 1.0f};
}

std::string_view material_name(Material m) noexcept {
    switch (m) {
        case Material::insulator: return "insulator";
        case Material::n_heavy: return "n+";
        case Material::n_light: return "n";
        case Material::p_channel: return "p";
        case Material::p_heavy: return "p+";
        case Material::substrate: return "substrate";
        case Material::metal: return "metal";
    }
    return "?";
}

// --- views ---------------------------------------------------------------------------------------

std::expected<FieldScene, base::Error> solid_view(const device::Device& d, const Geometry& g,
                                                  const ViewOptions& o) {
    if (auto c = check_view(d, g); !c) return std::unexpected(std::move(c.error()));
    FieldScene s;
    s.dimension = g.dimension();
    s.title = "Regions";
    if (g.dimension() == 2) {
        Batches batches(o.palette, false, false);
        for (std::size_t e = 0; e < g.element_count(); ++e) {
            for (const OwnedPolygon& piece : g.owner_pieces(g.element_polygon(e))) {
                const Material m = material_class(d, piece.owner, o.palette);
                Rgba c = material_color(m, o.palette);
                c.a = 1.0f;  // 2D maps are opaque
                fan(batches[m], piece.points, {}, c, {0.0, 0.0, 1.0});
            }
        }
        batches.move_into(s.triangles);
        if (o.region_outlines) outlines_2d(d, g, s.lines);
    } else {
        solid_3d(d, g, o, s, {}, nullptr);
    }
    if (o.contacts) contact_slabs(d, g, o.palette, s.triangles, o.clip);
    frame(s, o);
    return s;
}

std::expected<FieldScene, base::Error> mesh_view(const device::Device& d, const Geometry& g,
                                                 const ViewOptions& o) {
    auto s = solid_view(d, g, o);
    if (!s) return s;
    s->title = "Mesh";
    const auto points = d.mesh().points();
    LineBatch edges{"mesh", {0.15f, 0.17f, 0.2f, 0.75f}, 1.0f, g.dimension() == 2, {}};
    if (g.dimension() == 2) {
        for (const auto& [a, b] : element_edges(g)) {
            edges.points.push_back(um(points[static_cast<std::size_t>(a)]));
            edges.points.push_back(um(points[static_cast<std::size_t>(b)]));
        }
    } else {
        std::vector<std::array<NodeId, 2>> list;
        for (const Polygon& f : g.boundary_faces()) {
            for (std::size_t i = 0; i < 3; ++i) {
                const NodeId a = f.vertices[i].at.nodes[0], b = f.vertices[(i + 1) % 3].at.nodes[0];
                list.push_back({std::min(a, b), std::max(a, b)});
            }
        }
        std::ranges::sort(list);
        list.erase(std::unique(list.begin(), list.end()), list.end());
        for (const auto& [a, b] : list) {
            edges.points.push_back(um(points[static_cast<std::size_t>(a)]));
            edges.points.push_back(um(points[static_cast<std::size_t>(b)]));
        }
    }
    s->lines.push_back(std::move(edges));
    return s;
}

namespace {

std::expected<Colorbar, base::Error> colorbar_for(std::span<const double> values, const FieldOptions& o) {
    if (o.scale) return Colorbar{o.map, *o.scale, o.title};
    auto fitted = ColorScale::fit(ScaleKind::linear, values);
    if (!fitted) return std::unexpected(std::move(fitted.error()));
    return Colorbar{o.map, *fitted, o.title};
}

}  // namespace

std::expected<FieldScene, base::Error> field_view(const device::Device& d, const Geometry& g,
                                                  std::span<const double> values,
                                                  const FieldOptions& o) {
    if (auto c = check_view(d, g); !c) return std::unexpected(std::move(c.error()));
    if (values.size() != d.mesh().node_count()) {
        return std::unexpected(invalid("the field is not one value per mesh node", values.size()));
    }
    auto bar = colorbar_for(values, o);
    if (!bar) return std::unexpected(std::move(bar.error()));
    FieldScene s;
    s.dimension = g.dimension();
    s.title = o.title;
    if (g.dimension() == 2) {
        TriangleBatch field{"field", Shading::field, false, false, {}};
        for (std::size_t e = 0; e < g.element_count(); ++e) {
            const Polygon p = g.element_polygon(e);
            std::vector<Point> pts;
            std::vector<float> vals;
            for (const SurfaceVertex& v : p.vertices) {
                pts.push_back(v.position);
                vals.push_back(static_cast<float>(bar->scale.position(interpolate(v.at, values))));
            }
            fan(field, pts, vals, {}, {0.0, 0.0, 1.0});
        }
        if (o.contours) contours(field, *bar, s.lines);
        s.triangles.push_back(std::move(field));
        if (o.view.region_outlines) outlines_2d(d, g, s.lines);
    } else {
        solid_3d(d, g, o.view, s, values, &bar->scale);
        if (o.contours) {
            for (const TriangleBatch& t : s.triangles) {
                if (t.shading == Shading::field) contours(t, *bar, s.lines);
            }
        }
    }
    if (o.view.contacts) contact_slabs(d, g, o.view.palette, s.triangles, o.view.clip);
    if (o.view.mesh_edges && g.dimension() == 2) {
        LineBatch edges{"mesh", {0.15f, 0.17f, 0.2f, 0.5f}, 1.0f, true, {}};
        const auto points = d.mesh().points();
        for (const auto& [a, b] : element_edges(g)) {
            edges.points.push_back(um(points[static_cast<std::size_t>(a)]));
            edges.points.push_back(um(points[static_cast<std::size_t>(b)]));
        }
        s.lines.push_back(std::move(edges));
    }
    s.colorbar = std::move(*bar);
    frame(s, o.view);
    return s;
}

std::expected<FieldScene, base::Error> slice_view(const device::Device& d, const Geometry& g,
                                                  std::span<const double> values,
                                                  std::span<const Plane> planes,
                                                  const FieldOptions& o) {
    if (auto c = check_view(d, g); !c) return std::unexpected(std::move(c.error()));
    if (g.dimension() != 3) return std::unexpected(invalid("slices need a 3D geometry"));
    if (values.size() != d.mesh().node_count()) {
        return std::unexpected(invalid("the field is not one value per mesh node", values.size()));
    }
    auto bar = colorbar_for(values, o);
    if (!bar) return std::unexpected(std::move(bar.error()));
    FieldScene s;
    s.dimension = 3;
    s.title = o.title;
    TriangleBatch field{"slice", Shading::field, true, false, {}};
    for (std::size_t k = 0; k < planes.size(); ++k) {
        const std::vector<Polygon> polys = g.slice(planes[k]);
        if (polys.empty()) return std::unexpected(invalid("a slice plane misses the device", k));
        for (const Polygon& p : polys) {
            std::vector<Point> pts;
            std::vector<float> vals;
            for (const SurfaceVertex& v : p.vertices) {
                pts.push_back(v.position);
                vals.push_back(static_cast<float>(bar->scale.position(interpolate(v.at, values))));
            }
            fan(field, pts, vals, {}, planes[k].normal);
        }
    }
    if (o.contours) contours(field, *bar, s.lines);
    s.triangles.push_back(std::move(field));
    s.triangles.push_back(outline_shell(g));
    if (o.view.contacts) contact_slabs(d, g, o.view.palette, s.triangles, o.view.clip);
    s.colorbar = std::move(*bar);
    frame(s, o.view);
    return s;
}

std::expected<FieldScene, base::Error> iso_view(const device::Device& d, const Geometry& g,
                                                std::span<const double> values,
                                                std::span<const double> levels, const FieldOptions& o) {
    if (auto c = check_view(d, g); !c) return std::unexpected(std::move(c.error()));
    if (values.size() != d.mesh().node_count()) {
        return std::unexpected(invalid("the field is not one value per mesh node", values.size()));
    }
    if (levels.empty()) return std::unexpected(invalid("an isosurface needs a level"));
    auto bar = colorbar_for(values, o);
    if (!bar) return std::unexpected(std::move(bar.error()));
    for (std::size_t k = 0; k < levels.size(); ++k) {
        if (!std::isfinite(bar->scale.transform(levels[k]))) {
            return std::unexpected(invalid("a level is not finite on the field's scale", k));
        }
    }
    if (g.dimension() == 2) {
        FieldOptions plain = o;
        plain.contours = false;
        auto s = field_view(d, g, values, plain);
        if (!s) return s;
        std::vector<double> positions;
        for (const double level : levels) positions.push_back(bar->scale.position(level));
        for (const TriangleBatch& t : s->triangles) {
            if (t.shading == Shading::field) contours(t, *bar, s->lines, positions);
        }
        return s;
    }
    FieldScene s;
    s.dimension = 3;
    s.title = o.title;
    // On the scale's coordinate (BQ3), so a level lies where the colormap shows it.
    std::vector<double> scaled(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) scaled[i] = bar->scale.transform(values[i]);
    const bool one = levels.size() == 1;
    for (const double level : levels) {
        const IsoSurface iso = isosurface(g, scaled, bar->scale.transform(level));
        // Smooth shading: each vertex's normal is the area-weighted sum of its triangles'.
        std::vector<Point> normal(iso.vertices.size(), Point{});
        for (const auto& t : iso.triangles) {
            const Point& p0 = iso.vertices[t[0]].position;
            const Point n = cross(sub(iso.vertices[t[1]].position, p0), sub(iso.vertices[t[2]].position, p0));
            for (const std::uint32_t v : t) normal[v] = add(normal[v], n);
        }
        TriangleBatch batch{"isosurface", Shading::field, true, !one, {}};
        const float position = static_cast<float>(bar->scale.position(level));
        const Rgba tint{0.0f, 0.0f, 0.0f, one ? 1.0f : 0.55f};
        for (const auto& t : iso.triangles) {
            for (const std::uint32_t v : t) {
                batch.vertices.push_back({um(iso.vertices[v].position), to_f(unit(normal[v])), position, tint});
            }
        }
        if (!batch.vertices.empty()) s.triangles.push_back(std::move(batch));
    }
    s.triangles.push_back(outline_shell(g));
    if (o.view.contacts) contact_slabs(d, g, o.view.palette, s.triangles, o.view.clip);
    s.colorbar = std::move(*bar);
    frame(s, o.view);
    return s;
}

void add_streamlines(FieldScene& s, std::span<const Streamline> lines, const Rgba& color) {
    LineBatch batch{"streamlines", color, 1.5f, s.dimension == 2, {}};
    for (const Streamline& l : lines) {
        for (std::size_t i = 1; i < l.points.size(); ++i) {
            batch.points.push_back(um(l.points[i - 1]));
            batch.points.push_back(um(l.points[i]));
        }
    }
    if (!batch.points.empty()) s.lines.push_back(std::move(batch));
}

void add_cut_line(FieldScene& s, const CutLine& cut) {
    const Point d{cut.b[0] - cut.a[0], cut.b[1] - cut.a[1], cut.b[2] - cut.a[2]};
    const double length = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    LineBatch edge{"cut line", {0.0f, 0.0f, 0.0f, 1.0f}, 4.0f, true, {}};
    for (std::size_t i = 1; i < cut.samples.size(); ++i) {
        if (cut.samples[i].gap_before) continue;
        for (const std::size_t k : {i - 1, i}) {
            const double t = cut.samples[k].distance / length;
            edge.points.push_back(um({cut.a[0] + t * d[0], cut.a[1] + t * d[1], cut.a[2] + t * d[2]}));
        }
    }
    if (edge.points.empty()) return;
    LineBatch line = edge;
    line.color = {1.0f, 1.0f, 1.0f, 1.0f};
    line.width_px = 2.0f;
    s.lines.push_back(std::move(edge));
    s.lines.push_back(std::move(line));
}

void add_streamlines(FieldScene& s, std::span<const Streamline> lines, const Colorbar& colors) {
    constexpr std::size_t bins = 64;
    std::vector<LineBatch> batch(bins + 1);  // the last for segments without a position
    for (std::size_t b = 0; b <= bins; ++b) {
        const Rgb c = b < bins ? sample(colors.map, (static_cast<double>(b) + 0.5) / bins) : no_data;
        batch[b] = {"streamlines", {c.r, c.g, c.b, 1.0f}, 2.0f, s.dimension == 2, {}};
    }
    for (const Streamline& l : lines) {
        for (std::size_t i = 1; i < l.points.size(); ++i) {
            const double m = 0.5 * (l.magnitude[i - 1] + l.magnitude[i]);
            const double t = colors.scale.position(m);
            const std::size_t b = std::isfinite(t) ? std::min(bins - 1, static_cast<std::size_t>(t * bins)) : bins;
            batch[b].points.push_back(um(l.points[i - 1]));
            batch[b].points.push_back(um(l.points[i]));
        }
    }
    for (LineBatch& b : batch) {
        if (!b.points.empty()) s.lines.push_back(std::move(b));
    }
}

namespace {

// The triangles of one arrow, as polygons with their outward normals (cm).
void arrow_polygons(const FieldScene& s, const Arrow& a, double length,
                    std::vector<std::pair<std::vector<Point>, Point>>& out) {
    const Point d = a.direction;
    const double head = 0.35 * length;
    if (s.dimension == 2) {
        const double lift = 1e-3 * length;  // in front of the map (the camera looks down -z)
        const Point q{-d[1], d[0], 0.0};
        const Point c{a.position[0], a.position[1], lift};
        const Point tail = add(c, scale(d, -0.5 * length)), tip = add(c, scale(d, 0.5 * length));
        const Point neck = add(tip, scale(d, -head));
        const double w = 0.05 * length, hw = 0.18 * length;
        const Point up{0.0, 0.0, 1.0};
        out.push_back({{add(tail, scale(q, -w)), add(neck, scale(q, -w)), add(neck, scale(q, w)), add(tail, scale(q, w))}, up});
        out.push_back({{add(neck, scale(q, -hw)), tip, add(neck, scale(q, hw))}, up});
        return;
    }
    const Point c = add(a.position, scale(a.normal, 0.25 * length));
    const Point tail = add(c, scale(d, -0.5 * length)), tip = add(c, scale(d, 0.5 * length));
    const Point neck = add(tip, scale(d, -head));
    // A basis (d, u, w) around the axis.
    const Point helper = std::abs(d[0]) < 0.9 ? Point{1.0, 0.0, 0.0} : Point{0.0, 1.0, 0.0};
    const Point u = unit(cross(d, helper)), w = cross(d, u);
    constexpr int sides = 6;
    const double r = 0.04 * length, cone = 0.13 * length;
    const auto ring = [&](const Point& at, double radius, int k) {
        const double t = 2.0 * std::numbers::pi * k / sides;
        return add(at, add(scale(u, radius * std::cos(t)), scale(w, radius * std::sin(t))));
    };
    for (int k = 0; k < sides; ++k) {
        const double t = 2.0 * std::numbers::pi * (k + 0.5) / sides;
        const Point radial = add(scale(u, std::cos(t)), scale(w, std::sin(t)));
        out.push_back({{ring(tail, r, k), ring(tail, r, k + 1), ring(neck, r, k + 1), ring(neck, r, k)}, radial});
        out.push_back({{ring(neck, cone, k), ring(neck, cone, k + 1), tip}, add(radial, scale(d, cone / head))});
        out.push_back({{ring(neck, cone, k + 1), ring(neck, cone, k), neck}, scale(d, -1.0)});  // the cone's base
    }
}

}  // namespace

std::vector<Polygon> exposed_boundary(const device::Device& d, const Geometry& g) {
    std::vector<std::uint8_t> on_contact(d.mesh().node_count(), 0);
    for (const device::Contact& c : d.contacts()) {
        for (const NodeId v : c.nodes) on_contact[static_cast<std::size_t>(v)] = 1;
    }
    std::vector<Polygon> out;
    for (Polygon& f : g.boundary_faces()) {
        const bool covered = std::ranges::all_of(f.vertices, [&](const SurfaceVertex& v) {
            return on_contact[static_cast<std::size_t>(v.at.nodes[0])] != 0;
        });
        if (!covered) out.push_back(std::move(f));
    }
    return out;
}

void add_arrows(FieldScene& s, std::span<const Arrow> arrows, double length, const Rgba& color) {
    TriangleBatch batch{"arrows", Shading::flat, s.dimension == 3, false, {}};
    std::vector<std::pair<std::vector<Point>, Point>> polys;
    for (const Arrow& a : arrows) {
        polys.clear();
        arrow_polygons(s, a, length, polys);
        for (const auto& [pts, n] : polys) fan(batch, pts, {}, color, n);
    }
    if (!batch.vertices.empty()) s.triangles.push_back(std::move(batch));
}

void add_arrows(FieldScene& s, std::span<const Arrow> arrows, double length, const Colorbar& colors) {
    TriangleBatch batch{"arrows", Shading::field, s.dimension == 3, false, {}};
    std::vector<std::pair<std::vector<Point>, Point>> polys;
    for (const Arrow& a : arrows) {
        polys.clear();
        arrow_polygons(s, a, length, polys);
        const float t = static_cast<float>(colors.scale.position(a.magnitude));
        for (const auto& [pts, n] : polys) fan(batch, pts, std::vector<float>(pts.size(), t), {}, n);
    }
    if (!batch.vertices.empty()) s.triangles.push_back(std::move(batch));
}

std::expected<FieldScene, base::Error> outline_view(const device::Device& d, const Geometry& g,
                                                    const ViewOptions& o) {
    if (auto c = check_view(d, g); !c) return std::unexpected(std::move(c.error()));
    if (g.dimension() != 3) return std::unexpected(invalid("the outline view needs a 3D geometry"));
    FieldScene s;
    s.dimension = 3;
    s.title = "Streamlines";
    s.triangles.push_back(outline_shell(g));
    if (o.contacts) contact_slabs(d, g, o.palette, s.triangles, o.clip);
    frame(s, o);
    return s;
}

std::expected<std::vector<FieldScene>, base::Error> filmstrip_scenes(const device::Device& d, const Geometry& g,
                                                                    std::span<const FilmFrame> frames,
                                                                    ScaleKind kind, const FieldOptions& o) {
    if (frames.empty()) return std::unexpected(invalid("a filmstrip needs a frame"));
    FieldOptions shared = o;
    if (!shared.scale) {
        // CQ2: one scale over every frame, so a color means the same value in each.
        std::vector<double> all;
        for (const FilmFrame& f : frames) all.insert(all.end(), f.values.begin(), f.values.end());
        auto s = ColorScale::fit(kind, all, 1e10);
        if (!s) return std::unexpected(std::move(s.error()));
        shared.scale = *s;
    }
    std::vector<FieldScene> out;
    for (std::size_t k = 0; k < frames.size(); ++k) {
        auto s = field_view(d, g, frames[k].values, shared);
        if (!s) return std::unexpected(std::move(s.error()));
        s->title = o.title.empty() ? frames[k].label : o.title + ", " + frames[k].label;
        out.push_back(std::move(*s));
    }
    return out;
}

std::string time_label(double t) {
    if (t == 0.0) return "t = 0 s";
    struct Unit {
        double scale;
        const char* name;
    };
    constexpr std::array<Unit, 6> units{{{1.0, "s"}, {1e-3, "ms"}, {1e-6, "\xC2\xB5s"}, {1e-9, "ns"}, {1e-12, "ps"}, {1e-15, "fs"}}};
    const double a = std::abs(t);
    std::size_t u = 0;
    while (u + 1 < units.size() && a < units[u].scale * (1.0 - 5e-4)) ++u;
    const double v = t / units[u].scale;
    // Three significant digits.
    const int decimals = std::abs(v) >= 99.95 ? 0 : std::abs(v) >= 9.995 ? 1 : 2;
    char buffer[48];
    std::snprintf(buffer, sizeof buffer, "t = %.*f %s", decimals, v, units[u].name);
    return buffer;
}

double default_gap(const device::Device& d, int axis) noexcept {
    if (axis < 0 || axis > 2) return 0.0;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (const Point& p : d.mesh().points()) {
        lo = std::min(lo, p[static_cast<std::size_t>(axis)]);
        hi = std::max(hi, p[static_cast<std::size_t>(axis)]);
    }
    return 0.12 * (hi - lo);
}

std::expected<FieldScene, base::Error> exploded_view(const device::Device& d, const Geometry& g, int axis,
                                                     double gap, std::span<const double> values,
                                                     const FieldOptions& o) {
    if (auto c = check_view(d, g); !c) return std::unexpected(std::move(c.error()));
    if (g.dimension() != 3) return std::unexpected(invalid("the exploded view needs a 3D geometry"));
    if (axis < 0 || axis > 2) return std::unexpected(invalid("the axis is 0, 1 or 2", static_cast<std::size_t>(std::max(axis, 0))));
    if (!(gap >= 0.0) || !std::isfinite(gap)) return std::unexpected(invalid("the gap is a finite length, at least 0"));
    if (!values.empty() && values.size() != d.mesh().node_count()) {
        return std::unexpected(invalid("the field is not one value per mesh node", values.size()));
    }
    std::optional<Colorbar> bar;
    if (!values.empty()) {
        auto b = colorbar_for(values, o);
        if (!b) return std::unexpected(std::move(b.error()));
        bar = std::move(*b);
    }
    const auto regions = regions_of(d);
    const auto points = d.mesh().points();
    const auto ax = static_cast<std::size_t>(axis);
    // The regions' ranks along the axis, by their mean node coordinate.
    const std::size_t R = d.regions().size();
    std::vector<double> mean(R, 0.0);
    std::vector<std::size_t> count(R, 0);
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto r = static_cast<std::size_t>(regions[i]);
        mean[r] += points[i][ax];
        ++count[r];
    }
    std::vector<std::size_t> order;
    for (std::size_t r = 0; r < R; ++r) {
        if (count[r] > 0) {
            mean[r] /= static_cast<double>(count[r]);
            order.push_back(r);
        }
    }
    std::ranges::stable_sort(order, {}, [&](std::size_t r) { return mean[r]; });
    std::vector<Point> offset(R, Point{});
    for (std::size_t k = 0; k < order.size(); ++k) offset[order[k]][ax] = gap * static_cast<double>(k);

    FieldScene s;
    s.dimension = 3;
    s.title = bar ? o.title : "Exploded view";
    const ViewOptions& view = o.view;
    Batches batches(view.palette, true, view.transparent_insulators);
    TriangleBatch field_batch{"field", Shading::field, true, false, {}};
    const auto surfaces = g.region_surfaces(regions);
    for (std::size_t r = 0; r < surfaces.size(); ++r) {
        for (const Geometry::RegionFace& f : surfaces[r]) {
            std::vector<Point> pts;
            for (const Point& x : f.points) pts.push_back(add(x, offset[r]));
            const Point n = polygon_normal(f.points);
            if (bar && !d.is_insulator(f.owner)) {
                std::vector<float> vals;
                for (const Point& x : f.points) vals.push_back(position_at(g, f.element, x, values, bar->scale));
                fan(field_batch, pts, vals, {}, n);
            } else {
                const Material m = material_class(d, f.owner, view.palette);
                fan(batches[m], pts, {}, material_color(m, view.palette), n);
            }
        }
    }
    if (!field_batch.vertices.empty()) s.triangles.push_back(std::move(field_batch));
    batches.move_into(s.triangles);
    if (view.contacts) {
        // Each contact moves with the region of its first node.
        std::vector<Point> contact_offset;
        for (const device::Contact& c : d.contacts()) {
            contact_offset.push_back(c.nodes.empty() ? Point{} : offset[static_cast<std::size_t>(regions[static_cast<std::size_t>(c.nodes[0])])]);
        }
        contact_slabs(d, g, view.palette, s.triangles, std::nullopt, contact_offset);
    }
    if (bar) s.colorbar = std::move(*bar);
    frame(s, view);
    return s;
}

// --- plots ---------------------------------------------------------------------------------------

namespace {

// Tableau 10, the qualitative set the curves cycle through.
constexpr std::array<Rgb, 10> qualitative{{{0.122f, 0.467f, 0.706f}, {1.0f, 0.498f, 0.055f},
                                           {0.173f, 0.627f, 0.173f}, {0.839f, 0.153f, 0.157f},
                                           {0.580f, 0.404f, 0.741f}, {0.549f, 0.337f, 0.294f},
                                           {0.890f, 0.467f, 0.761f}, {0.498f, 0.498f, 0.498f},
                                           {0.737f, 0.741f, 0.133f}, {0.090f, 0.745f, 0.812f}}};

std::expected<Axis, base::Error> make_axis(std::string label, AxisScale kind,
                                           const std::vector<double>& data) {
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (const double v : data) {
        if (!std::isfinite(v) || (kind == AxisScale::log10 && !(v > 0.0))) continue;
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    if (!(lo <= hi)) {
        return std::unexpected(invalid(kind == AxisScale::log10 ? "a log axis needs a positive value"
                                                                : "an axis needs a finite value"));
    }
    Axis a{std::move(label), kind, lo, hi, {}};
    if (kind == AxisScale::log10) {
        const auto s = ColorScale::fit(ScaleKind::log10, data);
        a.lo = s->lo();
        a.hi = s->hi();
        a.ticks = log_ticks(a.lo, a.hi, 6);
    } else {
        if (lo == hi) {
            const double w = std::max(std::abs(lo) * 0.5, 1.0);
            lo -= w;
            hi += w;
        }
        const double pad = 0.04 * (hi - lo);
        a.lo = lo - pad;
        a.hi = hi + pad;
        a.ticks = linear_ticks(a.lo, a.hi, 6);
    }
    return a;
}

}  // namespace

std::expected<PlotScene, base::Error> curve_plot(std::span<const SeriesInput> input,
                                                 std::string x_label, std::string y_label,
                                                 AxisScale x_scale, AxisScale y_scale,
                                                 std::string title) {
    if (input.empty()) return std::unexpected(invalid("a plot needs a series"));
    PlotScene p;
    p.title = std::move(title);
    p.legend = input.size() > 1;
    std::vector<double> xs, ys;
    for (std::size_t k = 0; k < input.size(); ++k) {
        const analysis::Curve& c = input[k].curve;
        Series s;
        s.name = input[k].name;
        s.color = qualitative[k % qualitative.size()];
        s.x = c.x;
        s.y = c.y;
        bool all_negative = y_scale == AxisScale::log10;
        for (const double v : c.y) all_negative = all_negative && v < 0.0;
        if (all_negative) {
            for (double& v : s.y) v = -v;
        }
        s.markers = c.x.size() <= 60;
        s.flagged.assign(c.x.size(), 0);
        for (std::size_t i = 0; i < c.resolution.size() && i < c.y.size(); ++i) {
            s.flagged[i] = std::abs(c.y[i]) <= c.resolution[i] ? 1 : 0;
        }
        xs.insert(xs.end(), s.x.begin(), s.x.end());
        ys.insert(ys.end(), s.y.begin(), s.y.end());
        p.series.push_back(std::move(s));
    }
    auto x = make_axis(std::move(x_label), x_scale, xs);
    if (!x) return std::unexpected(std::move(x.error()));
    auto y = make_axis(std::move(y_label), y_scale, ys);
    if (!y) return std::unexpected(std::move(y.error()));
    p.x = std::move(*x);
    p.y = std::move(*y);
    return p;
}

namespace {

// One curve of y against x (um), as a profile draws it; a NaN breaks the line.
std::expected<PlotScene, base::Error> profile_against(std::vector<double> x, std::vector<double> y,
                                                      std::string x_label, std::string y_label, AxisScale y_scale) {
    PlotScene p;
    p.legend = false;
    Series s;
    s.color = qualitative[0];
    s.x = std::move(x);
    s.y = std::move(y);
    std::vector<double> xs = s.x, ys = s.y;
    p.series.push_back(std::move(s));
    auto xa = make_axis(std::move(x_label), AxisScale::linear, xs);
    if (!xa) return std::unexpected(std::move(xa.error()));
    auto ya = make_axis(std::move(y_label), y_scale, ys);
    if (!ya) return std::unexpected(std::move(ya.error()));
    p.x = std::move(*xa);
    p.y = std::move(*ya);
    return p;
}

// The four band edges against x (um), as a band diagram draws them.
std::expected<PlotScene, base::Error> bands_against(const std::vector<double>& x,
                                                    const std::array<std::vector<double>, 4>& bands,
                                                    std::string x_label) {
    PlotScene p;
    p.title = "Band diagram";
    constexpr std::array<std::string_view, 4> names{"Ec", "Ev", "EFn", "EFp"};
    const std::array<Rgb, 4> colors{{{0.12f, 0.47f, 0.71f}, {0.84f, 0.15f, 0.16f},
                                     {0.17f, 0.63f, 0.17f}, {1.0f, 0.5f, 0.06f}}};
    std::vector<double> xs, ys;
    for (std::size_t k = 0; k < bands.size(); ++k) {
        Series s;
        s.name = std::string(names[k]);
        s.color = colors[k];
        s.x = x;
        s.y = bands[k];  // NaN on insulators: a gap
        xs.insert(xs.end(), s.x.begin(), s.x.end());
        ys.insert(ys.end(), s.y.begin(), s.y.end());
        p.series.push_back(std::move(s));
    }
    auto xa = make_axis(std::move(x_label), AxisScale::linear, xs);
    if (!xa) return std::unexpected(std::move(xa.error()));
    auto ya = make_axis("E (eV)", AxisScale::linear, ys);
    if (!ya) return std::unexpected(std::move(ya.error()));
    p.x = std::move(*xa);
    p.y = std::move(*ya);
    return p;
}

std::expected<std::array<std::vector<double>, 4>, base::Error> band_parts(const results::BandDiagram& b, std::size_t n) {
    std::array<std::vector<double>, 4> parts{b.conduction_eV, b.valence_eV, b.electron_fermi_eV, b.hole_fermi_eV};
    for (const auto& v : parts) {
        if (v.size() != n) return std::unexpected(invalid("the bands are not one value per mesh node", v.size()));
    }
    return parts;
}

std::vector<double> node_x_um(const mesh::Mesh& m) {
    std::vector<double> x;
    for (const Point& p : m.points()) x.push_back(p[0] * um_per_cm);
    return x;
}

// The distance (um) at each sample, with a NaN point where the line runs outside the mesh, and the values
// to go with them.
std::vector<double> cut_x_um(const CutLine& cut) {
    std::vector<double> x;
    for (std::size_t i = 0; i < cut.samples.size(); ++i) {
        if (cut.samples[i].gap_before) x.push_back(0.5 * (cut.samples[i - 1].distance + cut.samples[i].distance) * um_per_cm);
        x.push_back(cut.samples[i].distance * um_per_cm);
    }
    return x;
}

std::vector<double> cut_y(const CutLine& cut, std::span<const double> field) {
    const std::vector<double> v = cut_values(cut, field);
    std::vector<double> y;
    for (std::size_t i = 0; i < cut.samples.size(); ++i) {
        if (cut.samples[i].gap_before) y.push_back(std::numeric_limits<double>::quiet_NaN());
        y.push_back(v[i]);
    }
    return y;
}

constexpr const char* distance_label = "distance (\xC2\xB5m)";

}  // namespace

std::expected<PlotScene, base::Error> profile_plot(const mesh::Mesh& m, std::span<const double> values,
                                                   std::string y_label, AxisScale y_scale) {
    if (m.dimension() != 1) return std::unexpected(invalid("a profile needs a 1D mesh"));
    if (values.size() != m.node_count()) {
        return std::unexpected(invalid("the field is not one value per mesh node", values.size()));
    }
    return profile_against(node_x_um(m), {values.begin(), values.end()}, "x (\xC2\xB5m)", std::move(y_label), y_scale);
}

std::expected<PlotScene, base::Error> band_plot(const mesh::Mesh& m, const results::BandDiagram& b) {
    if (m.dimension() != 1) return std::unexpected(invalid("a band diagram needs a 1D mesh"));
    auto parts = band_parts(b, m.node_count());
    if (!parts) return std::unexpected(std::move(parts.error()));
    return bands_against(node_x_um(m), *parts, "x (\xC2\xB5m)");
}

std::expected<PlotScene, base::Error> cut_plot(const mesh::Mesh& m, const CutLine& cut, std::span<const double> values,
                                               std::string y_label, AxisScale y_scale) {
    if (values.size() != m.node_count()) {
        return std::unexpected(invalid("the field is not one value per mesh node", values.size()));
    }
    if (cut.samples.empty()) return std::unexpected(invalid("the cut line misses the device"));
    return profile_against(cut_x_um(cut), cut_y(cut, values), distance_label, std::move(y_label), y_scale);
}

std::expected<PlotScene, base::Error> cut_band_plot(const mesh::Mesh& m, const CutLine& cut,
                                                    const results::BandDiagram& b) {
    auto parts = band_parts(b, m.node_count());
    if (!parts) return std::unexpected(std::move(parts.error()));
    if (cut.samples.empty()) return std::unexpected(invalid("the cut line misses the device"));
    std::array<std::vector<double>, 4> along;
    for (std::size_t k = 0; k < 4; ++k) along[k] = cut_y(cut, (*parts)[k]);
    return bands_against(cut_x_um(cut), along, distance_label);
}

}  // namespace NiTCAD::render
