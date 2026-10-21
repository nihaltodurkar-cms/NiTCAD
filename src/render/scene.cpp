#include "NiTCAD/render/scene.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
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
                   std::vector<TriangleBatch>& out) {
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
    for (const device::Contact& contact : d.contacts()) {
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
                const Point pa = points[static_cast<std::size_t>(a)], pb = points[static_cast<std::size_t>(b)];
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
            for (const SurfaceVertex& v : f.vertices) base.push_back(v.position);
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

void contours(const TriangleBatch& field, const Colorbar& bar, std::vector<LineBatch>& out) {
    LineBatch lines{"contours", {0.1f, 0.1f, 0.1f, 0.85f}, 1.0f, true, {}};
    std::vector<double> levels;
    for (const Tick& t : bar.scale.ticks()) {
        if (t.major && t.position > 0.0 && t.position < 1.0) levels.push_back(t.position);
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
    for (const Polygon& face : g.boundary_faces()) {
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
            if (field_scale) {
                std::vector<float> vals;
                for (const Point& x : p.points) vals.push_back(position_at(g, p.element, x, values, *field_scale));
                fan(field_batch, p.points, vals, {}, n);
            } else {
                const Material m = material_class(d, semi, o.palette);
                fan(batches[m], p.points, {}, material_color(m, o.palette), n);
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
    }
    return {"?", "?", ScaleKind::linear};
}

std::optional<Field> find_field(std::string_view key) noexcept {
    for (int f = 0; f <= static_cast<int>(Field::acceptors); ++f) {
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
    }
    return std::unexpected(invalid("unknown field"));
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
    if (o.contacts) contact_slabs(d, g, o.palette, s.triangles);
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
    if (o.view.contacts) contact_slabs(d, g, o.view.palette, s.triangles);
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
    // The device's boundary, faint and translucent, for context.
    TriangleBatch shell{"outline", Shading::flat, true, true, {}};
    for (const Polygon& f : g.boundary_faces()) {
        std::vector<Point> pts;
        for (const SurfaceVertex& v : f.vertices) pts.push_back(v.position);
        fan(shell, pts, {}, {0.82f, 0.85f, 0.88f, 0.12f}, polygon_normal(pts));
    }
    s.triangles.push_back(std::move(shell));
    if (o.view.contacts) contact_slabs(d, g, o.view.palette, s.triangles);
    s.colorbar = std::move(*bar);
    frame(s, o.view);
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

std::expected<PlotScene, base::Error> profile_plot(const mesh::Mesh& m, std::span<const double> values,
                                                   std::string y_label, AxisScale y_scale) {
    if (m.dimension() != 1) return std::unexpected(invalid("a profile needs a 1D mesh"));
    if (values.size() != m.node_count()) {
        return std::unexpected(invalid("the field is not one value per mesh node", values.size()));
    }
    PlotScene p;
    p.legend = false;
    Series s;
    s.color = qualitative[0];
    for (std::size_t i = 0; i < values.size(); ++i) {
        s.x.push_back(m.points()[i][0] * um_per_cm);
        s.y.push_back(values[i]);
    }
    std::vector<double> xs = s.x, ys = s.y;
    p.series.push_back(std::move(s));
    auto x = make_axis("x (\xC2\xB5m)", AxisScale::linear, xs);
    if (!x) return std::unexpected(std::move(x.error()));
    auto y = make_axis(std::move(y_label), y_scale, ys);
    if (!y) return std::unexpected(std::move(y.error()));
    p.x = std::move(*x);
    p.y = std::move(*y);
    return p;
}

std::expected<PlotScene, base::Error> band_plot(const mesh::Mesh& m, const results::BandDiagram& b) {
    if (m.dimension() != 1) return std::unexpected(invalid("a band diagram needs a 1D mesh"));
    const std::size_t n = m.node_count();
    for (const auto* v : {&b.conduction_eV, &b.valence_eV, &b.electron_fermi_eV, &b.hole_fermi_eV}) {
        if (v->size() != n) return std::unexpected(invalid("the bands are not one value per mesh node", v->size()));
    }
    PlotScene p;
    p.title = "Band diagram";
    const std::array<std::pair<const std::vector<double>*, std::string_view>, 4> parts{
        {{&b.conduction_eV, "Ec"}, {&b.valence_eV, "Ev"}, {&b.electron_fermi_eV, "EFn"},
         {&b.hole_fermi_eV, "EFp"}}};
    const std::array<Rgb, 4> colors{{{0.12f, 0.47f, 0.71f}, {0.84f, 0.15f, 0.16f},
                                     {0.17f, 0.63f, 0.17f}, {1.0f, 0.5f, 0.06f}}};
    std::vector<double> xs, ys;
    for (std::size_t k = 0; k < parts.size(); ++k) {
        Series s;
        s.name = std::string(parts[k].second);
        s.color = colors[k];
        for (std::size_t i = 0; i < n; ++i) {
            s.x.push_back(m.points()[i][0] * um_per_cm);
            s.y.push_back((*parts[k].first)[i]);  // NaN on insulators: a gap
        }
        xs.insert(xs.end(), s.x.begin(), s.x.end());
        ys.insert(ys.end(), s.y.begin(), s.y.end());
        p.series.push_back(std::move(s));
    }
    auto x = make_axis("x (\xC2\xB5m)", AxisScale::linear, xs);
    if (!x) return std::unexpected(std::move(x.error()));
    auto y = make_axis("E (eV)", AxisScale::linear, ys);
    if (!y) return std::unexpected(std::move(y.error()));
    p.x = std::move(*x);
    p.y = std::move(*y);
    return p;
}

}  // namespace NiTCAD::render
