// The 3D modeller's view (solid_view.hpp; DECISIONS.md Unit 26l, SD4).
#include "NiTCAD/app/solid_view.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string_view>

namespace NiTCAD::app::solid {

namespace {

constexpr double per_um_d = static_cast<double>(per_um);

double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
V3 add(const V3& a, const V3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
V3 mul(const V3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
V3 cross(const V3& a, const V3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
V3 unit(const V3& a) {
    const double n = std::sqrt(dot(a, a));
    return n > 0 ? mul(a, 1.0 / n) : V3{0, 0, 1};
}
V3 axis_vector(int a) {
    V3 e{};
    e[static_cast<std::size_t>(a)] = 1.0;
    return e;
}
render::Float3 f3(const V3& p) { return {static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2])}; }

// A material's colour: the 2D editor's palette (geometry_page.cpp).
render::Rgba material_color(std::string_view m, float alpha) {
    const auto c = [&](float r, float g, float b) { return render::Rgba{r, g, b, alpha}; };
    if (m == "silicon" || m == "silicon_legacy") return c(0.93f, 0.86f, 0.70f);
    if (m == "sio2") return c(0.78f, 0.87f, 0.97f);
    if (m == "insulator") return c(0.86f, 0.86f, 0.93f);
    if (m == "germanium") return c(0.80f, 0.86f, 0.64f);
    if (m == "gaas") return c(0.74f, 0.89f, 0.76f);
    if (m == "ingaas") return c(0.70f, 0.83f, 0.87f);
    if (m == "sic_4h") return c(0.86f, 0.80f, 0.92f);
    if (m == "algaas") return c(0.90f, 0.78f, 0.70f);
    std::size_t h = 1469598103934665603ull;
    for (const char ch : m) h = (h ^ static_cast<unsigned char>(ch)) * 1099511628211ull;
    const float hue = static_cast<float>(h % 360) / 60.0f;
    const auto f = [&](float n) {
        const float k = std::fmod(n + hue, 6.0f);
        return 0.88f - 0.18f * std::max(0.0f, std::min({k, 4.0f - k, 1.0f}));
    };
    return c(f(5), f(3), f(1));
}

render::Rgba contact_color(std::string_view kind) {
    if (kind == "gate") return {0.55f, 0.20f, 0.70f, 1.0f};
    if (kind == "schottky") return {0.10f, 0.50f, 0.25f, 1.0f};
    if (kind == "electrode") return {0.20f, 0.30f, 0.65f, 1.0f};
    return {0.85f, 0.50f, 0.05f, 1.0f};
}

// A face of a solid: a planar polygon (µm) with its outward normal and triangles.
struct Face {
    int face = -1;  // move_face's numbering
    V3 normal{};
    std::vector<V3> poly;
    std::vector<std::array<int, 3>> tris;
};

std::vector<Face> faces_of(const Solid& s) {
    std::vector<Face> out;
    if (s.kind == SolidKind::box) {
        const V3 lo = um_of(s.box.lo), hi = um_of(s.box.hi);
        for (int f = 0; f < 6; ++f) {
            const int a = f / 2, u = (a + 1) % 3, w = (a + 2) % 3;
            const double c = f % 2 ? hi[static_cast<std::size_t>(a)] : lo[static_cast<std::size_t>(a)];
            Face face;
            face.face = f;
            face.normal = mul(axis_vector(a), f % 2 ? 1.0 : -1.0);
            for (const auto& [i, j] : std::array<std::pair<int, int>, 4>{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}}) {
                V3 p{};
                p[static_cast<std::size_t>(a)] = c;
                p[static_cast<std::size_t>(u)] = i ? hi[static_cast<std::size_t>(u)] : lo[static_cast<std::size_t>(u)];
                p[static_cast<std::size_t>(w)] = j ? hi[static_cast<std::size_t>(w)] : lo[static_cast<std::size_t>(w)];
                face.poly.push_back(p);
            }
            face.tris = {{0, 1, 2}, {0, 2, 3}};
            out.push_back(std::move(face));
        }
        return out;
    }
    const int a = s.axis, u = (a + 1) % 3, w = (a + 2) % 3;
    const std::size_t n = s.polygon.size();
    std::int64_t area2 = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const Pt p = s.polygon[i], q = s.polygon[(i + 1) % n];
        area2 += p.x * q.y - q.x * p.y;
    }
    const double orient = area2 >= 0 ? 1.0 : -1.0;
    const auto tris = triangulate(s.polygon);
    for (const int cap : {0, 1}) {
        Face face;
        face.face = cap;
        face.normal = mul(axis_vector(a), cap ? 1.0 : -1.0);
        for (const Pt& p : s.polygon) face.poly.push_back(um_of(from_uv(a, p, cap ? s.to : s.from)));
        face.tris = tris;
        out.push_back(std::move(face));
    }
    for (std::size_t k = 0; k < n; ++k) {
        const Pt p = s.polygon[k], q = s.polygon[(k + 1) % n];
        Face face;
        face.face = 2 + static_cast<int>(k);
        V3 nrm{};
        nrm[static_cast<std::size_t>(u)] = orient * static_cast<double>(q.y - p.y);
        nrm[static_cast<std::size_t>(w)] = -orient * static_cast<double>(q.x - p.x);
        face.normal = unit(nrm);
        face.poly = {um_of(from_uv(a, p, s.from)), um_of(from_uv(a, q, s.from)), um_of(from_uv(a, q, s.to)), um_of(from_uv(a, p, s.to))};
        face.tris = {{0, 1, 2}, {0, 2, 3}};
        out.push_back(std::move(face));
    }
    return out;
}

Solid box_solid(const Box3& b) { return Solid{SolidKind::box, b, 2, 0, 0, {}}; }

// Whether a point lies inside a planar polygon (boundary included, to a relative 1e-9), in the plane's two
// coordinates other than its normal's largest.
bool inside_face(const Face& f, const V3& p) {
    std::size_t drop = 0;
    for (std::size_t k = 1; k < 3; ++k) {
        if (std::abs(f.normal[k]) > std::abs(f.normal[drop])) drop = k;
    }
    const std::size_t i0 = (drop + 1) % 3, i1 = (drop + 2) % 3;
    double scale = 0.0;
    for (const V3& q : f.poly) scale = std::max({scale, std::abs(q[i0]), std::abs(q[i1])});
    const double tol = 1e-9 * std::max(scale, 1e-6);
    bool in = false;
    const std::size_t n = f.poly.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const double ax = f.poly[i][i0], ay = f.poly[i][i1], bx = f.poly[j][i0], by = f.poly[j][i1];
        const double px = p[i0], py = p[i1];
        // On the edge?
        const double ex = bx - ax, ey = by - ay, len = std::hypot(ex, ey);
        if (len > 0) {
            const double d = std::abs((px - ax) * ey - (py - ay) * ex) / len;
            const double t = ((px - ax) * ex + (py - ay) * ey) / (len * len);
            if (d <= tol && t >= -tol / len && t <= 1 + tol / len) return true;
        }
        if ((ay > py) != (by > py) && px < (bx - ax) * (py - ay) / (by - ay) + ax) in = !in;
    }
    return in;
}

// The ray's hit with a face (t > 0), if any.
std::optional<double> hit(const render::ScreenRay& r, const Face& f) {
    const double dn = dot(r.direction, f.normal);
    if (std::abs(dn) < 1e-12) return std::nullopt;
    const double t = dot(sub(f.poly[0], r.origin), f.normal) / dn;
    if (!(t > 0)) return std::nullopt;
    const V3 p = add(r.origin, mul(r.direction, t));
    if (!inside_face(f, p)) return std::nullopt;
    return t;
}

bool cut_away(const ViewState& v, const V3& p) {
    if (v.section.axis < 0) return false;
    const double at = static_cast<double>(v.section.at) / per_um_d;
    return p[static_cast<std::size_t>(v.section.axis)] > at + 1e-9 * std::max(1.0, std::abs(at));
}

// A triangle clipped to the section's kept side: a polygon of 0, 3 or 4 points.
std::vector<V3> clip(const ViewState& v, std::vector<V3> poly) {
    if (v.section.axis < 0) return poly;
    const std::size_t a = static_cast<std::size_t>(v.section.axis);
    const double at = static_cast<double>(v.section.at) / per_um_d;
    std::vector<V3> out;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const V3& p = poly[i];
        const V3& q = poly[(i + 1) % poly.size()];
        const bool pin = p[a] <= at, qin = q[a] <= at;
        if (pin) out.push_back(p);
        if (pin != qin) out.push_back(add(p, mul(sub(q, p), (at - p[a]) / (q[a] - p[a]))));
    }
    return out;
}

void add_polygon(render::TriangleBatch& b, const std::vector<V3>& poly, const V3& normal, const render::Rgba& color) {
    for (std::size_t k = 1; k + 1 < poly.size(); ++k) {
        for (const V3* p : {&poly[0], &poly[k], &poly[k + 1]}) b.vertices.push_back({f3(*p), f3(normal), 0.0f, color});
    }
}

void add_face(render::TriangleBatch& b, const ViewState& v, const Face& f, double push, const render::Rgba& color) {
    const V3 off = mul(f.normal, push);
    for (const auto& t : f.tris) {
        std::vector<V3> tri{add(f.poly[static_cast<std::size_t>(t[0])], off), add(f.poly[static_cast<std::size_t>(t[1])], off),
                            add(f.poly[static_cast<std::size_t>(t[2])], off)};
        add_polygon(b, clip(v, std::move(tri)), f.normal, color);
    }
}

void add_segment(render::LineBatch& b, const ViewState& v, V3 p, V3 q) {
    if (v.section.axis >= 0) {
        const std::size_t a = static_cast<std::size_t>(v.section.axis);
        const double at = static_cast<double>(v.section.at) / per_um_d;
        if (p[a] > at && q[a] > at) return;
        if (p[a] > at) p = add(q, mul(sub(p, q), (at - q[a]) / (p[a] - q[a])));
        else if (q[a] > at) q = add(p, mul(sub(q, p), (at - p[a]) / (q[a] - p[a])));
    }
    b.points.push_back(f3(p));
    b.points.push_back(f3(q));
}

void add_edges(render::LineBatch& b, const ViewState& v, const Solid& s) {
    const auto vs = s.vertices();
    for (const auto& [i, j] : s.edges()) add_segment(b, v, um_of(vs[static_cast<std::size_t>(i)]), um_of(vs[static_cast<std::size_t>(j)]));
}

// A square of `px` pixels facing the camera at `p`.
void add_square(render::LineBatch& b, const render::FieldScene& scene, const ViewState& v, const V3& p, double px) {
    const render::ViewFrame f = render::view_frame(scene, v.width, v.height);
    const double h = 0.5 * px * render::um_per_pixel(scene, v.width, v.height, p);
    const V3 r = mul(f.right, h), u = mul(f.up, h);
    const std::array<V3, 4> c{sub(sub(p, r), u), sub(add(p, r), u), add(add(p, r), u), add(sub(p, r), u)};
    for (std::size_t k = 0; k < 4; ++k) {
        b.points.push_back(f3(c[k]));
        b.points.push_back(f3(c[(k + 1) % 4]));
    }
}

// The domain's faces a contact's patch box covers with area, as rectangles on the boundary.
std::vector<Face> contact_faces(const Box3& d, const Box3& b) {
    std::vector<Face> out;
    for (const Face& f : faces_of(box_solid(d))) {
        const int a = f.face / 2, u = (a + 1) % 3, w = (a + 2) % 3;
        const bool hi = f.face % 2;
        if (!(hi ? b.hi[a] == d.hi[a] : b.lo[a] == d.lo[a]) || !(b.lo[u] < b.hi[u] && b.lo[w] < b.hi[w])) continue;
        Face g;
        g.face = f.face;
        g.normal = f.normal;
        const double c = static_cast<double>(hi ? d.hi[a] : d.lo[a]) / per_um_d;
        for (const auto& [i, j] : std::array<std::pair<int, int>, 4>{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}}) {
            V3 p{};
            p[static_cast<std::size_t>(a)] = c;
            p[static_cast<std::size_t>(u)] = static_cast<double>(i ? b.hi[u] : b.lo[u]) / per_um_d;
            p[static_cast<std::size_t>(w)] = static_cast<double>(j ? b.hi[w] : b.lo[w]) / per_um_d;
            g.poly.push_back(p);
        }
        g.tris = {{0, 1, 2}, {0, 2, 3}};
        out.push_back(std::move(g));
    }
    return out;
}

double layer_push(const Model& m) {
    const V3 lo = um_of(m.domain.lo), hi = um_of(m.domain.hi);
    return 1e-5 * std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-9});
}

double point_segment_px(double x, double y, const std::array<double, 2>& a, const std::array<double, 2>& b) {
    const double ex = b[0] - a[0], ey = b[1] - a[1], l2 = ex * ex + ey * ey;
    const double t = l2 > 0 ? std::clamp(((x - a[0]) * ex + (y - a[1]) * ey) / l2, 0.0, 1.0) : 0.0;
    return std::hypot(x - (a[0] + t * ex), y - (a[1] + t * ey));
}

bool inside_quad(double x, double y, const std::array<std::array<double, 2>, 4>& q) {
    bool in = false;
    for (std::size_t i = 0, j = 3; i < 4; j = i++) {
        if ((q[i][1] > y) != (q[j][1] > y) && x < (q[j][0] - q[i][0]) * (y - q[i][1]) / (q[j][1] - q[i][1]) + q[i][0]) in = !in;
    }
    return in;
}

}  // namespace

V3 um_of(P3 p) noexcept { return {static_cast<double>(p.x) / per_um_d, static_cast<double>(p.y) / per_um_d, static_cast<double>(p.z) / per_um_d}; }

P3 quanta_of(const V3& p) noexcept {
    return {static_cast<Coord>(std::llround(p[0] * per_um_d)), static_cast<Coord>(std::llround(p[1] * per_um_d)),
            static_cast<Coord>(std::llround(p[2] * per_um_d))};
}

std::vector<std::array<int, 3>> triangulate(const std::vector<Pt>& p) {
    std::vector<std::array<int, 3>> out;
    const int n = static_cast<int>(p.size());
    if (n < 3) return out;
    const auto cross = [&](int a, int b, int c) {
        return (p[static_cast<std::size_t>(b)].x - p[static_cast<std::size_t>(a)].x) * (p[static_cast<std::size_t>(c)].y - p[static_cast<std::size_t>(a)].y) -
               (p[static_cast<std::size_t>(b)].y - p[static_cast<std::size_t>(a)].y) * (p[static_cast<std::size_t>(c)].x - p[static_cast<std::size_t>(a)].x);
    };
    std::int64_t area2 = 0;
    for (int i = 0; i < n; ++i) {
        const Pt a = p[static_cast<std::size_t>(i)], b = p[static_cast<std::size_t>((i + 1) % n)];
        area2 += a.x * b.y - b.x * a.y;
    }
    if (area2 == 0) return out;
    const int o = area2 > 0 ? 1 : -1;
    std::vector<int> idx(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) idx[static_cast<std::size_t>(i)] = i;
    while (idx.size() > 3) {
        bool progress = false;
        const std::size_t m = idx.size();
        for (std::size_t i = 0; i < m && !progress; ++i) {
            const int a = idx[(i + m - 1) % m], b = idx[i], c = idx[(i + 1) % m];
            const std::int64_t cr = cross(a, b, c);
            if (cr == 0) {  // a straight run: the vertex adds no triangle
                idx.erase(idx.begin() + static_cast<std::ptrdiff_t>(i));
                progress = true;
                break;
            }
            if ((cr > 0 ? 1 : -1) != o) continue;
            bool ear = true;
            for (const int j : idx) {
                if (j == a || j == b || j == c) continue;
                const Pt q = p[static_cast<std::size_t>(j)];
                if (q == p[static_cast<std::size_t>(a)] || q == p[static_cast<std::size_t>(b)] || q == p[static_cast<std::size_t>(c)]) continue;
                const std::int64_t c1 = cross(a, b, j), c2 = cross(b, c, j), c3 = cross(c, a, j);
                if (o > 0 ? (c1 >= 0 && c2 >= 0 && c3 >= 0) : (c1 <= 0 && c2 <= 0 && c3 <= 0)) {
                    ear = false;
                    break;
                }
            }
            if (!ear) continue;
            out.push_back({a, b, c});
            idx.erase(idx.begin() + static_cast<std::ptrdiff_t>(i));
            progress = true;
        }
        if (!progress) {  // not simple after all: a fan, so something is drawn
            for (std::size_t k = 1; k + 1 < idx.size(); ++k) out.push_back({idx[0], idx[k], idx[k + 1]});
            return out;
        }
    }
    if (cross(idx[0], idx[1], idx[2]) != 0) out.push_back({idx[0], idx[1], idx[2]});
    return out;
}

render::FieldScene frame_scene(const Model& m, const ViewState& v) {
    render::FieldScene s;
    s.dimension = 3;
    s.camera = v.camera;
    V3 lo = um_of(m.domain.lo), hi = um_of(m.domain.hi);
    if (lo == hi) hi = add(lo, {1, 1, 1});
    for (std::size_t k = 0; k < 3; ++k) {
        const double pad = 0.05 * std::max(hi[k] - lo[k], 1e-6);
        lo[k] -= pad;
        hi[k] += pad;
    }
    s.lo = f3(lo);
    s.hi = f3(hi);
    s.axis_triad = true;
    s.look = render::Look::classic;
    s.ground_shadow = false;
    s.show_axis_box = false;
    return s;
}

std::array<double, 2> screen_of(const Model& m, const ViewState& v, const V3& p) {
    const render::FieldScene s = frame_scene(m, v);
    const render::ViewFrame f = render::view_frame(s, v.width, v.height);
    const V3 d = sub(p, f.centre);
    const double vx = dot(d, f.right) - f.cx, vy = dot(d, f.up) - f.cy;
    double nx = vx / f.half_w, ny = vy / f.half_h;
    if (f.perspective) {
        const double w = dot(d, f.forward) + f.depth;
        const double aspect = f.half_w / f.half_h;
        nx = vx / (f.tan_half * aspect * w);
        ny = vy / (f.tan_half * w);
    }
    return {f.left + (nx * 0.5 + 0.5) * (f.right_px - f.left), f.top + (0.5 - ny * 0.5) * (f.bottom_px - f.top)};
}

render::ScreenRay ray_at(const Model& m, const ViewState& v, double x, double y) {
    return render::screen_ray(frame_scene(m, v), v.width, v.height, x, y);
}

std::vector<Handle> handles(const Model& m, const ViewState& v, SolidRef ref) {
    std::vector<Handle> out;
    if (ref.region < 0 || ref.region >= static_cast<int>(m.regions.size()) || ref.solid < 0 ||
        ref.solid >= static_cast<int>(m.regions[static_cast<std::size_t>(ref.region)].solids.size())) {
        return out;
    }
    const Solid& s = m.regions[static_cast<std::size_t>(ref.region)].solids[static_cast<std::size_t>(ref.solid)];
    const render::FieldScene scene = frame_scene(m, v);
    const Box3 b = s.bounds();
    const V3 o = mul(add(um_of(b.lo), um_of(b.hi)), 0.5);
    const double len = 70.0 * v.dpi / 96.0 * render::um_per_pixel(scene, v.width, v.height, o);
    for (int a = 0; a < 3; ++a) out.push_back({HandleKind::axis, a, {o, add(o, mul(axis_vector(a), len)), {}, {}}});
    for (int n = 0; n < 3; ++n) {
        const V3 ea = mul(axis_vector((n + 1) % 3), len), eb = mul(axis_vector((n + 2) % 3), len);
        out.push_back({HandleKind::plane, n,
                       {add(o, add(mul(ea, 0.25), mul(eb, 0.25))), add(o, add(mul(ea, 0.45), mul(eb, 0.25))),
                        add(o, add(mul(ea, 0.45), mul(eb, 0.45))), add(o, add(mul(ea, 0.25), mul(eb, 0.45)))}});
    }
    for (const Face& f : faces_of(s)) {
        V3 c{};
        if (s.kind == SolidKind::prism && f.face < 2) {  // a cap: its bounds' centre
            const V3 lo = um_of(b.lo), hi = um_of(b.hi);
            c = mul(add(lo, hi), 0.5);
            c[static_cast<std::size_t>(s.axis)] = f.poly[0][static_cast<std::size_t>(s.axis)];
        } else {
            for (const V3& p : f.poly) c = add(c, mul(p, 1.0 / static_cast<double>(f.poly.size())));
        }
        out.push_back({HandleKind::face, f.face, {c, {}, {}, {}}});
    }
    if (s.kind == SolidKind::prism) {
        for (std::size_t k = 0; k < s.polygon.size(); ++k) out.push_back({HandleKind::vertex, static_cast<int>(k), {um_of(from_uv(s.axis, s.polygon[k], s.to)), {}, {}, {}}});
    }
    return out;
}

render::FieldScene solid_scene(const Model& m, const ViewState& v, const Selection& sel, const Preview& pv, const Overlay& ov) {
    render::FieldScene scene = frame_scene(m, v);
    const double push = layer_push(m);
    const bool one = sel.solids.size() == 1 && sel.solids[0].solid >= 0;
    const auto selected = [&](int r, int k) {
        return std::ranges::any_of(sel.solids, [&](const SolidRef& s) { return s.region == r && (s.solid == k || s.solid < 0); });
    };
    render::TriangleBatch opaque{"solids", render::Shading::flat, true, false, {}, render::Finish::matte};
    render::TriangleBatch glass{"translucent", render::Shading::flat, true, true, {}, render::Finish::matte};
    render::LineBatch edges{"edges", {0.25f, 0.25f, 0.28f, 1.0f}, 1.0f, false, {}};
    render::LineBatch picked{"selection", {1.0f, 0.55f, 0.0f, 1.0f}, 2.5f, true, {}};
    // The background: the domain box, under every solid.
    const int bg = m.background();
    const float alpha = v.translucent ? 0.28f : 1.0f;
    if (bg >= 0) {
        const auto& r = m.regions[static_cast<std::size_t>(bg)];
        const bool on = selected(bg, -1);
        for (const Face& f : faces_of(box_solid(m.domain))) {
            add_face(v.translucent && !on ? glass : opaque, v, f, 0.0, material_color(r.material, v.translucent && !on ? alpha : 1.0f));
        }
        if (on) add_edges(picked, v, box_solid(m.domain));
    }
    add_edges(edges, v, box_solid(m.domain));
    for (std::size_t r = 0; r < m.regions.size(); ++r) {
        const Region& reg = m.regions[r];
        for (std::size_t k = 0; k < reg.solids.size(); ++k) {
            const bool on = selected(static_cast<int>(r), static_cast<int>(k));
            const bool see_through = v.translucent && !on;
            const render::Rgba c = material_color(reg.material, see_through ? alpha : 1.0f);
            for (const Face& f : faces_of(reg.solids[k])) add_face(see_through ? glass : opaque, v, f, push * static_cast<double>(r + 1), c);
            add_edges(on ? picked : edges, v, reg.solids[k]);
        }
    }
    // Contacts: their patch faces on the boundary, over every layer.
    render::TriangleBatch pads{"contacts", render::Shading::flat, true, false, {}, render::Finish::metal};
    render::LineBatch pad_edges{"contact edges", {0.15f, 0.15f, 0.15f, 1.0f}, 1.5f, false, {}};
    const double over = push * static_cast<double>(m.regions.size() + 2);
    for (std::size_t c = 0; c < m.contacts.size(); ++c) {
        const Contact& ct = m.contacts[c];
        if (!ct.editable) continue;
        const bool on = sel.contact == static_cast<int>(c);
        for (const Face& f : contact_faces(m.domain, ct.box)) {
            add_face(pads, v, f, over, contact_color(ct.kind));
            render::LineBatch& b = on ? picked : pad_edges;
            for (std::size_t i = 0; i < f.poly.size(); ++i) {
                add_segment(b, v, add(f.poly[i], mul(f.normal, over)), add(f.poly[(i + 1) % f.poly.size()], mul(f.normal, over)));
            }
        }
    }
    // The work plane's grid over the domain's extent.
    if (v.grid > 0) {
        render::LineBatch grid{"grid", {0.45f, 0.62f, 0.85f, 0.55f}, 1.0f, false, {}};
        const int n = v.plane.normal, u = (n + 1) % 3, w = (n + 2) % 3;
        const Coord u0 = m.domain.lo[u], u1 = m.domain.hi[u], w0 = m.domain.lo[w], w1 = m.domain.hi[w];
        Coord step = v.grid;
        while ((u1 - u0) / step > 200 || (w1 - w0) / step > 200) step *= 2;
        const auto pt = [&](Coord a, Coord b) { return um_of(from_uv(n, {a, b}, v.plane.offset)); };
        for (Coord a = (u0 / step) * step; a <= u1; a += step) {
            if (a >= u0) add_segment(grid, v, pt(a, w0), pt(a, w1));
        }
        for (Coord b = (w0 / step) * step; b <= w1; b += step) {
            if (b >= w0) add_segment(grid, v, pt(u0, b), pt(u1, b));
        }
        scene.lines.push_back(std::move(grid));
    }
    // What a tool draws.
    render::TriangleBatch ghost{"preview", render::Shading::flat, true, true, {}, render::Finish::matte};
    render::LineBatch ghost_edges{"preview edges", {0.0f, 0.45f, 0.95f, 1.0f}, 2.0f, true, {}};
    for (const Solid& s : pv.solids) {
        for (const Face& f : faces_of(s)) add_face(ghost, v, f, 0.0, {0.3f, 0.6f, 1.0f, 0.35f});
        add_edges(ghost_edges, v, s);
    }
    for (std::size_t k = 0; k + 1 < pv.polyline.size(); ++k) add_segment(ghost_edges, v, um_of(pv.polyline[k]), um_of(pv.polyline[k + 1]));
    if (pv.measure) add_segment(ghost_edges, v, um_of(pv.measure->first), um_of(pv.measure->second));
    // The selection's handles.
    render::LineBatch gizmo[3] = {{"gizmo x", {0.85f, 0.15f, 0.15f, 1.0f}, 3.0f, true, {}},
                                  {"gizmo y", {0.15f, 0.65f, 0.15f, 1.0f}, 3.0f, true, {}},
                                  {"gizmo z", {0.15f, 0.30f, 0.90f, 1.0f}, 3.0f, true, {}}};
    render::LineBatch knobs{"handles", {0.05f, 0.05f, 0.05f, 1.0f}, 2.0f, true, {}};
    render::LineBatch corners{"vertex handles", {0.75f, 0.10f, 0.65f, 1.0f}, 2.0f, true, {}};
    if (one) {
        const double px = 7.0 * v.dpi / 96.0;
        for (const Handle& h : handles(m, v, sel.solids[0])) {
            switch (h.kind) {
            case HandleKind::axis:
                gizmo[h.index].points.push_back(f3(h.pts[0]));
                gizmo[h.index].points.push_back(f3(h.pts[1]));
                break;
            case HandleKind::plane:
                for (std::size_t i = 0; i < 4; ++i) {
                    gizmo[(h.index + 2) % 3].points.push_back(f3(h.pts[i]));
                    gizmo[(h.index + 2) % 3].points.push_back(f3(h.pts[(i + 1) % 4]));
                }
                break;
            case HandleKind::face: add_square(knobs, scene, v, h.pts[0], px); break;
            case HandleKind::vertex: add_square(corners, scene, v, h.pts[0], px); break;
            }
        }
    }
    // The live check's mesh preview and problems.
    render::LineBatch mesh_lines{"mesh", {0.10f, 0.10f, 0.12f, 0.8f}, 1.0f, false, {}};
    for (std::size_t k = 0; k + 1 < ov.mesh.size(); k += 2) add_segment(mesh_lines, v, ov.mesh[k], ov.mesh[k + 1]);
    render::LineBatch marks{"problems", {0.90f, 0.05f, 0.05f, 1.0f}, 3.0f, true, {}};
    for (const Marker& mk : ov.markers) {
        switch (mk.kind) {
        case Marker::Kind::solid: {
            if (mk.solid.region < 0 || mk.solid.region >= static_cast<int>(m.regions.size())) break;
            const Region& reg = m.regions[static_cast<std::size_t>(mk.solid.region)];
            for (std::size_t k = 0; k < reg.solids.size(); ++k) {
                if (mk.solid.solid < 0 || mk.solid.solid == static_cast<int>(k)) add_edges(marks, v, reg.solids[k]);
            }
            if (reg.solids.empty()) add_edges(marks, v, box_solid(m.domain));
            break;
        }
        case Marker::Kind::contact:
            if (mk.contact >= 0 && mk.contact < static_cast<int>(m.contacts.size())) {
                const Box3& b = m.contacts[static_cast<std::size_t>(mk.contact)].box;
                const auto fs = contact_faces(m.domain, b);
                for (const Face& f : fs) {
                    for (std::size_t i = 0; i < f.poly.size(); ++i) add_segment(marks, v, f.poly[i], f.poly[(i + 1) % f.poly.size()]);
                }
                if (fs.empty()) add_edges(marks, v, box_solid(b));  // a patch box on no boundary face: the box itself
            }
            break;
        case Marker::Kind::point: add_square(marks, scene, v, um_of(mk.at), 9.0 * v.dpi / 96.0); break;
        case Marker::Kind::domain: add_edges(marks, v, box_solid(m.domain)); break;
        }
    }
    for (render::TriangleBatch* b : {&opaque, &pads, &glass, &ghost}) {
        if (!b->vertices.empty()) scene.triangles.push_back(std::move(*b));
    }
    for (render::LineBatch* b : {&mesh_lines, &edges, &pad_edges, &picked, &ghost_edges, &marks, &gizmo[0], &gizmo[1], &gizmo[2], &knobs, &corners}) {
        if (!b->points.empty()) scene.lines.push_back(std::move(*b));
    }
    return scene;
}

Pick pick(const Model& m, const ViewState& v, double x, double y, const Selection& sel) {
    Pick out;
    const double scale = v.dpi / 96.0;
    const bool one = sel.solids.size() == 1 && sel.solids[0].solid >= 0;
    if (one) {
        const auto hs = handles(m, v, sel.solids[0]);
        const auto px = [&](const V3& p) { return screen_of(m, v, p); };
        for (const HandleKind kind : {HandleKind::vertex, HandleKind::face, HandleKind::axis, HandleKind::plane}) {
            double best = std::numeric_limits<double>::infinity();
            const Handle* chosen = nullptr;
            for (const Handle& h : hs) {
                if (h.kind != kind) continue;
                double d = std::numeric_limits<double>::infinity();
                if (kind == HandleKind::axis) {
                    const double e = point_segment_px(x, y, px(h.pts[0]), px(h.pts[1]));
                    if (e <= 5.0 * scale) d = e;
                } else if (kind == HandleKind::plane) d = inside_quad(x, y, {px(h.pts[0]), px(h.pts[1]), px(h.pts[2]), px(h.pts[3])}) ? 0.0 : d;
                else {
                    const auto p = px(h.pts[0]);
                    const double e = std::hypot(x - p[0], y - p[1]);
                    if (e <= 7.0 * scale) d = e;
                }
                if (d < best) {
                    best = d;
                    chosen = &h;
                }
            }
            if (chosen) {
                out.kind = PickKind::handle;
                out.handle = chosen->kind;
                out.index = chosen->index;
                out.solid = sel.solids[0];
                out.at = chosen->pts[0];
                return out;
            }
        }
        // The selected solid's edges.
        const Solid& s = m.regions[static_cast<std::size_t>(sel.solids[0].region)].solids[static_cast<std::size_t>(sel.solids[0].solid)];
        const auto vs = s.vertices();
        const auto es = s.edges();
        double best = 4.0 * scale;
        int edge = -1;
        for (std::size_t k = 0; k < es.size(); ++k) {
            const double d = point_segment_px(x, y, px(um_of(vs[static_cast<std::size_t>(es[k].first)])), px(um_of(vs[static_cast<std::size_t>(es[k].second)])));
            if (d <= best) {
                best = d;
                edge = static_cast<int>(k);
            }
        }
        if (edge >= 0) {
            out.kind = PickKind::edge;
            out.index = edge;
            out.solid = sel.solids[0];
            return out;
        }
    }
    // Faces along the ray: the nearest wins; a tie (one plane) goes to a contact, then the later layer.
    const render::ScreenRay r = ray_at(m, v, x, y);
    const V3 dlo = um_of(m.domain.lo), dhi = um_of(m.domain.hi);
    const double tol = 1e-9 * std::max({1.0, std::abs(dlo[0]), std::abs(dhi[0]), std::abs(dlo[1]), std::abs(dhi[1]), std::abs(dlo[2]), std::abs(dhi[2])});
    double best_t = std::numeric_limits<double>::infinity();
    int best_rank = -1;
    const auto offer = [&](double t, int rank, const Pick& p) {
        const V3 at = add(r.origin, mul(r.direction, t));
        if (cut_away(v, at)) return;
        if (t < best_t - tol || (std::abs(t - best_t) <= tol && rank > best_rank)) {
            best_t = std::min(t, best_t);
            best_rank = rank;
            out = p;
            out.at = at;
        }
    };
    const int bg = m.background();
    for (const Face& f : faces_of(box_solid(m.domain))) {
        if (const auto t = hit(r, f)) offer(*t, 0, Pick{PickKind::face, HandleKind::axis, -1, {bg, -1}, -1, -1, {}});
    }
    for (std::size_t g = 0; g < m.regions.size(); ++g) {
        for (std::size_t k = 0; k < m.regions[g].solids.size(); ++k) {
            for (const Face& f : faces_of(m.regions[g].solids[k])) {
                if (const auto t = hit(r, f)) offer(*t, 1 + static_cast<int>(g), Pick{PickKind::face, HandleKind::axis, f.face, {static_cast<int>(g), static_cast<int>(k)}, -1, -1, {}});
            }
        }
    }
    for (std::size_t c = 0; c < m.contacts.size(); ++c) {
        if (!m.contacts[c].editable) continue;
        for (const Face& f : contact_faces(m.domain, m.contacts[c].box)) {
            if (const auto t = hit(r, f)) offer(*t, 1 << 20, Pick{PickKind::contact, HandleKind::axis, -1, {}, static_cast<int>(c), -1, {}});
        }
    }
    if (out.kind != PickKind::none) {
        for (int f = 0; f < 6; ++f) {
            const std::size_t a = static_cast<std::size_t>(f / 2);
            const double c = f % 2 ? dhi[a] : dlo[a];
            if (std::abs(out.at[a] - c) <= tol) {
                out.domain_face = f;
                break;
            }
        }
    }
    return out;
}

std::optional<V3> on_plane(const Model& m, const ViewState& v, double x, double y, int normal, double at) {
    const render::ScreenRay r = ray_at(m, v, x, y);
    const double d = r.direction[static_cast<std::size_t>(normal)];
    if (std::abs(d) < 1e-9) return std::nullopt;
    const double t = (at - r.origin[static_cast<std::size_t>(normal)]) / d;
    V3 p = add(r.origin, mul(r.direction, t));
    p[static_cast<std::size_t>(normal)] = at;
    return p;
}

double along_axis(const Model& m, const ViewState& v, double x, double y, const V3& origin, int axis) {
    // The nearest points of two lines: o + s e and r.o + t r.d.
    const render::ScreenRay r = ray_at(m, v, x, y);
    const V3 e = axis_vector(axis), w0 = sub(origin, r.origin);
    const double b = dot(e, r.direction), d = dot(e, w0), ee = dot(r.direction, w0);
    const double den = 1.0 - b * b;
    const double s = den > 1e-12 ? (b * ee - d) / den : 0.0;
    return origin[static_cast<std::size_t>(axis)] + s;
}

Coord snap_tolerance(const Model& m, const ViewState& v, const V3& p) {
    const double um = 8.0 * v.dpi / 96.0 * render::um_per_pixel(frame_scene(m, v), v.width, v.height, p);
    return std::max<Coord>(1, static_cast<Coord>(std::llround(um * per_um_d)));
}

Coord view_grid(const Model& m, const ViewState& v, const V3& p) {
    const double um = render::um_per_pixel(frame_scene(m, v), v.width, v.height, p);
    return grid_step(1.0 / std::max(um, 1e-12), v.dpi);
}

}  // namespace NiTCAD::app::solid
