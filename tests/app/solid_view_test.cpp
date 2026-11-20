// Unit 26l slice 2 (SD4): the 3D modeller's view without Win32. SE6: picking solids, faces, edges, vertices, gizmo
// handles and contacts at 96, 144 and 192 dpi, zooms from 1e-2 to 1e4 px/µm, the six named views and two oblique
// ones, as a table; the picking frame against the renderer's own projection; the drags' inverses; triangulation.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <format>
#include <string>
#include <vector>

#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/app/solid_view.hpp"
#include "NiTCAD/render/axes3d.hpp"
#include "NiTCAD/render/image.hpp"

namespace {

using namespace NiTCAD;
using namespace NiTCAD::app::solid;

// A cube domain of side L with, on its face `face` (0..5: x_min, x_max, y_min, ...): a box B (layer 1) and an
// L-shaped prism P along the face's normal (layer 2), each reaching a quarter of L into the domain, and a contact C;
// in the face's (u, v), in eighths of L: B over [1, 3]², P the L of (5,1) (7,1) (7,2) (6,2) (6,3) (5,3), C [5, 7]
// × [5, 7]. The background is layer 0.
struct Scene {
    Model m;
    int axis = 0;
    bool high = false;
    Coord L = 0, c = 0, depth = 0;
    [[nodiscard]] P3 at(double u8, double v8, Coord into = 0) const {
        const Coord w = high ? c - into : c + into;
        return from_uv(axis, {static_cast<Coord>(std::llround(u8 * static_cast<double>(L) / 8)), static_cast<Coord>(std::llround(v8 * static_cast<double>(L) / 8))}, w);
    }
};

Scene scene_on(int face, Coord L) {
    Scene s;
    s.axis = face / 2;
    s.high = face % 2 == 1;
    s.L = L;
    s.c = s.high ? L : 0;
    s.depth = L / 4;
    s.m.mesh = MeshKind::layered;
    s.m.domain = {{0, 0, 0}, {L, L, L}};
    s.m.regions.push_back({"bulk", "silicon", {}, -1});
    Solid b;
    b.kind = SolidKind::box;
    b.box = normalized(s.at(1, 1), s.at(3, 3, s.depth));
    s.m.regions.push_back({"well", "silicon", {b}, -1});
    Solid p;
    p.kind = SolidKind::prism;
    p.axis = s.axis;
    p.from = s.high ? L - s.depth : 0;
    p.to = s.high ? L : s.depth;
    for (const auto [u, v] : std::vector<std::pair<int, int>>{{5, 1}, {7, 1}, {7, 2}, {6, 2}, {6, 3}, {5, 3}}) p.polygon.push_back({u * L / 8, v * L / 8});
    s.m.regions.push_back({"fin", "sio2", {p}, -1});
    s.m.contacts.push_back({"gate", "ohmic", normalized(s.at(5, 5), s.at(7, 7)), true, -1});
    return s;
}

V3 add(V3 a, V3 b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
V3 mul(V3 a, double k) { return {a[0] * k, a[1] * k, a[2] * k}; }
double dist(std::array<double, 2> a, std::array<double, 2> b) { return std::hypot(a[0] - b[0], a[1] - b[1]); }
double seg(std::array<double, 2> p, std::array<double, 2> a, std::array<double, 2> b) {
    const double dx = b[0] - a[0], dy = b[1] - a[1], l2 = dx * dx + dy * dy;
    const double t = l2 > 0 ? std::clamp(((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2, 0.0, 1.0) : 0.0;
    return std::hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy);
}

struct View {
    std::string name;
    float azimuth, elevation;
};

std::vector<View> views() {
    std::vector<View> out;
    for (const auto& p : render::named_views()) out.push_back({p.name, p.azimuth_deg, p.elevation_deg});  // six and the isometric
    out.push_back({"oblique", -35.0f, 20.0f});
    return out;
}

}  // namespace

TEST_CASE("solid view: picking at three dpi, six decades of zoom, eight views (SE6)", "[app][solid][view]") {
    int faces = 0, contacts = 0, backgrounds = 0, nones = 0, sections = 0, vertices = 0, face_handles = 0, axes = 0, planes = 0, edges = 0;
    for (const int dpi : {96, 144, 192}) {
        for (const double target : {1e-2, 1.0, 1e2, 1e4}) {  // px/µm at 96 dpi
            const Coord L = static_cast<Coord>(std::llround(300.0 / target * static_cast<double>(per_um) / 80.0)) * 80;
            for (const View& vw : views()) {
                ViewState v;
                v.dpi = dpi;
                v.width = 800 * dpi / 96;
                v.height = 600 * dpi / 96;
                v.camera.azimuth_deg = vw.azimuth;
                v.camera.elevation_deg = vw.elevation;
                // The face most toward the viewer.
                const Scene probe = scene_on(0, L);
                const render::ViewFrame f0 = render::view_frame(frame_scene(probe.m, v), v.width, v.height);
                int face = 0;
                double best = -2;
                for (int k = 0; k < 6; ++k) {
                    const double d = (k % 2 ? 1.0 : -1.0) * -f0.forward[static_cast<std::size_t>(k / 2)];
                    if (d > best + 1e-9) {
                        best = d;
                        face = k;
                    }
                }
                const Scene s = scene_on(face, L);
                // The zoom giving `target` (scaled by dpi) px/µm at the centre.
                const V3 centre = mul(um_of(s.m.domain.hi), 0.5);
                const double upp = render::um_per_pixel(frame_scene(s.m, v), v.width, v.height, centre);
                v.camera.zoom = static_cast<float>(target * dpi / 96.0 * upp);
                const std::string where = std::format("{} dpi, {} px/um, {} view, face {}", dpi, target, vw.name, face);
                INFO(where);
                const auto px = [&](P3 p) { return screen_of(s.m, v, um_of(p)); };
                const auto at = [&](std::array<double, 2> p, const Selection& sel = {}) { return pick(s.m, v, p[0], p[1], sel); };

                // Faces on the front face: B, P's leg (exact polygon test), the L's notch (the background), the contact.
                const Pick on_b = at(px(s.at(2, 2)));
                CHECK(on_b.kind == PickKind::face);
                CHECK(on_b.solid == SolidRef{1, 0});
                CHECK(on_b.index == face);  // a box's face, as move_face numbers them
                CHECK(on_b.domain_face == face);
                const Pick on_p = at(px(s.at(5.5, 2.5)));
                CHECK(on_p.kind == PickKind::face);
                CHECK(on_p.solid == SolidRef{2, 0});
                CHECK(on_p.index == (s.high ? 1 : 0));  // the cap on the face
                const Pick notch = at(px(s.at(6.5, 2.5)));
                CHECK(notch.kind == PickKind::face);
                CHECK(notch.solid == SolidRef{0, -1});
                CHECK(notch.domain_face == face);
                const Pick on_c = at(px(s.at(6, 6)));
                CHECK(on_c.kind == PickKind::contact);
                CHECK(on_c.contact == 0);
                CHECK(on_c.domain_face == face);
                const Pick bg = at(px(s.at(2, 6)));
                CHECK(bg.kind == PickKind::face);
                CHECK(bg.solid == SolidRef{0, -1});
                faces += 2;
                contacts += 1;
                backgrounds += 2;
                // Off the domain.
                CHECK(at({2.0, 2.0}).kind == PickKind::none);
                ++nones;
                // A section plane cutting the nearer half of B away: B's far face, and the domain's opposite face.
                if (s.high && best > 0.999) {  // seen straight on: past the cut, the far faces
                    ViewState cut = v;
                    cut.section = {s.axis, s.c - s.depth / 2};
                    const auto q = screen_of(s.m, cut, um_of(s.at(2, 2)));
                    const Pick far = pick(s.m, cut, q[0], q[1], {});
                    CHECK(far.kind == PickKind::face);
                    CHECK(far.solid == SolidRef{1, 0});
                    CHECK(far.index == face - 1);
                    const auto r = screen_of(s.m, cut, um_of(s.at(2, 6)));
                    const Pick back = pick(s.m, cut, r[0], r[1], {});
                    CHECK(back.solid == SolidRef{0, -1});
                    CHECK(back.domain_face == face - 1);
                    sections += 2;
                }

                // P selected: its handles and edges, where nothing else of an earlier kind lies near on the screen.
                const Selection sel{{{2, 0}}, -1};
                const double scale = dpi / 96.0;
                const auto hs = handles(s.m, v, {2, 0});
                std::vector<std::array<double, 2>> points;  // vertex and face handles' centres
                for (const Handle& h : hs) {
                    if (h.kind == HandleKind::vertex || h.kind == HandleKind::face) points.push_back(screen_of(s.m, v, h.pts[0]));
                }
                const auto clear_of = [&](std::array<double, 2> p, double r, int skip) {
                    for (std::size_t k = 0; k < points.size(); ++k) {
                        if (static_cast<int>(k) != skip && dist(points[k], p) <= r) return false;
                    }
                    return true;
                };
                int point = 0;
                for (const Handle& h : hs) {
                    if (h.kind != HandleKind::vertex && h.kind != HandleKind::face) continue;
                    const auto p = screen_of(s.m, v, h.pts[0]);
                    const int me = point++;
                    if (!clear_of(p, 14.0 * scale, me)) continue;
                    // 5 px off the centre, still inside the 7 px.
                    const Pick k = at({p[0] + 3.0 * scale, p[1] - 4.0 * scale}, sel);
                    CHECK(k.kind == PickKind::handle);
                    CHECK(k.handle == h.kind);
                    CHECK(k.index == h.index);
                    (h.kind == HandleKind::vertex ? vertices : face_handles) += 1;
                    // And 9 px off, outside it: not this handle.
                    const Pick off = at({p[0] + 9.0 * scale, p[1]}, sel);
                    if (clear_of({p[0] + 9.0 * scale, p[1]}, 14.0 * scale, -1)) CHECK_FALSE((off.kind == PickKind::handle && off.handle == h.kind && off.index == h.index));
                }
                for (const Handle& h : hs) {
                    if (h.kind == HandleKind::axis) {
                        const auto a = screen_of(s.m, v, h.pts[0]), b = screen_of(s.m, v, h.pts[1]);
                        if (dist(a, b) < 30.0 * scale) continue;  // seen end on
                        const std::array<double, 2> p{a[0] + 0.75 * (b[0] - a[0]), a[1] + 0.75 * (b[1] - a[1])};
                        if (!clear_of(p, 12.0 * scale, -1)) continue;
                        const Pick k = at(p, sel);
                        CHECK(k.kind == PickKind::handle);
                        CHECK(k.handle == HandleKind::axis);
                        CHECK(k.index == h.index);
                        ++axes;
                    } else if (h.kind == HandleKind::plane) {
                        std::array<std::array<double, 2>, 4> q{};
                        for (int k = 0; k < 4; ++k) q[static_cast<std::size_t>(k)] = screen_of(s.m, v, h.pts[static_cast<std::size_t>(k)]);
                        const std::array<double, 2> mid{(q[0][0] + q[2][0]) / 2, (q[0][1] + q[2][1]) / 2};
                        const double area = std::abs((q[2][0] - q[0][0]) * (q[3][1] - q[1][1]) - (q[3][0] - q[1][0]) * (q[2][1] - q[0][1])) / 2;
                        if (area < 60.0 * scale * scale || !clear_of(mid, 12.0 * scale, -1)) continue;
                        bool near_axis = false;
                        for (const Handle& g : hs) {
                            if (g.kind == HandleKind::axis && seg(mid, screen_of(s.m, v, g.pts[0]), screen_of(s.m, v, g.pts[1])) <= 8.0 * scale) near_axis = true;
                        }
                        if (near_axis) continue;
                        const Pick k = at(mid, sel);
                        CHECK(k.kind == PickKind::handle);
                        CHECK(k.handle == HandleKind::plane);
                        CHECK(k.index == h.index);
                        ++planes;
                    }
                }
                // Edges: 30 % along each, clear of the handles; the edge picked runs within 4 px.
                const Solid& sp = s.m.regions[2].solids[0];
                const auto vs = sp.vertices();
                for (const auto& [i, j] : sp.edges()) {
                    const auto a = px(vs[static_cast<std::size_t>(i)]), b = px(vs[static_cast<std::size_t>(j)]);
                    if (dist(a, b) < 20.0 * scale) continue;
                    const std::array<double, 2> p{a[0] + 0.3 * (b[0] - a[0]), a[1] + 0.3 * (b[1] - a[1])};
                    if (!clear_of(p, 12.0 * scale, -1)) continue;
                    bool handle_near = false;
                    for (const Handle& g : hs) {
                        if (g.kind == HandleKind::axis && seg(p, screen_of(s.m, v, g.pts[0]), screen_of(s.m, v, g.pts[1])) <= 8.0 * scale) handle_near = true;
                        if (g.kind == HandleKind::plane) {
                            for (int k = 0; k < 4; ++k) {
                                if (dist(p, screen_of(s.m, v, g.pts[static_cast<std::size_t>(k)])) <= 40.0 * scale) handle_near = true;
                            }
                        }
                    }
                    if (handle_near) continue;
                    const Pick k = at(p, sel);
                    CHECK(k.kind == PickKind::edge);
                    if (k.kind == PickKind::edge) {
                        const auto [ki, kj] = sp.edges()[static_cast<std::size_t>(k.index)];
                        CHECK(seg(p, px(vs[static_cast<std::size_t>(ki)]), px(vs[static_cast<std::size_t>(kj)])) <= 4.0 * scale + 1e-9);
                    }
                    ++edges;
                }
            }
        }
    }
    // The table reached every kind many times.
    CHECK(faces >= 96);
    CHECK(contacts >= 96);
    CHECK(sections >= 24);
    CHECK(vertices >= 96);
    CHECK(face_handles >= 48);
    CHECK(axes >= 96);
    CHECK(planes >= 48);
    CHECK(edges >= 96);
    (void)backgrounds;
    (void)nones;
}

TEST_CASE("solid view: the picking frame is the renderer's projection; drags invert it (SE6)", "[app][solid][view]") {
    const Scene s = scene_on(5, 4 * per_um);
    for (const bool perspective : {false, true}) {
        for (const View& vw : views()) {
            ViewState v;
            v.camera.azimuth_deg = vw.azimuth;
            v.camera.elevation_deg = vw.elevation;
            v.camera.zoom = 1.7f;
            v.camera.pan = {0.13f, -0.21f};
            if (perspective) v.camera.projection = render::Projection::perspective;
            INFO(vw.name << (perspective ? " perspective" : " orthographic"));
            const render::FieldScene fs = frame_scene(s.m, v);
            for (const P3 corner : s.m.regions[2].solids[0].vertices()) {
                const V3 p = um_of(corner);
                const auto mine = screen_of(s.m, v, p);
                const auto theirs = render::screen_position(fs, v.width, v.height, {static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2])});
                CHECK(std::abs(mine[0] - theirs.x) < 0.02);
                CHECK(std::abs(mine[1] - theirs.y) < 0.02);
                // The ray under the pixel passes through the point.
                const auto r = ray_at(s.m, v, mine[0], mine[1]);
                const V3 d{p[0] - r.origin[0], p[1] - r.origin[1], p[2] - r.origin[2]};
                const double t = d[0] * r.direction[0] + d[1] * r.direction[1] + d[2] * r.direction[2];
                CHECK(std::hypot(d[0] - t * r.direction[0], d[1] - t * r.direction[1], d[2] - t * r.direction[2]) < 1e-8);
                // on_plane and along_axis give the point back.
                for (int a = 0; a < 3; ++a) {
                    if (const auto q = on_plane(s.m, v, mine[0], mine[1], a, p[static_cast<std::size_t>(a)])) {
                        CHECK(std::hypot((*q)[0] - p[0], (*q)[1] - p[1], (*q)[2] - p[2]) < 1e-7);
                    }
                    const V3 origin = add(p, {a == 0 ? -0.5 : 0.0, a == 1 ? -0.5 : 0.0, a == 2 ? -0.5 : 0.0});
                    const double front = std::abs(render::view_frame(fs, v.width, v.height).forward[static_cast<std::size_t>(a)]);
                    if (front < 0.99) CHECK(std::abs(along_axis(s.m, v, mine[0], mine[1], origin, a) - p[static_cast<std::size_t>(a)]) < 1e-7);
                }
            }
        }
    }
}

TEST_CASE("solid view: triangulation of rectilinear and collinear polygons", "[app][solid][view]") {
    const auto area2 = [](const std::vector<Pt>& p) {
        long double a = 0;
        for (std::size_t k = 0; k < p.size(); ++k) {
            const Pt q = p[k], r = p[(k + 1) % p.size()];
            a += static_cast<long double>(q.x) * r.y - static_cast<long double>(r.x) * q.y;
        }
        return std::abs(a);
    };
    const std::vector<std::vector<Pt>> shapes{
        {{0, 0}, {4, 0}, {4, 1}, {1, 1}, {1, 3}, {0, 3}},                                  // L
        {{0, 0}, {3, 0}, {3, 3}, {2, 3}, {2, 1}, {1, 1}, {1, 3}, {0, 3}},                  // U
        {{0, 0}, {2, 0}, {4, 0}, {4, 2}, {4, 4}, {0, 4}},                                  // collinear vertices
        {{0, 0}, {0, 3}, {1, 3}, {1, 1}, {2, 1}, {2, 3}, {3, 3}, {3, 0}},                  // U, clockwise
        {{0, 0}, {6, 0}, {6, 1}, {5, 1}, {5, 2}, {4, 2}, {4, 3}, {0, 3}},                  // stairs
    };
    for (const auto& p : shapes) {
        const auto tris = triangulate(p);
        CHECK(tris.size() == p.size() - 2);
        long double sum = 0;
        for (const auto& t : tris) {
            const Pt a = p[static_cast<std::size_t>(t[0])], b = p[static_cast<std::size_t>(t[1])], c = p[static_cast<std::size_t>(t[2])];
            const long double cr = static_cast<long double>(b.x - a.x) * (c.y - a.y) - static_cast<long double>(b.y - a.y) * (c.x - a.x);
            sum += std::abs(cr);
        }
        CHECK(sum == area2(p));
    }
}
