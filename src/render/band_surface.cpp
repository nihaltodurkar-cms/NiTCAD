// The band surface (band_surface.hpp).
#include "NiTCAD/render/band_surface.hpp"
#include "NiTCAD/render/streamlines.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <variant>

#include "scene_detail.hpp"

namespace NiTCAD::render {

std::vector<std::int32_t> band_material_keys(const device::Device& device) {
    // A key per region: insulators -1; a semiconductor the key of the first earlier region with equal
    // parameters (the solver's heterointerface test), else a new one.
    const auto regions = device.regions();
    std::vector<std::int32_t> region_key(regions.size(), -1);
    std::int32_t next = 0;
    for (std::size_t r = 0; r < regions.size(); ++r) {
        if (device::is_insulator(regions[r])) continue;
        const auto& params = std::get<physics::Semiconductor>(regions[r].material).parameters();
        for (std::size_t q = 0; q < r && region_key[r] < 0; ++q) {
            if (region_key[q] >= 0 && std::get<physics::Semiconductor>(regions[q].material).parameters() == params) {
                region_key[r] = region_key[q];
            }
        }
        if (region_key[r] < 0) region_key[r] = next++;
    }
    std::vector<std::int32_t> keys;
    keys.reserve(device.node_region().size());
    for (const auto r : device.node_region()) keys.push_back(region_key[static_cast<std::size_t>(r)]);
    return keys;
}

double band_value(const Geometry& geometry, std::size_t element, const Point& x, NodeId owner,
                  std::span<const std::int32_t> keys, std::span<const double> level) {
    const Interpolant at = geometry.interpolant(element, x);
    const std::int32_t key = keys[static_cast<std::size_t>(owner)];
    bool one = true;
    double weight = 0.0, sum = 0.0;
    for (std::size_t k = 0; k < at.nodes.size(); ++k) {
        if (at.weights[k] == 0.0) continue;
        const auto n = static_cast<std::size_t>(at.nodes[k]);
        if (keys[n] != key) {
            one = false;
            continue;
        }
        weight += at.weights[k];
        sum += at.weights[k] * level[n];
    }
    if (one) return interpolate(at, level);
    if (weight > 0.0) return sum / weight;
    return level[static_cast<std::size_t>(owner)];
}

namespace {

constexpr double um_per_cm = 1e4;
constexpr std::array<std::string_view, 4> level_names{"Ec", "Ev", "EFn", "EFp"};
// The band diagram's colours (scene.cpp bands_against).
constexpr std::array<Rgb, 4> level_colors{{{0.12f, 0.47f, 0.71f}, {0.84f, 0.15f, 0.16f}, {0.17f, 0.63f, 0.17f},
                                           {1.0f, 0.5f, 0.06f}}};

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}

Point sub(const Point& a, const Point& b) noexcept { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const Point& a, const Point& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Point cross(const Point& a, const Point& b) noexcept {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Point unit(const Point& a) noexcept {
    const double n = std::sqrt(dot(a, a));
    return n > 0.0 ? Point{a[0] / n, a[1] / n, a[2] / n} : a;
}

const std::vector<double>& level_values(const results::BandDiagram& b, std::size_t k) {
    switch (k) {
        case 0: return b.conduction_eV;
        case 1: return b.valence_eV;
        case 2: return b.electron_fermi_eV;
        default: return b.hole_fermi_eV;
    }
}

// The base plane's axes: u, v (unit, world) and their titles.
struct Base {
    Point u{1.0, 0.0, 0.0};
    Point v{0.0, 1.0, 0.0};
    std::string u_title = "x (\xC2\xB5m)";
    std::string v_title = "y (\xC2\xB5m)";
};

Base base_of(const std::optional<Plane>& plane) {
    if (!plane) return {};
    static constexpr std::array<const char*, 3> names{"x (\xC2\xB5m)", "y (\xC2\xB5m)", "z (\xC2\xB5m)"};
    const Point n = unit(plane->normal);
    for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t a = k == 0 ? 1 : 0, b = k == 2 ? 1 : 2;
        if (std::abs(n[a]) < 1e-12 && std::abs(n[b]) < 1e-12) {
            Base out;
            out.u = {0.0, 0.0, 0.0};
            out.v = {0.0, 0.0, 0.0};
            out.u[a] = 1.0;
            out.v[b] = 1.0;
            out.u_title = names[a];
            out.v_title = names[b];
            return out;
        }
    }
    // Oblique: u the projection of the first world axis not along the normal, v = n x u.
    Base out;
    for (std::size_t j = 0; j < 3; ++j) {
        if (std::abs(n[j]) > 0.9) continue;
        Point e{0.0, 0.0, 0.0};
        e[j] = 1.0;
        const double d = dot(e, n);
        out.u = unit({e[0] - d * n[0], e[1] - d * n[1], e[2] - d * n[2]});
        break;
    }
    out.v = cross(n, out.u);
    out.u_title = "u (\xC2\xB5m)";
    out.v_title = "v (\xC2\xB5m)";
    return out;
}

// A base polygon split for one control-volume owner, its points in cm and each shown level's energy there.
struct Piece {
    std::vector<Point> points;
    NodeId owner = 0;
    std::size_t element = 0;
    std::int32_t key = 0;
    std::array<std::vector<double>, 4> energy;  // per shown level, one per point
};

// A shared edge between two pieces of one element in different materials: its two points and each side's
// energies there.
struct Wall {
    std::array<Point, 2> points;
    std::array<std::array<double, 2>, 4> a{}, b{};  // per level: side a's, side b's, at the two points
};

// The overlap of segment q0-q1 with segment p0-p1 when both lie on one line (within tol), as parameters on p.
std::optional<std::pair<double, double>> overlap(const Point& p0, const Point& p1, const Point& q0, const Point& q1,
                                                 double tol) {
    const Point d = sub(p1, p0);
    const double len2 = dot(d, d);
    if (len2 <= tol * tol) return std::nullopt;
    const auto off_line = [&](const Point& q) {
        const Point r = sub(q, p0);
        const Point c = cross(r, d);
        return std::sqrt(dot(c, c) / len2) > tol;
    };
    if (off_line(q0) || off_line(q1)) return std::nullopt;
    const double t0 = dot(sub(q0, p0), d) / len2, t1 = dot(sub(q1, p0), d) / len2;
    const double lo = std::max(0.0, std::min(t0, t1)), hi = std::min(1.0, std::max(t0, t1));
    if ((hi - lo) * std::sqrt(len2) <= tol) return std::nullopt;
    return std::pair{lo, hi};
}

Point lerp(const Point& a, const Point& b, double t) noexcept {
    return {a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t};
}

// A lifted polygon: base points and an energy each; clipped to lo <= E <= hi (Sutherland-Hodgman, E linear).
struct Lifted {
    std::vector<Point> points;
    std::vector<double> energy;
};

Lifted clip_energy(Lifted in, double lo, double hi) {
    for (const int side : {0, 1}) {
        Lifted out;
        const std::size_t n = in.points.size();
        const auto inside = [&](double e) { return side == 0 ? e >= lo : e <= hi; };
        const double bound = side == 0 ? lo : hi;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t j = (i + 1) % n;
            const double ei = in.energy[i], ej = in.energy[j];
            if (inside(ei)) {
                out.points.push_back(in.points[i]);
                out.energy.push_back(ei);
            }
            if (inside(ei) != inside(ej)) {
                const double t = (bound - ei) / (ej - ei);
                out.points.push_back(lerp(in.points[i], in.points[j], t));
                out.energy.push_back(bound);
            }
        }
        in = std::move(out);
        if (in.points.size() < 3) return {};
    }
    return in;
}

}  // namespace

NodeId band_owner(const Geometry& geometry, std::size_t element, const Point& x) {
    if (geometry.kind() == GeometryKind::tensor) {
        const std::size_t c = geometry.element_cell(element);
        const auto [lo, hi] = geometry.cell_box(c);
        std::size_t q = 0;
        for (int a = 0; a < geometry.dimension(); ++a) {
            const auto ax = static_cast<std::size_t>(a);
            if (x[ax] > 0.5 * (lo[ax] + hi[ax])) q |= std::size_t{1} << a;
        }
        return geometry.cell_corners(c)[q];
    }
    const Interpolant at = geometry.interpolant(element, x);
    std::size_t best = 0;
    for (std::size_t k = 1; k < at.nodes.size(); ++k) {
        if (at.weights[k] > at.weights[best]) best = k;
    }
    return at.nodes[best];
}

std::optional<double> band_energy_at(const device::Device& device, const Geometry& geometry,
                                     const results::BandDiagram& bands, BandLevel level, std::size_t element,
                                     const Point& x) {
    const auto& values = level_values(bands, static_cast<std::size_t>(level));
    if (values.size() != device.mesh().node_count() || element >= geometry.element_count()) return std::nullopt;
    const auto keys = band_material_keys(device);
    const NodeId owner = band_owner(geometry, element, x);
    if (keys[static_cast<std::size_t>(owner)] < 0) return std::nullopt;
    const double e = band_value(geometry, element, x, owner, keys, values);
    if (!std::isfinite(e)) return std::nullopt;
    return e;
}

std::expected<BandSurfaceScene, base::Error> band_surface_view(const device::Device& device,
                                                              const Geometry& geometry,
                                                              const results::BandDiagram& bands,
                                                              const BandSurfaceOptions& o) {
    const int dim = device.mesh().dimension();
    if (dim < 2) return std::unexpected(invalid("a band surface needs a 2D or 3D device"));
    if (auto c = detail::check_view(device, geometry, o.view); !c) return std::unexpected(std::move(c.error()));
    if (dim == 2 && o.plane) return std::unexpected(invalid("a 2D band surface takes no plane: it is the whole device"));
    if (dim == 3 && !o.plane) return std::unexpected(invalid("a 3D band surface needs a plane"));
    const std::size_t nodes = device.mesh().node_count();
    for (std::size_t k = 0; k < 4; ++k) {
        const auto& v = level_values(bands, k);
        if (o.levels[k] && v.size() != nodes) {
            return std::unexpected(invalid(std::string(level_names[k]) + " is not one value per mesh node", v.size()));
        }
    }
    std::vector<std::size_t> shown;
    for (std::size_t k = 0; k < 4; ++k) {
        if (o.levels[k]) shown.push_back(k);
    }
    if (shown.empty()) return std::unexpected(invalid("a band surface needs a level to show"));
    if (!(o.exaggeration >= 0.1 && o.exaggeration <= 10.0)) {
        return std::unexpected(invalid("the vertical exaggeration is outside 0.1 to 10"));
    }
    if (!(o.opacity >= 0.2f && o.opacity <= 1.0f)) return std::unexpected(invalid("the opacity is outside 0.2 to 1"));
    if (o.range_eV && !(std::isfinite(o.range_eV->first) && std::isfinite(o.range_eV->second) &&
                        o.range_eV->first < o.range_eV->second)) {
        return std::unexpected(invalid("the energy range needs finite lo < hi"));
    }

    // The base polygons, nudged off a mesh plane between materials (Review Focus 4).
    std::vector<Polygon> base_polygons;
    if (dim == 2) {
        for (std::size_t e = 0; e < geometry.element_count(); ++e) base_polygons.push_back(geometry.element_polygon(e));
    } else {
        base_polygons = geometry.slice(*o.plane);
        if (base_polygons.empty()) return std::unexpected(invalid("the band-surface plane misses the device"));
    }
    const Base base = base_of(o.plane);
    const auto keys = band_material_keys(device);

    BandSurfaceScene out;
    std::vector<Piece> pieces;
    std::vector<Wall> walls;
    const auto values = [&](std::size_t k) -> std::span<const double> { return level_values(bands, k); };
    for (const Polygon& p : base_polygons) {
        std::vector<Point> pts;
        for (const SurfaceVertex& s : p.vertices) pts.push_back(s.position);
        if (pts.size() < 3 || polygon_area(pts) == 0.0) {
            ++out.dropped;
            continue;
        }
        const std::size_t first = pieces.size();
        for (OwnedPolygon& owned : geometry.owner_pieces(p)) {
            const std::int32_t key = keys[static_cast<std::size_t>(owned.owner)];
            if (key < 0 || owned.points.size() < 3) continue;  // insulator: a gap
            Piece piece{std::move(owned.points), owned.owner, p.element, key, {}};
            bool finite = true;
            for (const std::size_t k : shown) {
                for (const Point& x : piece.points) {
                    const double e = band_value(geometry, p.element, x, piece.owner, keys, values(k));
                    finite = finite && std::isfinite(e);
                    piece.energy[k].push_back(e);
                }
            }
            if (!finite) {
                ++out.dropped;
                continue;
            }
            pieces.push_back(std::move(piece));
        }
        // Walls: edges shared by pieces of this polygon in different materials, once per pair.
        double size = 0.0;
        for (const Point& q : pts) size = std::max(size, std::sqrt(dot(sub(q, pts[0]), sub(q, pts[0]))));
        const double tol = 1e-9 * size;
        for (std::size_t i = first; i < pieces.size(); ++i) {
            for (std::size_t j = i + 1; j < pieces.size(); ++j) {
                const Piece& a = pieces[i];
                const Piece& b = pieces[j];
                if (a.key == b.key) continue;
                for (std::size_t ia = 0; ia < a.points.size(); ++ia) {
                    const Point& a0 = a.points[ia];
                    const Point& a1 = a.points[(ia + 1) % a.points.size()];
                    for (std::size_t ib = 0; ib < b.points.size(); ++ib) {
                        const auto t = overlap(a0, a1, b.points[ib], b.points[(ib + 1) % b.points.size()], tol);
                        if (!t) continue;
                        Wall w{{lerp(a0, a1, t->first), lerp(a0, a1, t->second)}, {}, {}};
                        for (const std::size_t k : shown) {
                            for (std::size_t e = 0; e < 2; ++e) {
                                w.a[k][e] = band_value(geometry, a.element, w.points[e], a.owner, keys, values(k));
                                w.b[k][e] = band_value(geometry, b.element, w.points[e], b.owner, keys, values(k));
                            }
                        }
                        walls.push_back(w);
                    }
                }
            }
        }
    }
    if (pieces.empty()) {
        return std::unexpected(invalid(dim == 2 ? "no band values on the device" : "no band values on the plane"));
    }

    // The range, the height scale and the colours.
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    double ulo = lo, uhi = hi, vlo = lo, vhi = hi;
    for (const Piece& p : pieces) {
        for (const std::size_t k : shown) {
            for (const double e : p.energy[k]) {
                lo = std::min(lo, e);
                hi = std::max(hi, e);
            }
        }
        for (const Point& x : p.points) {
            ulo = std::min(ulo, dot(x, base.u));
            uhi = std::max(uhi, dot(x, base.u));
            vlo = std::min(vlo, dot(x, base.v));
            vhi = std::max(vhi, dot(x, base.v));
        }
    }
    if (o.range_eV) {
        lo = o.range_eV->first;
        hi = o.range_eV->second;
    } else if (!(hi - lo > 1e-9)) {
        // One value (to rounding: renormalised weights may differ from 1 by an ulp): widened, so the heights stay
        // on the scale of the base.
        const double mid = 0.5 * (lo + hi);
        lo = mid - 0.5;
        hi = mid + 0.5;
    }
    out.range_eV = {lo, hi};
    const double extent_um = std::max(uhi - ulo, vhi - vlo) * um_per_cm;
    out.um_per_eV = 0.6 * (extent_um > 0.0 ? extent_um : 1.0) / (hi - lo) * o.exaggeration;
    auto scale = ColorScale::create(ScaleKind::linear, lo, hi);
    if (!scale) return std::unexpected(std::move(scale.error()));

    FieldScene& s = out.scene;
    s.dimension = 3;
    s.title = "Band surface";
    const double cm_per_eV = out.um_per_eV / um_per_cm;
    // Scene point (cm, as fan takes it): x = u, y = height, z = v.
    const auto lift = [&](const Point& x, double e) -> Point { return {dot(x, base.u), cm_per_eV * e, dot(x, base.v)}; };
    for (const std::size_t k : shown) {
        TriangleBatch b{std::string(level_names[k]), o.by_level ? Shading::flat : Shading::field, true,
                        o.opacity < 1.0f, {}};
        const Rgb c = level_colors[k];
        const Rgba color = o.by_level ? Rgba{c.r, c.g, c.b, o.opacity} : Rgba{0.0f, 0.0f, 0.0f, o.opacity};
        const auto emit = [&](const Lifted& l, const Point& normal) {
            if (l.points.size() < 3) return;
            std::vector<Point> pts;
            std::vector<float> vals;
            for (std::size_t i = 0; i < l.points.size(); ++i) {
                pts.push_back(lift(l.points[i], l.energy[i]));
                if (!o.by_level) vals.push_back(static_cast<float>(scale->position(l.energy[i])));
            }
            detail::fan(b, std::move(pts), vals, color, normal);
        };
        for (const Piece& p : pieces) emit(clip_energy({p.points, p.energy[k]}, lo, hi), {0.0, 1.0, 0.0});
        for (const Wall& w : walls) {
            const Lifted quad{{w.points[0], w.points[1], w.points[1], w.points[0]},
                              {w.a[k][0], w.a[k][1], w.b[k][1], w.b[k][0]}};
            if (std::abs(w.a[k][0] - w.b[k][0]) + std::abs(w.a[k][1] - w.b[k][1]) == 0.0) continue;
            // The wall stands on the base edge: its quad, drawn as two triangles in the lifted vertical plane.
            const Point du = sub(lift(w.points[1], 0.0), lift(w.points[0], 0.0));
            const Point normal = unit({du[2], 0.0, -du[0]});
            const Lifted clipped = clip_energy(quad, lo, hi);
            if (clipped.points.size() < 3) continue;
            std::vector<Point> pts;
            std::vector<float> vals;
            for (std::size_t i = 0; i < clipped.points.size(); ++i) {
                pts.push_back(lift(clipped.points[i], clipped.energy[i]));
                if (!o.by_level) vals.push_back(static_cast<float>(scale->position(clipped.energy[i])));
            }
            detail::fan(b, std::move(pts), vals, color, normal);
        }
        if (!b.vertices.empty()) s.triangles.push_back(std::move(b));
    }
    // BSD9: a cut line in the base plane, drawn on each shown level at its heights; a step where the owner's
    // material changes (each segment split in 16, its owner at the part's middle: within 1/32 of an element).
    if (o.cut && o.cut->samples.size() >= 2) {
        const CutLine& c = *o.cut;
        bool in_plane = true;
        if (o.plane) {
            const double nn = std::sqrt(dot(o.plane->normal, o.plane->normal));
            const Point n = unit(o.plane->normal);
            const double off = o.plane->offset / nn;
            const double tol = 1e-9 * std::max(extent_um / um_per_cm, 1e-30);
            in_plane = std::abs(dot(n, c.a) - off) <= tol && std::abs(dot(n, c.b) - off) <= tol;
        }
        const Point d = sub(c.b, c.a);
        const double len = std::sqrt(dot(d, d));
        if (in_plane && len > 0.0) {
            const Locator locator(geometry);
            LineBatch edge{"band cut", {0.0f, 0.0f, 0.0f, 1.0f}, 4.0f, true, {}};
            const auto at = [&](double dist) {
                const double t = dist / len;
                return Point{c.a[0] + t * d[0], c.a[1] + t * d[1], c.a[2] + t * d[2]};
            };
            const auto put = [&](const Point& x, double e) {
                const Point q = lift(x, e);
                edge.points.push_back(Float3{static_cast<float>(q[0] * um_per_cm), static_cast<float>(q[1] * um_per_cm),
                                             static_cast<float>(q[2] * um_per_cm)});
            };
            for (const std::size_t k : shown) {
                std::optional<std::pair<Point, double>> last;
                for (std::size_t i = 1; i < c.samples.size(); ++i) {
                    if (c.samples[i].gap_before) {
                        last.reset();
                        continue;
                    }
                    constexpr int parts = 16;
                    for (int j = 0; j < parts; ++j) {
                        const double d0 = c.samples[i - 1].distance, d1 = c.samples[i].distance;
                        const Point p0 = at(d0 + (d1 - d0) * j / parts), p1 = at(d0 + (d1 - d0) * (j + 1) / parts);
                        const Point pm = at(d0 + (d1 - d0) * (j + 0.5) / parts);
                        const auto hit = locator.locate(pm, c.samples[i].element);
                        if (!hit) {
                            last.reset();
                            continue;
                        }
                        const NodeId owner = band_owner(geometry, hit->element, pm);
                        if (keys[static_cast<std::size_t>(owner)] < 0) {
                            last.reset();
                            continue;
                        }
                        const double e0 = band_value(geometry, hit->element, p0, owner, keys, values(k));
                        const double e1 = band_value(geometry, hit->element, p1, owner, keys, values(k));
                        if (!std::isfinite(e0) || !std::isfinite(e1) || e0 < lo || e0 > hi || e1 < lo || e1 > hi) {
                            last.reset();
                            continue;
                        }
                        if (last && last->second != e0) {  // a step: straight up or down
                            put(p0, last->second);
                            put(p0, e0);
                        }
                        put(p0, e0);
                        put(p1, e1);
                        last = std::pair{p1, e1};
                    }
                }
            }
            if (!edge.points.empty()) {
                LineBatch inner = edge;
                inner.color = {1.0f, 1.0f, 1.0f, 1.0f};
                inner.width_px = 2.0f;
                s.lines.push_back(std::move(edge));
                s.lines.push_back(std::move(inner));
            }
        }
    }
    if (!o.by_level) s.colorbar = Colorbar{o.map, *scale, "E (eV)"};
    s.axes = std::array<SceneAxis, 3>{SceneAxis{base.u_title, 0.0, 1.0}, SceneAxis{"E (eV)", 0.0, out.um_per_eV},
                                      SceneAxis{base.v_title, 0.0, 1.0}};
    ViewOptions view = o.view;
    view.clip.reset();
    detail::frame(s, view);
    s.ground_shadow = false;
    // BSD6: an axis-aligned plane's place in the device, for the corner inset.
    if (o.plane) {
        const Point& n = o.plane->normal;
        for (std::size_t a = 0; a < 3; ++a) {
            if (n[a] == 0.0 || n[(a + 1) % 3] != 0.0 || n[(a + 2) % 3] != 0.0) continue;
            PlaneInset inset;
            inset.axis = static_cast<int>(a);
            inset.at_um = o.plane->offset / n[a] * 1e4;
            inset.lo_um = {1e300, 1e300, 1e300};
            inset.hi_um = {-1e300, -1e300, -1e300};
            for (const auto& q : device.mesh().points()) {
                for (std::size_t k = 0; k < 3; ++k) {
                    inset.lo_um[k] = std::min(inset.lo_um[k], q[k] * 1e4);
                    inset.hi_um[k] = std::max(inset.hi_um[k], q[k] * 1e4);
                }
            }
            s.plane_inset = inset;
        }
    }
    return out;
}

std::optional<std::pair<double, double>> band_energy_range(const device::Device& device,
                                                           std::span<const results::BiasPoint> points,
                                                           const std::array<bool, 4>& levels) {
    const auto keys = band_material_keys(device);
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (const results::BiasPoint& p : points) {
        for (std::size_t k = 0; k < 4; ++k) {
            if (!levels[k]) continue;
            const auto& v = level_values(p.bands, k);
            if (v.size() != keys.size()) continue;
            for (std::size_t i = 0; i < v.size(); ++i) {
                if (keys[i] < 0 || !std::isfinite(v[i])) continue;
                lo = std::min(lo, v[i]);
                hi = std::max(hi, v[i]);
            }
        }
    }
    if (!(lo <= hi)) return std::nullopt;
    if (!(hi - lo > 1e-9)) {
        const double mid = 0.5 * (lo + hi);
        return std::pair{mid - 0.5, mid + 0.5};
    }
    return std::pair{lo, hi};
}

void add_plane_outline(FieldScene& scene, const Geometry& geometry, const Plane& plane) {
    // The section's edges keyed by their ends rounded to 1e-9 of the device's size: an edge two polygons share is
    // inside the section.
    const auto polygons = geometry.slice(plane);
    if (polygons.empty()) return;
    double size = 0.0;
    for (const Polygon& p : polygons) {
        for (const SurfaceVertex& sv : p.vertices) size = std::max({size, std::abs(sv.position[0]), std::abs(sv.position[1]), std::abs(sv.position[2])});
    }
    const double q = std::max(size, 1e-30) * 1e-9;
    using Key = std::array<std::int64_t, 3>;
    const auto key = [&](const Point& x) {
        return Key{std::llround(x[0] / q), std::llround(x[1] / q), std::llround(x[2] / q)};
    };
    std::map<std::pair<Key, Key>, std::pair<std::size_t, std::pair<Point, Point>>> edges;
    for (const Polygon& p : polygons) {
        for (std::size_t i = 0; i < p.vertices.size(); ++i) {
            const Point& a = p.vertices[i].position;
            const Point& b = p.vertices[(i + 1) % p.vertices.size()].position;
            Key ka = key(a), kb = key(b);
            if (ka == kb) continue;
            if (kb < ka) std::swap(ka, kb);
            auto [it, fresh] = edges.try_emplace({ka, kb}, 0, std::pair{a, b});
            ++it->second.first;
        }
    }
    LineBatch line{"band plane", {0.80f, 0.38f, 0.05f, 1.0f}, 2.5f, true, {}};
    for (const auto& [k, e] : edges) {
        if (e.first != 1) continue;
        for (const Point& x : {e.second.first, e.second.second}) {
            line.points.push_back(Float3{static_cast<float>(x[0] * 1e4), static_cast<float>(x[1] * 1e4),
                                         static_cast<float>(x[2] * 1e4)});
        }
    }
    if (!line.points.empty()) scene.lines.push_back(std::move(line));
}

}  // namespace NiTCAD::render
