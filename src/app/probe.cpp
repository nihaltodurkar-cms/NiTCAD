// Probe and measure (probe.hpp).
#include "NiTCAD/app/probe.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>

#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/render/band_surface.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/streamlines.hpp"

namespace NiTCAD::app {

Prober::Prober(const Workspace& w) {
    if (!w.device || !w.geometry) return;
    mesh_ = &w.device->mesh();
    dimension_ = mesh_->dimension();
    if (dimension_ == 2) {
        candidates_.resize(mesh_->node_count());
        for (std::size_t i = 0; i < candidates_.size(); ++i) candidates_[i] = static_cast<mesh::NodeId>(i);
    } else if (dimension_ == 3) {
        // The boundary faces' corners are nodes: each vertex's interpolant is one node at weight 1.
        std::vector<char> on(mesh_->node_count(), 0);
        for (const render::Polygon& f : w.geometry->boundary_faces()) {
            for (const render::SurfaceVertex& v : f.vertices) {
                for (std::size_t k = 0; k < v.at.nodes.size(); ++k) {
                    if (v.at.weights[k] > 0.0) on[static_cast<std::size_t>(v.at.nodes[k])] = 1;
                }
            }
        }
        for (std::size_t i = 0; i < on.size(); ++i) {
            if (on[i]) candidates_.push_back(static_cast<mesh::NodeId>(i));
        }
    }
    workspace_ = &w;
    all_ = candidates_;
}

void Prober::hide(std::span<const std::int32_t> regions) {
    candidates_ = all_;
    if (regions.empty() || !workspace_) return;
    const std::vector<std::int32_t> region(workspace_->device->node_region().begin(), workspace_->device->node_region().end());
    const auto hidden = [&](mesh::NodeId n) {
        return std::ranges::find(regions, region[static_cast<std::size_t>(n)]) != regions.end();
    };
    std::erase_if(candidates_, hidden);
    if (dimension_ != 3) return;
    // The faces the hidden regions opened: their shown nodes.
    std::vector<char> on(mesh_->node_count(), 0);
    for (const mesh::NodeId n : candidates_) on[static_cast<std::size_t>(n)] = 1;
    for (const render::InterfacePiece& p : workspace_->geometry->interfaces(region)) {
        const bool ha = hidden(p.node_a), hb = hidden(p.node_b);
        if (ha != hb) on[static_cast<std::size_t>(ha ? p.node_b : p.node_a)] = 1;
    }
    candidates_.clear();
    for (std::size_t i = 0; i < on.size(); ++i) {
        if (on[i]) candidates_.push_back(static_cast<mesh::NodeId>(i));
    }
}

std::optional<ProbeHit> Prober::pick(const render::FieldScene& scene, int width, int height, float x, float y,
                                     std::span<const double> values, float radius_px) const {
    std::optional<ProbeHit> best;
    float best_depth = std::numeric_limits<float>::infinity();
    for (const mesh::NodeId n : candidates_) {
        const auto& p = mesh_->points()[static_cast<std::size_t>(n)];
        const render::Float3 um{static_cast<float>(p[0] * 1e4), static_cast<float>(p[1] * 1e4),
                                dimension_ == 3 ? static_cast<float>(p[2] * 1e4) : 0.0f};
        const render::ScreenPoint s = render::screen_position(scene, width, height, um);
        const float d = std::hypot(s.x - x, s.y - y);
        if (!(d <= radius_px)) continue;
        // 2D: the nearest; 3D: the front-most, the nearer of equal depths.
        const bool better = !best || (dimension_ == 3 ? (s.depth < best_depth - 1e-6f ||
                                                         (std::abs(s.depth - best_depth) <= 1e-6f && d < best->distance_px))
                                                      : d < best->distance_px);
        if (!better) continue;
        ProbeHit h;
        h.node = n;
        h.position_um = {p[0] * 1e4, p[1] * 1e4, dimension_ == 3 ? p[2] * 1e4 : 0.0};
        h.value = static_cast<std::size_t>(n) < values.size() ? values[static_cast<std::size_t>(n)]
                                                               : std::numeric_limits<double>::quiet_NaN();
        h.distance_px = d;
        best = h;
        best_depth = s.depth;
    }
    return best;
}

std::optional<std::array<double, 3>> Prober::point_at(const render::FieldScene& scene, int width, int height, float x,
                                                      float y, float radius_px) const {
    if (dimension_ == 3) {
        const auto hit = pick(scene, width, height, x, y, {}, radius_px);
        if (!hit) return std::nullopt;
        return hit->position_um;
    }
    if (dimension_ != 2) return std::nullopt;
    // The view of the plane is affine: invert it from where the origin and the unit points land.
    const render::ScreenPoint o = render::screen_position(scene, width, height, {0.0f, 0.0f, 0.0f});
    const render::ScreenPoint ex = render::screen_position(scene, width, height, {1.0f, 0.0f, 0.0f});
    const render::ScreenPoint ey = render::screen_position(scene, width, height, {0.0f, 1.0f, 0.0f});
    const double a = ex.x - o.x, b = ey.x - o.x, c = ex.y - o.y, d = ey.y - o.y;
    const double det = a * d - b * c;
    if (det == 0.0) return std::nullopt;
    const double rx = x - o.x, ry = y - o.y;
    return std::array<double, 3>{(rx * d - b * ry) / det, (a * ry - rx * c) / det, 0.0};
}

double distance_um(const ProbeHit& a, const ProbeHit& b) noexcept {
    const double dx = a.position_um[0] - b.position_um[0], dy = a.position_um[1] - b.position_um[1],
                 dz = a.position_um[2] - b.position_um[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::string hover_text(std::string_view title, const ProbeHit& hit, int dimension) {
    const auto [symbol, unit] = title_parts(title);
    std::string text(symbol);
    if (std::isfinite(hit.value)) {
        text += std::format(" = {:.4g}", hit.value);
        if (!unit.empty()) text += " " + std::string(unit);
    } else {
        text += ": no value";
    }
    text += std::format("\nx {:.3g}  y {:.3g}", hit.position_um[0], hit.position_um[1]);
    if (dimension == 3) text += std::format("  z {:.3g}", hit.position_um[2]);
    text += std::format(" µm\nnode {}", hit.node);
    return text;
}

std::optional<CurveHit> pick_curve(const render::PlotScene& scene, int width, int height, float x, float y,
                                   float radius_px) {
    std::optional<CurveHit> best;
    for (std::size_t k = 0; k < scene.series.size(); ++k) {
        const render::Series& s = scene.series[k];
        const std::size_t n = std::min(s.x.size(), s.y.size());
        for (std::size_t i = 0; i < n; ++i) {
            const render::ScreenPoint p = render::plot_position(scene, width, height, s.x[i], s.y[i]);
            if (!std::isfinite(p.x) || !std::isfinite(p.y)) continue;
            const float d = std::hypot(p.x - x, p.y - y);
            if (!(d <= radius_px) || (best && !(d < best->distance_px))) continue;
            best = CurveHit{k, i, s.x[i], s.y[i], p.x, p.y, d};
        }
    }
    return best;
}

std::string curve_text(const render::PlotScene& scene, const CurveHit& hit) {
    std::string text;
    if (scene.series.size() > 1 && hit.series < scene.series.size()) text = scene.series[hit.series].name + "\n";
    text += std::format("{} = {:.4g}\n{} = {:.4g}", scene.x.label, hit.x, scene.y.label, hit.y);
    return text;
}

// --- band surface readout (Unit 26o, BSD7) ---------------------------------------------------------------

std::optional<BandHit> band_pick(const Workspace& w, const ViewSettings& s, const render::FieldScene& scene,
                                 int width, int height, float x, float y) {
    if (s.kind != ViewKind::band_surface || !scene.axes || !w.device || !w.geometry || w.points.empty()) {
        return std::nullopt;
    }
    const std::size_t point = s.point.value_or(w.points.size() - 1);
    if (point >= w.points.size()) return std::nullopt;
    const int dimension = w.device->mesh().dimension();
    if (dimension == 3 && !s.band_plane) return std::nullopt;
    constexpr std::array<std::string_view, 4> names{"Ec", "Ev", "EFn", "EFp"};
    // The nearest triangle of a shown level the ray meets (Moller-Trumbore, in double).
    const render::ScreenRay ray = render::screen_ray(scene, width, height, x, y);
    const auto sub = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };
    const auto cross = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return std::array<double, 3>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    };
    const auto dot = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    double best = std::numeric_limits<double>::infinity();
    int level = -1;
    for (const render::TriangleBatch& b : scene.triangles) {
        const auto n = std::ranges::find(names, b.name);
        if (n == names.end()) continue;
        const auto k = static_cast<std::size_t>(n - names.begin());
        if (!s.band_levels[k]) continue;
        for (std::size_t i = 0; i + 2 < b.vertices.size(); i += 3) {
            std::array<std::array<double, 3>, 3> v{};
            for (std::size_t j = 0; j < 3; ++j) {
                for (std::size_t c = 0; c < 3; ++c) v[j][c] = b.vertices[i + j].position[c];
            }
            const auto e1 = sub(v[1], v[0]), e2 = sub(v[2], v[0]);
            const auto p = cross(ray.direction, e2);
            const double det = dot(e1, p);
            if (std::abs(det) < 1e-300) continue;
            const auto t0 = sub(ray.origin, v[0]);
            const double u = dot(t0, p) / det;
            if (u < 0.0 || u > 1.0) continue;
            const auto q = cross(t0, e1);
            const double w2 = dot(ray.direction, q) / det;
            if (w2 < 0.0 || u + w2 > 1.0) continue;
            const double t = dot(e2, q) / det;
            if (t > 0.0 && t < best) {
                best = t;
                level = static_cast<int>(k);
            }
        }
    }
    if (level < 0) return std::nullopt;
    // The device point under the hit: scene x and z are the base's u and v (BSD1).
    const std::array<double, 3> at{ray.origin[0] + best * ray.direction[0], ray.origin[1] + best * ray.direction[1],
                                   ray.origin[2] + best * ray.direction[2]};
    BandHit hit;
    hit.level = level;
    if (dimension == 2) {
        hit.position_um = {at[0], at[2], 0.0};
    } else {
        const auto a = static_cast<std::size_t>(s.band_plane->axis);
        const std::size_t b0 = a == 0 ? 1 : 0, b1 = a == 2 ? 1 : 2;
        hit.position_um[a] = s.band_plane->at_um;
        hit.position_um[b0] = at[0];
        hit.position_um[b1] = at[2];
    }
    const render::Point x_cm{hit.position_um[0] * 1e-4, hit.position_um[1] * 1e-4, hit.position_um[2] * 1e-4};
    const render::Locator locator(*w.geometry);
    const auto located = locator.locate(x_cm);
    if (!located) return std::nullopt;
    const auto& bands = w.points[point].bands;
    for (std::size_t k = 0; k < 4; ++k) {
        if (!s.band_levels[k]) continue;
        hit.energy_eV[k] = render::band_energy_at(*w.device, *w.geometry, bands, static_cast<render::BandLevel>(k),
                                                  located->element, x_cm);
    }
    const mesh::NodeId owner = render::band_owner(*w.geometry, located->element, x_cm);
    const auto regions = w.device->regions();
    const auto r = w.device->node_region()[static_cast<std::size_t>(owner)];
    if (static_cast<std::size_t>(r) < regions.size()) hit.region = regions[static_cast<std::size_t>(r)].name;
    return hit;
}

std::string band_hover_text(const BandHit& hit, int dimension) {
    constexpr std::array<std::string_view, 4> names{"Ec", "Ev", "EFn", "EFp"};
    std::string text;
    const auto line = [&](std::size_t k) {
        if (!text.empty()) text += '\n';
        text += hit.energy_eV[k] ? std::format("{} = {:.4g} eV", names[k], *hit.energy_eV[k])
                                 : std::format("{}: no value", names[k]);
    };
    const auto first = static_cast<std::size_t>(std::clamp(hit.level, 0, 3));
    line(first);
    for (std::size_t k = 0; k < 4; ++k) {
        if (k != first && hit.energy_eV[k]) line(k);
    }
    text += std::format("\nx {:.3g}  y {:.3g}", hit.position_um[0], hit.position_um[1]);
    if (dimension == 3) text += std::format("  z {:.3g}", hit.position_um[2]);
    text += " \xC2\xB5m";
    if (!hit.region.empty()) text += "\n" + hit.region;
    return text;
}

}  // namespace NiTCAD::app
