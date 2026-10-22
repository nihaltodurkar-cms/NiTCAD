// Probe and measure (probe.hpp).
#include "NiTCAD/app/probe.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "NiTCAD/render/image.hpp"

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

}  // namespace NiTCAD::app
