// Mesh quality (quality.hpp).
#include "NiTCAD/render/quality.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <string>
#include <utility>

namespace NiTCAD::render {

namespace {

Point sub(const Point& a, const Point& b) noexcept { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const Point& a, const Point& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Point cross(const Point& a, const Point& b) noexcept {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Point& a) noexcept { return std::sqrt(dot(a, a)); }

// The angle between two vectors, radians, from its sine and cosine parts.
double angle(const Point& u, const Point& v) noexcept { return std::atan2(norm(cross(u, v)), dot(u, v)); }

constexpr double degrees = 180.0 / std::numbers::pi;

double triangle(const Point& a, const Point& b, const Point& c) noexcept {
    const Point ab = sub(b, a), ac = sub(c, a), bc = sub(c, b);
    const double size = std::max({norm(ab), norm(ac), norm(bc)});
    // Coincident or collinear points.
    if (!(size > 0.0) || norm(cross(ab, ac)) <= 1e-14 * size * size) return 0.0;
    const double angles[3] = {angle(ab, ac), angle(sub(a, b), bc), angle(sub(a, c), sub(b, c))};
    return std::min({angles[0], angles[1], angles[2]}) * degrees;
}

double tetrahedron(std::span<const Point> p) noexcept {
    double size = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = i + 1; j < 4; ++j) size = std::max(size, norm(sub(p[j], p[i])));
    }
    const double six_volume = dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0])));
    if (!(size > 0.0) || std::abs(six_volume) <= 1e-14 * size * size * size) return 0.0;
    // The outward normal of the face opposite each vertex.
    std::array<Point, 4> n{};
    for (std::size_t k = 0; k < 4; ++k) {
        std::array<std::size_t, 3> f{};
        for (std::size_t i = 0, j = 0; i < 4; ++i) {
            if (i != k) f[j++] = i;
        }
        n[k] = cross(sub(p[f[1]], p[f[0]]), sub(p[f[2]], p[f[0]]));
        if (dot(n[k], sub(p[k], p[f[0]])) > 0.0) n[k] = {-n[k][0], -n[k][1], -n[k][2]};
    }
    // Every two faces meet along one edge; the dihedral angle there is pi less the normals' angle.
    double smallest = std::numbers::pi;
    for (std::size_t k = 0; k < 4; ++k) {
        for (std::size_t l = k + 1; l < 4; ++l) smallest = std::min(smallest, std::numbers::pi - angle(n[k], n[l]));
    }
    return smallest * degrees;
}

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

}  // namespace

double min_angle_deg(std::span<const Point> v) noexcept {
    if (v.size() == 3) return triangle(v[0], v[1], v[2]);
    if (v.size() == 4) return tetrahedron(v);
    return 0.0;
}

std::expected<std::vector<double>, base::Error> element_quality(const Geometry& g) {
    if (g.dimension() < 2) return std::unexpected(invalid("mesh quality needs a 2D or 3D mesh"));
    if (g.element_count() == 0) return std::unexpected(invalid("mesh quality needs a mesh with elements"));
    const auto points = g.mesh().points();
    std::vector<double> out;
    out.reserve(g.element_count());
    std::vector<Point> v;
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        v.clear();
        for (const NodeId n : g.element(e)) v.push_back(points[static_cast<std::size_t>(n)]);
        out.push_back(min_angle_deg(v));
    }
    return out;
}

double quality_threshold_deg(int dimension) noexcept { return dimension == 3 ? 10.0 : 20.0; }

double best_quality_deg(int dimension) noexcept { return dimension == 3 ? std::acos(1.0 / 3.0) * degrees : 60.0; }

QualitySummary quality_summary(std::span<const double> angles, int dimension) {
    QualitySummary s;
    s.threshold_deg = quality_threshold_deg(dimension);
    s.elements = angles.size();
    if (angles.empty()) return s;
    std::vector<double> sorted(angles.begin(), angles.end());
    std::ranges::sort(sorted);
    s.min_deg = sorted.front();
    const std::size_t h = sorted.size() / 2;
    s.median_deg = sorted.size() % 2 == 1 ? sorted[h] : 0.5 * (sorted[h - 1] + sorted[h]);
    s.below = static_cast<std::size_t>(std::ranges::lower_bound(sorted, s.threshold_deg) - sorted.begin());
    return s;
}

QualityHistogram quality_histogram(std::span<const double> angles) noexcept {
    QualityHistogram h;
    const std::size_t last = h.counts.size() - 1;
    for (const double a : angles) {
        std::size_t k = 0;  // a NaN (no element shape) in the first bin
        if (a >= 90.0) k = last;
        else if (a > 0.0) k = std::min(static_cast<std::size_t>(a / QualityHistogram::bin_deg), last);
        ++h.counts[k];
    }
    return h;
}

}  // namespace NiTCAD::render
