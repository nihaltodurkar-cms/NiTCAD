#include "NiTCAD/render/cut_line.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <utility>

namespace NiTCAD::render {

namespace {

base::Error invalid(std::string message) {
    return base::Error{base::ErrorCode::invalid_input, std::move(message), std::nullopt};
}

Point sub(const Point& a, const Point& b) noexcept { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const Point& a, const Point& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Point cross(const Point& a, const Point& b) noexcept {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// A weight changing by at most this much over the whole line runs parallel to its face; a line on
// the face (weight 0 to within this) reads the element.
constexpr double parallel = 1e-10;
// Line parameters closer than this are the same point.
constexpr double same_t = 1e-12;
// Weights below this at a sample are zero (the sample is on the face away from the node).
constexpr double zero_weight = 1e-12;

// The element's barycentric weights at x, as affine functions extended beyond it (Cramer, as
// Locator::inside), or nothing for a degenerate element.
std::optional<std::array<double, 4>> barycentric(const Geometry& g, std::size_t e, const Point& x) {
    const auto n = g.element(e);
    const auto points = g.mesh().points();
    const int D = g.dimension();
    const Point& x0 = points[static_cast<std::size_t>(n[0])];
    std::array<Point, 3> col{};
    for (int k = 0; k < D; ++k) {
        col[static_cast<std::size_t>(k)] = sub(points[static_cast<std::size_t>(n[static_cast<std::size_t>(k + 1)])], x0);
    }
    const Point r = sub(x, x0);
    std::array<double, 4> w{};
    if (D == 1) {
        if (col[0][0] == 0.0) return std::nullopt;
        w[1] = r[0] / col[0][0];
    } else if (D == 2) {
        const double det = col[0][0] * col[1][1] - col[1][0] * col[0][1];
        if (det == 0.0) return std::nullopt;
        w[1] = (r[0] * col[1][1] - col[1][0] * r[1]) / det;
        w[2] = (col[0][0] * r[1] - r[0] * col[0][1]) / det;
    } else {
        const double det = dot(col[0], cross(col[1], col[2]));
        if (det == 0.0) return std::nullopt;
        w[1] = dot(r, cross(col[1], col[2])) / det;
        w[2] = dot(col[0], cross(r, col[2])) / det;
        w[3] = dot(col[0], cross(col[1], r)) / det;
    }
    w[0] = 1.0;
    for (int k = 1; k <= D; ++k) w[0] -= w[static_cast<std::size_t>(k)];
    return w;
}

// The part of the line in one element: t from t0 to t1, weights wa + t dw.
struct Piece {
    double t0 = 0.0;
    double t1 = 0.0;
    std::size_t element = 0;
    std::array<double, 4> wa{};
    std::array<double, 4> dw{};
};

CutSample sample_at(const Geometry& g, const Piece& p, double t, double length, bool gap) {
    CutSample s;
    s.distance = t * length;
    s.element = p.element;
    s.gap_before = gap;
    const auto n = g.element(p.element);
    double sum = 0.0;
    for (std::size_t k = 0; k < n.size(); ++k) {
        double w = p.wa[k] + t * p.dw[k];
        if (w < zero_weight) w = 0.0;
        s.at.nodes[k] = n[k];
        s.at.weights[k] = w;
        sum += w;
    }
    for (std::size_t k = 0; k < n.size(); ++k) s.at.weights[k] /= sum;
    return s;
}

}  // namespace

std::expected<CutLine, base::Error> cut_line(const Geometry& g, const Point& a, const Point& b) {
    if (g.element_count() == 0) return std::unexpected(invalid("a cut line needs a mesh with elements"));
    const auto D = static_cast<std::size_t>(g.dimension());
    CutLine c;
    for (std::size_t k = 0; k < D; ++k) {
        c.a[k] = a[k];
        c.b[k] = b[k];
    }
    const Point ab = sub(c.b, c.a);
    const double length = std::sqrt(dot(ab, ab));
    if (!std::isfinite(length) || !std::isfinite(dot(c.a, c.a))) {
        return std::unexpected(invalid("a cut line's end points must be finite"));
    }
    if (length == 0.0) return std::unexpected(invalid("a cut line needs two different end points"));

    // The line's bounding box, widened by a little more than rounding, rules out most elements.
    Point lo{}, hi{};
    for (std::size_t k = 0; k < D; ++k) {
        const double margin = 1e-9 * length;
        lo[k] = std::min(c.a[k], c.b[k]) - margin;
        hi[k] = std::max(c.a[k], c.b[k]) + margin;
    }
    const auto points = g.mesh().points();
    std::vector<Piece> pieces;
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        const auto n = g.element(e);
        bool apart = false;
        for (std::size_t k = 0; k < D && !apart; ++k) {
            double emin = points[static_cast<std::size_t>(n[0])][k], emax = emin;
            for (const NodeId v : n) {
                emin = std::min(emin, points[static_cast<std::size_t>(v)][k]);
                emax = std::max(emax, points[static_cast<std::size_t>(v)][k]);
            }
            apart = emax < lo[k] || emin > hi[k];
        }
        if (apart) continue;
        const auto wa = barycentric(g, e, c.a);
        const auto wb = barycentric(g, e, c.b);
        if (!wa || !wb) continue;
        Piece p{0.0, 1.0, e, *wa, {}};
        bool empty = false;
        for (std::size_t k = 0; k < n.size() && !empty; ++k) {
            p.dw[k] = (*wb)[k] - (*wa)[k];
            if (std::abs(p.dw[k]) <= parallel) {
                empty = p.wa[k] < -parallel;
            } else if (const double t = -p.wa[k] / p.dw[k]; p.dw[k] > 0.0) {
                p.t0 = std::max(p.t0, t);
            } else {
                p.t1 = std::min(p.t1, t);
            }
        }
        if (!empty && p.t1 - p.t0 > same_t) pieces.push_back(p);
    }
    std::ranges::sort(pieces, [](const Piece& x, const Piece& y) {
        return x.t0 != y.t0 ? x.t0 < y.t0 : x.t1 > y.t1;
    });

    // The union of the pieces, in order: an entry sample where the line enters the mesh (at the start
    // or after a gap), then each piece's exit; a piece within the part covered adds nothing, and one
    // overlapping it (a line along a shared face) continues from its end.
    double covered = 0.0;
    for (const Piece& p : pieces) {
        if (!c.samples.empty() && p.t1 <= covered + same_t) continue;
        const bool gap = !c.samples.empty() && p.t0 > covered + same_t;
        if (c.samples.empty() || gap) c.samples.push_back(sample_at(g, p, p.t0, length, gap));
        c.samples.push_back(sample_at(g, p, p.t1, length, false));
        covered = p.t1;
    }
    return c;
}

std::vector<double> cut_values(const CutLine& cut, std::span<const double> field) {
    std::vector<double> out;
    out.reserve(cut.samples.size());
    for (const CutSample& s : cut.samples) {
        double v = 0.0;
        for (std::size_t k = 0; k < s.at.nodes.size(); ++k) {
            if (s.at.weights[k] > 0.0) v += s.at.weights[k] * field[static_cast<std::size_t>(s.at.nodes[k])];
        }
        out.push_back(v);
    }
    return out;
}

}  // namespace NiTCAD::render
