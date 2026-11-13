// Point location for state transfer (Unit 17, slice 3; DECISIONS.md Unit 17 AD5).
#include "NiTCAD/mesh/locate.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "predicates.hpp"

namespace NiTCAD::mesh {
namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}

double det3(const Point& a, const Point& b, const Point& c) {
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}

Point sub(const Point& a, const Point& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

}  // namespace

std::expected<PointLocator, base::Error> PointLocator::create(const SimplexMesh& sm) {
    if (sm.dimension != 2 && sm.dimension != 3) return std::unexpected(invalid("a locator needs a 2D or 3D simplex mesh"));
    const std::size_t k = static_cast<std::size_t>(sm.dimension) + 1;
    if (sm.cells.empty() || sm.cells.size() % k != 0 || sm.points.empty()) {
        return std::unexpected(invalid("a locator needs points and whole cells"));
    }
    if (sm.cells.size() / k > std::numeric_limits<std::uint32_t>::max() / 8) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "a mesh too large to locate in", std::nullopt});
    }
    PointLocator L;
    L.dim_ = sm.dimension;
    L.points_ = sm.points;
    L.nodes_ = sm.points.size();
    L.simplices_ = sm.cells;
    const std::size_t nc = sm.cells.size() / k;
    for (std::size_t c = 0; c < nc; ++c) {
        NodeId* v = &L.simplices_[k * c];
        for (std::size_t i = 0; i < k; ++i) {
            if (v[i] < 0 || static_cast<std::size_t>(v[i]) >= sm.points.size()) {
                return std::unexpected(invalid("a cell names a point out of range", c));
            }
        }
        const auto P = [&](std::size_t i) -> const Point& { return L.points_[static_cast<std::size_t>(v[i])]; };
        const int o = L.dim_ == 2 ? predicates::orient2d(P(0), P(1), P(2)) : predicates::orient3d(P(0), P(1), P(2), P(3));
        if (o == 0) return std::unexpected(invalid("a flat cell", c));
        if (o < 0) std::swap(v[1], v[2]);
    }
    Point lo = L.points_[0], hi = L.points_[0];
    for (const Point& p : L.points_) {
        for (int a = 0; a < 3; ++a) {
            lo[static_cast<std::size_t>(a)] = std::min(lo[static_cast<std::size_t>(a)], p[static_cast<std::size_t>(a)]);
            hi[static_cast<std::size_t>(a)] = std::max(hi[static_cast<std::size_t>(a)], p[static_cast<std::size_t>(a)]);
        }
    }
    double extent = std::numeric_limits<double>::min();
    for (int a = 0; a < L.dim_; ++a) extent = std::max(extent, hi[static_cast<std::size_t>(a)] - lo[static_cast<std::size_t>(a)]);
    const double per_axis = std::max(1.0, std::pow(static_cast<double>(nc) / 2.0, 1.0 / L.dim_));
    L.cell_ = extent / per_axis;
    L.lo_ = lo;
    for (int a = 0; a < L.dim_; ++a) {
        L.n_[static_cast<std::size_t>(a)] =
            std::min<std::size_t>(512, static_cast<std::size_t>((hi[static_cast<std::size_t>(a)] - lo[static_cast<std::size_t>(a)]) / L.cell_) + 1);
    }
    const auto bucket = [&](int a, double x) {
        const double c = std::floor((x - L.lo_[static_cast<std::size_t>(a)]) / L.cell_);
        return static_cast<std::size_t>(std::clamp(c, 0.0, static_cast<double>(L.n_[static_cast<std::size_t>(a)] - 1)));
    };
    const auto for_buckets = [&](std::size_t c, auto&& fn) {
        std::array<std::size_t, 3> b0{0, 0, 0}, b1{0, 0, 0};
        for (int a = 0; a < L.dim_; ++a) {
            double x0 = std::numeric_limits<double>::infinity(), x1 = -x0;
            for (std::size_t i = 0; i < k; ++i) {
                const double x = L.points_[static_cast<std::size_t>(L.simplices_[k * c + i])][static_cast<std::size_t>(a)];
                x0 = std::min(x0, x);
                x1 = std::max(x1, x);
            }
            b0[static_cast<std::size_t>(a)] = bucket(a, x0);
            b1[static_cast<std::size_t>(a)] = bucket(a, x1);
        }
        for (std::size_t z = b0[2]; z <= b1[2]; ++z) {
            for (std::size_t y = b0[1]; y <= b1[1]; ++y) {
                for (std::size_t x = b0[0]; x <= b1[0]; ++x) fn((z * L.n_[1] + y) * L.n_[0] + x);
            }
        }
    };
    const std::size_t nb = L.n_[0] * L.n_[1] * L.n_[2];
    std::vector<std::uint32_t> count(nb + 1, 0);
    for (std::size_t c = 0; c < nc; ++c) for_buckets(c, [&](std::size_t b) { ++count[b + 1]; });
    for (std::size_t b = 1; b < count.size(); ++b) count[b] += count[b - 1];
    L.start_ = count;
    L.items_.assign(count.back(), 0);
    std::vector<std::uint32_t> fill(count.begin(), count.end() - 1);
    for (std::size_t c = 0; c < nc; ++c) {
        for_buckets(c, [&](std::size_t b) { L.items_[fill[b]++] = static_cast<std::uint32_t>(c); });
    }
    return L;
}

PointLocator PointLocator::create(const TensorCells& cells) {
    PointLocator L;
    L.cells_ = cells;
    L.dim_ = cells.dimension();
    L.nodes_ = 1;
    for (int a = 0; a < L.dim_; ++a) L.nodes_ *= cells.size(a);
    return L;
}

bool PointLocator::contains(std::size_t c, const Point& p) const {
    const std::size_t k = static_cast<std::size_t>(dim_) + 1;
    const auto P = [&](std::size_t i) -> const Point& { return points_[static_cast<std::size_t>(simplices_[k * c + i])]; };
    if (dim_ == 2) {
        return predicates::orient2d(P(0), P(1), p) >= 0 && predicates::orient2d(P(1), P(2), p) >= 0 &&
               predicates::orient2d(P(2), P(0), p) >= 0;
    }
    // A positively oriented tetrahedron (a, b, c, d): p replacing each vertex keeps the sign.
    return predicates::orient3d(p, P(1), P(2), P(3)) >= 0 && predicates::orient3d(P(0), p, P(2), P(3)) >= 0 &&
           predicates::orient3d(P(0), P(1), p, P(3)) >= 0 && predicates::orient3d(P(0), P(1), P(2), p) >= 0;
}

Interpolant PointLocator::weights(std::size_t c, const Point& p) const {
    const std::size_t k = static_cast<std::size_t>(dim_) + 1;
    Interpolant r;
    r.count = static_cast<int>(k);
    for (std::size_t i = 0; i < k; ++i) r.nodes[i] = simplices_[k * c + i];
    const auto P = [&](std::size_t i) -> const Point& { return points_[static_cast<std::size_t>(simplices_[k * c + i])]; };
    if (dim_ == 2) {
        const Point& a = P(0);
        const Point& b = P(1);
        const Point& d = P(2);
        const double det = (b[0] - a[0]) * (d[1] - a[1]) - (d[0] - a[0]) * (b[1] - a[1]);
        double lb = ((p[0] - a[0]) * (d[1] - a[1]) - (d[0] - a[0]) * (p[1] - a[1])) / det;
        double lc = ((b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1])) / det;
        lb = std::clamp(lb, 0.0, 1.0);
        lc = std::clamp(lc, 0.0, 1.0 - lb);
        r.weights[0] = 1.0 - lb - lc;
        r.weights[1] = lb;
        r.weights[2] = lc;
        return r;
    }
    const Point e1 = sub(P(1), P(0)), e2 = sub(P(2), P(0)), e3 = sub(P(3), P(0)), q = sub(p, P(0));
    const double det = det3(e1, e2, e3);
    double l1 = std::clamp(det3(q, e2, e3) / det, 0.0, 1.0);
    double l2 = std::clamp(det3(e1, q, e3) / det, 0.0, 1.0 - l1);
    double l3 = std::clamp(det3(e1, e2, q) / det, 0.0, 1.0 - l1 - l2);
    r.weights[0] = 1.0 - l1 - l2 - l3;
    r.weights[1] = l1;
    r.weights[2] = l2;
    r.weights[3] = l3;
    return r;
}

Interpolant PointLocator::nearest(const Point& p) const {
    std::size_t best = 0;
    double bd = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < points_.size(); ++i) {
        const double dx = points_[i][0] - p[0], dy = points_[i][1] - p[1], dz = points_[i][2] - p[2];
        const double d = dx * dx + dy * dy + dz * dz;
        if (d < bd) {
            bd = d;
            best = i;
        }
    }
    Interpolant r;
    r.count = 1;
    r.nodes[0] = static_cast<NodeId>(best);
    r.weights[0] = 1.0;
    r.inside = false;
    return r;
}

Interpolant PointLocator::locate(const Point& p) const {
    if (cells_) {
        const TensorCells& tc = *cells_;
        TensorCells::Index cell{0, 0, 0};
        for (int a = 0; a < tc.dimension(); ++a) {
            const auto ax = tc.axis(a);
            const auto it = std::upper_bound(ax.begin(), ax.end(), p[static_cast<std::size_t>(a)]);
            const auto i = static_cast<std::size_t>(std::distance(ax.begin(), it));
            cell[static_cast<std::size_t>(a)] = std::clamp<std::size_t>(i == 0 ? 0 : i - 1, 0, ax.size() - 2);
        }
        const auto s = tc.stencil(cell, p);
        Interpolant r;
        r.count = s.count;
        r.nodes = s.nodes;
        r.weights = s.weights;
        for (int a = 0; a < tc.dimension(); ++a) {
            const auto ax = tc.axis(a);
            const double x = p[static_cast<std::size_t>(a)];
            if (x < ax.front() || x > ax.back()) r.inside = false;
        }
        return r;
    }
    std::array<std::size_t, 3> b{0, 0, 0};
    for (int a = 0; a < dim_; ++a) {
        const double c = std::floor((p[static_cast<std::size_t>(a)] - lo_[static_cast<std::size_t>(a)]) / cell_);
        if (!(c >= 0.0 && c < static_cast<double>(n_[static_cast<std::size_t>(a)]))) return nearest(p);
        b[static_cast<std::size_t>(a)] = static_cast<std::size_t>(c);
    }
    const std::size_t bk = (b[2] * n_[1] + b[1]) * n_[0] + b[0];
    for (std::uint32_t i = start_[bk]; i < start_[bk + 1]; ++i) {
        if (contains(items_[i], p)) return weights(items_[i], p);
    }
    return nearest(p);
}

}  // namespace NiTCAD::mesh
