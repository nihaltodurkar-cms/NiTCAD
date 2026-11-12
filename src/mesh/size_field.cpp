// Size fields for adaptive refinement (Unit 17, slice 1; DECISIONS.md Unit 17 AD4, AD10).
#include "NiTCAD/mesh/size_field.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <utility>

#include "predicates.hpp"

namespace NiTCAD::mesh {
namespace {

using base::ErrorCode;

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}};
}

// The point of segment [a, b] nearest p, as a parameter in [0, 1].
double segment_parameter(const Point& p, const Point& a, const Point& b) {
    const double dx = b[0] - a[0], dy = b[1] - a[1];
    const double len2 = dx * dx + dy * dy;
    if (len2 == 0.0) return 0.0;
    return std::clamp(((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2, 0.0, 1.0);
}

double distance2(const Point& p, const Point& q) {
    const double dx = q[0] - p[0], dy = q[1] - p[1];
    return dx * dx + dy * dy;
}

}  // namespace

std::expected<std::shared_ptr<const SizeField>, base::Error> SizeField::create(std::vector<Point> points,
                                                                             std::vector<NodeId> triangles,
                                                                             std::vector<double> sizes) {
    if (points.empty()) return std::unexpected(invalid("a size field needs points"));
    if (triangles.empty() || triangles.size() % 3 != 0) {
        return std::unexpected(invalid("a size field needs triangles, three node ids each"));
    }
    if (sizes.size() != points.size()) return std::unexpected(invalid("a size field needs one size per point"));
    if (points.size() > std::numeric_limits<std::uint32_t>::max() / 2 ||
        triangles.size() / 3 > std::numeric_limits<std::uint32_t>::max() / 4) {
        return std::unexpected(base::Error{ErrorCode::resource_exhausted, "a size field is too large", std::nullopt});
    }
    for (std::size_t i = 0; i < points.size(); ++i) {
        const Point& p = points[i];
        if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || p[2] != 0.0) {
            return std::unexpected(invalid("a size field point must be finite with z = 0", i));
        }
        if (!(std::isfinite(sizes[i]) && sizes[i] > 0.0)) {
            return std::unexpected(invalid("a size field size must be finite and positive", i));
        }
    }
    const std::size_t nt = triangles.size() / 3;
    for (std::size_t t = 0; t < nt; ++t) {
        NodeId* v = &triangles[3 * t];
        for (int k = 0; k < 3; ++k) {
            if (v[k] < 0 || static_cast<std::size_t>(v[k]) >= points.size()) {
                return std::unexpected(invalid("a size field triangle names a point out of range", t));
            }
        }
        if (v[0] == v[1] || v[1] == v[2] || v[0] == v[2]) {
            return std::unexpected(invalid("a size field triangle repeats a point", t));
        }
        const int o = predicates::orient2d(points[static_cast<std::size_t>(v[0])], points[static_cast<std::size_t>(v[1])],
                                           points[static_cast<std::size_t>(v[2])]);
        if (o == 0) return std::unexpected(invalid("a size field triangle is flat", t));
        if (o < 0) std::swap(v[1], v[2]);
    }

    std::shared_ptr<SizeField> f(new SizeField());
    f->points_ = std::move(points);
    f->triangles_ = std::move(triangles);
    f->sizes_ = std::move(sizes);

    Point lo = f->points_[0], hi = f->points_[0];
    for (const Point& p : f->points_) {
        lo = {std::min(lo[0], p[0]), std::min(lo[1], p[1]), 0.0};
        hi = {std::max(hi[0], p[0]), std::max(hi[1], p[1]), 0.0};
    }
    const double wx = hi[0] - lo[0], wy = hi[1] - lo[1];
    const double extent = std::max({wx, wy, std::numeric_limits<double>::min()});
    // About two triangles per bucket.
    const double per_axis = std::max(1.0, std::sqrt(static_cast<double>(nt) / 2.0));
    f->cell_ = extent / per_axis;
    f->lo_ = lo;
    f->nx_ = std::min<std::size_t>(4096, static_cast<std::size_t>(wx / f->cell_) + 1);
    f->ny_ = std::min<std::size_t>(4096, static_cast<std::size_t>(wy / f->cell_) + 1);
    const auto bucket_x = [&](double x) {
        const double c = std::floor((x - f->lo_[0]) / f->cell_);
        return static_cast<std::size_t>(std::clamp(c, 0.0, static_cast<double>(f->nx_ - 1)));
    };
    const auto bucket_y = [&](double y) {
        const double c = std::floor((y - f->lo_[1]) / f->cell_);
        return static_cast<std::size_t>(std::clamp(c, 0.0, static_cast<double>(f->ny_ - 1)));
    };
    std::vector<std::uint32_t> count(f->nx_ * f->ny_ + 1, 0);
    const auto for_buckets = [&](std::size_t t, auto&& fn) {
        const NodeId* v = &f->triangles_[3 * t];
        double x0 = f->points_[static_cast<std::size_t>(v[0])][0], x1 = x0;
        double y0 = f->points_[static_cast<std::size_t>(v[0])][1], y1 = y0;
        for (int k = 1; k < 3; ++k) {
            const Point& p = f->points_[static_cast<std::size_t>(v[k])];
            x0 = std::min(x0, p[0]);
            x1 = std::max(x1, p[0]);
            y0 = std::min(y0, p[1]);
            y1 = std::max(y1, p[1]);
        }
        for (std::size_t by = bucket_y(y0); by <= bucket_y(y1); ++by) {
            for (std::size_t bx = bucket_x(x0); bx <= bucket_x(x1); ++bx) fn(by * f->nx_ + bx);
        }
    };
    for (std::size_t t = 0; t < nt; ++t) for_buckets(t, [&](std::size_t b) { ++count[b + 1]; });
    for (std::size_t b = 1; b < count.size(); ++b) count[b] += count[b - 1];
    f->bucket_start_ = count;
    f->bucket_items_.assign(count.back(), 0);
    std::vector<std::uint32_t> fill(count.begin(), count.end() - 1);
    for (std::size_t t = 0; t < nt; ++t) {
        for_buckets(t, [&](std::size_t b) { f->bucket_items_[fill[b]++] = static_cast<std::uint32_t>(t); });
    }
    return std::shared_ptr<const SizeField>(std::move(f));
}

double SizeField::interpolate(std::size_t t, const Point& p) const {
    const Point& a = points_[static_cast<std::size_t>(triangles_[3 * t])];
    const Point& b = points_[static_cast<std::size_t>(triangles_[3 * t + 1])];
    const Point& c = points_[static_cast<std::size_t>(triangles_[3 * t + 2])];
    const double det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
    double lb = ((p[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (p[1] - a[1])) / det;
    double lc = ((b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1])) / det;
    lb = std::clamp(lb, 0.0, 1.0);
    lc = std::clamp(lc, 0.0, 1.0 - lb);
    const double la = 1.0 - lb - lc;
    const double h = la * sizes_[static_cast<std::size_t>(triangles_[3 * t])] +
                     lb * sizes_[static_cast<std::size_t>(triangles_[3 * t + 1])] +
                     lc * sizes_[static_cast<std::size_t>(triangles_[3 * t + 2])];
    // Rounding cannot take a convex combination outside its corners' range.
    const double lo = std::min({sizes_[static_cast<std::size_t>(triangles_[3 * t])],
                                sizes_[static_cast<std::size_t>(triangles_[3 * t + 1])],
                                sizes_[static_cast<std::size_t>(triangles_[3 * t + 2])]});
    const double hi = std::max({sizes_[static_cast<std::size_t>(triangles_[3 * t])],
                                sizes_[static_cast<std::size_t>(triangles_[3 * t + 1])],
                                sizes_[static_cast<std::size_t>(triangles_[3 * t + 2])]});
    return std::clamp(h, lo, hi);
}

double SizeField::at(const Point& p) const {
    const double cx = std::floor((p[0] - lo_[0]) / cell_);
    const double cy = std::floor((p[1] - lo_[1]) / cell_);
    const bool inside_grid = cx >= 0.0 && cy >= 0.0 && cx < static_cast<double>(nx_) && cy < static_cast<double>(ny_);
    if (inside_grid) {
        const std::size_t b = static_cast<std::size_t>(cy) * nx_ + static_cast<std::size_t>(cx);
        for (std::uint32_t k = bucket_start_[b]; k < bucket_start_[b + 1]; ++k) {
            const std::size_t t = bucket_items_[k];
            const Point& a = points_[static_cast<std::size_t>(triangles_[3 * t])];
            const Point& q = points_[static_cast<std::size_t>(triangles_[3 * t + 1])];
            const Point& c = points_[static_cast<std::size_t>(triangles_[3 * t + 2])];
            if (predicates::orient2d(a, q, p) >= 0 && predicates::orient2d(q, c, p) >= 0 &&
                predicates::orient2d(c, a, p) >= 0) {
                return interpolate(t, p);
            }
        }
    }
    // Outside every triangle: the nearest triangle (lowest index on ties).
    std::size_t best = 0;
    Point best_q = p;
    double best_d = std::numeric_limits<double>::infinity();
    const std::size_t nt = triangles_.size() / 3;
    for (std::size_t t = 0; t < nt; ++t) {
        for (int k = 0; k < 3; ++k) {
            const Point& a = points_[static_cast<std::size_t>(triangles_[3 * t + static_cast<std::size_t>(k)])];
            const Point& b = points_[static_cast<std::size_t>(triangles_[3 * t + static_cast<std::size_t>((k + 1) % 3)])];
            const double s = segment_parameter(p, a, b);
            const Point q{a[0] + s * (b[0] - a[0]), a[1] + s * (b[1] - a[1]), 0.0};
            const double d = distance2(p, q);
            if (d < best_d) {
                best_d = d;
                best = t;
                best_q = q;
            }
        }
    }
    return interpolate(best, best_q);
}

double SizeField::min_size() const { return *std::ranges::min_element(sizes_); }

std::uint64_t SizeField::digest() const {
    std::uint64_t h = 0xcbf29ce484222325ULL;
    const auto word = [&](std::uint64_t v) {
        for (int i = 0; i < 8; ++i) {
            h ^= (v >> (8 * i)) & 0xffU;
            h *= 0x100000001b3ULL;
        }
    };
    word(0x53495a45464cULL);  // "SIZEFL"
    word(points_.size());
    for (const Point& p : points_) {
        word(std::bit_cast<std::uint64_t>(p[0]));
        word(std::bit_cast<std::uint64_t>(p[1]));
    }
    word(triangles_.size());
    for (const NodeId v : triangles_) word(static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
    for (const double s : sizes_) word(std::bit_cast<std::uint64_t>(s));
    return h;
}

namespace {

// The unique edges of a simplex mesh's cells, sorted.
std::vector<std::pair<NodeId, NodeId>> cell_edges(const SimplexMesh& sm) {
    const std::size_t k = static_cast<std::size_t>(sm.dimension) + 1;
    std::vector<std::pair<NodeId, NodeId>> e;
    if (k < 2) return e;
    e.reserve(sm.cells.size() / k * (k * (k - 1) / 2));
    for (std::size_t c = 0; c + k <= sm.cells.size(); c += k) {
        for (std::size_t i = 0; i < k; ++i) {
            for (std::size_t j = i + 1; j < k; ++j) {
                const NodeId a = sm.cells[c + i], b = sm.cells[c + j];
                e.emplace_back(std::min(a, b), std::max(a, b));
            }
        }
    }
    std::ranges::sort(e);
    e.erase(std::unique(e.begin(), e.end()), e.end());
    return e;
}

}  // namespace

std::vector<double> local_sizes(const SimplexMesh& simplices) {
    std::vector<double> h(simplices.points.size(), 0.0);
    for (const auto& [a, b] : cell_edges(simplices)) {
        const Point& p = simplices.points[static_cast<std::size_t>(a)];
        const Point& q = simplices.points[static_cast<std::size_t>(b)];
        const double len = std::sqrt((q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1]) +
                                     (q[2] - p[2]) * (q[2] - p[2]));
        h[static_cast<std::size_t>(a)] = std::max(h[static_cast<std::size_t>(a)], len);
        h[static_cast<std::size_t>(b)] = std::max(h[static_cast<std::size_t>(b)], len);
    }
    return h;
}

std::expected<void, base::Error> grade_sizes(const SimplexMesh& simplices, std::vector<double>& sizes, double g,
                                             const std::vector<char>& seeds) {
    if (!(std::isfinite(g) && g >= 1.0)) return std::unexpected(invalid("the grading factor must be finite and >= 1"));
    if (sizes.size() != simplices.points.size()) return std::unexpected(invalid("grading needs one size per point"));
    if (!seeds.empty() && seeds.size() != sizes.size()) return std::unexpected(invalid("grading needs one seed flag per point"));
    for (std::size_t i = 0; i < sizes.size(); ++i) {
        if (!(std::isfinite(sizes[i]) && sizes[i] > 0.0)) {
            return std::unexpected(invalid("a size must be finite and positive", i));
        }
    }
    // Dijkstra on log h with edge weight log g: each node takes min over k of g^(dist) h_k. Ties and
    // order are deterministic (the set orders by size, then node).
    const auto edges = cell_edges(simplices);
    std::vector<std::vector<NodeId>> adj(sizes.size());
    for (const auto& [a, b] : edges) {
        adj[static_cast<std::size_t>(a)].push_back(b);
        adj[static_cast<std::size_t>(b)].push_back(a);
    }
    std::set<std::pair<double, std::size_t>> queue;
    for (std::size_t i = 0; i < sizes.size(); ++i) {
        if (seeds.empty() || seeds[i] != 0) queue.emplace(sizes[i], i);
    }
    std::vector<char> done(sizes.size(), 0);
    while (!queue.empty()) {
        const auto [h, i] = *queue.begin();
        queue.erase(queue.begin());
        if (done[i] != 0) continue;
        done[i] = 1;
        const double limit = g * h;
        for (const NodeId jn : adj[i]) {
            const auto j = static_cast<std::size_t>(jn);
            if (done[j] == 0 && limit < sizes[j]) {
                queue.erase({sizes[j], j});
                sizes[j] = limit;
                queue.emplace(limit, j);
            }
        }
    }
    return {};
}

std::expected<AxisRefinement, base::Error> refine_axis(const std::vector<double>& lines,
                                                       const std::vector<char>& marked, double g,
                                                       double min_interval,
                                                       const std::vector<AxisRange>& protected_ranges) {
    if (lines.size() < 2) return std::unexpected(invalid("an axis needs at least two lines"));
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (!std::isfinite(lines[i])) return std::unexpected(invalid("an axis line must be finite", i));
        if (i > 0 && !(lines[i] > lines[i - 1])) return std::unexpected(invalid("axis lines must increase strictly", i));
    }
    if (marked.size() + 1 != lines.size()) return std::unexpected(invalid("one mark per axis interval"));
    if (!(std::isfinite(g) && g >= 1.0)) return std::unexpected(invalid("the grading factor must be finite and >= 1"));
    if (!(std::isfinite(min_interval) && min_interval >= 0.0)) {
        return std::unexpected(invalid("the smallest interval must be finite and >= 0"));
    }
    AxisRefinement out;
    const auto is_protected = [&](double a, double b) {
        for (const AxisRange& r : protected_ranges) {
            if (std::min(b, r.hi) > std::max(a, r.lo)) return true;
        }
        return false;
    };
    // May [a, b] be split? Its halves must stay at least min_interval long.
    const auto splittable = [&](double a, double b) {
        const double m = 0.5 * (a + b);
        return m > a && m < b && m - a >= min_interval && b - m >= min_interval;
    };
    std::vector<double> l;
    l.reserve(2 * lines.size());
    l.push_back(lines[0]);
    for (std::size_t i = 0; i + 1 < lines.size(); ++i) {
        const double a = lines[i], b = lines[i + 1];
        if (marked[i] != 0) {
            if (is_protected(a, b)) {
                ++out.protected_skips;
            } else if (splittable(a, b)) {
                l.push_back(0.5 * (a + b));
                ++out.split;
            }
        }
        l.push_back(b);
    }
    // Grading to a fixed point: split an interval longer than g times a neighbour made by a split
    // in this call (fresh). Ratios already in the given axis (an interface straddle beside coarse
    // intervals) are kept: grading follows the refinement only. Each pass goes left to right; a
    // split never lengthens an interval, so this terminates.
    std::vector<char> fresh;  // per interval of l
    {
        std::size_t j = 0;
        for (std::size_t i = 0; i + 1 < l.size(); ++i) {
            const bool given = j + 1 < lines.size() && l[i] == lines[j] && l[i + 1] == lines[j + 1];
            fresh.push_back(given ? 0 : 1);
            if (l[i + 1] == lines[j + 1]) ++j;
        }
    }
    // A midpoint split always leaves a ratio of 2 beside it, so on an axis the grading factor is at
    // least 2; below that, every split would cascade over the whole axis.
    const double g_axis = std::max(g, 2.0);
    std::set<double> skipped;  // protected intervals counted once, by their left end
    for (bool changed = true; changed;) {
        changed = false;
        std::vector<double> next;
        std::vector<char> next_fresh;
        next.reserve(l.size() * 2);
        next.push_back(l[0]);
        for (std::size_t i = 0; i + 1 < l.size(); ++i) {
            const double len = l[i + 1] - l[i];
            double nb = std::numeric_limits<double>::infinity();
            if (i > 0 && fresh[i - 1] != 0) nb = std::min(nb, l[i] - l[i - 1]);
            if (i + 2 < l.size() && fresh[i + 1] != 0) nb = std::min(nb, l[i + 2] - l[i + 1]);
            if (len > g_axis * nb) {
                if (is_protected(l[i], l[i + 1])) {
                    if (skipped.insert(l[i]).second) ++out.protected_skips;
                } else if (splittable(l[i], l[i + 1])) {
                    next.push_back(0.5 * (l[i] + l[i + 1]));
                    next_fresh.push_back(1);
                    ++out.graded;
                    changed = true;
                    next.push_back(l[i + 1]);
                    next_fresh.push_back(1);
                    continue;
                }
            }
            next.push_back(l[i + 1]);
            next_fresh.push_back(fresh[i]);
        }
        l = std::move(next);
        fresh = std::move(next_fresh);
    }
    out.lines = std::move(l);
    return out;
}

}  // namespace NiTCAD::mesh
