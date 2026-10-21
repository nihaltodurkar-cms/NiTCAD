#include "NiTCAD/render/streamlines.hpp"
#include "NiTCAD/render/vectors.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>

namespace NiTCAD::render {

namespace {

Point sub(const Point& a, const Point& b) noexcept { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Point add(const Point& a, const Point& b) noexcept { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Point scale(const Point& a, double s) noexcept { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Point& a, const Point& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm(const Point& a) noexcept { return std::sqrt(dot(a, a)); }
Point lerp(const Point& a, const Point& b, double t) noexcept { return add(a, scale(sub(b, a), t)); }
Point cross(const Point& a, const Point& b) noexcept {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Places `count` points at the midpoints of equal shares of a distribution over pieces: piece k has
// weight w[k], and at(k, f) is the point a fraction f into it.
template <typename At>
std::vector<Point> quantiles(std::span<const double> w, std::size_t count, At&& at) {
    std::vector<Point> out;
    double total = 0.0;
    for (const double x : w) total += x;
    if (!(total > 0.0) || count == 0) return out;
    std::size_t k = 0;
    double before = 0.0;  // weight of the pieces before k
    for (std::size_t i = 0; i < count; ++i) {
        const double q = (static_cast<double>(i) + 0.5) / static_cast<double>(count) * total;
        while (k + 1 < w.size() && before + w[k] <= q) before += w[k++];
        const double f = w[k] > 0.0 ? std::clamp((q - before) / w[k], 0.0, 1.0) : 0.5;
        out.push_back(at(k, f));
    }
    return out;
}

// Dormand-Prince 5(4) (Hairer, Norsett, Wanner, table 5.2).
constexpr double a21 = 1.0 / 5.0;
constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0,
                 a54 = -212.0 / 729.0;
constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0,
                 a65 = -5103.0 / 18656.0;
constexpr double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0,
                 b6 = 11.0 / 84.0;
// b - b* (the 4th-order weights 5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40).
constexpr double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0,
                 e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;

enum class Sample : std::uint8_t { ok, outside, no_field };

}  // namespace

// --- point location ---------------------------------------------------------------------------

Locator::Locator(const Geometry& g) : g_(&g) {
    const auto points = g.mesh().points();
    const int D = g.dimension();
    const std::size_t E = g.element_count();
    Point lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
    Point hi = scale(lo, -1.0);
    for (const Point& p : points) {
        for (std::size_t a = 0; a < 3; ++a) lo[a] = std::min(lo[a], p[a]), hi[a] = std::max(hi[a], p[a]);
    }
    if (points.empty() || E == 0) return;
    lo_ = lo;
    double volume = 1.0;
    for (int a = 0; a < D; ++a) {
        size_ = std::max(size_, hi[static_cast<std::size_t>(a)] - lo[static_cast<std::size_t>(a)]);
        volume *= std::max(hi[static_cast<std::size_t>(a)] - lo[static_cast<std::size_t>(a)], 1e-300);
    }
    // About one element per bucket.
    cell_ = std::pow(volume / static_cast<double>(E), 1.0 / D);
    for (int a = 0; a < D; ++a) {
        const auto ua = static_cast<std::size_t>(a);
        n_[ua] = std::clamp<std::size_t>(static_cast<std::size_t>((hi[ua] - lo[ua]) / cell_) + 1, 1, 4096);
    }
    const auto bucket = [&](const Point& p, std::size_t a) {
        const double t = (p[a] - lo_[a]) / cell_;
        return std::min(n_[a] - 1, static_cast<std::size_t>(std::max(0.0, t)));
    };
    std::vector<std::uint32_t> count(n_[0] * n_[1] * n_[2] + 1, 0);
    std::vector<std::array<std::size_t, 6>> range(E);
    for (std::size_t e = 0; e < E; ++e) {
        Point blo = points[static_cast<std::size_t>(g.element(e)[0])], bhi = blo;
        for (const NodeId v : g.element(e)) {
            const Point& p = points[static_cast<std::size_t>(v)];
            for (std::size_t a = 0; a < 3; ++a) blo[a] = std::min(blo[a], p[a]), bhi[a] = std::max(bhi[a], p[a]);
        }
        auto& r = range[e];
        for (std::size_t a = 0; a < 3; ++a) r[a] = bucket(blo, a), r[a + 3] = bucket(bhi, a);
        for (std::size_t k = r[2]; k <= r[5]; ++k)
            for (std::size_t j = r[1]; j <= r[4]; ++j)
                for (std::size_t i = r[0]; i <= r[3]; ++i) ++count[(k * n_[1] + j) * n_[0] + i + 1];
    }
    for (std::size_t b = 1; b < count.size(); ++b) count[b] += count[b - 1];
    bucket_start_ = count;
    bucket_items_.resize(count.back());
    for (std::size_t e = 0; e < E; ++e) {
        const auto& r = range[e];
        for (std::size_t k = r[2]; k <= r[5]; ++k)
            for (std::size_t j = r[1]; j <= r[4]; ++j)
                for (std::size_t i = r[0]; i <= r[3]; ++i) bucket_items_[count[(k * n_[1] + j) * n_[0] + i]++] = static_cast<std::uint32_t>(e);
    }
    // Node to elements.
    node_start_.assign(points.size() + 1, 0);
    for (std::size_t e = 0; e < E; ++e) for (const NodeId v : g.element(e)) ++node_start_[static_cast<std::size_t>(v) + 1];
    for (std::size_t i = 1; i < node_start_.size(); ++i) node_start_[i] += node_start_[i - 1];
    node_items_.resize(node_start_.back());
    std::vector<std::uint32_t> fill(node_start_.begin(), node_start_.end() - 1);
    for (std::size_t e = 0; e < E; ++e) {
        for (const NodeId v : g.element(e)) node_items_[fill[static_cast<std::size_t>(v)]++] = static_cast<std::uint32_t>(e);
    }
}

std::optional<Locator::Hit> Locator::inside(std::size_t e, const Point& x) const noexcept {
    const auto n = g_->element(e);
    const auto points = g_->mesh().points();
    const int D = g_->dimension();
    const Point& x0 = points[static_cast<std::size_t>(n[0])];
    // Solve sum_k l_k (x_k - x0) = x - x0 for l_1..l_D (Cramer).
    std::array<Point, 3> col{};
    for (int k = 0; k < D; ++k) col[static_cast<std::size_t>(k)] = sub(points[static_cast<std::size_t>(n[static_cast<std::size_t>(k + 1)])], x0);
    const Point r = sub(x, x0);
    Hit h;
    h.element = e;
    if (D == 1) {
        if (col[0][0] == 0.0) return std::nullopt;
        h.weights[1] = r[0] / col[0][0];
    } else if (D == 2) {
        const double det = col[0][0] * col[1][1] - col[1][0] * col[0][1];
        if (det == 0.0) return std::nullopt;
        h.weights[1] = (r[0] * col[1][1] - col[1][0] * r[1]) / det;
        h.weights[2] = (col[0][0] * r[1] - r[0] * col[0][1]) / det;
    } else {
        const double det = dot(col[0], cross(col[1], col[2]));
        if (det == 0.0) return std::nullopt;
        h.weights[1] = dot(r, cross(col[1], col[2])) / det;
        h.weights[2] = dot(col[0], cross(r, col[2])) / det;
        h.weights[3] = dot(col[0], cross(col[1], r)) / det;
    }
    h.weights[0] = 1.0;
    for (int k = 1; k <= D; ++k) h.weights[0] -= h.weights[static_cast<std::size_t>(k)];
    for (int k = 0; k <= D; ++k) {
        if (h.weights[static_cast<std::size_t>(k)] < -1e-10) return std::nullopt;
    }
    return h;
}

std::optional<Locator::Hit> Locator::locate(const Point& x, std::optional<std::size_t> hint) const {
    if (g_->element_count() == 0) return std::nullopt;
    if (hint) {
        if (auto h = inside(*hint, x)) return h;
        for (const NodeId v : g_->element(*hint)) {
            const auto i = static_cast<std::size_t>(v);
            for (std::uint32_t k = node_start_[i]; k < node_start_[i + 1]; ++k) {
                if (auto h = inside(node_items_[k], x)) return h;
            }
        }
    }
    std::array<std::size_t, 3> b{};
    for (std::size_t a = 0; a < 3; ++a) {
        const double t = (x[a] - lo_[a]) / cell_;
        if (a < static_cast<std::size_t>(g_->dimension()) &&
            (t < -1e-9 || t > static_cast<double>(n_[a]) + 1e-9)) {
            return std::nullopt;
        }
        b[a] = std::min(n_[a] - 1, static_cast<std::size_t>(std::max(0.0, t)));
    }
    const std::size_t cell = (b[2] * n_[1] + b[1]) * n_[0] + b[0];
    for (std::uint32_t k = bucket_start_[cell]; k < bucket_start_[cell + 1]; ++k) {
        if (auto h = inside(bucket_items_[k], x)) return h;
    }
    return std::nullopt;
}

// --- tracing ----------------------------------------------------------------------------------

Streamline trace_streamline(const Locator& loc, std::span<const Point> field, const Point& seed,
                            const StreamlineOptions& o) {
    const Geometry& g = loc.geometry();
    double largest = 0.0;
    for (const Point& v : field) {
        const double m = norm(v);
        if (std::isfinite(m)) largest = std::max(largest, m);
    }
    const double vanish = 1e-12 * largest;
    const double size = loc.size();
    const double tol = o.tolerance * size;
    const double h_min = 1e-7 * size;
    const double max_length = o.max_length > 0.0 ? o.max_length : 20.0 * size;

    std::optional<std::size_t> hint;
    double last = std::numeric_limits<double>::quiet_NaN();  // |v| at the last point sampled
    // The unit direction at x, times sgn.
    const auto direction = [&](const Point& x, double sgn, Point& out) -> Sample {
        const auto h = loc.locate(x, hint);
        if (!h) return Sample::outside;
        hint = h->element;
        Point v{};
        const auto n = g.element(h->element);
        for (std::size_t k = 0; k < n.size(); ++k) {
            v = add(v, scale(field[static_cast<std::size_t>(n[k])], h->weights[k]));
        }
        const double m = norm(v);
        last = m;
        if (!std::isfinite(m) || !(m > vanish)) return Sample::no_field;
        out = scale(v, sgn / m);
        return Sample::ok;
    };

    Streamline line;
    std::vector<Point> back;
    std::vector<double> back_m;
    double seed_m = std::numeric_limits<double>::quiet_NaN();
    for (const double sgn : {1.0, -1.0}) {
        if ((sgn > 0.0 && !o.forward) || (sgn < 0.0 && !o.backward)) continue;
        std::vector<Point>& pts = sgn > 0.0 ? line.points : back;
        std::vector<double>& mags = sgn > 0.0 ? line.magnitude : back_m;
        StreamlineStop& stop = sgn > 0.0 ? line.forward : line.backward;
        hint.reset();
        Point x = seed;
        Point k1{};
        const Sample first = direction(x, sgn, k1);
        seed_m = last;
        double x_m = last;
        if (first != Sample::ok) {
            stop = first == Sample::outside ? StreamlineStop::boundary : StreamlineStop::no_field;
            continue;
        }
        // The first step: a twentieth of the local element's size.
        double h = 0.05 * std::pow(element_volume(g, *hint), 1.0 / g.dimension());
        double s = 0.0;
        std::size_t steps = 0;
        for (;;) {
            if (steps >= o.max_steps) { stop = StreamlineStop::steps; break; }
            if (s >= max_length * (1.0 - 1e-15)) { stop = StreamlineStop::length; break; }
            h = std::min(h, max_length - s);
            Point k2{}, k3{}, k4{}, k5{}, k6{}, k7{};
            const auto at = [&](std::initializer_list<std::pair<double, const Point*>> terms) {
                Point y = x;
                for (const auto& [c, k] : terms) y = add(y, scale(*k, h * c));
                return y;
            };
            Sample fail = Sample::ok;
            if ((fail = direction(at({{a21, &k1}}), sgn, k2)) == Sample::ok &&
                (fail = direction(at({{a31, &k1}, {a32, &k2}}), sgn, k3)) == Sample::ok &&
                (fail = direction(at({{a41, &k1}, {a42, &k2}, {a43, &k3}}), sgn, k4)) == Sample::ok &&
                (fail = direction(at({{a51, &k1}, {a52, &k2}, {a53, &k3}, {a54, &k4}}), sgn, k5)) == Sample::ok &&
                (fail = direction(at({{a61, &k1}, {a62, &k2}, {a63, &k3}, {a64, &k4}, {a65, &k5}}), sgn, k6)) ==
                    Sample::ok) {
                const Point y = at({{b1, &k1}, {b3, &k3}, {b4, &k4}, {b5, &k5}, {b6, &k6}});
                fail = direction(y, sgn, k7);
                const double y_m = last;
                if (fail == Sample::ok) {
                    Point err{};
                    for (const auto& [c, k] : {std::pair{e1, &k1}, std::pair{e3, &k3}, std::pair{e4, &k4},
                                               std::pair{e5, &k5}, std::pair{e6, &k6}, std::pair{e7, &k7}}) {
                        err = add(err, scale(*k, h * c));
                    }
                    const double e = norm(err);
                    const double factor = e > 0.0 ? 0.9 * std::pow(tol / e, 0.2) : 5.0;
                    if (e <= tol || h <= h_min) {
                        if (pts.empty()) {
                            pts.push_back(x);
                            mags.push_back(x_m);
                        }
                        x = y;
                        x_m = y_m;
                        s += h;
                        ++steps;
                        pts.push_back(x);
                        mags.push_back(x_m);
                        // A step that turns the direction back has crossed a sink of the field (where
                        // it reverses: its zero): the line ends there.
                        if (dot(k1, k7) < 0.0) {
                            stop = StreamlineStop::no_field;
                            break;
                        }
                        k1 = k7;  // first same as last
                        // At most half the local element's size: the drawn polyline follows the
                        // field from element to element.
                        h = std::min(h * std::clamp(factor, 0.2, 5.0),
                                     0.5 * std::pow(element_volume(g, *hint), 1.0 / g.dimension()));
                    } else {
                        h *= std::clamp(factor, 0.2, 1.0);
                    }
                    continue;
                }
            }
            // A stage left the mesh or the field: shorten the step, or stop there.
            if (h <= h_min) {
                stop = fail == Sample::outside ? StreamlineStop::boundary : StreamlineStop::no_field;
                break;
            }
            h = std::max(0.25 * h, h_min);
            // The failed stage moved the hint; restore it to x's element.
            hint = loc.locate(x, std::nullopt).transform([](const Locator::Hit& hh) { return hh.element; });
        }
    }
    // Backward part reversed, then the forward part (both start at the seed).
    std::vector<Point> joined(back.rbegin(), back.rend());
    std::vector<double> joined_m(back_m.rbegin(), back_m.rend());
    if (!joined.empty() && !line.points.empty()) {  // the seed once
        joined.pop_back();
        joined_m.pop_back();
    }
    joined.insert(joined.end(), line.points.begin(), line.points.end());
    joined_m.insert(joined_m.end(), line.magnitude.begin(), line.magnitude.end());
    if (joined.empty()) {
        joined.push_back(seed);
        joined_m.push_back(seed_m);
    }
    line.points = std::move(joined);
    line.magnitude = std::move(joined_m);
    return line;
}

// --- seeds ------------------------------------------------------------------------------------

std::vector<double> contact_node_currents(const device::Device& d, const results::BiasPoint& p,
                                          std::size_t contact) {
    const auto edges = d.mesh().edges();
    if (p.edge_current_n.size() != edges.size() || p.edge_current_p.size() != edges.size()) return {};
    const device::Contact& c = d.contacts()[contact];
    std::map<NodeId, std::size_t> index;
    for (std::size_t k = 0; k < c.nodes.size(); ++k) index.emplace(c.nodes[k], k);
    std::vector<double> current(c.nodes.size(), 0.0);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto a = index.find(edges[e].first), b = index.find(edges[e].second);
        if ((a == index.end()) == (b == index.end())) continue;  // both or neither in the contact
        const double i = p.edge_current_n[e] + p.edge_current_p[e];  // first to second
        if (a != index.end()) current[a->second] += i;
        else current[b->second] -= i;
    }
    return current;
}

std::vector<Point> contact_seeds(const device::Device& d, const Locator& loc, std::size_t contact,
                                 std::size_t count, std::span<const double> node_current) {
    const Geometry& g = loc.geometry();
    const auto points = d.mesh().points();
    const device::Contact& c = d.contacts()[contact];
    std::map<NodeId, std::size_t> index;
    for (std::size_t k = 0; k < c.nodes.size(); ++k) index.emplace(c.nodes[k], k);
    const auto in_contact = [&](NodeId v) { return index.contains(v); };

    // The contact's boundary pieces: 2D segments, 3D triangles (node lists).
    std::vector<std::vector<NodeId>> faces;
    if (g.dimension() == 2) {
        std::map<std::pair<NodeId, NodeId>, int> uses;
        for (std::size_t e = 0; e < g.element_count(); ++e) {
            const auto n = g.element(e);
            for (std::size_t i = 0; i < 3; ++i) {
                const NodeId a = n[i], b = n[(i + 1) % 3];
                ++uses[{std::min(a, b), std::max(a, b)}];
            }
        }
        for (const auto& [ab, u] : uses) {
            if (u == 1 && in_contact(ab.first) && in_contact(ab.second)) faces.push_back({ab.first, ab.second});
        }
    } else if (g.dimension() == 3) {
        for (const Polygon& f : g.boundary_faces()) {
            std::vector<NodeId> nodes;
            for (const SurfaceVertex& v : f.vertices) nodes.push_back(v.at.nodes[0]);
            if (std::ranges::all_of(nodes, in_contact)) faces.push_back(std::move(nodes));
        }
    }
    if (faces.empty()) {
        // No boundary piece: the nodes themselves, evenly chosen.
        std::vector<Point> out;
        for (std::size_t i = 0; i < count && !c.nodes.empty(); ++i) {
            out.push_back(points[static_cast<std::size_t>(c.nodes[(i * c.nodes.size()) / count])]);
        }
        return out;
    }
    const auto measure = [&](const std::vector<NodeId>& f) {
        const Point& a = points[static_cast<std::size_t>(f[0])];
        if (f.size() == 2) return norm(sub(points[static_cast<std::size_t>(f[1])], a));
        return 0.5 * norm(cross(sub(points[static_cast<std::size_t>(f[1])], a), sub(points[static_cast<std::size_t>(f[2])], a)));
    };
    // Each node's part of the contact, and its current density over it.
    std::vector<double> part(c.nodes.size(), 0.0);
    for (const auto& f : faces) {
        for (const NodeId v : f) part[index.at(v)] += measure(f) / static_cast<double>(f.size());
    }
    const bool by_current = node_current.size() == c.nodes.size() &&
                            std::ranges::any_of(node_current, [](double i) { return i != 0.0; });
    const auto density = [&](NodeId v) {
        const std::size_t i = index.at(v);
        return by_current ? (part[i] > 0.0 ? std::abs(node_current[i]) / part[i] : 0.0) : 1.0;
    };
    // Pieces: the part of each face nearest each of its nodes, as a path a -> b a seed may lie on.
    struct Piece {
        Point a, b;
    };
    std::vector<Piece> pieces;
    std::vector<double> weight;
    const auto at = [&](NodeId v) { return points[static_cast<std::size_t>(v)]; };
    if (g.dimension() == 2) {
        // Along each chain of segments in order, so the seeds are spaced along the contact: each
        // segment's two halves, the first from its start node to the midpoint, the second on to its end.
        std::map<NodeId, std::vector<std::size_t>> touching;
        for (std::size_t s = 0; s < faces.size(); ++s) for (const NodeId v : faces[s]) touching[v].push_back(s);
        std::vector<std::uint8_t> used(faces.size(), 0);
        const auto walk = [&](NodeId from) {
            for (;;) {
                const auto next = std::ranges::find_if(touching[from], [&](std::size_t s) { return !used[s]; });
                if (next == touching[from].end()) return;
                used[*next] = 1;
                const auto& f = faces[*next];
                const NodeId to = f[0] == from ? f[1] : f[0];
                const Point mid = lerp(at(from), at(to), 0.5);
                const double half = 0.5 * measure(f);
                pieces.push_back({at(from), mid});
                weight.push_back(half * density(from));
                pieces.push_back({mid, at(to)});
                weight.push_back(half * density(to));
                from = to;
            }
        };
        for (const auto& [v, segments] : touching) if (segments.size() == 1) walk(v);  // open chains
        for (const auto& [v, segments] : touching) walk(v);                           // closed loops
    } else {
        // From each node toward the face's centroid, through the node's third of the face.
        for (const auto& f : faces) {
            Point centroid{};
            for (const NodeId v : f) centroid = add(centroid, scale(at(v), 1.0 / 3.0));
            for (const NodeId v : f) {
                pieces.push_back({at(v), centroid});
                weight.push_back(measure(f) / 3.0 * density(v));
            }
        }
    }
    return quantiles(weight, count, [&](std::size_t k, double f) { return lerp(pieces[k].a, pieces[k].b, f); });
}

// --- arrows -----------------------------------------------------------------------------------

namespace {

// The field at x as an arrow, or nothing outside the mesh or where the field is undefined or zero.
std::optional<Arrow> arrow_at(const Locator& loc, std::span<const Point> field, const Point& x,
                              const Point& normal) {
    const auto h = loc.locate(x);
    if (!h) return std::nullopt;
    const auto n = loc.geometry().element(h->element);
    Arrow a;
    a.position = x;
    a.normal = normal;
    Point v{};
    for (std::size_t k = 0; k < n.size(); ++k) {
        v = add(v, scale(field[static_cast<std::size_t>(n[k])], h->weights[k]));
        a.at.nodes[k] = n[k];
        a.at.weights[k] = h->weights[k];
    }
    a.magnitude = norm(v);
    if (!std::isfinite(a.magnitude) || !(a.magnitude > 0.0)) return std::nullopt;
    a.direction = scale(v, 1.0 / a.magnitude);
    return a;
}

}  // namespace

std::vector<Arrow> grid_arrows(const Locator& loc, std::span<const Point> field, std::size_t count,
                               double& spacing) {
    std::vector<Arrow> out;
    spacing = 0.0;
    const Geometry& g = loc.geometry();
    if (g.dimension() != 2 || g.element_count() == 0 || count == 0) return out;
    double area = 0.0;
    for (std::size_t e = 0; e < g.element_count(); ++e) area += element_volume(g, e);
    Point lo = g.mesh().points()[0], hi = lo;
    for (const Point& p : g.mesh().points()) {
        for (std::size_t c = 0; c < 2; ++c) lo[c] = std::min(lo[c], p[c]), hi[c] = std::max(hi[c], p[c]);
    }
    const double h = std::sqrt(area / static_cast<double>(count));
    spacing = h;
    const auto nx = static_cast<std::size_t>(std::max(1.0, std::floor((hi[0] - lo[0]) / h)));
    const auto ny = static_cast<std::size_t>(std::max(1.0, std::floor((hi[1] - lo[1]) / h)));
    // Centred in the box.
    const double x0 = lo[0] + 0.5 * ((hi[0] - lo[0]) - static_cast<double>(nx) * h);
    const double y0 = lo[1] + 0.5 * ((hi[1] - lo[1]) - static_cast<double>(ny) * h);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const Point x{x0 + (static_cast<double>(i) + 0.5) * h, y0 + (static_cast<double>(j) + 0.5) * h, 0.0};
            if (auto a = arrow_at(loc, field, x, Point{})) out.push_back(*a);
        }
    }
    return out;
}

std::vector<Arrow> surface_arrows(const Locator& loc, std::span<const Point> field,
                                  std::span<const Polygon> surface, std::size_t count, double& spacing) {
    std::vector<std::array<Point, 3>> triangles;
    std::vector<Point> normals;
    std::vector<double> area;
    for (const Polygon& p : surface) {
        for (std::size_t i = 1; i + 1 < p.vertices.size(); ++i) {
            const std::array<Point, 3> t{p.vertices[0].position, p.vertices[i].position, p.vertices[i + 1].position};
            const Point n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
            const double a = 0.5 * norm(n);
            if (!(a > 0.0)) continue;
            triangles.push_back(t);
            normals.push_back(scale(n, 0.5 / a));
            area.push_back(a);
        }
    }
    double total = 0.0;
    for (const double a : area) total += a;
    spacing = count > 0 ? std::sqrt(total / static_cast<double>(count)) : 0.0;
    std::vector<std::size_t> which;
    std::size_t i = 0;
    const std::vector<Point> sites = quantiles(area, count, [&](std::size_t k, double f) {
        which.push_back(k);
        const double gold = std::fmod(static_cast<double>(i++) * 0.6180339887498949, 1.0);
        const auto& t = triangles[k];
        return lerp(t[0], lerp(t[1], t[2], gold), std::sqrt(f));
    });
    std::vector<Arrow> out;
    for (std::size_t s = 0; s < sites.size(); ++s) {
        if (auto a = arrow_at(loc, field, sites[s], normals[which[s]])) out.push_back(*a);
    }
    return out;
}

std::vector<Point> segment_seeds(const Point& a, const Point& b, std::size_t count) {
    std::vector<Point> out;
    for (std::size_t i = 0; i < count; ++i) out.push_back(lerp(a, b, (static_cast<double>(i) + 0.5) / static_cast<double>(count)));
    return out;
}

std::vector<Point> plane_seeds(const Geometry& g, const Plane& plane, std::size_t count) {
    std::vector<std::array<Point, 3>> triangles;
    std::vector<double> area;
    for (const Polygon& p : g.slice(plane)) {
        for (std::size_t i = 1; i + 1 < p.vertices.size(); ++i) {
            const std::array<Point, 3> t{p.vertices[0].position, p.vertices[i].position, p.vertices[i + 1].position};
            triangles.push_back(t);
            area.push_back(0.5 * norm(cross(sub(t[1], t[0]), sub(t[2], t[0]))));
        }
    }
    std::size_t i = 0;
    return quantiles(area, count, [&](std::size_t k, double f) {
        // Within the triangle: sqrt(f) toward the far side, at a golden-ratio fraction along it.
        const double gold = std::fmod(static_cast<double>(i++) * 0.6180339887498949, 1.0);
        const auto& t = triangles[k];
        return lerp(t[0], lerp(t[1], t[2], gold), std::sqrt(f));
    });
}

}  // namespace NiTCAD::render
