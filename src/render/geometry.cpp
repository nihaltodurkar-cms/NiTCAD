#include "NiTCAD/render/geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <unordered_map>
#include <utility>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::render {

namespace {

using Vec = Point;

Vec sub(const Vec& a, const Vec& b) noexcept { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec add(const Vec& a, const Vec& b) noexcept { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec scale(const Vec& a, double s) noexcept { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec& a, const Vec& b) noexcept { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec cross(const Vec& a, const Vec& b) noexcept {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Vec lerp(const Vec& a, const Vec& b, double t) noexcept { return add(a, scale(sub(b, a), t)); }
double norm(const Vec& a) noexcept { return std::sqrt(dot(a, a)); }

base::Error invalid(std::string message, std::size_t index) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

// Barycentric coordinates are kept as four weights (unused ones zero).
using Bary = std::array<double, 4>;

// Sutherland-Hodgman: the part of a convex polygon (or, with two vertices, a segment) where f >= 0.
// V is the vertex type; f evaluates a vertex; mix interpolates two vertices.
template <class V, class F, class Mix>
std::vector<V> clip(const std::vector<V>& in, F&& f, Mix&& mix) {
    std::vector<V> out;
    const std::size_t n = in.size();
    if (n == 2) {
        const double a = f(in[0]), b = f(in[1]);
        if (a >= 0.0 && b >= 0.0) return in;
        if (a < 0.0 && b < 0.0) return out;
        const V cut = mix(in[0], in[1], a / (a - b));
        if (a >= 0.0) return {in[0], cut};
        return {cut, in[1]};
    }
    for (std::size_t i = 0; i < n; ++i) {
        const V& p = in[i];
        const V& q = in[(i + 1) % n];
        const double a = f(p), b = f(q);
        if (a >= 0.0) out.push_back(p);
        if ((a >= 0.0) != (b >= 0.0)) out.push_back(mix(p, q, a / (a - b)));
    }
    return out;
}

// The measure of a polygon's point list: length for two points, area for more.
double measure(std::span<const Point> pts) noexcept {
    if (pts.size() < 2) return 0.0;
    if (pts.size() == 2) return norm(sub(pts[1], pts[0]));
    return polygon_area(pts);
}

// Orders the points of a convex planar polygon around their centroid (deterministic).
void order_convex(std::vector<SurfaceVertex>& v, const Vec& normal) {
    if (v.size() < 4) {
        if (v.size() == 3) {
            const Vec n = cross(sub(v[1].position, v[0].position), sub(v[2].position, v[0].position));
            if (dot(n, normal) < 0.0) std::swap(v[1], v[2]);
        }
        return;
    }
    Vec c{};
    for (const SurfaceVertex& s : v) c = add(c, s.position);
    c = scale(c, 1.0 / static_cast<double>(v.size()));
    Vec u = sub(v[0].position, c);
    const double lu = norm(u);
    if (lu == 0.0) return;
    u = scale(u, 1.0 / lu);
    const Vec w = cross(normal, u);
    std::vector<std::pair<double, SurfaceVertex>> keyed;
    for (const SurfaceVertex& s : v) {
        const Vec d = sub(s.position, c);
        keyed.emplace_back(std::atan2(dot(d, w), dot(d, u)), s);
    }
    std::ranges::stable_sort(keyed, {}, &std::pair<double, SurfaceVertex>::first);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = keyed[i].second;
}

}  // namespace

// --- construction -----------------------------------------------------------------------------

std::expected<Geometry, base::Error> Geometry::tensor(const mesh::Mesh& m,
                                                      const mesh::TensorCells& cells) {
    if (!cells.matches(m)) {
        return std::unexpected(invalid("the tensor cells are not those of the mesh", 0));
    }
    Geometry g;
    g.mesh_ = &m;
    g.kind_ = GeometryKind::tensor;
    g.split_tensor(cells);
    return g;
}

void Geometry::split_tensor(const mesh::TensorCells& cells) {
    const int D = cells.dimension();
    per_element_ = D + 1;
    const std::size_t nx = cells.size(0) - 1;
    const std::size_t ny = D >= 2 ? cells.size(1) - 1 : 1;
    const std::size_t nz = D >= 3 ? cells.size(2) - 1 : 1;
    const int corners = 1 << D;
    const auto points = mesh_->points();
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t c = cell_lo_.size();
                std::array<NodeId, 8> corner{};
                for (int q = 0; q < corners; ++q) {
                    corner[static_cast<std::size_t>(q)] = cells.node(
                        {i + static_cast<std::size_t>(q & 1),
                         D >= 2 ? j + static_cast<std::size_t>((q >> 1) & 1) : 0,
                         D >= 3 ? k + static_cast<std::size_t>((q >> 2) & 1) : 0});
                    cell_corners_.push_back(corner[static_cast<std::size_t>(q)]);
                }
                cell_lo_.push_back(points[static_cast<std::size_t>(corner[0])]);
                cell_hi_.push_back(points[static_cast<std::size_t>(corner[static_cast<std::size_t>(corners - 1)])]);
                const auto emit = [&](std::initializer_list<int> q) {
                    for (const int x : q) elements_.push_back(corner[static_cast<std::size_t>(x)]);
                    element_cell_.push_back(c);
                };
                if (D == 1) {
                    emit({0, 1});
                } else if (D == 2) {
                    emit({0, 1, 3});
                    emit({0, 3, 2});
                } else {
                    // Kuhn: 000 -> e_a -> e_a + e_b -> 111 for each order (a, b, c) of the axes.
                    constexpr std::array<std::array<int, 3>, 6> orders{
                        {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}};
                    for (const auto& o : orders) {
                        const int v1 = 1 << o[0];
                        const int v2 = v1 | (1 << o[1]);
                        std::array<int, 4> t{0, v1, v2, 7};
                        const Vec p0 = points[static_cast<std::size_t>(corner[0])];
                        const Vec a = sub(points[static_cast<std::size_t>(corner[static_cast<std::size_t>(t[1])])], p0);
                        const Vec b = sub(points[static_cast<std::size_t>(corner[static_cast<std::size_t>(t[2])])], p0);
                        const Vec d = sub(points[static_cast<std::size_t>(corner[7])], p0);
                        if (dot(cross(a, b), d) < 0.0) std::swap(t[1], t[2]);
                        emit({t[0], t[1], t[2], t[3]});
                    }
                }
            }
        }
    }
}

std::expected<Geometry, base::Error> Geometry::simplex(const mesh::Mesh& m,
                                                       const mesh::SimplexMesh& s) {
    if (s.dimension != m.dimension()) {
        return std::unexpected(invalid("the simplices' dimension is not the mesh's",
                                       static_cast<std::size_t>(s.dimension)));
    }
    if (s.points.size() != m.node_count()) {
        return std::unexpected(invalid("the simplices' point count is not the mesh's node count",
                                       s.points.size()));
    }
    const auto points = m.points();
    for (std::size_t i = 0; i < s.points.size(); ++i) {
        if (s.points[i] != points[i]) {
            return std::unexpected(invalid("a simplex point is not the mesh node's", i));
        }
    }
    const std::size_t per = static_cast<std::size_t>(s.dimension) + 1;
    if (s.cells.size() % per != 0) {
        return std::unexpected(invalid("the cell list is not whole cells", s.cells.size()));
    }
    Geometry g;
    g.mesh_ = &m;
    g.kind_ = GeometryKind::simplex;
    g.per_element_ = static_cast<int>(per);
    g.elements_ = s.cells;
    for (std::size_t e = 0; e * per < g.elements_.size(); ++e) {
        NodeId* c = g.elements_.data() + e * per;
        for (std::size_t k = 0; k < per; ++k) {
            if (c[k] < 0 || static_cast<std::size_t>(c[k]) >= m.node_count()) {
                return std::unexpected(invalid("a cell node is out of range", e));
            }
        }
        const Vec p0 = points[static_cast<std::size_t>(c[0])];
        const Vec a = sub(points[static_cast<std::size_t>(c[1])], p0);
        const Vec b = sub(points[static_cast<std::size_t>(c[2])], p0);
        const double orientation = per == 3 ? a[0] * b[1] - a[1] * b[0]
                                            : dot(cross(a, b), sub(points[static_cast<std::size_t>(c[3])], p0));
        if (orientation < 0.0) std::swap(c[1], c[2]);
    }
    return g;
}

Geometry Geometry::graph(const mesh::Mesh& m) {
    Geometry g;
    g.mesh_ = &m;
    g.kind_ = GeometryKind::graph;
    if (m.dimension() == 1) {
        g.per_element_ = 2;
        for (const mesh::Edge& e : m.edges()) {
            g.elements_.push_back(e.first);
            g.elements_.push_back(e.second);
        }
    }
    return g;
}

// --- access -----------------------------------------------------------------------------------

std::size_t Geometry::element_count() const noexcept {
    return per_element_ == 0 ? 0 : elements_.size() / static_cast<std::size_t>(per_element_);
}

std::span<const NodeId> Geometry::element(std::size_t e) const noexcept {
    const auto per = static_cast<std::size_t>(per_element_);
    return std::span<const NodeId>(elements_).subspan(e * per, per);
}

std::pair<Point, Point> Geometry::cell_box(std::size_t c) const noexcept {
    return {cell_lo_[c], cell_hi_[c]};
}

std::span<const NodeId> Geometry::cell_corners(std::size_t c) const noexcept {
    const auto per = std::size_t{1} << dimension();
    return std::span<const NodeId>(cell_corners_).subspan(c * per, per);
}

Polygon Geometry::element_polygon(std::size_t e) const {
    Polygon p;
    p.element = e;
    const auto nodes = element(e);
    for (const NodeId n : nodes) {
        SurfaceVertex v;
        v.position = mesh_->points()[static_cast<std::size_t>(n)];
        v.at.nodes[0] = n;
        v.at.weights[0] = 1.0;
        p.vertices.push_back(v);
    }
    return p;
}

// --- boundary ---------------------------------------------------------------------------------

std::vector<Polygon> Geometry::boundary_faces() const {
    std::vector<Polygon> out;
    const std::size_t count = element_count();
    if (dimension() == 2) {
        for (std::size_t e = 0; e < count; ++e) out.push_back(element_polygon(e));
        return out;
    }
    if (dimension() != 3) return out;
    struct Face {
        std::array<NodeId, 3> key;      // sorted
        std::array<NodeId, 3> oriented; // outward from its element
        std::size_t element;
    };
    std::vector<Face> faces;
    faces.reserve(count * 4);
    for (std::size_t e = 0; e < count; ++e) {
        const auto v = element(e);
        const std::array<std::array<NodeId, 3>, 4> f{
            {{v[1], v[2], v[3]}, {v[0], v[3], v[2]}, {v[0], v[1], v[3]}, {v[0], v[2], v[1]}}};
        for (const auto& o : f) {
            std::array<NodeId, 3> k = o;
            std::ranges::sort(k);
            faces.push_back({k, o, e});
        }
    }
    std::ranges::sort(faces, [](const Face& a, const Face& b) {
        return a.key != b.key ? a.key < b.key : a.element < b.element;
    });
    const auto points = mesh_->points();
    for (std::size_t i = 0; i < faces.size();) {
        std::size_t j = i + 1;
        while (j < faces.size() && faces[j].key == faces[i].key) ++j;
        if (j - i == 1) {
            Polygon p;
            p.element = faces[i].element;
            for (const NodeId n : faces[i].oriented) {
                SurfaceVertex sv;
                sv.position = points[static_cast<std::size_t>(n)];
                sv.at.nodes[0] = n;
                sv.at.weights[0] = 1.0;
                p.vertices.push_back(sv);
            }
            out.push_back(std::move(p));
        }
        i = j;
    }
    return out;
}

// --- slices -----------------------------------------------------------------------------------

std::vector<Polygon> Geometry::slice(const Plane& plane) const {
    std::vector<Polygon> out;
    const int D = dimension();
    if (D < 2) return out;
    const auto points = mesh_->points();
    const std::size_t count = element_count();
    const double nn = norm(plane.normal);
    if (nn == 0.0) return out;
    for (std::size_t e = 0; e < count; ++e) {
        const auto v = element(e);
        const std::size_t k = v.size();
        std::array<double, 4> s{};
        std::array<Vec, 4> x{};
        int zeros = 0, positive = 0, negative = 0;
        double scale_s = 0.0;
        for (std::size_t i = 0; i < k; ++i) {
            x[i] = points[static_cast<std::size_t>(v[i])];
            s[i] = dot(plane.normal, x[i]) - plane.offset;
            scale_s = std::max(scale_s, std::abs(s[i]));
        }
        for (std::size_t i = 0; i < k; ++i) {
            if (s[i] == 0.0) ++zeros;
            else if (s[i] > 0.0) ++positive;
            else ++negative;
        }
        // A facet (3D face, 2D edge) lying in the plane is taken from the element on its positive
        // side only, so it appears once.
        const int facet = D;  // vertices of a facet
        if (zeros == facet && negative == 1) continue;
        if (zeros < facet && (positive == 0 || negative == 0)) continue;
        Polygon p;
        p.element = e;
        for (std::size_t i = 0; i < k; ++i) {
            if (s[i] == 0.0) {
                SurfaceVertex sv;
                sv.position = x[i];
                sv.at.nodes[0] = v[i];
                sv.at.weights[0] = 1.0;
                p.vertices.push_back(sv);
            }
        }
        for (std::size_t i = 0; i < k; ++i) {
            for (std::size_t j = i + 1; j < k; ++j) {
                if (!((s[i] > 0.0 && s[j] < 0.0) || (s[i] < 0.0 && s[j] > 0.0))) continue;
                const double t = s[i] / (s[i] - s[j]);
                SurfaceVertex sv;
                sv.position = lerp(x[i], x[j], t);
                sv.at.nodes = {v[i], v[j], 0, 0};
                sv.at.weights = {1.0 - t, t, 0.0, 0.0};
                p.vertices.push_back(sv);
            }
        }
        if (static_cast<int>(p.vertices.size()) < D) continue;
        if (D == 3) order_convex(p.vertices, plane.normal);
        out.push_back(std::move(p));
    }
    return out;
}

// --- clipping box -----------------------------------------------------------------------------

std::vector<std::vector<Point>> outside_box(std::span<const Point> polygon, const Box& box) {
    std::vector<std::vector<Point>> out;
    std::vector<Point> rest(polygon.begin(), polygon.end());
    const double whole = measure(rest);
    if (!(whole > 0.0)) return out;
    const auto mix = [](const Point& a, const Point& b, double t) { return lerp(a, b, t); };
    for (std::size_t a = 0; a < 3 && rest.size() >= 3; ++a) {
        for (const bool upper : {false, true}) {
            // Beyond this face: x_a < lo_a, or x_a > hi_a.
            const auto beyond = [&](const Point& p) { return upper ? p[a] - box.hi[a] : box.lo[a] - p[a]; };
            std::vector<Point> part = clip(rest, beyond, mix);
            // A part lying in the face's plane (nothing strictly beyond it) belongs to the box.
            double depth = 0.0;
            for (const Point& p : part) depth = std::max(depth, beyond(p));
            if (part.size() >= 3 && depth > 0.0 && measure(part) > 1e-12 * whole) out.push_back(std::move(part));
            rest = clip(rest, [&](const Point& p) { return -beyond(p); }, mix);
            if (rest.size() < 3) break;
        }
    }
    return out;
}

std::vector<Polygon> Geometry::clipped_surface(const Box& box) const {
    std::vector<Polygon> out;
    if (dimension() != 3 || element_count() == 0) return out;
    const auto with_interpolants = [&](std::vector<Point> pts, std::size_t element) {
        Polygon p;
        p.element = element;
        for (const Point& x : pts) p.vertices.push_back({x, interpolant(element, x)});
        return p;
    };
    // The boundary's parts outside the box.
    for (const Polygon& f : boundary_faces()) {
        std::vector<Point> pts;
        for (const SurfaceVertex& v : f.vertices) pts.push_back(v.position);
        for (auto& part : outside_box(pts, box)) out.push_back(with_interpolants(std::move(part), f.element));
    }
    // The cut faces: on each box face strictly inside the mesh's bounding box.
    const auto points = mesh_->points();
    Point lo = points[0], hi = points[0];
    for (const Point& p : points) {
        for (std::size_t a = 0; a < 3; ++a) lo[a] = std::min(lo[a], p[a]), hi[a] = std::max(hi[a], p[a]);
    }
    const auto mix = [](const Point& a, const Point& b, double t) { return lerp(a, b, t); };
    for (std::size_t a = 0; a < 3; ++a) {
        for (const bool upper : {false, true}) {
            const double at = upper ? box.hi[a] : box.lo[a];
            if (!(at > lo[a] && at < hi[a])) continue;
            // What is left lies beyond the face; its outward normal points into the box.
            Plane plane;
            plane.normal[a] = upper ? -1.0 : 1.0;
            plane.offset = upper ? -at : at;
            for (const Polygon& s : slice(plane)) {
                std::vector<Point> pts;
                for (const SurfaceVertex& v : s.vertices) pts.push_back(v.position);
                const double whole = measure(pts);
                for (std::size_t b = 0; b < 3 && pts.size() >= 3; ++b) {
                    if (b == a) continue;
                    pts = clip(pts, [&](const Point& p) { return p[b] - box.lo[b]; }, mix);
                    pts = clip(pts, [&](const Point& p) { return box.hi[b] - p[b]; }, mix);
                }
                if (pts.size() >= 3 && measure(pts) > 1e-12 * whole) out.push_back(with_interpolants(std::move(pts), s.element));
            }
        }
    }
    return out;
}

// --- control-volume pieces --------------------------------------------------------------------

std::vector<OwnedPolygon> Geometry::owner_pieces(const Polygon& polygon) const {
    std::vector<OwnedPolygon> out;
    std::vector<Point> pts;
    for (const SurfaceVertex& s : polygon.vertices) pts.push_back(s.position);
    const double whole = measure(pts);
    if (whole == 0.0) return out;
    const double tiny = whole * 1e-12;
    const auto mix = [](const Point& a, const Point& b, double t) { return lerp(a, b, t); };
    const int D = dimension();
    if (kind_ == GeometryKind::tensor) {
        const std::size_t c = element_cell_[polygon.element];
        const Point lo = cell_lo_[c], hi = cell_hi_[c];
        const auto corners = cell_corners(c);
        for (std::size_t q = 0; q < corners.size(); ++q) {
            std::vector<Point> piece = pts;
            for (int a = 0; a < D && !piece.empty(); ++a) {
                const auto ax = static_cast<std::size_t>(a);
                const double mid = 0.5 * (lo[ax] + hi[ax]);
                const bool upper = ((q >> a) & 1) != 0;
                piece = clip(piece, [&](const Point& p) { return upper ? p[ax] - mid : mid - p[ax]; }, mix);
            }
            if (piece.size() >= 2 && measure(piece) > tiny) out.push_back({std::move(piece), corners[q]});
        }
        return out;
    }
    // Barycentric: lambda = M^-1 (x - v0) for the element's vertices.
    const auto v = element(polygon.element);
    const auto points = mesh_->points();
    const Vec p0 = points[static_cast<std::size_t>(v[0])];
    const std::size_t k = v.size();
    std::array<Vec, 3> col{};
    for (std::size_t i = 1; i < k; ++i) col[i - 1] = sub(points[static_cast<std::size_t>(v[i])], p0);
    const auto bary = [&](const Point& x) {
        const Vec r = sub(x, p0);
        Bary l{};
        if (k == 2) {
            const double t = dot(r, col[0]) / dot(col[0], col[0]);
            l = {1.0 - t, t, 0.0, 0.0};
        } else if (k == 3) {
            const double det = col[0][0] * col[1][1] - col[0][1] * col[1][0];
            const double a = (r[0] * col[1][1] - r[1] * col[1][0]) / det;
            const double b = (col[0][0] * r[1] - col[0][1] * r[0]) / det;
            l = {1.0 - a - b, a, b, 0.0};
        } else {
            const double det = dot(col[0], cross(col[1], col[2]));
            const double a = dot(r, cross(col[1], col[2])) / det;
            const double b = dot(col[0], cross(r, col[2])) / det;
            const double c = dot(col[0], cross(col[1], r)) / det;
            l = {1.0 - a - b - c, a, b, c};
        }
        return l;
    };
    struct BV {
        Point x;
        Bary l;
    };
    std::vector<BV> start;
    for (const Point& x : pts) start.push_back({x, bary(x)});
    const auto bmix = [](const BV& a, const BV& b, double t) {
        BV r{lerp(a.x, b.x, t), {}};
        for (std::size_t i = 0; i < 4; ++i) r.l[i] = a.l[i] + t * (b.l[i] - a.l[i]);
        return r;
    };
    for (std::size_t q = 0; q < k; ++q) {
        std::vector<BV> piece = start;
        for (std::size_t j = 0; j < k && !piece.empty(); ++j) {
            if (j == q) continue;
            piece = clip(piece, [&](const BV& b) { return b.l[q] - b.l[j]; }, bmix);
        }
        std::vector<Point> xs;
        for (const BV& b : piece) xs.push_back(b.x);
        if (xs.size() >= 2 && measure(xs) > tiny) out.push_back({std::move(xs), v[q]});
    }
    return out;
}

std::vector<InterfacePiece> Geometry::interfaces(std::span<const std::int32_t> node_region) const {
    NITCAD_EXPECTS(node_region.size() == mesh_->node_count());
    std::vector<InterfacePiece> out;
    const int D = dimension();
    if (D < 2 || element_count() == 0) return out;
    const auto region = [&](NodeId n) { return node_region[static_cast<std::size_t>(n)]; };
    if (kind_ == GeometryKind::tensor) {
        for (std::size_t c = 0; c < cell_count(); ++c) {
            const auto corners = cell_corners(c);
            const Point lo = cell_lo_[c], hi = cell_hi_[c];
            Point mid{};
            for (std::size_t a = 0; a < 3; ++a) mid[a] = 0.5 * (lo[a] + hi[a]);
            for (std::size_t q = 0; q < corners.size(); ++q) {
                for (int a = 0; a < D; ++a) {
                    if ((q >> a) & 1) continue;
                    const std::size_t q2 = q | (std::size_t{1} << a);
                    std::int32_t ra = region(corners[q]), rb = region(corners[q2]);
                    if (ra == rb) continue;
                    if (ra > rb) std::swap(ra, rb);
                    // The face x_a = mid between the two sub-boxes, over corner q's half of the
                    // other axes.
                    std::array<std::pair<double, double>, 3> span{};
                    for (int b = 0; b < 3; ++b) {
                        const auto bx = static_cast<std::size_t>(b);
                        if (b == a || b >= D) {
                            span[bx] = {mid[bx], mid[bx]};
                            if (b >= D) span[bx] = {lo[bx], lo[bx]};
                        } else if ((q >> b) & 1) {
                            span[bx] = {mid[bx], hi[bx]};
                        } else {
                            span[bx] = {lo[bx], mid[bx]};
                        }
                    }
                    const bool flip = region(corners[q]) != ra;
                    const std::size_t per_cell = element_count() / cell_count();
                    InterfacePiece piece{{}, ra, rb, flip ? corners[q2] : corners[q],
                                         flip ? corners[q] : corners[q2], c * per_cell};
                    const auto ax = static_cast<std::size_t>(a);
                    if (D == 2) {
                        const std::size_t b = 1 - ax;
                        Point p{}, r{};
                        p[ax] = r[ax] = mid[ax];
                        p[b] = span[b].first;
                        r[b] = span[b].second;
                        piece.points = {p, r};
                    } else {
                        const std::size_t b = (ax + 1) % 3, d = (ax + 2) % 3;
                        for (const auto& [u, w] : {std::pair{0, 0}, std::pair{1, 0}, std::pair{1, 1},
                                                   std::pair{0, 1}}) {
                            Point p{};
                            p[ax] = mid[ax];
                            p[b] = u ? span[b].second : span[b].first;
                            p[d] = w ? span[d].second : span[d].first;
                            piece.points.push_back(p);
                        }
                    }
                    out.push_back(std::move(piece));
                }
            }
        }
        return out;
    }
    // Simplices: the plane lambda_i = lambda_j within the element, where i and j weigh most.
    const auto points = mesh_->points();
    for (std::size_t e = 0; e < element_count(); ++e) {
        const auto v = element(e);
        const std::size_t k = v.size();
        bool mixed = false;
        for (std::size_t i = 1; i < k; ++i) mixed = mixed || region(v[i]) != region(v[0]);
        if (!mixed) continue;
        for (std::size_t i = 0; i < k; ++i) {
            for (std::size_t j = i + 1; j < k; ++j) {
                std::int32_t ra = region(v[i]), rb = region(v[j]);
                if (ra == rb) continue;
                if (ra > rb) std::swap(ra, rb);
                // In barycentric coordinates: the midpoint of edge ij and the other vertices.
                std::vector<Bary> section;
                Bary m{};
                m[i] = m[j] = 0.5;
                section.push_back(m);
                for (std::size_t q = 0; q < k; ++q) {
                    if (q == i || q == j) continue;
                    Bary b{};
                    b[q] = 1.0;
                    section.push_back(b);
                }
                const auto mix = [](const Bary& a, const Bary& b, double t) {
                    Bary r{};
                    for (std::size_t z = 0; z < 4; ++z) r[z] = a[z] + t * (b[z] - a[z]);
                    return r;
                };
                for (std::size_t q = 0; q < k && section.size() >= 2; ++q) {
                    if (q == i || q == j) continue;
                    section = clip(section, [&](const Bary& b) { return b[i] - b[q]; }, mix);
                }
                std::vector<Point> xs;
                for (const Bary& b : section) {
                    Point x{};
                    for (std::size_t z = 0; z < k; ++z) {
                        x = add(x, scale(points[static_cast<std::size_t>(v[z])], b[z]));
                    }
                    xs.push_back(x);
                }
                if (static_cast<int>(xs.size()) >= D && measure(xs) > 0.0) {
                    const bool flip = region(v[i]) != ra;
                    out.push_back({std::move(xs), ra, rb, flip ? v[j] : v[i], flip ? v[i] : v[j], e});
                }
            }
        }
    }
    return out;
}

Interpolant Geometry::interpolant(std::size_t e, const Point& x) const {
    std::size_t first = e, last = e + 1;
    if (kind_ == GeometryKind::tensor) {
        const std::size_t per_cell = element_count() / cell_count();
        first = element_cell_[e] * per_cell;
        last = first + per_cell;
    }
    const auto points = mesh_->points();
    Interpolant best;
    double best_min = -std::numeric_limits<double>::infinity();
    for (std::size_t q = first; q < last; ++q) {
        const auto v = element(q);
        const std::size_t k = v.size();
        const Vec p0 = points[static_cast<std::size_t>(v[0])];
        const Vec r = sub(x, p0);
        std::array<Vec, 3> col{};
        for (std::size_t i = 1; i < k; ++i) col[i - 1] = sub(points[static_cast<std::size_t>(v[i])], p0);
        Bary l{};
        if (k == 2) {
            const double t = dot(r, col[0]) / dot(col[0], col[0]);
            l = {1.0 - t, t, 0.0, 0.0};
        } else if (k == 3) {
            const double det = col[0][0] * col[1][1] - col[0][1] * col[1][0];
            const double a = (r[0] * col[1][1] - r[1] * col[1][0]) / det;
            const double b = (col[0][0] * r[1] - col[0][1] * r[0]) / det;
            l = {1.0 - a - b, a, b, 0.0};
        } else {
            const double det = dot(col[0], cross(col[1], col[2]));
            const double a = dot(r, cross(col[1], col[2])) / det;
            const double b = dot(col[0], cross(r, col[2])) / det;
            const double c = dot(col[0], cross(col[1], r)) / det;
            l = {1.0 - a - b - c, a, b, c};
        }
        double lowest = l[0];
        for (std::size_t i = 1; i < k; ++i) lowest = std::min(lowest, l[i]);
        if (lowest > best_min) {
            best_min = lowest;
            best = Interpolant{};
            for (std::size_t i = 0; i < k; ++i) {
                best.nodes[i] = v[i];
                best.weights[i] = l[i];
            }
        }
    }
    return best;
}

// --- free functions ---------------------------------------------------------------------------

double interpolate(const Interpolant& at, std::span<const double> field) noexcept {
    double v = 0.0;
    for (std::size_t k = 0; k < 4; ++k) {
        if (at.weights[k] != 0.0) v += at.weights[k] * field[static_cast<std::size_t>(at.nodes[k])];
    }
    return v;
}

std::optional<std::array<Point, 2>> contour_segment(const std::array<Point, 3>& p,
                                                    const std::array<double, 3>& v,
                                                    double level) noexcept {
    std::array<bool, 3> above{};
    for (std::size_t i = 0; i < 3; ++i) {
        if (!std::isfinite(v[i])) return std::nullopt;
        above[i] = v[i] >= level;
    }
    std::array<Point, 2> seg{};
    int found = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        const std::size_t j = (i + 1) % 3;
        if (above[i] == above[j]) continue;
        const double t = (level - v[i]) / (v[j] - v[i]);
        seg[static_cast<std::size_t>(found++)] = lerp(p[i], p[j], t);
    }
    if (found != 2) return std::nullopt;
    return seg;
}

IsoSurface isosurface(const Geometry& g, std::span<const double> values, double level) {
    IsoSurface out;
    if (g.dimension() != 3 || g.element_count() == 0) return out;
    const auto points = g.mesh().points();
    // A vertex per cut edge (lo, hi), or per node (n, n) for a crossing at the node.
    std::unordered_map<std::uint64_t, std::uint32_t> vertex_of;
    const auto vertex = [&](NodeId i, NodeId j) -> std::uint32_t {
        const NodeId a = std::min(i, j), b = std::max(i, j);
        const double va = values[static_cast<std::size_t>(a)], vb = values[static_cast<std::size_t>(b)];
        // Exactly at a node: that node.
        NodeId at = -1;
        if (va == level) at = a;
        else if (vb == level) at = b;
        const NodeId ka = at >= 0 ? at : a, kb = at >= 0 ? at : b;
        const std::uint64_t key = (static_cast<std::uint64_t>(kb) << 32) | static_cast<std::uint32_t>(ka);
        const auto [it, fresh] = vertex_of.try_emplace(key, static_cast<std::uint32_t>(out.vertices.size()));
        if (!fresh) return it->second;
        SurfaceVertex v;
        if (at >= 0) {
            v.position = points[static_cast<std::size_t>(at)];
            v.at.nodes = {at, 0, 0, 0};
            v.at.weights = {1.0, 0.0, 0.0, 0.0};
        } else {
            // From the lower node, so every element around the edge computes the same point.
            const double t = (level - va) / (vb - va);
            v.position = lerp(points[static_cast<std::size_t>(a)], points[static_cast<std::size_t>(b)], t);
            v.at.nodes = {a, b, 0, 0};
            v.at.weights = {1.0 - t, t, 0.0, 0.0};
        }
        out.vertices.push_back(v);
        return it->second;
    };
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        const auto n = g.element(e);
        std::array<double, 4> v{};
        bool finite = true;
        for (std::size_t k = 0; k < 4; ++k) {
            v[k] = values[static_cast<std::size_t>(n[k])];
            finite = finite && std::isfinite(v[k]);
        }
        if (!finite) continue;
        std::array<std::size_t, 4> above{}, below{};
        std::size_t na = 0, nb = 0;
        for (std::size_t k = 0; k < 4; ++k) (v[k] >= level ? above[na++] : below[nb++]) = k;
        if (na == 0 || nb == 0) continue;
        std::vector<std::uint32_t> poly;
        if (na == 1 || nb == 1) {
            const std::size_t odd = na == 1 ? above[0] : below[0];
            for (std::size_t k = 0; k < 4; ++k) {
                if (k != odd) poly.push_back(vertex(n[odd], n[k]));
            }
        } else {
            // Around the quadrilateral: each consecutive pair of cut edges shares a node.
            poly = {vertex(n[above[0]], n[below[0]]), vertex(n[above[0]], n[below[1]]),
                    vertex(n[above[1]], n[below[1]]), vertex(n[above[1]], n[below[0]])};
        }
        for (std::size_t i = 1; i + 1 < poly.size(); ++i) {
            std::array<std::uint32_t, 3> t{poly[0], poly[i], poly[i + 1]};
            if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) continue;
            const Vec& p0 = out.vertices[t[0]].position;
            const Vec normal = cross(sub(out.vertices[t[1]].position, p0), sub(out.vertices[t[2]].position, p0));
            // Toward higher values: sum_k (v_k - level)(x_k - p0) has the gradient's direction's sign
            // (for a linear field its dot with the gradient is sum_k (v_k - level)^2).
            Vec up{};
            for (std::size_t k = 0; k < 4; ++k) {
                up = add(up, scale(sub(points[static_cast<std::size_t>(n[k])], p0), v[k] - level));
            }
            if (dot(normal, up) < 0.0) std::swap(t[1], t[2]);
            out.triangles.push_back(t);
        }
    }
    return out;
}

double polygon_area(std::span<const Point> points) noexcept {
    Vec n{};
    const std::size_t k = points.size();
    for (std::size_t i = 0; i < k; ++i) n = add(n, cross(points[i], points[(i + 1) % k]));
    return 0.5 * norm(n);
}

double enclosed_volume(std::span<const Polygon> triangles) noexcept {
    double v = 0.0;
    for (const Polygon& t : triangles) {
        for (std::size_t i = 1; i + 1 < t.vertices.size(); ++i) {
            v += dot(t.vertices[0].position, cross(t.vertices[i].position, t.vertices[i + 1].position));
        }
    }
    return v / 6.0;
}

}  // namespace NiTCAD::render
