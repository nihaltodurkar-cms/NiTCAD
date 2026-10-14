#include "triangulation2d.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "NiTCAD/base/contract.hpp"
#include "predicates.hpp"

namespace NiTCAD::mesh::detail {

Triangulation2D::Triangulation2D(const Point& lo, const Point& hi) {
    const double cx = 0.5 * (lo[0] + hi[0]);
    const double cy = 0.5 * (lo[1] + hi[1]);
    const double extent = std::max({hi[0] - lo[0], hi[1] - lo[1], std::abs(cx) + std::abs(cy), 1e-300});
    const double d = 64.0 * extent;
    points_ = {{cx, cy + 2.0 * d, 0.0}, {cx - 1.7320508075688772 * d, cy - d, 0.0},
               {cx + 1.7320508075688772 * d, cy - d, 0.0}};
    vertex_triangle_ = {0, 0, 0};
    tris_.push_back(Triangle{{0, 1, 2}, {-1, -1, -1}, true});
    last_ = 0;
}

NodeId Triangulation2D::tie_id(int v) const noexcept {
    return v < 3 ? std::numeric_limits<NodeId>::max() - v : static_cast<NodeId>(v - 3);
}

bool Triangulation2D::in_circle(int t, const Point& p, int pid) const {
    const Triangle& tr = tris_[static_cast<std::size_t>(t)];
    return predicates::incircle_perturbed(
               {point(tr.v[0]), point(tr.v[1]), point(tr.v[2]), p},
               {tie_id(tr.v[0]), tie_id(tr.v[1]), tie_id(tr.v[2]), tie_id(pid)}) > 0;
}

int Triangulation2D::locate(const Point& p, int start) const {
    int t = start;
    if (t < 0 || static_cast<std::size_t>(t) >= tris_.size() || !tris_[static_cast<std::size_t>(t)].alive) t = last_;
    const std::size_t cap = 4 * tris_.size() + 16;
    for (std::size_t step = 0; step < cap; ++step) {
        const Triangle& tr = tris_[static_cast<std::size_t>(t)];
        bool moved = false;
        for (int i = 0; i < 3 && !moved; ++i) {
            const Point& a = point(tr.v[static_cast<std::size_t>((i + 1) % 3)]);
            const Point& b = point(tr.v[static_cast<std::size_t>((i + 2) % 3)]);
            if (predicates::orient2d(a, b, p) < 0) {
                const int next = tr.n[static_cast<std::size_t>(i)];
                NITCAD_EXPECTS(next >= 0);  // p is inside the super triangle
                t = next;
                moved = true;
            }
        }
        if (!moved) return t;
    }
    // The walk did not settle (it should on a regular triangulation): scan.
    for (std::size_t k = 0; k < tris_.size(); ++k) {
        const Triangle& tr = tris_[k];
        if (!tr.alive) continue;
        bool inside = true;
        for (int i = 0; i < 3 && inside; ++i) {
            inside = predicates::orient2d(point(tr.v[static_cast<std::size_t>((i + 1) % 3)]),
                                          point(tr.v[static_cast<std::size_t>((i + 2) % 3)]), p) >= 0;
        }
        if (inside) return static_cast<int>(k);
    }
    NITCAD_EXPECTS(false);
    return -1;
}

Triangulation2D::Cavity Triangulation2D::cavity(const Point& p, int containing) const {
    const int pid = static_cast<int>(points_.size());
    Cavity c;
    std::vector<int> stack{containing};
    std::vector<int> seen{containing};
    NITCAD_EXPECTS(in_circle(containing, p, pid));
    c.triangles.push_back(containing);
    while (!stack.empty()) {
        const int t = stack.back();
        stack.pop_back();
        for (const int nb : tris_[static_cast<std::size_t>(t)].n) {
            if (nb < 0 || std::ranges::find(seen, nb) != seen.end()) continue;
            seen.push_back(nb);
            if (in_circle(nb, p, pid)) {
                c.triangles.push_back(nb);
                stack.push_back(nb);
            }
        }
    }
    for (const int t : c.triangles) {
        const Triangle& tr = tris_[static_cast<std::size_t>(t)];
        for (std::size_t i = 0; i < 3; ++i) {
            const int nb = tr.n[i];
            if (nb >= 0 && std::ranges::find(c.triangles, nb) != c.triangles.end()) continue;
            c.boundary.push_back({tr.v[(i + 1) % 3], tr.v[(i + 2) % 3], nb});
        }
    }
    return c;
}

int Triangulation2D::new_triangle(const std::array<int, 3>& v) {
    int t = 0;
    if (!free_.empty()) {
        t = free_.back();
        free_.pop_back();
        tris_[static_cast<std::size_t>(t)] = Triangle{v, {-1, -1, -1}, true};
    } else {
        t = static_cast<int>(tris_.size());
        tris_.push_back(Triangle{v, {-1, -1, -1}, true});
    }
    return t;
}

int Triangulation2D::insert(const Point& p, const Cavity& c) {
    const int id = static_cast<int>(points_.size());
    points_.push_back(p);
    vertex_triangle_.push_back(-1);
    for (const int t : c.triangles) {
        tris_[static_cast<std::size_t>(t)].alive = false;
        free_.push_back(t);
    }
    std::vector<std::pair<int, int>> by_first;   // (a, new triangle (a, b, p))
    std::vector<std::pair<int, int>> by_second;  // (b, new triangle (a, b, p))
    for (const Cavity::Edge& e : c.boundary) {
        const int t = new_triangle({e.a, e.b, id});
        tris_[static_cast<std::size_t>(t)].n[2] = e.outside;
        if (e.outside >= 0) {
            Triangle& out = tris_[static_cast<std::size_t>(e.outside)];
            for (std::size_t j = 0; j < 3; ++j) {
                if (out.v[(j + 1) % 3] == e.b && out.v[(j + 2) % 3] == e.a) out.n[j] = t;
            }
        }
        by_first.emplace_back(e.a, t);
        by_second.emplace_back(e.b, t);
        vertex_triangle_[static_cast<std::size_t>(e.a)] = t;
        vertex_triangle_[static_cast<std::size_t>(e.b)] = t;
        vertex_triangle_[static_cast<std::size_t>(id)] = t;
        last_ = t;
    }
    for (const auto& [a, t] : by_first) {
        Triangle& tr = tris_[static_cast<std::size_t>(t)];
        const int b = tr.v[1];
        // Opposite a: the edge b -> p, shared with the new triangle (b, c, p).
        for (const auto& [x, u] : by_first) {
            if (x == b) tr.n[0] = u;
        }
        // Opposite b: the edge p -> a, shared with the new triangle (z, a, p).
        for (const auto& [y, u] : by_second) {
            if (y == a) tr.n[1] = u;
        }
    }
    return id;
}

int Triangulation2D::insert(const Point& p) {
    const int t = locate(p, last_);
    return insert(p, cavity(p, t));
}

std::optional<int> Triangulation2D::triangle_with_edge(int a, int b) const {
    const int start = vertex_triangle_[static_cast<std::size_t>(a)];
    if (start < 0) return std::nullopt;
    // Rotate around a both ways from start.
    for (const int direction : {0, 1}) {
        int t = start;
        for (std::size_t guard = 0; guard < tris_.size() + 1; ++guard) {
            const Triangle& tr = tris_[static_cast<std::size_t>(t)];
            std::size_t i = 0;
            while (tr.v[i] != a) ++i;
            if (tr.v[(i + 1) % 3] == b) return t;
            // Across the edge (a, v[i+1]) (opposite v[i+2]) or the edge (v[i+2], a) (opposite v[i+1]).
            const int next = direction == 0 ? tr.n[(i + 2) % 3] : tr.n[(i + 1) % 3];
            if (next < 0 || next == start) break;
            t = next;
        }
    }
    return std::nullopt;
}

}  // namespace NiTCAD::mesh::detail
