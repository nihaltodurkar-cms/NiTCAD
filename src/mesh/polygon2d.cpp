#include "polygon2d.hpp"

#include <algorithm>

#include "predicates.hpp"

namespace NiTCAD::mesh::detail {

using predicates::orient2d;

Point project(const Point& p, int drop) {
    return {p[static_cast<std::size_t>((drop + 1) % 3)], p[static_cast<std::size_t>((drop + 2) % 3)], 0.0};
}

bool collinear(const Point& a, const Point& b, const Point& c) {
    for (int k = 0; k < 3; ++k) {
        if (orient2d(project(a, k), project(b, k), project(c, k)) != 0) return false;
    }
    return true;
}

bool in_triangle_2d(const Point& a, const Point& b, const Point& c, const Point& x) {
    const int o1 = orient2d(a, b, x);
    const int o2 = orient2d(b, c, x);
    const int o3 = orient2d(c, a, x);
    return (o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0);
}

bool on_segment_2d(const Point& a, const Point& b, const Point& x) {
    if (orient2d(a, b, x) != 0) return false;
    for (std::size_t k = 0; k < 2; ++k) {
        if (x[k] < std::min(a[k], b[k]) || x[k] > std::max(a[k], b[k])) return false;
    }
    return true;
}

bool cross_properly_2d(const Point& a, const Point& b, const Point& c, const Point& d) {
    const int o1 = orient2d(a, b, c);
    const int o2 = orient2d(a, b, d);
    const int o3 = orient2d(c, d, a);
    const int o4 = orient2d(c, d, b);
    return o1 * o2 < 0 && o3 * o4 < 0;
}

bool segments_meet_2d(const Point& a, const Point& b, const Point& c, const Point& d) {
    return cross_properly_2d(a, b, c, d) || on_segment_2d(a, b, c) || on_segment_2d(a, b, d) ||
           on_segment_2d(c, d, a) || on_segment_2d(c, d, b);
}

int polygon_turn(const std::vector<Point>& q) {
    std::size_t m = 0;
    for (std::size_t i = 1; i < q.size(); ++i) {
        if (q[i][0] < q[m][0] || (q[i][0] == q[m][0] && q[i][1] < q[m][1])) m = i;
    }
    const std::size_t n = q.size();
    return orient2d(q[(m + n - 1) % n], q[m], q[(m + 1) % n]);
}

bool simple_2d(const std::vector<Point>& q) {
    const std::size_t n = q.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = q[i];
        const Point& b = q[(i + 1) % n];
        const Point& c = q[(i + 2) % n];
        if (orient2d(a, b, c) == 0 && predicates::dot_sign(b, a, c) >= 0) return false;  // a reversal
        for (std::size_t j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) continue;  // consecutive through vertex 0
            if (segments_meet_2d(a, b, q[j], q[(j + 1) % n])) return false;
        }
    }
    return true;
}

bool strictly_inside_2d(const std::vector<Point>& q, const Point& p) {
    const std::size_t n = q.size();
    bool in = false;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = q[i];
        const Point& b = q[(i + 1) % n];
        if (on_segment_2d(a, b, p)) return false;
        if ((a[1] > p[1]) != (b[1] > p[1])) {
            const int o = orient2d(a, b, p);
            if (b[1] > a[1] ? o > 0 : o < 0) in = !in;
        }
    }
    return in;
}

bool inside_closed_2d(const std::vector<Point>& q, const Point& p) {
    for (std::size_t i = 0; i < q.size(); ++i) {
        if (on_segment_2d(q[i], q[(i + 1) % q.size()], p)) return true;
    }
    return strictly_inside_2d(q, p);
}

std::optional<std::vector<std::array<std::size_t, 3>>> ear_clip(const std::vector<Point>& q, int turn) {
    std::vector<std::size_t> left(q.size());
    for (std::size_t i = 0; i < q.size(); ++i) left[i] = i;
    std::vector<std::array<std::size_t, 3>> out;
    while (left.size() > 3) {
        bool found = false;
        const std::size_t n = left.size();
        for (std::size_t k = 0; k < n && !found; ++k) {
            const std::size_t a = left[(k + n - 1) % n];
            const std::size_t b = left[k];
            const std::size_t c = left[(k + 1) % n];
            if (orient2d(q[a], q[b], q[c]) != turn) continue;
            bool empty = true;
            for (const std::size_t v : left) {
                if (v != a && v != b && v != c && in_triangle_2d(q[a], q[b], q[c], q[v])) {
                    empty = false;
                    break;
                }
            }
            if (!empty) continue;
            out.push_back({a, b, c});
            left.erase(left.begin() + static_cast<std::ptrdiff_t>(k));
            found = true;
        }
        if (!found) return std::nullopt;
    }
    if (orient2d(q[left[0]], q[left[1]], q[left[2]]) != turn) return std::nullopt;
    out.push_back({left[0], left[1], left[2]});
    return out;
}

}  // namespace NiTCAD::mesh::detail
