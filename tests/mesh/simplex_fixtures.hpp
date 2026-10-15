// Shared fixtures for the simplex-mesh tests (Unit 16): grid points, box boundary tagging, a
// brute-force Delaunay triangulation under the perturbed predicates, and small helpers.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"
#include "predicates.hpp"

namespace mesh_fixtures {

using NiTCAD::mesh::NodeId;
using NiTCAD::mesh::Point;
using NiTCAD::mesh::SimplexMesh;

inline bool close(double a, double b, double rel) {
    return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}

struct SplitMix {
    std::uint64_t state;
    double uniform() {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return static_cast<double>((z ^ (z >> 31)) >> 11) * 0x1p-53;
    }
};

// Grid points in make_tensor_grid's order (x fastest).
inline std::vector<Point> grid_points(const std::vector<double>& x, const std::vector<double>& y,
                               const std::vector<double>& z = {}) {
    std::vector<Point> p;
    const std::size_t nz = z.empty() ? 1 : z.size();
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < y.size(); ++j) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                p.push_back({x[i], y[j], z.empty() ? 0.0 : z[k]});
            }
        }
    }
    return p;
}

// Lists every topological boundary facet whose nodes lie on one face of the box [lo, hi], under
// the tensor grid's patch names and order: x_min, x_max, y_min, y_max, z_min, z_max.
inline void tag_box_boundary(SimplexMesh& s, const Point& lo, const Point& hi) {
    const auto dim = static_cast<std::size_t>(s.dimension);
    const std::size_t per = dim + 1;
    std::map<std::vector<NodeId>, int> count;
    for (std::size_t c = 0; c < s.cells.size() / per; ++c) {
        for (std::size_t l = 0; l < per; ++l) {
            std::vector<NodeId> key;
            for (std::size_t q = 0; q < per; ++q) {
                if (q != l) key.push_back(s.cells[c * per + q]);
            }
            std::ranges::sort(key);
            ++count[key];
        }
    }
    const std::array<const char*, 6> names{"x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
    s.patch_names.assign(names.begin(), names.begin() + static_cast<std::ptrdiff_t>(2 * dim));
    for (const auto& [key, c] : count) {
        if (c != 1) continue;
        for (std::size_t a = 0; a < dim; ++a) {
            for (const int side : {0, 1}) {
                const double plane = side == 0 ? lo[a] : hi[a];
                const bool on = std::ranges::all_of(
                    key, [&](NodeId v) { return s.points[static_cast<std::size_t>(v)][a] == plane; });
                if (!on) continue;
                s.facets.insert(s.facets.end(), key.begin(), key.end());
                s.facet_patch.push_back(static_cast<std::int32_t>(2 * a) + side);
            }
        }
    }
}


// The Delaunay triangulation of the points under the perturbed predicates, by brute force: every
// non-degenerate simplex whose perturbed circumsphere holds no other point.
inline std::vector<NodeId> brute_delaunay(int dim, const std::vector<Point>& p) {
    const auto n = static_cast<NodeId>(p.size());
    const auto P = [&](NodeId i) -> const Point& { return p[static_cast<std::size_t>(i)]; };
    std::vector<NodeId> cells;
    if (dim == 2) {
        for (NodeId i = 0; i < n; ++i) {
            for (NodeId j = i + 1; j < n; ++j) {
                for (NodeId k = j + 1; k < n; ++k) {
                    const int o = NiTCAD::mesh::predicates::orient2d(P(i), P(j), P(k));
                    if (o == 0) continue;
                    const NodeId a = o > 0 ? i : j, b = o > 0 ? j : i;
                    bool empty = true;
                    for (NodeId l = 0; l < n && empty; ++l) {
                        if (l == i || l == j || l == k) continue;
                        empty = NiTCAD::mesh::predicates::incircle_perturbed({P(a), P(b), P(k), P(l)}, {a, b, k, l}) < 0;
                    }
                    if (empty) cells.insert(cells.end(), {a, b, k});
                }
            }
        }
        return cells;
    }
    for (NodeId i = 0; i < n; ++i) {
        for (NodeId j = i + 1; j < n; ++j) {
            for (NodeId k = j + 1; k < n; ++k) {
                for (NodeId m = k + 1; m < n; ++m) {
                    const int o = NiTCAD::mesh::predicates::orient3d(P(i), P(j), P(k), P(m));
                    if (o == 0) continue;
                    const NodeId a = o > 0 ? i : j, b = o > 0 ? j : i;
                    bool empty = true;
                    for (NodeId l = 0; l < n && empty; ++l) {
                        if (l == i || l == j || l == k || l == m) continue;
                        empty = NiTCAD::mesh::predicates::insphere_perturbed({P(a), P(b), P(k), P(m), P(l)},
                                                       {a, b, k, m, l}) < 0;
                    }
                    if (empty) cells.insert(cells.end(), {a, b, k, m});
                }
            }
        }
    }
    return cells;
}

// The total measure of the cells (positive orientation assumed), in double.
inline double cell_measure(int dim, const SimplexMesh& s) {
    const std::size_t per = static_cast<std::size_t>(dim) + 1;
    double total = 0.0;
    for (std::size_t c = 0; c < s.cells.size() / per; ++c) {
        const auto P = [&](std::size_t k) { return s.points[static_cast<std::size_t>(s.cells[c * per + k])]; };
        const Point a = P(0), b = P(1), d = P(2);
        const std::array<double, 3> u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const std::array<double, 3> v{d[0] - a[0], d[1] - a[1], d[2] - a[2]};
        if (dim == 2) {
            total += 0.5 * std::abs(u[0] * v[1] - u[1] * v[0]);
        } else {
            const Point e = P(3);
            const std::array<double, 3> w{e[0] - a[0], e[1] - a[1], e[2] - a[2]};
            total += std::abs(u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) +
                              u[2] * (v[0] * w[1] - v[1] * w[0])) /
                     6.0;
        }
    }
    return total;
}

inline SimplexMesh two_d(std::vector<Point> points, std::vector<NodeId> cells) {
    SimplexMesh s;
    s.dimension = 2;
    s.points = std::move(points);
    s.cells = std::move(cells);
    return s;
}

}  // namespace mesh_fixtures
