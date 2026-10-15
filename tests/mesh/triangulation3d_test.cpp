// The 3D Delaunay tetrahedralization of the generator (Unit 16, commit 6a): against the brute-force
// Delaunay of the perturbed predicates, its invariants, determinism and the insertion order.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <set>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "dual_geometry.hpp"
#include "predicates.hpp"
#include "simplex_fixtures.hpp"
#include "triangulation3d.hpp"

using NiTCAD::mesh::NodeId;
using NiTCAD::mesh::Point;
using NiTCAD::mesh::detail::brio_order;
using NiTCAD::mesh::detail::Triangulation3D;
namespace predicates = NiTCAD::mesh::predicates;

namespace {

Point lower(const std::vector<Point>& p) {
    Point lo = p[0];
    for (const Point& q : p) {
        for (std::size_t k = 0; k < 3; ++k) lo[k] = std::min(lo[k], q[k]);
    }
    return lo;
}

Point upper(const std::vector<Point>& p) {
    Point hi = p[0];
    for (const Point& q : p) {
        for (std::size_t k = 0; k < 3; ++k) hi[k] = std::max(hi[k], q[k]);
    }
    return hi;
}

// Inserts p in index order, so the tie-break ids are the indices.
Triangulation3D build(const std::vector<Point>& p) {
    Triangulation3D t(lower(p), upper(p));
    for (const Point& q : p) t.insert(q);
    return t;
}

// The tetrahedra with no super vertex, as sorted point indices.
std::set<std::array<NodeId, 4>> real_tets(const Triangulation3D& t) {
    std::set<std::array<NodeId, 4>> out;
    for (std::size_t k = 0; k < t.tet_slots(); ++k) {
        const auto& tt = t.tet(static_cast<int>(k));
        if (!tt.alive) continue;
        if (std::ranges::any_of(tt.v, [&](int v) { return t.is_super(v); })) continue;
        std::array<NodeId, 4> a{};
        for (std::size_t i = 0; i < 4; ++i) a[i] = t.tie_id(tt.v[i]);
        std::ranges::sort(a);
        out.insert(a);
    }
    return out;
}

std::set<std::array<NodeId, 4>> brute(const std::vector<Point>& p) {
    const std::vector<NodeId> cells = mesh_fixtures::brute_delaunay(3, p);
    std::set<std::array<NodeId, 4>> out;
    for (std::size_t c = 0; c < cells.size(); c += 4) {
        std::array<NodeId, 4> a{cells[c], cells[c + 1], cells[c + 2], cells[c + 3]};
        std::ranges::sort(a);
        out.insert(a);
    }
    return out;
}

// Orientation, neighbour symmetry (each shared face has the same three vertices) and the empty
// perturbed circumsphere of every alive tetrahedron against every vertex.
void check_invariants(const Triangulation3D& t) {
    std::size_t alive = 0;
    for (std::size_t k = 0; k < t.tet_slots(); ++k) {
        const auto& tt = t.tet(static_cast<int>(k));
        if (!tt.alive) continue;
        ++alive;
        REQUIRE(predicates::orient3d(t.point(tt.v[0]), t.point(tt.v[1]), t.point(tt.v[2]), t.point(tt.v[3])) > 0);
        for (std::size_t i = 0; i < 4; ++i) {
            const int nb = tt.n[i];
            if (nb < 0) continue;
            const auto& o = t.tet(nb);
            REQUIRE(o.alive);
            std::size_t back = 4;
            for (std::size_t j = 0; j < 4; ++j) back = o.n[j] == static_cast<int>(k) ? j : back;
            REQUIRE(back < 4);
            std::array<int, 3> f{};
            std::array<int, 3> g{};
            for (std::size_t a = 0, m = 0, l = 0; a < 4; ++a) {
                if (a != i) f[m++] = tt.v[a];
                if (a != back) g[l++] = o.v[a];
            }
            std::ranges::sort(f);
            std::ranges::sort(g);
            REQUIRE(f == g);
        }
        for (std::size_t v = 0; v < t.vertex_count(); ++v) {
            const int vi = static_cast<int>(v);
            if (std::ranges::find(tt.v, vi) != tt.v.end()) continue;
            REQUIRE(predicates::insphere_perturbed(
                        {t.point(tt.v[0]), t.point(tt.v[1]), t.point(tt.v[2]), t.point(tt.v[3]), t.point(vi)},
                        {t.tie_id(tt.v[0]), t.tie_id(tt.v[1]), t.tie_id(tt.v[2]), t.tie_id(tt.v[3]), t.tie_id(vi)}) < 0);
        }
    }
    CHECK(alive > 0);
}

std::vector<Point> random_points(std::size_t n, std::uint64_t seed) {
    mesh_fixtures::SplitMix r{seed};
    std::vector<Point> p(n);
    for (Point& q : p) q = {r.uniform(), r.uniform(), r.uniform()};
    return p;
}

}  // namespace

TEST_CASE("3D Delaunay: a lattice with every tie equals the brute-force perturbed Delaunay", "[mesh][delaunay3d]") {
    // A 3 x 3 x 3 lattice: every cube is cospherical, so every cell is decided by the tie-break.
    const std::vector<Point> p = mesh_fixtures::grid_points({0, 0.5, 1}, {0, 0.25, 1}, {0, 0.5, 1.5});
    const Triangulation3D t = build(p);
    check_invariants(t);
    // The hull's tetrahedra are well shaped, so none reaches a super vertex: equal sets.
    CHECK(real_tets(t) == brute(p));
}

TEST_CASE("3D Delaunay: random points give a subset of the brute-force Delaunay covering the well-shaped part",
          "[mesh][delaunay3d]") {
    for (const std::uint64_t seed : {1ull, 2ull, 3ull}) {
        const std::vector<Point> p = random_points(18, seed);
        const Triangulation3D t = build(p);
        check_invariants(t);
        const auto ours = real_tets(t);
        const auto all = brute(p);
        // Real tetrahedra Delaunay among all vertices are Delaunay among the real ones. A brute
        // tetrahedron may be missing only when its circumsphere reaches a super vertex, which needs
        // a circumradius far beyond the points' extent.
        for (const auto& a : ours) CHECK(all.contains(a));
        for (const auto& a : all) {
            if (ours.contains(a)) continue;
            const auto& [i, j, k, l] = a;
            const auto P = [&](NodeId n) { return p[static_cast<std::size_t>(n)]; };
            const Point c = NiTCAD::mesh::detail::tetrahedron_circumcentre({P(i), P(j), P(k), P(l)});
            CHECK(NiTCAD::mesh::detail::distance(c, P(i)) > 4.0);
        }
    }
}

TEST_CASE("3D Delaunay: points on a plane and on a sphere (exact ties beyond cubes)", "[mesh][delaunay3d]") {
    // Eight cospherical points (a cube's corners rotated about z by a rational rotation, so the
    // coordinates are exact) plus the centre, and a planar 4 x 4 layer under them.
    std::vector<Point> p;
    for (const double x : {-0.6, 0.6}) {
        for (const double y : {-0.8, 0.8}) {
            for (const double z : {-1.0, 1.0}) p.push_back({x, y, z});
        }
    }
    p.push_back({0, 0, 0});
    for (const double x : {-1.0, -0.5, 0.5, 1.0}) {
        for (const double y : {-1.0, -0.25, 0.25, 1.0}) p.push_back({x, y, -2.0});
    }
    const Triangulation3D t = build(p);
    check_invariants(t);
    const auto all = brute(p);
    for (const auto& a : real_tets(t)) CHECK(all.contains(a));
}

TEST_CASE("3D Delaunay: the same input gives the same tetrahedra; the order changes only tie decisions",
          "[mesh][delaunay3d]") {
    const std::vector<Point> p = random_points(200, 7);
    const Triangulation3D a = build(p);
    const Triangulation3D b = build(p);
    REQUIRE(a.tet_slots() == b.tet_slots());
    for (std::size_t k = 0; k < a.tet_slots(); ++k) {
        CHECK(a.tet(static_cast<int>(k)).v == b.tet(static_cast<int>(k)).v);
        CHECK(a.tet(static_cast<int>(k)).n == b.tet(static_cast<int>(k)).n);
        CHECK(a.tet(static_cast<int>(k)).alive == b.tet(static_cast<int>(k)).alive);
    }
    // Random points have no ties, so another insertion order gives the same Delaunay set.
    const std::vector<std::size_t> order = brio_order(p);
    std::vector<Point> q;
    for (const std::size_t i : order) q.push_back(p[i]);
    const Triangulation3D c = build(q);
    check_invariants(c);
    std::set<std::array<NodeId, 4>> mapped;
    for (auto tet : real_tets(c)) {
        for (NodeId& v : tet) v = static_cast<NodeId>(order[static_cast<std::size_t>(v)]);
        std::ranges::sort(tet);
        mapped.insert(tet);
    }
    CHECK(mapped == real_tets(a));
}

TEST_CASE("3D Delaunay: incident tetrahedra and the BRIO order", "[mesh][delaunay3d]") {
    const std::vector<Point> p = random_points(200, 11);
    const Triangulation3D t = build(p);
    for (int v = 4; v < static_cast<int>(t.vertex_count()); v += 17) {
        std::vector<int> expected;
        for (std::size_t k = 0; k < t.tet_slots(); ++k) {
            const auto& tt = t.tet(static_cast<int>(k));
            if (tt.alive && std::ranges::find(tt.v, v) != tt.v.end()) expected.push_back(static_cast<int>(k));
        }
        CHECK(t.incident(v) == expected);
    }
    const std::vector<std::size_t> order = brio_order(p);
    std::vector<std::size_t> sorted = order;
    std::ranges::sort(sorted);
    for (std::size_t i = 0; i < sorted.size(); ++i) REQUIRE(sorted[i] == i);
    CHECK(brio_order(p) == order);
    CHECK(order != sorted);
}

TEST_CASE("3D Delaunay: a walk across many points", "[mesh][delaunay3d]") {
    const std::vector<Point> p = random_points(5000, 13);
    const std::vector<std::size_t> order = brio_order(p);
    Triangulation3D t(lower(p), upper(p));
    for (const std::size_t i : order) t.insert(p[i]);
    std::size_t alive = 0;
    for (std::size_t k = 0; k < t.tet_slots(); ++k) alive += t.tet(static_cast<int>(k)).alive ? 1 : 0;
    // About 6.7 tetrahedra per point for uniform random points, plus the hull's.
    CHECK(alive > 5 * p.size());
    CHECK(alive < 9 * p.size());
    for (std::size_t probe = 0; probe < p.size(); probe += 97) {
        const int at = t.locate(p[probe], 0);
        const auto& tt = t.tet(at);
        CHECK(std::ranges::find(tt.v, static_cast<int>(std::ranges::find(order, probe) - order.begin()) + 4) !=
              tt.v.end());
    }
}

TEST_CASE("3D Delaunay: speed (1e5 points in BRIO order)", "[.benchmark][mesh][delaunay3d]") {
    const std::vector<Point> p = random_points(100'000, 17);
    const auto start = std::chrono::steady_clock::now();
    const std::vector<std::size_t> order = brio_order(p);
    Triangulation3D t(lower(p), upper(p));
    for (const std::size_t i : order) t.insert(p[i]);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("3D Delaunay: %zu points in %.3f s\n", p.size(), seconds);
    CHECK(t.vertex_count() == p.size() + 4);
}
