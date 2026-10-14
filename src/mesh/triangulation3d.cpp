#include "triangulation3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <tuple>

#include "NiTCAD/base/contract.hpp"
#include "predicates.hpp"

namespace NiTCAD::mesh::detail {

Triangulation3D::Triangulation3D(const Point& lo, const Point& hi) {
    const double cx = 0.5 * (lo[0] + hi[0]);
    const double cy = 0.5 * (lo[1] + hi[1]);
    const double cz = 0.5 * (lo[2] + hi[2]);
    const double extent = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2],
                                    std::abs(cx) + std::abs(cy) + std::abs(cz), 1e-300});
    // A regular tetrahedron on alternate corners of a cube of half-side d: its faces lie d / sqrt 3
    // from the centre, far outside the box.
    const double d = 64.0 * extent;
    points_ = {{cx + d, cy + d, cz + d}, {cx - d, cy + d, cz - d}, {cx + d, cy - d, cz - d}, {cx - d, cy - d, cz + d}};
    NITCAD_EXPECTS(predicates::orient3d(points_[0], points_[1], points_[2], points_[3]) > 0);
    vertex_tet_ = {0, 0, 0, 0};
    tets_.push_back(Tetrahedron{{0, 1, 2, 3}, {-1, -1, -1, -1}, true});
    mark_.push_back(0);
    last_ = 0;
}

NodeId Triangulation3D::tie_id(int v) const noexcept {
    return v < 4 ? std::numeric_limits<NodeId>::max() - v : static_cast<NodeId>(v - 4);
}

bool Triangulation3D::in_sphere(int t, const Point& p, int pid) const {
    const Tetrahedron& tt = tets_[static_cast<std::size_t>(t)];
    return predicates::insphere_perturbed(
               {point(tt.v[0]), point(tt.v[1]), point(tt.v[2]), point(tt.v[3]), p},
               {tie_id(tt.v[0]), tie_id(tt.v[1]), tie_id(tt.v[2]), tie_id(tt.v[3]), tie_id(pid)}) > 0;
}

namespace {

// orient3d of tetrahedron v with v[i] replaced by p: negative when p is beyond face i.
int side(const std::array<Point, 4>& q, std::size_t i, const Point& p) {
    std::array<Point, 4> r = q;
    r[i] = p;
    return predicates::orient3d(r[0], r[1], r[2], r[3]);
}

}  // namespace

int Triangulation3D::locate(const Point& p, int start) const {
    int t = start;
    if (t < 0 || static_cast<std::size_t>(t) >= tets_.size() || !tets_[static_cast<std::size_t>(t)].alive) t = last_;
    const std::size_t cap = 4 * tets_.size() + 16;
    for (std::size_t step = 0; step < cap; ++step) {
        const Tetrahedron& tt = tets_[static_cast<std::size_t>(t)];
        const std::array<Point, 4> q{point(tt.v[0]), point(tt.v[1]), point(tt.v[2]), point(tt.v[3])};
        bool moved = false;
        for (std::size_t i = 0; i < 4 && !moved; ++i) {
            if (side(q, i, p) < 0) {
                const int next = tt.n[i];
                NITCAD_EXPECTS(next >= 0);  // p is inside the super tetrahedron
                t = next;
                moved = true;
            }
        }
        if (!moved) return t;
    }
    // The walk did not settle (it should on a regular tetrahedralization): scan.
    for (std::size_t k = 0; k < tets_.size(); ++k) {
        const Tetrahedron& tt = tets_[k];
        if (!tt.alive) continue;
        const std::array<Point, 4> q{point(tt.v[0]), point(tt.v[1]), point(tt.v[2]), point(tt.v[3])};
        bool inside = true;
        for (std::size_t i = 0; i < 4 && inside; ++i) inside = side(q, i, p) >= 0;
        if (inside) return static_cast<int>(k);
    }
    NITCAD_EXPECTS(false);
    return -1;
}

Triangulation3D::Cavity Triangulation3D::cavity(const Point& p, int containing) const {
    const int pid = static_cast<int>(points_.size());
    if (++stamp_ == 0) {  // wrapped: clear the stamps
        std::ranges::fill(mark_, 0u);
        stamp_ = 1;
    }
    // While searching, mark_ == stamp_ means visited (tested once).
    Cavity c;
    NITCAD_EXPECTS(in_sphere(containing, p, pid));
    c.tetrahedra.push_back(containing);
    const auto mark = [&](int t) -> std::uint32_t& { return mark_[static_cast<std::size_t>(t)]; };
    std::vector<int> stack{containing};
    std::vector<int> outside_seen;
    mark(containing) = stamp_;
    while (!stack.empty()) {
        const int t = stack.back();
        stack.pop_back();
        for (const int nb : tets_[static_cast<std::size_t>(t)].n) {
            if (nb < 0 || mark(nb) == stamp_) continue;
            mark(nb) = stamp_;
            if (in_sphere(nb, p, pid)) {
                c.tetrahedra.push_back(nb);
                stack.push_back(nb);
            } else {
                outside_seen.push_back(nb);
            }
        }
    }
    // Re-mark: cavity members keep stamp_, the outside neighbours seen get stamp_ - 1 so that the
    // boundary scan below can tell them apart. (Stamps only grow, so stamp_ - 1 is never a later
    // search's stamp.)
    for (const int t : outside_seen) mark(t) = stamp_ - 1;
    for (const int t : c.tetrahedra) {
        const Tetrahedron& tt = tets_[static_cast<std::size_t>(t)];
        for (std::size_t i = 0; i < 4; ++i) {
            const int nb = tt.n[i];
            if (nb >= 0 && mark(nb) == stamp_) continue;
            c.boundary.push_back({t, static_cast<int>(i), nb});
        }
    }
    return c;
}

int Triangulation3D::new_tet(const std::array<int, 4>& v) {
    int t = 0;
    if (!free_.empty()) {
        t = free_.back();
        free_.pop_back();
        tets_[static_cast<std::size_t>(t)] = Tetrahedron{v, {-1, -1, -1, -1}, true};
    } else {
        t = static_cast<int>(tets_.size());
        tets_.push_back(Tetrahedron{v, {-1, -1, -1, -1}, true});
        mark_.push_back(0);
    }
    return t;
}

int Triangulation3D::insert(const Point& p, const Cavity& c) {
    const int id = static_cast<int>(points_.size());
    // Read the new tetrahedra off the cavity before its slots are reused.
    struct Pending {
        std::array<int, 4> v;
        int slot;
        int outside;
        int old;
    };
    std::vector<Pending> pending;
    pending.reserve(c.boundary.size());
    for (const Cavity::Face& f : c.boundary) {
        std::array<int, 4> v = tets_[static_cast<std::size_t>(f.tet)].v;
        v[static_cast<std::size_t>(f.slot)] = id;
        pending.push_back({v, f.slot, f.outside, f.tet});
    }
    points_.push_back(p);
    vertex_tet_.push_back(-1);
    // The cavity's slots are freed only after the outside neighbours are relinked below: they
    // find their old neighbour by slot, so a slot must not be reused while that link is pending.
    for (const int t : c.tetrahedra) tets_[static_cast<std::size_t>(t)].alive = false;
    // The faces among the new tetrahedra: the face opposite v[j] (j != slot) holds p and two old
    // vertices; each such pair of vertices is shared by exactly two new tetrahedra.
    std::vector<std::tuple<int, int, int, int>> shared;  // (lo, hi, new tet, j)
    shared.reserve(3 * pending.size());
    for (const Pending& q : pending) {
        const int t = new_tet(q.v);
        Tetrahedron& nt = tets_[static_cast<std::size_t>(t)];
#ifndef NDEBUG
        NITCAD_EXPECTS(predicates::orient3d(point(nt.v[0]), point(nt.v[1]), point(nt.v[2]), point(nt.v[3])) > 0);
#endif
        nt.n[static_cast<std::size_t>(q.slot)] = q.outside;
        if (q.outside >= 0) {
            Tetrahedron& out = tets_[static_cast<std::size_t>(q.outside)];
            for (std::size_t j = 0; j < 4; ++j) {
                if (out.n[j] == q.old) out.n[j] = t;
            }
        }
        for (std::size_t j = 0; j < 4; ++j) {
            if (static_cast<int>(j) == q.slot) continue;
            int a = -1;
            int b = -1;
            for (std::size_t k = 0; k < 4; ++k) {
                if (k == j || static_cast<int>(k) == q.slot) continue;
                (a < 0 ? a : b) = nt.v[k];
            }
            shared.emplace_back(std::min(a, b), std::max(a, b), t, static_cast<int>(j));
        }
        for (const int v : nt.v) vertex_tet_[static_cast<std::size_t>(v)] = t;
        last_ = t;
    }
    std::ranges::sort(shared);
    for (std::size_t k = 0; k + 1 < shared.size(); k += 2) {
        const auto& [a0, b0, t0, j0] = shared[k];
        const auto& [a1, b1, t1, j1] = shared[k + 1];
        NITCAD_EXPECTS(a0 == a1 && b0 == b1);
        tets_[static_cast<std::size_t>(t0)].n[static_cast<std::size_t>(j0)] = t1;
        tets_[static_cast<std::size_t>(t1)].n[static_cast<std::size_t>(j1)] = t0;
    }
    for (const int t : c.tetrahedra) free_.push_back(t);
    return id;
}

int Triangulation3D::insert(const Point& p) {
    const int t = locate(p, last_);
    return insert(p, cavity(p, t));
}

std::vector<int> Triangulation3D::incident(int v) const {
    std::vector<int> out;
    const int start = vertex_tet_[static_cast<std::size_t>(v)];
    if (start < 0) return out;
    if (++stamp_ == 0) {
        std::ranges::fill(mark_, 0u);
        stamp_ = 1;
    }
    std::vector<int> stack{start};
    mark_[static_cast<std::size_t>(start)] = stamp_;
    while (!stack.empty()) {
        const int t = stack.back();
        stack.pop_back();
        out.push_back(t);
        const Tetrahedron& tt = tets_[static_cast<std::size_t>(t)];
        for (std::size_t i = 0; i < 4; ++i) {
            if (tt.v[i] == v) continue;  // the face opposite v does not hold v
            const int nb = tt.n[i];
            if (nb < 0 || mark_[static_cast<std::size_t>(nb)] == stamp_) continue;
            mark_[static_cast<std::size_t>(nb)] = stamp_;
            stack.push_back(nb);
        }
    }
    std::ranges::sort(out);
    return out;
}

namespace {

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

// The Hilbert index of a cell (Skilling 2004, "Programming the Hilbert curve": axes to transpose,
// then the transposed bits interleaved), 16 bits per axis.
std::uint64_t hilbert3(std::array<std::uint32_t, 3> x) {
    constexpr std::uint32_t m = 1u << 15;
    for (std::uint32_t q = m; q > 1; q >>= 1) {
        const std::uint32_t pm = q - 1;
        for (std::size_t i = 0; i < 3; ++i) {
            if ((x[i] & q) != 0) {
                x[0] ^= pm;
            } else {
                const std::uint32_t t = (x[0] ^ x[i]) & pm;
                x[0] ^= t;
                x[i] ^= t;
            }
        }
    }
    for (std::size_t i = 1; i < 3; ++i) x[i] ^= x[i - 1];
    std::uint32_t t = 0;
    for (std::uint32_t q = m; q > 1; q >>= 1) {
        if ((x[2] & q) != 0) t ^= q - 1;
    }
    for (std::size_t i = 0; i < 3; ++i) x[i] ^= t;
    std::uint64_t key = 0;
    for (int bit = 15; bit >= 0; --bit) {
        for (std::size_t i = 0; i < 3; ++i) key = (key << 1) | ((x[i] >> bit) & 1u);
    }
    return key;
}

}  // namespace

std::vector<std::size_t> brio_order(const std::vector<Point>& p) {
    if (p.empty()) return {};
    Point lo = p[0];
    Point hi = p[0];
    for (const Point& q : p) {
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], q[k]);
            hi[k] = std::max(hi[k], q[k]);
        }
    }
    const double span = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-300});
    constexpr std::uint32_t n = 1u << 16;
    // Round: the number of leading one bits of the hash, so half the points are in the last round,
    // a quarter in the one before, ...; at most 24 rounds.
    std::vector<std::tuple<int, std::uint64_t, std::size_t>> key(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) {
        std::array<std::uint32_t, 3> c{};
        for (std::size_t k = 0; k < 3; ++k) {
            c[k] = static_cast<std::uint32_t>(std::min((p[i][k] - lo[k]) / span * (n - 1), double(n - 1)));
        }
        const std::uint64_t h = splitmix64(0x4E695443414Dull ^ static_cast<std::uint64_t>(i));
        int level = 0;
        while (level < 24 && ((h >> (63 - level)) & 1u) != 0) ++level;
        key[i] = {-level, hilbert3(c), i};
    }
    std::ranges::sort(key);
    std::vector<std::size_t> order(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) order[i] = std::get<2>(key[i]);
    return order;
}

}  // namespace NiTCAD::mesh::detail
