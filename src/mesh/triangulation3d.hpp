// An incremental 3D Delaunay tetrahedralization (Bowyer-Watson) for the mesh generator (Unit 16;
// DECISIONS.md M3, M4, M13). Internal to the mesh layer.
//
// The 3D counterpart of Triangulation2D. Every topological decision is an exact predicate:
// orientation for point location, and the perturbed insphere test (the node-id tie-break) for the
// cavity, so ties resolve one consistent way and the result is the Delaunay tetrahedralization of
// the perturbed lifting. The points lie inside a large super tetrahedron whose four vertices are
// vertices 0-3; the caller's points are vertices 4, 5, ... in insertion order (deterministic: no
// hashing, no randomness). For the tie-break a super vertex ranks above every real one.
//
// Every tetrahedron is positively oriented: orient3d(v[0], v[1], v[2], v[3]) > 0. Face i is the
// face opposite v[i]. Point location walks from a given tetrahedron across faces the point is
// beyond (a visibility walk; it terminates on regular, hence Delaunay, tetrahedralizations; a step
// cap falls back to a scan).
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh::detail {

class Triangulation3D {
public:
    struct Tetrahedron {
        std::array<int, 4> v{};  // positively oriented
        std::array<int, 4> n{};  // n[i]: the neighbour across the face opposite v[i], or -1
        bool alive = false;
    };

    // The cavity of a point: the tetrahedra whose perturbed circumspheres contain it, and its
    // boundary faces, each as the cavity tetrahedron and the slot of the face (the face opposite
    // v[slot]), with the tetrahedron outside it (-1 on the super tetrahedron's hull).
    struct Cavity {
        std::vector<int> tetrahedra;
        struct Face {
            int tet;
            int slot;
            int outside;
        };
        std::vector<Face> boundary;
    };

    // A super tetrahedron enclosing the box [lo, hi] with a wide margin.
    Triangulation3D(const Point& lo, const Point& hi);

    [[nodiscard]] std::size_t vertex_count() const noexcept { return points_.size(); }
    [[nodiscard]] const Point& point(int v) const { return points_[static_cast<std::size_t>(v)]; }
    [[nodiscard]] bool is_super(int v) const noexcept { return v < 4; }
    [[nodiscard]] std::size_t tet_slots() const noexcept { return tets_.size(); }
    [[nodiscard]] const Tetrahedron& tet(int t) const { return tets_[static_cast<std::size_t>(t)]; }
    // A tetrahedron created most recently: a good start for a walk.
    [[nodiscard]] int last_tet() const noexcept { return last_; }
    // An alive tetrahedron incident to vertex v.
    [[nodiscard]] int vertex_tet(int v) const { return vertex_tet_[static_cast<std::size_t>(v)]; }

    // The alive tetrahedron containing p (on its boundary counts), walking from `start`.
    [[nodiscard]] int locate(const Point& p, int start) const;
    // The cavity of p, found from the tetrahedron containing it. Does not modify the mesh.
    [[nodiscard]] Cavity cavity(const Point& p, int containing) const;
    // Inserts p with this cavity (from cavity()); returns its vertex id.
    int insert(const Point& p, const Cavity& c);
    // locate, cavity and insert. Precondition: p strictly inside the super tetrahedron, not a vertex.
    int insert(const Point& p);

    // The tetrahedra created by the last insertion.
    [[nodiscard]] const std::vector<int>& created() const noexcept { return created_; }
    // The alive tetrahedra incident to vertex v.
    [[nodiscard]] std::vector<int> incident(int v) const;
    // The alive tetrahedra holding the triangle abc (at most two; -1 for none), found by walking the
    // tetrahedra around a and stopping at the second.
    [[nodiscard]] std::array<int, 2> tets_with_face(int a, int b, int c) const;
    // The alive tetrahedra holding the edge ab, into out (cleared first), in walk order.
    void tets_with_edge(int a, int b, std::vector<int>& out) const;
    // The id the tie-break uses for vertex v.
    [[nodiscard]] NodeId tie_id(int v) const noexcept;

private:
    [[nodiscard]] bool in_sphere(int t, const Point& p, int pid) const;
    int new_tet(const std::array<int, 4>& v);
    // Visits the alive tetrahedra around vertex v (a walk across faces holding v) until visit(t)
    // returns false.
    template <class Visit>
    void walk_around(int v, Visit&& visit) const;

    std::vector<Point> points_;
    std::vector<Tetrahedron> tets_;
    std::vector<int> free_;
    std::vector<int> vertex_tet_;  // an alive tetrahedron incident to each vertex
    mutable std::vector<std::uint32_t> mark_;  // per tetrahedron slot: the stamp of the last search
    mutable std::uint32_t stamp_ = 0;
    mutable std::vector<int> walk_;  // walk_around's stack, kept between calls
    std::vector<int> created_;
    int last_ = 0;
};

// The insertion order of a point set (DECISIONS.md M13): biased randomized insertion (BRIO) rounds,
// a point entering round r with probability 1/2 per round from the last (a splitmix64 hash of its
// index with a fixed seed, so the order is the same on every machine), each round in Hilbert-curve
// order of the points' bounding box (2^16 cells a side).
[[nodiscard]] std::vector<std::size_t> brio_order(const std::vector<Point>& p);

}  // namespace NiTCAD::mesh::detail
