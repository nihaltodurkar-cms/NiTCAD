// An incremental 2D Delaunay triangulation (Bowyer-Watson) for the mesh generator (Unit 16;
// DECISIONS.md M3, M4, M13). Internal to the mesh layer.
//
// Every topological decision is an exact predicate: orientation for point location, and the
// perturbed incircle test (the node-id tie-break) for the cavity, so ties resolve one consistent
// way and the result is the Delaunay triangulation of the perturbed lifting. The points lie inside
// a large super triangle whose three vertices are the triangulation's vertices 0, 1, 2; the
// caller's points are vertices 3, 4, ... in insertion order (deterministic: no hashing, no
// randomness). For the tie-break a super vertex ranks above every real one.
//
// Point location walks from a given triangle across edges the point is beyond (a visibility walk;
// it terminates on regular, hence Delaunay, triangulations; a step cap falls back to a scan).
#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh::detail {

class Triangulation2D {
public:
    struct Triangle {
        std::array<int, 3> v{};  // counterclockwise
        std::array<int, 3> n{};  // n[i]: the neighbour across the edge opposite v[i], or -1
        bool alive = false;
    };

    // The cavity of a point: the triangles whose perturbed circumcircles contain it, and the
    // cavity's boundary edges (a, b) as in their triangle (counterclockwise), each with the
    // triangle outside it (-1 on the super triangle's hull).
    struct Cavity {
        std::vector<int> triangles;
        struct Edge {
            int a;
            int b;
            int outside;
        };
        std::vector<Edge> boundary;
    };

    // A super triangle enclosing the box [lo, hi] (2D) with a wide margin.
    Triangulation2D(const Point& lo, const Point& hi);

    [[nodiscard]] std::size_t vertex_count() const noexcept { return points_.size(); }
    [[nodiscard]] const Point& point(int v) const { return points_[static_cast<std::size_t>(v)]; }
    [[nodiscard]] bool is_super(int v) const noexcept { return v < 3; }
    [[nodiscard]] std::size_t triangle_slots() const noexcept { return tris_.size(); }
    [[nodiscard]] const Triangle& triangle(int t) const { return tris_[static_cast<std::size_t>(t)]; }
    // A triangle created most recently: a good start for a walk.
    [[nodiscard]] int last_triangle() const noexcept { return last_; }

    // The alive triangle containing p (on its boundary counts), walking from `start`.
    [[nodiscard]] int locate(const Point& p, int start) const;
    // The cavity of p, found from the triangle containing it. Does not modify the triangulation.
    [[nodiscard]] Cavity cavity(const Point& p, int containing) const;
    // Inserts p with this cavity (from cavity()); returns its vertex id.
    int insert(const Point& p, const Cavity& c);
    // locate, cavity and insert. Precondition: p strictly inside the super triangle, not a vertex.
    int insert(const Point& p);

    // The alive triangle with directed edge a -> b (counterclockwise in it), if any.
    [[nodiscard]] std::optional<int> triangle_with_edge(int a, int b) const;
    // The id the tie-break uses for vertex v.
    [[nodiscard]] NodeId tie_id(int v) const noexcept;

private:
    [[nodiscard]] bool in_circle(int t, const Point& p, int pid) const;
    int new_triangle(const std::array<int, 3>& v);

    std::vector<Point> points_;
    std::vector<Triangle> tris_;
    std::vector<int> free_;
    std::vector<int> vertex_triangle_;  // an alive triangle incident to each vertex
    int last_ = 0;
};

}  // namespace NiTCAD::mesh::detail
