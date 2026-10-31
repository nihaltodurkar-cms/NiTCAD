// The validated polyhedral complex (Unit 53; DECISIONS.md FD1, FD6): the facets' kinds, planes and
// triangulations, and exact point location. Internal to the mesh layer.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/complex.hpp"

namespace NiTCAD::mesh::detail {

enum class FacetKind : std::uint8_t { boundary, interface, internal };

struct IndexedFacet {
    FacetKind kind = FacetKind::boundary;
    std::array<int, 2> volumes{-1, -1};  // the volumes beside it (the second -1 on the boundary)
    std::array<std::size_t, 3> plane{};  // three non-collinear vertices, in facet order
    int drop = 0;                        // the coordinate dropped to project the facet one to one
    std::vector<std::array<std::size_t, 3>> triangles;  // vertex indices, in the facet's orientation
    int patch = -1;                      // a boundary facet's index in ComplexReport::patches
    Point lo{}, hi{};                    // its bounding box (a filter before the exact tests)
};

class ComplexIndex {
public:
    // Validates (check_complex's rules) and indexes. Errors: invalid_input.
    [[nodiscard]] static std::expected<ComplexIndex, base::Error> create(PolyhedralGeometry geometry);

    [[nodiscard]] const PolyhedralGeometry& geometry() const noexcept { return g_; }
    [[nodiscard]] const std::vector<IndexedFacet>& facets() const noexcept { return facets_; }
    [[nodiscard]] const ComplexReport& report() const noexcept { return report_; }

    struct Location {
        enum class Kind : std::uint8_t { outside, inside, on_facet, undetermined } kind = Kind::outside;
        int index = -1;  // the volume (inside) or the lowest-index facet containing the point (on_facet)
    };
    // Where p lies, exactly (FD6): on a facet (closed), strictly inside a volume, or outside the
    // domain. Inside is decided by the parity of a segment's crossings of the facet triangles,
    // from p to a point outside the bounding box; a segment meeting a triangle's edge or vertex,
    // or lying in its plane, is degenerate and replaced by the next of a fixed sequence of
    // segments, so the answer is exact and deterministic. undetermined only if 64 segments in a
    // row are degenerate.
    [[nodiscard]] Location locate(const Point& p) const;

private:
    // locate among the listed facets only, as the surface of one volume (index 0), or among all
    // facets (only null). The volume-overlap check asks this of each volume's own facets.
    [[nodiscard]] Location locate(const Point& p, const std::vector<std::size_t>* only) const;

    PolyhedralGeometry g_;
    std::vector<IndexedFacet> facets_;
    ComplexReport report_;
    Point lo_{}, hi_{};  // the bounding box
};

}  // namespace NiTCAD::mesh::detail
