// Feature edges and the enhanced 3D look (Unit 26f; DECISIONS.md KD1, KD3).
#include "NiTCAD/render/edges.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numbers>
#include <string_view>
#include <unordered_map>

namespace NiTCAD::render {

namespace {

// -0 and +0 are one position: both compare equal, so they must hash alike.
Float3 canonical(const Float3& p) noexcept { return {p[0] + 0.0f, p[1] + 0.0f, p[2] + 0.0f}; }

bool before(const Float3& p, const Float3& q) noexcept {
    return std::lexicographical_compare(p.begin(), p.end(), q.begin(), q.end());
}

struct EdgeHash {
    std::size_t operator()(const Segment& s) const noexcept {
        std::uint64_t h = 1469598103934665603ull;
        for (const Float3& p : s) {
            for (const float f : p) {
                std::uint32_t u = 0;
                std::memcpy(&u, &f, sizeof u);
                h = (h ^ u) * 1099511628211ull;
            }
        }
        return static_cast<std::size_t>(h);
    }
};

struct EdgeFaces {
    std::array<double, 3> normal{};  // of the first face
    int faces = 0;
    bool sharp = false;
};

constexpr std::string_view edges_name = "feature edges";

}  // namespace

std::vector<Segment> feature_edges(const TriangleBatch& batch, float angle_deg) {
    const double cos_limit = std::cos(static_cast<double>(angle_deg) * std::numbers::pi / 180.0);
    std::unordered_map<Segment, EdgeFaces, EdgeHash> edges;
    edges.reserve(batch.vertices.size());
    for (std::size_t t = 0; t + 2 < batch.vertices.size(); t += 3) {
        const Float3 p[3] = {canonical(batch.vertices[t].position), canonical(batch.vertices[t + 1].position),
                             canonical(batch.vertices[t + 2].position)};
        const std::array<double, 3> u{double(p[1][0]) - p[0][0], double(p[1][1]) - p[0][1], double(p[1][2]) - p[0][2]};
        const std::array<double, 3> v{double(p[2][0]) - p[0][0], double(p[2][1]) - p[0][1], double(p[2][2]) - p[0][2]};
        std::array<double, 3> n{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
        const double length = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (!(length > 0.0) || !std::isfinite(length)) continue;  // degenerate
        for (double& c : n) c /= length;
        for (int e = 0; e < 3; ++e) {
            const Float3& a = p[e];
            const Float3& b = p[(e + 1) % 3];
            const Segment key = before(b, a) ? Segment{b, a} : Segment{a, b};
            auto [it, fresh] = edges.try_emplace(key, EdgeFaces{n, 1, false});
            if (fresh) continue;
            EdgeFaces& f = it->second;
            ++f.faces;
            // Either winding of the neighbour: the angle between the faces' planes.
            const double d = std::abs(f.normal[0] * n[0] + f.normal[1] * n[1] + f.normal[2] * n[2]);
            if (d < cos_limit) f.sharp = true;
        }
    }
    std::vector<Segment> out;
    for (const auto& [segment, f] : edges) {
        if (f.faces != 2 || f.sharp) out.push_back(segment);
    }
    std::sort(out.begin(), out.end(), [](const Segment& x, const Segment& y) {
        if (x[0] != y[0]) return before(x[0], y[0]);
        return before(x[1], y[1]);
    });
    return out;
}

bool takes_feature_edges(const TriangleBatch& batch) noexcept {
    return batch.lit && batch.name != "arrows" && batch.name != "isosurface" && batch.name != "outline";
}

void set_look(FieldScene& scene, Look look) {
    scene.look = look;
    std::erase_if(scene.lines, [](const LineBatch& l) { return l.name == edges_name; });
    if (look != Look::enhanced || scene.dimension != 3) return;
    LineBatch lines{std::string(edges_name), {0.08f, 0.09f, 0.11f, 0.85f}, 1.7f, false, {}};
    for (const TriangleBatch& b : scene.triangles) {
        if (!takes_feature_edges(b)) continue;
        for (const Segment& s : feature_edges(b)) {
            lines.points.push_back(s[0]);
            lines.points.push_back(s[1]);
        }
    }
    if (!lines.points.empty()) scene.lines.push_back(std::move(lines));
}

}  // namespace NiTCAD::render
