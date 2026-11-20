// Layered 2D geometry resolved into a tiling (Unit 26k, GD1-GD2; generate.hpp).
//
// The overlay: every shape's x and y breakpoints (cut to the domain) make a grid of cells; each cell takes the
// region of the last shape covering it (an even-odd scan of the shape's vertical edges per row, exact: only
// comparisons of input coordinates), else the background. Same-region cells sharing an edge form parts; each part's
// boundary is traced with the part on its left (counter-clockwise). A part whose boundary touches itself at a grid
// vertex (two outgoing edges there) or has more than one loop (a hole) is cut along the vertical breakpoint line
// through that vertex or the hole's leftmost line, and its pieces are traced again: each cut strictly divides the
// part, so this ends. Last, every polygon takes on its edges the corners of all polygons (T-junctions).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "polygon2d.hpp"

namespace NiTCAD::mesh {

namespace {

std::unexpected<base::Error> invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return std::unexpected(base::Error{base::ErrorCode::invalid_input, std::move(message), base::ErrorContext{.index = index, .value = std::nullopt}});
}

// The breakpoints along one axis: the domain's ends and every vertex coordinate inside them.
std::vector<double> breakpoints(const LayeredGeometry2D& g, int axis) {
    std::vector<double> v{g.lo[static_cast<std::size_t>(axis)], g.hi[static_cast<std::size_t>(axis)]};
    for (const LayeredShape& s : g.shapes) {
        for (const Point& p : s.vertices) {
            const double c = p[static_cast<std::size_t>(axis)];
            if (c > v[0] && c < v[1]) v.push_back(c);
        }
    }
    std::ranges::sort(v);
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

// The index of the breakpoint at or nearest inside the domain for a coordinate (cut to the domain).
std::size_t index_of(const std::vector<double>& b, double c) {
    if (c <= b.front()) return 0;
    if (c >= b.back()) return b.size() - 1;
    return static_cast<std::size_t>(std::ranges::lower_bound(b, c) - b.begin());  // c is a breakpoint
}

}  // namespace

std::expected<Geometry2D, base::Error> layered_partition2d(const LayeredGeometry2D& g) {
    for (int a = 0; a < 2; ++a) {
        const auto k = static_cast<std::size_t>(a);
        if (!(std::isfinite(g.lo[k]) && std::isfinite(g.hi[k]) && g.lo[k] < g.hi[k])) {
            return invalid("the layered domain must be finite with lo < hi");
        }
    }
    for (std::size_t s = 0; s < g.shapes.size(); ++s) {
        const auto& v = g.shapes[s].vertices;
        if (g.shapes[s].region < 0) return invalid("a layered shape's region must be >= 0", s);
        if (v.size() < 4) return invalid("a layered shape needs at least four vertices", s);
        for (std::size_t k = 0; k < v.size(); ++k) {
            const Point& p = v[k];
            const Point& q = v[(k + 1) % v.size()];
            if (!(std::isfinite(p[0]) && std::isfinite(p[1]))) return invalid("a layered shape's vertices must be finite", s);
            const bool horizontal = p[1] == q[1], vertical = p[0] == q[0];
            if (horizontal == vertical) {  // slanted, or a repeated vertex
                return invalid(horizontal ? "a layered shape repeats a vertex"
                                          : "a layered shape's edges must be horizontal or vertical (slanted edges: a later slice)",
                               s);
            }
        }
        if (!detail::simple_2d(v) || detail::polygon_turn(v) == 0) return invalid("a layered shape must be simple, of non-zero area", s);
    }

    std::vector<double> xs = breakpoints(g, 0), ys = breakpoints(g, 1);
    std::size_t nx = xs.size() - 1, ny = ys.size() - 1;
    // Each row painted from the top shape down, every cell once: `free_from` skips the cells already taken
    // (path-halving union-find over the row's columns).
    std::vector<std::int32_t> owner(nx * ny, g.background);
    std::vector<std::array<double, 2>> yrange(g.shapes.size());
    for (std::size_t s = 0; s < g.shapes.size(); ++s) {
        double lo = g.shapes[s].vertices[0][1], hi = lo;
        for (const Point& p : g.shapes[s].vertices) {
            lo = std::min(lo, p[1]);
            hi = std::max(hi, p[1]);
        }
        yrange[s] = {lo, hi};
    }
    std::vector<std::size_t> cross, free_from(nx + 1);
    const auto find = [&](std::size_t i) {
        while (free_from[i] != i) {
            free_from[i] = free_from[free_from[i]];
            i = free_from[i];
        }
        return i;
    };
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) free_from[i] = i;
        for (std::size_t s = g.shapes.size(); s-- > 0;) {
            if (yrange[s][0] > ys[j] || yrange[s][1] < ys[j + 1]) continue;
            const auto& v = g.shapes[s].vertices;
            cross.clear();
            for (std::size_t k = 0; k < v.size(); ++k) {
                const Point& p = v[k];
                const Point& q = v[(k + 1) % v.size()];
                if (p[0] != q[0]) continue;  // horizontal
                const double lo = std::min(p[1], q[1]), hi = std::max(p[1], q[1]);
                if (lo <= ys[j] && hi >= ys[j + 1]) cross.push_back(index_of(xs, p[0]));
            }
            std::ranges::sort(cross);
            for (std::size_t c = 0; c + 1 < cross.size(); c += 2) {
                for (std::size_t i = find(cross[c]); i < cross[c + 1]; i = find(i + 1)) {
                    owner[j * nx + i] = g.shapes[s].region;
                    free_from[i] = i + 1;
                }
            }
        }
    }
    if (const auto it = std::ranges::find(owner, -1); it != owner.end()) {
        const auto c = static_cast<std::size_t>(it - owner.begin());
        return invalid("the layers leave part of the domain uncovered (near (" + std::to_string(xs[c % nx]) + ", " +
                       std::to_string(ys[c / nx]) + ")) and no region is the background");
    }
    // Breakpoint lines with the same owners on both sides along their whole length (hidden under later shapes)
    // carry no edge of any polygon: merged away, exactly, so the later passes see a smaller grid.
    {
        std::vector<char> keep_x(nx + 1, 0), keep_y(ny + 1, 0);
        keep_x[0] = keep_x[nx] = keep_y[0] = keep_y[ny] = 1;
        for (std::size_t j = 0; j < ny; ++j) {
            const std::int32_t* r = owner.data() + j * nx;
            for (std::size_t i = 1; i < nx; ++i) keep_x[i] |= static_cast<char>(r[i - 1] != r[i]);
            if (j > 0) {
                const std::int32_t* below = r - nx;
                for (std::size_t i = 0; i < nx && keep_y[j] == 0; ++i) keep_y[j] = static_cast<char>(below[i] != r[i]);
            }
        }
        std::vector<std::size_t> cx, cy;  // the first old column (row) of each new one
        std::vector<double> x2, y2;
        for (std::size_t i = 0; i <= nx; ++i) {
            if (keep_x[i] == 0) continue;
            x2.push_back(xs[i]);
            if (i < nx) cx.push_back(i);
        }
        for (std::size_t j = 0; j <= ny; ++j) {
            if (keep_y[j] == 0) continue;
            y2.push_back(ys[j]);
            if (j < ny) cy.push_back(j);
        }
        if (cx.size() < nx || cy.size() < ny) {
            std::vector<std::int32_t> o2(cx.size() * cy.size());
            for (std::size_t j = 0; j < cy.size(); ++j) {
                for (std::size_t i = 0; i < cx.size(); ++i) o2[j * cx.size() + i] = owner[cy[j] * nx + cx[i]];
            }
            owner = std::move(o2);
            xs = std::move(x2);
            ys = std::move(y2);
            nx = cx.size();
            ny = cy.size();
        }
    }
    // Each cell's column and row, looked up (a division per cell per pass costs more than the pass).
    std::vector<std::uint32_t> col(nx * ny), row(nx * ny);
    for (std::size_t j = 0, c = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i, ++c) {
            col[c] = static_cast<std::uint32_t>(i);
            row[c] = static_cast<std::uint32_t>(j);
        }
    }

    // Parts: same-region cells sharing an edge, labelled by a flood from each unlabelled cell in row-major order.
    // A part's boundary is its directed unit edges with the part on the left, found for every part in one pass.
    constexpr std::uint32_t none = UINT32_MAX;
    std::vector<std::uint32_t> label(nx * ny, none);
    std::vector<std::size_t> stack;
    struct PartInfo {
        std::int32_t region;
        std::size_t seed;                                       // a cell of the part
        std::vector<std::pair<std::size_t, std::size_t>> edges;  // grid vertex keys
    };
    std::vector<PartInfo> parts;
    const std::size_t vx = nx + 1;
    const auto key = [&](std::size_t i, std::size_t j) { return j * vx + i; };
    // Flood from c0 over its region's cells with label `from` whose column lies in [lo, hi): they take `to`. A span
    // fill: each popped seed's row run is taken whole, and one seed is pushed per run of matches above and below.
    const auto flood = [&](std::size_t c0, std::uint32_t from, std::uint32_t to, std::size_t lo, std::size_t hi,
                           std::vector<std::size_t>* cells) {
        const std::int32_t region = owner[c0];
        const auto match = [&](std::size_t n) { return label[n] == from && owner[n] == region; };
        stack.assign(1, c0);
        while (!stack.empty()) {
            std::size_t c = stack.back();
            stack.pop_back();
            if (!match(c)) continue;
            const std::size_t j = row[c], base = j * nx;
            std::size_t a = col[c];
            while (a > lo && match(base + a - 1)) --a;
            std::size_t b = a;
            while (b < hi && match(base + b)) {
                label[base + b] = to;
                if (cells) cells->push_back(base + b);
                ++b;
            }
            for (const std::size_t r : {j - 1, j + 1}) {  // j - 1 wraps past 0 when j == 0: skipped as >= ny
                if (r >= ny) continue;
                bool run = false;
                for (std::size_t x = a; x < b; ++x) {
                    const bool m = match(r * nx + x);
                    if (m && !run) stack.push_back(r * nx + x);
                    run = m;
                }
            }
        }
    };
    for (std::size_t c = 0; c < nx * ny; ++c) {
        if (label[c] != none) continue;
        parts.push_back({owner[c], c, {}});
        flood(c, none, static_cast<std::uint32_t>(parts.size() - 1), 0, nx, nullptr);
    }
    // The boundary edges of the cells given (all of them, or a cut part's).
    const auto boundary = [&](std::size_t c) {
        const std::uint32_t l = label[c];
        auto& e = parts[l].edges;
        const std::size_t i = col[c], j = row[c];
        if (j == 0 || label[c - nx] != l) e.emplace_back(key(i, j), key(i + 1, j));
        if (i + 1 == nx || label[c + 1] != l) e.emplace_back(key(i + 1, j), key(i + 1, j + 1));
        if (j + 1 == ny || label[c + nx] != l) e.emplace_back(key(i + 1, j + 1), key(i, j + 1));
        if (i == 0 || label[c - 1] != l) e.emplace_back(key(i, j + 1), key(i, j));
    };
    for (std::size_t c = 0; c < nx * ny; ++c) boundary(c);

    // A hole's cut (a ring has no simple outline, Unit 26k slice 1): two rays along grid lines from the hole's
    // corners out to the part's outer loop. A straight cut across the whole part would also meet the domain's
    // boundary, and a vertex inside a patch there breaks generate2d's rule that a patch lies inside one boundary
    // edge, and one ending on an interface makes a corner zone generate2d's P4b check can refuse; so a ray ends, best
    // first: on the domain's boundary off every patch's inside (a corner, then not), on an interface's corner, on an
    // interface, inside a patch; then the shorter. The
    // part's cells, walled along the rays, fall into pieces, each a new part (traced again: other holes, pinches).
    // 1: cut; 0: fewer than two rays (nothing changed); -1: the rays did not divide the part.
    std::vector<char> on_outer((nx + 1) * (ny + 1), 0), vwall(nx * ny, 0), hwall(nx * ny, 0);
    std::vector<std::size_t> cells;
    const auto in_patch = [&](std::size_t v) {
        const std::size_t i = v % vx, j = v / vx;
        if (i != 0 && i != nx && j != 0 && j != ny) return false;
        const double x = xs[i], y = ys[j];
        for (const PatchSegment& p : g.patches) {
            if (p.a[0] == x && p.b[0] == x && y > std::min(p.a[1], p.b[1]) && y < std::max(p.a[1], p.b[1])) return true;
            if (p.a[1] == y && p.b[1] == y && x > std::min(p.a[0], p.b[0]) && x < std::max(p.a[0], p.b[0])) return true;
        }
        return false;
    };
    const auto corner_at = [](const std::vector<std::size_t>& l, std::size_t k) {
        const std::size_t n = l.size();
        const auto d1 = static_cast<std::ptrdiff_t>(l[k]) - static_cast<std::ptrdiff_t>(l[(k + n - 1) % n]);
        const auto d2 = static_cast<std::ptrdiff_t>(l[(k + 1) % n]) - static_cast<std::ptrdiff_t>(l[k]);
        return d1 != d2;
    };
    const auto ray_cut = [&](std::size_t q, const std::vector<std::vector<std::size_t>>& ls) -> int {
        const auto lq = static_cast<std::uint32_t>(q);
        const auto mine_cell = [&](std::ptrdiff_t i, std::ptrdiff_t j) {
            return i >= 0 && j >= 0 && i < static_cast<std::ptrdiff_t>(nx) && j < static_cast<std::ptrdiff_t>(ny) &&
                   label[static_cast<std::size_t>(j) * nx + static_cast<std::size_t>(i)] == lq;
        };
        // One unit step from (i, j) in direction d (0 up, 1 down, 2 right, 3 left): whether both cells along it are
        // the part's; `wall` the wall it sets (`vertical`: a wall on the cell's left side, else on its bottom).
        const auto step = [&](std::ptrdiff_t& i, std::ptrdiff_t& j, int d, bool set) {
            std::ptrdiff_t ai = i - 1, aj = j, bi = i, bj = j, ni = i, nj = j + 1;
            if (d == 1) { aj = j - 1; bj = j - 1; nj = j - 1; }
            if (d == 2) { ai = i; aj = j - 1; nj = j; ni = i + 1; }
            if (d == 3) { ai = i - 1; aj = j - 1; bi = i - 1; bj = j; ni = i - 1; nj = j; }
            if (!mine_cell(ai, aj) || !mine_cell(bi, bj)) return false;
            if (set) {
                const std::size_t wc = static_cast<std::size_t>(bj) * nx + static_cast<std::size_t>(bi);
                (d < 2 ? vwall : hwall)[wc] = 1;
            }
            i = ni;
            j = nj;
            return true;
        };
        const auto& outer = ls[0];
        const auto& hole = ls[1];
        for (std::size_t k = 0; k < outer.size(); ++k) on_outer[outer[k]] = corner_at(outer, k) ? 2 : 1;
        struct Ray {
            std::size_t start;
            int dir, score;
            std::size_t steps;
        };
        std::vector<Ray> found;
        for (std::size_t k = 0; k < hole.size(); ++k) {
            if (!corner_at(hole, k)) continue;
            for (int d = 0; d < 4; ++d) {
                auto i = static_cast<std::ptrdiff_t>(hole[k] % vx), j = static_cast<std::ptrdiff_t>(hole[k] / vx);
                std::size_t n = 0;
                while (step(i, j, d, false)) ++n;
                const std::size_t end = static_cast<std::size_t>(j) * vx + static_cast<std::size_t>(i);
                if (n == 0 || on_outer[end] == 0) continue;
                const std::size_t ei = end % vx, ej = end / vx;
                const bool wall = ei == 0 || ei == nx || ej == 0 || ej == ny, corner = on_outer[end] == 2;
                const int score = in_patch(end) ? 16 : wall ? (corner ? 0 : 1) : (corner ? 4 : 8);
                found.push_back({hole[k], d, score, n});
            }
        }
        for (const std::size_t v : outer) on_outer[v] = 0;
        std::ranges::sort(found, [](const Ray& a, const Ray& b) { return std::tie(a.score, a.steps, a.start, a.dir) < std::tie(b.score, b.steps, b.start, b.dir); });
        // The best ray, and the best from another corner (two from one corner leave the rest touching itself there).
        const auto second = std::ranges::find_if(found, [&](const Ray& r) { return found.empty() || r.start != found.front().start; });
        if (found.empty() || second == found.end()) return 0;
        for (const Ray* r : {&found.front(), &*second}) {
            auto i = static_cast<std::ptrdiff_t>(r->start % vx), j = static_cast<std::ptrdiff_t>(r->start / vx);
            while (step(i, j, r->dir, true)) {}
        }
        // The pieces: a plain flood that stops at walls.
        const std::uint32_t work = none - 1;
        cells.clear();
        flood(parts[q].seed, lq, work, 0, nx, &cells);
        const auto [lo_c, hi_c] = std::ranges::minmax(cells);
        const std::size_t first = parts.size();
        for (std::size_t c0 = lo_c; c0 <= hi_c; ++c0) {
            if (label[c0] != work) continue;
            const auto to = static_cast<std::uint32_t>(parts.size());
            parts.push_back({parts[q].region, c0, {}});
            stack.assign(1, c0);
            label[c0] = to;
            while (!stack.empty()) {
                const std::size_t c = stack.back();
                stack.pop_back();
                const std::size_t i = col[c], j = row[c];
                const auto visit = [&](std::size_t n) {
                    if (label[n] == work) {
                        label[n] = to;
                        stack.push_back(n);
                    }
                };
                if (i > 0 && vwall[c] == 0) visit(c - 1);
                if (i + 1 < nx && vwall[c + 1] == 0) visit(c + 1);
                if (j > 0 && hwall[c] == 0) visit(c - nx);
                if (j + 1 < ny && hwall[c + nx] == 0) visit(c + nx);
            }
        }
        for (const std::size_t c : cells) vwall[c] = hwall[c] = 0;
        if (parts.size() - first < 2) return -1;
        for (const std::size_t c : cells) boundary(c);
        return 1;
    };

    // Grid vertex (i, j) as one key; a loop as its unit steps' vertices.
    std::vector<std::pair<std::int32_t, std::vector<std::size_t>>> loops;  // region, vertices (every unit step)
    constexpr std::size_t npos = SIZE_MAX;
    std::vector<std::size_t> next((nx + 1) * (ny + 1), npos);  // a boundary vertex's outgoing edge's end
    std::vector<char> seen((nx + 1) * (ny + 1), 0);
    std::vector<std::size_t> touched, cuts;
    for (std::size_t q = 0; q < parts.size(); ++q) {
        // A vertex with two outgoing edges is a pinch: a cut along its column line.
        touched.clear();
        cuts.clear();
        for (const auto& [a, b] : parts[q].edges) {
            if (next[a] != npos) {
                cuts.push_back(a % vx);
                continue;
            }
            next[a] = b;
            touched.push_back(a);
        }
        std::vector<std::vector<std::size_t>> mine;
        int rays = 0;
        if (cuts.empty()) {
            // Loops from their lowest-then-leftmost vertex (keys order by row, then column); the first is the outer
            // one, each other a hole: the first hole cut by two rays (ray_cut), or, without two, every hole by a cut
            // along its leftmost line.
            std::ranges::sort(touched);
            for (const std::size_t s : touched) {
                if (seen[s] != 0) continue;
                std::vector<std::size_t> loop;
                for (std::size_t v = s; seen[v] == 0; v = next[v]) {
                    seen[v] = 1;
                    loop.push_back(v);
                }
                mine.push_back(std::move(loop));
            }
            if (mine.size() > 1) {
                rays = ray_cut(q, mine);
                if (rays < 0) return invalid("internal: a layered part's cut did not divide it");
                for (std::size_t k = 1; rays == 0 && k < mine.size(); ++k) {
                    std::size_t left = npos;
                    for (const std::size_t v : mine[k]) left = std::min(left, v % vx);
                    cuts.push_back(left);
                }
            }
        }
        for (const std::size_t v : touched) {
            next[v] = npos;
            seen[v] = 0;
        }
        parts[q].edges.clear();
        parts[q].edges.shrink_to_fit();
        if (rays > 0) continue;
        if (cuts.empty()) {
            loops.emplace_back(parts[q].region, std::move(mine.front()));
            continue;
        }
        // Every cut at once: the part's cells split into the bands between cut lines, each band's pieces new parts
        // (cutting adds no hole or pinch), their edges found again.
        std::ranges::sort(cuts);
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        const std::uint32_t old = static_cast<std::uint32_t>(q), work = none - 1;
        cells.clear();
        flood(parts[q].seed, old, work, 0, nx, &cells);
        // In row-major order without a sort: a scan of the part's index range.
        const auto [lo_c, hi_c] = std::ranges::minmax(cells);
        const std::size_t first = parts.size();
        for (std::size_t c = lo_c; c <= hi_c; ++c) {
            if (label[c] != work) continue;
            const auto band = static_cast<std::size_t>(std::ranges::upper_bound(cuts, col[c]) - cuts.begin());
            const std::size_t lo = band == 0 ? 0 : cuts[band - 1], hi = band == cuts.size() ? nx : cuts[band];
            parts.push_back({parts[q].region, c, {}});
            flood(c, work, static_cast<std::uint32_t>(parts.size() - 1), lo, hi, nullptr);
        }
        if (parts.size() - first < 2) return invalid("internal: a layered part's cut did not divide it");
        for (const std::size_t c : cells) boundary(c);
    }
    // The corners of every loop; each polygon keeps the vertices along its boundary that are any loop's corner.
    std::vector<char> corners((nx + 1) * (ny + 1), 0);
    const auto is_corner = [&](const std::vector<std::size_t>& l, std::size_t k) {
        const std::size_t a = l[(k + l.size() - 1) % l.size()], b = l[k], c = l[(k + 1) % l.size()];
        const bool ab_h = a / vx == b / vx, bc_h = b / vx == c / vx;
        return ab_h != bc_h;
    };
    for (const auto& [region, l] : loops) {
        for (std::size_t k = 0; k < l.size(); ++k) {
            if (is_corner(l, k)) corners[l[k]] = 1;
        }
    }
    Geometry2D out;
    out.patches = g.patches;
    for (const auto& [region, l] : loops) {
        RegionPolygon poly;
        poly.region = region;
        for (const std::size_t v : l) {
            if (corners[v] != 0) poly.vertices.push_back({xs[v % vx], ys[v / vx], 0.0});
        }
        out.polygons.push_back(std::move(poly));
    }
    return out;
}

}  // namespace NiTCAD::mesh
