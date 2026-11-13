// The adaptive refinement loop (Unit 17, slices 3-5; DECISIONS.md Unit 17 AD1, AD4-AD6, AD9, AD11-AD13,
// ID1-ID7, AR5, AR7).
#include "NiTCAD/solve/adapt.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/mesh/locate.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/solve/adapt_adjoint.hpp"

namespace NiTCAD::solve {
namespace {

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

using clock_type = std::chrono::steady_clock;
double seconds_since(clock_type::time_point t0) {
    return std::chrono::duration<double>(clock_type::now() - t0).count();
}

// Per edge, the largest excess over the kinds (AD2's normalization).
std::vector<double> edge_excess(const EdgeIndicators& e, const IndicatorOptions& o) {
    const std::array<double, indicator_count> t{o.potential_vt, o.density_decades, o.density_decades, o.doping_decades,
                                                o.generation_share};
    std::vector<double> x(e.raw[0].size(), 0.0);
    for (std::size_t k = 0; k < indicator_count; ++k) {
        for (std::size_t j = 0; j < x.size(); ++j) x[j] = std::max(x[j], std::max(0.0, e.raw[k][j] / t[k] - 1.0));
    }
    return x;
}

// The 2D refinement of a triangulation's nodal sizes (AD4): marked nodes halved, or (targets given,
// ID2) each node's size lowered to its target; `changed` counts the nodes lowered.
std::expected<std::shared_ptr<const mesh::SizeField>, base::Error> refine_field(
    const mesh::SimplexMesh& sm, const std::shared_ptr<const mesh::SizeField>& previous, const std::vector<char>& marked,
    double min_size, double grading, bool uniform, const std::vector<double>* targets = nullptr,
    std::size_t* changed = nullptr) {
    std::vector<double> h = mesh::local_sizes(sm);
    if (previous) {
        for (std::size_t i = 0; i < h.size(); ++i) h[i] = std::min(h[i], previous->at(sm.points[i]));
    }
    std::size_t count = 0;
    if (uniform) {
        for (double& x : h) x *= 0.5;
        count = h.size();
    } else {
        std::vector<char> seeds(h.size(), 0);
        for (std::size_t i = 0; i < h.size(); ++i) {
            if (!(h[i] > min_size)) continue;
            if (targets != nullptr) {
                // Halved at most per call (as marks are), toward the target (AR5).
                const double t = std::max({(*targets)[i], 0.5 * h[i], min_size});
                if (!(t < h[i])) continue;
                h[i] = t;
            } else {
                if (marked[i] == 0) continue;
                h[i] = std::max(0.5 * h[i], min_size);
            }
            seeds[i] = 1;
            ++count;
        }
        auto g = mesh::grade_sizes(sm, h, grading, seeds);
        if (!g) return std::unexpected(g.error());
    }
    if (changed != nullptr) *changed = count;
    return mesh::SizeField::create(sm.points, sm.cells, std::move(h));
}

double smallest_interval(const std::vector<double>& a) {
    double m = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i + 1 < a.size(); ++i) m = std::min(m, a[i + 1] - a[i]);
    return m;
}

// An axis's smallest interval for refinement: the option, else a quarter of the narrowest straddle,
// else 1/64 of the axis's smallest interval.
double axis_floor(const AdaptOptions& o, const std::vector<double>& axis, const std::vector<mesh::AxisRange>& prot) {
    if (o.min_size > 0.0) return o.min_size;
    double w = std::numeric_limits<double>::infinity();
    for (const mesh::AxisRange& r : prot) w = std::min(w, r.hi - r.lo);
    if (std::isfinite(w)) return 0.25 * w;
    return smallest_interval(axis) / 64.0;
}

std::expected<std::vector<double>, base::Error> refine_lines(const std::vector<double>& lines,
                                                             const std::vector<char>& interval_marks,
                                                             const std::vector<mesh::AxisRange>& prot,
                                                             const AdaptOptions& o, bool uniform) {
    std::vector<char> m = uniform ? std::vector<char>(lines.size() - 1, 1) : interval_marks;
    auto r = mesh::refine_axis(lines, m, o.grading, uniform ? 0.0 : axis_floor(o, lines, prot), prot);
    if (!r) return std::unexpected(r.error());
    return std::move(r->lines);
}

// The refinement of `cur` from node marks (the 2D size field) and an edge predicate (axis
// intervals), or from per-node targets (ID2: the size field lowered to them; axis intervals by the
// predicate). `changed`: whether anything was refined.
std::expected<MeshRefinement, base::Error> refine_with(const AdaptableMesh& cur, const MeshRefinement& previous,
                                                       const std::vector<char>& mk,
                                                       const std::function<bool(std::size_t)>& edge_marked,
                                                       const std::vector<double>* targets, const AdaptOptions& o,
                                                       bool uniform, bool* changed);

}  // namespace

std::expected<MeshRefinement, base::Error> next_refinement(const AdaptableMesh& cur, const MeshRefinement& previous,
                                                           const std::vector<char>& marked, const EdgeIndicators* edges,
                                                           const IndicatorOptions& indicators, const AdaptOptions& o,
                                                           bool uniform) {
    const mesh::Mesh& m = cur.device.mesh();
    const std::size_t nn = m.node_count();
    if (!uniform && marked.size() != nn) return std::unexpected(invalid("marks must be one per node"));
    const std::vector<char> none(nn, 0);
    const std::vector<char>& mk = uniform ? none : marked;
    std::vector<double> excess;
    if (edges != nullptr) excess = edge_excess(*edges, indicators);
    const auto edge_marked = [&](std::size_t e) {
        if (uniform) return true;
        const auto a = static_cast<std::size_t>(m.edges()[e].first), b = static_cast<std::size_t>(m.edges()[e].second);
        const bool has = excess.empty() || excess[e] > 0.0;
        return has && (mk[a] != 0 || mk[b] != 0);
    };
    return refine_with(cur, previous, mk, edge_marked, nullptr, o, uniform, nullptr);
}

namespace {

std::expected<MeshRefinement, base::Error> refine_with(const AdaptableMesh& cur, const MeshRefinement& previous,
                                                       const std::vector<char>& mk,
                                                       const std::function<bool(std::size_t)>& edge_marked,
                                                       const std::vector<double>* targets, const AdaptOptions& o,
                                                       bool uniform, bool* changed) {
    const mesh::Mesh& m = cur.device.mesh();
    const std::size_t nn = m.node_count();
    bool any = false;
    const double floor2d = o.min_size > 0.0             ? o.min_size
                           : cur.interface_resolution > 0 ? 0.25 * cur.interface_resolution
                                                          : cur.max_size / 64.0;
    MeshRefinement next = previous;
    switch (cur.kind) {
        case AdaptMeshKind::generated2d: {
            if (!cur.generated || cur.generated->simplices.points.size() != nn) {
                return std::unexpected(invalid("a generated 2D mesh needs its simplices, one point per node"));
            }
            std::size_t count = 0;
            auto f = refine_field(cur.generated->simplices, previous.size_field, mk, floor2d, o.grading, uniform, targets,
                                  &count);
            if (!f) return std::unexpected(f.error());
            next.size_field = std::move(*f);
            if (changed != nullptr) *changed = count > 0;
            return next;
        }
        case AdaptMeshKind::tensor: {
            const mesh::TensorCells* tc = cur.device.cells();
            if (tc == nullptr) return std::unexpected(invalid("a tensor mesh needs its cells"));
            for (int a = 0; a < tc->dimension(); ++a) {
                const auto ax = tc->axis(a);
                std::vector<double> lines(ax.begin(), ax.end());
                std::vector<char> im(lines.size() - 1, 0);
                for (std::size_t e = 0; e < m.edges().size(); ++e) {
                    const auto& ed = m.edges()[e];
                    const auto ia = tc->index(ed.first), ib = tc->index(ed.second);
                    bool along = true;
                    for (int b = 0; b < tc->dimension(); ++b) {
                        if (b != a && ia[static_cast<std::size_t>(b)] != ib[static_cast<std::size_t>(b)]) along = false;
                    }
                    if (!along || ia[static_cast<std::size_t>(a)] == ib[static_cast<std::size_t>(a)]) continue;
                    if (edge_marked(e)) im[std::min(ia[static_cast<std::size_t>(a)], ib[static_cast<std::size_t>(a)])] = 1;
                }
                auto r = refine_lines(lines, im, cur.axis_protected[static_cast<std::size_t>(a)], o, uniform);
                if (!r) return std::unexpected(r.error());
                any = any || r->size() != lines.size();
                next.axes[static_cast<std::size_t>(a)] = std::move(*r);
            }
            if (changed != nullptr) *changed = any;
            return next;
        }
        case AdaptMeshKind::swept: {
            if (!cur.section) return std::unexpected(invalid("a swept mesh needs its cross-section"));
            const std::size_t n2 = cur.section->points.size();
            if (n2 == 0 || nn % n2 != 0) return std::unexpected(invalid("a swept mesh's nodes are not whole lines"));
            const std::size_t nl = nn / n2;
            const auto pts = m.points();
            const auto ax = static_cast<std::size_t>(cur.sweep_axis);
            std::vector<double> lines(nl);
            for (std::size_t j = 0; j < nl; ++j) lines[j] = pts[j * n2][ax];
            std::vector<char> section_marks(n2, 0);
            for (std::size_t i = 0; i < nn; ++i) {
                if (mk[i] != 0) section_marks[i % n2] = 1;
            }
            std::vector<char> im(nl - 1, 0);
            for (std::size_t e = 0; e < m.edges().size(); ++e) {
                const auto a = static_cast<std::size_t>(m.edges()[e].first), b = static_cast<std::size_t>(m.edges()[e].second);
                if (b != a + n2) continue;  // along the axis: the same section node on consecutive lines
                if (edge_marked(e)) im[a / n2] = 1;
            }
            // Targets on the cross-section: per section node, the smallest over the lines.
            std::vector<double> section_targets;
            if (targets != nullptr) {
                section_targets.assign(n2, std::numeric_limits<double>::infinity());
                for (std::size_t i = 0; i < nn; ++i) section_targets[i % n2] = std::min(section_targets[i % n2], (*targets)[i]);
            }
            std::size_t count = 0;
            auto f = refine_field(*cur.section, previous.size_field, section_marks, floor2d, o.grading, uniform,
                                  targets != nullptr ? &section_targets : nullptr, &count);
            if (!f) return std::unexpected(f.error());
            next.size_field = std::move(*f);
            auto r = refine_lines(lines, im, cur.sweep_protected, o, uniform);
            if (!r) return std::unexpected(r.error());
            if (changed != nullptr) *changed = count > 0 || r->size() != lines.size();
            // The generator adds the straddles' lines itself: drop them (the protected ranges' ends).
            std::vector<double> given;
            for (const double x : *r) {
                bool straddle_line = false;
                for (const mesh::AxisRange& p : cur.sweep_protected) straddle_line = straddle_line || x == p.lo || x == p.hi;
                if (!straddle_line) given.push_back(x);
            }
            next.sweep_lines = std::move(given);
            return next;
        }
    }
    return std::unexpected(invalid("unknown mesh kind"));
}

// ID6's DOF growth cap per pass.
constexpr double growth_cap = 2.5;

// The equation set's unknowns per node (ID4): Poisson alone 1; drift-diffusion 3 per semiconductor
// node, the electrothermal temperature one more and energy balance two more; insulator nodes the
// potential (and the temperature).
std::size_t unknowns(const BiasOptions& b, bool insulator) {
    if (b.equations == Equations::equilibrium_poisson) return 1;
    const std::size_t thermal = b.models.electrothermal ? 1 : 0;
    if (insulator) return 1 + thermal;
    return 3 + thermal + (b.models.energy_balance ? 2 : 0);
}

std::size_t dofs_of(const device::Device& d, const BiasOptions& b) {
    std::size_t n = 0;
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) n += unknowns(b, d.is_insulator(static_cast<mesh::NodeId>(i)));
    return n;
}

// Per kind and node, the largest over the node's edges of the raw value over its threshold (ID5's
// local error level, AR7: the excess vanishes below the threshold, the level does not).
std::array<std::vector<double>, indicator_count> node_levels(const mesh::Mesh& m, const EdgeIndicators& e,
                                                             const IndicatorOptions& o) {
    const std::array<double, indicator_count> t{o.potential_vt, o.density_decades, o.density_decades, o.doping_decades,
                                                o.generation_share};
    std::array<std::vector<double>, indicator_count> level;
    for (std::size_t k = 0; k < indicator_count; ++k) {
        level[k].assign(m.node_count(), 0.0);
        for (std::size_t j = 0; j < m.edges().size(); ++j) {
            const double r = e.raw[k][j] / t[k];
            for (const mesh::NodeId q : {m.edges()[j].first, m.edges()[j].second}) {
                auto& x = level[k][static_cast<std::size_t>(q)];
                x = std::max(x, r);
            }
        }
    }
    return level;
}

// Nodes sorted by each coordinate, for the nodes of a slab or a disc.
class AxisIndex {
public:
    explicit AxisIndex(std::span<const mesh::Point> pts) : pts_(pts) {
        for (std::size_t a = 0; a < 3; ++a) {
            order_[a].resize(pts.size());
            std::iota(order_[a].begin(), order_[a].end(), std::size_t{0});
            std::ranges::stable_sort(order_[a], [&](std::size_t i, std::size_t j) { return pts[i][a] < pts[j][a]; });
            coord_[a].resize(pts.size());
            for (std::size_t q = 0; q < pts.size(); ++q) coord_[a][q] = pts[order_[a][q]][a];
        }
    }
    template <typename F>
    void each_between(std::size_t a, double lo, double hi, F&& f) const {
        auto it = std::ranges::lower_bound(coord_[a], lo);
        for (auto q = static_cast<std::size_t>(it - coord_[a].begin()); q < coord_[a].size() && coord_[a][q] <= hi; ++q) {
            f(order_[a][q]);
        }
    }
    [[nodiscard]] std::span<const mesh::Point> points() const { return pts_; }

private:
    std::span<const mesh::Point> pts_;
    std::array<std::vector<std::size_t>, 3> order_;
    std::array<std::vector<double>, 3> coord_;
};

// Where a candidate refines (ID5's realized ratio): a slab [lo, hi] on an axis, or a disc of a
// radius about a centre over the nodes of one region (the skipped axis, a swept mesh's, ignored).
struct Place {
    int axis = -1;
    double lo = 0.0, hi = 0.0;
    mesh::Point centre{};
    double radius = 0.0;
    int skip_axis = -1;
    std::size_t region = 0;
};

// The weighted mean of a level over a place's nodes, weights omega (1 when empty); NaN without weight.
double local_level(const device::Device& d, const AxisIndex& index, const Place& w, const std::vector<double>& level,
                   std::span<const double> omega) {
    double sum = 0.0, weight = 0.0;
    const auto add = [&](std::size_t i) {
        const double om = omega.empty() ? 1.0 : omega[i];
        sum += om * level[i];
        weight += om;
    };
    if (w.axis >= 0) {
        index.each_between(static_cast<std::size_t>(w.axis), w.lo, w.hi, add);
    } else {
        const std::size_t a = w.skip_axis == 0 ? 1 : 0;
        const auto regions = d.node_region();
        const auto pts = index.points();
        index.each_between(a, w.centre[a] - w.radius, w.centre[a] + w.radius, [&](std::size_t i) {
            if (static_cast<std::size_t>(regions[i]) != w.region) return;
            double s = 0.0;
            for (std::size_t b = 0; b < 3; ++b) {
                if (static_cast<int>(b) == w.skip_axis) continue;
                s += (pts[i][b] - w.centre[b]) * (pts[i][b] - w.centre[b]);
            }
            if (std::sqrt(s) <= w.radius) add(i);
        });
    }
    return weight > 0.0 ? sum / weight : std::numeric_limits<double>::quiet_NaN();
}

// A candidate taken by ID6, waiting for the next pass's ratio.
struct Pending {
    Place where;
    Indicator kind = Indicator::potential;
    std::size_t region = 0;
    double before = 0.0;  // its local level on this pass
    double error = 0.0;   // its share of the estimate on this pass, relative to the allowance
    int levels = 1;       // halvings taken
    std::size_t listed = std::numeric_limits<std::size_t>::max();  // its index in the pass's list
};

struct Allocation {
    std::optional<MeshRefinement> next;
    std::vector<AdaptCandidate> listed;
    std::vector<Pending> pending;
    std::size_t count = 0;
    std::size_t taken = 0;
    double predicted_relative = std::numeric_limits<double>::quiet_NaN();
    std::size_t predicted_dofs = 0;
    bool budget_spent = false;  // budget mode: nothing fits
};

// What a pass gives the allocation.
struct PassState {
    double estimate = 0.0;  // the relative estimate E
    double newton = 0.0;
    double nnz = 0.0;
    double fill = 0.0;
    double dofs = 0.0;
    double nodes = 0.0;
    double seconds_so_far = 0.0;  // measured, this pass included
    double added_scale = 1.0;     // budget mode: the remaining DOF budget over this (grading's inflation)
};

constexpr std::size_t listed_candidates = 20;   // ID7
constexpr std::size_t measured_candidates = 200;  // AR7: the realized ratios measured per pass
constexpr int max_levels = 3;                     // AR7: an interval split at most into 8 per pass

bool inside_protected(double lo, double hi, const std::vector<mesh::AxisRange>& prot) {
    for (const mesh::AxisRange& r : prot) {
        if (std::min(hi, r.hi) - std::max(lo, r.lo) > 0.0) return true;
    }
    return false;
}

// An axis's lines with interval k of `lines` split to level[k] (1 halves, 2 quarters, 3 eighths):
// repeated midpoint splits, each graded, floored and protected as AD4's (refine_lines).
std::expected<std::vector<double>, base::Error> split_levels(const std::vector<double>& lines, const std::vector<int>& level,
                                                             const std::vector<mesh::AxisRange>& prot, const AdaptOptions& o) {
    std::vector<double> cur = lines;
    const int deepest = level.empty() ? 0 : *std::ranges::max_element(level);
    for (int l = 1; l <= deepest; ++l) {
        std::vector<char> m(cur.size() - 1, 0);
        for (std::size_t j = 0; j + 1 < cur.size(); ++j) {
            const double mid = 0.5 * (cur[j] + cur[j + 1]);
            const auto k = static_cast<std::size_t>(std::ranges::upper_bound(lines, mid) - lines.begin()) - 1;
            if (k < level.size() && level[k] >= l) m[j] = 1;
        }
        auto r = refine_lines(cur, m, prot, o, false);
        if (!r) return std::unexpected(r.error());
        cur = std::move(*r);
    }
    return cur;
}

// The deepest level an interval of length len can take above a floor.
int levels_within(double len, double floor) {
    int l = 0;
    while (l < max_levels && len * std::exp2(-(l + 1)) >= floor) ++l;
    return l;
}

// ID4 and ID6 (AR7): the candidates of `cur`, their benefits and costs, and the refinement taking
// them in decreasing benefit over work until the predicted estimate is at most min(1/2, E/2) (relative
// to the allowance: tol / 2 at most), the pass's DOF growth cap, or a budget.
std::expected<Allocation, base::Error> allocate(const AdaptableMesh& cur, const MeshRefinement& previous,
                                                const IndicatorField& united, const EdgeIndicators& edges,
                                                const std::array<std::vector<double>, indicator_count>& levels,
                                                std::span<const double> omega, const PassState& st,
                                                const CostModel& model, std::span<const OrderFit> fits,
                                                const BiasOptions& bias, const AdaptOptions& o) {
    Allocation out;
    const device::Device& dev = cur.device;
    const mesh::Mesh& m = dev.mesh();
    const std::size_t nn = m.node_count();
    const auto pts = m.points();
    const auto regions = dev.node_region();
    const auto w = [&](std::size_t i) { return omega.empty() ? 1.0 : omega[i]; };
    const auto u = [&](std::size_t i) { return unknowns(bias, dev.is_insulator(static_cast<mesh::NodeId>(i))); };
    const std::vector<double> ee = edge_excess(edges, o.indicators);
    const double degree = nn > 0 ? 2.0 * static_cast<double>(m.edges().size()) / static_cast<double>(nn) : 0.0;
    const double floor2d = o.min_size > 0.0             ? o.min_size
                           : cur.interface_resolution > 0 ? 0.25 * cur.interface_resolution
                                                          : cur.max_size / 64.0;
    struct Cand {
        AdaptCandidate c;
        double share = 0.0;
        std::size_t node = 0;  // carries the largest weighted excess
        Place where;
        int max_level = 1;  // an interval: the deepest split above the floor (AR7)
    };
    std::vector<Cand> cands;
    // The node of an edge with the larger weighted excess (the lower index on a tie).
    const auto heavier = [&](std::size_t a, std::size_t b) {
        const double x = w(a) * united.excess[a], y = w(b) * united.excess[b];
        return y > x || (y == x && b < a) ? b : a;
    };
    std::size_t sweep_n2 = 0;
    switch (cur.kind) {
        case AdaptMeshKind::tensor: {
            const mesh::TensorCells* tc = dev.cells();
            if (tc == nullptr) return std::unexpected(invalid("a tensor mesh needs its cells"));
            const int dim = tc->dimension();
            std::array<std::vector<double>, 3> lines, line_dofs;
            std::array<std::vector<double>, 3> share;
            std::array<std::vector<std::size_t>, 3> best;
            std::array<std::vector<char>, 3> seen;
            for (int a = 0; a < dim; ++a) {
                const auto ax = tc->axis(a);
                lines[static_cast<std::size_t>(a)].assign(ax.begin(), ax.end());
                line_dofs[static_cast<std::size_t>(a)].assign(ax.size(), 0.0);
                share[static_cast<std::size_t>(a)].assign(ax.size() - 1, 0.0);
                best[static_cast<std::size_t>(a)].assign(ax.size() - 1, 0);
                seen[static_cast<std::size_t>(a)].assign(ax.size() - 1, 0);
            }
            for (std::size_t i = 0; i < nn; ++i) {
                const auto idx = tc->index(static_cast<mesh::NodeId>(i));
                for (int a = 0; a < dim; ++a) line_dofs[static_cast<std::size_t>(a)][idx[static_cast<std::size_t>(a)]] += static_cast<double>(u(i));
            }
            for (std::size_t e = 0; e < m.edges().size(); ++e) {
                const auto& ed = m.edges()[e];
                const auto ia = tc->index(ed.first), ib = tc->index(ed.second);
                int axis = -1;
                bool aligned = true;
                for (int b = 0; b < dim; ++b) {
                    if (ia[static_cast<std::size_t>(b)] == ib[static_cast<std::size_t>(b)]) continue;
                    if (axis >= 0) aligned = false;
                    axis = b;
                }
                if (!aligned || axis < 0) continue;
                const auto A = static_cast<std::size_t>(axis);
                const std::size_t k = std::min(ia[A], ib[A]);
                const auto a = static_cast<std::size_t>(ed.first), b = static_cast<std::size_t>(ed.second);
                share[A][k] += std::max(w(a), w(b)) * ee[e];
                const std::size_t h = heavier(a, b);
                if (seen[A][k] == 0 || heavier(best[A][k], h) == h) best[A][k] = h;
                seen[A][k] = 1;
            }
            for (int a = 0; a < dim; ++a) {
                const auto A = static_cast<std::size_t>(a);
                const auto& L = lines[A];
                const double floor = axis_floor(o, L, cur.axis_protected[A]);
                const std::size_t per_line = nn / L.size();
                for (std::size_t k = 0; k + 1 < L.size(); ++k) {
                    if (!(share[A][k] > 0.0)) continue;
                    const int lv = levels_within(L[k + 1] - L[k], floor);
                    if (lv == 0 || inside_protected(L[k], L[k + 1], cur.axis_protected[A])) continue;
                    Cand c;
                    c.share = share[A][k];
                    c.node = best[A][k];
                    c.c.axis = a;
                    c.c.index = k;
                    c.c.location = pts[c.node];
                    c.c.location[A] = 0.5 * (L[k] + L[k + 1]);
                    c.c.added_nodes = per_line;
                    c.c.added_dofs = static_cast<std::size_t>(std::llround(0.5 * (line_dofs[A][k] + line_dofs[A][k + 1])));
                    c.where = Place{.axis = a, .lo = L[k], .hi = L[k + 1]};
                    c.max_level = lv;
                    cands.push_back(c);
                }
            }
            break;
        }
        case AdaptMeshKind::generated2d: {
            if (!cur.generated || cur.generated->simplices.points.size() != nn) {
                return std::unexpected(invalid("a generated 2D mesh needs its simplices, one point per node"));
            }
            std::vector<double> h = mesh::local_sizes(cur.generated->simplices);
            if (previous.size_field) {
                for (std::size_t i = 0; i < nn; ++i) h[i] = std::min(h[i], previous.size_field->at(cur.generated->simplices.points[i]));
            }
            for (std::size_t i = 0; i < nn; ++i) {
                if (!(united.excess[i] > 0.0 && w(i) > 0.0 && h[i] > floor2d)) continue;
                Cand c;
                c.share = w(i) * united.excess[i];
                c.node = i;
                c.c.index = i;
                c.c.location = pts[i];
                c.c.added_nodes = 3;  // a halved size quadruples the local density in 2D
                c.c.added_dofs = 3 * u(i);
                c.where = Place{.centre = pts[i], .radius = h[i], .region = static_cast<std::size_t>(regions[i])};
                cands.push_back(c);
            }
            break;
        }
        case AdaptMeshKind::swept: {
            if (!cur.section) return std::unexpected(invalid("a swept mesh needs its cross-section"));
            const std::size_t n2 = cur.section->points.size();
            if (n2 == 0 || nn % n2 != 0) return std::unexpected(invalid("a swept mesh's nodes are not whole lines"));
            sweep_n2 = n2;
            const std::size_t nl = nn / n2;
            const auto ax = static_cast<std::size_t>(cur.sweep_axis);
            std::vector<double> L(nl);
            for (std::size_t j = 0; j < nl; ++j) L[j] = pts[j * n2][ax];
            std::vector<double> h = mesh::local_sizes(*cur.section);
            if (previous.size_field) {
                for (std::size_t j = 0; j < n2; ++j) h[j] = std::min(h[j], previous.size_field->at(cur.section->points[j]));
            }
            for (std::size_t j = 0; j < n2; ++j) {
                if (!(h[j] > floor2d)) continue;
                double s = 0.0, dofs = 0.0;
                std::size_t heavy = j;
                for (std::size_t l = 0; l < nl; ++l) {
                    const std::size_t i = j + l * n2;
                    s += w(i) * united.excess[i];
                    dofs += static_cast<double>(u(i));
                    heavy = heavier(heavy, i);
                }
                if (!(s > 0.0)) continue;
                Cand c;
                c.share = s;
                c.node = heavy;
                c.c.index = j;
                c.c.location = pts[heavy];
                c.c.added_nodes = 3 * nl;
                c.c.added_dofs = static_cast<std::size_t>(3.0 * dofs);
                c.where = Place{.centre = pts[heavy], .radius = h[j], .skip_axis = cur.sweep_axis, .region = static_cast<std::size_t>(regions[heavy])};
                cands.push_back(c);
            }
            std::vector<double> share(nl - 1, 0.0);
            std::vector<std::size_t> best(nl - 1, 0);
            std::vector<char> seen(nl - 1, 0);
            for (std::size_t e = 0; e < m.edges().size(); ++e) {
                const auto a = static_cast<std::size_t>(m.edges()[e].first), b = static_cast<std::size_t>(m.edges()[e].second);
                if (b != a + n2) continue;
                const std::size_t l = a / n2;
                share[l] += std::max(w(a), w(b)) * ee[e];
                const std::size_t hv = heavier(a, b);
                if (seen[l] == 0 || heavier(best[l], hv) == hv) best[l] = hv;
                seen[l] = 1;
            }
            const double floor = axis_floor(o, L, cur.sweep_protected);
            for (std::size_t l = 0; l + 1 < nl; ++l) {
                if (!(share[l] > 0.0)) continue;
                const int lv = levels_within(L[l + 1] - L[l], floor);
                if (lv == 0 || inside_protected(L[l], L[l + 1], cur.sweep_protected)) continue;
                double dofs = 0.0;
                for (std::size_t j = 0; j < n2; ++j) dofs += 0.5 * static_cast<double>(u(j + l * n2) + u(j + (l + 1) * n2));
                Cand c;
                c.share = share[l];
                c.node = best[l];
                c.c.axis = cur.sweep_axis;
                c.c.index = l;
                c.c.location = pts[best[l]];
                c.c.location[ax] = 0.5 * (L[l] + L[l + 1]);
                c.c.added_nodes = n2;
                c.c.added_dofs = static_cast<std::size_t>(std::llround(dofs));
                c.where = Place{.axis = cur.sweep_axis, .lo = L[l], .hi = L[l + 1]};
                c.max_level = lv;
                cands.push_back(c);
            }
            break;
        }
    }
    out.count = cands.size();
    double total = 0.0;
    for (const Cand& c : cands) total += c.share;
    const double E = st.estimate;
    if (cands.empty() || !(total > 0.0) || !(std::isfinite(E) && E > 0.0)) return out;
    // Benefits and costs (ID4): the share of the estimate, the fitted order, the added nonzeros at
    // the mean degree and the candidate's unknowns per node, their fill at the current ratio.
    const double newton = st.newton > 0.0 ? st.newton : 1.0;
    const double rho = st.nnz > 0.0 ? st.fill / st.nnz : 0.0;
    const double W = st.nnz + st.fill;
    // One entry per candidate and level (AR7): level l leaves 2^-(l-1)p of its error to reduce by
    // (1 - 2^-p) and adds 2^(l-1) times the first level's nodes (an interval's next lines).
    struct Entry {
        std::size_t cand = 0;
        int level = 1;
        double benefit = 0.0, work = 0.0, cost = 0.0, dnnz = 0.0;
        double dnodes = 0.0, ddofs = 0.0;
    };
    std::vector<Entry> entries;
    for (std::size_t q = 0; q < cands.size(); ++q) {
        Cand& c = cands[q];
        c.c.region = static_cast<std::size_t>(regions[c.node]);
        c.c.indicator = united.dominant[c.node];
        c.c.order = order_of(fits, c.c.indicator, c.c.region);
        const double error = E * c.share / total;
        const double dn1 = static_cast<double>(c.c.added_nodes);
        const double s = static_cast<double>(c.c.added_dofs) / dn1;
        for (int l = 1; l <= c.max_level; ++l) {
            Entry e{.cand = q, .level = l};
            const double mult = std::exp2(l - 1);
            e.benefit = error * std::exp2(-(l - 1) * c.c.order) * (1.0 - std::exp2(-c.c.order));
            e.dnodes = dn1 * mult;
            e.ddofs = static_cast<double>(c.c.added_dofs) * mult;
            e.dnnz = (degree * e.dnodes + e.dnodes) * s * s;
            e.work = newton * (std::pow(W + e.dnnz * (1.0 + rho), 1.5) - std::pow(W, 1.5));
            e.cost = model.samples == 0
                         ? std::numeric_limits<double>::quiet_NaN()
                         : newton * (model.a * (std::pow(st.nnz + e.dnnz, model.b) - std::pow(st.nnz, model.b)) +
                                     model.c * rho * e.dnnz + model.d * e.ddofs) +
                               (model.mesh_per_node + model.transfer_per_node) * e.dnodes;
            entries.push_back(e);
        }
    }
    std::vector<std::size_t> order(entries.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::stable_sort(order, [&](std::size_t i, std::size_t j) {
        return entries[i].benefit * entries[j].work > entries[j].benefit * entries[i].work;
    });
    const double target = std::min(0.5, 0.5 * E);
    // The caps on the predicted DOFs, each's headroom over the measured inflation (AR7).
    double cap = st.dofs + (growth_cap - 1.0) * st.dofs / st.added_scale;
    const bool budget = o.intelligent.dof_budget > 0 || o.intelligent.seconds_budget > 0.0;
    if (o.intelligent.dof_budget > 0) {
        const double remaining = std::max(0.0, static_cast<double>(o.intelligent.dof_budget) - st.dofs);
        cap = std::min(cap, st.dofs + remaining / st.added_scale);
    }
    double predicted = E, pdofs = st.dofs, pnodes = st.nodes, pnnz = st.nnz;
    std::vector<int> taken_level(cands.size(), 0);
    std::vector<std::size_t> first_taken;  // candidates in the order of their first level
    for (const std::size_t q : order) {
        const Entry& e = entries[q];
        if (!budget && predicted <= target) break;
        if (!(e.benefit > 0.0)) break;
        if (taken_level[e.cand] != e.level - 1) continue;  // a level after the one before it
        if (pdofs + e.ddofs > cap) break;
        if (o.intelligent.seconds_budget > 0.0 && model.samples > 0) {
            const double nnz_next = pnnz + e.dnnz;
            const double t = model.pass_seconds(newton, nnz_next, rho * nnz_next, pdofs + e.ddofs, pnodes + e.dnodes);
            if (st.seconds_so_far + t > o.intelligent.seconds_budget) break;
        }
        taken_level[e.cand] = e.level;
        if (e.level == 1) first_taken.push_back(e.cand);
        AdaptCandidate& c = cands[e.cand].c;
        if (e.level == 1) {
            c.benefit = 0.0;
            c.work = 0.0;
            c.cost = 0.0;
            c.added_nodes = 0;
            c.added_dofs = 0;
            c.added_nonzeros = 0;
        }
        c.levels = e.level;
        c.benefit += e.benefit;
        c.work += e.work;
        c.cost += e.cost;
        c.added_nodes += static_cast<std::size_t>(std::llround(e.dnodes));
        c.added_dofs += static_cast<std::size_t>(std::llround(e.ddofs));
        c.added_nonzeros += static_cast<std::size_t>(std::llround(e.dnnz));
        predicted -= e.benefit;
        pdofs += e.ddofs;
        pnodes += e.dnodes;
        pnnz += e.dnnz;
    }
    out.taken = first_taken.size();
    if (out.taken == 0) {
        out.budget_spent = budget;
        return out;
    }
    out.predicted_relative = predicted;
    out.predicted_dofs = static_cast<std::size_t>(std::llround(pdofs));
    // The first taken are the listed (ID7) and the measured (AR7), with their local levels before
    // refinement for the ratios measured on the next pass.
    const AxisIndex index(pts);
    for (const std::size_t q : first_taken) {
        const Cand& c = cands[q];
        const bool listed = out.listed.size() < listed_candidates;
        if (listed) out.listed.push_back(c.c);
        if (out.pending.size() >= measured_candidates) continue;
        Pending p;
        p.where = c.where;
        p.kind = c.c.indicator;
        p.region = c.c.region;
        p.error = E * c.share / total;
        p.levels = c.c.levels;
        if (listed) p.listed = out.listed.size() - 1;
        p.before = local_level(dev, index, p.where, levels[static_cast<std::size_t>(p.kind)], omega);
        out.pending.push_back(p);
    }
    // The refinement: taken nodes marked (their size halved), taken intervals split to their levels.
    std::vector<char> mk(nn, 0);
    std::array<std::vector<int>, 3> interval_level;
    for (const std::size_t q : first_taken) {
        const AdaptCandidate& c = cands[q].c;
        if (c.axis >= 0) {
            auto& v = interval_level[static_cast<std::size_t>(c.axis)];
            if (v.size() <= c.index) v.resize(c.index + 1, 0);
            v[c.index] = c.levels;
        } else if (cur.kind == AdaptMeshKind::swept) {
            for (std::size_t i = c.index; i < nn; i += sweep_n2) mk[i] = 1;
        } else {
            mk[c.index] = 1;
        }
    }
    bool changed = false;
    auto r = refine_with(cur, previous, mk, [](std::size_t) { return false; }, nullptr, o, false, &changed);
    if (!r) return std::unexpected(r.error());
    MeshRefinement next = std::move(*r);
    if (cur.kind == AdaptMeshKind::tensor) {
        const mesh::TensorCells* tc = dev.cells();
        for (int a = 0; a < tc->dimension(); ++a) {
            const auto A = static_cast<std::size_t>(a);
            if (interval_level[A].empty()) continue;
            const auto ax = tc->axis(a);
            const std::vector<double> lines(ax.begin(), ax.end());
            auto split = split_levels(lines, interval_level[A], cur.axis_protected[A], o);
            if (!split) return std::unexpected(split.error());
            changed = changed || split->size() != lines.size();
            next.axes[A] = std::move(*split);
        }
    } else if (cur.kind == AdaptMeshKind::swept) {
        const auto A = static_cast<std::size_t>(cur.sweep_axis);
        if (!interval_level[A].empty()) {
            const std::size_t nl = nn / sweep_n2;
            std::vector<double> lines(nl);
            for (std::size_t j = 0; j < nl; ++j) lines[j] = pts[j * sweep_n2][A];
            auto split = split_levels(lines, interval_level[A], cur.sweep_protected, o);
            if (!split) return std::unexpected(split.error());
            changed = changed || split->size() != lines.size();
            std::vector<double> given;  // the generator adds the straddles' lines itself
            for (const double x : *split) {
                bool straddle_line = false;
                for (const mesh::AxisRange& p : cur.sweep_protected) straddle_line = straddle_line || x == p.lo || x == p.hi;
                if (!straddle_line) given.push_back(x);
            }
            next.sweep_lines = std::move(given);
        }
    }
    if (changed) out.next = std::move(next);
    return out;
}

struct Goals {
    std::vector<double> value;
    std::vector<double> resolution;
};

Goals read_goals(const results::Sweep& s, const AdaptOptions& o, const std::vector<std::size_t>& at) {
    Goals g;
    for (const AdaptGoal& goal : o.goals) {
        for (const std::size_t k : at) {
            const results::BiasPoint& p = s.points[k];
            const bool weighted = !p.weighted_current.empty();
            g.value.push_back(weighted ? p.weighted_current[goal.contact] : p.terminal_current[goal.contact]);
            const auto& res = weighted ? p.weighted_current_resolution : p.terminal_current_resolution;
            g.resolution.push_back(goal.contact < res.size() ? res[goal.contact] : 0.0);
        }
    }
    return g;
}

std::expected<mesh::PointLocator, base::Error> locator_of(const AdaptableMesh& m) {
    if (m.kind == AdaptMeshKind::tensor) {
        if (m.device.cells() == nullptr) return std::unexpected(invalid("a tensor mesh needs its cells"));
        return mesh::PointLocator::create(*m.device.cells());
    }
    if (!m.generated) return std::unexpected(invalid("a generated mesh needs its simplices"));
    return mesh::PointLocator::create(m.generated->simplices);
}

struct Solved {
    results::Sweep sweep;
    TransferReport transfer;
    double seconds_transfer = 0.0;
    double seconds_solve = 0.0;
};

// Transfers `from`'s first point onto `to` (when given) and solves the sweep there.
std::expected<Solved, base::Error> solve_on(const AdaptableMesh& to, const AdaptableMesh* from,
                                            const results::Sweep* from_sweep, std::span<const std::vector<double>> points,
                                            const BiasOptions& bias, const RunControl& control) {
    Solved out;
    std::optional<results::NodeFields> initial;
    BiasOptions opts = bias;
    if (from != nullptr && from_sweep != nullptr && !from_sweep->points.empty()) {
        const auto t0 = clock_type::now();
        auto loc = locator_of(*from);
        if (!loc) return std::unexpected(loc.error());
        const results::BiasPoint& p0 = from_sweep->points[0];
        auto st = transfer_state(from->device, *loc, p0.fields, p0.bands, to.device, &out.transfer);
        if (!st) return std::unexpected(st.error());
        initial = std::move(*st);
        opts.step_control.initial_bias_V = p0.bias_V;
        out.seconds_transfer = seconds_since(t0);
    }
    const auto t1 = clock_type::now();
    auto s = sweep_bias(to.device, points, opts, initial ? &*initial : nullptr, control);
    out.seconds_solve = seconds_since(t1);
    if (!s) return std::unexpected(s.error());
    if (s->stopped) return std::unexpected(*s->stopped);
    out.sweep = std::move(*s);
    return out;
}

// The verified error's distribution (AD12 c, AR2): per node of the coarse mesh M, the largest over
// the adaptation points of |psi' - psi| / V_T, |ln n' - ln n| and |ln p' - ln p|, the primes M''s
// solution interpolated at the node (over M' nodes of its own region). Solves only: no indicator,
// weight or prediction.
std::expected<IndicatorField, base::Error> difference_field(const AdaptableMesh& coarse, const results::Sweep& cs,
                                                            const AdaptableMesh& fine, const results::Sweep& fs,
                                                            const std::vector<std::size_t>& at) {
    auto loc = locator_of(fine);
    if (!loc) return std::unexpected(loc.error());
    const auto pts = coarse.device.mesh().points();
    const auto cr = coarse.device.node_region();
    const auto fr = fine.device.node_region();
    const double vt = base::thermal_voltage(coarse.device.temperature_K());
    IndicatorField f;
    f.excess.assign(pts.size(), 0.0);
    f.dominant.assign(pts.size(), Indicator::potential);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const mesh::Interpolant it = loc->locate(pts[i]);
        double sum = 0.0;
        for (int j = 0; j < it.count; ++j) {
            if (fr[static_cast<std::size_t>(it.nodes[static_cast<std::size_t>(j)])] == cr[i]) sum += it.weights[static_cast<std::size_t>(j)];
        }
        if (!(sum > 0.0)) continue;
        const bool insulator = coarse.device.is_insulator(static_cast<mesh::NodeId>(i));
        for (const std::size_t k : at) {
            const results::NodeFields& a = cs.points[k].fields;
            const results::NodeFields& b = fs.points[k].fields;
            double psi = 0.0, ln_n = 0.0, ln_p = 0.0;
            for (int j = 0; j < it.count; ++j) {
                const auto q = static_cast<std::size_t>(it.nodes[static_cast<std::size_t>(j)]);
                if (fr[q] != cr[i]) continue;
                const double w = it.weights[static_cast<std::size_t>(j)] / sum;
                psi += w * b.potential_V[q];
                if (!insulator) {
                    ln_n += w * std::log(b.n_cm3[q]);
                    ln_p += w * std::log(b.p_cm3[q]);
                }
            }
            const std::array<double, 3> d{std::abs(psi - a.potential_V[i]) / vt,
                                          insulator ? 0.0 : std::abs(ln_n - std::log(a.n_cm3[i])),
                                          insulator ? 0.0 : std::abs(ln_p - std::log(a.p_cm3[i]))};
            const std::array<Indicator, 3> kind{Indicator::potential, Indicator::electron_density, Indicator::hole_density};
            for (std::size_t m = 0; m < 3; ++m) {
                if (d[m] > f.excess[i]) {
                    f.excess[i] = d[m];
                    f.dominant[i] = kind[m];
                }
            }
        }
    }
    return f;
}

}  // namespace

std::expected<AdaptResult, base::Error> adapt(const DeviceBuilder& build, std::span<const std::vector<double>> points,
                                              const BiasOptions& bias, const AdaptOptions& o, const RunControl& control) {
    if (points.empty()) return std::unexpected(invalid("adaptation needs bias points"));
    if (o.goals.empty()) return std::unexpected(invalid("adaptation needs at least one goal"));
    if (!(o.tolerance > 0.0 && o.tolerance < 1.0)) return std::unexpected(invalid("the tolerance must be in (0, 1)"));
    if (!(o.theta > 0.0 && o.theta <= 1.0)) return std::unexpected(invalid("theta must be in (0, 1]"));
    if (!(std::isfinite(o.grading) && o.grading >= 1.0)) return std::unexpected(invalid("the grading must be >= 1"));
    if (!(std::isfinite(o.min_size) && o.min_size >= 0.0)) return std::unexpected(invalid("min_size must be >= 0"));
    if (o.max_passes < 1) return std::unexpected(invalid("max_passes must be at least 1"));
    const bool adjoint = o.weighting == Weighting::adjoint;
    if (adjoint) {
        if (auto e = check_adjoint(bias)) return std::unexpected(std::move(*e));
    }
    std::vector<std::size_t> at = o.at;
    if (at.empty()) {
        at.push_back(0);
        if (points.size() > 1) at.push_back(points.size() - 1);
    }
    for (const std::size_t k : at) {
        if (k >= points.size()) return std::unexpected(invalid("an adaptation point is outside the sweep"));
    }

    AdaptResult result;
    MeshRefinement refinement;
    auto t0 = clock_type::now();
    auto built = build(refinement);
    if (!built) return std::unexpected(built.error());
    for (const AdaptGoal& g : o.goals) {
        if (g.contact >= built->device.contacts().size()) return std::unexpected(invalid("a goal names no contact"));
        const auto kind = built->device.contacts()[g.contact].kind;
        if (kind == device::ContactKind::gate || kind == device::ContactKind::electrode) {
            return std::unexpected(invalid("a goal's contact is a gate or an electrode, which carries no current"));
        }
        if (adjoint && kind != device::ContactKind::ohmic) {
            return std::unexpected(invalid("adjoint goal weighting needs ohmic goal contacts"));
        }
    }
    std::vector<std::size_t> goal_contacts;
    for (const AdaptGoal& g : o.goals) goal_contacts.push_back(g.contact);
    // AD6's estimate of goal j on a pass (goal-major, as the goals): the adjoint's (AR3), or the
    // last change without weighting (infinity on the first pass).
    const auto estimate_of = [&](const AdaptPass& rep, std::size_t j) {
        if (adjoint) return rep.estimate[j];
        return rep.change.empty() ? std::numeric_limits<double>::infinity() : rep.change[j];
    };
    // M' for a pass: the uniform refinement of its mesh, built (the estimate's and AD11's mesh).
    const auto build_uniform = [&](const AdaptableMesh& m,
                                   const MeshRefinement& r) -> std::expected<AdaptableMesh, base::Error> {
        auto uniform = next_refinement(m, r, {}, nullptr, o.indicators, o, true);
        if (!uniform) return std::unexpected(uniform.error());
        return build(*uniform);
    };
    // The adjoint estimate per goal and point and AD3's weights omega (per node, the largest over
    // goals and points of the semiconductor-normalized weight), on M' (AR3).
    struct AdjointPart {
        std::optional<AdaptableMesh> fine;
        std::vector<double> estimate;
        std::vector<double> omega;
        bool fell_back = false;
        double seconds_factorize = 0.0;
        double seconds_transposed = 0.0;
        double seconds = 0.0;
    };
    const auto adjoint_part = [&](const AdaptableMesh& m, const MeshRefinement& r,
                                  const results::Sweep& s) -> std::expected<AdjointPart, base::Error> {
        AdjointPart out;
        const auto te = clock_type::now();
        auto f = build_uniform(m, r);
        if (!f) return std::unexpected(f.error());
        out.fine = std::move(*f);
        auto loc = locator_of(m);
        if (!loc) return std::unexpected(loc.error());
        const std::size_t nn = m.device.mesh().node_count();
        out.omega.assign(nn, 0.0);
        out.estimate.assign(goal_contacts.size() * at.size(), 0.0);
        for (std::size_t a = 0; a < at.size(); ++a) {
            auto e = estimate_goals(m.device, *loc, s.points[at[a]], out.fine->device, bias, goal_contacts);
            if (!e) return std::unexpected(e.error());
            out.fell_back = out.fell_back || e->fell_back;
            out.seconds_factorize += e->seconds_factorize;
            out.seconds_transposed = std::max(out.seconds_transposed, e->seconds_transposed);
            for (std::size_t gi = 0; gi < goal_contacts.size(); ++gi) {
                out.estimate[gi * at.size() + a] = e->estimate[gi];
                // Normalized so the semiconductor nodes' weights sum to 1 (AD3).
                double total = 0.0;
                for (std::size_t i = 0; i < nn; ++i) {
                    if (!m.device.is_insulator(static_cast<mesh::NodeId>(i))) total += e->node_weight[gi][i];
                }
                if (!(total > 0.0)) continue;
                for (std::size_t i = 0; i < nn; ++i) {
                    if (m.device.is_insulator(static_cast<mesh::NodeId>(i))) continue;
                    out.omega[i] = std::max(out.omega[i], e->node_weight[gi][i] / total);
                }
            }
        }
        out.seconds = seconds_since(te);
        return out;
    };
    // The largest over goals and points of a per-goal quantity over its allowance tol max(|G|, r).
    const auto relative = [&](const std::vector<double>& x, const Goals& g) {
        double r = 0.0;
        for (std::size_t j = 0; j < x.size(); ++j) {
            r = std::max(r, x[j] / (o.tolerance * std::max(std::abs(g.value[j]), g.resolution[j])));
        }
        return r;
    };

    // The decision layer before pass 0 (ID3, ID2): the adaptation points from a coarse sweep on the
    // start mesh, and the start refined toward the a priori prediction.
    const bool intelligent = o.policy == Policy::intelligent;
    std::optional<Solved> coarse;
    double coarse_relative = std::numeric_limits<double>::infinity();
    MeshSource source = MeshSource::start;
    if (intelligent) {
        if (!(o.intelligent.miss_margin >= 0.0 && o.intelligent.miss_margin < 1.0)) {
            return std::unexpected(invalid("the miss margin must be in [0, 1)"));
        }
        if (o.at.empty() && o.intelligent.max_points > 0) {
            const auto tc = clock_type::now();
            auto c = solve_on(*built, nullptr, nullptr, points, bias, control);
            if (!c) return std::unexpected(c.error());
            std::vector<IndicatorField> fields;
            for (const results::BiasPoint& p : c->sweep.points) {
                auto f = compute_indicators(built->device, &p.fields, &p.rates, o.indicators);
                if (!f) return std::unexpected(f.error());
                fields.push_back(std::move(*f));
            }
            // Eligible: every goal above its current resolution at the point (AR5).
            std::vector<char> eligible(c->sweep.points.size(), 1);
            for (std::size_t k = 0; k < eligible.size(); ++k) {
                const results::BiasPoint& p = c->sweep.points[k];
                const bool weighted = !p.weighted_current.empty();
                const auto& cur = weighted ? p.weighted_current : p.terminal_current;
                const auto& res = weighted ? p.weighted_current_resolution : p.terminal_current_resolution;
                for (const AdaptGoal& goal : o.goals) {
                    const double r = goal.contact < res.size() ? res[goal.contact] : 0.0;
                    if (!(std::abs(cur[goal.contact]) > r)) eligible[k] = 0;
                }
            }
            auto chosen = choose_points(fields, o.intelligent.max_points, eligible);
            if (!chosen) return std::unexpected(chosen.error());
            at = std::move(*chosen);
            if (adjoint) {
                auto part = adjoint_part(*built, refinement, c->sweep);
                if (!part) return std::unexpected(part.error());
                coarse_relative = relative(part->estimate, read_goals(c->sweep, o, at));
            }
            result.coarse_nodes = built->device.mesh().node_count();
            result.coarse_seconds = seconds_since(tc);
            coarse = std::move(*c);
        }
        const std::size_t start_dofs = dofs_of(built->device, bias);
        for (int round = 0; round < 12; ++round) {
            SizeTargets t;
            if (o.intelligent.prediction) {
                t = o.intelligent.prediction(built->device);
            } else {
                auto p = predict_sizes(built->device, points, bias);
                if (!p) return std::unexpected(p.error());
                if (round == 0) {
                    result.junction_points = p->junction_points;
                    result.interface_points = p->interface_points;
                    result.predicted_nodes = p->predicted_nodes;
                }
                t = std::move(p->targets);
            }
            const std::size_t nn = built->device.mesh().node_count();
            if (t.size.size() != nn || t.normal.size() != nn) {
                return std::unexpected(invalid("a prediction must give one target and normal per node"));
            }
            // Per edge, how far it exceeds its target (projected length over the smaller target of its
            // ends; 0 when it does not); per node, the largest over its edges.
            const device::Device& dev = built->device;
            const auto& edges = dev.mesh().edges();
            const auto dpts = dev.mesh().points();
            std::vector<double> ratio(edges.size(), 0.0), node_ratio(nn, 0.0);
            for (std::size_t e = 0; e < edges.size(); ++e) {
                if (!edge_exceeds(dev, t, e)) continue;
                const auto a = static_cast<std::size_t>(edges[e].first), b = static_cast<std::size_t>(edges[e].second);
                double best = 0.0;
                for (const std::size_t k : {a, b}) {
                    if (!std::isfinite(t.size[k])) continue;
                    const auto& n = t.normal[k];
                    double along = edges[e].length;
                    if (n[0] != 0.0 || n[1] != 0.0 || n[2] != 0.0) {
                        along = std::abs((dpts[b][0] - dpts[a][0]) * n[0] + (dpts[b][1] - dpts[a][1]) * n[1] +
                                         (dpts[b][2] - dpts[a][2]) * n[2]);
                    }
                    best = std::max(best, along / t.size[k]);
                }
                ratio[e] = best;
                node_ratio[a] = std::max(node_ratio[a], best);
                node_ratio[b] = std::max(node_ratio[b], best);
            }
            // The refinement taking the edges (and nodes) at or above rho, built; nullopt when it
            // refines nothing or exceeds a budget: max_nodes or ID6's growth cap, the prediction being
            // pass 0's proposal (ID1, AR5).
            struct Round {
                MeshRefinement refinement;
                AdaptableMesh mesh;
            };
            bool over = false;
            const auto try_round = [&](double rho) -> std::expected<std::optional<Round>, base::Error> {
                std::vector<double> sizes = t.size;
                for (std::size_t i = 0; i < nn; ++i) {
                    if (!(node_ratio[i] >= rho && node_ratio[i] > 0.0)) sizes[i] = std::numeric_limits<double>::infinity();
                }
                bool changed = false;
                auto r = refine_with(*built, refinement, std::vector<char>(nn, 0),
                                     [&](std::size_t e) { return ratio[e] > 0.0 && ratio[e] >= rho; }, &sizes, o, false,
                                     &changed);
                if (!r) return std::unexpected(r.error());
                over = false;
                if (!changed) return std::optional<Round>{};
                auto nb = build(*r);
                if (!nb) return std::unexpected(nb.error());
                const std::size_t nd = dofs_of(nb->device, bias);
                if ((o.max_nodes > 0 && nb->device.mesh().node_count() > o.max_nodes) ||
                    static_cast<double>(nd) > growth_cap * static_cast<double>(start_dofs) ||
                    (o.intelligent.dof_budget > 0 && nd > o.intelligent.dof_budget)) {
                    over = true;
                    return std::optional<Round>{};
                }
                return std::optional<Round>{Round{std::move(*r), std::move(*nb)}};
            };
            auto full = try_round(0.0);
            if (!full) return std::unexpected(full.error());
            std::optional<Round> taken = std::move(*full);
            bool last = false;
            if (!taken && over) {
                // Over a budget: the edges that exceed their targets most, as many as fit (bisection on
                // the sorted ratios, the largest set within the budget); this round is the last.
                last = true;
                std::vector<double> sorted;
                for (const double x : ratio) {
                    if (x > 0.0) sorted.push_back(x);
                }
                std::ranges::sort(sorted, std::greater<>());
                sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
                std::size_t lo = 0, hi = sorted.size();  // the answer takes sorted[0 .. k), k in [lo, hi)
                while (lo + 1 < hi) {
                    const std::size_t mid = (lo + hi) / 2;
                    auto attempt = try_round(sorted[mid - 1]);
                    if (!attempt) return std::unexpected(attempt.error());
                    if (*attempt) {
                        taken = std::move(*attempt);
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
            }
            if (!taken) break;
            built = std::move(taken->mesh);
            refinement = std::move(taken->refinement);
            ++result.prediction_rounds;
            source = MeshSource::prediction;
            if (last) break;
        }
    }
    bool fallen_back = false;
    int consecutive_misses = 0;
    double previous_relative = coarse_relative;
    double previous_change = std::numeric_limits<double>::infinity();
    const auto fall_back = [&](FallbackTrigger t, int pass) {
        if (!intelligent || fallen_back) return;
        fallen_back = true;
        result.fallback = t;
        result.fallback_pass = pass;
    };

    // Pass 0's mesh time: the start's build and the prediction's (not the coarse sweep).
    double seconds_mesh = seconds_since(t0) - result.coarse_seconds;
    std::optional<AdaptableMesh> previous_mesh;
    results::Sweep previous_sweep;
    std::vector<double> previous_goal;
    // ID5: the passes' cost samples and the realized ratios so far; the candidates the last allocation
    // took, measured on the next pass.
    std::vector<CostSample> cost_samples;
    std::vector<OrderSample> order_samples;
    std::vector<OrderFit> fits;
    std::vector<Pending> pending;
    double seconds_so_far = result.coarse_seconds;
    double added_scale = 1.0;  // the last allocation's built over predicted added DOFs (>= 1)

    // AD11 on the mesh of a pass: M' solved from the transferred state; the verified error per goal.
    struct Verification {
        std::vector<double> goal;
        std::vector<double> error;
        std::vector<double> order;
        std::size_t nodes = 0;
        double seconds = 0.0;
        bool accepts = true;  // every goal's estimate and verified error within its allowance
        IndicatorField difference;  // the verified error's distribution on M (AR2)
    };
    std::optional<Verification> verification;  // of the current mesh, when run
    // `prebuilt`: M' when the pass built it for the estimate (the same refinement, so the same mesh).
    const auto verify_pass = [&](const AdaptableMesh& m, const MeshRefinement& r, const results::Sweep& s,
                                 const AdaptPass& rep,
                                 const AdaptableMesh* prebuilt) -> std::expected<Verification, base::Error> {
        Verification v;
        const auto tv = clock_type::now();
        std::optional<AdaptableMesh> own;
        if (prebuilt == nullptr) {
            auto f = build_uniform(m, r);
            if (!f) return std::unexpected(f.error());
            own = std::move(*f);
        }
        const AdaptableMesh* fine = prebuilt != nullptr ? prebuilt : &*own;
        auto solved = solve_on(*fine, &m, &s, points, bias, control);
        if (!solved) return std::unexpected(solved.error());
        v.seconds = seconds_since(tv);
        v.nodes = fine->device.mesh().node_count();
        v.goal = read_goals(solved->sweep, o, at).value;
        auto diff = difference_field(m, s, *fine, solved->sweep, at);
        if (!diff) return std::unexpected(diff.error());
        v.difference = std::move(*diff);
        for (std::size_t j = 0; j < v.goal.size(); ++j) {
            const double d_fine = std::abs(v.goal[j] - rep.goal[j]);
            double p = 1.0;
            if (!rep.change.empty() && d_fine > 0.0 && rep.change[j] > 0.0) {
                p = std::clamp(std::log2(rep.change[j] / d_fine), 0.5, 3.0);
            }
            v.order.push_back(p);
            v.error.push_back(d_fine / (std::exp2(p) - 1.0));
            const double allowed = o.tolerance * std::max(std::abs(rep.goal[j]), rep.resolution[j]);
            if (!(estimate_of(rep, j) <= allowed && v.error.back() <= allowed)) v.accepts = false;
        }
        return v;
    };

    for (int pass = 0;; ++pass) {
        auto solved = [&]() -> std::expected<Solved, base::Error> {
            // ID3's coarse sweep is pass 0's solve when the prediction refined nothing.
            if (pass == 0 && coarse && result.prediction_rounds == 0) return std::move(*coarse);
            return solve_on(*built, previous_mesh ? &*previous_mesh : nullptr, previous_mesh ? &previous_sweep : nullptr,
                            points, bias, control);
        }();
        if (!solved) return std::unexpected(solved.error());
        AdaptPass rep;
        rep.nodes = built->device.mesh().node_count();
        rep.edges = built->device.mesh().edges().size();
        rep.transfer = solved->transfer;
        rep.seconds_mesh = seconds_mesh;
        rep.seconds_transfer = solved->seconds_transfer;
        rep.seconds_solve = solved->seconds_solve;
        for (const results::BiasPoint& p : solved->sweep.points) {
            rep.newton_iterations += static_cast<int>(p.convergence.iterations.size());
            rep.matrix_nonzeros = std::max(rep.matrix_nonzeros, p.convergence.matrix_nonzeros);
            rep.factor_fill = std::max(rep.factor_fill, p.convergence.factor_fill.value_or(0));
        }
        rep.dofs = dofs_of(built->device, bias);
        rep.mesh_digest = mesh::graph_digest(built->device.mesh());
        if (o.on_mesh) {
            if (auto r = o.on_mesh(static_cast<std::size_t>(pass), *built); !r) return std::unexpected(r.error());
        }
        // ID5's cost model: this pass predicted from the earlier passes (their last Newton count, this
        // mesh's size), then refitted with it.
        const CostSample sample{.newton_iterations = static_cast<double>(rep.newton_iterations),
                                .matrix_nonzeros = static_cast<double>(rep.matrix_nonzeros),
                                .fill = static_cast<double>(rep.factor_fill),
                                .dofs = static_cast<double>(rep.dofs),
                                .nodes = static_cast<double>(rep.nodes),
                                .seconds_solve = rep.seconds_solve,
                                .seconds_mesh = rep.seconds_mesh,
                                .seconds_transfer = rep.seconds_transfer};
        if (intelligent && !cost_samples.empty()) {
            const CostModel before = fit_cost_model(cost_samples);
            rep.predicted_seconds = before.pass_seconds(cost_samples.back().newton_iterations, sample.matrix_nonzeros,
                                                        sample.fill, sample.dofs, sample.nodes);
            rep.predicted_seconds_nodes = before.pass_seconds_nodes(sample.nodes);
        }
        cost_samples.push_back(sample);
        if (intelligent) rep.cost_model = fit_cost_model(cost_samples);
        const Goals g = read_goals(solved->sweep, o, at);
        rep.goal = g.value;
        rep.resolution = g.resolution;
        // The adjoint estimate and weights (AD3, AD6; AR3): on M', per adaptation point.
        std::optional<AdaptableMesh> fine_mesh;
        std::vector<double> omega;  // AD3's node weight: per node, the largest over goals and points
        if (adjoint) {
            auto part = adjoint_part(*built, refinement, solved->sweep);
            if (!part) return std::unexpected(part.error());
            fine_mesh = std::move(part->fine);
            omega = std::move(part->omega);
            rep.estimate = std::move(part->estimate);
            rep.estimate_fell_back = part->fell_back;
            rep.seconds_estimate_factorize = part->seconds_factorize;
            rep.seconds_transposed = part->seconds_transposed;
            rep.estimate_nodes = fine_mesh->device.mesh().node_count();
            rep.seconds_estimate = part->seconds;
        }
        seconds_so_far += rep.seconds_mesh + rep.seconds_transfer + rep.seconds_solve + rep.seconds_estimate;
        // ID6's state on this pass, kept for a re-allocation when the built proposal is over a cap.
        bool allocated = false;
        std::optional<PassState> alloc_state;
        std::optional<std::array<std::vector<double>, indicator_count>> alloc_levels;
        // Indicators at every adaptation point, united (AQ3): computed once per pass when needed.
        std::optional<EdgeIndicators> edges;
        std::optional<IndicatorField> united;
        const auto indicators = [&]() -> std::expected<void, base::Error> {
            if (united) return {};
            std::vector<IndicatorField> fields;
            for (const std::size_t k : at) {
                const results::BiasPoint& p = solved->sweep.points[k];
                auto f = compute_indicators(built->device, &p.fields, &p.rates, o.indicators);
                if (!f) return std::unexpected(f.error());
                fields.push_back(std::move(*f));
                auto e = edge_indicators(built->device, &p.fields, &p.rates, o.indicators);
                if (!e) return std::unexpected(e.error());
                if (!edges) {
                    edges = std::move(*e);
                } else {
                    for (std::size_t kind = 0; kind < indicator_count; ++kind) {
                        for (std::size_t j = 0; j < edges->raw[kind].size(); ++j) {
                            edges->raw[kind][j] = std::max(edges->raw[kind][j], e->raw[kind][j]);
                        }
                    }
                }
            }
            auto u = unite(fields);
            if (!u) return std::unexpected(u.error());
            united = std::move(*u);
            return {};
        };
        // ID5: the last allocation's candidates measured on this mesh (their local level after over
        // before), and the orders refitted. Decisions only (OA6).
        if (!pending.empty()) {
            if (auto ok = indicators(); !ok) return std::unexpected(ok.error());
            const auto levels = node_levels(built->device.mesh(), *edges, o.indicators);
            const AxisIndex index(built->device.mesh().points());
            AdaptPass& last = result.passes.back();
            for (const Pending& p : pending) {
                const double after = local_level(built->device, index, p.where, levels[static_cast<std::size_t>(p.kind)], omega);
                if (!(std::isfinite(after) && std::isfinite(p.before) && p.before > 0.0)) continue;
                const double ratio = after / p.before;
                // Per halving: the ratio's levels-th root.
                order_samples.push_back(
                    OrderSample{.kind = p.kind, .region = p.region, .ratio = std::pow(ratio, 1.0 / p.levels)});
                if (p.listed < last.candidates.size()) {
                    last.candidates[p.listed].realized_ratio = ratio;
                    last.candidates[p.listed].realized = p.error * (1.0 - ratio);
                }
            }
            pending.clear();
            fits = fit_orders(order_samples);
        }
        rep.orders = fits;
        bool converged = pass > 0;
        if (pass > 0) {
            rep.change.resize(g.value.size());
            for (std::size_t j = 0; j < g.value.size(); ++j) {
                rep.change[j] = std::abs(g.value[j] - previous_goal[j]);
                const double allowed = o.tolerance * std::max(std::abs(g.value[j]), g.resolution[j]);
                if (rep.change[j] > allowed) converged = false;
                if (adjoint && !(rep.estimate[j] <= allowed)) converged = false;  // AD6: and the estimate
            }
        }
        // The guard (ID1) and AD12's triggers a and b: decisions only, after the solve; nothing
        // here enters the stop rule, AD11 or the certificate (OA6).
        const double inf = std::numeric_limits<double>::infinity();
        const double change_relative = rep.change.empty() ? inf : relative(rep.change, g);
        rep.source = source;
        rep.relative_estimate = adjoint ? relative(rep.estimate, g) : change_relative;
        if (intelligent && !fallen_back) {
            const bool proposal = source == MeshSource::prediction || source == MeshSource::proposal;
            if (proposal && std::isfinite(previous_relative)) {
                rep.miss = !(rep.relative_estimate <= (1.0 - o.intelligent.miss_margin) * previous_relative);
                if (rep.miss) {
                    ++result.misses;
                    if (++consecutive_misses >= 2) fall_back(FallbackTrigger::two_misses, pass);
                } else {
                    consecutive_misses = 0;
                }
            }
            if (pass >= 2 && std::isfinite(previous_relative) && std::isfinite(previous_change) &&
                !(rep.relative_estimate < previous_relative) && !(change_relative < previous_change)) {
                fall_back(FallbackTrigger::no_progress, pass);
            }
        }
        std::optional<AdaptStop> stop;
        verification.reset();
        if (converged && o.verify) {
            // A proposed stop, subject to AD11; a rejected stop goes on refining (AR2).
            auto v = verify_pass(*built, refinement, solved->sweep, rep, fine_mesh ? &*fine_mesh : nullptr);
            if (!v) return std::unexpected(v.error());
            verification = std::move(*v);
            if (verification->accepts) {
                stop = AdaptStop::converged;
            } else {
                rep.rejected = true;
                ++result.rejected_stops;
                fall_back(FallbackTrigger::rejected_stop, pass);  // AD12 c
            }
        } else if (converged) {
            stop = AdaptStop::converged;
        }
        if (!stop && pass + 1 >= o.max_passes) stop = AdaptStop::max_passes;
        std::optional<MeshRefinement> next;
        MeshSource next_source = MeshSource::marks;
        if (!stop) {
            // The indicators marked by AD3: weighted by the adjoint's omega, or feature-only (w = 1).
            if (auto ok = indicators(); !ok) return std::unexpected(ok.error());
            rep.nodes_with_excess = united->nodes_with_excess;
            // After a rejected stop the marks come from the verified error's distribution (AD12 c,
            // AR2): the indicators found nothing more to refine where the error is. That distribution
            // is marked unweighted: it is already the solves' own difference (OA6).
            const bool by_difference = rep.rejected;
            const std::span<const double> weights =
                adjoint && !by_difference ? std::span<const double>(omega) : std::span<const double>();
            auto mk = mark_nodes(built->device, by_difference ? verification->difference : *united, weights, o.theta);
            if (!mk) return std::unexpected(mk.error());
            rep.marking = *mk;
            rep.marked_by_difference = by_difference;
            // The decision layer's proposal (ID1), unless the run fell back, this pass missed or a stop
            // was rejected (AD3's marks then). A proposal that refines nothing is a miss itself.
            const bool proposing = intelligent && !fallen_back && !rep.miss && !rep.rejected;
            if (proposing && adjoint && !o.intelligent.prediction && std::isfinite(rep.relative_estimate)) {
                // ID6's allocation of ID4's candidates (AR7), with the adjoint's estimate and weights
                // (ID4's benefit is in the goal's units through them; without weighting the last
                // change is no error to reduce, and AD3's marks are used); none taken is an empty
                // proposal (a miss), or in budget mode the budget spent (the run stops).
                const PassState st{.estimate = rep.relative_estimate,
                                   .newton = static_cast<double>(rep.newton_iterations),
                                   .nnz = static_cast<double>(rep.matrix_nonzeros),
                                   .fill = static_cast<double>(rep.factor_fill),
                                   .dofs = static_cast<double>(rep.dofs),
                                   .nodes = static_cast<double>(rep.nodes),
                                   .seconds_so_far = seconds_so_far,
                                   .added_scale = added_scale};
                alloc_levels = node_levels(built->device.mesh(), *edges, o.indicators);
                alloc_state = st;
                auto al = allocate(*built, refinement, *united, *edges, *alloc_levels, omega, st, rep.cost_model, fits, bias, o);
                if (!al) return std::unexpected(al.error());
                rep.candidate_count = al->count;
                rep.candidates_taken = al->taken;
                rep.candidates = std::move(al->listed);
                rep.predicted_relative = al->predicted_relative;
                rep.predicted_dofs = al->predicted_dofs;
                if (al->next) {
                    next = std::move(al->next);
                    next_source = MeshSource::proposal;
                    pending = std::move(al->pending);
                    allocated = true;
                } else if (al->budget_spent) {
                    stop = AdaptStop::max_nodes;
                } else {
                    rep.empty_proposal = true;
                    ++result.misses;
                    if (++consecutive_misses >= 2) fall_back(FallbackTrigger::two_misses, pass);
                }
            }
            if (proposing && o.intelligent.prediction) {
                const SizeTargets t = o.intelligent.prediction(built->device);
                const std::size_t nn = built->device.mesh().node_count();
                if (t.size.size() != nn || t.normal.size() != nn) {
                    return std::unexpected(invalid("a prediction must give one target and normal per node"));
                }
                bool changed = false;
                const device::Device& dev = built->device;
                auto r = refine_with(*built, refinement, std::vector<char>(nn, 0),
                                     [&](std::size_t e) { return edge_exceeds(dev, t, e); }, &t.size, o, false, &changed);
                if (!r) return std::unexpected(r.error());
                if (changed) {
                    next = std::move(*r);
                    next_source = MeshSource::proposal;
                } else {
                    rep.empty_proposal = true;
                    ++result.misses;
                    if (++consecutive_misses >= 2) fall_back(FallbackTrigger::two_misses, pass);
                }
            }
            if (next || stop) {
                // The proposal, or the budget spent.
            } else if (by_difference && mk->count > 0) {
                auto r = next_refinement(*built, refinement, mk->marked, nullptr, o.indicators, o, false);
                if (!r) return std::unexpected(r.error());
                next = std::move(*r);
                next_source = MeshSource::difference;
            } else if (mk->count == 0 && rep.rejected) {
                // A rejected stop with nothing marked: the uniform refinement (M') is the next mesh.
                auto r = next_refinement(*built, refinement, {}, nullptr, o.indicators, o, true);
                if (!r) return std::unexpected(r.error());
                next = std::move(*r);
                next_source = MeshSource::uniform;
            } else if (mk->count == 0) {
                stop = AdaptStop::nothing_marked;
            } else {
                auto r = next_refinement(*built, refinement, mk->marked, &*edges, o.indicators, o, false);
                if (!r) return std::unexpected(r.error());
                next = std::move(*r);
            }
        }
        result.passes.push_back(std::move(rep));
        if (next) {
            t0 = clock_type::now();
            auto nb = build(*next);
            if (!nb) return std::unexpected(nb.error());
            std::size_t next_dofs = dofs_of(nb->device, bias);
            // ID6's caps hold on the built mesh (AR7): the allocation's predicted DOFs miss grading's
            // added lines and sizes, so a proposal built over the pass's growth cap or the budget is
            // allocated again with the measured inflation, at most three times.
            AdaptPass& made = result.passes.back();
            const auto inflation = [&](std::size_t built_dofs) {
                return made.predicted_dofs > made.dofs
                           ? std::max(1.0, static_cast<double>(built_dofs - std::min(built_dofs, made.dofs)) /
                                               static_cast<double>(made.predicted_dofs - made.dofs))
                           : 1.0;
            };
            const auto over = [&](std::size_t d) {
                return static_cast<double>(d) > growth_cap * static_cast<double>(made.dofs) ||
                       (o.intelligent.dof_budget > 0 && d > o.intelligent.dof_budget);
            };
            for (int retry = 0; allocated && retry < 3 && over(next_dofs); ++retry) {
                alloc_state->added_scale *= inflation(next_dofs);
                auto al = allocate(*built, refinement, *united, *edges, *alloc_levels, omega, *alloc_state, made.cost_model,
                                   fits, bias, o);
                if (!al) return std::unexpected(al.error());
                if (!al->next) break;
                made.candidate_count = al->count;
                made.candidates_taken = al->taken;
                made.candidates = std::move(al->listed);
                made.predicted_relative = al->predicted_relative;
                made.predicted_dofs = al->predicted_dofs;
                pending = std::move(al->pending);
                next = std::move(al->next);
                auto again = build(*next);
                if (!again) return std::unexpected(again.error());
                nb = std::move(again);
                next_dofs = dofs_of(nb->device, bias);
            }
            seconds_mesh = seconds_since(t0);
            if (allocated) added_scale = std::max(added_scale, inflation(next_dofs));  // for the next allocation
            if ((o.max_nodes > 0 && nb->device.mesh().node_count() > o.max_nodes) ||
                (intelligent && o.intelligent.dof_budget > 0 && next_dofs > o.intelligent.dof_budget)) {
                stop = AdaptStop::max_nodes;
            } else {
                previous_mesh = std::move(*built);
                previous_sweep = std::move(solved->sweep);
                previous_goal = g.value;
                previous_relative = result.passes.back().relative_estimate;
                previous_change = change_relative;
                built = std::move(*nb);
                refinement = std::move(*next);
                source = next_source;
                continue;
            }
        }
        // Stopped: the result is this pass.
        result.stop = *stop;
        result.sweep = std::move(solved->sweep);
        break;
    }

    // The certificate (AD13): from the final mesh's verification when one was run (a goal stop, or a
    // rejected stop at the last pass), else unverified.
    const AdaptPass& last = result.passes.back();
    const std::size_t ng = last.goal.size();
    std::vector<double> verified(ng, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> order(ng, 1.0);
    const bool verify = verification.has_value();
    if (verify) {
        result.verification_seconds = verification->seconds;
        result.verification_nodes = verification->nodes;
        result.verification_goal = verification->goal;
        verified = verification->error;
        order = verification->order;
    }
    for (std::size_t gi = 0; gi < o.goals.size(); ++gi) {
        for (std::size_t a = 0; a < at.size(); ++a) {
            const std::size_t j = gi * at.size() + a;
            GoalCertificate c;
            c.goal = gi;
            c.point = at[a];
            c.value = last.goal[j];
            c.estimate = estimate_of(last, j);
            c.verified_error = verified[j];
            c.order = order[j];
            c.allowed = o.tolerance * std::max(std::abs(c.value), last.resolution[j]);
            if (!verify) {
                c.status = CertificateStatus::unverified;
            } else {
                c.status = c.estimate <= c.allowed && c.verified_error <= c.allowed ? CertificateStatus::certified
                                                                                    : CertificateStatus::not_certified;
            }
            result.certificates.push_back(c);
        }
    }
    result.points = at;
    result.policy = o.policy;
    result.mesh = std::move(*built);
    return result;
}

}  // namespace NiTCAD::solve
