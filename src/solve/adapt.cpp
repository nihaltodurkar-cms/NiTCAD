// The adaptive refinement loop (Unit 17, slice 3; DECISIONS.md Unit 17 AD1, AD4-AD6, AD9, AD11, AD13).
#include "NiTCAD/solve/adapt.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/mesh/locate.hpp"
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

// The 2D refinement of a triangulation's nodal sizes (AD4).
std::expected<std::shared_ptr<const mesh::SizeField>, base::Error> refine_field(
    const mesh::SimplexMesh& sm, const std::shared_ptr<const mesh::SizeField>& previous, const std::vector<char>& marked,
    double min_size, double grading, bool uniform) {
    std::vector<double> h = mesh::local_sizes(sm);
    if (previous) {
        for (std::size_t i = 0; i < h.size(); ++i) h[i] = std::min(h[i], previous->at(sm.points[i]));
    }
    if (uniform) {
        for (double& x : h) x *= 0.5;
    } else {
        std::vector<char> seeds(h.size(), 0);
        for (std::size_t i = 0; i < h.size(); ++i) {
            if (marked[i] == 0 || !(h[i] > min_size)) continue;
            h[i] = std::max(0.5 * h[i], min_size);
            seeds[i] = 1;
        }
        auto g = mesh::grade_sizes(sm, h, grading, seeds);
        if (!g) return std::unexpected(g.error());
    }
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
    const double floor2d = o.min_size > 0.0             ? o.min_size
                           : cur.interface_resolution > 0 ? 0.25 * cur.interface_resolution
                                                          : cur.max_size / 64.0;
    MeshRefinement next = previous;
    switch (cur.kind) {
        case AdaptMeshKind::generated2d: {
            if (!cur.generated || cur.generated->simplices.points.size() != nn) {
                return std::unexpected(invalid("a generated 2D mesh needs its simplices, one point per node"));
            }
            auto f = refine_field(cur.generated->simplices, previous.size_field, mk, floor2d, o.grading, uniform);
            if (!f) return std::unexpected(f.error());
            next.size_field = std::move(*f);
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
                next.axes[static_cast<std::size_t>(a)] = std::move(*r);
            }
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
            auto f = refine_field(*cur.section, previous.size_field, section_marks, floor2d, o.grading, uniform);
            if (!f) return std::unexpected(f.error());
            next.size_field = std::move(*f);
            auto r = refine_lines(lines, im, cur.sweep_protected, o, uniform);
            if (!r) return std::unexpected(r.error());
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

namespace {

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
    double seconds_mesh = seconds_since(t0);
    std::optional<AdaptableMesh> previous_mesh;
    results::Sweep previous_sweep;
    std::vector<double> previous_goal;

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
    // M' for a pass: the uniform refinement of its mesh, built (the estimate's and AD11's mesh).
    const auto build_uniform = [&](const AdaptableMesh& m,
                                   const MeshRefinement& r) -> std::expected<AdaptableMesh, base::Error> {
        auto uniform = next_refinement(m, r, {}, nullptr, o.indicators, o, true);
        if (!uniform) return std::unexpected(uniform.error());
        return build(*uniform);
    };
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
        auto solved = solve_on(*built, previous_mesh ? &*previous_mesh : nullptr, previous_mesh ? &previous_sweep : nullptr,
                               points, bias, control);
        if (!solved) return std::unexpected(solved.error());
        AdaptPass rep;
        rep.nodes = built->device.mesh().node_count();
        rep.edges = built->device.mesh().edges().size();
        rep.transfer = solved->transfer;
        rep.seconds_mesh = seconds_mesh;
        rep.seconds_transfer = solved->seconds_transfer;
        rep.seconds_solve = solved->seconds_solve;
        for (const results::BiasPoint& p : solved->sweep.points) rep.newton_iterations += static_cast<int>(p.convergence.iterations.size());
        const Goals g = read_goals(solved->sweep, o, at);
        rep.goal = g.value;
        rep.resolution = g.resolution;
        // The adjoint estimate and weights (AD3, AD6; AR3): on M', per adaptation point.
        std::optional<AdaptableMesh> fine_mesh;
        std::vector<double> omega;  // AD3's node weight: per node, the largest over goals and points
        if (adjoint) {
            const auto te = clock_type::now();
            auto f = build_uniform(*built, refinement);
            if (!f) return std::unexpected(f.error());
            fine_mesh = std::move(*f);
            auto loc = locator_of(*built);
            if (!loc) return std::unexpected(loc.error());
            const std::size_t nn = rep.nodes;
            omega.assign(nn, 0.0);
            rep.estimate.assign(g.value.size(), 0.0);
            for (std::size_t a = 0; a < at.size(); ++a) {
                auto e = estimate_goals(built->device, *loc, solved->sweep.points[at[a]], fine_mesh->device, bias,
                                        goal_contacts);
                if (!e) return std::unexpected(e.error());
                rep.estimate_fell_back = rep.estimate_fell_back || e->fell_back;
                for (std::size_t gi = 0; gi < goal_contacts.size(); ++gi) {
                    rep.estimate[gi * at.size() + a] = e->estimate[gi];
                    // Normalized so the semiconductor nodes' weights sum to 1 (AD3).
                    double total = 0.0;
                    for (std::size_t i = 0; i < nn; ++i) {
                        if (!built->device.is_insulator(static_cast<mesh::NodeId>(i))) total += e->node_weight[gi][i];
                    }
                    if (!(total > 0.0)) continue;
                    for (std::size_t i = 0; i < nn; ++i) {
                        if (built->device.is_insulator(static_cast<mesh::NodeId>(i))) continue;
                        omega[i] = std::max(omega[i], e->node_weight[gi][i] / total);
                    }
                }
            }
            rep.estimate_nodes = fine_mesh->device.mesh().node_count();
            rep.seconds_estimate = seconds_since(te);
        }
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
            }
        } else if (converged) {
            stop = AdaptStop::converged;
        }
        if (!stop && pass + 1 >= o.max_passes) stop = AdaptStop::max_passes;
        std::optional<MeshRefinement> next;
        if (!stop) {
            // Indicators at every adaptation point, united (AQ3), marked by AD3: weighted by the
            // adjoint's omega, or feature-only (w = 1).
            std::vector<IndicatorField> fields;
            std::optional<EdgeIndicators> edges;
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
            auto united = unite(fields);
            if (!united) return std::unexpected(united.error());
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
            if (by_difference && mk->count > 0) {
                auto r = next_refinement(*built, refinement, mk->marked, nullptr, o.indicators, o, false);
                if (!r) return std::unexpected(r.error());
                next = std::move(*r);
            } else if (mk->count == 0 && rep.rejected) {
                // A rejected stop with nothing marked: the uniform refinement (M') is the next mesh.
                auto r = next_refinement(*built, refinement, {}, nullptr, o.indicators, o, true);
                if (!r) return std::unexpected(r.error());
                next = std::move(*r);
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
            seconds_mesh = seconds_since(t0);
            if (o.max_nodes > 0 && nb->device.mesh().node_count() > o.max_nodes) {
                stop = AdaptStop::max_nodes;
            } else {
                previous_mesh = std::move(*built);
                previous_sweep = std::move(solved->sweep);
                previous_goal = g.value;
                built = std::move(*nb);
                refinement = std::move(*next);
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
    result.mesh = std::move(*built);
    return result;
}

}  // namespace NiTCAD::solve
