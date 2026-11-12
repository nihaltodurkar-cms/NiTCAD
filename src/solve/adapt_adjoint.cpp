// Goal adjoints and the adjoint error estimate (Unit 17, slice 4; DECISIONS.md Unit 17 AQ4, AD3, AR3).
#include "NiTCAD/solve/adapt_adjoint.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/solve/adapt_transfer.hpp"
#include "fields.hpp"
#include "weighted_current.hpp"

namespace NiTCAD::solve {
namespace {

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

// The scaled system of `device` at a bias, with a given concentration scale when given.
struct System {
    assemble::Scaling scaling;
    std::optional<assemble::DriftDiffusion> dd;
};

std::expected<System, base::Error> make_system(const device::Device& device, const BiasOptions& options,
                                               std::span<const double> bias_V, std::optional<double> Ns) {
    System s;
    auto scaling = assemble::make_scaling(device, Ns ? Ns : options.Ns_override);
    if (!scaling) return std::unexpected(scaling.error());
    s.scaling = *scaling;
    auto dd = assemble::DriftDiffusion::create(device, s.scaling, options.models);
    if (!dd) return std::unexpected(dd.error());
    s.dd.emplace(std::move(*dd));
    s.dd->set_minority_floor(options.minority_floor);
    if (auto ok = s.dd->set_bias(bias_V); !ok) return std::unexpected(ok.error());
    return s;
}

// Per node: its ohmic contact, or -1.
std::vector<std::int32_t> ohmic_of(const device::Device& device) {
    std::vector<std::int32_t> o(device.mesh().node_count(), -1);
    const auto contacts = device.contacts();
    for (std::size_t c = 0; c < contacts.size(); ++c) {
        if (contacts[c].kind != device::ContactKind::ohmic) continue;
        for (const mesh::NodeId k : contacts[c].nodes) o[static_cast<std::size_t>(k)] = static_cast<std::int32_t>(c);
    }
    return o;
}

double weighted_goal(const assemble::DriftDiffusion& dd, const device::Device& device, std::span<const double> x,
                     std::span<const double> w) {
    const auto edges = device.mesh().edges();
    const auto j = dd.edge_currents(x);
    double g = 0.0;
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const double dw = w[static_cast<std::size_t>(edges[e].first)] - w[static_cast<std::size_t>(edges[e].second)];
        if (dw != 0.0) g += (j[e].first + j[e].second) * dw;
    }
    return g;
}

std::optional<base::Error> check_contacts(const device::Device& device, std::span<const std::size_t> contacts) {
    for (const std::size_t c : contacts) {
        if (c >= device.contacts().size() || device.contacts()[c].kind != device::ContactKind::ohmic) {
            return invalid("an adjoint's goal contact is out of range or not ohmic");
        }
    }
    return std::nullopt;
}

// The adjoints at state x of `dd` (its Jacobian `jac` factorized in `solver`).
std::expected<PointAdjoint, base::Error> adjoints_at(const device::Device& device, const assemble::DriftDiffusion& dd,
                                                     const System& sys, std::vector<double> x,
                                                     const linalg::SparseMatrix& jac, linalg::LinearSolver& solver,
                                                     const BiasOptions& options,
                                                     std::span<const std::size_t> contacts) {
    const std::size_t nodes = device.mesh().node_count();
    const std::size_t n = dd.unknowns();
    PointAdjoint out;
    out.state = std::move(x);
    out.current_scale = sys.scaling.J0 * std::pow(sys.scaling.L_D, device.mesh().dimension() - 1);
    std::optional<std::vector<std::vector<double>>> weights;
    if (options.weighted_current) {
        auto wc = detail::WeightedCurrents::create(device, options.linear);
        if (!wc) return std::unexpected(wc.error());
        auto w = wc->weights(dd, out.state);
        if (!w) return std::unexpected(w.error());
        weights = std::move(*w);
    }
    const auto rows = dd.small_signal_currents(out.state, 0.0);
    const detail::Slots slot = detail::slots(dd.dofs());
    const std::vector<std::int32_t> ohmic = ohmic_of(device);
    const auto offsets = jac.row_offsets();
    const auto columns = jac.col_indices();
    const auto values = jac.values();
    const std::vector<double> terminal = dd.terminal_currents(out.state);
    for (const std::size_t c : contacts) {
        std::vector<double> h(n, 0.0);
        for (std::size_t k = 0; k < rows[c].columns.size(); ++k) h[rows[c].columns[k]] += rows[c].values[k].real();
        std::vector<double> y(n, 0.0);
        auto s = solver.solve_transposed(h, y);
        if (!s) return std::unexpected(s.error());
        out.fell_back = out.fell_back || s->fell_back;
        std::vector<double> z = y;
        std::vector<double> g = h;
        std::vector<double> w;
        double goal = terminal[c];
        if (weights) {
            w = (*weights)[c];
            for (std::size_t i = 0; i < nodes; ++i) {
                if (ohmic[i] >= 0 || device.is_insulator(static_cast<mesh::NodeId>(i)) || w[i] == 0.0) continue;
                for (const std::size_t r : {slot.block * i + slot.electrons, slot.block * i + slot.holes}) {
                    z[r] += w[i];
                    for (auto k = static_cast<std::size_t>(offsets[r]); k < static_cast<std::size_t>(offsets[r + 1]); ++k) {
                        g[static_cast<std::size_t>(columns[k])] += w[i] * values[k];
                    }
                }
            }
            goal = weighted_goal(dd, device, out.state, w);
        }
        out.weights.push_back(std::move(w));
        out.goal.push_back(goal);
        out.gradient.push_back(std::move(g));
        out.adjoint.push_back(std::move(z));
        out.terminal.push_back(std::move(y));
    }
    return out;
}

}  // namespace

std::optional<base::Error> check_adjoint(const BiasOptions& options) {
    if (options.equations != Equations::drift_diffusion) {
        return invalid("adjoint goal weighting needs the drift-diffusion equations");
    }
    if (options.models.electrothermal) return invalid("adjoint goal weighting does not support the electrothermal model");
    if (options.models.energy_balance) return invalid("adjoint goal weighting does not support the energy-balance model");
    if (options.models.btbt_nonlocal != assemble::NonlocalTunnelling::off) {
        return invalid("adjoint goal weighting does not support nonlocal tunnelling");
    }
    if (options.linear.backend == linalg::SolverBackend::gmres_ilut &&
        options.linear.fallback == linalg::LinearFallback::none) {
        return invalid("adjoint goal weighting needs a transposed solve: gmres_ilut without the sparse_lu "
                       "fallback cannot give one");
    }
    return std::nullopt;
}

std::expected<PointAdjoint, base::Error> point_adjoint(const device::Device& device, const BiasOptions& options,
                                                       const results::BiasPoint& point,
                                                       std::span<const std::size_t> contacts) {
    if (auto e = check_adjoint(options)) return std::unexpected(std::move(*e));
    if (auto e = check_contacts(device, contacts)) return std::unexpected(std::move(*e));
    const std::size_t nodes = device.mesh().node_count();
    if (point.fields.potential_V.size() != nodes || point.fields.n_cm3.size() != nodes ||
        point.fields.p_cm3.size() != nodes) {
        return std::unexpected(invalid("an adjoint's point needs its fields, one per node"));
    }
    auto sys = make_system(device, options, point.bias_V, std::nullopt);
    if (!sys) return std::unexpected(sys.error());
    const assemble::DriftDiffusion& dd = *sys->dd;
    std::vector<double> x(dd.unknowns(), 0.0);
    detail::scaled_state(dd, point.fields, sys->scaling, x);
    dd.stamp_contacts(x);
    linalg::SparseMatrix jac = dd.make_jacobian();
    std::vector<double> f(dd.unknowns());
    dd.evaluate(x, f, jac);
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(solver.error());
    if (auto fz = solver->factorize(jac); !fz) return std::unexpected(fz.error());
    return adjoints_at(device, dd, *sys, std::move(x), jac, *solver, options, contacts);
}

std::expected<double, base::Error> goal_at(const device::Device& device, const BiasOptions& options,
                                           std::span<const double> bias_V, std::span<const double> x,
                                           std::size_t contact, std::span<const double> weights) {
    auto sys = make_system(device, options, bias_V, std::nullopt);
    if (!sys) return std::unexpected(sys.error());
    if (weights.empty()) return sys->dd->terminal_currents(x)[contact];
    return weighted_goal(*sys->dd, device, x, weights);
}

std::expected<GoalEstimate, base::Error> estimate_goals(const device::Device& coarse,
                                                        const mesh::PointLocator& coarse_locator,
                                                        const results::BiasPoint& point, const device::Device& fine,
                                                        const BiasOptions& options,
                                                        std::span<const std::size_t> contacts) {
    if (auto e = check_adjoint(options)) return std::unexpected(std::move(*e));
    if (auto e = check_contacts(coarse, contacts)) return std::unexpected(std::move(*e));
    if (auto e = check_contacts(fine, contacts)) return std::unexpected(std::move(*e));
    for (const std::size_t c : contacts) {
        const auto& g = options.weighted_current ? point.weighted_current : point.terminal_current;
        if (c >= g.size()) return std::unexpected(invalid("an estimate's point carries no current for its goal"));
    }
    // M's concentration scale on M', so both residuals and goals are in the same units.
    auto coarse_scaling = assemble::make_scaling(coarse, options.Ns_override);
    if (!coarse_scaling) return std::unexpected(coarse_scaling.error());
    auto fine_sys = make_system(fine, options, point.bias_V, coarse_scaling->Ns);
    if (!fine_sys) return std::unexpected(fine_sys.error());
    const assemble::DriftDiffusion& fd = *fine_sys->dd;
    auto moved = transfer_state(coarse, coarse_locator, point.fields, point.bands, fine, nullptr,
                                TransferMode::exponential_fit);
    if (!moved) return std::unexpected(moved.error());
    const std::size_t n = fd.unknowns();
    std::vector<double> u(n, 0.0);
    detail::scaled_state(fd, *moved, fine_sys->scaling, u);
    fd.stamp_contacts(u);
    linalg::SparseMatrix jac = fd.make_jacobian();
    std::vector<double> r(n, 0.0);
    fd.evaluate(u, r, jac);
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(solver.error());
    if (auto fz = solver->factorize(jac); !fz) return std::unexpected(fz.error());
    auto adj = adjoints_at(fine, fd, *fine_sys, u, jac, *solver, options, contacts);
    if (!adj) return std::unexpected(adj.error());

    const std::size_t nf = fine.mesh().node_count();
    const std::size_t nc = coarse.mesh().node_count();
    const auto fine_points = fine.mesh().points();
    const auto cr = coarse.node_region();
    const auto fr = fine.node_region();
    const detail::Slots fs = detail::slots(fd.dofs());
    // Each fine node's owner: the coarse node of its largest interpolation weight in its own
    // region (any region when none is; lowest index on ties).
    std::vector<std::size_t> owner(nf, 0);
    for (std::size_t q = 0; q < nf; ++q) {
        const mesh::Interpolant it = coarse_locator.locate(fine_points[q]);
        bool own = false;
        for (int j = 0; j < it.count; ++j) own = own || cr[static_cast<std::size_t>(it.nodes[static_cast<std::size_t>(j)])] == fr[q];
        double best = -1.0;
        for (int j = 0; j < it.count; ++j) {
            const auto k = static_cast<std::size_t>(it.nodes[static_cast<std::size_t>(j)]);
            const double wt = it.weights[static_cast<std::size_t>(j)];
            if (own && cr[k] != fr[q]) continue;
            if (wt > best || (wt == best && k < owner[q])) {
                best = wt;
                owner[q] = k;
            }
        }
    }

    GoalEstimate out;
    out.fell_back = adj->fell_back;
    const double scale = adj->current_scale;
    for (std::size_t gi = 0; gi < contacts.size(); ++gi) {
        const std::size_t c = contacts[gi];
        const double coarse_goal = options.weighted_current ? point.weighted_current[c] : point.terminal_current[c];
        std::vector<double> local(nc, 0.0);
        const std::vector<double>& z = adj->adjoint[gi];
        double dual = 0.0;
        for (std::size_t q = 0; q < nf; ++q) {
            double v = 0.0;
            for (const std::size_t k : {fs.potential, fs.electrons, fs.holes}) v += z[fs.block * q + k] * r[fs.block * q + k];
            local[owner[q]] += v;
            dual += v;
        }
        const double prolongated = adj->goal[gi] * scale - coarse_goal;
        const double d = prolongated - dual * scale;
        out.difference.push_back(d);
        out.estimate.push_back(2.0 * std::abs(d));
        out.prolongated_part.push_back(prolongated);
        for (double& v : local) v = std::abs(v) * scale;
        out.node_weight.push_back(std::move(local));
    }
    return out;
}

}  // namespace NiTCAD::solve
