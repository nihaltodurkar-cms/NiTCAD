#include "weighted_current.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace NiTCAD::solve::detail {

std::expected<WeightedCurrents, base::Error> WeightedCurrents::create(
    const device::Device& device, const linalg::SolverConfig& linear) {
    // Unit 48 (DECISIONS.md SD7): the weights are built on ohmic contacts only.
    if (std::ranges::any_of(device.contacts(), [](const device::Contact& c) {
            return c.kind == device::ContactKind::schottky;
        })) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input,
            "the weighted terminal current with a Schottky contact is not supported", {}});
    }
    // The Laplace problem is an M-matrix whose conductances span many decades: its pivot ratio says
    // nothing about singularity (every part of it is tied to an ohmic contact by the device's
    // topology check), so the pivot test is off for it.
    linalg::SolverConfig laplace_config = linear;
    laplace_config.min_pivot_ratio = 0.0;
    auto laplace = linalg::LinearSolver::create(laplace_config);
    if (!laplace) return std::unexpected(std::move(laplace.error()));
    auto newton = linalg::LinearSolver::create(linear);
    if (!newton) return std::unexpected(std::move(newton.error()));
    WeightedCurrents w(std::move(*laplace), std::move(*newton));
    w.device_ = &device;
    const std::size_t n = device.mesh().node_count();
    w.ohmic_.assign(n, -1);
    w.row_.assign(n, -1);
    w.carrier_.assign(n, 0);
    const auto contacts = device.contacts();
    for (std::size_t c = 0; c < contacts.size(); ++c) {
        if (contacts[c].kind != device::ContactKind::ohmic) continue;
        for (const mesh::NodeId k : contacts[c].nodes) {
            w.ohmic_[static_cast<std::size_t>(k)] = static_cast<std::int32_t>(c);
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        w.carrier_[i] = device.is_insulator(static_cast<mesh::NodeId>(i)) ? 0 : 1;
        if (w.carrier_[i] != 0 && w.ohmic_[i] < 0) {
            w.row_[i] = static_cast<std::int32_t>(w.rows_++);
        }
    }
    return w;
}

std::expected<std::vector<std::vector<double>>, base::Error> WeightedCurrents::weights(
    const assemble::DriftDiffusion& system, std::span<const double> x) {
    const auto contacts = device_->contacts();
    const auto edges = device_->mesh().edges();
    const std::vector<double> t = system.edge_flux_terms(x);
    double largest = 0.0;
    for (const double v : t) largest = std::max(largest, v);
    const double floor = 1e-30 * largest;
    // Each edge between two semiconductor nodes couples them with g = t / 2 (+ the floor); an edge
    // to an ohmic contact node moves its share to the right-hand side of that contact's problem.
    std::vector<linalg::Triplet> entries;
    std::vector<std::vector<double>> rhs(contacts.size());
    for (std::size_t c = 0; c < contacts.size(); ++c) {
        if (contacts[c].kind == device::ContactKind::ohmic) rhs[c].assign(rows_, 0.0);
    }
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto a = static_cast<std::size_t>(edges[e].first);
        const auto b = static_cast<std::size_t>(edges[e].second);
        if (carrier_[a] == 0 || carrier_[b] == 0) continue;
        const double g = 0.5 * t[e] + floor;
        const std::int32_t ra = row_[a], rb = row_[b];
        if (ra >= 0) entries.push_back({ra, ra, g});
        if (rb >= 0) entries.push_back({rb, rb, g});
        if (ra >= 0 && rb >= 0) {
            entries.push_back({ra, rb, -g});
            entries.push_back({rb, ra, -g});
        } else if (ra >= 0 && ohmic_[b] >= 0) {
            rhs[static_cast<std::size_t>(ohmic_[b])][static_cast<std::size_t>(ra)] += g;
        } else if (rb >= 0 && ohmic_[a] >= 0) {
            rhs[static_cast<std::size_t>(ohmic_[a])][static_cast<std::size_t>(rb)] += g;
        }
    }
    std::vector<std::vector<double>> w(contacts.size());
    const std::size_t n = device_->mesh().node_count();
    std::vector<double> solution(rows_);
    if (rows_ > 0) {
        const auto size = static_cast<linalg::Index>(rows_);
        auto matrix = linalg::SparseMatrix::from_triplets(size, size, entries);
        if (!matrix) return std::unexpected(std::move(matrix.error()));
        if (auto f = laplace_.factorize(*matrix); !f) return std::unexpected(std::move(f.error()));
    }
    for (std::size_t c = 0; c < contacts.size(); ++c) {
        if (contacts[c].kind != device::ContactKind::ohmic) continue;
        if (rows_ > 0) {
            if (auto s = laplace_.solve(rhs[c], solution); !s) {
                return std::unexpected(std::move(s.error()));
            }
        }
        w[c].assign(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (row_[i] >= 0) w[c][i] = solution[static_cast<std::size_t>(row_[i])];
            else if (ohmic_[i] == static_cast<std::int32_t>(c)) w[c][i] = 1.0;
        }
    }
    return w;
}

std::vector<double> WeightedCurrents::currents(const assemble::DriftDiffusion& system,
                                               std::span<const double> x,
                                               const std::vector<std::vector<double>>& w) const {
    const auto edges = device_->mesh().edges();
    const auto j = system.edge_currents(x);
    std::vector<double> current(w.size(), 0.0);
    for (std::size_t c = 0; c < w.size(); ++c) {
        if (w[c].empty()) continue;
        for (std::size_t e = 0; e < edges.size(); ++e) {
            const double dw = w[c][static_cast<std::size_t>(edges[e].first)] -
                              w[c][static_cast<std::size_t>(edges[e].second)];
            if (dw != 0.0) current[c] += (j[e].first + j[e].second) * dw;
        }
    }
    return current;
}

std::expected<void, base::Error> WeightedCurrents::fill(const assemble::DriftDiffusion& system,
                                                        std::span<const double> x,
                                                        double max_update,
                                                        const assemble::Scaling& scaling,
                                                        int dimension, results::BiasPoint& point) {
    auto w = weights(system, x);
    if (!w) return std::unexpected(std::move(w.error()));
    const std::vector<double> I0 = currents(system, x, *w);

    // The rounding of the weighted edge currents.
    const auto edges = device_->mesh().edges();
    const std::vector<double> t = system.edge_flux_terms(x);
    const double eps8 = 8.0 * std::numeric_limits<double>::epsilon();
    std::vector<double> rounding(w->size(), 0.0);
    for (std::size_t c = 0; c < w->size(); ++c) {
        if ((*w)[c].empty()) continue;
        for (std::size_t e = 0; e < edges.size(); ++e) {
            const double dw = (*w)[c][static_cast<std::size_t>(edges[e].first)] -
                              (*w)[c][static_cast<std::size_t>(edges[e].second)];
            rounding[c] += eps8 * t[e] * std::abs(dw);
        }
    }

    // The adjoint estimate: the change over one more Newton correction, w held fixed.
    const double current_scale = scaling.J0 * std::pow(scaling.L_D, dimension - 1);
    std::vector<double> estimate(w->size(), 0.0);
    bool estimated = false;
    {
        const std::size_t n = system.unknowns();
        linalg::SparseMatrix jacobian = system.make_jacobian();
        std::vector<double> f(n), rhs(n), dx(n);
        system.evaluate(x, f, jacobian);
        bool finite = true;
        for (std::size_t i = 0; i < n; ++i) {
            finite = finite && std::isfinite(f[i]);
            rhs[i] = -f[i];
        }
        if (finite && newton_.factorize(jacobian) && newton_.solve(rhs, dx)) {
            std::vector<double> next(x.begin(), x.end());
            system.apply_update(next, dx, max_update);
            const std::vector<double> I1 = currents(system, next, *w);
            for (std::size_t c = 0; c < w->size(); ++c) estimate[c] = std::abs(I1[c] - I0[c]);
            estimated = true;
        }
    }
    point.weighted_current.assign(w->size(), 0.0);
    point.weighted_current_resolution.assign(w->size(), 0.0);
    for (std::size_t c = 0; c < w->size(); ++c) {
        if ((*w)[c].empty()) continue;
        point.weighted_current[c] = I0[c] * current_scale;
        const double fallback = c < point.terminal_current_resolution.size()
                                    ? point.terminal_current_resolution[c]
                                    : 0.0;
        point.weighted_current_resolution[c] =
            estimated ? (estimate[c] + rounding[c]) * current_scale
                      : std::max(rounding[c] * current_scale, fallback);
    }
    return {};
}

}  // namespace NiTCAD::solve::detail
