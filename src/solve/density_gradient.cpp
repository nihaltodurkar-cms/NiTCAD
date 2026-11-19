#include "density_gradient.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>

#include "NiTCAD/assemble/equilibrium.hpp"
#include "NiTCAD/base/constants.hpp"

namespace NiTCAD::solve::detail {

std::expected<std::vector<double>, base::Error> classical_equilibrium(
    const device::Device& device, const assemble::Scaling& scaling,
    const assemble::PhysicsModels& models, std::span<const double> bias, const NewtonOptions& newton,
    linalg::LinearSolver& solver, results::ConvergenceRecord& record, const std::stop_token& stop,
    const IterationObserver& observe) {
    assemble::PhysicsModels plain = models;
    plain.quantum_confinement = assemble::QuantumConfinement::off;
    auto system = assemble::System::create(device, scaling, plain, assemble::Formulation::equilibrium);
    if (!system) return std::unexpected(std::move(system.error()));
    if (auto ok = system->set_bias(bias); !ok) return std::unexpected(std::move(ok.error()));
    std::vector<double> psi = assemble::charge_neutral_potential(*system);
    if (auto ok = newton_solve(*system, psi, newton, solver, record, stop, observe); !ok) {
        return std::unexpected(std::move(ok.error()));
    }
    return psi;
}

std::expected<QuantumEquilibrium, base::Error> solve_quantum_equilibrium(
    const device::Device& device, const assemble::Scaling& scaling,
    const assemble::PhysicsModels& models, std::span<const double> bias, std::span<const double> psi,
    const NewtonOptions& newton, linalg::LinearSolver& solver, results::ConvergenceRecord& record,
    const std::stop_token& stop, const IterationObserver& observe) {
    const std::size_t n = psi.size();
    // The system at beta k: beta is linear in gamma (physics::dg_coefficient_cm2), and the barrier
    // terms are beta B, so gamma k scales the whole quantum term.
    const auto solve_at = [&](double k, std::vector<double>& x) -> std::expected<void, base::Error> {
        assemble::PhysicsModels m = models;
        m.dg_gamma_n *= k;
        m.dg_gamma_p *= k;
        auto system = assemble::System::create(device, scaling, m, assemble::Formulation::equilibrium);
        if (!system) return std::unexpected(std::move(system.error()));
        if (auto ok = system->set_bias(bias); !ok) return ok;
        return newton_solve(*system, x, newton, solver, record, stop, observe);
    };
    std::vector<double> start(3 * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) start[3 * i] = psi[i];
    QuantumEquilibrium out;
    std::vector<double> x = start;
    auto first = solve_at(1.0, x);
    if (!first && first.error().code != base::ErrorCode::non_convergence) {
        return std::unexpected(std::move(first.error()));
    }
    if (!first) {  // QP9's ladder from the classical state
        std::vector<double> last = start;
        double done = 0.0;
        for (const double target : {0.125, 0.25, 0.5, 1.0}) {
            double k = target;
            int halvings = 0;
            while (done < target) {
                x = last;
                auto ok = solve_at(k, x);
                if (ok) {
                    last = x;
                    done = k;
                    if (k < 1.0) out.rungs.push_back(k);
                    k = target;
                    halvings = 0;
                    continue;
                }
                if (ok.error().code != base::ErrorCode::non_convergence || ++halvings > 4) {
                    std::string message = "density gradient: the coefficient continuation failed at k = ";
                    char buffer[64];
                    std::snprintf(buffer, sizeof buffer, "%.6g", k);
                    message += buffer;
                    message += " (converged:";
                    for (const double r : out.rungs) {
                        std::snprintf(buffer, sizeof buffer, " %.6g", r);
                        message += buffer;
                    }
                    message += out.rungs.empty() ? " none)" : ")";
                    message += "; ";
                    message += ok.error().message;
                    return std::unexpected(base::Error{ok.error().code, std::move(message), ok.error().context});
                }
                k = done + 0.5 * (k - done);
            }
        }
        x = last;
    }
    out.psi.resize(n);
    out.lambda_n.resize(n);
    out.lambda_p.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        out.psi[i] = x[3 * i];
        out.lambda_n[i] = x[3 * i + 1];
        out.lambda_p[i] = x[3 * i + 2];
    }
    return out;
}

std::string continuation_warning(std::span<const double> rungs) {
    if (rungs.empty()) return {};
    std::string s = "density gradient: the equilibrium start used the coefficient continuation (k =";
    char buffer[32];
    for (std::size_t i = 0; i < rungs.size(); ++i) {
        std::snprintf(buffer, sizeof buffer, "%s %.6g", i == 0 ? "" : ",", rungs[i]);
        s += buffer;
    }
    return s + ")";
}

CoarseScan::CoarseScan(const device::Device& device)
    : bulk_(device.mesh().node_count(), 0),
      kT_eV_(base::thermal_voltage(device.temperature_K())),
      flagged_(device.mesh().node_count(), 0) {
    const std::size_t n = bulk_.size();
    for (std::size_t i = 0; i < n; ++i) bulk_[i] = device.is_insulator(static_cast<mesh::NodeId>(i)) ? 0 : 1;
    for (const mesh::Edge& e : device.mesh().edges()) {
        const auto a = static_cast<std::size_t>(e.first), b = static_cast<std::size_t>(e.second);
        if (device.is_insulator(e.first) != device.is_insulator(e.second)) bulk_[a] = bulk_[b] = 0;
    }
    for (const device::Contact& c : device.contacts()) {
        if (c.kind != device::ContactKind::gate) continue;
        for (const mesh::NodeId i : c.nodes) bulk_[static_cast<std::size_t>(i)] = 0;
    }
}

void CoarseScan::add(const results::NodeFields& f) {
    if (f.quantum_potential_n_eV.size() != bulk_.size()) return;
    ++results_;
    bool any = false;
    for (std::size_t i = 0; i < bulk_.size(); ++i) {
        if (bulk_[i] == 0) continue;
        const double l = std::max(std::abs(f.quantum_potential_n_eV[i]), std::abs(f.quantum_potential_p_eV[i])) / kT_eV_;
        if (!(l > 10.0)) continue;
        any = true;
        flagged_[i] = 1;
        if (l > largest_) {
            largest_ = l;
            where_ = i;
        }
    }
    if (any) ++flagged_results_;
}

std::string CoarseScan::warning() const {
    if (flagged_results_ == 0) return {};
    std::size_t nodes = 0;
    for (const char c : flagged_) nodes += c != 0 ? 1 : 0;
    char buffer[320];
    std::snprintf(buffer, sizeof buffer,
                  "density gradient: |Lambda| above 10 kT at %zu node(s) away from the interfaces in %zu of %zu "
                  "result(s), largest %.3g kT (node %zu); the mesh may not resolve the densities there (QP15)",
                  nodes, flagged_results_, results_, largest_, where_);
    return buffer;
}

}  // namespace NiTCAD::solve::detail
