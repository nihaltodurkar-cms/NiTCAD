#include "density_gradient.hpp"

#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>

#include "NiTCAD/assemble/equilibrium.hpp"

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

}  // namespace NiTCAD::solve::detail
