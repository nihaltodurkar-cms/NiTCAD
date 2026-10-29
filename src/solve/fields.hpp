// Private to the solve layer: the physical fields of an equilibrium-Poisson state, shared by
// solve_equilibrium and the quasi-static bias sweep, and of a drift-diffusion state.
#pragma once

#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::solve::detail {

// Potential in V and the Boltzmann carriers in cm^-3 of the scaled potential psi.
inline results::NodeFields equilibrium_fields(const assemble::EquilibriumPoisson& system,
                                              std::span<const double> psi,
                                              const assemble::Scaling& scaling) {
    const std::size_t n = system.unknowns();
    results::NodeFields f{std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    system.carriers(psi, f.n_cm3, f.p_cm3);
    for (std::size_t i = 0; i < n; ++i) {
        f.potential_V[i] = psi[i] * scaling.V_T;
        f.n_cm3[i] *= scaling.Ns;
        f.p_cm3[i] *= scaling.Ns;
    }
    return f;
}

// The band diagram in eV from the assemblers' scaled one (units of V_T).
inline results::BandDiagram band_diagram(const assemble::BandEdges& b, double V_T) {
    results::BandDiagram d{b.conduction, b.valence, b.electron_fermi, b.hole_fermi};
    for (auto* v : {&d.conduction_eV, &d.valence_eV, &d.electron_fermi_eV, &d.hole_fermi_eV}) {
        for (double& e : *v) e *= V_T;
    }
    return d;
}

// Scale from scaled gate charges to C / cm^(3-D): q Ns L_D^D.
inline double charge_scale(const assemble::Scaling& scaling, int dimension) {
    return base::q_C * scaling.Ns * std::pow(scaling.L_D, dimension);
}

// The positions of the drift-diffusion fields in each node's block of unknowns (the system's
// assemble::DofMap, Unit 40): field f of node i is unknown block i + f. temperature only with the
// electrothermal model's field.
struct Slots {
    std::size_t block, potential, electrons, holes, temperature;
};

inline Slots slots(const assemble::DofMap& dofs) {
    using assemble::FieldId;
    return {dofs.block(), *dofs.slot(FieldId::potential), *dofs.slot(FieldId::electrons),
            *dofs.slot(FieldId::holes), dofs.slot(FieldId::temperature).value_or(0)};
}

// The scaled state of physical node fields into x (the inverse of drift_diffusion_point):
// psi / V_T, n / Ns, p / Ns and, with the temperature field, the rise tau = (T - T0) / T0 (0 when
// the fields carry no temperature). Precondition: x has the system's unknowns (or more: the rest
// is not written).
inline void scaled_state(const assemble::DriftDiffusion& system, const results::NodeFields& f,
                         const assemble::Scaling& scaling, std::span<double> x) {
    const Slots k = slots(system.dofs());
    const bool thermal = system.dofs().has(assemble::FieldId::temperature);
    const double T0 = scaling.temperature_K;
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[k.block * i + k.potential] = f.potential_V[i] / scaling.V_T;
        x[k.block * i + k.electrons] = f.n_cm3[i] / scaling.Ns;
        x[k.block * i + k.holes] = f.p_cm3[i] / scaling.Ns;
        if (thermal) {
            x[k.block * i + k.temperature] =
                f.temperature_K.empty() ? 0.0 : (f.temperature_K[i] - T0) / T0;
        }
    }
}

// The generation-recombination rates of a drift-diffusion state (Unit 45): the system's mechanisms
// (DriftDiffusion::rates, already in cm^-3 s^-1) into RateFields, and their sum.
inline results::RateFields rate_fields(const assemble::DriftDiffusion& system,
                                       std::span<const double> x) {
    using M = assemble::DriftDiffusion::RateMechanism;
    results::RateFields r;
    const std::size_t nodes = system.node_count();
    r.net_n.assign(nodes, 0.0);
    r.net_p.assign(nodes, 0.0);
    for (assemble::DriftDiffusion::MechanismRate& m : system.rates(x)) {
        for (std::size_t i = 0; i < nodes; ++i) {
            r.net_n[i] += m.electron[i];
            r.net_p[i] += m.hole[i];
        }
        switch (m.mechanism) {
            case M::srh: r.srh = std::move(m.electron); break;
            case M::auger: r.auger = std::move(m.electron); break;
            case M::radiative: r.radiative = std::move(m.electron); break;
            case M::impact: r.impact = std::move(m.electron); break;
            case M::interfaces:
                r.interface_n = std::move(m.electron);
                r.interface_p = std::move(m.hole);
                break;
            case M::band_to_band:
                r.band_to_band_n = std::move(m.electron);
                r.band_to_band_p = std::move(m.hole);
                break;
        }
    }
    return r;
}

// The fields and terminal quantities of a drift-diffusion state x into `point` (its bias and
// convergence are the caller's), in the public units; shared by sweep_bias and trace_bias. With
// `rates`, its generation-recombination rates too (Unit 45).
inline void drift_diffusion_point(const assemble::DriftDiffusion& system,
                                  std::span<const double> x, const assemble::Scaling& scaling,
                                  int dimension, results::BiasPoint& point, bool rates = true) {
    if (rates) point.rates = rate_fields(system, x);
    const std::size_t nodes = system.node_count();
    const double current_scale = scaling.J0 * std::pow(scaling.L_D, dimension - 1);
    const double charges = charge_scale(scaling, dimension);
    point.fields = {std::vector<double>(nodes), std::vector<double>(nodes),
                    std::vector<double>(nodes)};
    const Slots k = slots(system.dofs());
    for (std::size_t i = 0; i < nodes; ++i) {
        point.fields.potential_V[i] = x[k.block * i + k.potential] * scaling.V_T;
        point.fields.n_cm3[i] = x[k.block * i + k.electrons] * scaling.Ns;
        point.fields.p_cm3[i] = x[k.block * i + k.holes] * scaling.Ns;
    }
    if (system.electrothermal()) {  // Unit 23: the temperature and the heat to each sink
        for (std::size_t i = 0; i < nodes; ++i) {
            const double T0 = scaling.temperature_K;  // tau = (T - T0) / T0
            point.fields.temperature_K.push_back(T0 + T0 * x[k.block * i + k.temperature]);
        }
        point.thermal_contact_heat = system.thermal_contact_heat(x);
        for (double& P : point.thermal_contact_heat) P *= scaling.V_T * current_scale;
    }
    point.bands = band_diagram(system.band_edges(x), scaling.V_T);
    point.terminal_current = system.terminal_currents(x);
    for (double& I : point.terminal_current) I *= current_scale;
    point.terminal_current_resolution = system.terminal_current_resolution(x);
    for (double& I : point.terminal_current_resolution) I *= current_scale;
    point.gate_charge = system.gate_charges(x);
    for (double& Q : point.gate_charge) Q *= charges;
    point.interface_trap_charge = system.interface_trap_charges(x);
    for (double& Q : point.interface_trap_charge) Q *= charges;
    for (const auto& [jn, jp] : system.edge_currents(x)) {
        point.edge_current_n.push_back(jn * current_scale);
        point.edge_current_p.push_back(jp * current_scale);
    }
    for (const assemble::TunnelPathState& s : system.path_states(x)) {
        point.tunnel_paths.push_back(
            {s.start, s.length_cm, s.field_V_per_cm, s.rate_cm3_s, s.reached});
    }
    if (system.gate_tunnelling()) {  // Unit 51
        auto t = system.gate_tunnelling_currents(x);
        for (std::vector<double>* v : {&t.ecb, &t.hvb, &t.evb}) {
            for (double& I : *v) I *= current_scale;
        }
        point.tunnelling_ecb = std::move(t.ecb);
        point.tunnelling_hvb = std::move(t.hvb);
        point.tunnelling_evb = std::move(t.evb);
    }
}

}  // namespace NiTCAD::solve::detail
