// Private to the solve layer: the physical fields of an equilibrium-Poisson state, shared by
// solve_equilibrium and the quasi-static bias sweep, and of a drift-diffusion state.
#pragma once

#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/charges.hpp"
#include "NiTCAD/assemble/currents.hpp"
#include "NiTCAD/assemble/energy.hpp"
#include "NiTCAD/assemble/equilibrium.hpp"
#include "NiTCAD/assemble/derived.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/assemble/thermal.hpp"
#include "NiTCAD/assemble/tunnelling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::solve::detail {

// Potential in V and the Boltzmann carriers in cm^-3 of the scaled potential psi.
inline results::NodeFields equilibrium_fields(const assemble::System& system,
                                              std::span<const double> psi,
                                              const assemble::Scaling& scaling) {
    const std::size_t n = system.node_count();
    results::NodeFields f{std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    equilibrium_carriers(system, psi, f.n_cm3, f.p_cm3);
    for (std::size_t i = 0; i < n; ++i) {
        f.potential_V[i] = psi[i] * scaling.V_T;
        f.n_cm3[i] *= scaling.Ns;
        f.p_cm3[i] *= scaling.Ns;
    }
    return f;
}

// A state's band diagram in eV, through the derived-field table (Unit 57 GE4; thermal
// equilibrium's too, GE5).
inline results::BandDiagram band_diagram(const assemble::System& system,
                                         std::span<const double> x, double V_T) {
    using assemble::FieldId;
    const assemble::Diagnostic& bands = assemble::diagnostic(assemble::DiagnosticId::bands);
    assemble::DerivedValues v = bands.evaluate(system, x);
    results::BandDiagram d{take(v, bands, FieldId::conduction_band),
                           take(v, bands, FieldId::valence_band),
                           take(v, bands, FieldId::electron_fermi),
                           take(v, bands, FieldId::hole_fermi)};
    for (auto* e : {&d.conduction_eV, &d.valence_eV, &d.electron_fermi_eV, &d.hole_fermi_eV}) {
        for (double& level : *e) level *= V_T;
    }
    return d;
}

// Appends the electron and hole current through each mesh edge, times current_scale, through
// the table.
inline void append_edge_currents(const assemble::System& system, std::span<const double> x,
                                 double current_scale, std::vector<double>& jn,
                                 std::vector<double>& jp) {
    using assemble::FieldId;
    const assemble::Diagnostic& currents =
        assemble::diagnostic(assemble::DiagnosticId::edge_currents);
    assemble::DerivedValues v = currents.evaluate(system, x);
    for (double J : take(v, currents, FieldId::electron_edge_current)) {
        jn.push_back(J * current_scale);
    }
    for (double J : take(v, currents, FieldId::hole_edge_current)) {
        jp.push_back(J * current_scale);
    }
}

// The lattice and carrier temperatures in K (empty without their models), through the table.
inline void temperatures(const assemble::System& system, std::span<const double> x,
                         results::NodeFields& f) {
    using assemble::FieldId;
    const assemble::Diagnostic& t = assemble::diagnostic(assemble::DiagnosticId::temperatures);
    assemble::DerivedValues v = t.evaluate(system, x);
    f.temperature_K = take(v, t, FieldId::lattice_temperature);
    f.electron_temperature_K = take(v, t, FieldId::electron_temperature_K);
    f.hole_temperature_K = take(v, t, FieldId::hole_temperature_K);
}

// Unit 62 (QP12): the quantum potentials in eV of a density-gradient drift-diffusion state into f (left
// empty without the model).
inline void quantum_potentials(const assemble::System& system, std::span<const double> x, double V_T,
                               results::NodeFields& f) {
    const auto sn = system.dofs().slot(assemble::FieldId::quantum_potential_n);
    if (!sn) return;
    const std::size_t sp = *system.dofs().slot(assemble::FieldId::quantum_potential_p);
    const std::size_t m = system.dofs().block(), n = system.node_count();
    f.quantum_potential_n_eV.resize(n);
    f.quantum_potential_p_eV.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        f.quantum_potential_n_eV[i] = x[m * i + *sn] * V_T;
        f.quantum_potential_p_eV[i] = x[m * i + sp] * V_T;
    }
}

// The drift-diffusion state of an equilibrium result's fields (solve_equilibrium): the classical
// state of its potential, or with density gradient the state its potential and quantum potentials
// slave (Unit 62, QP8).
inline std::vector<double> state_from_equilibrium(const assemble::System& dd,
                                                  const results::NodeFields& f, double V_T) {
    const std::size_t n = dd.node_count();
    std::vector<double> psi(n);
    for (std::size_t i = 0; i < n; ++i) psi[i] = f.potential_V[i] / V_T;
    if (f.quantum_potential_n_eV.empty()) return dd.state_from_potential(psi);
    std::vector<double> ln(n), lp(n);
    for (std::size_t i = 0; i < n; ++i) {
        ln[i] = f.quantum_potential_n_eV[i] / V_T;
        lp[i] = f.quantum_potential_p_eV[i] / V_T;
    }
    return dd.state_from_quantum_potential(psi, ln, lp);
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
inline void scaled_state(const assemble::System& system, const results::NodeFields& f,
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
    // Unit 62: Lambda_n, Lambda_p in kT (0 when the fields carry none).
    if (const auto sn = system.dofs().slot(assemble::FieldId::quantum_potential_n)) {
        const std::size_t sp = *system.dofs().slot(assemble::FieldId::quantum_potential_p);
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            x[k.block * i + *sn] =
                f.quantum_potential_n_eV.empty() ? 0.0 : f.quantum_potential_n_eV[i] / scaling.V_T;
            x[k.block * i + sp] =
                f.quantum_potential_p_eV.empty() ? 0.0 : f.quantum_potential_p_eV[i] / scaling.V_T;
        }
    }
    if (energy_balance(system)) {  // Unit 52: tau_n, tau_p (0 when the fields carry none)
        const std::size_t sn = *system.dofs().slot(assemble::FieldId::electron_temperature);
        const std::size_t sp = *system.dofs().slot(assemble::FieldId::hole_temperature);
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            x[k.block * i + sn] = f.electron_temperature_K.empty()
                                      ? 0.0
                                      : (f.electron_temperature_K[i] - T0) / T0;
            x[k.block * i + sp] =
                f.hole_temperature_K.empty() ? 0.0 : (f.hole_temperature_K[i] - T0) / T0;
        }
    }
}

// The generation-recombination rates of a drift-diffusion state (Unit 45): the system's mechanisms
// (assemble::rates, already in cm^-3 s^-1) and their sum, through the derived-field table.
inline results::RateFields rate_fields(const assemble::System& system,
                                       std::span<const double> x) {
    using assemble::FieldId;
    const assemble::Diagnostic& d = assemble::diagnostic(assemble::DiagnosticId::rates);
    assemble::DerivedValues v = d.evaluate(system, x);
    results::RateFields r;
    r.srh = take(v, d, FieldId::srh_rate);
    r.auger = take(v, d, FieldId::auger_rate);
    r.radiative = take(v, d, FieldId::radiative_rate);
    r.impact = take(v, d, FieldId::impact_rate);
    r.interface_n = take(v, d, FieldId::interface_rate_n);
    r.interface_p = take(v, d, FieldId::interface_rate_p);
    r.band_to_band_n = take(v, d, FieldId::band_to_band_rate_n);
    r.band_to_band_p = take(v, d, FieldId::band_to_band_rate_p);
    r.net_n = take(v, d, FieldId::net_rate_n);
    r.net_p = take(v, d, FieldId::net_rate_p);
    return r;
}

// The fields and terminal quantities of a drift-diffusion state x into `point` (its bias and
// convergence are the caller's), in the public units; shared by sweep_bias and trace_bias. With
// `rates`, its generation-recombination rates too (Unit 45).
inline void drift_diffusion_point(const assemble::System& system,
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
    quantum_potentials(system, x, scaling.V_T, point.fields);
    // Units 23 and 52: the lattice and carrier temperatures (T0 (1 + tau), T0 the scaling's).
    temperatures(system, x, point.fields);
    if (electrothermal(system)) {  // the heat to each sink
        point.thermal_contact_heat = thermal_contact_heat(system, x);
        for (double& P : point.thermal_contact_heat) P *= scaling.V_T * current_scale;
    }
    point.bands = band_diagram(system, x, scaling.V_T);
    point.terminal_current = terminal_currents(system, x);
    for (double& I : point.terminal_current) I *= current_scale;
    point.terminal_current_resolution = terminal_current_resolution(system, x);
    for (double& I : point.terminal_current_resolution) I *= current_scale;
    point.gate_charge = gate_charges(system, x);
    for (double& Q : point.gate_charge) Q *= charges;
    point.interface_trap_charge = interface_trap_charges(system, x);
    for (double& Q : point.interface_trap_charge) Q *= charges;
    append_edge_currents(system, x, current_scale, point.edge_current_n, point.edge_current_p);
    for (const assemble::TunnelPathState& s : path_states(system, x)) {
        point.tunnel_paths.push_back(
            {s.start, s.length_cm, s.field_V_per_cm, s.rate_cm3_s, s.reached});
    }
    if (gate_tunnelling(system)) {  // Unit 51
        auto t = gate_tunnelling_currents(system, x);
        for (std::vector<double>* v : {&t.ecb, &t.hvb, &t.evb}) {
            for (double& I : *v) I *= current_scale;
        }
        point.tunnelling_ecb = std::move(t.ecb);
        point.tunnelling_hvb = std::move(t.hvb);
        point.tunnelling_evb = std::move(t.evb);
    }
}

}  // namespace NiTCAD::solve::detail
