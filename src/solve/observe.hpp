// Private to the solve layer: the results of each public solve, handed to the assembler's test
// observer (assemble/observer.hpp; Unit 40, DECISIONS.md FG2) when one is set, so the assembly
// trace also records every field, current, charge, admittance and convergence record a run
// returns. Without an observer nothing is copied.
#pragma once

#include <bit>
#include <complex>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "NiTCAD/assemble/observer.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/small_signal.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/results/transient.hpp"

namespace NiTCAD::solve::detail {

using assemble::detail::observe;

inline void observe_count(std::string_view tag, std::size_t n) {
    observe(tag, static_cast<double>(n));
}

inline void observe(std::string_view tag, const results::ConvergenceRecord& c) {
    std::vector<double> v;
    for (const results::IterationRecord& r : c.iterations) {
        v.insert(v.end(), {static_cast<double>(r.iteration), r.update, r.residual, r.step});
    }
    v.push_back(c.converged ? 1.0 : 0.0);
    v.push_back(c.smallest_pivot_ratio.value_or(-1.0));
    observe(tag, v);
}

inline void observe(const results::NodeFields& f) {
    observe("result.potential", f.potential_V);
    observe("result.n", f.n_cm3);
    observe("result.p", f.p_cm3);
    observe("result.temperature", f.temperature_K);
}

inline void observe(const results::BandDiagram& b) {
    observe("result.conduction", b.conduction_eV);
    observe("result.valence", b.valence_eV);
    observe("result.electron_fermi", b.electron_fermi_eV);
    observe("result.hole_fermi", b.hole_fermi_eV);
}

inline void observe(const results::RunRecord& r) {
    observe("result.identity", std::bit_cast<double>(r.input_identity));
    std::vector<double> v;
    for (const auto& s : r.settings) v.push_back(s.second);
    observe("result.settings", v);
}

inline void observe_result(const results::EquilibriumResult& e) {
    if (assemble::detail::observer == nullptr) return;
    observe(e.fields);
    observe(e.bands);
    observe("result.gate_charge", e.gate_charge);
    observe("result.interface_trap_charge", e.interface_trap_charge);
    observe("result.convergence", e.convergence);
}

inline void observe_result(const results::BiasPoint& p) {
    if (assemble::detail::observer == nullptr) return;
    observe("result.bias", p.bias_V);
    observe(p.fields);
    observe(p.bands);
    observe("result.terminal_current", p.terminal_current);
    observe("result.terminal_current_resolution", p.terminal_current_resolution);
    observe("result.gate_charge", p.gate_charge);
    observe("result.interface_trap_charge", p.interface_trap_charge);
    observe("result.edge_current_n", p.edge_current_n);
    observe("result.edge_current_p", p.edge_current_p);
    observe("result.convergence", p.convergence);
    std::vector<double> paths;
    for (const results::TunnelPath& t : p.tunnel_paths) {
        paths.insert(paths.end(), {static_cast<double>(t.start_node), t.length_cm,
                                   t.field_V_per_cm, t.generation_cm3_s, t.reached ? 1.0 : 0.0});
    }
    observe("result.tunnel_paths", paths);
    observe_count("result.path_relocations", static_cast<std::size_t>(p.path_relocations));
    observe("result.thermal_bias", p.thermal_bias_K);
    observe("result.thermal_contact_heat", p.thermal_contact_heat);
    observe_count("result.bias_substeps", static_cast<std::size_t>(p.bias_substeps));
    observe("result.weighted_current", p.weighted_current);
    observe("result.weighted_current_resolution", p.weighted_current_resolution);
}

inline void observe_result(const results::Sweep& s) {
    if (assemble::detail::observer == nullptr) return;
    observe(s.run);
    for (const results::BiasPoint& p : s.points) observe_result(p);
    observe_count("result.stopped", s.stopped ? 1 : 0);
    if (s.unfinished) observe("result.unfinished", *s.unfinished);
}

inline void observe_result(const results::SmallSignal& s) {
    if (assemble::detail::observer == nullptr) return;
    observe(s.run);
    for (const results::SmallSignalPoint& p : s.points) {
        observe_result(p.dc);
        for (const auto& y : p.admittance) observe("result.admittance", std::span(y));
        observe("result.pivot_ratio", p.pivot_ratio);
        observe("result.backward_error", p.backward_error);
        for (const auto& r : p.resolution) observe("result.ac_resolution", r);
        for (const results::SmallSignalFields& f : p.fields) {
            observe("result.ac_frequency", f.frequency_Hz);
            observe_count("result.ac_contact", f.contact);
            observe("result.ac_potential", std::span(f.potential));
            observe("result.ac_n", std::span(f.n_cm3));
            observe("result.ac_p", std::span(f.p_cm3));
            observe("result.ac_temperature", std::span(f.temperature_K));
        }
    }
    observe_count("result.stopped", s.stopped ? 1 : 0);
    if (s.unfinished) observe("result.unfinished", *s.unfinished);
}

inline void observe_result(const results::Transient& t) {
    if (assemble::detail::observer == nullptr) return;
    observe(t.run);
    for (const results::TimePoint& p : t.points) {
        observe("result.time", std::vector<double>{p.time_s, p.step_s,
                                                   static_cast<double>(p.order), p.error_ratio});
        observe("result.bias", p.bias_V);
        observe("result.terminal_current", p.terminal_current);
        observe("result.conduction_current", p.conduction_current);
        observe("result.displacement_current", p.displacement_current);
        observe("result.contact_charge", p.contact_charge);
        observe("result.interface_trap_charge", p.interface_trap_charge);
        observe("result.thermal_contact_heat", p.thermal_contact_heat);
        observe("result.convergence", p.convergence);
    }
    for (const results::TransientSnapshot& s : t.snapshots) {
        observe("result.snapshot_time", s.time_s);
        observe_count("result.snapshot_point", s.point);
        observe(s.fields);
        observe(s.bands);
        observe("result.edge_current_n", s.edge_current_n);
        observe("result.edge_current_p", s.edge_current_p);
        observe("result.trap_occupancy", s.trap_occupancy);
    }
    observe("result.steps", std::vector<double>{static_cast<double>(t.rejected_steps),
                                                static_cast<double>(t.retraced_steps),
                                                static_cast<double>(t.lagged_steps)});
    observe_count("result.stopped", t.stopped ? 1 : 0);
    if (t.unfinished) observe("result.unfinished", *t.unfinished);
}

}  // namespace NiTCAD::solve::detail
