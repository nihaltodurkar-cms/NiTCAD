// A transient run (ARCHITECTURE.md 6.6, Unit 21): its record, one time point per accepted step and
// field snapshots at the output times; plain data in the public unit convention (s, V, cm^-3,
// A / cm^(3-D), C / cm^(3-D)).
//
// A run that is cancelled or fails part-way keeps every accepted step; `stopped` then says why and
// `unfinished` holds the convergence history of the step that was being solved (as Sweep).
#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/convergence.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::results {

struct TimePoint {
    double time_s;
    double step_s;  // the step that reached it (0 for the starting state)
    int order;      // of that step: 1 backward Euler, 2 BDF2 (0 for the starting state)
    // The step's estimated local error over the tolerance (accepted at <= 1); NaN when the step
    // had no estimate (the first step after a start or a breakpoint, or fixed steps).
    double error_ratio;
    std::vector<double> bias_V;  // per contact
    // Per contact, the current entering the device through it, in A / cm^(3-D): the total current
    // is the conduction current (carriers crossing into the device) plus the displacement current
    // (the rate of change of the contact's charge). The total currents of all contacts sum to
    // zero; a gate or electrode carries displacement current only. In the starting (steady) state
    // the displacement current is zero.
    std::vector<double> terminal_current;
    std::vector<double> conduction_current;
    std::vector<double> displacement_current;
    // Per contact, the charge on it (the displacement flux leaving its nodes; for a gate as
    // BiasPoint::gate_charge), in C / cm^(3-D).
    std::vector<double> contact_charge;
    // Per declared interface, as BiasPoint::interface_trap_charge, with the traps' occupancy of
    // this time.
    std::vector<double> interface_trap_charge;
    ConvergenceRecord convergence;  // empty for the starting state
};

// The fields at one time point.
struct TransientSnapshot {
    double time_s;
    std::size_t point;  // index in Transient::points
    NodeFields fields;
    BandDiagram bands;
    // Conduction current through each mesh edge (as BiasPoint).
    std::vector<double> edge_current_n;
    std::vector<double> edge_current_p;
    // Electron occupancy of each interface trap level on each interface edge (empty without
    // traps), in the assembler's slot order.
    std::vector<double> trap_occupancy;
};

struct Transient {
    RunRecord run;
    std::vector<TimePoint> points;
    std::vector<TransientSnapshot> snapshots;
    std::size_t rejected_steps = 0;  // by the error estimate or by Newton
    // Nonlocal tunnelling (Unit 20): steps repeated because the paths re-traced after them
    // differed, and steps accepted with the paths lagging (they changed again after the repeat;
    // the next step uses the new ones).
    std::size_t retraced_steps = 0;
    std::size_t lagged_steps = 0;
    std::optional<base::Error> stopped;
    std::optional<ConvergenceRecord> unfinished;
};

}  // namespace NiTCAD::results
