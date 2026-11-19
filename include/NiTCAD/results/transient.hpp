// A transient run (ARCHITECTURE.md 6.6, Unit 21): its record, one time point per accepted step and
// field snapshots at the output times; plain data in the public unit convention (s, V, cm^-3,
// A / cm^(3-D), C / cm^(3-D)).
//
// A run that is cancelled or fails part-way keeps every accepted step; `stopped` then says why and
// `unfinished` holds the convergence history of the step that was being solved (as Sweep).
#pragma once

#include <cstddef>
#include <cstdint>
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
    // Electrothermal (Unit 23; empty with the model off): the heat leaving through each thermal
    // contact, as BiasPoint::thermal_contact_heat, at this time.
    std::vector<double> thermal_contact_heat;
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
    // Generation-recombination rates (Unit 45; steady traps); empty when not stored.
    RateFields rates;
};

// Where a transient run can be continued from (Unit 61, RD3; DECISIONS.md RQ1): the integrator's own
// state after the last accepted step, bit for bit, so a continued run repeats the uninterrupted
// run's arithmetic. Solver state, not a result: the vectors are the assembler's scaled unknowns
// and stored quantities, the history's times scaled by the assembler's time scale.
struct TransientResume {
    struct State {
        double t;  // scaled time
        std::vector<double> x, storage, traps, charges;
    };
    std::vector<State> history;  // the accepted states since the last (re)start, newest last (1 to 3)
    TimePoint point;             // the last accepted step, as Transient::points has it
    double time_s = 0.0;         // its time
    double next_step_s = 0.0;    // the step to try next
    std::size_t landing = 0;     // the next landing time (breakpoints, output times, t_end) to reach
    std::size_t attempt = 0;     // step attempts so far (the progress and error index)
    // The run's identity without its initial state (RunRecord::input_identity of the same device,
    // options and waveforms with no initial state): a continuation must be of this run.
    std::uint64_t identity = 0;
    // The last accepted state's fields, as a snapshot of it has them (what a checkpoint file shows).
    NodeFields fields;
};

struct Transient {
    RunRecord run;
    std::vector<TimePoint> points;
    std::vector<TransientSnapshot> snapshots;
    // Set once a step has been accepted, also in a run that stopped (Unit 61, RD3); a run that
    // reached t_end holds it too, with nothing left to continue.
    std::optional<TransientResume> resume;
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
