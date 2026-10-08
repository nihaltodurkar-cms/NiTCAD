// A small-signal (AC) run (ARCHITECTURE.md 6.6, Unit 22): its record, and per bias point the
// operating point and the admittance matrix at each frequency; plain data in the public unit
// convention (Hz, V, cm^-3, S / cm^(3-D)).
//
// The admittance Y_ij = dI_i / dV_j is the small-signal current entering through contact i per
// volt on contact j, every other contact held at its bias (the short-circuit admittance
// parameters). I is the total current, conduction plus displacement (the rate of change of the
// contact's charge), so every column sums to zero (conservation) and so does every row (a common
// shift of all biases changes nothing). A run that is cancelled or fails part-way keeps every
// completed point; `stopped` then says why and `unfinished` holds the convergence history of the
// operating point that was being solved, if one was.
#pragma once

#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/convergence.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::results {

// The small-signal fields at one frequency for a unit signal on one contact: per mesh node, the
// complex amplitude per volt of the potential (V/V), and of the electron and hole densities
// (cm^-3 / V; 0 on insulator nodes).
struct SmallSignalFields {
    double frequency_Hz;
    std::size_t contact;  // the driven contact
    std::vector<std::complex<double>> potential;
    std::vector<std::complex<double>> n_cm3;
    std::vector<std::complex<double>> p_cm3;
};

struct SmallSignalPoint {
    BiasPoint dc;  // the operating point, as sweep_bias reports it
    // admittance[k][i * contacts + j] = Y_ij at frequency k, in S / cm^(3-D) (S/cm^2 in 1D, S/cm
    // in 2D, S in 3D); contacts in the device's order.
    std::vector<std::vector<std::complex<double>>> admittance;
    // Per frequency: the linear solve's smallest pivot ratio and its largest backward error over
    // the contacts' solves.
    std::vector<double> pivot_ratio;
    std::vector<double> backward_error;
    // resolution[k][j]: per frequency and driven contact j, a bound on how far the column's
    // currents can miss conservation (|sum_i Y_ij|), in S / cm^(3-D): the residual the linear
    // solve leaves in the continuity rows and, times omega, in the Poisson rows, plus the rounding
    // of the large terms that cancel in each current. As BiasPoint::terminal_current_resolution
    // for the DC currents, an admittance far below it is not resolved by the (psi, n, p) state.
    std::vector<std::vector<double>> resolution;
    std::vector<SmallSignalFields> fields;  // at the requested field frequencies
};

struct SmallSignal {
    RunRecord run;
    std::size_t contacts = 0;
    std::vector<double> frequency_Hz;
    std::vector<SmallSignalPoint> points;
    std::optional<base::Error> stopped;
    std::optional<ConvergenceRecord> unfinished;

    [[nodiscard]] std::complex<double> admittance(std::size_t point, std::size_t frequency,
                                                  std::size_t i, std::size_t j) const {
        return points[point].admittance[frequency][i * contacts + j];
    }
    // G_ij = Re Y_ij [S / cm^(3-D)].
    [[nodiscard]] double conductance(std::size_t point, std::size_t frequency, std::size_t i,
                                     std::size_t j) const {
        return admittance(point, frequency, i, j).real();
    }
    // C_ij = Im Y_ij / (2 pi f) [F / cm^(3-D)]; NaN at f = 0. (The capacitance of a contact to the
    // rest is C_ii; C_ij for i != j is negative for a plain capacitor between them.)
    [[nodiscard]] double capacitance(std::size_t point, std::size_t frequency, std::size_t i,
                                     std::size_t j) const {
        const double f = frequency_Hz[frequency];
        if (f == 0.0) return std::numeric_limits<double>::quiet_NaN();
        return admittance(point, frequency, i, j).imag() / (2.0 * std::numbers::pi * f);
    }
};

}  // namespace NiTCAD::results
