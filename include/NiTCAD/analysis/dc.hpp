// DC extraction from current-voltage curves: transistor, diode and breakdown figures
// (ARCHITECTURE.md section 11, Unit 24; legacy gui/services/characterization.py,
// gui/services/sweep_derived.py and tests/test_validation_2d.py _extract_vth_max_gm).
//
// The transistor extractors take a transfer curve (gate bias against drain current) whose current
// rises with the gate, as an n-channel device's drain current does; scale a p-channel curve by -1
// (curve.hpp scaled) and read its gate biases with the opposite sign.
// A window from_V to to_V holds the points curve.hpp slice would take; they must be consecutive.
#pragma once

#include <expected>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/base/error.hpp"

namespace NiTCAD::analysis {

// The largest dy/dx (curve.hpp derivative), in A cm^(D-3) / V; its window the point and its
// neighbours. Errors: x not strictly monotone; the largest not positive.
[[nodiscard]] std::expected<Extraction, base::Error> max_transconductance(const Curve& transfer);

// The threshold voltage by linear extrapolation at the largest transconductance: the tangent there
// reaches zero current at V - I / g_m, less drain_V / 2 (the linear-region correction), as the
// legacy _extract_vth_max_gm. Errors: those of max_transconductance.
[[nodiscard]] std::expected<Extraction, base::Error> threshold_max_gm(const Curve& transfer,
                                                                      double drain_V);

// The gate bias where the current first reaches criterion (positive), interpolated in ln I.
// OLD / NEW / REASON: OLD, the legacy interpolated linearly in I. NEW, in ln I. REASON, below
// threshold the current is exponential in the gate bias, which ln I interpolates exactly.
// Errors: criterion not finite and positive; the current never reaches it.
[[nodiscard]] std::expected<Extraction, base::Error> threshold_constant_current(
    const Curve& transfer, double criterion);

// The subthreshold swing in mV/decade: the steepest segment between consecutive points where both
// currents are positive, above their resolution, and at most ceiling_fraction of the largest
// current; the window the segment's two points.
// OLD / NEW / REASON: OLD, the legacy took the 25th percentile of the local inverse slopes over the
// lowest decades. NEW, the steepest segment. REASON, the swing is defined by the steepest part of
// the curve; a percentile depends on how the points are spaced.
// Errors: ceiling_fraction not in (0, 1]; no such segment with a rising current.
[[nodiscard]] std::expected<Extraction, base::Error> subthreshold_swing(
    const Curve& transfer, double ceiling_fraction = 1e-2);

// Drain-induced barrier lowering in V/V: the constant-current thresholds at two drain biases,
// (V_th(low) - V_th(high)) / (high_drain_V - low_drain_V); the window that of the high-drain
// threshold. Errors: the two drain biases equal; those of threshold_constant_current.
[[nodiscard]] std::expected<Extraction, base::Error> dibl(const Curve& low_drain,
                                                          const Curve& high_drain,
                                                          double low_drain_V,
                                                          double high_drain_V, double criterion);

// I(on_V) / I(off_V), each interpolated in ln |I|.
// OLD / NEW / REASON: OLD, the legacy divided the largest current of the sweep by its smallest.
// NEW, the currents at named biases. REASON, the extremes of a sweep depend on its range, not on
// the device. Errors: those of value_at; the two currents of opposite signs.
[[nodiscard]] std::expected<Extraction, base::Error> on_off_ratio(const Curve& transfer,
                                                                  double on_V, double off_V);

// The output conductance g_ds in A cm^(D-3) / V: the least-squares slope over the last
// fit_fraction of the curve's bias span (the saturation tail). Errors: fit_fraction not in (0, 1];
// fewer than 3 points in the tail; x not strictly monotone.
[[nodiscard]] std::expected<Extraction, base::Error> output_conductance(
    const Curve& output, double fit_fraction = 0.25);

struct DiodeFit {
    Extraction ideality;            // dimensionless
    Extraction saturation_current;  // A cm^(D-3)
};

// The least-squares line ln I = ln I_s + V / (n V_T) over the points with from_V <= V <= to_V.
// Errors: temperature not positive; fewer than 2 points in the window; a current there not
// positive; the slope not positive.
[[nodiscard]] std::expected<DiodeFit, base::Error> diode_fit(const Curve& forward, double from_V,
                                                             double to_V, double temperature_K);

// n = 1 / (V_T d ln I / dV) at every point. Errors: a current not positive; x not strictly
// monotone; temperature not positive.
[[nodiscard]] std::expected<std::vector<double>, base::Error> local_ideality(
    const Curve& forward, double temperature_K);

struct SeriesResistanceFit {
    Extraction ideality;    // dimensionless
    Extraction resistance;  // Ohm cm^(3-D)
};

// With V = n V_T ln(I / I_s) + I R_s, dV / d ln I = n V_T + R_s I: the least-squares line of
// dV / d ln I (curve.hpp derivative of V against ln I) against I over the points with
// from_V <= V <= to_V, leaving out the curve's two end points (their derivative is one-sided).
// Errors: temperature not positive; fewer than 3 points left in the window; a current there or at
// a neighbour not positive; ln I not strictly monotone there.
[[nodiscard]] std::expected<SeriesResistanceFit, base::Error> series_resistance(
    const Curve& forward, double from_V, double to_V, double temperature_K);

// The bias where |I| first reaches limit (positive) along the curve, interpolated in ln |I|.
// Errors: limit not finite and positive; |I| never reaches it.
[[nodiscard]] std::expected<Extraction, base::Error> breakdown_at_current(const Curve& curve,
                                                                          double limit);

// The first turning point of the bias along a traced curve (trace_bias): where the bias, after
// moving away from its start, first turns back (snapback). Refined by the parabola V(ln |I|)
// through the extreme point and its neighbours when ln |I| rises (or falls) strictly through them
// and the vertex lies between them; else the extreme point itself. Errors: the bias never turns.
[[nodiscard]] std::expected<Extraction, base::Error> turning_point(const Curve& curve);

}  // namespace NiTCAD::analysis
