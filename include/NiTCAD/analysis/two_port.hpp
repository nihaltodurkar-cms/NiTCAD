// Two-port parameters and RF figures of merit from the small-signal admittance (ARCHITECTURE.md
// section 11, Unit 24; deferred from Unit 22; legacy gui/services/characterization.py
// estimate_ft).
//
// A two-port is two contacts of the admittance matrix with every other contact at AC ground (the
// common terminal): its Y parameters are the 2 x 2 block of those two rows and columns, exactly,
// since the grounded contacts' voltages are zero. Y, Z and h are in the results' units (Y in
// S cm^(D-3)); gains and the figures of merit are ratios and need no size. S parameters relate to
// a reference impedance in ohm, so they take Y in siemens: scale first by the device's size in
// cm^(3-D) (its width in 2D, area in 1D; 1 in 3D).
#pragma once

#include <complex>
#include <cstddef>
#include <expected>
#include <span>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/small_signal.hpp"

namespace NiTCAD::analysis {

using Complex = std::complex<double>;

struct TwoPort {
    Complex p11, p12, p21, p22;
};

// Y of the two-port input (port 1) and output (port 2) at one operating point and frequency.
// Errors: point, frequency or contact out of range; input equal to output.
[[nodiscard]] std::expected<TwoPort, base::Error> admittance_two_port(
    const results::SmallSignal& small_signal, std::size_t point, std::size_t frequency,
    std::size_t input, std::size_t output);

[[nodiscard]] TwoPort scaled(const TwoPort& m, double factor) noexcept;

// Z = Y^-1 and back. Errors (singular_system): a zero determinant.
[[nodiscard]] std::expected<TwoPort, base::Error> inverse(const TwoPort& m);

// h11 = 1/y11, h12 = -y12/y11, h21 = y21/y11, h22 = det Y / y11, and back (the same map with h and
// y exchanged). Errors (singular_system): y11 (h11) zero.
[[nodiscard]] std::expected<TwoPort, base::Error> y_to_h(const TwoPort& y);
[[nodiscard]] std::expected<TwoPort, base::Error> h_to_y(const TwoPort& h);

// S = (1 - z0 Y)(1 + z0 Y)^-1 with Y in siemens, and Y = (1 - S)(1 + S)^-1 / z0. Errors:
// reference_ohm not finite and positive (invalid_input); the matrix inverted singular
// (singular_system).
[[nodiscard]] std::expected<TwoPort, base::Error> y_to_s(const TwoPort& y_siemens,
                                                         double reference_ohm = 50.0);
[[nodiscard]] std::expected<TwoPort, base::Error> s_to_y(const TwoPort& s,
                                                         double reference_ohm = 50.0);

// h21 = y21 / y11, the short-circuit current gain. Errors (singular_system): y11 zero.
[[nodiscard]] std::expected<Complex, base::Error> current_gain(const TwoPort& y);

// Mason's unilateral power gain U = |y21 - y12|^2 / (4 (Re y11 Re y22 - Re y12 Re y21)), invariant
// under lossless reciprocal embedding. Errors (invalid_input): the denominator not positive.
[[nodiscard]] std::expected<double, base::Error> unilateral_gain(const TwoPort& y);

// Rollett's stability factor k = (2 Re y11 Re y22 - Re(y12 y21)) / |y12 y21|; infinite when
// y12 y21 = 0.
[[nodiscard]] double stability_factor(const TwoPort& y) noexcept;

// The maximum available gain |y21 / y12| (k - sqrt(k^2 - 1)) when k >= 1, else the maximum stable
// gain |y21 / y12|. Errors (invalid_input): y12 zero (unilateral: use unilateral_gain).
[[nodiscard]] std::expected<double, base::Error> maximum_gain(const TwoPort& y);

// The frequency where an amplitude gain (|h21|, sqrt(U)) falls through 1, interpolated in
// ln f and ln gain between the samples around it; the window those frequency indices. Samples with
// f or the gain not finite and positive are skipped. If the gain stays above 1, extrapolated from
// the last sample at -20 dB/decade (f |gain|, the single-pole roll-off) and flagged.
// OLD / NEW / REASON: OLD, the legacy extrapolated along the last segment's slope. NEW, at
// -20 dB/decade. REASON, the last segment's slope is not yet the asymptotic one below f_T and
// throws the estimate far off; -20 dB/decade is the standard extrapolation of f_T and f_max.
// Errors (invalid_input): sizes differ; fewer than 1 usable sample; frequencies not increasing; the
// gain already below 1 at the first.
[[nodiscard]] std::expected<Extraction, base::Error> unity_gain_frequency(
    std::span<const double> frequency_Hz, std::span<const double> gain);

// f_T: unity_gain_frequency of |h21| over the run's frequencies at one point. Errors: those of
// admittance_two_port and unity_gain_frequency.
[[nodiscard]] std::expected<Extraction, base::Error> transition_frequency(
    const results::SmallSignal& small_signal, std::size_t point, std::size_t input,
    std::size_t output);

// f_max: unity_gain_frequency of sqrt(U) (frequencies where U is undefined skipped). Errors: as
// transition_frequency.
[[nodiscard]] std::expected<Extraction, base::Error> maximum_oscillation_frequency(
    const results::SmallSignal& small_signal, std::size_t point, std::size_t input,
    std::size_t output);

}  // namespace NiTCAD::analysis
