// Contact bias as a function of time for transient runs (ARCHITECTURE.md section 11, Unit 21;
// legacy transient.py StepWaveform, RampWaveform, PulseWaveform, ConstantWaveform).
//
// A waveform is piecewise linear or a sine:
// - piecewise linear: corners (t_k, v_k) with non-decreasing times, linear between them, held at
//   the first value before the first corner and at the last after the last one. Two corners at one
//   time are a jump: the value just before it is the first, from it on the second (the waveform is
//   right-continuous). The legacy step, ramp and pulse are piecewise linear (step, ramp, pulse).
// - sine (SPICE SIN): v0 + va sin(phase) before the delay, v0 + va sin(2 pi f (t - delay) + phase)
//   from it on.
// The breakpoints are the corner times and the sine's delay: there the waveform or its slope is
// not smooth, so a transient run lands a step on each and restarts its integrator there.
#pragma once

#include <cstdint>
#include <expected>
#include <utility>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::solve {

class Waveform {
public:
    enum class Kind : std::uint8_t { piecewise_linear, sine };

    // A constant bias (one corner). Errors: none; a value not finite is rejected when the run
    // checks the waveform.
    [[nodiscard]] static Waveform constant(double value_V);
    // Errors (invalid_input): no corners; a time or value not finite; times decreasing; more than
    // two corners at one time (context index: the corner).
    [[nodiscard]] static std::expected<Waveform, base::Error> piecewise_linear(
        std::vector<std::pair<double, double>> corners);
    // v0 before t, v1 from t on (legacy StepWaveform). Errors: as piecewise_linear.
    [[nodiscard]] static std::expected<Waveform, base::Error> step(double v0, double v1, double t);
    // v0 up to t0, linear to v1 at t1, v1 after (legacy RampWaveform). Errors: t1 <= t0, and as
    // piecewise_linear.
    [[nodiscard]] static std::expected<Waveform, base::Error> ramp(double v0, double v1, double t0,
                                                                   double t1);
    // v_base, rising linearly from t_start over `rise` to v_pulse, held for `width`, falling over
    // `fall` back to v_base; a zero rise or fall is a jump (legacy PulseWaveform: both zero).
    // Errors: width not positive, rise or fall negative, and as piecewise_linear.
    [[nodiscard]] static std::expected<Waveform, base::Error> pulse(double v_base, double v_pulse,
                                                                    double t_start, double width,
                                                                    double rise = 0.0,
                                                                    double fall = 0.0);
    // Errors (invalid_input): a parameter not finite, frequency not positive.
    [[nodiscard]] static std::expected<Waveform, base::Error> sine(double offset_V,
                                                                   double amplitude_V,
                                                                   double frequency_Hz,
                                                                   double delay_s = 0.0,
                                                                   double phase_rad = 0.0);

    [[nodiscard]] Kind kind() const noexcept { return kind_; }
    // The value at t (right-continuous) and the limit from below.
    [[nodiscard]] double value(double t_s) const noexcept;
    [[nodiscard]] double left_value(double t_s) const noexcept;
    // The breakpoints in (t0, t1], ascending, each once.
    [[nodiscard]] std::vector<double> breakpoints(double t0_s, double t1_s) const;
    // The defining numbers, for the run record: the corners as (t, v) pairs, or the sine's
    // (offset, amplitude, frequency, delay, phase).
    [[nodiscard]] std::vector<double> parameters() const;

private:
    Waveform() = default;
    Kind kind_ = Kind::piecewise_linear;
    std::vector<std::pair<double, double>> corners_;
    double offset_ = 0.0, amplitude_ = 0.0, frequency_ = 0.0, delay_ = 0.0, phase_ = 0.0;
};

}  // namespace NiTCAD::solve
