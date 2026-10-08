#include "NiTCAD/analysis/two_port.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "common.hpp"

namespace NiTCAD::analysis {

using detail::invalid;

namespace {

base::Error singular(const char* message) {
    return {base::ErrorCode::singular_system, message, std::nullopt};
}

TwoPort multiply(const TwoPort& a, const TwoPort& b) noexcept {
    return {a.p11 * b.p11 + a.p12 * b.p21, a.p11 * b.p12 + a.p12 * b.p22,
            a.p21 * b.p11 + a.p22 * b.p21, a.p21 * b.p12 + a.p22 * b.p22};
}

// (1 - m)(1 + m)^-1: the Cayley map, its own inverse.
std::expected<TwoPort, base::Error> cayley(const TwoPort& m) {
    const TwoPort plus{1.0 + m.p11, m.p12, m.p21, 1.0 + m.p22};
    const TwoPort minus{1.0 - m.p11, -m.p12, -m.p21, 1.0 - m.p22};
    const auto inv = inverse(plus);
    if (!inv) return inv;
    return multiply(minus, *inv);
}

}  // namespace

std::expected<TwoPort, base::Error> admittance_two_port(const results::SmallSignal& small_signal,
                                                        std::size_t point, std::size_t frequency,
                                                        std::size_t input, std::size_t output) {
    if (point >= small_signal.points.size()) return std::unexpected(invalid("point out of range"));
    if (frequency >= small_signal.frequency_Hz.size() ||
        frequency >= small_signal.points[point].admittance.size()) {
        return std::unexpected(invalid("frequency out of range"));
    }
    if (input >= small_signal.contacts || output >= small_signal.contacts) {
        return std::unexpected(invalid("contact out of range"));
    }
    if (input == output) return std::unexpected(invalid("the input and output contacts are equal"));
    const auto Y = [&](std::size_t i, std::size_t j) {
        return small_signal.admittance(point, frequency, i, j);
    };
    return TwoPort{Y(input, input), Y(input, output), Y(output, input), Y(output, output)};
}

TwoPort scaled(const TwoPort& m, double factor) noexcept {
    return {m.p11 * factor, m.p12 * factor, m.p21 * factor, m.p22 * factor};
}

std::expected<TwoPort, base::Error> inverse(const TwoPort& m) {
    const Complex det = m.p11 * m.p22 - m.p12 * m.p21;
    if (det == 0.0) return std::unexpected(singular("the determinant is zero"));
    return TwoPort{m.p22 / det, -m.p12 / det, -m.p21 / det, m.p11 / det};
}

std::expected<TwoPort, base::Error> y_to_h(const TwoPort& y) {
    if (y.p11 == 0.0) return std::unexpected(singular("y11 is zero"));
    const Complex det = y.p11 * y.p22 - y.p12 * y.p21;
    return TwoPort{1.0 / y.p11, -y.p12 / y.p11, y.p21 / y.p11, det / y.p11};
}

std::expected<TwoPort, base::Error> h_to_y(const TwoPort& h) {
    if (h.p11 == 0.0) return std::unexpected(singular("h11 is zero"));
    return y_to_h(h);  // the map is an involution
}

std::expected<TwoPort, base::Error> y_to_s(const TwoPort& y_siemens, double reference_ohm) {
    if (!(std::isfinite(reference_ohm) && reference_ohm > 0.0)) {
        return std::unexpected(invalid("the reference impedance is not finite and positive"));
    }
    return cayley(scaled(y_siemens, reference_ohm));
}

std::expected<TwoPort, base::Error> s_to_y(const TwoPort& s, double reference_ohm) {
    if (!(std::isfinite(reference_ohm) && reference_ohm > 0.0)) {
        return std::unexpected(invalid("the reference impedance is not finite and positive"));
    }
    const auto y = cayley(s);
    if (!y) return y;
    return scaled(*y, 1.0 / reference_ohm);
}

std::expected<Complex, base::Error> current_gain(const TwoPort& y) {
    if (y.p11 == 0.0) return std::unexpected(singular("y11 is zero"));
    return y.p21 / y.p11;
}

std::expected<double, base::Error> unilateral_gain(const TwoPort& y) {
    const double d = 4.0 * (y.p11.real() * y.p22.real() - y.p12.real() * y.p21.real());
    if (!(d > 0.0)) return std::unexpected(invalid("U is undefined: its denominator is not positive"));
    return std::norm(y.p21 - y.p12) / d;
}

double stability_factor(const TwoPort& y) noexcept {
    const double m = std::abs(y.p12 * y.p21);
    if (m == 0.0) return std::numeric_limits<double>::infinity();
    return (2.0 * y.p11.real() * y.p22.real() - (y.p12 * y.p21).real()) / m;
}

std::expected<double, base::Error> maximum_gain(const TwoPort& y) {
    if (y.p12 == 0.0) return std::unexpected(invalid("y12 is zero: the two-port is unilateral"));
    const double msg = std::abs(y.p21 / y.p12);
    const double k = stability_factor(y);
    if (k < 1.0) return msg;
    return msg * (k - std::sqrt(k * k - 1.0));
}

std::expected<Extraction, base::Error> unity_gain_frequency(std::span<const double> frequency_Hz,
                                                            std::span<const double> gain) {
    if (frequency_Hz.size() != gain.size()) {
        return std::unexpected(invalid("the frequencies and gains differ in size"));
    }
    std::vector<std::size_t> index;
    for (std::size_t k = 0; k < gain.size(); ++k) {
        const double f = frequency_Hz[k], g = gain[k];
        if (!(std::isfinite(f) && f > 0.0 && std::isfinite(g) && g > 0.0)) continue;
        if (!index.empty() && !(f > frequency_Hz[index.back()])) {
            return std::unexpected(invalid("the frequencies do not increase", k));
        }
        index.push_back(k);
    }
    if (index.empty()) return std::unexpected(invalid("no usable sample"));
    if (gain[index.front()] < 1.0) {
        return std::unexpected(invalid("the gain is below 1 at the lowest frequency"));
    }
    for (std::size_t j = 0; j + 1 < index.size(); ++j) {
        const std::size_t a = index[j], b = index[j + 1];
        if (gain[a] == 1.0) {
            return Extraction{.value = frequency_Hz[a], .unit = "Hz", .window = {a, a}};
        }
        if (!(gain[b] <= 1.0)) continue;
        const double la = std::log(gain[a]), lb = std::log(gain[b]);
        const double t = la / (la - lb);
        const double lf = std::log(frequency_Hz[a]) +
                          t * (std::log(frequency_Hz[b]) - std::log(frequency_Hz[a]));
        return Extraction{.value = std::exp(lf), .unit = "Hz", .window = {a, b}};
    }
    const std::size_t last = index.back();
    return Extraction{.value = frequency_Hz[last] * gain[last],
                      .unit = "Hz",
                      .window = {last, last},
                      .extrapolated = gain[last] != 1.0};
}

namespace {

template <class Gain>
std::expected<Extraction, base::Error> unity_of(const results::SmallSignal& small_signal,
                                                std::size_t point, std::size_t input,
                                                std::size_t output, Gain gain) {
    std::vector<double> g(small_signal.frequency_Hz.size(),
                          std::numeric_limits<double>::quiet_NaN());
    for (std::size_t k = 0; k < g.size(); ++k) {
        const auto y = admittance_two_port(small_signal, point, k, input, output);
        if (!y) return std::unexpected(y.error());
        if (small_signal.frequency_Hz[k] > 0.0) g[k] = gain(*y);
    }
    return unity_gain_frequency(small_signal.frequency_Hz, g);
}

}  // namespace

std::expected<Extraction, base::Error> transition_frequency(
    const results::SmallSignal& small_signal, std::size_t point, std::size_t input,
    std::size_t output) {
    return unity_of(small_signal, point, input, output, [](const TwoPort& y) {
        const auto h = current_gain(y);
        return h ? std::abs(*h) : std::numeric_limits<double>::quiet_NaN();
    });
}

std::expected<Extraction, base::Error> maximum_oscillation_frequency(
    const results::SmallSignal& small_signal, std::size_t point, std::size_t input,
    std::size_t output) {
    return unity_of(small_signal, point, input, output, [](const TwoPort& y) {
        const auto u = unilateral_gain(y);
        return u ? std::sqrt(*u) : std::numeric_limits<double>::quiet_NaN();
    });
}

}  // namespace NiTCAD::analysis
