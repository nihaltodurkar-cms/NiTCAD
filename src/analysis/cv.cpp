#include "NiTCAD/analysis/cv.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "common.hpp"

namespace NiTCAD::analysis {

using detail::below;
using detail::invalid;

std::expected<Curve, base::Error> quasi_static_capacitance(const Curve& charge) {
    auto d = derivative(charge);
    if (!d) return std::unexpected(d.error());
    return Curve{charge.x, std::move(*d), {}};
}

std::expected<Extraction, base::Error> accumulation_capacitance(const Curve& capacitance) {
    const auto it = std::max_element(capacitance.y.begin(), capacitance.y.end());
    const auto k = static_cast<std::size_t>(it - capacitance.y.begin());
    if (!(*it > 0.0)) return std::unexpected(invalid("the capacitance is nowhere positive"));
    return Extraction{.value = *it,
                      .unit = "F cm^(D-3)",
                      .window = {k, k},
                      .extrapolated = false,
                      .below_resolution = below(capacitance, k, k)};
}

std::expected<DopingProfile, base::Error> doping_profile(const Curve& capacitance,
                                                         double permittivity_F_per_cm,
                                                         double series_capacitance) {
    if (!(permittivity_F_per_cm > 0.0)) {
        return std::unexpected(invalid("the permittivity is not positive"));
    }
    if (!(series_capacitance >= 0.0)) {
        return std::unexpected(invalid("the series capacitance is negative"));
    }
    Curve inverse_square{capacitance.x, capacitance.y, {}};
    std::vector<double> depth;
    for (std::size_t k = 0; k < capacitance.y.size(); ++k) {
        const double C = capacitance.y[k];
        if (!(C > 0.0) || (series_capacitance > 0.0 && !(C < series_capacitance))) {
            return std::unexpected(
                invalid("a capacitance is not positive or not below the series capacitance", k));
        }
        const double C_d = series_capacitance > 0.0 ? 1.0 / (1.0 / C - 1.0 / series_capacitance) : C;
        inverse_square.y[k] = 1.0 / (C_d * C_d);
        depth.push_back(permittivity_F_per_cm / C_d);
    }
    const auto d = derivative(inverse_square);
    if (!d) return std::unexpected(d.error());
    DopingProfile p{std::move(depth), {}};
    for (std::size_t k = 0; k < d->size(); ++k) {
        if ((*d)[k] == 0.0) return std::unexpected(invalid("d(1/C^2)/dV is zero", k));
        p.doping_cm3.push_back(2.0 / (base::q_C * permittivity_F_per_cm * std::abs((*d)[k])));
    }
    return p;
}

double flat_band_capacitance(double oxide_capacitance, double permittivity_F_per_cm,
                             double doping_cm3, double temperature_K) {
    NITCAD_EXPECTS(oxide_capacitance > 0.0 && permittivity_F_per_cm > 0.0 && doping_cm3 > 0.0 &&
                   temperature_K > 0.0);
    const double L_D = std::sqrt(permittivity_F_per_cm * base::thermal_voltage(temperature_K) /
                                 (base::q_C * doping_cm3));
    return 1.0 / (1.0 / oxide_capacitance + L_D / permittivity_F_per_cm);
}

std::expected<Extraction, base::Error> flat_band_voltage(const Curve& capacitance,
                                                         double flat_band) {
    return crossing(capacitance, flat_band);
}

std::expected<ConductancePeak, base::Error> conductance_peak(
    const results::SmallSignal& small_signal, std::size_t point, std::size_t contact,
    double oxide_capacitance) {
    if (point >= small_signal.points.size()) return std::unexpected(invalid("point out of range"));
    if (contact >= small_signal.contacts) return std::unexpected(invalid("contact out of range"));
    if (!(oxide_capacitance > 0.0)) {
        return std::unexpected(invalid("the oxide capacitance is not positive"));
    }
    std::vector<std::size_t> index;
    std::vector<double> log_f, g;
    for (std::size_t k = 0; k < small_signal.frequency_Hz.size(); ++k) {
        const double f = small_signal.frequency_Hz[k];
        if (!(f > 0.0)) continue;
        if (!log_f.empty() && !(std::log(f) > log_f.back())) {
            return std::unexpected(invalid("the frequencies do not increase", k));
        }
        const double w = 2.0 * std::numbers::pi * f;
        const std::complex<double> y = small_signal.admittance(point, k, contact, contact);
        const std::complex<double> ys =
            1.0 / (1.0 / y - 1.0 / std::complex<double>{0.0, w * oxide_capacitance});
        index.push_back(k);
        log_f.push_back(std::log(f));
        g.push_back(ys.real() / w);
    }
    if (g.size() < 3) return std::unexpected(invalid("fewer than 3 positive frequencies"));
    const auto j = static_cast<std::size_t>(std::max_element(g.begin(), g.end()) - g.begin());
    if (j == 0 || j + 1 == g.size()) {
        return std::unexpected(invalid("the peak is not inside the frequency range", index[j]));
    }
    // The parabola ln g(t), t = ln f, through the three samples: its vertex. (A single level's
    // ln g = ln(C_it / 2) - ln cosh(t - t_peak) is a parabola to fourth order in t - t_peak; g
    // itself only to second.)
    const bool positive = g[j - 1] > 0.0 && g[j + 1] > 0.0;
    const auto level = [&](std::size_t i) { return positive ? std::log(g[i]) : g[i]; };
    const double t0 = log_f[j - 1], t1 = log_f[j], t2 = log_f[j + 1];
    const double d01 = (level(j) - level(j - 1)) / (t1 - t0);
    const double d12 = (level(j + 1) - level(j)) / (t2 - t1);
    const double a = (d12 - d01) / (t2 - t0);
    double t = t1, value = g[j];
    if (a < 0.0) {
        t = 0.5 * (t0 + t1) - d01 / (2.0 * a);
        const double v = level(j - 1) + (t - t0) * (d01 + a * (t - t1));
        value = positive ? std::exp(v) : v;
    }
    const Window w{index[j - 1], index[j + 1]};
    return ConductancePeak{
        .conductance_over_omega = {.value = value, .unit = "F cm^(D-3)", .window = w,
                                   .extrapolated = false, .below_resolution = false},
        .frequency = {.value = std::exp(t), .unit = "Hz", .window = w,
                      .extrapolated = false, .below_resolution = false}};
}

double interface_trap_density(double peak_conductance_over_omega) {
    return 2.5 * peak_conductance_over_omega / base::q_C;
}

}  // namespace NiTCAD::analysis
