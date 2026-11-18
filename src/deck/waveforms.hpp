// Waveforms of the deck's [transient] section (private; deck/build.hpp, Unit 43, DECISIONS.md AN2,
// AQ1): a solve::Waveform factory written as a call, its arguments numbers in V, s, Hz and rad;
// times and frequencies may carry a unit (2ns, 10MHz; Unit 61, RD1).
// Each reader returns the waveform or a message without the location, as values.hpp.
#pragma once

#include <expected>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "NiTCAD/solve/waveform.hpp"
#include "values.hpp"

namespace NiTCAD::deck::detail {

// constant(v), step(v0, v1, t), ramp(v0, v1, t0, t1), pulse(base, pulse, t_start, width[, rise,
// fall]), pwl(t v, t v, ...), sine(offset, amplitude, frequency[, delay, phase]).
inline std::expected<solve::Waveform, std::string> waveform(std::string_view s) {
    s = trimmed(s);
    const std::size_t open = s.find('(');
    if (open == std::string_view::npos || !s.ends_with(")")) {
        return std::unexpected("'" + std::string(s) +
                               "' is not a waveform: constant(...), step(...), ramp(...), "
                               "pulse(...), pwl(...) or sine(...)");
    }
    const std::string_view name = trimmed(s.substr(0, open));
    const auto args = items(s.substr(open + 1, s.size() - open - 2));
    const auto lifted = [](std::expected<solve::Waveform, base::Error> w)
        -> std::expected<solve::Waveform, std::string> {
        if (!w) return std::unexpected(w.error().message);
        return std::move(*w);
    };
    if (name == "pwl") {
        std::vector<std::pair<double, double>> corners;
        for (const std::string_view a : args) {
            const auto t = tokens(a);
            if (t.size() != 2) return std::unexpected("pwl takes corners 't v', separated by commas");
            const auto time = seconds(t[0]);
            const auto value = voltage(t[1]);
            if (!time) return std::unexpected(time.error());
            if (!value) return std::unexpected(value.error());
            corners.emplace_back(*time, *value);
        }
        return lifted(solve::Waveform::piecewise_linear(std::move(corners)));
    }
    // Each argument by its position's kind (Unit 61, RD1): V a voltage, T a time (s, ms, ..., fs),
    // F a frequency (Hz, ..., THz), R a phase in rad (a plain number).
    const std::string_view kinds = name == "constant" ? "V"
                                   : name == "step"   ? "VVT"
                                   : name == "ramp"   ? "VVTT"
                                   : name == "pulse"  ? "VVTTTT"
                                   : name == "sine"   ? "VVFTR"
                                                      : "";
    std::vector<double> v;
    for (std::size_t k = 0; k < args.size(); ++k) {
        const char kind = k < kinds.size() ? kinds[k] : 'R';  // too many arguments: arity's message below
        const auto x = kind == 'V' ? voltage(args[k])
                       : kind == 'T' ? seconds(args[k])
                       : kind == 'F' ? hertz(args[k])
                                     : number(args[k]);
        if (!x) return std::unexpected(x.error());
        v.push_back(*x);
    }
    const auto arity = [&](std::size_t lo, std::size_t hi, const char* form)
        -> std::expected<void, std::string> {
        if (v.size() < lo || v.size() > hi) return std::unexpected(std::string(name) + " takes " + form);
        return {};
    };
    if (name == "constant") {
        if (auto a = arity(1, 1, "(v)"); !a) return std::unexpected(a.error());
        return solve::Waveform::constant(v[0]);
    }
    if (name == "step") {
        if (auto a = arity(3, 3, "(v0, v1, t)"); !a) return std::unexpected(a.error());
        return lifted(solve::Waveform::step(v[0], v[1], v[2]));
    }
    if (name == "ramp") {
        if (auto a = arity(4, 4, "(v0, v1, t0, t1)"); !a) return std::unexpected(a.error());
        return lifted(solve::Waveform::ramp(v[0], v[1], v[2], v[3]));
    }
    if (name == "pulse") {
        if (auto a = arity(4, 6, "(base, pulse, t_start, width[, rise, fall])"); !a) {
            return std::unexpected(a.error());
        }
        if (v.size() == 5) return std::unexpected("pulse takes both rise and fall, or neither");
        return lifted(v.size() == 6 ? solve::Waveform::pulse(v[0], v[1], v[2], v[3], v[4], v[5])
                                    : solve::Waveform::pulse(v[0], v[1], v[2], v[3]));
    }
    if (name == "sine") {
        if (auto a = arity(3, 5, "(offset, amplitude, frequency[, delay, phase])"); !a) {
            return std::unexpected(a.error());
        }
        return lifted(solve::Waveform::sine(v[0], v[1], v[2], v.size() > 3 ? v[3] : 0.0,
                                            v.size() > 4 ? v[4] : 0.0));
    }
    return std::unexpected("unknown waveform '" + std::string(name) + "'");
}

}  // namespace NiTCAD::deck::detail
