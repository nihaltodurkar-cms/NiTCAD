// A test hook for refactors of the assembler (Unit 40, DECISIONS.md FE7, FG1, FG2): when an
// observer is set, DriftDiffusion's public evaluations hand it their inputs and outputs, each
// under a short tag, so a test run can record their bit patterns and a later run compare them.
// No production path sets one; unset, each evaluation pays one test of a null pointer.
#pragma once

#include <complex>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace NiTCAD::assemble {

using Observer = void (*)(std::string_view tag, std::span<const double> data);

// Sets the observer (nullptr: none). Not synchronized: set it before any evaluation runs.
void set_observer(Observer observer) noexcept;

namespace detail {

extern Observer observer;

inline void observe(std::string_view tag, std::span<const double> data) {
    if (observer != nullptr) observer(tag, data);
}
inline void observe(std::string_view tag, double value) {
    if (observer != nullptr) observer(tag, std::span<const double>(&value, 1));
}
inline void observe(std::string_view tag, std::span<const std::pair<double, double>> data) {
    if (observer != nullptr) {
        std::vector<double> flat;
        flat.reserve(2 * data.size());
        for (const auto& [a, b] : data) {
            flat.push_back(a);
            flat.push_back(b);
        }
        observer(tag, flat);
    }
}
// std::complex<double> is layout-compatible with double[2].
inline void observe(std::string_view tag, std::span<const std::complex<double>> data) {
    if (observer != nullptr) {
        observer(tag, std::span<const double>(reinterpret_cast<const double*>(data.data()),
                                              2 * data.size()));
    }
}

}  // namespace detail

}  // namespace NiTCAD::assemble
