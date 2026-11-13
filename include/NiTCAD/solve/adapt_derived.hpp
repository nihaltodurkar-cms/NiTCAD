// Derived adaptation goals (Unit 17c, DECISIONS.md DG1-DG11): the threshold voltage, the largest
// transconductance and the subthreshold swing of a transfer curve, each with its sensitivities to
// the currents of its window (DG3), so the loop combines the per-point estimates and weights (DG4).
#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::solve {

// A goal's kind: a terminal current at the adaptation points (Unit 17), or a figure of the transfer
// curve (DG1).
enum class GoalKind : unsigned char { current, vt, gm, swing };

[[nodiscard]] constexpr const char* to_string(GoalKind k) noexcept {
    switch (k) {
        case GoalKind::current: return "current";
        case GoalKind::vt: return "vt";
        case GoalKind::gm: return "gm";
        case GoalKind::swing: return "swing";
    }
    return "?";
}

// A derived goal on one pass's sweep: its value (V, A cm^(D-3)/V or mV/decade), its resolution
// (DG5: sum of |c_k| r_k), its window (consecutive sweep points first..last), the point its
// certificate names (the largest transconductance's, or the swing segment's first) and
// dG/dI_k for k = first..last, the window held fixed (DG3), with respect to the goal contact's
// current as the sweep reports it (the sign of a p-channel current undone).
struct DerivedGoal {
    double value = 0.0;
    double resolution = 0.0;
    std::size_t first = 0;
    std::size_t last = 0;
    std::size_t point = 0;
    std::vector<double> sensitivity;
};

// The transfer curve of a derived goal: gate (or electrode) biases, the goal contact's currents and
// their resolutions, one per sweep point; drain_V the goal contact's fixed bias.
struct TransferCurve {
    std::vector<double> gate_V;
    std::vector<double> current;
    std::vector<double> resolution;
    double drain_V = 0.0;
};

// DG1: the figure by analysis/dc.hpp's extractors (threshold_max_gm, max_transconductance,
// subthreshold_swing with its default ceiling) and DG3's sensitivities. A curve whose largest
// |current| is negative (p-channel) is read as dc.hpp documents: currents, gate and drain biases
// times -1, so V_T is reported in the device's own sign and g_m, swing positive.
// Errors (invalid_input): kind current; those of the extractors.
[[nodiscard]] std::expected<DerivedGoal, base::Error> derived_goal(GoalKind kind, const TransferCurve& curve);

}  // namespace NiTCAD::solve
