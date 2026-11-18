// Tunnelling in a system (Unit 57, DECISIONS.md GE3): the nonlocal band-to-band paths (Unit 20;
// system.hpp's header comment, "Band-to-band tunnelling") and the gate-tunnelling currents
// (Unit 51).
#pragma once

#include <span>
#include <vector>

#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/assemble/tunnel_paths.hpp"

namespace NiTCAD::assemble {

// Nonlocal tunnelling: whether the model is on; the paths traced at state x (tunnel_paths.hpp)
// whose rate at x is at least 1e-10 of the largest (the others would only fill the Jacobian);
// the frozen set the rows use (empty at creation), whose setting rebuilds the Jacobian pattern
// (make_jacobian and make_small_signal_matrix again); and each path's state at x.
// Preconditions (NITCAD_EXPECTS): x has unknowns() entries; the paths of set_paths are
// trace_paths's of this system.
[[nodiscard]] bool tunnelling(const System& system) noexcept;
[[nodiscard]] TunnelPaths trace_paths(const System& system, std::span<const double> x);
void set_paths(System& system, TunnelPaths paths);
[[nodiscard]] const TunnelPaths& paths(const System& system) noexcept;
[[nodiscard]] std::vector<TunnelPathState> path_states(const System& system,
                                                       std::span<const double> x);
// Whether two path sets of this system act alike at state x: every path whose rate at x is at
// least 1e-10 of the largest rate in either set has a path from the same start in the other with
// the same geometry (tunnel_paths.hpp same_geometry). Paths below that are left out: a start at
// the edge of the start criterion, or a crossing near the longest length, flips in and out of
// the set with the rounding of the state, at a rate near exp(-690) of the scale.
[[nodiscard]] bool paths_agree(const System& system, const TunnelPaths& a, const TunnelPaths& b,
                               std::span<const double> x);

// Gate tunnelling (Unit 51, GD10): whether any contact tunnels, and per contact the part of its
// terminal current each process carries at state x (scaled as terminal_currents; zero off a
// tunnelling contact). They sum to the tunnelling part of terminal_currents. Empty vectors
// without gate tunnelling. Precondition (NITCAD_EXPECTS): x has unknowns() entries.
[[nodiscard]] bool gate_tunnelling(const System& system) noexcept;
struct GateTunnellingCurrents {
    std::vector<double> ecb, hvb, evb;
};
[[nodiscard]] GateTunnellingCurrents gate_tunnelling_currents(const System& system,
                                                              std::span<const double> x);

}  // namespace NiTCAD::assemble
