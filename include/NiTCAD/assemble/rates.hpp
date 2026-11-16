// The generation-recombination rates of a system's state, per mechanism (Unit 45, DECISIONS.md
// RS1; Unit 57 GE3).
#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "NiTCAD/assemble/system.hpp"

namespace NiTCAD::assemble {

// Each generation-recombination mechanism at state x, per node, in cm^-3 s^-1, read from its own
// equation term: the term's contribution to the node's electron and hole continuity rows divided
// by the node's control volume, so the rates are the residual's, operand for operand. Sign: net
// recombination R - G (a generation is negative). A mechanism is listed only when the system has
// its term; the local and nonlocal band-to-band terms add into one entry. Interfaces: their
// edges' carrier exchange (surface recombination and trap capture) over the node's control
// volume, an effective volume rate at interface nodes; read with steady traps (no time step).
// Precondition (NITCAD_EXPECTS): x has unknowns() entries.
enum class RateMechanism : std::uint8_t { srh, auger, radiative, interfaces, impact, band_to_band };
struct MechanismRate {
    RateMechanism mechanism;
    std::vector<double> electron;  // per node, as the electron row sees it
    std::vector<double> hole;
};
[[nodiscard]] std::vector<MechanismRate> rates(const System& system, std::span<const double> x);

}  // namespace NiTCAD::assemble
