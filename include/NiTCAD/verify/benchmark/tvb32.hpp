// TVB-32, the consistency suite C1-C12 on real solves of the three TVR-1 decks (DECISIONS Unit 55 11.6, TVR1-Q5,
// the TVB-32 rulings and RC-Q5 of 11.13). The decks are embedded at build time (verify_decks.hpp, generated from
// examples/pn_diode_1d.deck, decks/schottky_diode_1d.deck, examples/moscap_cv.deck) and read by deck::parse, outputs
// cleared. Every check whose inputs the run supplies is applied at every operating point: C1-C5 and C12 at each DC
// point, C1 and C3 on the equilibrium, C7 and C8 at each C-V point and at the zero-bias AC point added through the
// API. C6, C9, C10 and C11 need inputs no deck produces (listed in the benchmark's notes).
#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/verify/consistency/checks.hpp"

namespace NiTCAD::verify::benchmark {

struct Tvb32Entry {
    std::string deck;  // "pn_diode_1d", "schottky_diode_1d", "moscap_cv"
    std::string at;    // "equilibrium", "bias (0.1) V", "added AC point, bias (0) V"
    bool zero_bias = false;  // every bias exactly 0 and one thermal temperature (the C1 predicates)
    std::size_t contacts = 0;
    consistency::CheckResult result;
};

struct Tvb32Run {
    std::vector<Tvb32Entry> entries;
    std::optional<std::string> run_error;
};

[[nodiscard]] Tvb32Run tvb32_entries();

}  // namespace NiTCAD::verify::benchmark
