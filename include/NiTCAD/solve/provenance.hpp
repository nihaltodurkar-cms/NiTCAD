// The provenance of a run (ARCHITECTURE.md 6.6, Unit 42; DECISIONS.md PV1-PV3): which build, on
// which platform and when. Every make_run_record attaches it; the public runs (sweep_bias,
// solve_transient, solve_small_signal, trace_bias) add their wall time when they return.
#pragma once

#include "NiTCAD/results/run.hpp"

namespace NiTCAD::solve {

// The build's facts (version, git revision and dirty state, build type, compiler; written at build
// time by cmake/build_info.cmake), Eigen's version and fixed cache sizes, the loaded MKL runtime's
// version (empty without one), the OS version, the processor's brand string and logical CPU count,
// and the time now (UTC) as the start; wall_s 0.
[[nodiscard]] results::Provenance provenance();

}  // namespace NiTCAD::solve
