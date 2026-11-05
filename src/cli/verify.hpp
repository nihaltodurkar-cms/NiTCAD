// `nitcad verify [--tier fast|full|nightly] [--id TVB-...] [--report <file>]` (Unit 55, TVR-1; DECISIONS 11.7,
// rulings RC-Q1-Q7 in 11.13).
//
// Runs the registry's benchmarks of the tier (every TVR-1 benchmark is fast; full and nightly run the same set,
// RC-Q6), or the one named by --id, prints a line per benchmark (its overall outcome, counts and notes) and the
// summary line, and with --report writes the report (io's container, io/verify_report.hpp) and its CSV summary at
// "<file>.csv".
// Exit status (TVD-23): 0 when no benchmark fails (finding, not_applicable and not_evaluated are counted, not
// failed; an all-not_applicable benchmark is printed "unverified"); 1 on any fail, a benchmark that could not run, or
// a report that could not be written (the reason on stderr); 64 usage (an unknown option, tier or id; the known ids
// are listed).
#pragma once

#include <span>

#include "NiTCAD/verify/benchmark/registry.hpp"

namespace NiTCAD::cli {

// argv[1] is "verify". `registry` is the shipped one (benchmark::benchmarks()) in nitcad; the VR8 test build passes
// it with a test-only failing entry (RC-Q7).
int verify(int argc, char** argv, std::span<const verify::benchmark::Benchmark> registry);

}  // namespace NiTCAD::cli
