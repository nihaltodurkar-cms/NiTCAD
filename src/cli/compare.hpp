// nitcad compare (Unit 60, XC3): two sweep files compared by tolerances (analysis/compare.hpp),
// called by main with the whole argument list (argv[1] is "compare"). Returns the exit status:
// 0 match, 4 differ, 1 a file that does not load or sweeps that are not comparable, 64 usage.
#pragma once

namespace NiTCAD::cli {

[[nodiscard]] int compare(int argc, char** argv);

}  // namespace NiTCAD::cli
