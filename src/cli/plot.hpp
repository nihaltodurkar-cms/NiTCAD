// nitcad plot (Unit 25a, VR6): the command's entry point, called by main with the whole argument
// list (argv[1] is "plot"). Returns the exit status: 0 written, 1 a failure (an input that does not
// load, a result that does not belong to the device's mesh, a drawing error), 64 usage.
#pragma once

namespace NiTCAD::cli {

[[nodiscard]] int plot(int argc, char** argv);
[[nodiscard]] const char* plot_usage_text();

}  // namespace NiTCAD::cli
