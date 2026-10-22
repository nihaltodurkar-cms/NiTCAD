// Run comparison in the application (Unit 26b, DECISIONS.md XD6): a second result, run B, beside the
// workspace shown (run A). The difference maps and the overlaid curves are views (view_settings.hpp,
// render/compare.hpp); here are what the results browser lists: the shift of Unit 24's figures (the
// threshold's is ΔV_th) and the two run records side by side.
#pragma once

#include <string>
#include <vector>

#include "NiTCAD/app/workspace.hpp"

namespace NiTCAD::app {

// Whether B's fields can be compared node by node with A's: both devices have the same mesh (the mesh
// digest).
[[nodiscard]] bool comparable(const Workspace& a, const Workspace& b);

// The figures (extracted_figures) both runs have, by name and curve, in A's order: B - A.
struct FigureShift {
    std::string name;
    std::string curve;
    std::string unit;
    double a = 0.0;
    double b = 0.0;
    [[nodiscard]] double delta() const noexcept { return b - a; }
};
[[nodiscard]] std::vector<FigureShift> figure_shifts(const Workspace& a, const Workspace& b);

// The two records (run_record) side by side: every line of A in its order, then the lines only B has;
// a line missing from one run is empty there. `differs` when the values are not equal.
struct RecordRow {
    std::string name;
    std::string a;
    std::string b;
    bool differs = false;
};
[[nodiscard]] std::vector<RecordRow> compare_records(const Workspace& a, const Workspace& b);

}  // namespace NiTCAD::app
