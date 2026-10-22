// Run comparison in the application (compare.hpp).
#include "NiTCAD/app/compare.hpp"

#include <algorithm>
#include <utility>

namespace NiTCAD::app {

bool comparable(const Workspace& a, const Workspace& b) {
    return a.device && b.device && same_mesh(a, b.device->mesh());
}

std::vector<FigureShift> figure_shifts(const Workspace& a, const Workspace& b) {
    std::vector<FigureShift> out;
    const std::vector<Figure> fb = extracted_figures(b);
    for (const Figure& f : extracted_figures(a)) {
        const auto match = std::ranges::find_if(fb, [&](const Figure& g) { return g.name == f.name && g.curve == f.curve; });
        if (match != fb.end()) out.push_back({f.name, f.curve, f.unit, f.value, match->value});
    }
    return out;
}

std::vector<RecordRow> compare_records(const Workspace& a, const Workspace& b) {
    const std::vector<Line> la = run_record(a), lb = run_record(b);
    std::vector<RecordRow> out;
    std::vector<bool> used(lb.size(), false);
    for (const Line& l : la) {
        RecordRow r{l.name, l.value, {}, true};
        for (std::size_t k = 0; k < lb.size(); ++k) {
            if (!used[k] && lb[k].name == l.name) {
                used[k] = true;
                r.b = lb[k].value;
                r.differs = r.a != r.b;
                break;
            }
        }
        out.push_back(std::move(r));
    }
    for (std::size_t k = 0; k < lb.size(); ++k) {
        if (!used[k]) out.push_back({lb[k].name, {}, lb[k].value, true});
    }
    return out;
}

}  // namespace NiTCAD::app
