// Private to the solve layer: a steady solve with nonlocal tunnel paths (Unit 20). The paths are
// frozen geometry (assemble/tunnel_paths.hpp): after Newton converges they are re-traced at the
// solution, and while the traced set does not agree with the frozen one at the solution
// (assemble::DriftDiffusion::paths_agree) the system takes it and Newton runs again from the
// solution. At most max_relocations re-tracings; a set that agrees with one used before (a cycle)
// or the limit ends the solve with non_convergence. Without the model it is the Newton run alone.
#pragma once

#include <expected>
#include <optional>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::solve::detail {

inline constexpr int max_relocations = 8;

// `newton()` runs Newton on `system` from x (in place) and returns its result. A system without
// paths yet takes those traced at x first. `relocations` is the number of re-tracings. System is
// assemble::DriftDiffusion (or a test's stand-in with its tunnelling, paths, set_paths,
// trace_paths and paths_agree).
template <class System, class Newton>
std::expected<void, base::Error> solve_relocating(System& system, std::span<double> x,
                                                  Newton newton, int& relocations) {
    relocations = 0;
    if (!system.tunnelling()) return newton();
    if (system.paths().paths.empty()) system.set_paths(system.trace_paths(x));
    std::vector<std::decay_t<decltype(system.paths())>> used;
    while (true) {
        if (auto ok = newton(); !ok) return ok;
        auto traced = system.trace_paths(x);
        if (system.paths_agree(traced, system.paths(), x)) return {};
        for (const auto& p : used) {
            if (system.paths_agree(traced, p, x)) {
                return std::unexpected(base::Error{
                    base::ErrorCode::non_convergence,
                    "the tunnel paths cycle between geometries after " +
                        std::to_string(relocations) + " re-tracings",
                    base::ErrorContext{.index = static_cast<std::size_t>(relocations),
                                       .value = std::nullopt}});
            }
        }
        if (relocations == max_relocations) {
            return std::unexpected(base::Error{
                base::ErrorCode::non_convergence,
                "the tunnel paths did not settle in " + std::to_string(max_relocations) +
                    " re-tracings",
                base::ErrorContext{.index = static_cast<std::size_t>(relocations),
                                   .value = std::nullopt}});
        }
        used.push_back(system.paths());
        system.set_paths(std::move(traced));
        ++relocations;
    }
}

}  // namespace NiTCAD::solve::detail
