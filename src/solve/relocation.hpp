// Private to the solve layer: a steady solve with nonlocal tunnel paths (Unit 20). The paths are
// frozen geometry (assemble/tunnel_paths.hpp): after Newton converges they are re-traced at the
// solution, and while the traced set does not agree with the frozen one at the solution
// (assemble::paths_agree) the system takes it and Newton runs again from the
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

#include "NiTCAD/assemble/tunnelling.hpp"
#include "NiTCAD/base/error.hpp"

namespace NiTCAD::solve::detail {

inline constexpr int max_relocations = 8;

// `newton()` runs Newton on `system` from x (in place) and returns its result. A system without
// paths yet takes those traced at x first. `relocations` is the number of re-tracings. System is
// assemble::System (or a test's stand-in), with the free functions tunnelling, paths, set_paths,
// trace_paths and paths_agree that argument-dependent lookup finds (assemble/tunnelling.hpp).
template <class System, class Newton>
std::expected<void, base::Error> solve_relocating(System& system, std::span<double> x,
                                                  Newton newton, int& relocations) {
    relocations = 0;
    if (!tunnelling(system)) return newton();
    if (paths(system).paths.empty()) set_paths(system, trace_paths(system, x));
    std::vector<std::decay_t<decltype(paths(system))>> used;
    while (true) {
        if (auto ok = newton(); !ok) return ok;
        auto traced = trace_paths(system, x);
        if (paths_agree(system, traced, paths(system), x)) return {};
        for (const auto& p : used) {
            if (paths_agree(system, traced, p, x)) {
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
        used.push_back(paths(system));
        set_paths(system, std::move(traced));
        ++relocations;
    }
}

}  // namespace NiTCAD::solve::detail
