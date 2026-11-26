// Material registry (Unit 68, MAD1): one table of every material a deck can name, its kind, the
// composition arguments it takes and the sources of its parameters. `deck::build` looks names up
// here; `nitcad materials` and the GUI region picker list the same table. The existing names map
// to the existing constexpr sets unchanged (MAC1: results bit-identical).
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

enum class MaterialKind : std::uint8_t { semiconductor, insulator };

// Where a parameter set's values come from (MAD3). The existing sets carry the sources already
// named in their header comments.
struct Citation {
    std::string_view source;  // author, journal, year
    std::string_view where;   // which values, table or equation
};

using MaterialParameters = std::variant<SemiconductorParameters, InsulatorParameters>;

struct MaterialEntry {
    std::string_view name;
    MaterialKind kind;
    std::size_t arguments_min;     // composition arguments after the name
    std::size_t arguments_max;
    std::string_view usage;        // the refusal when the argument count is wrong ("" for none)
    std::string_view summary;      // one line for listings
    std::span<const Citation> sources;
    // Builds the set from the arguments (their count already checked). Errors: as the material's
    // own factory (for example algaas_parameters).
    std::expected<MaterialParameters, base::Error> (*make)(std::span<const double> arguments);
};

// Every registry entry, in listing order.
[[nodiscard]] std::span<const MaterialEntry> material_registry() noexcept;

// The entry named `name`, or nullptr.
[[nodiscard]] const MaterialEntry* find_material(std::string_view name) noexcept;

// A parameter of a set by its C++ member path (DK7: eps_r, electron_mobility.mu_max, ...), or
// nullptr. The density-of-states mass fits are not reachable this way.
[[nodiscard]] double* semiconductor_parameter(SemiconductorParameters& p, std::string_view path) noexcept;
[[nodiscard]] double* insulator_parameter(InsulatorParameters& p, std::string_view path) noexcept;

// Every member path of the set a deck can set, with its value, in that order (`nitcad materials`).
[[nodiscard]] std::vector<std::pair<std::string_view, double>> parameter_values(const MaterialParameters& p);

}  // namespace NiTCAD::physics
