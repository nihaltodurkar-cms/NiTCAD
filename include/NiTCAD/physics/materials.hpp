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

// metal (Unit 68, MAD6): a metal region, an insulator set flagged metal (InsulatorParameters::metal).
enum class MaterialKind : std::uint8_t { semiconductor, insulator, metal };

// Where a parameter set's values come from (MAD3). The existing sets carry the sources already
// named in their header comments.
struct Citation {
    std::string_view source;  // author, journal, year
    std::string_view where;   // which values, table or equation
};

// The source of one parameter of a set (MAD3): its member path, the index of its Citation in the
// entry's sources, and how the value was got: "measured", "fit", "interpolated", "derived",
// "secondary" (through a review, not the primary paper) or "default" (carried, not sourced).
struct ParameterCitation {
    std::string_view path;
    std::size_t source;
    std::string_view basis;
};

// The bases a ParameterCitation may name.
// "approximate": a general value where no source in hand gives one (owner ruling 2026-11-26, Unit 68 slice 5b).
inline constexpr std::string_view citation_bases[] = {"measured", "fit", "interpolated", "derived", "secondary", "default", "approximate"};

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
    // Each parameter's source (MAD3). Empty for the sets before Unit 68, which cite `sources` for
    // the whole set (the citations already in their header comments, no new reading); every
    // material Unit 68 adds cites every parameter_values path here (MAC6).
    std::span<const ParameterCitation> parameter_sources = {};
    // Whether a run records a region of this material (RunRecord::materials). False for the names
    // before Unit 68, so their result files do not change (MAC1); true for every material it adds.
    bool recorded = false;
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

// The name of `names` nearest `word` by edit distance, when near enough to be a misspelling
// (distance at most a third of the word, at least 1; MAC13), or "" (the first on a tie).
[[nodiscard]] std::string_view nearest_name(std::string_view word, std::span<const std::string_view> names);

// A named gate or contact metal (Unit 68, MAD6, MAD14; slice 5a). A metal preset is a work
// function; `n_poly` and `p_poly` are the polysilicon electrodes, whose Fermi level sits at a band
// edge of silicon (device::GateElectrode), so they carry no number. A silicide preset (slice 5c)
// is a measured barrier on n-type silicon, for Schottky contacts on silicon only: silicide
// barriers do not follow Schottky-Mott.
enum class MetalPresetKind : std::uint8_t { metal, n_poly, p_poly, silicide };

struct MetalPreset {
    std::string_view name;
    MetalPresetKind kind;
    double work_function_eV;  // the metal's; 0 for a polysilicon or silicide preset
    Citation source;
    double barrier_eV = 0.0;  // a silicide's Phi_Bn on n-Si; 0 otherwise
};

// Every metal preset, in listing order.
[[nodiscard]] std::span<const MetalPreset> metal_presets() noexcept;

// The preset named `name`, or nullptr.
[[nodiscard]] const MetalPreset* find_metal_preset(std::string_view name) noexcept;

// A 64-bit FNV-1a digest of the set's kind, its parameter_values (paths and value bits), its
// direct-gap flag and its density-of-states mass fits: equal sets give equal hashes (MAD3).
[[nodiscard]] std::uint64_t parameter_hash(const MaterialParameters& p);

}  // namespace NiTCAD::physics
