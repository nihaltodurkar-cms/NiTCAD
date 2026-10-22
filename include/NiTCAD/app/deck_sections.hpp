// The Device, Mesh, Physics and Simulation menus (Unit 26, DECISIONS.md VQ13, WG7): each opens the deck
// at its sections and lists keys to insert. The deck stays the only source of truth: a menu edits its
// text, nothing else.
//
// The keys listed are those with a default (deck/build.hpp), inserted with that default, so an
// inserted key gives a deck that builds with the same meaning; the user then edits the value. Keys
// without a default (lengths, contacts, profiles, analysis points, outputs) are not listed: their
// sections are opened for typing.
#pragma once

#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace NiTCAD::app {

enum class DeckMenu : unsigned char { device, mesh, physics, simulation };

[[nodiscard]] std::string_view to_string(DeckMenu m) noexcept;  // "Device", "Mesh", ...

struct KeySpec {
    std::string_view key;
    std::string_view value;  // the default
    std::string_view help;
};

struct SectionSpec {
    std::string_view kind;   // the header's kind: "device", "region", "solver", ...
    std::string_view title;  // for the menu
    std::span<const KeySpec> keys;
};

// The section kinds a menu opens, in order.
[[nodiscard]] std::span<const SectionSpec> menu_sections(DeckMenu m) noexcept;
[[nodiscard]] const SectionSpec* find_section_spec(std::string_view kind) noexcept;

// A section header found in a deck's text (a line whose first non-blank character is '['; the text
// need not build).
struct SectionLocation {
    int line = 0;       // 1-based, the header's line
    int last_line = 0;  // the section's last entry line (the header's when it has none)
    std::string kind;
    std::string name;
    struct Entry {
        std::string key;
        std::string value;  // trimmed, without its comment (a continued value: its first line)
        int line = 0;
    };
    std::vector<Entry> entries;  // in order
};
[[nodiscard]] std::vector<SectionLocation> find_sections(std::string_view text);
// The first section of a kind, if the deck has one.
[[nodiscard]] std::optional<SectionLocation> find_section(std::string_view text, std::string_view kind);

// The keys of the section's spec it may take and does not have: a generated mesh's options by its type
// (min_angle in 2D, radius_edge in 3D, max_nodes both), none for a tensor mesh.
[[nodiscard]] std::vector<KeySpec> insertable_keys(const SectionLocation& section);

// The text with "key = value" inserted as the section's last entry (after its last entry line), and
// the inserted line's number. A deck without the section gets it at its end, except an analysis
// section ([transient], ...), which would make a second analysis: nothing then.
struct Insertion {
    std::string text;
    int line = 0;
};
[[nodiscard]] std::optional<Insertion> insert_key(std::string_view text, std::string_view kind, const KeySpec& key);

}  // namespace NiTCAD::app
