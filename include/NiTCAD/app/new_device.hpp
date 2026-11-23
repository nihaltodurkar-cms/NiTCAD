// New device (Unit 26n slice 1, DT4, DECISIONS.md): the device library's families (examples/devices), a deck's
// main dimensions (its first [parameters] group), a thumbnail of its regions for the dialog, and the copy Create
// writes. No Win32.
#pragma once

#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/examples.hpp"

namespace NiTCAD::app::devices {

struct Family {
    Example example;
    int dimension = 0;  // the built mesh's, 0 when the deck does not build
};

// The decks of `dir` in list_examples order, each with its dimension.
[[nodiscard]] std::vector<Family> list_families(const std::filesystem::path& dir);

// The names of the rows of the first group in [parameters] (the rows before the first `# group:` line when it has
// none of its own); empty without [parameters] or when they do not resolve.
[[nodiscard]] std::vector<std::string> main_dimensions(std::string_view deck_text);

// The regions as boxes in x-y, cm, in deck order (a later box drawn over an earlier one, as the deck assigns
// regions): a region without a box is the whole domain. A 3D deck is cut at its domain's middle z; a 1D deck is a
// bar from y 0 to 1. A deck that does not build: its message, no boxes.
struct ThumbBox {
    double x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;
    std::string material;
};
struct Thumbnail {
    double x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;  // the domain
    std::vector<ThumbBox> boxes;
    std::string error;
};
[[nodiscard]] Thumbnail thumbnail(std::string_view deck_text);

// The deck written to <folder>/<name>.deck (".deck" added when the name has none; the folder made when missing);
// its path; each [output] file is renamed <name>.<its extension>, so two devices of one family in one folder keep
// their own results. Refused, nothing written: an empty name, a name holding \ / : * ? " < > |, a file of that name already
// there ("<file> is already in <folder>: choose another name"), or a write that fails.
[[nodiscard]] std::expected<std::filesystem::path, std::string> create_device(std::string_view deck_text,
                                                                              const std::filesystem::path& folder,
                                                                              std::string name);

}  // namespace NiTCAD::app::devices
