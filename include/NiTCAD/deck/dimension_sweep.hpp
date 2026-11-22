// Dimension sweeps (Unit 26n, DT6, DECISIONS.md): a deck run with [parameters] values overridden, once
// (`nitcad run <deck> --set Lg=30nm`) or once per value of one parameter (`nitcad sweep <deck> Lg=20nm,30nm`).
//
// An override is applied after reading, as substitution is: the deck's text with the parameter's line rewritten
// (its key, spacing and comment kept) and parsed again, so a run with an override is the run of a copy of the deck
// with the value written in. The deck file is never edited. A sweep writes each value's outputs to its own folder,
// `<deck stem>_sweep_<name>/<name>=<value>/` beside the deck (NQ4 (a)), each an ordinary run's files, with the names
// the deck's [output] gives them.
#pragma once

#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"

namespace NiTCAD::deck {

struct Override {
    std::string name;
    std::string value;  // as it would be written on the parameter's line
};

// `name=value`, both trimmed. Errors: no '=', an empty name or value.
[[nodiscard]] std::expected<Override, std::string> parse_override(std::string_view text);

// The deck as if each override's value were written on its parameter's [parameters] line: the text with those lines
// rewritten, parsed again with the deck's source and directory. Errors (invalid_input, the message naming the
// parameter): the deck has no parameter of that name; a length parameter given a value that is not a length, a number
// parameter given one that is not a number; a value on a continued line; the parameters with the values do not
// resolve.
[[nodiscard]] std::expected<Deck, base::Error> with_overrides(const Deck& deck, const std::vector<Override>& overrides);

// A sweep's values as they will be written: `a, b, c` (each kept as written) or `first : last : step` (first + k step
// up to last, each written in first's unit, shortest form). Each a literal length (with nm, um or cm) when `length`,
// else a literal number. Errors: an empty list, a value of the other kind, a reference, a bad range, a value twice
// (two values would share a folder), more than 1000 values.
[[nodiscard]] std::expected<std::vector<std::string>, std::string> sweep_values(std::string_view list, bool length);

// Where a sweep of parameter `name` of the deck at `deck_path` writes: `<deck dir>/<deck stem>_sweep_<name>`.
[[nodiscard]] std::filesystem::path sweep_directory(const std::filesystem::path& deck_path, std::string_view name);
// One value's folder in it: `<name>=<value>`.
[[nodiscard]] std::filesystem::path sweep_point_directory(const std::filesystem::path& sweep_dir, std::string_view name,
                                                          std::string_view value);

// The outputs written under `to` instead: an output below `deck_dir` keeps its path relative to it, one
// elsewhere keeps its file name; an output not written stays not written.
[[nodiscard]] Outputs with_outputs_in(Outputs outputs, const std::filesystem::path& deck_dir, const std::filesystem::path& to);

}  // namespace NiTCAD::deck
