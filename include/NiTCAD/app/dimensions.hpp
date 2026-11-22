// A device's named dimensions (Unit 26n, DT1, DECISIONS.md): the deck's [parameters] as labelled rows for the
// Geometry tab's dimensions panel. No Win32.
//
// A row's label is its parameter line's trailing comment (`Lg = 20nm  # Gate length`), else the name; its group the
// text of the last `# group: <name>` comment line above it inside [parameters] (NQ1 (a)). A row's value is written by
// a minimal edit: only that line changes, its key and comment kept as written. The deck stays the only source of
// truth (VQ13).
#pragma once

#include <expected>
#include <string>
#include <string_view>
#include <vector>

namespace NiTCAD::app::dims {

struct Dimension {
    std::string name, label, group;
    std::string expr;        // as written: "20nm", "$Lsd + $Lg", "1e20"
    std::string value_text;  // resolved, for display: "20 nm", "1.5 µm", "1e+20"
    bool length = false;
    double value = 0.0;      // cm for a length
    int line = 0;            // 1-based deck line
};

// The deck's parameters as rows, in deck order; none without [parameters]. Errors: the deck does not parse, or its
// parameters do not resolve (the deck layer's message).
[[nodiscard]] std::expected<std::vector<Dimension>, std::string> read_dimensions(std::string_view deck_text);

// The deck with parameter `name`'s value replaced by `typed` (one line changed). A bare number typed for a length takes
// the row's unit. Errors (the deck unchanged; the message names the parameter): no such parameter, or the deck with the
// value does not parse, resolve or substitute.
[[nodiscard]] std::expected<std::string, std::string> set_dimension(std::string_view deck_text, std::string_view name,
                                                                    std::string_view typed);

// An entry line with its value replaced, its key, spacing and trailing comment kept ("Lg = 20nm  # Gate length", "25nm"
// -> "Lg = 25nm  # Gate length"). The 3D modeller writes its parameters with it too.
[[nodiscard]] std::string with_value(std::string_view line, std::string_view value);

// The expression one mouse-wheel move of `notches` gives: a number times 10^(notches/10); a length plus notches times
// the largest 1-2-5 multiple of the lattice step (0.1 nm) not above a tenth of its value, never below one lattice step,
// in the unit it was written in. A derived row (a sum or a reference) is not stepped: its expression back.
[[nodiscard]] std::string step_dimension(const Dimension& d, int notches);

}  // namespace NiTCAD::app::dims
