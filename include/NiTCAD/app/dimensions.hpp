// A device's named dimensions (Unit 26n, DT1, DECISIONS.md): the deck's [parameters] as labelled rows for the
// Geometry tab's dimensions panel. No Win32.
//
// A row's label is its parameter line's trailing comment (`Lg = 20nm  # Gate length`), else the name; its group the
// text of the last `# group: <name>` comment line above it inside [parameters] (NQ1 (a)). A row's value is written by
// a minimal edit: only that line changes, its key and comment kept as written. The deck stays the only source of
// truth (VQ13).
#pragma once

#include <array>
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

// Where a length dimension shows as a line (DT2): the deck substituted at p and at p + δ (δ a hundredth of |p|, at
// least one lattice step); among the extents of region and contact boxes (and a generate3d domain), [doping] windows
// (`x=a:b`) and the [mesh] axes' linspace pieces, a candidate changes its length by ±δ with exactly one end moving.
// One span per axis, the shortest candidate (ties: boxes, then doping windows, then mesh pieces, each in deck order).
// None for a number, an unknown name, or a deck that does not substitute; polygons and prisms give none (slice 0).
struct Span {
    int axis = 0;               // 0 x, 1 y, 2 z
    double lo = 0.0, hi = 0.0;  // cm along the axis, at the deck's value
    std::array<double, 3> at{};  // cm: the lo end; the other two coordinates place the line (a box's low corner)
    std::string source;         // "region oxide", "contact gate", "mesh x", "doping 2" (1-based profile index)
    int grows = 1;              // +1: the moving end goes toward +axis as the value grows; -1 toward -axis
};
[[nodiscard]] std::vector<Span> dimension_spans(std::string_view deck_text, std::string_view name);

// What the editors draw (DT2): a line per length row and span, from `a` (the span's lo end) to `b` (its hi end), cm,
// labelled "<label> = <value>"; dragging it toward `grows` along `axis` makes the value larger.
struct Mark {
    std::string name, label;
    int axis = 0;
    std::array<double, 3> a{}, b{};
    int grows = 1;
};
[[nodiscard]] std::vector<Mark> dimension_marks(std::string_view deck_text);

// A length row's value `cm` as a literal rounded to the lattice (0.1 nm) in the row's unit (a derived row: nm below
// 1 µm, else um); a positive row stays at least one lattice step.
[[nodiscard]] std::string length_literal(const Dimension& d, double cm);

// The deck's text with every reference in [mesh], [region], [contact] and [doping] written out (deck::substitute), each
// line rewritten in place, [parameters] and every other line kept: what the 2D editor draws a parametric deck from
// (view only). Errors: the deck does not parse or a reference does not resolve.
[[nodiscard]] std::expected<std::string, std::string> substituted_text(std::string_view deck_text);

}  // namespace NiTCAD::app::dims
