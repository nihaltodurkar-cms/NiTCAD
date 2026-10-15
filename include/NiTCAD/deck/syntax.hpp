// The input deck's syntax (ARCHITECTURE.md 6.5 and 11.3, Unit 31; DECISIONS.md DK1-DK10).
//
// A deck is a text file of sections, each a header line and key = value entries:
//
//     # a comment, to the end of the line (outside double quotes)
//     [region silicon]            <- a section: its kind ("region") and name ("silicon")
//     material = silicon          <- an entry: key "material", value "silicon"
//     polygon = 0 0, 1um 0, \     <- a trailing backslash continues the value on the next line
//               1um 1um, 0 1um
//
// - A header is `[kind]` or `[kind name]`: the kind one identifier, the name the rest of the
//   bracket with runs of blanks collapsed to one (an interface's name is its two regions).
// - A key is one or more identifiers separated by blanks (`contact anode`), collapsed likewise;
//   an identifier is letters, digits, '_', '.' and '-', starting with a letter.
// - A value is the text after the first '=' with the comment removed, trimmed; it may not be
//   empty. Double quotes keep '#' and blanks as written (paths); the quotes stay in the value
//   for the reader of that key to remove (unquote).
// - Keys may repeat; their order is kept. What each section and key means is not decided here
//   (deck/build.hpp): this layer only splits the text, so a syntax error and a meaning error are
//   reported apart, both with the line.
// Errors: invalid_input "<source>:<line>: <what>" with the line as the context index: an entry
// before the first header, a malformed header or key, a missing '=', an empty value, an
// unterminated quote, a continuation at the end of the file, a NUL byte.
#pragma once

#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::deck {

struct Entry {
    std::string key;
    std::string value;
    int line = 0;  // 1-based, the line the entry starts on
};

struct Section {
    std::string kind;
    std::string name;  // empty when the header has none
    int line = 0;
    std::vector<Entry> entries;
};

struct Deck {
    std::string source;              // the file name used in messages ("<text>" when parsed from text)
    std::filesystem::path directory;  // relative paths in the deck are relative to this
    std::vector<Section> sections;
};

[[nodiscard]] std::expected<Deck, base::Error> parse(std::string_view text,
                                                     std::string source = "<text>",
                                                     std::filesystem::path directory = {});

// The file's text parsed, with its name as the source and its directory as the base for relative
// paths. Errors: invalid_input if it cannot be read; those of parse.
[[nodiscard]] std::expected<Deck, base::Error> read(const std::filesystem::path& path);

// A value without its surrounding double quotes, if it has them.
[[nodiscard]] std::string unquote(std::string_view value);

}  // namespace NiTCAD::deck
