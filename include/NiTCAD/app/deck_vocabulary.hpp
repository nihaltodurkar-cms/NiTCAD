// The deck's vocabulary for the editor (Unit 26b, DECISIONS.md XD2, XD3): every section kind, the keys
// each takes, the choices of the values that have them and a line of help, as deck/build.hpp's
// reference lists them (a test reads that comment and asks deck::build about each key, so the two
// cannot drift). From it: completion at the caret and "did you mean" for a misspelled name.
//
// A material parameter of a region by its member path (eps_r, electron_mobility.mu_max, ...) is
// open-ended and not listed.
#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/document.hpp"

namespace NiTCAD::app {

struct KeyInfo {
    std::string_view key;       // the key's first word: "max_size", "contact", "patch"
    int names = 0;              // words after it naming something: "contact <name>" 1, "patch <face> <name>" 2
    bool repeatable = false;    // may appear more than once (named keys always may, with other names)
    std::string_view value;     // the default, or the value's form
    std::string_view help;
    std::span<const std::string_view> choices;   // the value's first word is one of these, if any
    std::span<const std::string_view> choices2;  // and its second word (a doping profile's shape)
    std::string_view contexts;  // the mesh types or contact kinds it belongs to, blank-separated ("" all)
};

struct SectionInfo {
    std::string_view kind;
    int names = 0;               // [region <name>] 1, [interface <a> <b>] 2
    bool repeatable = false;
    std::string_view help;
    std::span<const KeyInfo> keys;
    std::string_view context_key;  // the key whose value selects the keys: "type" ([mesh]), "kind" ([contact])
};

// Every section kind, in the order of deck/build.hpp's reference.
[[nodiscard]] std::span<const SectionInfo> deck_vocabulary() noexcept;
[[nodiscard]] const SectionInfo* find_section_info(std::string_view kind) noexcept;
// The keys of a section that apply in a context (the mesh's type, the contact's kind); all of them
// when the context is empty or unknown.
[[nodiscard]] std::vector<const KeyInfo*> section_keys(const SectionInfo& s, std::string_view context);
// The key an entry's key text starts with, in the context.
[[nodiscard]] const KeyInfo* find_key_info(const SectionInfo& s, std::string_view key, std::string_view context);

// Levenshtein distance (insertions, deletions, substitutions).
[[nodiscard]] std::size_t edit_distance(std::string_view a, std::string_view b);
// The name nearest the word, if within min(2, max(1, size / 3)) edits (the first in order on a tie);
// nothing for an exact match.
[[nodiscard]] std::optional<std::string> nearest(std::string_view word, std::span<const std::string_view> names);

// "did you mean 'x'?" for the diagnostic's line in the deck text: an unknown section kind, a key the
// section does not take, or a value word not among its key's choices, when a name is near enough.
[[nodiscard]] std::optional<std::string> suggestion(std::string_view text, const Diagnostic& d);

// Completion at a caret (line 1-based, column a byte offset into the line's UTF-8): the section
// kinds after '[', a section's keys at an entry's start (those present not offered again unless they
// repeat), a value's choices after '='. Nothing in a comment, a header's name, a continued value or
// before the first section.
struct Candidate {
    std::string label;   // as listed
    std::string insert;  // the text that replaces [from, to)
    std::string help;
};
struct Completion {
    int line = 0;
    std::size_t from = 0, to = 0;  // byte columns on the line: the word at the caret
    std::vector<Candidate> items;  // those starting with the typed part of the word, in vocabulary order
};
[[nodiscard]] std::optional<Completion> complete(std::string_view text, int line, std::size_t column);

}  // namespace NiTCAD::app
