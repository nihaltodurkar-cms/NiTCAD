// The deck document of the application (Unit 26, DECISIONS.md VA2, WG3): a deck's text, the file it
// came from and whether it has been modified, and a check that reports what `nitcad check` reports,
// each error at its line.
//
// The deck stays the single source of truth (VA2, VQ13): the application edits this text and runs
// what deck::build makes of it; it keeps no second model of the device.
#pragma once

#include <expected>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/deck/build.hpp"

namespace NiTCAD::app {

// One problem of a deck: its 1-based line (0 when the error names none) and the message as
// `nitcad check` prints it ("<source>:<line>: <what>").
struct Diagnostic {
    int line = 0;
    std::string message;
};

// The line an error names: the number after "<source>:" at the start of its message, else 0.
[[nodiscard]] int error_line(const base::Error& error, std::string_view source);

struct CheckResult {
    std::vector<Diagnostic> diagnostics;  // empty when the deck builds
    std::optional<deck::Plan> plan;      // the deck's plan when it builds
};

// The text parsed and built as the file `path` (its name the source of messages, its directory the
// base of the deck's relative paths); an untitled deck (empty path) is "<untitled>" relative to the
// current directory. Nothing is solved; generating a mesh may take a while.
[[nodiscard]] CheckResult check_deck(std::string_view text, const std::filesystem::path& path);

class DeckDocument {
public:
    // A new deck, not yet saved.
    [[nodiscard]] static DeckDocument untitled(std::string text = {});
    // The file's text. Errors (invalid_input): it cannot be read.
    [[nodiscard]] static std::expected<DeckDocument, base::Error> open(const std::filesystem::path& path);

    [[nodiscard]] const std::filesystem::path& path() const noexcept { return path_; }
    [[nodiscard]] const std::string& text() const noexcept { return text_; }
    [[nodiscard]] bool modified() const noexcept { return modified_; }
    [[nodiscard]] bool is_untitled() const noexcept { return path_.empty(); }
    // The file name, or "Untitled", with " *" when modified.
    [[nodiscard]] std::string title() const;

    // Replaces the text; the document is modified when the text changes.
    void set_text(std::string text);
    // Writes the text to its file (a temporary file next to it, then replaced, so a failed write
    // leaves the old file whole) and clears the modified flag. Errors (invalid_input): untitled;
    // (resource_exhausted): the file cannot be written.
    [[nodiscard]] std::expected<void, base::Error> save();
    // save to a new path, which becomes the document's.
    [[nodiscard]] std::expected<void, base::Error> save_as(const std::filesystem::path& path);

    [[nodiscard]] CheckResult check() const { return check_deck(text_, path_); }

private:
    std::filesystem::path path_;
    std::string text_;
    bool modified_ = false;
};

}  // namespace NiTCAD::app
