// The bundled example decks (Unit 26, DECISIONS.md VQ11, VA9, WG8): `examples/` beside the executable (or,
// in a build tree, the source's), listed on the start page. Each deck's first comment line is its title,
// the second its description. A user runs a copy: a run writes its outputs beside its deck, and the
// installed examples stay as shipped.
#pragma once

#include <expected>
#include <filesystem>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::app {

struct Example {
    std::filesystem::path path;
    std::string title;        // the first line's comment text, else the file's stem
    std::string description;  // the second line's, else empty
};

// The *.deck files of a directory, by file name. A missing directory has none.
[[nodiscard]] std::vector<Example> list_examples(const std::filesystem::path& directory);

// <executable_dir>/examples when it holds decks, else the fallback (the source tree's examples/).
[[nodiscard]] std::filesystem::path examples_directory(const std::filesystem::path& executable_dir,
                                                       const std::filesystem::path& fallback);

// The example copied into `to` (created if needed) under its own file name, unless a copy is already
// there, which is kept (it may hold the user's edits); the copy's path. Errors (resource_exhausted):
// the directory or the copy cannot be made.
[[nodiscard]] std::expected<std::filesystem::path, base::Error> copy_example(const Example& example,
                                                                             const std::filesystem::path& to);

}  // namespace NiTCAD::app
