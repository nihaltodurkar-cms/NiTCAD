// The adaptation report on disk (Unit 17, slice 6; DECISIONS.md AD8, AR8): results::AdaptReport in
// io's container (magic, format version, content kind adapt_report, length, keyed payload,
// FNV-1a 64), written atomically, doubles by bit pattern, so a loaded report equals the saved one
// field for field (NaN included). Errors as files.hpp's; a report whose per-pass goal arrays are not
// one value per goal and point, or whose certificates name a goal it does not have, is corrupt.
#pragma once

#include <expected>
#include <filesystem>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/adapt.hpp"

namespace NiTCAD::io {

[[nodiscard]] std::expected<void, base::Error> save_adapt_report(const std::filesystem::path& path,
                                                                 const results::AdaptReport& report);
[[nodiscard]] std::expected<results::AdaptReport, base::Error> load_adapt_report(const std::filesystem::path& path);

}  // namespace NiTCAD::io
