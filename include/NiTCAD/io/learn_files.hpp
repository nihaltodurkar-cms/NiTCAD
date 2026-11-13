// Cross-run learning files (Unit 17b; DECISIONS.md Unit 17 ID8, ID10, Unit 17b LB6, LB7).
//
// The adaptation history (.nithist) in io's container (content kind adapt_history; versioned,
// checksummed, atomic, doubles by bit pattern, so a loaded history equals the saved one bit for bit).
// The manifest is a text file:
//   nitcad-learn-manifest 1
//   schema <n>
//   <digest, 16 hexadecimal digits> <path>
//   ...
// one line per history, sorted by digest; the digest is FNV-1a 64 of the whole history file. A
// manifest's relative paths are relative to the manifest's directory.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <span>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/learn.hpp"

namespace NiTCAD::io {

// Errors as files.hpp's; a history whose passes' features are not rows x feature_count, or whose
// labels are not one per row, is corrupt; a schema other than this build's is refused (invalid_input,
// both versions named).
[[nodiscard]] std::expected<void, base::Error> save_adapt_history(const std::filesystem::path& path,
                                                                  const results::AdaptHistory& history);
[[nodiscard]] std::expected<results::AdaptHistory, base::Error> load_adapt_history(const std::filesystem::path& path);

// FNV-1a 64 of a file's bytes. Errors (invalid_input): the file cannot be read.
[[nodiscard]] std::expected<std::uint64_t, base::Error> file_digest(const std::filesystem::path& path);

// The manifest of these histories: each loaded (so a corrupt or other-schema file is refused here),
// digested, sorted by digest. Errors: a history's load error; two files with one digest (a duplicate).
[[nodiscard]] std::expected<results::Manifest, base::Error> make_manifest(std::span<const std::filesystem::path> histories);
[[nodiscard]] std::expected<void, base::Error> save_manifest(const std::filesystem::path& path,
                                                             const results::Manifest& manifest);
// Errors (invalid_input): not a manifest, another manifest or schema version, a malformed line,
// lines not sorted by digest or a duplicate.
[[nodiscard]] std::expected<results::Manifest, base::Error> load_manifest(const std::filesystem::path& path);

// LB6's per-history cap.
inline constexpr std::size_t max_samples_per_history = 20000;

// The training set of a manifest (LB6): each history loaded from its path (relative paths against
// `base`) and checked against its digest; its labelled edges (NaN labels dropped) in pass then edge
// order, every k-th taken (k = ceil(n / max_per_history)); feature 1 re-indexed to the sorted union
// of the histories' material keys. Errors: a load error; a digest that does not match its file
// (invalid_input, the path named).
[[nodiscard]] std::expected<results::TrainingSet, base::Error> build_training_set(
    const results::Manifest& manifest, const std::filesystem::path& base = {},
    std::size_t max_per_history = max_samples_per_history);

// The learned model (.nitlearn, content kind learned_model; ID10): written with its format and schema
// versions first; a model whose format or schema differs from this build's is refused (invalid_input,
// both versions named: "retrain"); trees whose nodes do not form a tree over valid features, ranges
// or categories not one per feature, are corrupt.
[[nodiscard]] std::expected<void, base::Error> save_learned_model(const std::filesystem::path& path,
                                                                  const results::LearnedModel& model);
[[nodiscard]] std::expected<results::LearnedModel, base::Error> load_learned_model(const std::filesystem::path& path);

}  // namespace NiTCAD::io
