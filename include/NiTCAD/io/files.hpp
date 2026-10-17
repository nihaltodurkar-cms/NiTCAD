// Files: meshes, sweeps and checkpoints on disk (ARCHITECTURE.md 6.5, 6.6 and 11.3, Unit 29; the
// io layer approved by the owner on 2026-10-09).
//
// One container format for every content:
//   magic "NiTCADio" (8 bytes) | format version (u32) | content kind (u32) |
//   payload length (u64) | payload | FNV-1a 64 of the payload (u64)
// all little-endian. The payload is a keyed record (src/io/record.hpp): unknown keys are skipped and
// absent ones take defaults, so files of a later minor revision load and earlier files keep
// loading; a different format version is refused. Doubles are stored by bit pattern, so every
// value read back equals the value written, bit for bit.
// Writes are atomic: the file is written beside the target as "<name>.tmp", flushed, closed and
// renamed over the target, so a crash or a full disk never leaves a half-written file under the
// target's name (an existing target is kept until the rename).
// The in-memory types know nothing of this format (R2): io depends on them, never the reverse.
//
// Errors (base::Error):
// - invalid_input: a file that cannot be opened; not a NiTCAD file; another format version or
//   content kind; a checksum mismatch, truncation or a malformed record ("corrupt or incompatible
//   file: ..."); a mesh the mesh layer rejects (its own errors, degenerate_mesh included); a
//   checkpoint of another mesh;
// - resource_exhausted: the file or its temporary could not be written, flushed or renamed.
#pragma once

#include <cstdint>
#include <expected>
#include <filesystem>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::io {

inline constexpr std::uint32_t format_version = 1;

enum class Content : std::uint32_t { mesh = 1, sweep = 2, checkpoint = 3, device = 4 };

// A 64-bit FNV-1a digest of the mesh graph (dimension, points, volumes, edges, boundary patches;
// doubles by bit pattern). Equal meshes give equal digests; it identifies a mesh on any machine,
// since it hashes the stored values, not computed ones.
[[nodiscard]] std::uint64_t mesh_digest(const mesh::Mesh& m);

[[nodiscard]] std::expected<void, base::Error> save_mesh(const std::filesystem::path& path,
                                                         const mesh::Mesh& m);
// Rebuilt through mesh::Mesh::from_parts, so a loaded mesh is validated like any other.
[[nodiscard]] std::expected<mesh::Mesh, base::Error> load_mesh(const std::filesystem::path& path);

// A generated mesh (Unit 16; DECISIONS.md M14): the mesh file with optional records for its
// simplices, node regions, generator report and P4 certificate, so a reloaded mesh can be rendered,
// refined or given to a multi-region device without regenerating it. load_mesh reads the same
// file as a plain mesh (unknown records are skipped). load_generated_mesh refuses (invalid_input)
// a file without those records, simplices whose dual is not the stored mesh, and a certificate
// that does not certify the stored mesh and node regions.
[[nodiscard]] std::expected<void, base::Error> save_generated_mesh(const std::filesystem::path& path,
                                                                   const mesh::GeneratedMesh& m);
[[nodiscard]] std::expected<mesh::GeneratedMesh, base::Error> load_generated_mesh(
    const std::filesystem::path& path);

// Every field of a sweep: the run record, every point (fields, bands, currents, resolutions,
// charges, edge currents, convergence, tunnel paths, thermal quantities, sub-steps), the stopping
// error and the unfinished history.
[[nodiscard]] std::expected<void, base::Error> save_sweep(const std::filesystem::path& path,
                                                          const results::Sweep& sweep);
[[nodiscard]] std::expected<results::Sweep, base::Error> load_sweep(
    const std::filesystem::path& path);

// A restart point: the state of one converged bias point and what it belongs to. A sweep restarts
// from it with solve::sweep_bias(device, remaining_points, options, &c.fields) and
// options.step_control.initial_bias_V = c.bias_V.
struct Checkpoint {
    std::uint64_t mesh_digest = 0;     // of the mesh the state lives on
    std::uint64_t run_identity = 0;    // RunRecord::input_identity of the run that produced it
    std::vector<double> bias_V;        // per contact
    std::vector<double> thermal_bias_K;  // per thermal contact; empty without the model
    results::NodeFields fields;
    // Of the run that wrote it (Unit 42); empty for a checkpoint written before Unit 42.
    results::Provenance provenance;
};

[[nodiscard]] std::expected<void, base::Error> save_checkpoint(const std::filesystem::path& path,
                                                               const mesh::Mesh& m,
                                                               const results::RunRecord& run,
                                                               const results::BiasPoint& point);
// Refused (invalid_input) unless the checkpoint was written for a mesh with this mesh's digest and
// its per-node fields have one value per node of it.
[[nodiscard]] std::expected<Checkpoint, base::Error> load_checkpoint(
    const std::filesystem::path& path, const mesh::Mesh& m);

}  // namespace NiTCAD::io
