// The Intel MKL runtime, loaded at run time for the mkl_pardiso backend (ARCHITECTURE.md 6.10,
// Unit 18).
//
// NiTCAD is MIT-licensed and does not need MKL to build or run: MKL is an optional external
// dependency under its own licence (the Intel Simplified Software License), which the user installs
// and names by path. Nothing is linked against it and no MKL header is used; the few entry points
// the backend calls are resolved from the loaded mkl_rt with GetProcAddress. Bundling MKL with a
// NiTCAD release would need its licence and notices handled explicitly (ARCHITECTURE.md N6).
//
// The runtime is process-wide and stays loaded. It is set up once, at the first successful load:
// the LP64 interface layer (32-bit indices, as SparseMatrix's) and the Intel OpenMP threading layer
// (libiomp5md.dll must sit next to mkl_rt, as in Intel's redistributable packages). The DLL search
// for MKL's own libraries starts in mkl_rt's directory.
#pragma once

#include <expected>
#include <filesystem>
#include <optional>
#include <string>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::linalg {

struct MklRuntime {
    std::filesystem::path path;  // the mkl_rt loaded
    std::string version;         // MKL_Get_Version_String
};

// Loads mkl_rt from path (for example .../mkl_rt.3.dll), resolves the entry points, sets the layers
// and checks the runtime with a 1 x 1 PARDISO solve. Thread-safe; after a success, later calls
// return the runtime already loaded whatever their path.
// Errors (invalid_input): the file cannot be loaded; an entry point is missing; the layers cannot
// be set; the check solve fails. A failed load leaves nothing loaded.
[[nodiscard]] std::expected<MklRuntime, base::Error> load_mkl_runtime(
    const std::filesystem::path& path);

// The runtime loaded, if any.
[[nodiscard]] std::optional<MklRuntime> loaded_mkl_runtime();

}  // namespace NiTCAD::linalg
