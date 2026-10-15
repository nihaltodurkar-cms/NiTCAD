// Device descriptions on disk (ARCHITECTURE.md 6.3, 6.5 and 11.3, Unit 31).
//
// A device file is the container of files.hpp with content kind `device`: every field of the
// description (device/device.hpp) as a keyed record, the mesh as the mesh file stores it, every
// material parameter by its own key, so a value read back is the value written, bit for bit, and a
// later parameter can be added without a format version change (absent keys take the parameter
// struct's default). Tensor cells are stored as their axes and rebuilt against the mesh; the P4
// certificate is stored as issued.
// A loaded description goes through device::Device::create, so a loaded device is validated like
// any other and a file cannot smuggle in a device the constructor would reject.
//
// Errors: those of files.hpp (invalid_input for a missing, foreign, corrupt or incompatible file;
// resource_exhausted for a failed write), "corrupt or incompatible file: ..." for an unknown
// region kind or trap type, and those of Semiconductor::create, Insulator::create,
// TensorCells::create and Device::create for a description they reject.
#pragma once

#include <expected>
#include <filesystem>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"

namespace NiTCAD::io {

[[nodiscard]] std::expected<void, base::Error> save_device(const std::filesystem::path& path,
                                                           const device::Device& d);
[[nodiscard]] std::expected<device::Device, base::Error> load_device(
    const std::filesystem::path& path);

}  // namespace NiTCAD::io
