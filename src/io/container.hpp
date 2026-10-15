// The container and the shared records of the io layer (private; files.hpp, device_file.hpp).
#pragma once

#include <cstddef>
#include <expected>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "record.hpp"

// Unwraps an expected into `var`, or returns its error from the enclosing function (every read is
// a chain of fallible field reads).
#define NITCAD_IO_TRY(var, expr)                                                  \
    auto var##_expected = (expr);                                                 \
    if (!var##_expected) return std::unexpected(std::move(var##_expected.error())); \
    auto var = std::move(*var##_expected)

namespace NiTCAD::io::detail {

base::Error invalid(std::string message);

// The container around a payload, written atomically (files.hpp).
std::expected<void, base::Error> write_file(const std::filesystem::path& path, Content kind,
                                            const Writer& payload);
// The payload of a container of this kind, checked (files.hpp).
std::expected<std::vector<std::byte>, base::Error> read_file(const std::filesystem::path& path,
                                                             Content kind);

Writer write_mesh(const mesh::Mesh& m);
// Rebuilt through mesh::Mesh::from_parts.
std::expected<mesh::Mesh, base::Error> read_mesh(const Reader& r);

Writer write_certificate(const mesh::InterfaceCertificate& c);
// The stored certificate, unchecked (the caller checks it against its mesh and regions).
std::expected<mesh::InterfaceCertificate, base::Error> read_certificate(const Reader& r);

}  // namespace NiTCAD::io::detail
