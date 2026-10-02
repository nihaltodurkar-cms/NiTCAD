// Private to the assemble layer: the scaled per-node and per-edge data every assembler starts from,
// and the checks they share (scaling temperature, heterojunctions).
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <vector>

#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace NiTCAD::assemble::detail {

struct ScaledEdge {
    std::size_t i, j;  // end nodes (mesh Edge::first, Edge::second)
    double geometry;   // coupling area / length / L_D^(D-2)
    double et;         // relative permittivity over the reference one (legacy et)
};

struct ScaledDevice {
    std::vector<double> volume;          // control volume / L_D^D
    std::vector<double> doping;          // (N_D - N_A) / Ns
    std::vector<double> n_ie;            // n_ie / Ns (effective n_ie with band-gap narrowing)
    std::vector<std::int32_t> contact;   // index of the node's ohmic contact, or -1
    GateNodes gates;                     // the gate nodes (gate.hpp)
    std::vector<device::ContactKind> kinds;  // per contact, in device.contacts() order
    std::vector<ScaledEdge> edges;       // in mesh edge order
};

// Errors (invalid_input): scaling.temperature_K differs from the device's; a heterojunction (an
// edge between regions whose material parameters differ; band offsets and permittivity steps are
// deferred), with the edge as the context index.
[[nodiscard]] std::expected<ScaledDevice, base::Error> make_scaled_device(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models);

// Position of (row, col) in a CSR matrix's values; the entry must exist (NITCAD_EXPECTS).
[[nodiscard]] std::size_t position(const linalg::SparseMatrix& m, std::size_t row,
                                   std::size_t col);

}  // namespace NiTCAD::assemble::detail
