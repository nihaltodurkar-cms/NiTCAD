// State transfer between meshes of one device (Unit 17, slice 3; DECISIONS.md Unit 17 AD5).
//
// transfer_state carries a solved state from a device onto a remeshed copy of it (the same deck
// with another mesh: the same regions in the same order) as an initial state for Newton. It
// interpolates, at every new node, by the old mesh's interpolant (mesh/locate.hpp) over the old
// nodes of the new node's own region (renormalized; the old interpolant's nodes of another region
// are dropped, so an interface straddled by a cell takes each side's values; with none, the nearest
// old node of the region):
// - the potential psi and the electron and hole quasi-Fermi levels E_Fn, E_Fp (the band diagram's,
//   eV), not the densities, which interpolate badly across decades;
// - the per-node density offsets c_n = ln n - (psi + E_Fn) / V_T and c_p = ln p + (psi + E_Fp) / V_T,
//   from which n = exp(c_n + (psi + E_Fn) / V_T) and p = exp(c_p - (psi + E_Fp) / V_T) at the new
//   node. Under Boltzmann statistics c is ln n_ie plus the band shift, constant inside one material
//   without band-gap narrowing, so a linear psi and linear quasi-Fermi levels transfer exactly; with
//   doping-dependent n_ie (narrowing) or Fermi-Dirac statistics it is an approximation, which only
//   costs Newton iterations: the transferred state is an initial state, never a result;
// - the lattice and carrier temperatures when present.
// Insulator nodes get psi only (densities 0). V_T is the node's lattice temperature's (the device's
// without the electrothermal model).
#pragma once

#include <cstddef>
#include <expected>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/locate.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::solve {

struct TransferReport {
    std::size_t nodes = 0;          // new nodes
    std::size_t outside = 0;        // located outside every old cell (nearest node)
    std::size_t other_region = 0;   // no old interpolant node of their region: nearest of the region
};

// Errors (invalid_input): the fields or bands not one per old node; the locator not over the old
// mesh (node count); devices with different region lists (count or insulator kind); a new node's
// region with no old node.
[[nodiscard]] std::expected<results::NodeFields, base::Error> transfer_state(
    const device::Device& from, const mesh::PointLocator& from_locator, const results::NodeFields& state,
    const results::BandDiagram& bands, const device::Device& to, TransferReport* report = nullptr);

}  // namespace NiTCAD::solve
