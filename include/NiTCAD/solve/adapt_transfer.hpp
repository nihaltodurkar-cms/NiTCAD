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

// How the densities are carried (Unit 17 AR3). quasi_fermi: as above (Newton's initial state).
// exponential_fit: the Scharfetter-Gummel profile instead, for the adjoint estimate's prolongation
// (adapt_adjoint.hpp): with psi linear between two nodes and a constant current, the electrons'
// Slotboom variable n e^(-psi/V_T) follows (1 - e^(-d t)) / (1 - e^(-d)) along the edge
// (d = (psi_b - psi_a) / V_T), which at the midpoint is the weighted mean with weights
// e^(psi_j / 2 V_T); generalized to any interpolant with its weights lambda_j,
//     n = sum lambda_j n_j e^((psi - psi_j) / 2 V_T) / sum lambda_k e^((psi_k - psi) / 2 V_T),
// and for holes with the signs of the exponents reversed. It is exact at an edge's midpoint (so
// the coarse edge's flux is that of its two halves: no residual from the interpolation itself)
// and in equilibrium (n_j = n_i e^(psi_j / V_T) gives n_i e^(psi / V_T)); elsewhere an
// approximation. Quasi-Fermi levels linear across an edge put the drop of the level where the
// edge's average is, not where the density is small, so on a coarse edge through a depletion
// region they make the new halves' currents wrong by orders (measured: 10 times the current at
// 0.6 V on the 41-node diode), which a linearization cannot correct.
enum class TransferMode : unsigned char { quasi_fermi, exponential_fit };

// Errors (invalid_input): the fields or bands not one per old node; the locator not over the old
// mesh (node count); devices with different region lists (count or insulator kind); a new node's
// region with no old node.
[[nodiscard]] std::expected<results::NodeFields, base::Error> transfer_state(
    const device::Device& from, const mesh::PointLocator& from_locator, const results::NodeFields& state,
    const results::BandDiagram& bands, const device::Device& to, TransferReport* report = nullptr,
    TransferMode mode = TransferMode::quasi_fermi);

}  // namespace NiTCAD::solve
