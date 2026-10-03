// Private to the assemble layer: the scaled per-node and per-edge data every assembler starts from,
// and the check they share (scaling temperature).
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
#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::assemble::detail {

struct ScaledEdge {
    std::size_t i, j;  // end nodes (mesh Edge::first, Edge::second)
    double geometry;   // coupling area / length / L_D^(D-2)
    double length_cm;  // physical edge length
    double et;         // relative permittivity over the reference one, harmonic mean of the ends
                       // (legacy et)
    bool interface;    // the ends are in materials whose parameters differ (a heterointerface)
};

struct ScaledDevice {
    std::vector<double> volume;          // control volume / L_D^D
    std::vector<double> doping;          // (N_D - N_A) / Ns
    std::vector<double> n_ie;            // n_ie / Ns (effective n_ie with band-gap narrowing)
    std::vector<double> log_dos_n;       // ln(Nc / n_ie), the statistics' g for electrons
    std::vector<double> log_dos_p;       // ln(Nv / n_ie), for holes
    // The band shift s (Unit 15): the node material's intrinsic-level depth below the vacuum
    // level minus that of node 0's material, over V_T. The carriers see psi + s
    // (n = n_ie e^(psi + s) under Boltzmann statistics), so psi is the electrostatic potential,
    // continuous across a heterointerface, and s carries the band offsets; 0 on every node of a
    // single material.
    std::vector<double> band_shift;
    std::vector<std::int32_t> contact;   // index of the node's ohmic contact, or -1
    GateNodes gates;                     // the gate nodes (gate.hpp)
    std::vector<device::ContactKind> kinds;  // per contact, in device.contacts() order
    std::vector<ScaledEdge> edges;       // in mesh edge order
};

// Errors (invalid_input): scaling.temperature_K differs from the device's. An edge between
// materials whose parameters differ is a heterointerface (ScaledEdge::interface): the interface
// lies at the edge's midpoint, each node's control volume is of its own material, and the edge's
// permittivity is the harmonic mean of the ends' (legacy Device1D, _eps_tilde_edge), which is the
// series permittivity of the two half-edges.
[[nodiscard]] std::expected<ScaledDevice, base::Error> make_scaled_device(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models);

// The statistics the assemblers call (physics/statistics.hpp), Boltzmann or Fermi-Dirac, on one
// node's scaled n_ie and its g = ln(N / n_ie) for the carrier's band (unused by Boltzmann). At
// equilibrium eta is psi + s for electrons and -(psi + s) for holes.
[[nodiscard]] inline physics::DensityResult density(bool fermi_dirac, double n_ie, double log_dos,
                                                    double eta) noexcept {
    return fermi_dirac ? physics::fermi_dirac_density(n_ie, log_dos, eta)
                       : physics::boltzmann_density(n_ie, eta);
}

[[nodiscard]] inline physics::NeutralEquilibrium neutral_equilibrium(bool fermi_dirac,
                                                                     double doping, double n_ie,
                                                                     double log_dos_n,
                                                                     double log_dos_p) {
    return fermi_dirac
               ? physics::fermi_dirac_neutral_equilibrium(doping, n_ie, log_dos_n, log_dos_p)
               : physics::boltzmann_neutral_equilibrium(doping, n_ie);
}

// Position of (row, col) in a CSR matrix's values; the entry must exist (NITCAD_EXPECTS).
[[nodiscard]] std::size_t position(const linalg::SparseMatrix& m, std::size_t row,
                                   std::size_t col);

}  // namespace NiTCAD::assemble::detail
