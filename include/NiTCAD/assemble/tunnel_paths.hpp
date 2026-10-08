// Nonlocal band-to-band tunnelling paths (Unit 20; legacy nonlocal_path.py and paths.cpp;
// ARCHITECTURE.md 6.2): frozen geometry traced through a tensor mesh's cells, evaluated live at
// each state.
//
// A path is a polyline of samples, each with the multilinear stencil of the cell it lies in, so a
// sample's potential and band edges are linear in the nodal potentials. It is traced once at a
// state and then frozen; everything that depends on the potential is evaluated live with exact
// partials (the geometry's own dependence on the potential is not differentiated: a relocation
// re-traces it, drift_diffusion.hpp).
//
// Tracing (trace_tunnel_paths), in the scaled potential psi (units of V_T) and cm:
// - Only cells whose corners are all semiconductor nodes are traced; a face to any other cell or
//   the device boundary is a Neumann boundary for the path: the direction's outward component is
//   dropped there.
// - A node i starts a path if and only if (1) it is a semiconductor node and not an ohmic-contact
//   node, (2) its material has the model's parameters (max_length_cm[i] > 0), (3) some
//   semiconductor node j lies within max_length_cm[i] plus the largest cell diagonal of i with
//   end_band[j] <= valence[i] (a superset of the reachable crossings: a crossing within the length
//   is a weighted mean of its stencil's end_band, each node of which lies within a diagonal of it),
//   and (4) the nodal gradient of psi at i is not zero.
// - The path runs along +grad(psi) (electrons tunnel towards higher potential), the direction at
//   a point the multilinear interpolation of the nodal gradients (second-order central differences,
//   first order at the mesh's ends). Each step runs straight to the nearer of the cell's next face
//   (landing exactly on it) and half the cell's smallest width, so each segment lies in one cell.
// - With Delta(s) = valence[start] - end_band(s), the crossing is the first segment on which Delta
//   reaches 0. A path stops when Delta reaches half the start's gap beyond it (room for the
//   crossing to move during Newton; the legacy margin), before a sample with weight on a contact
//   node, when its direction vanishes, after 64 steps without a new maximum of psi, or after
//   100000 steps; one that has not crossed within max_length_cm is dropped (it would carry no
//   generation: the length is where the rate falls below exp(-690) of its scale).
// - `extend` (the direct-gap WKB rate, which needs the band extrema around the path): the path
//   continues past the crossing, and a second polyline runs from the start along -grad(psi), each
//   until the field falls below 1e3 V/cm (the end of the band bending) or the length reaches ten
//   times max_length_cm.
#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"

namespace NiTCAD::assemble {

struct TunnelSample {
    mesh::Point position;            // [cm]
    std::array<std::size_t, 8> nodes{};
    std::array<double, 8> weights{};
    int count = 0;                   // 2^D
};

struct TunnelPath {
    std::size_t start = 0;
    std::vector<TunnelSample> forward;  // forward[0] is the start node
    std::vector<double> length_cm;      // forward segment k joins samples k and k + 1
    // The segment on which the crossing lay when the path was traced: the live evaluation keeps
    // it, extrapolating along it while the crossing moves (relocation then moves it).
    std::size_t crossing = 0;
    std::vector<TunnelSample> backward;  // `extend` only: from the start along -grad(psi)
};

struct TunnelPaths {
    std::vector<TunnelPath> paths;  // in increasing start node
};

// Whether two paths from the same start have the same geometry: crossing segments at most one
// apart, and the samples up to the end of the later crossing segment at the same stencil nodes and
// within tolerance_cm of each other. The samples beyond, and the backward polyline, are not
// compared (they enter the rate only through the band extrema). Which paths must agree is the
// assembler's (DriftDiffusion::paths_agree: those with a rate that matters at the state).
[[nodiscard]] bool same_geometry(const TunnelPath& a, const TunnelPath& b,
                                 double tolerance_cm = 1e-9);

// What tracing reads, per node (scaled potentials and band edges in units of V_T).
struct TunnelTraceInput {
    const mesh::TensorCells* cells;
    std::span<const double> psi;
    std::span<const double> valence;        // E_v
    std::span<const double> end_band;       // the band the electron ends in (see the models)
    std::span<const char> semiconductor;    // 1 on semiconductor nodes
    std::span<const char> contact;          // 1 on ohmic-contact nodes
    std::span<const double> max_length_cm;  // 0 where no path may start
    double V_T = 0.0;                       // [V]
    bool extend = false;
};

// Preconditions (NITCAD_EXPECTS): cells set; every span has the cells' node count; V_T positive.
[[nodiscard]] TunnelPaths trace_tunnel_paths(const TunnelTraceInput& input);

// One path at a state, in physical units (DriftDiffusion::path_states).
struct TunnelPathState {
    std::size_t start;
    double length_cm;       // start to the crossing, along the path
    double field_V_per_cm;  // the potential drop over that length divided by it
    double rate_cm3_s;      // the path's generation rate (per volume of its start node)
    double fraction;        // t: where the crossing lies along the frozen segment
    bool reached;           // the crossing lies on the frozen segment (t in [0, 1] to 1e-6)
};

}  // namespace NiTCAD::assemble
