// The input deck's meaning (ARCHITECTURE.md 6.5 and 11.3, Unit 31; DECISIONS.md DK1-DK10): a parsed
// deck (deck/syntax.hpp) turned into a device, solver options and a run plan.
//
// Sections (each at most once unless noted; an unknown section or key is an error, DK3):
//   [device]          temperature = <K>                                   (default 300)
//   [mesh]            type = tensor | generate2d | generate3d | file, and per type:
//     tensor:         x = <axis>, y = <axis>, z = <axis>   (one per dimension, in order)
//                       an axis: lengths separated by commas, or pieces linspace(a, b, n) and
//                       lists joined by '+' (a repeated join point is kept once);
//                     straddle x = <position> <spacing>       (repeatable, any axis; mesh/tensor_grid
//                       straddle_interface: a material interface on the axis)
//     generate2d:     patch <name> = <x y>, <x y>          (repeatable, in order: boundary patches)
//                     max_size, interface_resolution (lengths), min_angle, max_nodes,
//                     box = <x0 y0 x1 y1 h>                (repeatable: a smaller size inside)
//                     region polygons are given in the region sections (polygon = ...)
//     generate3d:     vertex = <x y z>                     (repeatable, in order, numbered from 0)
//                     face <name> = <i j k ...>            (repeatable; a domain face, also its patch)
//                     cuboid = <x0 x1 y0 y1 z0 z1>         (instead of vertices and faces: a box with
//                       faces x_min ... z_max)
//                     patch <face> <name> = <x y z>, ...   (repeatable: a polygon of a tiling of the
//                       face into patches; the face's own name is then not a patch)
//                     plane = <x y z>, <x y z>, <x y z>    (repeatable: a planar interface)
//                     max_size, interface_resolution, radius_edge, max_nodes,
//                     box = <x0 y0 z0 x1 y1 z1 h>          (repeatable)
//                     region seeds are given in the region sections (seed = x y z)
//     file:           path = <file>                        (a generated mesh file carries its node
//                       regions and certificate; a plain mesh file is a one-region mesh)
//   [region <name>]   (repeatable, in order: region k is the k-th section, also for a generated
//                     mesh and a mesh file, whose node regions are indices)
//                     material = silicon | silicon_legacy | germanium | gaas | ingaas | sic_4h
//                                | algaas <x> [<share>] | sio2 | insulator
//                       (silicon: Couderc et al.'s band set, Unit 35; silicon_legacy: the legacy's)
//                     <parameter> = <value>   (a material parameter by its C++ member path, e.g.
//                       eps_r, electron_mobility.mu_max, thermal.conductivity_W_cmK,
//                       passler_gap.Theta_K; DK7; the density-of-states mass fits are not
//                       overridable here)
//                     tensor meshes: box = <lo hi> per dimension (repeatable; a region without boxes
//                       is the background, at most one; later regions override earlier ones, DK8)
//                     generate2d: polygon = <x y>, ... (repeatable)
//                     generate3d with planes: seed = <x y z>, a point inside its slab
//   [doping]          profile = <donor|acceptor> <constant|gaussian|erfc> <peak> [axis=x|y|z]
//                       [center=<l>] [width=<l>] [x=<lo>:<hi>] [y=...] [z=...] [lateral=<l>]
//                       (repeatable; device/doping_profiles.hpp)
//   [contact <name>]  (repeatable) kind = ohmic | gate | electrode; patch = <boundary patch>;
//                     box = <lo hi> per dimension (optional: only the patch nodes inside);
//                     gate: oxide_thickness, oxide_permittivity, electrode = n_poly|p_poly|metal,
//                           work_function, fixed_charge;
//                     electrode: electrode, work_function
//                     The contact's nodes are the patch's semiconductor nodes (ohmic, gate) or its
//                     insulator nodes (electrode).
//   [interface <a> <b>] (repeatable) transport = drift_diffusion | thermionic; fixed_charge;
//                     recombination_velocity_n, recombination_velocity_p;
//                     trap_level = <donor|acceptor> <N_t> <E - E_i> [<sigma_n> <sigma_p>];
//                     trap_band = <donor|acceptor> <D_it> <E_low> <E_high> [<sigma_n> <sigma_p>];
//                     thermal_velocity_n, thermal_velocity_p
//   [thermal_contact <name>] (repeatable) patch; kind = isothermal | resistance; temperature;
//                     resistance
//   [models]          each assemble::PhysicsModels switch by its name = on | off;
//                     btbt_nonlocal = off | kane | direct_wkb; impact_current_resolution = <x>;
//                     hole_mobility = automatic | caughey_thomas | masetti  (Units 34, 35;
//                       automatic, the default: Masetti in silicon regions);
//                     electron_mobility = caughey_thomas | masetti  (Unit 38; opt-in)
//                     ionization_levels = fixed | altermatt  (Unit 37; with incomplete_ionization)
//   [solver]          equations = drift_diffusion | equilibrium_poisson; max_iterations; tol_update;
//                     max_update; line_search = off | residual | natural; line_search.max_backtracks;
//                     line_search.sufficient_decrease; max_cuts; growth; cut_max_iterations;
//                     backend = eigen | pardiso; mkl = <path to mkl_rt>; threads; equilibrate;
//                     max_backward_error; max_refinement_steps; min_pivot_ratio; Ns = <cm^-3>;
//                     weighted_current = on | off  (Unit 32: the weighted terminal current)
//   [equilibrium]     (empty) the analysis is the equilibrium solve, or
//   [sweep]           contact <name> = <series>  (a series: one value held at every point, a comma
//                       list or first : last : step; all series longer than one value have the
//                       same length; a contact not named is held at 0 V);
//                     thermal <name> = <series>  (electrothermal: a thermal contact's temperature);
//                     restart = <checkpoint file>  (continue after the point it holds, DK9)
//   [output]          sweep, device, mesh, checkpoint, table = <file>  (relative to the deck)
// Exactly one of [equilibrium] and [sweep]. Lengths take the suffixes cm, um, nm (DK5); every
// other quantity is in the units of the C++ field it sets (V, K, cm^-3, eV, cm^-2, cm/s, s).
//
// Errors: invalid_input "<source>:<line>: <what>" for anything the deck states wrongly, with the
// line as the context index; the errors of the mesh generators, the materials, the doping profiles
// and Device::create, prefixed with the section they come from.
#pragma once

#include <expected>
#include <filesystem>
#include <optional>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/solve/bias.hpp"

namespace NiTCAD::deck {

enum class Analysis : unsigned char { equilibrium, sweep };

// Files the run writes; an empty path is not written.
struct Outputs {
    std::filesystem::path sweep;       // io::save_sweep
    std::filesystem::path device;      // io::save_device
    std::filesystem::path mesh;        // io::save_generated_mesh, or io::save_mesh
    std::filesystem::path checkpoint;  // the last converged point (equilibrium: the equilibrium)
    std::filesystem::path table;       // a CSV of the bias and terminal currents per point
};

struct Plan {
    device::Device device;
    std::optional<mesh::GeneratedMesh> generated;  // the mesh the deck generated or loaded, if any
    Analysis analysis = Analysis::sweep;
    solve::BiasOptions options;
    std::vector<std::vector<double>> points;  // one bias per contact, device contact order
    std::filesystem::path restart;            // a checkpoint to continue from, or empty
    std::filesystem::path mkl_runtime;        // mkl_rt to load for backend = pardiso, or empty
    Outputs outputs;
};

[[nodiscard]] std::expected<Plan, base::Error> build(const Deck& deck);

}  // namespace NiTCAD::deck
