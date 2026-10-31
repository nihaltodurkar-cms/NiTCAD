// The input deck's meaning (ARCHITECTURE.md 6.5 and 11.3, Unit 31; DECISIONS.md DK1-DK10): a parsed
// deck (deck/syntax.hpp) turned into a device, solver options and a run plan.
//
// Sections (each at most once unless noted; an unknown section or key is an error, DK3):
//   [device]          temperature = <K>                                   (default 300)
//   [mesh]            type = tensor | generate2d | generate3d | extrude3d | file, and per type:
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
//     extrude3d:      axis = x | y | z; stations = <lengths>      (increasing; Unit 53, FD13)
//                     cell <name> = <u v>, <u v>, ...      (repeatable: a conforming tiling of the cross-section;
//                       u, v follow the axis cyclically: x gives (y, z), y (z, x), z (x, y))
//                     patch <name> = <x0 y0 z0 x1 y1 z1>   (repeatable: the whole boundary facets inside
//                       the closed box, to rounding, form the patch; a box selecting none is an error; the rest are "boundary")
//                     max_size, interface_resolution, radius_edge, max_nodes, box as generate3d
//                     sweep = on | off                     (default off; Unit 54, UD10) on: the swept mesh
//                       (mesh::generate_swept): the cross-section meshed in 2D and swept along the axis, with
//                       lines = <axis>                     (required: as a tensor axis; the first and last
//                         the end stations exactly)
//                       straddle = <station> <width>       (repeatable: lines at station +- width/2, required at
//                         every interior station where a cell changes region)
//                       layer = <u v>, <u v> along=<axis> lines=<axis> [straddle=<length>]
//                                                          (repeatable: a structured layer on the base from the
//                         first point to the second, its lines signed, positive to the base's left; straddle the
//                         width about each interface it covers, default interface_resolution; two layers whose
//                         bases share an end point form a corner block)
//                       max_size, interface_resolution, min_angle as generate2d; max_nodes over the 3D nodes;
//                       box = <u0 v0 u1 v1 h> in the cross-section; radius_edge refused. Off: lines, straddle,
//                       layer and min_angle refused.
//                     region extrusions are given in the region sections (extrude = <cells> <t_lo>:<t_hi>,
//                       the ends being stations, to rounding; a cell and section assigned twice is an error, an
//                       unassigned one is void)
//     file:           path = <file>                        (a generated mesh file carries its node
//                       regions and certificate; a plain mesh file is a one-region mesh)
//   [region <name>]   (repeatable, in order: region k is the k-th section, also for a generated
//                     mesh and a mesh file, whose node regions are indices)
//                     material = silicon | silicon_legacy | germanium | gaas | ingaas | sic_4h
//                                | algaas <x> [<share>] | sio2 | insulator
//                       (silicon: Couderc et al.'s band set, Unit 35; silicon_legacy: the legacy's)
//                     <parameter> = <value>   (a material parameter by its C++ member path, e.g.
//                       eps_r, electron_mobility.mu_max, thermal.conductivity_W_cmK,
//                       passler_gap.Theta_K, energy_relaxation.tau_n_s (Unit 52); DK7; the
//                       density-of-states mass fits are not overridable here)
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
//                     surface_mobility = off | lombardi  (Unit 46; opt-in)
//                     quantum_confinement = off | mlda  (Unit 47; opt-in); mlda_mass_n,
//                       mlda_mass_p = <m / m0>  (defaults 0.98, 0.49)
//                     energy_balance = on | off  (Unit 52; opt-in: carrier temperatures)
//                     ionization_levels = fixed | altermatt  (Unit 37; with incomplete_ionization)
//   [solver]          equations = drift_diffusion | equilibrium_poisson; max_iterations; tol_update;
//                     max_update; line_search = off | residual | natural; line_search.max_backtracks;
//                     line_search.sufficient_decrease; max_cuts; growth; cut_max_iterations;
//                     backend = eigen | pardiso; mkl = <path to mkl_rt>; threads; equilibrate;
//                     max_backward_error; max_refinement_steps; min_pivot_ratio; Ns = <cm^-3>;
//                     weighted_current = on | off  (Unit 32: the weighted terminal current)
//                     rates = on | off  (Unit 45: store the generation-recombination rates; on)
//                     minority_floor = on | off  (Unit 50: densities against eps of the node's
//                       majority in the update size; off)
//   [equilibrium]     (empty) the analysis is the equilibrium solve, or
//   [sweep]           contact <name> = <series>  (a series: one value held at every point, a comma
//                       list or first : last : step; all series longer than one value have the
//                       same length; a contact not named is held at 0 V);
//                     thermal <name> = <series>  (electrothermal: a thermal contact's temperature);
//                     restart = <checkpoint file>  (continue after the point it holds, DK9)
//   [transient]       (Unit 43; solve/transient.hpp) contact <name> = <waveform>  (a contact not
//                       named is held at 0 V): constant(v), step(v0, v1, t),
//                       ramp(v0, v1, t0, t1), pulse(base, pulse, t_start, width[, rise, fall]),
//                       pwl(t v, t v, ...), sine(offset, amplitude, frequency[, delay, phase]);
//                     t_end (required), dt_initial, dt_min, dt_max, integrator = bdf2 |
//                     backward_euler, adaptive = on | off, rtol, density_ref,
//                     output_times = <series>, fields_every_step = on | off;
//                     thermal <name> = <K>  (electrothermal: held for the whole run)
//   [small_signal]    (Unit 43; solve/small_signal.hpp) contact <name> = <series> and
//                     thermal <name> = <series>, the operating points as in [sweep];
//                     frequencies = <series> (required), field_frequencies = <series>
//   [trace]           (Unit 43; solve/trace.hpp) contact = <name> and end = <V> (required);
//                     start <name> = <V>  (each contact's bias at the start, default 0);
//                     step, max_step, min_step, current_limit, max_points, density_floor;
//                     thermal <name> = <K>
//   [output]          sweep, device, mesh, checkpoint, table = <file>  (relative to the deck);
//                     transient, small_signal = <file>  (Unit 43). Each analysis takes only what
//                     it writes: equilibrium device, mesh, checkpoint; sweep and trace sweep,
//                     device, mesh, checkpoint, table; transient transient, device, mesh, table
//                     (no checkpoint: its end state is not steady); small_signal small_signal,
//                     device, mesh, checkpoint (the last operating point), table.
// Exactly one analysis section: [equilibrium], [sweep], [transient], [small_signal] or [trace].
// A series is a comma list, first : last : step (DK6) or logspace(first, last, count) (Unit 43:
// count values evenly spaced in log10, first and last exact, both positive). Lengths take the
// suffixes cm, um, nm (DK5); every other quantity is in the units of the C++ field it sets (V, K,
// cm^-3, eV, cm^-2, cm/s, s, Hz).
//
// Errors: invalid_input "<source>:<line>: <what>" for anything the deck states wrongly, with the
// line as the context index; the errors of the mesh generators, the materials, the doping profiles
// and Device::create, prefixed with the section they come from.
#pragma once

#include <expected>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "NiTCAD/solve/waveform.hpp"

namespace NiTCAD::deck {

enum class Analysis : unsigned char { equilibrium, sweep, transient, small_signal, trace };

// Files the run writes; an empty path is not written.
struct Outputs {
    std::filesystem::path sweep;         // io::save_sweep (a sweep or a trace)
    std::filesystem::path device;        // io::save_device
    std::filesystem::path mesh;          // io::save_generated_mesh, or io::save_mesh
    std::filesystem::path checkpoint;    // the last converged point (equilibrium: the equilibrium)
    std::filesystem::path table;         // a CSV of the analysis's terminal quantities
    std::filesystem::path transient;     // io::save_transient (Unit 43)
    std::filesystem::path small_signal;  // io::save_small_signal (Unit 43)
};

struct Plan {
    device::Device device;
    std::optional<mesh::GeneratedMesh> generated;  // the mesh the deck generated or loaded, if any
    std::string mesh_type;                         // the [mesh] type, as written; extrude3d for a file mesh with a generation digest (Unit 53, FD14)
    Analysis analysis = Analysis::sweep;
    // Models, solver and, for a sweep or a small-signal run, the thermal contacts' temperature per
    // point (for a transient or a trace, one list held for the run).
    solve::BiasOptions options;
    // Sweep and small-signal: the points, one bias per contact, device contact order.
    std::vector<std::vector<double>> points;
    std::filesystem::path restart;            // a checkpoint to continue from, or empty
    std::filesystem::path mkl_runtime;        // mkl_rt to load for backend = pardiso, or empty
    Outputs outputs;
    // The analysis's own options (Unit 43). Their `steady` is left default: the run uses `options`.
    solve::TransientOptions transient;
    std::vector<solve::Waveform> waveforms;  // one per contact, device contact order
    solve::SmallSignalOptions small_signal;
    solve::TraceOptions trace;
};

[[nodiscard]] std::expected<Plan, base::Error> build(const Deck& deck);

}  // namespace NiTCAD::deck
