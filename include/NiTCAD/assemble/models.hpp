// Physical model switches shared by the assemblers (the legacy Models flags of the baseline
// silicon models; the defaults are the legacy ones).
#pragma once

namespace NiTCAD::assemble {

struct PhysicsModels {
    // Caughey-Thomas mobility; false: the lattice mobility (Caughey-Thomas at N = 0).
    bool doping_mobility = true;
    // SRH recombination with Scharfetter lifetimes.
    bool srh = true;
    // Auger recombination (physics::auger_recombination).
    bool auger = true;
    // Slotboom band-gap narrowing: the effective n_ie of every node (statistics, recombination,
    // contacts) and the ln(n_ie) term of the Scharfetter-Gummel driving force.
    bool bgn = true;
    // Canali velocity saturation of each edge's mobility in the edge's own field
    // (physics::canali_mobility; drift-diffusion only). Off by default, as in the legacy.
    bool field_mobility = false;
    // Fermi-Dirac statistics for parabolic bands (physics::fermi_dirac_density and the rest of
    // statistics.hpp): carrier densities, contacts and the charge-neutral guess, the degeneracy
    // term of the Scharfetter-Gummel driving force, and the equilibrium product of SRH and Auger.
    // Off by default (Boltzmann), as the legacy fd flag.
    bool fermi_dirac = false;
};

}  // namespace NiTCAD::assemble
