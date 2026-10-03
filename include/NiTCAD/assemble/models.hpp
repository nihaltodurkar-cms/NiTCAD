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
};

}  // namespace NiTCAD::assemble
