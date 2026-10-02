// Physical model switches shared by the assemblers (the legacy Models flags of the baseline
// silicon models; all on by default, as in the legacy).
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
};

}  // namespace NiTCAD::assemble
