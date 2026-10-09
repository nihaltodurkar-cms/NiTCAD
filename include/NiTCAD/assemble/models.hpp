// Physical model switches shared by the assemblers (the legacy Models flags of the baseline
// silicon models; the defaults are the legacy ones).
#pragma once

#include <cstdint>

namespace NiTCAD::assemble {

enum class NonlocalTunnelling : std::uint8_t { off, kane, direct_wkb };

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
    // Radiative recombination (physics::radiative_recombination) with each material's coefficient
    // (drift-diffusion only). On by default: silicon's coefficient is 0, so silicon results do not
    // change (Unit 15).
    bool radiative = true;
    // Incomplete dopant ionization (physics/ionization.hpp) with each material's levels: the
    // ionized doping in Poisson, the contacts and the neutral guess. Off by default, as the legacy
    // incomplete_ion flag. Thermionic emission is not a switch here: it is chosen per interface
    // (device::Interface; Unit 15).
    bool incomplete_ionization = false;
    // Impact ionization (physics/impact_ionization.hpp; drift-diffusion only, Unit 19): generation
    // (alpha_n |J_n| + alpha_p |J_p|) / q at each node with each material's coefficients, the field
    // taken along each carrier's current (drift_diffusion.hpp). Off by default, as the legacy.
    bool impact_ionization = false;
    // The resolution of a current for impact ionization, relative to the opposing flux terms it is
    // the difference of (drift_diffusion.hpp): a current at or below it does not ionize. A
    // numerical option, not physics, off (0) by default: 1e-12 (50 times the measured rounding of
    // the spill-over layer at a junction) removes the generation of that layer's unresolved
    // currents and the corrector cycles they cause near breakdown, but stalls the trace of an
    // open-base transistor at its start (ARCHITECTURE.md 6.2, Unit 19).
    double impact_current_resolution = 0.0;
    // Local band-to-band tunnelling (physics/band_to_band.hpp; drift-diffusion only, Unit 20):
    // Kane generation A F^2 exp(-B / F) at each semiconductor node off the ohmic contacts, F the
    // magnitude of the node's reconstructed field, with each material's (A, B). Pure generation
    // (no Hurkx D factor). Off by default, as the legacy btbt flag.
    bool btbt_local = false;
    // Nonlocal band-to-band tunnelling along traced paths (assemble/tunnel_paths.hpp; Unit 20),
    // off by default: `kane` the calibrated rate A F^2 exp(-B / F) at each path's mean field (the
    // silicon model), `direct_wkb` the direct-gap WKB rate of Esseni et al. (2017) for materials
    // with a direct gap and cited tunnelling masses (never silicon). Needs the device's tensor
    // cells; excludes btbt_local (each would count the same tunnelling).
    NonlocalTunnelling btbt_nonlocal = NonlocalTunnelling::off;
    // Electrothermal coupling (Unit 23; DECISIONS.md T1-T14; drift_diffusion.hpp): the lattice
    // temperature is a fourth unknown per node, with the heat equation, the thermopower in the
    // carrier fluxes and the temperature dependence of the band, mobility and contact models; heat
    // leaves through the device's thermal contacts. Drift-diffusion only. Off by default; off,
    // every result is that of the isothermal device at its temperature.
    bool electrothermal = false;
};

}  // namespace NiTCAD::assemble
