// Gate contact boundary terms in scaled variables (ARCHITECTURE.md 6.4; legacy
// moscap.MOSCapacitor.solve_psi row 0, Device2D/3D GateBC, unstructured_dd3d gate_trans).
//
// The oxide under a gate is a lumped capacitor (device/contact.hpp). A gate node i's Poisson row
// gains the displacement flux through it and the fixed oxide charge on its face:
//     F_i += G_i (psi_G,i - psi_i) + S_i
//     G_i     = (eps_ox / eps) (A_i / L_D^(D-1)) / (t_ox / L_D)   (an edge of length t_ox)
//     psi_G,i = (V_G - Phi_i) / V_T,   Phi_i = phi_m - chi_i - Eg_i(T) / 2
//     S_i     = Q_f A_i / (Ns L_D^D)
// with A_i the node's face area on the gate's boundary patch, eps the scaling's reference
// permittivity, and chi_i, Eg_i of the node's material. psi_G is the potential, referenced like psi
// to the intrinsic level, that the electrode's Fermi level imposes on the oxide side: psi is the
// intrinsic level, taken at midgap as in the legacy, so the vacuum level is chi + Eg / 2 above it.
// In 1D with a uniform substrate this is the legacy row
//     kappa (V_G - V_FB - (psi_0 - psi_b) V_T) / V_T,  kappa = eps_ox L_D / (eps t_ox),
//     V_FB = phi_m - (chi + Eg / 2 - psi_b V_T) - q Q_f / C_ox,
// in which psi_b cancels. The legacy subtracted the neutral potential psi_b of the doping it was
// given for V_FB and added that of the node; the two differ only where the doping under the gate is
// not that doping, and there the legacy made the electrode's potential depend on the local doping.
// Here it does not (ARCHITECTURE.md 6.4, OLD / NEW / REASON).
//
// The charge on the gate electrode (per node, in units of q Ns L_D^D) is G_i (psi_G,i - psi_i):
// the displacement flux leaving the electrode. With the fixed charge it balances the
// semiconductor's charge (Gauss).
#pragma once

#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/contact.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::assemble {

// phi_m [eV] of the gate electrode over material m at lattice temperature T. Preconditions as
// physics::band_gap_eV.
[[nodiscard]] double gate_work_function_eV(const device::GateStack& gate,
                                           const physics::Semiconductor& m, double temperature_K);

// Phi = phi_m - chi - Eg(T) / 2 [V]: the electrode's potential at zero bias, referenced to the
// semiconductor's intrinsic (midgap) level.
[[nodiscard]] double gate_midgap_offset_V(const device::GateStack& gate,
                                          const physics::Semiconductor& m, double temperature_K);

struct GateTerm {
    double coupling;      // G_i
    double offset;        // Phi_i / V_T, so psi_G,i = V_G / V_T - offset
    double sheet_charge;  // S_i
};

// The scaled gate term of a node with face area area_cm (cm^(D-1)) in a D-dimensional mesh.
[[nodiscard]] GateTerm gate_term(const device::GateStack& gate, const physics::Semiconductor& m,
                                 double area_cm, int dimension, const Scaling& scaling);

}  // namespace NiTCAD::assemble
