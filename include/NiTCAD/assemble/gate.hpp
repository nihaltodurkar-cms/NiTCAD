// Gate contact boundary terms in scaled variables (ARCHITECTURE.md 6.4; legacy
// moscap.MOSCapacitor.solve_psi row 0, Device2D/3D GateBC, unstructured_dd3d gate_trans).
//
// The oxide under a gate is a lumped capacitor (device/contact.hpp). A gate node i's Poisson row
// gains the displacement flux through it and the fixed oxide charge on its face:
//     F_i += G_i (psi_G,i - psi_i) + S_i
//     G_i     = (eps_ox / eps) (A_i / L_D^(D-1)) / (t_ox / L_D)   (an edge of length t_ox)
//     psi_G,i = (V_G - Phi_i) / V_T
//     Phi_i   = phi_m - chi_i - (E_c - E_i),   E_c - E_i = Eg_i(T) / 2 + (k T / 2) ln(Nc_i / Nv_i)
//     S_i     = Q_f A_i / (Ns L_D^D)
// with A_i the node's face area on the gate's boundary patch, eps the scaling's reference
// permittivity, and chi_i, Eg_i, Nc_i, Nv_i of the node's material. psi is referenced to the
// intrinsic level E_i (n = n_ie e^psi), which lies (k T / 2) ln(Nv / Nc) above midgap; the vacuum
// level is chi + (E_c - E_i) above it, and psi_G is the potential, so referenced, that the
// electrode's Fermi level imposes on the oxide side. Band-gap narrowing (Slotboom) moves both band
// edges by half the narrowing and leaves E_i in place, so it does not enter.
//
// The legacy row is kappa (V_G - V_FB - (psi_0 - psi_b) V_T) / V_T with
// kappa = eps_ox L_D / (eps t_ox) and V_FB = phi_m - (chi + Eg / 2 - psi_b V_T) - q Q_f / C_ox
// (ARCHITECTURE.md 6.4, OLD / NEW / REASON):
// - it put the intrinsic level at midgap, which for silicon shifts the whole curve by
//   (k T / 2) ln(Nv / Nc) = 1.04 mV at 300 K (more for materials with unequal Nc and Nv);
// - it subtracted the neutral potential psi_b of the doping given for V_FB and added that of the
//   node, so the electrode's potential depended on the doping under it wherever the two differ.
// Neither is carried over; with a uniform substrate the rows differ only by the first.
//
// The charge on the gate electrode (per node, in units of q Ns L_D^D) is G_i (psi_G,i - psi_i):
// the displacement flux leaving the electrode. With the fixed charge it balances the
// semiconductor's charge (Gauss).
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/contact.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::assemble {

// phi_m [eV] of the gate electrode over material m at lattice temperature T. Preconditions as
// physics::band_gap_eV.
[[nodiscard]] double gate_work_function_eV(const device::GateStack& gate,
                                           const physics::Semiconductor& m, double temperature_K);

// Phi = phi_m - chi - (E_c - E_i) [V]: the electrode's potential at zero bias, referenced to the
// semiconductor's intrinsic level.
[[nodiscard]] double gate_intrinsic_offset_V(const device::GateStack& gate,
                                             const physics::Semiconductor& m,
                                             double temperature_K);

struct GateTerm {
    double coupling;      // G_i
    double offset;        // Phi_i / V_T, so psi_G,i = V_G / V_T - offset
    double sheet_charge;  // S_i
};

// The scaled gate term of a node with face area area_cm (cm^(D-1)) in a D-dimensional mesh.
[[nodiscard]] GateTerm gate_term(const device::GateStack& gate, const physics::Semiconductor& m,
                                 double area_cm, int dimension, const Scaling& scaling);

// phi_m [eV] of an `electrode` contact (device/contact.hpp: polysilicon is silicon's).
// Precondition: temperature_K finite and positive.
[[nodiscard]] double electrode_work_function_eV(const device::Electrode& electrode,
                                                double temperature_K);

// The gate nodes of a device, shared by both assemblers: per node, its gate contact (or -1), its
// term and the electrode potential psi_G at the gate's current bias.
class GateNodes {
public:
    GateNodes() = default;
    GateNodes(std::vector<std::int32_t> contact, std::vector<GateTerm> term, double V_T);

    [[nodiscard]] bool on_gate(std::size_t node) const noexcept { return contact_[node] >= 0; }
    // The gate contact of a gate node, in device.contacts() order.
    [[nodiscard]] std::size_t contact(std::size_t node) const noexcept {
        return static_cast<std::size_t>(contact_[node]);
    }
    [[nodiscard]] const GateTerm& term(std::size_t node) const noexcept { return term_[node]; }
    // psi_G of a gate node.
    [[nodiscard]] double potential(std::size_t node) const noexcept { return psi_gate_[node]; }

    // Sets psi_G of every gate node from the bias of its contact (bias_V in device.contacts()
    // order, already checked by check_contact_bias). Initially every bias is 0.
    void set_bias(std::span<const double> bias_V);

    // The node's oxide term G (psi_G - psi) + S; its derivative in psi is -G.
    [[nodiscard]] double row_term(std::size_t node, double psi) const noexcept {
        const GateTerm& g = term_[node];
        return g.coupling * (psi_gate_[node] - psi) + g.sheet_charge;
    }

    // Scaled charge on each gate electrode, the sum of G (psi_G - psi) over its nodes, in a vector
    // of `contacts` entries (zero for the other contacts). The potential of node i is
    // state[stride * i].
    [[nodiscard]] std::vector<double> charges(std::span<const double> state, std::size_t stride,
                                              std::size_t contacts) const;

private:
    std::vector<std::int32_t> contact_;
    std::vector<GateTerm> term_;
    std::vector<double> psi_gate_;
    double V_T_ = 1.0;
};

}  // namespace NiTCAD::assemble
