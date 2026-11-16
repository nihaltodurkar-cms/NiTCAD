// The carriers' energy balance of a system (Unit 52, models.energy_balance; Unit 57,
// DECISIONS.md GE3): tau_n and tau_p are the fourth and fifth unknowns of every node.
#pragma once

#include <span>

#include "NiTCAD/assemble/system.hpp"

namespace NiTCAD::assemble {

// Whether the model is on.
[[nodiscard]] bool energy_balance(const System& system) noexcept;

// The carriers' energy budget of state x (Unit 52, EC7 as ruled 2026-10-29), in the units of
// electrical_power (k T_L per scaled particle flux): in a steady state
//     input = relaxation + recombination + contact_flux + contact_joule - peltier + residual
// to rounding, summed over the semiconductor nodes off the ohmic contacts:
// - input: sum_c V_c I_c (electrical_power);
// - relaxation: V k (3/2) (n tau_n / tau_wn + p tau_p / tau_wp), the power to the lattice;
// - recombination: V k U (E_c - E_v + (3/2)(theta_n + theta_p)), each recombined pair's gap and
//   thermal energy (U the continuity rows' net rate);
// - contact_flux: the energy flux the edges carry into the ohmic contact nodes (where tau = 0);
// - contact_joule: the half of the Joule heat of the contact edges no energy row keeps;
// - peltier: sum over the contact nodes of the carriers' band-edge energy above their Fermi
//   level times their flow into the device, (E_c + V) D_n + (E_v + V) D_p;
// - residual: the Newton residual's power, sum of the energy rows' residuals and of E_c and E_v
//   times the continuity rows' (zero for an exact solution).
// After a time step (Unit 52, EC13) the balance is the step's: input is the conduction power,
// recombination counts the rates alone, and two terms store energy:
//     input = ... + residual + thermal_storage + carrier_storage
// - thermal_storage: V rate (W_n + W_p - their history), the carriers' thermal energy (ED9);
// - carrier_storage: V rate (E_c (S_n - c_n) - E_v (S_p - c_p)), the band-edge energy of the
//   carriers stored (the continuity rows' storage terms, at the step's band edges).
// Both are 0 without a step. Preconditions (NITCAD_EXPECTS): the model is on and the
// temperatures are not held; x has unknowns() entries; step's spans of evaluate's sizes.
struct EnergyPowers {
    double input, relaxation, recombination, contact_flux, contact_joule, peltier, residual;
    double thermal_storage = 0.0, carrier_storage = 0.0;
};
[[nodiscard]] EnergyPowers energy_powers(const System& system, std::span<const double> x,
                                         const System::TimeStep* step = nullptr);

// The two-stage start of a bias point (Unit 52, ED12 as ruled 2026-10-29): hold_temperatures
// makes every energy row tau - tau_held = 0, tau_held the rises of state x (so psi, n and p are
// solved at the given carrier temperatures), until release_temperatures. Precondition
// (NITCAD_EXPECTS): the model is on; x has unknowns() entries.
void hold_temperatures(System& system, std::span<const double> x);
void release_temperatures(System& system) noexcept;
[[nodiscard]] bool temperatures_held(const System& system) noexcept;

}  // namespace NiTCAD::assemble
