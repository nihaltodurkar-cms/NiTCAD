// The lattice heat of a system (Unit 23, models.electrothermal; Unit 57, DECISIONS.md GE3;
// system.hpp's header comment, "Electrothermal coupling").
#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/base/error.hpp"

namespace NiTCAD::assemble {

// Whether the model is on.
[[nodiscard]] bool electrothermal(const System& system) noexcept;

// The temperature of each thermal contact in K (isothermal: its own; resistance: the ambient
// behind it), in the order of device.thermal_contacts(); initially the device's values. Errors
// (invalid_input, the index the thermal contact, the value it): the size is not the thermal
// contact count; a value not finite and positive. (The solve layer checks a temperature against
// the materials' range.) Precondition (NITCAD_EXPECTS): the model is on.
[[nodiscard]] std::expected<void, base::Error> set_thermal_bias(
    System& system, std::span<const double> temperature_K);
// Whether every connected part of the mesh (every edge conducts heat) has a node of a thermal
// contact. A steady state needs one (DECISIONS.md T7); a transient does not. False without the
// model.
[[nodiscard]] bool heat_sinks_complete(const System& system) noexcept;
// Whether the connected thermal domain of a node has a node of a thermal contact (DECISIONS.md
// U5): a region without a sink of its own but joined to one has. Preconditions (NITCAD_EXPECTS):
// the model is on; node < node_count().
[[nodiscard]] bool reaches_heat_sink(const System& system, std::size_t node);
// K0 = kappa_ref T0 / (V_T J0 L_D) (system.hpp's header comment): the heat rows are the nodes'
// heat balances in units of V_T J0 L_D^(D-1) divided by it. 0 without the model.
[[nodiscard]] double heat_scale(const System& system) noexcept;
// Scaled heat leaving through each thermal contact at state x (the heat its isothermal nodes'
// rows would balance without their Dirichlet terms, or its R_th flux), in units of
// V_T J0 L_D^(D-1) (a power in W / cm^(3-D) after multiplying by it), in
// device.thermal_contacts() order. In a steady state their sum is electrical_power(x).
// Preconditions (NITCAD_EXPECTS): the model is on; x has unknowns() entries. With a time step,
// the heat rows' storage terms are in the balance.
[[nodiscard]] std::vector<double> thermal_contact_heat(const System& system,
                                                       std::span<const double> x,
                                                       const System::TimeStep* step = nullptr);
// The power the ohmic contacts deliver, sum_c V_c I_c, in the same units (also with energy
// balance, Unit 52).
[[nodiscard]] double electrical_power(const System& system, std::span<const double> x);
// The Joule part of each carrier edge's heat (Jn, Jp), in the same units: each flux times the
// drop of its quasi-Fermi level along it, theta_e (ln(w_b / w_a) - delta_n) for electrons and
// theta_e (ln(w_a / w_b) - delta_p) for holes; each >= 0, since a flux runs down that drop. Zero
// on edges without carriers. Preconditions as thermal_contact_heat.
[[nodiscard]] std::vector<std::pair<double, double>> edge_joule_heat(const System& system,
                                                                     std::span<const double> x);

}  // namespace NiTCAD::assemble
