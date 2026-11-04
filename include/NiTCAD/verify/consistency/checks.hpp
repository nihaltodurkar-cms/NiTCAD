// Physical-consistency checks on NiTCAD results (Unit 55, TVR-1, DECISIONS 5.1, 11.6; TVD-16, TVD-20).
//
// Each check is a function of (result, mesh, device data) returning a CheckResult: the outcome, every part's
// measured residual and frozen tolerance, and the applicability predicates. Tolerances use only quantities the
// result carries (11.6); none is set after a run. A check never calls NiTCAD's physics: where it would need it,
// the part is NOT EVALUATED.
//
// The checks in this file (the DC ones):
//   C1  equilibrium: at zero bias with the thermal contacts at one temperature, every terminal current within its
//       own terminal_current_resolution; every edge's electron and hole current within the smallest nonzero
//       contact resolution; |F_n - F_p| <= 64 eps max(|E_c|, |E_v|, V_T) + V_T u_last at each semiconductor
//       node, u_last the point's last Newton update (owner ruling after measurement, DECISIONS 11.13); the
//       equilibrium solve's result is judged on the quasi-Fermi part alone, with u_last = 0.
//       (NiTCAD has no optical or prescribed generation (solve/bias.hpp), so that predicate of 11.6 holds by
//       construction; a BiasPoint is a steady state, so "no transient" holds by type.)
//   C2  at C1's predicates, each stored mechanism's q sum_i |R_k,i| V_i <= the smallest nonzero terminal
//       resolution (owner ruling TVR1-Q3); NOT EVALUATED when no rate is stored.
//   C3  dual-cell Gauss law, at a converged state: for each cell sum_f C_f (psi_j - psi_i) + q V_i (p - n + N_D -
//       N_A) = 0, C_f = eps_0 eps_f A_f / d_f with eps_f the harmonic mean of the ends' permittivities; tolerance
//       64 eps s_i + N_i u_last (s_i the row's operand scale, N_i its change over the last Newton update; owner
//       ruling, DECISIONS 11.13). Ohmic and Schottky contact cells are skipped (their charge is not stored);
//       each lumped gate, meshed electrode and charged interface is one group whose sum includes its stored
//       charge (gate_charge, q Q_f A, q sigma_f A, interface_trap_charge); global neutrality is the sum of
//       those residuals with each contact's charge from its cells' flux. NOT EVALUATED with incomplete
//       ionization (the ionized densities are not stored) or without a Newton record.
//   C5  steady energy balance with the electrothermal model: |sum_c V_c I_c - sum_k H_k| <= E_max sum_c R_c +
//       64 eps (sum_c |V_c I_c| + sum_k |H_k|), H_k the heat leaving through thermal contact k, R_c the terminal
//       resolutions, E_max = max(|V_c|; |E_c|, |E_v| on semiconductor nodes) + max |r + 5/2| k T_max (Unit 23's
//       U8 scale restated; owner ruling, DECISIONS 11.13). The global balance only: the per-node heat is not in
//       the results (a listed gap, not a part). NOT APPLICABLE with the model off or no thermal_contact_heat; a
//       BiasPoint is steady by type.
//   C6  transient charge bookkeeping, per accepted step k and contact c: beta h_k I_disp = Q_k - a1 Q_k-1 + a2 Q_k-2
//       by the step's own formula (backward Euler; variable-step BDF2 with w = h_k / h_k-1: beta = (1 + w)/(1 + 2w),
//       a1 = (1 + w)^2/(1 + 2w), a2 = w^2/(1 + 2w)), tolerance 64 eps times the operands' magnitudes; each gate and
//       meshed electrode without gate tunnelling carries conduction current 0 and terminal = displacement current,
//       exactly (owner rulings, DECISIONS 11.13). Diagnostics, not judged: the trapezoidal integral of I_disp
//       against Q_N - Q_0 per contact with sum_k error_ratio_k rtol |dQ_k| beside it, and the cumulative identity.
//       The identity checks bookkeeping, not charge physics (the solver forms I_disp from the same stored Q).
//   C7  reciprocity (Onsager) of the admittance at equilibrium: |Y_ij - Y_ji| <= resolution[k][i] +
//       resolution[k][j] at every frequency k of the point (f = 0 included), unnormalized; NOT APPLICABLE at a
//       bias not exactly 0 or thermal contacts at more than one temperature (C1's predicates). Every model is
//       judged (owner rulings, DECISIONS 11.13). With two contacts the symmetry follows from zero column and row
//       sums, so a diagnostic says the pass shows conservation, not reciprocity. No magnetic field exists.
//   C9  time-frequency consistency at low frequency, from four DC points about the AC point (DcStencil):
//       conductance, every row: Re Y_ij(0) against the fourth-order central difference D4 of the DC currents,
//       tolerance the difference's resolution term + |D4 - D2| (the stated truncation estimate) + 64 eps of its
//       terms + the column's AC resolution; capacitance, gate and electrode rows only (the DC results store no
//       ohmic contact's charge): Im Y_ij(f1)/omega1 against D4 of gate_charge, tolerance as the conductance's plus
//       resolution/omega1 and the omega^2 term from the two lowest nonzero frequencies. Kramers-Kronig on
//       (Y_ii - Y_ii(0))/(j omega) is a diagnostic, never judged (owner rulings, DECISIONS 11.13).
//   C10 limit continuity: over a caller-supplied x_k -> 0 (>= 3 values, else NOT APPLICABLE), the difference to
//       the limit model (a terminal current or the potential's L-infinity) above its floor at every value (else
//       NOT EVALUATED), and each observed rate between successive values within 0.2 of the stated one (owner
//       ruling, DECISIONS 11.13). The differences are reported as diagnostics.
//   C11 symmetry and invariance: two results related by a caller-supplied node and contact map, checked first to
//       be a symmetry of the device (bijections; doping, permittivity, regions equal; edges, lengths, couplings
//       and volumes to 64 eps; contacts and biases corresponding; else NOT APPLICABLE with the first failing item);
//       terminal currents to res_A + res_B + 64 eps max |I|, gate charges to 64 eps max |Q| plus the C3-style
//       Newton term (owner ruling, DECISIONS 11.13).
//   C12 n, p > 0 on semiconductor nodes; every stored temperature > 0; trap occupancies are not stored, so that
//       part is NOT EVALUATED when the device has interface traps (and vacuous without).
// Each outcome is from contract.hpp's judge: NOT APPLICABLE when a predicate is false (named, with its value),
// NOT EVALUATED when an input is missing, else pass or fail.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/results/small_signal.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/results/transient.hpp"
#include "NiTCAD/verify/consistency/device_data.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

namespace NiTCAD::verify::consistency {

inline constexpr double rounding_factor = 64.0;  // 64 eps (Unit 50's convention; DECISIONS 4.4, 11.6)

struct CheckPart {
    std::string name;
    Outcome outcome = Outcome::not_evaluated;  // pass, fail or NOT EVALUATED
    double residual = 0.0;                      // the worst measured residual of the part
    double tolerance = 0.0;                     // its frozen tolerance at that location
    std::string where;                          // the node, edge or contact of the worst case
    std::string reason;                         // for NOT EVALUATED, the missing input
};

// A measured value reported beside a check and never judged (DECISIONS 11.13, C6; TVD-21).
struct Diagnostic {
    std::string name;
    double value = 0.0;  // the measured difference
    double scale = 0.0;  // the quantity it is read against (NaN when not available)
    std::string where;
    std::string note;
};

struct CheckResult {
    std::string id;  // "C1", ...
    Comparison comparison;
    std::vector<CheckPart> parts;
    std::vector<Diagnostic> diagnostics;
};

// The contract of each check (class D, I2: a separate restatement of a discrete property, same arithmetic).
[[nodiscard]] ReferenceInfo check_info(std::string_view id);

[[nodiscard]] CheckResult check_c1(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d);
// The equilibrium solve's form: the quasi-Fermi part only, 64 eps (no Newton term).
[[nodiscard]] CheckResult check_c1(const results::EquilibriumResult& eq, const mesh::Mesh& m, const DeviceData& d);
[[nodiscard]] CheckResult check_c2(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d);
[[nodiscard]] CheckResult check_c3(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d);
[[nodiscard]] CheckResult check_c3(const results::EquilibriumResult& eq, const mesh::Mesh& m, const DeviceData& d);
[[nodiscard]] CheckResult check_c4(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d);
[[nodiscard]] CheckResult check_c5(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d);
[[nodiscard]] CheckResult check_c6(const results::Transient& run, const mesh::Mesh& m, const DeviceData& d);
// Point `point` of a small-signal run, all its frequencies.
[[nodiscard]] CheckResult check_c7(const results::SmallSignal& run, std::size_t point, const mesh::Mesh& m,
                                   const DeviceData& d);

// C9's DC points about the AC point: contact `contact` at V - 2h, V - h, V + h, V + 2h, every other bias equal.
struct DcStencil {
    std::size_t contact = 0;
    double h_V = 0.0;
    std::array<const results::BiasPoint*, 4> points{};
};
[[nodiscard]] CheckResult check_c9(const results::SmallSignal& run, std::size_t point, const DcStencil& stencil,
                                   const mesh::Mesh& m, const DeviceData& d);

// C10's series: the limit parameter x_k -> 0, a model result per value, and one limit-model result (or one per
// value, when the limit model's device moves with x, as the doping in Fermi-Dirac -> Boltzmann). The results live
// in the caller's storage.
enum class LimitQuantity : std::uint8_t {
    terminal_current,  // |I - I_limit| at `contact`; floor: the two resolutions
    potential,         // max over nodes |psi - psi_limit|; floor: 64 eps max |psi|
};
struct LimitSeries {
    std::string parameter;  // its name and unit, for the report
    std::vector<double> x;
    std::vector<const results::BiasPoint*> model;
    std::vector<const results::BiasPoint*> limit;
    LimitQuantity quantity = LimitQuantity::terminal_current;
    std::size_t contact = 0;
    double stated_rate = 1.0;  // derived before measurement (DECISIONS 11.13)
};
[[nodiscard]] CheckResult check_c10(const LimitSeries& series, const mesh::Mesh& m, const DeviceData& d);

// C11's maps from result A's device to result B's: node i of A is node[i] of B, contact c of A is contact[c] of B.
struct SymmetryMap {
    std::vector<mesh::NodeId> node;
    std::vector<std::size_t> contact;
};
[[nodiscard]] CheckResult check_c11(const results::BiasPoint& a, const mesh::Mesh& ma, const DeviceData& da,
                                    const results::BiasPoint& b, const mesh::Mesh& mb, const DeviceData& db,
                                    const SymmetryMap& map);
[[nodiscard]] CheckResult check_c12(const results::NodeFields& fields, const results::BiasPoint* point,
                                    const mesh::Mesh& m, const DeviceData& d);

}  // namespace NiTCAD::verify::consistency
