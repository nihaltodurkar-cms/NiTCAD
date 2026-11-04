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
//   C12 n, p > 0 on semiconductor nodes; every stored temperature > 0; trap occupancies are not stored, so that
//       part is NOT EVALUATED when the device has interface traps (and vacuous without).
// Each outcome is from contract.hpp's judge: NOT APPLICABLE when a predicate is false (named, with its value),
// NOT EVALUATED when an input is missing, else pass or fail.
#pragma once

#include <string>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/results/solution.hpp"
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

struct CheckResult {
    std::string id;  // "C1", ...
    Comparison comparison;
    std::vector<CheckPart> parts;
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
[[nodiscard]] CheckResult check_c12(const results::NodeFields& fields, const results::BiasPoint* point,
                                    const mesh::Mesh& m, const DeviceData& d);

}  // namespace NiTCAD::verify::consistency
