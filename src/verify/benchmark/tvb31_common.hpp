// Shared machinery of the TVB-31 model-implementation entries (Unit 55, TVR-3; DECISIONS 15.2, 15.5, 15.9).
// Internal to verify/benchmark.
//
//   - value_check: one comparison under 15.9's value rule, |q_sim - q_ref| <= max(1e-13 |q_ref|, 64 eps S(q)), with
//     q_ref the double-double reference f(x) rounded and S(q) = sum_i |x_i dq/dx_i| from condition_sum over every
//     input. The same function judges a partial: f is then the reference's analytic derivative (15.9, "Partials").
//     U_ref = 2^-80 (|q_ref| + S(q)): double-double evaluation (u_dd = 2^-104 per operation, a few dozen operations)
//     with margin for exp's subnormal low part near an underflow edge. By TB-Q4 a point where U_ref exceeds 1e-3 of
//     the bound is not_evaluated.
//   - table_check: one register row (verify/reference/model_parameters.hpp) against the shipped field it binds,
//     bit-equal to the correctly rounded parse of the printed decimal (15.9, "Parameter register"). A mismatch is a
//     finding recorded with its source; grade L and PC are findings at best (the contract); tag C is pinned (I0).
//   - register_parameters: a shipped set with the register's printed values of the given entries written over the
//     fields they bind, so the formula check runs NiTCAD's function on the source's numbers (15.2).
#pragma once

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <initializer_list>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/condition.hpp"
#include "NiTCAD/verify/reference/decimal.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/model_parameters.hpp"

namespace NiTCAD::verify::benchmark::tvb31 {

namespace ref = verify::reference;
namespace reg = verify::reference::model_parameters;
using ref::DD;

inline constexpr double eps = DBL_EPSILON;
inline constexpr std::string_view value_criterion = "|q_sim - q_ref| <= max(1e-13 |q_ref|, 64 eps S(q)) (15.9)";

// The value rule on a reference value q and its condition sum S (computed by the caller).
[[nodiscard]] inline ComparisonRecord value_record(const ReferenceInfo& info, std::string quantity, std::string unit,
                                                   double parameter, double q_sim, const DD& q, double S,
                                                   const std::vector<Predicate>& predicates) {
    const double q_ref = q.to_double();
    const double bound = std::fmax(1e-13 * std::abs(q_ref), 64.0 * eps * S);
    const double u_ref = 0x1p-80 * (std::abs(q_ref) + S);
    const double diff = std::abs((DD(q_sim) - q).hi);
    ComparisonRecord c{std::move(quantity), std::move(unit), parameter, q_sim, q_ref, u_ref,
                       bound > 0.0 ? diff / bound : (diff == 0.0 ? 0.0 : INFINITY), std::string(value_criterion), {}};
    const bool applicable = std::all_of(predicates.begin(), predicates.end(), [](const Predicate& p) { return p.holds; });
    if (applicable && u_ref > 1e-3 * bound) {
        c.comparison.outcome = Outcome::not_evaluated;
        c.comparison.predicates = predicates;
        c.comparison.reason = "U_ref " + shortest(u_ref) + " above 1e-3 of the bound " + shortest(bound) + " (TB-Q4)";
        return c;
    }
    c.comparison = judge(info, predicates, std::isfinite(q_sim) && diff <= bound, c.criterion);
    return c;
}

template <std::size_t N, class F>
[[nodiscard]] ComparisonRecord value_check(const ReferenceInfo& info, std::string quantity, std::string unit,
                                           double parameter, double q_sim, F&& f, const std::array<DD, N>& x,
                                           const std::vector<Predicate>& predicates) {
    return value_record(info, std::move(quantity), std::move(unit), parameter, q_sim, f(x), ref::condition_sum(f, x),
                        predicates);
}

// The NiTCAD field a register binding names ("legacy_silicon.Eg0_eV", "silicon.dos_mass.valence[3]"), in `p`;
// nullptr when the binding names another set or a field this resolver does not know.
[[nodiscard]] double* bound_field(physics::SemiconductorParameters& p, std::string_view binding,
                                  std::string_view set);

// The shipped value of a binding: legacy_silicon_parameters or silicon_parameters.
[[nodiscard]] std::optional<double> shipped_value(std::string_view binding);

// `base` (the shipped set named `set`) with the printed register values of `entries` written over their fields.
[[nodiscard]] physics::SemiconductorParameters register_parameters(physics::SemiconductorParameters base,
                                                                   std::string_view set,
                                                                   std::initializer_list<std::string_view> entries);

// The register's printed value of a binding as a double-double (exact decimal); throws if absent.
[[nodiscard]] DD register_dd(std::string_view binding);

// The table checks of every register row of `entry`.
void table_checks(BenchmarkResult& r, std::string_view entry);

// T from lo to hi on 64 points a decade, with hi, 300 K when inside, and each of `beyond` (reported
// not_applicable by the range predicate).
[[nodiscard]] std::vector<double> temperature_grid(double lo, double hi, std::initializer_list<double> beyond);

// lo, lo 10^(1/per), ..., hi (both ends included), per points a decade.
[[nodiscard]] std::vector<double> log_grid(double lo, double hi, int per);

// Points a decade on the dense grids: 64 (15.5) in every tier unless the close-out's Debug timing moves them
// (MI-Q9).
[[nodiscard]] int per_decade(Tier tier);

[[nodiscard]] Predicate in_range(std::string name, double value, double lo, double hi);

// An M reference contract for an entry (I1, G-equations `equations`, the register as version).
[[nodiscard]] inline ReferenceInfo entry_info(std::string id, Provenance g, std::string citation, std::string form,
                                              std::vector<std::string> applicability, std::vector<int> equations,
                                              std::string derivation, std::string verification) {
    ReferenceInfo r;
    r.id = std::move(id);
    r.cls = ReferenceClass::M;
    r.independence = Independence::I1;
    r.provenance = g;
    r.citation = std::move(citation);
    r.version = "register rows of verify/reference/model_parameters.hpp (DECISIONS 15.10)";
    r.equations = std::move(equations);
    r.form = std::move(form);
    r.assumptions = {"parameters: the register's printed decimals, exact in double-double"};
    r.applicability = std::move(applicability);
    r.derivation = std::move(derivation);
    r.verification = std::move(verification);
    return r;
}

// A predicate that holds at exactly 0 or inside [lo, hi].
[[nodiscard]] inline Predicate zero_or_in_range(std::string name, double value, double lo, double hi) {
    Predicate p = in_range(std::move(name), value, lo, hi);
    if (value == 0.0) p.holds = true;
    return p;
}

}  // namespace NiTCAD::verify::benchmark::tvb31
