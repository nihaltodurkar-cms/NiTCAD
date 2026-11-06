// Shared machinery of the TVB-31 model-implementation entries (Unit 55, TVR-3; DECISIONS 15.2, 15.5, 15.9).
// Internal to verify/benchmark.
//
//   - value_check: one comparison under 15.9's value rule, |q_sim - q_ref| <= max(1e-13 |q_ref|, 64 eps S(q)), with
//     q_ref the double-double reference f(x) rounded and S(q) = sum_i |x_i dq/dx_i| from condition_sum over every
//     input. The same function judges a partial: f is then the reference's analytic derivative (15.9, "Partials").
//     U_ref = 2^-80 (|q_ref| + S(q)): double-double evaluation (u_dd = 2^-104 per operation, a few dozen operations)
//     with margin for exp's subnormal low part near an underflow edge. By TB-Q4 a point where U_ref exceeds 1e-3 of
//     the bound is not_evaluated. Below 2^-969 the low part is subnormal: U_ref is floored at its spacing 2^-1074
//     (MI-Q19).
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
    // Below 2^-969 the double-double's low part is subnormal: its spacing 2^-1074 (and a subnormal q_ref's own
    // rounding) floors U_ref (DECISIONS Unit 55 15.18, MI-Q19).
    const double u_floor = q_ref != 0.0 && std::abs(q_ref) < 0x1p-969 ? 0x1p-1074 : 0.0;
    const double u_ref = std::fmax(0x1p-80 * (std::abs(q_ref) + S), u_floor);
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

// A documented cut (15.9, "Documented cuts"): NiTCAD returns exactly 0 where the code documents a cut, and that passes
// if |q_ref| there is at most |q_ref| at the cut's threshold. q is the uncut reference value at the point, q_threshold
// its value at the threshold.
[[nodiscard]] inline ComparisonRecord cut_record(const ReferenceInfo& info, std::string quantity, std::string unit,
                                                 double parameter, double q_sim, const DD& q, const DD& q_threshold,
                                                 const std::vector<Predicate>& predicates) {
    const double q_ref = q.to_double();
    ComparisonRecord c{std::move(quantity), std::move(unit), parameter, q_sim, q_ref, 0x1p-80 * std::abs(q_ref),
                       q_sim == 0.0 ? 0.0 : INFINITY,
                       "documented cut: q_sim = 0 and |q_ref| <= |q_ref at the threshold| (15.9)", {}};
    const bool ok = q_sim == 0.0 && std::abs(q.hi) <= std::abs(q_threshold.hi);
    c.comparison = judge(info, predicates, ok, c.criterion);
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

// The shipped value of a binding: legacy_silicon_parameters, silicon_parameters, a non-silicon set, SiO2 or
// another shipped constant.
[[nodiscard]] std::optional<double> shipped_value(std::string_view binding);

// The non-silicon constant sets by register name (15.19): germanium, gallium_arsenide, indium_gallium_arsenide,
// silicon_carbide_4h.
struct NamedSet {
    std::string_view name;
    physics::SemiconductorParameters params;
};
[[nodiscard]] std::array<NamedSet, 4> non_silicon_sets();

// Every field name bound_field resolves (the coverage rule of 15.19, MI-Q25).
[[nodiscard]] std::vector<std::string_view> field_names();

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

// A known finding defined by a rounding floor (DECISIONS Unit 55 15.15): a comparison that misses the value rule is
// the finding only where its error is within 64 eps T, T the size of the shipped form's terms; a larger miss stays a
// fail. classify_floor applies it to the last K comparisons (T[k] = 0 leaves k alone); note_if_stale adds the note
// when nothing was converted.
struct FloorFinding {
    std::string_view label, ruling;
    int converted = 0;
};

// The last K comparisons of r; T[k] = 0 for a quantity outside the finding.
template <std::size_t K>
void classify_floor(BenchmarkResult& r, const std::array<double, K>& T, FloorFinding& ff) {
    for (std::size_t k = 0; k < K; ++k) {
        ComparisonRecord& c = r.comparisons[r.comparisons.size() - K + k];
        if (c.comparison.outcome != Outcome::fail || !(T[k] > 0.0)) continue;
        if (!(std::abs(c.q_sim - c.q_ref) <= 64.0 * eps * T[k])) continue;
        c.comparison.outcome = Outcome::finding;
        c.comparison.reason = "known finding " + std::string(ff.label) + " (" + std::string(ff.ruling) +
                              "): " + c.comparison.reason + "; within 64 eps T, T " + shortest(T[k]);
        ++ff.converted;
    }
}

inline void note_if_stale(BenchmarkResult& r, const FloorFinding& ff) {
    if (ff.converted == 0) {
        r.notes.push_back("known finding " + std::string(ff.label) + " no longer reproduces: remove the entry");
    }
}

// A predicate that holds at exactly 0 or inside [lo, hi].
[[nodiscard]] inline Predicate zero_or_in_range(std::string name, double value, double lo, double hi) {
    Predicate p = in_range(std::move(name), value, lo, hi);
    if (value == 0.0) p.holds = true;
    return p;
}

}  // namespace NiTCAD::verify::benchmark::tvb31
