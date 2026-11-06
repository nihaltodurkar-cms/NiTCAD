// The reference-model contract and comparison outcomes (Unit 55, TVR-1, DECISIONS 3.3, 4.4, 11.4; TVD-16, TVD-18).
//
// Plain data with validation of their own fields, and one judging rule that turns a criterion's result into an
// outcome the report prints. The rules enforced (owner rulings, DECISIONS 11.13, contract readings):
//   - a contract is invalid when a required field is empty (id, citation, version, equations, form, derivation,
//     verification), an equation number is outside G1-G15 or repeated, an E reference is not PA or PD with a second
//     method named, an E reference graded I1 names no second I1/I2 cross-check, or an A reference states no known
//     approximation error (3.3, 2.1);
//   - a comparison against an invalid contract is not_evaluated, the problems named;
//   - every applicability predicate the contract lists must be supplied (else not_evaluated, the missing one named)
//     and hold (else not_applicable, the predicate and its measured value named) before the criterion is judged;
//     a comparison outside its applicability is refused, not failed (TVD-6);
//   - an I0 reference cannot pass or fail: finding, "regression pin, not evidence" (2.1);
//   - a PC reference cannot pass or fail: finding, "provenance insufficient" (3.3, TVD-18); nor can an L reference
//     (legacy value or form, no source read; MI-Q5 of Unit 55 15.8, treated as PC).
// not_applicable and not_evaluated never count as a pass (TVD-16).
#pragma once

#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace NiTCAD::verify {

enum class ReferenceClass { E, D, A, M, H, P, X };  // DECISIONS 1.2
enum class Independence { I0, I1, I2, I3 };          // DECISIONS 2.1
enum class Provenance { PA, PB, PC, PD, L };         // DECISIONS 3.3; L: legacy, no source read (Unit 55 15.8, MI-Q5)
enum class Outcome { pass, finding, fail, not_applicable, not_evaluated };  // DECISIONS 4.4

inline constexpr int governing_equation_count = 15;  // G1-G15 (DECISIONS 1.3)

[[nodiscard]] std::string_view to_string(ReferenceClass c);
[[nodiscard]] std::string_view to_string(Independence g);
[[nodiscard]] std::string_view to_string(Provenance g);
[[nodiscard]] std::string_view to_string(Outcome o);

// An applicability predicate evaluated on the benchmark's parameters or the solved state.
struct Predicate {
    std::string name;
    double value = 0.0;      // the measured value
    double threshold = 0.0;  // the value it is compared with
    bool holds = false;
};

// A named output with its unit and the reference's own numerical uncertainty.
struct Quantity {
    std::string name;
    std::string unit;
    double value = 0.0;
    double u_ref = 0.0;  // >= 0
};

struct ReferenceInfo {
    std::string id;  // e.g. "TVB-02"
    ReferenceClass cls = ReferenceClass::E;
    Independence independence = Independence::I0;
    Provenance provenance = Provenance::PC;
    std::string citation;  // authors, journal, year, equation and table numbers, how accessed
    std::string version;
    std::vector<int> equations;  // the subset of G1-G15 (1-based)
    std::string form;            // the exact form solved
    std::vector<std::string> assumptions;
    std::vector<std::string> applicability;  // names of the predicates a comparison must supply
    std::string derivation;           // where the independent derivation is written
    std::string verification;         // the L0 evidence
    std::string approximation_error;  // A class: the bound or estimate and its source, or "measured"
    std::string second_method;        // E with PD: the second method; E with I1: the second I1/I2 cross-check
};

// Field problems, empty when valid.
[[nodiscard]] std::vector<std::string> validate(const ReferenceInfo& info);
[[nodiscard]] std::vector<std::string> validate(const Predicate& p);
[[nodiscard]] std::vector<std::string> validate(const Quantity& q);

struct Comparison {
    Outcome outcome = Outcome::not_evaluated;
    std::string reason;                // why (always set except for a plain pass)
    std::vector<Predicate> predicates;  // as supplied
};

// Judges one comparison against a reference: the contract, then applicability, then the criterion, then the
// independence and provenance rules above. `criterion` names the criterion for the reason text.
[[nodiscard]] Comparison judge(const ReferenceInfo& info, std::span<const Predicate> predicates, bool criterion_met,
                               std::string_view criterion);

// A comparison whose inputs are missing (no ladder, no stored rates, ...): not_evaluated, the input named.
[[nodiscard]] Comparison missing_input(std::string_view input);

// Field problems of an outcome: finding, fail, not_applicable and not_evaluated need a reason; not_applicable
// needs a supplied predicate that does not hold.
[[nodiscard]] std::vector<std::string> validate(const Comparison& c);

}  // namespace NiTCAD::verify
