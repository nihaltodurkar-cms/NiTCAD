// The reference-model contract and comparison outcomes (Unit 55, TVR-1, DECISIONS 11.4). See contract.hpp.
#include "NiTCAD/verify/reference/contract.hpp"

#include <algorithm>
#include <cmath>
#include <charconv>
#include <string>

namespace NiTCAD::verify {

std::string_view to_string(ReferenceClass c) {
    switch (c) {
        case ReferenceClass::E: return "E";
        case ReferenceClass::D: return "D";
        case ReferenceClass::A: return "A";
        case ReferenceClass::M: return "M";
        case ReferenceClass::H: return "H";
        case ReferenceClass::P: return "P";
        case ReferenceClass::X: return "X";
    }
    return "?";
}

std::string_view to_string(Independence g) {
    switch (g) {
        case Independence::I0: return "I0";
        case Independence::I1: return "I1";
        case Independence::I2: return "I2";
        case Independence::I3: return "I3";
    }
    return "?";
}

std::string_view to_string(Provenance g) {
    switch (g) {
        case Provenance::PA: return "PA";
        case Provenance::PB: return "PB";
        case Provenance::PC: return "PC";
        case Provenance::PD: return "PD";
    }
    return "?";
}

std::string_view to_string(Outcome o) {
    switch (o) {
        case Outcome::pass: return "pass";
        case Outcome::finding: return "finding";
        case Outcome::fail: return "fail";
        case Outcome::not_applicable: return "NOT APPLICABLE";
        case Outcome::not_evaluated: return "NOT EVALUATED";
    }
    return "?";
}

std::vector<std::string> validate(const ReferenceInfo& info) {
    std::vector<std::string> p;
    auto need = [&p](const std::string& field, const char* name) {
        if (field.empty()) p.push_back(std::string(name) + " empty");
    };
    need(info.id, "id");
    need(info.citation, "citation");
    need(info.version, "version");
    need(info.form, "form");
    need(info.derivation, "derivation");
    need(info.verification, "verification");
    if (info.equations.empty()) p.push_back("equations empty");
    for (std::size_t i = 0; i < info.equations.size(); ++i) {
        const int g = info.equations[i];
        if (g < 1 || g > governing_equation_count) p.push_back("equation G" + std::to_string(g) + " outside G1-G15");
        if (std::count(info.equations.begin(), info.equations.begin() + static_cast<std::ptrdiff_t>(i), g) > 0) {
            p.push_back("equation G" + std::to_string(g) + " repeated");
        }
    }
    for (const auto& a : info.assumptions) {
        if (a.empty()) p.push_back("empty assumption");
    }
    for (const auto& a : info.applicability) {
        if (a.empty()) p.push_back("empty applicability name");
    }
    if (info.cls == ReferenceClass::E) {
        if (info.provenance != Provenance::PA && info.provenance != Provenance::PD) {
            p.push_back("E reference graded " + std::string(to_string(info.provenance)) + " (needs PA or PD)");
        }
        if (info.provenance == Provenance::PD && info.second_method.empty()) {
            p.push_back("E reference graded PD names no second method");
        }
        if (info.independence == Independence::I1 && info.second_method.empty()) {
            p.push_back("E reference graded I1 names no second I1/I2 cross-check");
        }
    }
    if (info.cls == ReferenceClass::A && info.approximation_error.empty()) {
        p.push_back("A reference states no known approximation error");
    }
    return p;
}

std::vector<std::string> validate(const Predicate& pr) {
    std::vector<std::string> p;
    if (pr.name.empty()) p.push_back("predicate name empty");
    if (!std::isfinite(pr.value) || !std::isfinite(pr.threshold)) p.push_back("predicate '" + pr.name + "' not finite");
    return p;
}

std::vector<std::string> validate(const Quantity& q) {
    std::vector<std::string> p;
    if (q.name.empty()) p.push_back("quantity name empty");
    if (q.unit.empty()) p.push_back("quantity '" + q.name + "' unit empty");
    if (!std::isfinite(q.value)) p.push_back("quantity '" + q.name + "' value not finite");
    if (!std::isfinite(q.u_ref) || q.u_ref < 0.0) p.push_back("quantity '" + q.name + "' U_ref not finite and >= 0");
    return p;
}

namespace {

std::string join(const std::vector<std::string>& v) {
    std::string s;
    for (const auto& x : v) {
        if (!s.empty()) s += "; ";
        s += x;
    }
    return s;
}

// The shortest round-trip form of a double, for reason texts.
std::string number(double x) {
    char buf[32];
    const auto r = std::to_chars(buf, buf + sizeof buf, x);
    return std::string(buf, r.ptr);
}

}  // namespace

Comparison judge(const ReferenceInfo& info, std::span<const Predicate> predicates, bool criterion_met,
                 std::string_view criterion) {
    Comparison c;
    c.predicates.assign(predicates.begin(), predicates.end());
    if (const auto problems = validate(info); !problems.empty()) {
        c.outcome = Outcome::not_evaluated;
        c.reason = "contract invalid: " + join(problems);
        return c;
    }
    for (const auto& pr : predicates) {
        if (const auto problems = validate(pr); !problems.empty()) {
            c.outcome = Outcome::not_evaluated;
            c.reason = "predicate invalid: " + join(problems);
            return c;
        }
    }
    for (const auto& name : info.applicability) {
        const bool supplied = std::any_of(predicates.begin(), predicates.end(),
                                          [&name](const Predicate& pr) { return pr.name == name; });
        if (!supplied) {
            c.outcome = Outcome::not_evaluated;
            c.reason = "predicate '" + name + "' not supplied";
            return c;
        }
    }
    for (const auto& pr : predicates) {
        if (!pr.holds) {
            c.outcome = Outcome::not_applicable;
            c.reason = pr.name + " = " + number(pr.value) + " (threshold " + number(pr.threshold) + ")";
            return c;
        }
    }
    const std::string verdict = std::string(criterion) + (criterion_met ? " met" : " missed");
    if (info.independence == Independence::I0) {
        c.outcome = Outcome::finding;
        c.reason = "regression pin, not evidence (I0): " + verdict;
        return c;
    }
    if (info.provenance == Provenance::PC) {
        c.outcome = Outcome::finding;
        c.reason = "provenance insufficient (PC): " + verdict;
        return c;
    }
    c.outcome = criterion_met ? Outcome::pass : Outcome::fail;
    if (!criterion_met) c.reason = verdict;
    return c;
}

Comparison missing_input(std::string_view input) {
    Comparison c;
    c.outcome = Outcome::not_evaluated;
    c.reason = input.empty() ? std::string() : "missing input: " + std::string(input);
    return c;
}

std::vector<std::string> validate(const Comparison& c) {
    std::vector<std::string> p;
    if (c.outcome != Outcome::pass && c.reason.empty()) {
        p.push_back(std::string(to_string(c.outcome)) + " without a reason");
    }
    if (c.outcome == Outcome::not_applicable &&
        std::none_of(c.predicates.begin(), c.predicates.end(), [](const Predicate& pr) { return !pr.holds; })) {
        p.push_back("NOT APPLICABLE without a failing predicate");
    }
    for (const auto& pr : c.predicates) {
        for (auto& x : validate(pr)) p.push_back(std::move(x));
    }
    return p;
}

}  // namespace NiTCAD::verify
