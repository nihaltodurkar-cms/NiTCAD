// Contract and outcome types (DECISIONS Unit 55 sections 3.3, 4.4, 11.4; owner rulings 11.13, contract readings).
//   - a filled contract validates; each required field, an out-of-range or repeated equation, an E reference not
//     PA or PD-with-two-methods, an E reference at I1 without a cross-check and an A reference without an
//     approximation error each make it invalid, and a comparison against it is NOT EVALUATED, the problem named;
//   - applicability: a listed predicate not supplied is NOT EVALUATED, one that does not hold is NOT APPLICABLE
//     with its measured value named, never fail;
//   - a PC reference cannot pass or fail (VR3), nor can an I0 one: both report finding;
//   - outcome field validation.
#include "NiTCAD/verify/reference/contract.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace v = NiTCAD::verify;
using v::Outcome;

namespace {

v::ReferenceInfo filled() {
    v::ReferenceInfo r;
    r.id = "TVB-03";
    r.cls = v::ReferenceClass::E;
    r.independence = v::Independence::I3;
    r.provenance = v::Provenance::PA;
    r.citation = "test citation, eq. (1)";
    r.version = "1";
    r.equations = {2, 3};
    r.form = "B(x) = x / (e^x - 1)";
    r.assumptions = {"none"};
    r.applicability = {"|x| bound"};
    r.derivation = "DECISIONS 11.8";
    r.verification = "series and closed form agree";
    return r;
}

bool contains(const std::string& s, const char* part) { return s.find(part) != std::string::npos; }

const std::array<v::Predicate, 1> holds{{{"|x| bound", 3.0, 700.0, true}}};

}  // namespace

TEST_CASE("a filled contract validates and passes on a met criterion", "[verify][contract]") {
    const auto r = filled();
    CHECK(v::validate(r).empty());
    const auto c = v::judge(r, holds, true, "4 eps relative");
    CHECK(c.outcome == Outcome::pass);
    CHECK(v::validate(c).empty());
    const auto f = v::judge(r, holds, false, "4 eps relative");
    CHECK(f.outcome == Outcome::fail);
    CHECK(f.reason == "4 eps relative missed");
    CHECK(v::validate(f).empty());
}

TEST_CASE("invalid contracts are not evaluated, the problem named", "[verify][contract]") {
    struct Case {
        const char* what;
        void (*mutate)(v::ReferenceInfo&);
        const char* expected;
    };
    const std::vector<Case> cases{
        {"id", [](v::ReferenceInfo& r) { r.id.clear(); }, "id empty"},
        {"citation", [](v::ReferenceInfo& r) { r.citation.clear(); }, "citation empty"},
        {"version", [](v::ReferenceInfo& r) { r.version.clear(); }, "version empty"},
        {"form", [](v::ReferenceInfo& r) { r.form.clear(); }, "form empty"},
        {"derivation", [](v::ReferenceInfo& r) { r.derivation.clear(); }, "derivation empty"},
        {"verification", [](v::ReferenceInfo& r) { r.verification.clear(); }, "verification empty"},
        {"equations", [](v::ReferenceInfo& r) { r.equations.clear(); }, "equations empty"},
        {"G0", [](v::ReferenceInfo& r) { r.equations = {0}; }, "G0 outside"},
        {"G16", [](v::ReferenceInfo& r) { r.equations = {16}; }, "G16 outside"},
        {"repeat", [](v::ReferenceInfo& r) { r.equations = {2, 5, 2}; }, "G2 repeated"},
        {"assumption", [](v::ReferenceInfo& r) { r.assumptions = {""}; }, "empty assumption"},
        {"applicability", [](v::ReferenceInfo& r) { r.applicability = {""}; }, "empty applicability"},
        {"E PB", [](v::ReferenceInfo& r) { r.provenance = v::Provenance::PB; }, "graded PB (needs PA or PD)"},
        {"E PC", [](v::ReferenceInfo& r) { r.provenance = v::Provenance::PC; }, "graded PC (needs PA or PD)"},
        {"E PD alone", [](v::ReferenceInfo& r) { r.provenance = v::Provenance::PD; }, "PD names no second method"},
        {"E I1 alone", [](v::ReferenceInfo& r) { r.independence = v::Independence::I1; }, "I1 names no second"},
        {"A no error",
         [](v::ReferenceInfo& r) {
             r.cls = v::ReferenceClass::A;
             r.provenance = v::Provenance::PB;
         },
         "no known approximation error"},
    };
    for (const auto& k : cases) {
        INFO(k.what);
        auto r = filled();
        k.mutate(r);
        const auto problems = v::validate(r);
        REQUIRE(problems.size() == 1);
        CHECK(contains(problems[0], k.expected));
        const auto c = v::judge(r, holds, true, "criterion");
        CHECK(c.outcome == Outcome::not_evaluated);
        CHECK(contains(c.reason, "contract invalid"));
        CHECK(contains(c.reason, k.expected));
    }
}

TEST_CASE("the E and A rules are satisfied by what they ask for", "[verify][contract]") {
    auto pd = filled();
    pd.provenance = v::Provenance::PD;
    pd.second_method = "shooting in double-double";
    CHECK(v::validate(pd).empty());
    auto i1 = filled();
    i1.independence = v::Independence::I1;
    i1.second_method = "I2 cross-check by quadrature";
    CHECK(v::validate(i1).empty());
    auto a = filled();
    a.cls = v::ReferenceClass::A;
    a.provenance = v::Provenance::PB;
    a.approximation_error = "(L_D / W)^2, stated bound";
    CHECK(v::validate(a).empty());
    // PB is allowed for A and M (printed in the report line), not for E.
    auto m = filled();
    m.cls = v::ReferenceClass::M;
    m.provenance = v::Provenance::PB;
    CHECK(v::validate(m).empty());
}

TEST_CASE("applicability: refused, not failed", "[verify][contract]") {
    const auto r = filled();
    SECTION("a listed predicate not supplied") {
        const auto c = v::judge(r, std::span<const v::Predicate>{}, false, "criterion");
        CHECK(c.outcome == Outcome::not_evaluated);
        CHECK(c.reason == "predicate '|x| bound' not supplied");
        CHECK(v::validate(c).empty());
    }
    SECTION("a predicate that does not hold") {
        const std::array<v::Predicate, 1> off{{{"|x| bound", 750.5, 700.0, false}}};
        for (const bool met : {true, false}) {
            const auto c = v::judge(r, off, met, "criterion");
            CHECK(c.outcome == Outcome::not_applicable);
            CHECK(c.reason == "|x| bound = 750.5 (threshold 700)");
            CHECK(v::validate(c).empty());
        }
    }
    SECTION("an extra caller predicate that does not hold also refuses") {
        const std::array<v::Predicate, 2> two{{{"|x| bound", 3.0, 700.0, true}, {"bias", 0.1, 0.0, false}}};
        CHECK(v::judge(r, two, true, "criterion").outcome == Outcome::not_applicable);
    }
    SECTION("a non-finite predicate is not evaluated") {
        const std::array<v::Predicate, 1> nan{{{"|x| bound", NAN, 700.0, true}}};
        const auto c = v::judge(r, nan, true, "criterion");
        CHECK(c.outcome == Outcome::not_evaluated);
        CHECK(contains(c.reason, "predicate invalid"));
    }
}

TEST_CASE("a PC reference cannot pass or fail (VR3)", "[verify][contract]") {
    auto r = filled();
    r.cls = v::ReferenceClass::M;  // PC is a valid grade to hold; it is the outcome that is limited
    r.provenance = v::Provenance::PC;
    REQUIRE(v::validate(r).empty());
    for (const bool met : {true, false}) {
        const auto c = v::judge(r, holds, met, "4 eps relative");
        CHECK(c.outcome == Outcome::finding);
        CHECK(c.reason == std::string("provenance insufficient (PC): 4 eps relative ") + (met ? "met" : "missed"));
        CHECK(v::validate(c).empty());
    }
}

TEST_CASE("an L reference cannot pass or fail, and is no E grade (Unit 55 15.8, MI-Q5)", "[verify][contract]") {
    auto r = filled();
    r.cls = v::ReferenceClass::M;
    r.provenance = v::Provenance::L;
    REQUIRE(v::validate(r).empty());
    CHECK(v::to_string(v::Provenance::L) == "L");
    for (const bool met : {true, false}) {
        const auto c = v::judge(r, holds, met, "bit-equal");
        CHECK(c.outcome == Outcome::finding);
        CHECK(c.reason ==
              std::string("provenance insufficient (L: legacy, no source read): bit-equal ") + (met ? "met" : "missed"));
        CHECK(v::validate(c).empty());
    }
    auto e = filled();
    e.provenance = v::Provenance::L;
    const auto problems = v::validate(e);
    REQUIRE(problems.size() == 1);
    CHECK(contains(problems[0], "graded L (needs PA or PD)"));
}

TEST_CASE("an I0 reference cannot pass or fail", "[verify][contract]") {
    auto r = filled();
    r.cls = v::ReferenceClass::M;
    r.independence = v::Independence::I0;
    r.provenance = v::Provenance::PA;
    REQUIRE(v::validate(r).empty());
    for (const bool met : {true, false}) {
        const auto c = v::judge(r, holds, met, "pin");
        CHECK(c.outcome == Outcome::finding);
        CHECK(contains(c.reason, "regression pin, not evidence (I0)"));
    }
    // Applicability still comes first.
    const std::array<v::Predicate, 1> off{{{"|x| bound", 800.0, 700.0, false}}};
    CHECK(v::judge(r, off, true, "pin").outcome == Outcome::not_applicable);
}

TEST_CASE("quantities, missing inputs and outcome validation", "[verify][contract]") {
    CHECK(v::validate(v::Quantity{"I_d", "A", 1e-3, 1e-15}).empty());
    CHECK(v::validate(v::Quantity{"", "A", 1e-3, 0.0}).size() == 1);
    CHECK(v::validate(v::Quantity{"I_d", "", 1e-3, 0.0}).size() == 1);
    CHECK(v::validate(v::Quantity{"I_d", "A", NAN, 0.0}).size() == 1);
    CHECK(v::validate(v::Quantity{"I_d", "A", 1.0, -1e-20}).size() == 1);
    CHECK(v::validate(v::Quantity{"I_d", "A", 1.0, INFINITY}).size() == 1);

    const auto m = v::missing_input("no tolerance ladder");
    CHECK(m.outcome == Outcome::not_evaluated);
    CHECK(m.reason == "missing input: no tolerance ladder");
    CHECK(v::validate(m).empty());
    CHECK_FALSE(v::validate(v::missing_input("")).empty());

    v::Comparison bad;
    bad.outcome = Outcome::not_applicable;
    bad.reason = "something";
    CHECK_FALSE(v::validate(bad).empty());  // no failing predicate
    bad.outcome = Outcome::finding;
    bad.reason.clear();
    CHECK_FALSE(v::validate(bad).empty());  // no reason

    CHECK(v::to_string(Outcome::not_applicable) == "NOT APPLICABLE");
    CHECK(v::to_string(Outcome::not_evaluated) == "NOT EVALUATED");
    CHECK(v::to_string(v::Provenance::PD) == "PD");
    CHECK(v::to_string(v::Independence::I3) == "I3");
    CHECK(v::to_string(v::ReferenceClass::H) == "H");
}
