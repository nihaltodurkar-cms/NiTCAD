// TVR-4 device-theory benchmarks (DECISIONS Unit 55 16.4; criteria 16.8 and 16.19, frozen before code): their
// frozen expected outcomes.
// TVB-06 (dielectric stack, insulator-only device, 16.22 DT-Q19 A) is D-class: C from both admittances and the
// field ratio pass to 64 eps N^2 at every level; fast runs 1D N = 32, 128 (16.7), full adds 2D 32 x 32 and runs in
// the Release-only ctest verify_benchmark.tvr4_full with TVB-05.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/benchmark/tvb32.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

namespace bench = NiTCAD::verify::benchmark;
using NiTCAD::verify::Outcome;

namespace {

void require_all_pass(const bench::BenchmarkResult& r, std::size_t expected_comparisons) {
    INFO(r.id);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    REQUIRE(r.comparisons.size() == expected_comparisons);
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.q_sim, c.q_ref, c.metric, c.criterion, c.comparison.reason);
        std::printf("%s %s: %.3e (%s)\n", r.id.c_str(), c.quantity.c_str(), c.metric, c.criterion.c_str());
        CHECK(c.comparison.outcome == Outcome::pass);
    }
    CHECK(bench::overall(r) == Outcome::pass);
}

// A frozen outcome recorded as found: per comparison, in order, its quantity, outcome, a part of its reason ("" for
// none) and, where the criterion ends in "= <bound>", that bound to 1e-6 (0: not pinned).
struct Found {
    std::string quantity;
    Outcome outcome;
    std::string reason;
    double bound = 0.0;
};

void require_as_found(const bench::BenchmarkResult& r, const std::vector<Found>& expected) {
    INFO(r.id);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("%s note: %s\n", r.id.c_str(), n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("%s [%s] %s: %.4g (%s) %s\n", r.id.c_str(),
                    std::string(to_string(c.comparison.outcome)).c_str(), c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const auto& c = r.comparisons[i];
        CAPTURE(c.quantity, c.metric, c.criterion, c.comparison.reason);
        CHECK(c.quantity == expected[i].quantity);
        CHECK(c.comparison.outcome == expected[i].outcome);
        CHECK(c.comparison.reason.find(expected[i].reason) != std::string::npos);
        if (expected[i].bound > 0.0) {
            const double bound = std::stod(c.criterion.substr(c.criterion.rfind("= ") + 2));
            CHECK(std::abs(bound - expected[i].bound) <= 1e-6 * expected[i].bound);
        }
    }
}

}  // namespace

TEST_CASE("TVB-06: the dielectric stack, fast tier (frozen outcome: pass)", "[tvr4]") {
    require_all_pass(bench::run_tvb06(bench::Tier::fast), 6);  // C twice and the field ratio, N = 32 and 128
}

TEST_CASE("TVB-06: the dielectric stack, full tier (frozen outcome: pass)", "[.tvr4_full]") {
    require_all_pass(bench::run_tvb06(bench::Tier::full), 9);  // + 2D 32 x 32
}

// TVB-13 (the uniform resistor, D; 16.2 (f), 16.8 frozen; construction 16.42). Frozen outcome, set before its first
// run: every comparison passes. Fast: the Canali current at four fields; full: + T_n - T_L at two fields, N = 128
// and 512, on the nodes beyond the layer's tail. Found on the first run (16.42): the drift passes; T_n - T_L fails at
// both fields and both N (the reference taken at V/L, not the interior's own field). DT-Q37 A (16.48): the reference
// is taken at each judged node's own field (16.2 (f)), and all 8 pass.
TEST_CASE("TVB-13: Canali drift, fast tier (frozen outcome: pass)", "[tvr4]") {
    require_all_pass(bench::run_tvb13(bench::Tier::fast), 4);
}

TEST_CASE("TVB-13: drift and energy balance, full tier (frozen outcome: pass)", "[.tvr4_full]") {
    require_all_pass(bench::run_tvb13(bench::Tier::full), 8);
}

// TVB-05 (coaxial capacitor, 2D generated, M = 32 ... 256 sectors) is E-class and full tier only (16.7). Release-only
// ctest verify_benchmark.tvr4_full (hidden tag), as tvb31_full. Recorded as found (16.24, DT-Q20 A): the finest
// accuracy passes; the order is not_evaluated, its successive orders 2.0545, 2.0825, 2.1846 spreading 0.102 > 0.1
// (the M = 32 pair), with p-hat 2.054 in the band. Neither the criteria nor the ladder changed (TVM2).
// The (h, tol) factorial (16.2 (p), DT-Q23 A) needs M = 512, about 4 minutes a solve, so it runs only in the nightly
// tier (verify_full.tvb05, local; 16.30). In full the finest accuracy is therefore held at not_evaluated (4.3).
TEST_CASE("TVB-05: the coaxial capacitor, full tier (frozen outcome: order and accuracy not_evaluated)",
          "[.tvr4_full]") {
    const auto r = bench::run_tvb05(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    REQUIRE(r.comparisons.size() == 2);
    const auto& order = r.comparisons[0];
    const auto& finest = r.comparisons[1];
    CAPTURE(order.metric, order.comparison.reason, finest.metric);
    std::printf("TVB-05 %s: %.4f (%s)\n", order.quantity.c_str(), order.metric, order.comparison.reason.c_str());
    std::printf("TVB-05 %s: %.3e (%s)\n", finest.quantity.c_str(), finest.metric, finest.criterion.c_str());
    CHECK(order.quantity == "C order");
    CHECK(order.comparison.outcome == Outcome::not_evaluated);
    CHECK(order.comparison.reason == "not in the asymptotic range");
    CHECK(order.metric >= 1.8);  // p-hat itself is in the band
    CHECK(order.metric <= 2.3);
    CHECK(finest.metric <= 1e-3);  // the accuracy itself (3.7e-5), not judged without the factorial
    CHECK(finest.comparison.outcome == Outcome::not_evaluated);
    CHECK(finest.comparison.reason == "the (h, tol) factorial was not evaluated");
    CHECK(bench::overall(r) == Outcome::not_evaluated);
}

TEST_CASE("TVB-05: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb05(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-15 (minority diffusion with SRH, the TVQ3 hook's G) is E (linear model, low injection) + A (injection ladder),
// full tier only (16.7): the order and finest accuracy of Delta n pass, and the departure's order in G and its
// G -> 0 limit pass (16.8, 16.19).
TEST_CASE("TVB-15: minority diffusion, full tier (frozen outcome: pass)", "[.tvr4_full]") {
    require_all_pass(bench::run_tvb15(bench::Tier::full), 4);
}

TEST_CASE("TVB-15: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb15(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-17 (depletion capacitance against the QE model, DT-Q4's predicate) is E, full tier only (16.7). Recorded as
// found (16.28, DT-Q22 A): at each of V = -0.5, -1, -2 V the finest accuracy of C_j passes, and the order is
// not_evaluated (spread 0.13-0.23 from the N = 64 pair, h about the Debye length) with p-hat in the band. Neither
// criteria, ladder nor construction changed (TVM2).
TEST_CASE("TVB-17: depletion capacitance, full tier (frozen outcome: orders not_evaluated, accuracy pass)",
          "[.tvr4_full]") {
    const auto r = bench::run_tvb17(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    REQUIRE(r.comparisons.size() == 6);
    for (std::size_t k = 0; k < 6; k += 2) {
        const auto& order = r.comparisons[k];
        const auto& finest = r.comparisons[k + 1];
        CAPTURE(order.quantity, order.metric, order.comparison.reason, finest.quantity, finest.metric);
        std::printf("%s: %.4f (%s)\n", order.quantity.c_str(), order.metric, order.comparison.reason.c_str());
        std::printf("%s: %.3e (%s)\n", finest.quantity.c_str(), finest.metric, finest.criterion.c_str());
        CHECK(order.comparison.outcome == Outcome::not_evaluated);
        CHECK(order.comparison.reason == "not in the asymptotic range");
        CHECK(order.metric >= 1.8);
        CHECK(order.metric <= 2.3);
        CHECK(finest.comparison.outcome == Outcome::pass);
    }
    CHECK(bench::overall(r) == Outcome::not_evaluated);
}

TEST_CASE("TVB-17: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb17(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-16 (ideal diode, three configurations; A, 16.19 frozen; DT-Q24 A, DT-Q25 A; construction 16.31) is full tier
// only (16.7). Frozen outcome: as found (16.35), 13 pass and 28 not_evaluated.
TEST_CASE("TVB-16: the ideal diode, full tier (frozen outcome: as found, 16.37)", "[.tvr4_full]") {
    const auto r = bench::run_tvb16(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("TVB-16 note: %s\n", n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("TVB-16 [%s] %s: %.4g (%s) %s\n", std::string(to_string(c.comparison.outcome)).c_str(),
                    c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == 41);
    // Frozen as found after DT-Q32 A (DECISIONS Unit 55 16.37): these pass; the long base at 0.5 and 0.6 V are known
    // finding (x) (DT-Q27 A); every other comparison is not_evaluated (U_disc not estimated: outside the asymptotic
    // range and the finest-triple order outside 2 +- 0.5 or oscillating; or a factorial state not reached).
    // The passes' bounds as found, to 1e-6: they pin U_num's wiring (U_disc, the factorial, eps).
    const std::vector<std::string> findings{"J (ii) SRH, long base, V = 0.5 V", "J (ii) SRH, long base, V = 0.6 V",
                                            "ideality (ii) SRH, long base, V = 0.5 .. 0.6 V"};
    const std::vector<std::pair<std::string, double>> passing{
        {"J (i) no G-R, short base, V = 0.2 V", 0.012969149188851229},
        {"J (i) no G-R, short base, V = 0.3 V", 0.012900968947520106},
        {"J (i) no G-R, short base, V = 0.4 V", 0.01284897206434447},
        {"J (i) no G-R, short base, V = 0.5 V", 0.012984518235497101},
        {"J (i) no G-R, short base, V = 0.6 V", 0.02062468841401274},
        {"ideality (i) no G-R, short base, V = 0.2 .. 0.3 V", 0.015638331230131353},
        {"ideality (i) no G-R, short base, V = 0.3 .. 0.4 V", 0.015576529070217962},
        {"ideality (i) no G-R, short base, V = 0.4 .. 0.5 V", 0.015640656838250177},
        {"ideality (i) no G-R, short base, V = 0.5 .. 0.6 V", 0.02048675753687821},
        {"J (ii) SRH, short base, V = 0.5 V", 0.012982402586626456},
        {"J (ii) SRH, short base, V = 0.6 V", 0.02062463318097297},
        {"J (ii) SRH, short base, V = -1 V", 0.012000282544006202},
        {"J (ii) SRH, short base, V = -2 V", 0.011961397873277196},
        {"ideality (ii) SRH, short base, V = 0.5 .. 0.6 V", 0.020490327411534195},
        {"J (ii) SRH, long base, V = 0.4 V", 0.0003116398018498728},
        {"J (ii) SRH, long base, V = -0.5 V", 0.0008452570169339834},
        {"J (ii) SRH, long base, V = -1 V", 0.0004328477970241768},
        {"J (ii) SRH, long base, V = -2 V", 0.00022225934829510192},
        {"reverse exponent (ii) SRH, long base, V = -0.5 .. -1 V", 0.0027973850856005764}};
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.metric, c.criterion, c.comparison.reason);
        const auto it = std::find_if(passing.begin(), passing.end(), [&](const auto& p) { return p.first == c.quantity; });
        const bool pass = it != passing.end();
        const bool finding = std::find(findings.begin(), findings.end(), c.quantity) != findings.end();
        CHECK(c.comparison.outcome == (pass ? Outcome::pass : finding ? Outcome::finding : Outcome::not_evaluated));
        if (pass) {
            const double bound = std::stod(c.criterion.substr(c.criterion.rfind("= ") + 2));
            CHECK(std::abs(bound - it->second) <= 1e-6 * it->second);
        }
    }
}

TEST_CASE("TVB-16: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb16(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-18 (the diode's admittance at forward bias, 5 um bases; A, 16.19 frozen; DT-Q29-Q31 A; construction 16.34) is
// full tier only (16.7). 63 comparisons: |Y|, Re Y, arg Y at 3 biases x 7 frequencies. Frozen outcome: as found after the
// DT-Q34 A reference correction (i w C_SCR, DECISIONS Unit 55 16.38): no fail; at 0.35 V |Y|, Re Y and arg Y at
// w tau_p = 0.001 .. 0.1 and Re Y at 1 and 10 are not_evaluated (U_num not estimated: differences changing sign, or
// arg Y's bound taking |Y|'s U_num); every other comparison passes. (Set before its first run: every comparison
// passes; that run's one fail was the reference's omission, DT-Q33 A.)
TEST_CASE("TVB-18: the diode admittance, full tier (frozen outcome: as found, 16.38)", "[.tvr4_full]") {
    const auto r = bench::run_tvb18(bench::Tier::full);
    if (r.run_error) FAIL(*r.run_error);
    for (const auto& n : r.notes) std::printf("TVB-18 note: %s\n", n.c_str());
    for (const auto& c : r.comparisons) {
        std::printf("TVB-18 [%s] %s: %.4g (%s) %s\n", std::string(to_string(c.comparison.outcome)).c_str(),
                    c.quantity.c_str(), c.metric, c.criterion.c_str(), c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == 63);
    const auto not_evaluated = [](const std::string& q) {
        if (q.find("V = 0.35 V") == std::string::npos) return false;
        for (const char* w : {"w tau_p = 0.001", "w tau_p = 0.01", "w tau_p = 0.1"}) {
            if (q.ends_with(w)) return true;
        }
        return q.starts_with("Re Y") && (q.ends_with("w tau_p = 1") || q.ends_with("w tau_p = 10"));
    };
    // The passes' bounds as found, to 1e-6 (as TVB-16's, 16.36): they pin U_num's wiring (U_disc, the factorial,
    // eps_Y, U_ref).
    const std::vector<std::pair<std::string, double>> passing{
        {"|Y| V = 0.35 V, w tau_p = 1", 8.808578696416149e-05},
        {"arg Y V = 0.35 V, w tau_p = 1", 0.014832918846696491},
        {"|Y| V = 0.35 V, w tau_p = 10", 0.0008383929445536035},
        {"arg Y V = 0.35 V, w tau_p = 10", 0.017068339655162565},
        {"|Y| V = 0.35 V, w tau_p = 100", 0.012940053191018656},
        {"Re Y V = 0.35 V, w tau_p = 100", 0.01218796097561107},
        {"arg Y V = 0.35 V, w tau_p = 100", 0.026406683468167724},
        {"|Y| V = 0.35 V, w tau_p = 1000", 0.2763215294156682},
        {"Re Y V = 0.35 V, w tau_p = 1000", 0.2688630994753438},
        {"arg Y V = 0.35 V, w tau_p = 1000", 0.05641374563828519},
        {"|Y| V = 0.4 V, w tau_p = 0.001", 0.0002339523335332707},
        {"Re Y V = 0.4 V, w tau_p = 0.001", 0.00023395233147367267},
        {"arg Y V = 0.4 V, w tau_p = 0.001", 0.010109126866218723},
        {"|Y| V = 0.4 V, w tau_p = 0.01", 0.0002339514578328694},
        {"Re Y V = 0.4 V, w tau_p = 0.01", 0.00023395125187305386},
        {"arg Y V = 0.4 V, w tau_p = 0.01", 0.010109062490821455},
        {"|Y| V = 0.4 V, w tau_p = 0.1", 0.00023404158423036216},
        {"Re Y V = 0.4 V, w tau_p = 0.1", 0.00023402099349257191},
        {"arg Y V = 0.4 V, w tau_p = 0.1", 0.010110303533754503},
        {"|Y| V = 0.4 V, w tau_p = 1", 0.00024311714675549916},
        {"Re Y V = 0.4 V, w tau_p = 1", 0.00024110992216568887},
        {"arg Y V = 0.4 V, w tau_p = 1", 0.010237231943646619},
        {"|Y| V = 0.4 V, w tau_p = 10", 0.0007738890490044848},
        {"Re Y V = 0.4 V, w tau_p = 10", 0.0006915658285820066},
        {"arg Y V = 0.4 V, w tau_p = 10", 0.013320211768472834},
        {"|Y| V = 0.4 V, w tau_p = 100", 0.012221447409990777},
        {"Re Y V = 0.4 V, w tau_p = 100", 0.011322861583325365},
        {"arg Y V = 0.4 V, w tau_p = 100", 0.022913256758529853},
        {"|Y| V = 0.4 V, w tau_p = 1000", 0.28187994411004363},
        {"Re Y V = 0.4 V, w tau_p = 1000", 0.27296328488109567},
        {"arg Y V = 0.4 V, w tau_p = 1000", 0.05291925488939189},
        {"|Y| V = 0.45 V, w tau_p = 0.001", 0.0013796196410596588},
        {"Re Y V = 0.45 V, w tau_p = 0.001", 0.00137961964066477},
        {"arg Y V = 0.45 V, w tau_p = 0.001", 0.008642224599368788},
        {"|Y| V = 0.45 V, w tau_p = 0.01", 0.0013796235262199786},
        {"Re Y V = 0.45 V, w tau_p = 0.01", 0.0013796234867318779},
        {"arg Y V = 0.45 V, w tau_p = 0.01", 0.008642248344870075},
        {"|Y| V = 0.45 V, w tau_p = 0.1", 0.0013796937420468888},
        {"Re Y V = 0.45 V, w tau_p = 0.1", 0.0013796897931653415},
        {"arg Y V = 0.45 V, w tau_p = 0.1", 0.008642628944188506},
        {"|Y| V = 0.45 V, w tau_p = 1", 0.0013880049117046638},
        {"Re Y V = 0.45 V, w tau_p = 1", 0.0013876102942143344},
        {"arg Y V = 0.45 V, w tau_p = 1", 0.008688737458754295},
        {"|Y| V = 0.45 V, w tau_p = 10", 0.0018167462270400497},
        {"Re Y V = 0.45 V, w tau_p = 10", 0.0017797293975423622},
        {"arg Y V = 0.45 V, w tau_p = 10", 0.010665824791732422},
        {"|Y| V = 0.45 V, w tau_p = 100", 0.013186925602938966},
        {"Re Y V = 0.45 V, w tau_p = 100", 0.012159244195724592},
        {"arg Y V = 0.45 V, w tau_p = 100", 0.021436684973540406},
        {"|Y| V = 0.45 V, w tau_p = 1000", 0.30621993551079657},
        {"Re Y V = 0.45 V, w tau_p = 1000", 0.295663707607053},
        {"arg Y V = 0.45 V, w tau_p = 1000", 0.05155669754630884},
    };
    int evaluated = 0;
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.metric, c.criterion, c.comparison.reason);
        const bool ne = not_evaluated(c.quantity);
        evaluated += ne ? 0 : 1;
        CHECK(c.comparison.outcome == (ne ? Outcome::not_evaluated : Outcome::pass));
        if (!ne) {
            const auto it =
                std::find_if(passing.begin(), passing.end(), [&](const auto& p) { return p.first == c.quantity; });
            REQUIRE(it != passing.end());
            const double bound = std::stod(c.criterion.substr(c.criterion.rfind("= ") + 2));
            CHECK(std::abs(bound - it->second) <= 1e-6 * it->second);
        }
    }
    CHECK(evaluated == 52);
    CHECK(passing.size() == 52);
}

TEST_CASE("TVB-18: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb18(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-10 (the MOS capacitor's HF C-V, A, and deep depletion, E; criteria 16.8 and 16.19, frozen; construction 16.40) is
// full tier only (16.7). 14 comparisons: C_HF at 6 gate biases; Q_s order and finest accuracy at 4. Frozen outcome: as
// found (DT-Q35 A, 16.48): 13 pass; C_HF at 2 V is NOT EVALUATED (U_num not estimated: the h/2 state's tolerance
// sensitivity; its |C_sim - C_ref| is 5.3e-13, as its neighbours'). The C_HF passes' bounds are pinned to 1e-6.
TEST_CASE("TVB-10: the MOS capacitor HF C-V and deep depletion, full tier (frozen outcome: as found, 16.48)",
          "[.tvr4_full]") {
    const Outcome pass = Outcome::pass;
    std::vector<Found> expected{
        {"C_HF, V_G - V_FB = 0.2 V", pass, "", 3.5149300759295247e-12},
        {"C_HF, V_G - V_FB = 0.5 V", pass, "", 2.302502110943355e-12},
        {"C_HF, V_G - V_FB = 1 V", pass, "", 1.6855634414557072e-12},
        {"C_HF, V_G - V_FB = 1.5 V", pass, "", 3.22814423343859e-10},
        {"C_HF, V_G - V_FB = 2 V", Outcome::not_evaluated, "U_num not estimated", 0.0},
        {"C_HF, V_G - V_FB = 2.5 V", pass, "", 4.740650305197262e-09}};
    for (const char* v : {"1", "2", "3", "4"}) {
        expected.push_back({std::string("Q_s order (deep depletion), V_G - V_FB = ") + v + " V", pass, "", 0.0});
        expected.push_back({std::string("Q_s finest (deep depletion, N = 512), V_G - V_FB = ") + v + " V", pass, "", 0.0});
    }
    require_as_found(bench::run_tvb10(bench::Tier::full), expected);
}

TEST_CASE("TVB-10: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb10(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

TEST_CASE("TVB-19: nothing in the fast tier", "[tvr4]") {
    const auto r = bench::run_tvb19(bench::Tier::fast);
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.empty());
}

// TVB-27 (lattice heat, E; 16.8 frozen; construction 16.43, the resistor's changed by DT-Q38 B, 16.48) and TVB-28
// (thermoelectric identities, E and D) are full tier only (16.7). TVB-27: the order and the finest accuracy of T for
// the Kirchhoff slab and the self-heated resistor (4). TVB-28: the same for S and Pi, and C5's energy residual on four
// states (8). Frozen outcome: as found (DT-Q38 B, 16.48).
// TVB-27: both accuracies pass; both orders are NOT EVALUATED (the slab is exact in h; the resistor, 100 um on the
// legacy set, has p-hat 2.26 with its successive orders spreading more than 0.1).
// TVB-28: 6 pass; S's order is NOT EVALUATED (its finest pair 0.01 below the window); Pi's order fails (Pi is exact
// in h: DT-Q36's class).
TEST_CASE("TVB-27: lattice heat, full tier (frozen outcome: as found, 16.48)", "[.tvr4_full]") {
    const char* range = "not in the asymptotic range";
    require_as_found(bench::run_tvb27(bench::Tier::full),
                     {{"T (Kirchhoff slab) order in h", Outcome::not_evaluated, range},
                      {"T (Kirchhoff slab) finest accuracy", Outcome::pass, ""},
                      {"T (self-heated resistor) order in h", Outcome::not_evaluated, range},
                      {"T (self-heated resistor) finest accuracy", Outcome::pass, ""}});
}

TEST_CASE("TVB-28: thermoelectric identities, full tier (frozen outcome: as found, 16.48)", "[.tvr4_full]") {
    require_as_found(bench::run_tvb28(bench::Tier::full),
                     {{"S (open-circuit rise) order in h", Outcome::not_evaluated, "not in the asymptotic range"},
                      {"S (open-circuit rise) finest accuracy", Outcome::pass, ""},
                      {"Pi (contact heat per current) order in h", Outcome::fail, "p-hat in [1.8, 2.3] (p_th 2) missed"},
                      {"Pi (contact heat per current) finest accuracy", Outcome::pass, ""},
                      {"energy residual (C5), self-heated resistor, 4 V", Outcome::pass, ""},
                      {"energy residual (C5), Seebeck bar", Outcome::pass, ""},
                      {"energy residual (C5), Peltier bar, +1 mV", Outcome::pass, ""},
                      {"energy residual (C5), Peltier bar, -1 mV", Outcome::pass, ""}});
}

TEST_CASE("TVB-27 and TVB-28: nothing in the fast tier", "[tvr4]") {
    for (const auto& r : {bench::run_tvb27(bench::Tier::fast), bench::run_tvb28(bench::Tier::fast)}) {
        REQUIRE_FALSE(r.run_error.has_value());
        REQUIRE(r.comparisons.empty());
    }
}

// TVB-32 on the TVR-4 benchmark devices (16.44; TVB-32's full tier). Frozen expectations, set before code and before
// any run (11.13's rule): C1, C2, C7, C8 pass where every bias is exactly 0 and the thermal contacts are at one
// temperature, else NOT APPLICABLE naming that predicate; C1 and C3 on the equilibrium, C3, C4, C12, C6 and C9 pass;
// C5 NOT APPLICABLE ("electrothermal on") with the model off and pass with it on (TVB-27, 28). Any other outcome,
// and any failed or NOT EVALUATED part, is a finding. Per device and check the worst residual/tolerance is printed.
// Frozen as found (DT-Q39 A, 16.48): 501 entries. Beside the rule above, the checks' own contracts (11.6) and one
// known finding are pinned by their counts and reasons:
//   C2 NOT EVALUATED, "no rates stored", at zero-bias points with G-R off (10);
//   C1 and C2 NOT EVALUATED, "no nonzero terminal resolution", at zero bias on the insulator-only devices and at
//   TVB-10's HF point (5); C12 NOT EVALUATED, "no semiconductor nodes", on the insulator-only devices (4);
//   C9 NOT EVALUATED at TVB-18 (no gate or electrode row for its capacitance part) (1);
//   C6 NOT APPLICABLE, "transient with a gate or electrode", on TVB-19's transients (3);
//   the finding (16.45; its fix deferred past TVR-4, DT-Q11): C1 and C2 fail at TVB-15's generated point (2).
TEST_CASE("TVB-32: the consistency suite on the TVR-4 devices (frozen outcomes)", "[.tvr4_full]") {
    namespace c = NiTCAD::verify::consistency;
    const auto name_of = [](Outcome o) {
        switch (o) {
            case Outcome::pass: return "pass";
            case Outcome::fail: return "FAIL";
            case Outcome::not_applicable: return "not applicable";
            case Outcome::not_evaluated: return "NOT EVALUATED";
            default: return "other";
        }
    };
    struct Tally {
        std::map<std::string, int> outcomes;
        double worst = 0.0;
        std::string where;
    };
    const bench::Tvb32Run run = bench::tvb32_tvr4_entries();
    CAPTURE(run.run_error.value_or(""));
    CHECK_FALSE(run.run_error.has_value());
    std::map<std::string, Tally> tally;  // "device check"
    std::map<std::string, int> pinned;   // the contract classes and the finding, by count
    const auto pinned_class = [](const bench::Tvb32Entry& e) -> std::string {
        const c::CheckResult& r = e.result;
        const auto has = [&](const char* s) { return r.comparison.reason.find(s) != std::string::npos; };
        const Outcome o = r.comparison.outcome;
        if (o == Outcome::not_evaluated && r.id == "C2" && has("no rates stored")) return "no rates stored";
        if (o == Outcome::not_evaluated && (r.id == "C1" || r.id == "C2") && has("no nonzero terminal resolution")) {
            return "no nonzero terminal resolution";
        }
        if (o == Outcome::not_evaluated && r.id == "C12" && has("no semiconductor nodes")) {
            return "no semiconductor nodes";
        }
        if (o == Outcome::not_evaluated && r.id == "C9" && e.deck.rfind("TVB-18", 0) == 0 &&
            has("capacitance (gate and electrode rows)")) {
            return "C9 without a gate or electrode row";
        }
        if (o == Outcome::not_applicable && r.id == "C6" && has("transient with a gate or electrode")) {
            return "C6 on a transient without a gate or electrode";
        }
        if (o == Outcome::fail && (r.id == "C1" || r.id == "C2") && e.deck.rfind("TVB-15", 0) == 0 &&
            e.at.find("generation") != std::string::npos) {
            return "finding: TVB-15's generated point";
        }
        return {};
    };
    for (const bench::Tvb32Entry& e : run.entries) {
        const c::CheckResult& r = e.result;
        const std::string& id = r.id;
        const bool thermal = e.deck.rfind("TVB-27", 0) == 0 || e.deck.rfind("TVB-28", 0) == 0;
        Outcome want = Outcome::pass;
        std::vector<std::string> predicates;
        if (e.at != "equilibrium" && (id == "C1" || id == "C2" || id == "C7" || id == "C8")) {
            want = e.zero_bias ? Outcome::pass : Outcome::not_applicable;
            predicates = {"all biases 0", "thermal contacts at one temperature"};
        } else if (id == "C5" && !thermal) {
            want = Outcome::not_applicable;
            predicates = {"electrothermal on"};
        }
        Tally& t = tally[e.deck + " " + id];
        ++t.outcomes[name_of(r.comparison.outcome)];
        for (const auto& p : r.parts) {
            if (p.outcome != Outcome::pass || !(p.tolerance > 0.0)) continue;
            const double q = p.residual / p.tolerance;
            if (q > t.worst) {
                t.worst = q;
                t.where = e.at + ", " + p.name + " at " + p.where;
            }
        }
        CAPTURE(e.deck, e.at, id, r.comparison.reason);
        if (const std::string cls = pinned_class(e); !cls.empty()) {
            ++pinned[cls];
            continue;
        }
        for (const auto& p : r.parts) {
            CAPTURE(p.name, name_of(p.outcome), p.residual, p.tolerance, p.where, p.reason);
            if (p.outcome == Outcome::fail || p.outcome == Outcome::not_evaluated) FAIL_CHECK("part " << p.name);
        }
        CHECK(r.comparison.outcome == want);
        if (want == Outcome::not_applicable) {
            CHECK(std::any_of(predicates.begin(), predicates.end(), [&](const std::string& s) {
                return r.comparison.reason.find(s) != std::string::npos;
            }));
        }
    }
    CHECK(run.entries.size() == 501);
    const std::map<std::string, int> expected{{"no rates stored", 10},
                                              {"no nonzero terminal resolution", 5},
                                              {"no semiconductor nodes", 4},
                                              {"C9 without a gate or electrode row", 1},
                                              {"C6 on a transient without a gate or electrode", 3},
                                              {"finding: TVB-15's generated point", 2}};
    for (const auto& [cls, n] : pinned) std::printf("TVB-32 pinned: %s: %d\n", cls.c_str(), n);
    CHECK(pinned == expected);
    for (const auto& [key, t] : tally) {
        std::string counts;
        for (const auto& [o, n] : t.outcomes) counts += (counts.empty() ? "" : ", ") + std::to_string(n) + " " + o;
        std::printf("TVB-32 %s: %s; worst %.3g of the bound (%s)\n", key.c_str(), counts.c_str(), t.worst,
                    t.where.c_str());
    }
}
