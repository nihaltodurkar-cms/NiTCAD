// Validation against published data (Units 33-35; DECISIONS.md V1-V6, Unit 34, Unit 35): the
// benchmark decks in benchmarks/ against values taken from their cited sources. A disagreement is
// reported, not tuned: each gap is pinned at its measured value (V5). Units 33 and 34 ran on legacy
// silicon, as their pins still do; Unit 35 adds the defaults' own pins.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

// ---------------------------------------------------------------------------------------------
// W. R. Thurber, R. L. Mattis, Y. M. Liu, J. J. Filliben, "The relationship between resistivity
// and dopant density for phosphorus- and boron-doped silicon", NBS Special Publication 400-64
// (1981), doi:10.6028/NBS.SP.400-64. Read from the Internet Archive's text of the report
// (archive.org/details/relationshipbetw4006thur); DECISIONS.md V1, V2.

constexpr double q_fit = 1.602e-19;  // the charge the fits were made with (the report's footnote)

// Phosphorus, eq (1): log10(P / P0) = (A0 + A1 X + A2 X^2 + A3 X^3) / (1 + B1 X + B2 X^2 + B3 X^3),
// P = q rho N, P0 = 1 V s / cm^2; table 9.
struct Rational {
    std::array<double, 4> A;
    std::array<double, 3> B;
    double operator()(double X) const {
        const double num = A[0] + X * (A[1] + X * (A[2] + X * A[3]));
        const double den = 1.0 + X * (B[0] + X * (B[1] + X * B[2]));
        return std::pow(10.0, num / den);
    }
};
// X = log10(N / 1e16 cm^-3).
constexpr Rational phosphorus_of_N_300K{{-3.0652, 2.1853, -0.61080, 0.056189},
                                       {-0.67642, 0.19542, -0.018100}};
constexpr Rational phosphorus_of_N_23C{{-3.0769, 2.2108, -0.62272, 0.057501},
                                      {-0.68157, 0.19833, -0.018376}};
// X = log10(rho / 1 ohm cm) (only to check the transcription, V2).
constexpr Rational phosphorus_of_rho_300K{{-3.0951, -3.2303, -1.2024, -0.13679},
                                         {1.0205, 0.38382, 0.041338}};
constexpr Rational phosphorus_of_rho_23C{{-3.1083, -3.2626, -1.2196, -0.13923},
                                        {1.0265, 0.38755, 0.041833}};

// Boron, eqs (3) and (4): qrhoN = min + (max - min) / (1 + (x / x_ref)^a); tables 11 and 12.
// Table 12's N_ref exponents are lost in the scan's text; 18 is the only decade that keeps the two
// fits within the report's stated 3 percent of each other (V2).
struct MinMax {
    double min, max, ref, a;
    double operator()(double x) const { return min + (max - min) / (1.0 + std::pow(x / ref, a)); }
};
constexpr MinMax boron_of_N_300K{0.00215, 0.02053, 4.09e18, -0.727};
constexpr MinMax boron_of_N_23C{0.00209, 0.02024, 3.88e18, -0.737};
constexpr MinMax boron_of_rho_300K{0.00220, 0.01973, 0.01782, 1.086};
constexpr MinMax boron_of_rho_23C{0.00213, 0.01947, 0.01833, 1.105};

enum class Dopant { phosphorus, boron };

// The report's resistivity at dopant density N [ohm cm], at 300 K or 23 C.
double thurber_rho(Dopant d, double N, bool at_300K) {
    const double P = d == Dopant::phosphorus
                         ? (at_300K ? phosphorus_of_N_300K : phosphorus_of_N_23C)(std::log10(N / 1e16))
                         : (at_300K ? boron_of_N_300K : boron_of_N_23C)(N);
    return P / (q_fit * N);
}

// The report's other fit, the product as a function of resistivity, at that resistivity.
double thurber_product_of_rho(Dopant d, double rho, bool at_300K) {
    return d == Dopant::phosphorus
               ? (at_300K ? phosphorus_of_rho_300K : phosphorus_of_rho_23C)(std::log10(rho))
               : (at_300K ? boron_of_rho_300K : boron_of_rho_23C)(rho);
}

// ---------------------------------------------------------------------------------------------

// A benchmark and its measured gaps, simulated / reference - 1 (V5: pinned, not tuned): at 300 K,
// at 23 C, and at 300 K with incomplete ionization.
struct Benchmark {
    Dopant dopant;
    const char* file;
    double N;
    double gap_300K, gap_23C, gap_incomplete;
};

constexpr Benchmark resistors[] = {
    {Dopant::phosphorus, "phosphorus_1e14.deck", 1e14, 0.039, 0.040, 0.039},
    {Dopant::phosphorus, "phosphorus_1e15.deck", 1e15, 0.005, 0.005, 0.006},
    {Dopant::phosphorus, "phosphorus_1e16.deck", 1e16, -0.069, -0.072, -0.065},
    {Dopant::phosphorus, "phosphorus_1e17.deck", 1e17, -0.088, -0.097, -0.053},
    {Dopant::phosphorus, "phosphorus_1e18.deck", 1e18, 0.006, -0.015, 0.314},
    {Dopant::phosphorus, "phosphorus_1e19.deck", 1e19, -0.071, -0.074, 1.377},
    {Dopant::boron, "boron_1e15.deck", 1e15, -0.044, -0.043, -0.043},
    {Dopant::boron, "boron_1e16.deck", 1e16, 0.034, 0.035, 0.042},
    {Dopant::boron, "boron_1e17.deck", 1e17, 0.300, 0.295, 0.389},
    {Dopant::boron, "boron_1e18.deck", 1e18, 0.480, 0.463, 1.209},
    {Dopant::boron, "boron_1e19.deck", 1e19, 0.233, 0.228, 3.017},
};

// A pinned gap is met within 0.2 percentage points (the printed values are rounded to 0.1).
constexpr double pin = 0.002;

fs::path deck_path(const char* file) { return fs::path(NITCAD_BENCHMARKS_DIR) / "resistivity" / file; }

// The material the decks' silicon is run with (Unit 35): the deck's own `silicon` (the published
// band set, the default) or `silicon_legacy`, on which Units 33 and 34 measured their pins.
enum class Silicon { published, legacy };

// The deck, with its [device] temperature replaced when given, and its silicon the legacy set when
// asked.
deck::Deck load(const char* file, const char* temperature = nullptr,
                Silicon silicon = Silicon::published) {
    auto d = deck::read(deck_path(file));
    if (!d) FAIL(d.error().message);
    for (deck::Section& s : d->sections) {
        for (deck::Entry& e : s.entries) {
            if (temperature && s.kind == "device" && e.key == "temperature") e.value = temperature;
            if (silicon == Silicon::legacy && s.kind == "region" && e.key == "material") {
                REQUIRE(e.value == "silicon");
                e.value = "silicon_legacy";
            }
        }
    }
    return std::move(*d);
}

// The simulated resistivity of a resistor deck: 0.01 V over the 10 um bar. The hole mobility is
// the deck's (the default, automatic) unless given.
double simulated_rho(const deck::Deck& d, bool incomplete_ionization = false,
                     std::optional<assemble::HoleMobility> holes = std::nullopt,
                     assemble::IonizationLevels levels = assemble::IonizationLevels::fixed) {
    auto plan = deck::build(d);
    if (!plan) FAIL(plan.error().message);
    plan->options.models.incomplete_ionization = incomplete_ionization;
    plan->options.models.ionization_levels = levels;
    if (holes) plan->options.models.hole_mobility = *holes;
    const auto report = deck::run(*plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    REQUIRE_FALSE(report->sweep->stopped.has_value());
    const auto& p = report->sweep->points;
    REQUIRE(p.size() == 2);
    const double J = p[1].terminal_current[0] - p[0].terminal_current[0];  // A/cm^2
    return 0.01 / (J * 10e-4);
}

}  // namespace

TEST_CASE("validation: the Thurber fits as transcribed agree with themselves (V2)", "[validation]") {
    // Each fit's resistivity, put into the report's other fit, gives back its product within the
    // self-consistency the report and ASTM F723 state: boron 3 percent, phosphorus 4.5 percent.
    for (const bool at_300K : {true, false}) {
        for (const Benchmark& b : resistors) {
            const double rho = thurber_rho(b.dopant, b.N, at_300K);
            const double product = q_fit * rho * b.N;
            const double back = thurber_product_of_rho(b.dopant, rho, at_300K);
            CAPTURE(b.file, at_300K, rho, product, back);
            CHECK(std::abs(back / product - 1.0) < (b.dopant == Dopant::boron ? 0.03 : 0.045));
        }
    }
}

TEST_CASE("validation: resistivity against dopant density, Thurber et al. (V3)", "[validation]") {
    // Unit 33's baseline, on the legacy options it was measured with (Unit 35, ruling 4): legacy
    // silicon and the Caughey-Thomas hole mobility.
    const auto ct = assemble::HoleMobility::caughey_thomas;
    const auto legacy = Silicon::legacy;
    std::printf("resistivity benchmark (NBS SP 400-64), legacy options: N, reference, simulated, "
                "gap; 300 K | 23 C | 300 K with incomplete ionization\n");
    for (const Benchmark& b : resistors) {
        const double ref300 = thurber_rho(b.dopant, b.N, true);
        const double ref23 = thurber_rho(b.dopant, b.N, false);
        const double sim300 = simulated_rho(load(b.file, nullptr, legacy), false, ct);
        const double sim23 = simulated_rho(load(b.file, "296.15", legacy), false, ct);
        const double simII = simulated_rho(load(b.file, nullptr, legacy), true, ct);
        std::printf("%-20s %8.3g ohm cm: %8.4g %+6.1f%% | %8.4g vs %8.4g %+6.1f%% | %8.4g %+6.1f%%\n",
                    b.file, ref300, sim300, 100.0 * (sim300 / ref300 - 1.0), sim23, ref23,
                    100.0 * (sim23 / ref23 - 1.0), simII, 100.0 * (simII / ref300 - 1.0));
        CAPTURE(b.file);
        CHECK(std::abs(sim300 / ref300 - 1.0 - b.gap_300K) < pin);
        CHECK(std::abs(sim23 / ref23 - 1.0 - b.gap_23C) < pin);
        CHECK(std::abs(simII / ref300 - 1.0 - b.gap_incomplete) < pin);
    }
}

TEST_CASE("validation: the intrinsic carrier density, Misiakos and Tsamakis (V4)", "[validation]") {
    // K. Misiakos, D. Tsamakis, "Accurate measurements of the silicon intrinsic carrier density from
    // 78 to 340 K", J. Appl. Phys. 74, 3293 (1993), doi:10.1063/1.354551: "At 300 K the value of
    // n(i) was found to be (9.7 +/- 0.1) X 10(9) cm-3" (abstract). The model's n_i is read from
    // the solved 1e14 phosphorus bar at 0 V, sqrt(n p) at its middle node (no band-gap narrowing
    // at that doping).
    const double measured = 9.7e9, uncertainty = 0.1e9;
    const auto solved_ni = [&](Silicon silicon) {
        auto plan = deck::build(load("phosphorus_1e14.deck", nullptr, silicon));
        if (!plan) FAIL(plan.error().message);
        const auto report = deck::run(*plan);
        REQUIRE(report.has_value());
        REQUIRE(report->sweep.has_value());
        const auto& f = report->sweep->points.front().fields;
        const std::size_t mid = f.n_cm3.size() / 2;
        const double n_i = std::sqrt(f.n_cm3[mid] * f.p_cm3[mid]);
        std::printf("intrinsic density at 300 K, %s silicon: model %.5g cm^-3, measured %.3g +/- "
                    "%.1g (gap %+.1f%%, %.1f times the uncertainty)\n",
                    silicon == Silicon::legacy ? "legacy" : "published", n_i, measured, uncertainty,
                    100.0 * (n_i / measured - 1.0), (n_i - measured) / uncertainty);
        return n_i;
    };
    // The legacy silicon value, 1.0674e10 (Unit 5), lies 10% above the measurement: a finding (V5).
    CHECK(std::abs(solved_ni(Silicon::legacy) / measured - 1.0 - 0.100) < pin);
    // The default, Couderc et al.'s band set (Unit 35, DC4): within the measurement's uncertainty.
    const double n_i = solved_ni(Silicon::published);
    CHECK(std::abs(n_i / measured - 1.0) <= 0.011);
    CHECK(std::abs(n_i - measured) <= uncertainty);
}

// Unit 34: the measured Masetti gaps (simulated / reference - 1) at 300 K and 23 C, pinned like
// Unit 33's (V5), per boron density 1e15 ... 1e19.
constexpr std::array<double, 5> masetti_gap_300K{-0.011, -0.011, -0.053, -0.110, -0.039};
constexpr std::array<double, 5> masetti_gap_23C{-0.011, -0.011, -0.058, -0.124, -0.049};

// The Masetti gaps of the boron decks at 300 K and 23 C, and the phosphorus decks' (which must not
// move: electrons keep Caughey-Thomas).
struct MasettiRun {
    std::vector<double> boron_300K, boron_23C, phosphorus_300K, phosphorus_23C;
};

// Unit 34's runs, on the legacy silicon they were measured with (Unit 35).
MasettiRun run_masetti() {
    MasettiRun r;
    const auto masetti = assemble::HoleMobility::masetti;
    const auto legacy = Silicon::legacy;
    for (const Benchmark& b : resistors) {
        const double gap300 = simulated_rho(load(b.file, nullptr, legacy), false, masetti) /
                                  thurber_rho(b.dopant, b.N, true) - 1.0;
        const double gap23 = simulated_rho(load(b.file, "296.15", legacy), false, masetti) /
                                 thurber_rho(b.dopant, b.N, false) - 1.0;
        (b.dopant == Dopant::boron ? r.boron_300K : r.phosphorus_300K).push_back(gap300);
        (b.dopant == Dopant::boron ? r.boron_23C : r.phosphorus_23C).push_back(gap23);
    }
    return r;
}

TEST_CASE("validation: Masetti hole mobility against the boron benchmarks (Unit 34)", "[validation]") {
    // The acceptance criteria were fixed before this ran (DECISIONS.md, Unit 34). Checked here:
    // AC4 (at least 10 points better than Caughey-Thomas at 1e17, 1e18, 1e19 and a smaller largest
    // |gap|), AC6 (the phosphorus gaps unchanged), and AC3 (|gap| <= 10% at 300 K) and AC5 (<= 12% at
    // 23 C) at every boron density but 1e18, where they are NOT met (the next test case). Every
    // measured gap is also pinned.
    const MasettiRun r = run_masetti();
    std::printf("Masetti hole mobility (NBS SP 400-64): gap Caughey-Thomas -> Masetti, 300 K | 23 C\n");
    double worst_masetti = 0.0, worst_ct = 0.0;
    std::size_t boron = 0, phosphorus = 0;
    for (const Benchmark& b : resistors) {
        CAPTURE(b.file);
        if (b.dopant == Dopant::phosphorus) {
            CHECK(std::abs(r.phosphorus_300K[phosphorus] - b.gap_300K) < pin);  // AC6
            CHECK(std::abs(r.phosphorus_23C[phosphorus] - b.gap_23C) < pin);
            ++phosphorus;
            continue;
        }
        const double g300 = r.boron_300K[boron], g23 = r.boron_23C[boron];
        std::printf("%-20s %+6.1f%% -> %+6.1f%% | %+6.1f%% -> %+6.1f%%\n", b.file, 100.0 * b.gap_300K,
                    100.0 * g300, 100.0 * b.gap_23C, 100.0 * g23);
        CHECK(std::abs(g300 - masetti_gap_300K[boron]) < pin);
        CHECK(std::abs(g23 - masetti_gap_23C[boron]) < pin);
        if (b.N != 1e18) {
            CHECK(std::abs(g300) <= 0.10);  // AC3
            CHECK(std::abs(g23) <= 0.12);   // AC5
        }
        if (b.N >= 1e17) CHECK(std::abs(g300) <= std::abs(b.gap_300K) - 0.10);  // AC4
        worst_masetti = std::max(worst_masetti, std::abs(g300));
        worst_ct = std::max(worst_ct, std::abs(b.gap_300K));
        ++boron;
    }
    CHECK(worst_masetti < worst_ct);  // AC4
}

TEST_CASE("validation: Masetti AC3 and AC5 at boron 1e18 as fixed in advance (not met, Unit 34)",
          "[validation][!shouldfail]") {
    // Expected to fail, and reported as a failure of the predefined criteria: at 1e18 cm^-3 the
    // Masetti gap is -11.0% at 300 K (AC3 asks <= 10%) and -12.4% at 23 C (AC5 asks <= 12%). The
    // criteria are not loosened after the fact; if this case starts passing, Catch2 reports it.
    const MasettiRun r = run_masetti();
    CHECK(std::abs(r.boron_300K[3]) <= 0.10);  // AC3 at 1e18
    CHECK(std::abs(r.boron_23C[3]) <= 0.12);   // AC5 at 1e18
}

// Unit 35: the defaults (Couderc et al.'s silicon band set, the automatic hole mobility, Masetti on
// silicon) against the same benchmarks, each gap pinned at its measured value (V5), per deck in
// the order of `resistors`: 300 K, 23 C, and 300 K with incomplete ionization.
struct DefaultGaps {
    double gap_300K, gap_23C, gap_incomplete;
};
constexpr DefaultGaps default_gaps[] = {
    {0.039, 0.040, 0.039},     // phosphorus 1e14
    {0.005, 0.005, 0.006},     // 1e15
    {-0.069, -0.072, -0.065},  // 1e16
    {-0.088, -0.097, -0.053},  // 1e17
    {0.006, -0.015, 0.314},    // 1e18
    {-0.071, -0.074, 1.378},   // 1e19
    {-0.011, -0.011, -0.010},  // boron 1e15
    {-0.011, -0.011, -0.002},  // 1e16
    {-0.053, -0.058, 0.026},   // 1e17
    {-0.110, -0.124, 0.403},   // 1e18
    {-0.039, -0.049, 2.414},   // 1e19
};

TEST_CASE("validation: the defaults against the resistivity benchmarks (Unit 35, DC5)", "[validation]") {
    // DC5, fixed before this ran (DECISIONS.md, Unit 35): the boron gaps equal Unit 34's Masetti
    // pins and the phosphorus gaps Unit 33's Caughey-Thomas pins, each within 0.2 points (a
    // majority-carrier resistivity does not depend on n_i), and the largest boron |gap| is at most
    // 11.2%. Every measured gap, including those with incomplete ionization (whose occupancy does
    // depend on Nc and Nv), is pinned.
    std::printf("resistivity benchmark (NBS SP 400-64), defaults: gap at 300 K | 23 C | 300 K with "
                "incomplete ionization (legacy options in brackets)\n");
    double worst_boron = 0.0;
    std::size_t boron = 0;
    for (std::size_t k = 0; k < std::size(resistors); ++k) {
        const Benchmark& b = resistors[k];
        const double g300 = simulated_rho(load(b.file)) / thurber_rho(b.dopant, b.N, true) - 1.0;
        const double g23 =
            simulated_rho(load(b.file, "296.15")) / thurber_rho(b.dopant, b.N, false) - 1.0;
        const double gII = simulated_rho(load(b.file), true) / thurber_rho(b.dopant, b.N, true) - 1.0;
        std::printf("%-20s %+6.1f%% (%+6.1f%%) | %+6.1f%% (%+6.1f%%) | %+6.1f%% (%+6.1f%%)\n", b.file,
                    100.0 * g300, 100.0 * b.gap_300K, 100.0 * g23, 100.0 * b.gap_23C, 100.0 * gII,
                    100.0 * b.gap_incomplete);
        CAPTURE(b.file, g300, g23, gII);
        if (b.dopant == Dopant::boron) {
            CHECK(std::abs(g300 - masetti_gap_300K[boron]) < pin);  // DC5
            CHECK(std::abs(g23 - masetti_gap_23C[boron]) < pin);
            worst_boron = std::max(worst_boron, std::abs(g300));
            ++boron;
        } else {
            CHECK(std::abs(g300 - b.gap_300K) < pin);  // DC5
            CHECK(std::abs(g23 - b.gap_23C) < pin);
        }
        CHECK(std::abs(g300 - default_gaps[k].gap_300K) < pin);
        CHECK(std::abs(g23 - default_gaps[k].gap_23C) < pin);
        CHECK(std::abs(gII - default_gaps[k].gap_incomplete) < pin);
    }
    CHECK(worst_boron <= 0.112);  // DC5
}

// Unit 37: the gaps with incomplete ionization at Altermatt et al.'s density-dependent levels, on
// the defaults, per deck in the order of `resistors` (pinned, V5).
constexpr double altermatt_gaps[] = {0.039, 0.006, -0.065, -0.052, 0.246, 0.028,    // phosphorus
                                     -0.010, -0.002, 0.021, 0.176, 0.043};  // boron

TEST_CASE("validation: incomplete ionization at density-dependent levels (Unit 37, II4)",
          "[validation]") {
    // II4, fixed before this ran (DECISIONS.md, Unit 37): at 1e18 and 1e19 cm^-3 the gaps with
    // Altermatt's levels are smaller in magnitude than the fixed 45 meV level's (Unit 35's pins),
    // for both dopants; every gap is pinned.
    std::printf("resistivity benchmark (NBS SP 400-64), defaults with incomplete ionization: gap with "
                "Altermatt's levels (the fixed level's in brackets)\n");
    for (std::size_t k = 0; k < std::size(resistors); ++k) {
        const Benchmark& b = resistors[k];
        const double gap = simulated_rho(load(b.file), true, std::nullopt,
                                         assemble::IonizationLevels::altermatt) /
                               thurber_rho(b.dopant, b.N, true) - 1.0;
        const double fixed = default_gaps[k].gap_incomplete;
        std::printf("%-20s %+7.1f%% (%+7.1f%%)\n", b.file, 100.0 * gap, 100.0 * fixed);
        CAPTURE(b.file, gap, fixed);
        if (b.N >= 1e18) CHECK(std::abs(gap) < std::abs(fixed));  // II4
        CHECK(std::abs(gap - altermatt_gaps[k]) < pin);
    }
}
