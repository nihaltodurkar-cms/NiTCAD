// C5, the steady energy balance with the electrothermal model (DECISIONS Unit 55 section 11.6 and the three C5
// rulings of 11.13): passes on real electrothermal solves of a 1D pn diode written here (an isothermal and an R_th
// sink, Unit 23's U8 device), at equilibrium and forward bias; is NOT APPLICABLE on examples/pn_diode_1d.deck
// (model off) and without stored heat; each injected defect fails it; each missing input makes it NOT EVALUATED;
// and the reported tolerance equals the ruled one restated here. A hidden case ([.finding]) reports the one real
// solve found to miss the ruled bound (a fully compensated bar, DECISIONS 11.14).
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace c = NiTCAD::verify::consistency;
namespace v = NiTCAD::verify;
using v::Outcome;

namespace {

struct Solved {
    NiTCAD::deck::Plan plan;
    NiTCAD::deck::RunReport report;
    c::DeviceData data;
    const NiTCAD::mesh::Mesh& mesh() const { return plan.device.mesh(); }
    const NiTCAD::results::BiasPoint& point(std::size_t k) const { return report.sweep->points.at(k); }
};

std::optional<Solved> solve(std::expected<NiTCAD::deck::Deck, NiTCAD::base::Error> deck,
                            std::vector<std::vector<double>> points) {
    if (!deck) return std::nullopt;
    auto plan = NiTCAD::deck::build(*deck);
    if (!plan) return std::nullopt;
    plan->outputs = {};
    plan->analysis = NiTCAD::deck::Analysis::sweep;
    plan->points = std::move(points);
    auto report = NiTCAD::deck::run(*plan);
    if (!report) return std::nullopt;
    Solved s{std::move(*plan), std::move(*report), {}};
    s.data = NiTCAD::verify::benchmark::device_data(s.plan.device, s.plan.options.models);
    return s;
}

// Unit 23's U8 device mirrored: 1 um silicon, n 1e17 (net, an erfc step at 0.5 um) | p 1e17, ohmic at both ends
// (the anode on the p side at x_max), isothermal at x_min and a thermal resistance of 1e-3 K cm^2/W at x_max.
constexpr const char* diode_doping = R"(
profile = acceptor constant 1e17
profile = donor erfc 2e17 axis=x center=0.5um width=0.01um
)";
// The same bar fully compensated (N_D = N_A = 1e17 everywhere): a resistor, used by the hidden finding case with
// the anode at x_min (next to the isothermal sink); with the anode at x_max it passes (largest 1.5e-3 of the bound).
constexpr const char* resistor_doping = R"(
profile = acceptor constant 1e17
profile = donor constant 1e17
)";

std::string thermal_deck(const char* doping, const char* anode_patch, const char* extra = "") {
    const std::string cathode_patch = std::string(anode_patch) == "x_max" ? "x_min" : "x_max";
    return std::string(R"(
[mesh]
type = tensor
x = linspace(0, 1um, 41)

[region silicon]
material = silicon

[doping])") + doping + R"(
[contact anode]
kind = ohmic
patch = )" + anode_patch + R"(

[contact cathode]
kind = ohmic
patch = )" + cathode_patch + R"(

[models]
electrothermal = on

[thermal_contact left]
patch = x_min
kind = isothermal
temperature = 300

[thermal_contact right]
patch = x_max
kind = resistance
temperature = 300
resistance = 1e-3

[sweep]
contact anode = 0 : 0.8 : 0.1
)" + extra;
}

const std::vector<double> anode_volts{0.0, 0.3, 0.5, 0.7, 0.8};

std::optional<Solved> solve_bar(const char* doping, const char* anode_patch, const char* extra = "") {
    std::vector<std::vector<double>> points;
    for (double V : anode_volts) points.push_back({V, 0.0});
    return solve(NiTCAD::deck::parse(thermal_deck(doping, anode_patch, extra)), points);
}

void require_complete(const std::optional<Solved>& s) {
    REQUIRE(s.has_value());
    REQUIRE(s->report.sweep.has_value());
    if (s->report.sweep->stopped) std::printf("thermal sweep stopped: %s\n", s->report.sweep->stopped->message.c_str());
    REQUIRE(s->report.sweep->points.size() == anode_volts.size());
}

const Solved& diode() {
    static const auto s = solve_bar(diode_doping, "x_max");
    require_complete(s);
    return *s;
}

std::filesystem::path source(const char* file) { return std::filesystem::path(NITCAD_SOURCE_DIR) / file; }

// The ruled tolerance restated here, independently of checks.cpp (DECISIONS 11.13, C5 tolerance).
double ruled_tolerance(const NiTCAD::results::BiasPoint& p, const c::DeviceData& d) {
    double e = 0.0, resolution = 0.0, magnitude = 0.0;
    for (std::size_t k = 0; k < p.bias_V.size(); ++k) {
        e = std::max(e, std::abs(p.bias_V[k]));
        resolution += p.terminal_current_resolution[k];
        magnitude += std::abs(p.bias_V[k] * p.terminal_current[k]);
    }
    for (double h : p.thermal_contact_heat) magnitude += std::abs(h);
    double r = 0.0;
    for (std::size_t i = 0; i < d.insulator.size(); ++i) {
        if (d.insulator[i] != 0) continue;
        e = std::max({e, std::abs(p.bands.conduction_eV[i]), std::abs(p.bands.valence_eV[i])});
        r = std::max({r, std::abs(d.thermopower_n[i] + 2.5), std::abs(d.thermopower_p[i] + 2.5)});
    }
    const double t = *std::ranges::max_element(p.fields.temperature_K);
    e += r * NiTCAD::base::k_B_eV_per_K * t;
    return e * resolution + 64.0 * DBL_EPSILON * magnitude;
}

double power(const NiTCAD::results::BiasPoint& p) {
    double P = 0.0;
    for (std::size_t k = 0; k < p.bias_V.size(); ++k) P += p.bias_V[k] * p.terminal_current[k];
    return P;
}

}  // namespace

TEST_CASE("C5 contract", "[verify][consistency]") {
    const auto info = c::check_info("C5");
    CHECK(v::validate(info).empty());
    CHECK(info.cls == v::ReferenceClass::D);
    CHECK(info.applicability.size() == 2);
}

TEST_CASE("C5 passes on electrothermal pn diode solves", "[verify][consistency]") {
    const Solved& s = diode();
    // The benchmark copy: the run's model flag and silicon's thermopower exponents (-1/2 for both carriers, acoustic-
    // phonon scattering, physics/semiconductor.hpp) on every node.
    CHECK(s.data.electrothermal);
    REQUIRE(s.data.thermopower_n.size() == s.mesh().node_count());
    CHECK(std::ranges::all_of(s.data.thermopower_n, [](double r) { return r == -0.5; }));
    CHECK(std::ranges::all_of(s.data.thermopower_p, [](double r) { return r == -0.5; }));
    for (std::size_t k = 0; k < anode_volts.size(); ++k) {
        const auto& p = s.point(k);
        const auto r = c::check_c5(p, s.mesh(), s.data);
        REQUIRE(r.parts.size() == 1);
        const auto& part = r.parts[0];
        std::printf("C5 %.1f V: power %.6e heat %.6e W/cm^2, residual %.3e tolerance %.3e (%.3g)\n", anode_volts[k],
                    power(p), p.thermal_contact_heat[0] + p.thermal_contact_heat[1], part.residual, part.tolerance,
                    part.residual / part.tolerance);
        CAPTURE(k, part.residual, part.tolerance);
        CHECK(r.comparison.outcome == Outcome::pass);
        // The reported tolerance is the ruled one.
        const double ruled = ruled_tolerance(p, s.data);
        CHECK(std::abs(part.tolerance - ruled) <= 1e-9 * ruled);
    }
    // The last point heats the device measurably: the check is not judged on a cold device alone.
    CHECK(power(s.point(anode_volts.size() - 1)) > 1.0);
    CHECK(*std::ranges::max_element(s.point(anode_volts.size() - 1).fields.temperature_K) > 300.005);
}

TEST_CASE("C5 is NOT APPLICABLE without the model or the heat", "[verify][consistency]") {
    static const auto iso = solve(NiTCAD::deck::read(source("examples/pn_diode_1d.deck")), {{0.0, 0.3}});
    REQUIRE(iso.has_value());
    const auto r = c::check_c5(iso->point(0), iso->mesh(), iso->data);
    CHECK(r.comparison.outcome == Outcome::not_applicable);
    CHECK(r.comparison.reason.find("electrothermal on") != std::string::npos);
    CHECK(r.parts.empty());

    const Solved& s = diode();
    auto bad = s.point(2);
    bad.thermal_contact_heat.clear();
    const auto h = c::check_c5(bad, s.mesh(), s.data);
    CHECK(h.comparison.outcome == Outcome::not_applicable);
    CHECK(h.comparison.reason.find("thermal_contact_heat present") != std::string::npos);

    auto off = s.data;
    off.electrothermal = false;
    CHECK(c::check_c5(s.point(2), s.mesh(), off).comparison.outcome == Outcome::not_applicable);
}

TEST_CASE("C5 reports NOT EVALUATED for missing inputs", "[verify][consistency]") {
    const Solved& s = diode();
    const auto& p = s.point(3);
    auto expect = [&](const NiTCAD::results::BiasPoint& q, const c::DeviceData& d, const char* what) {
        const auto r = c::check_c5(q, s.mesh(), d);
        CAPTURE(what, r.comparison.reason);
        CHECK(r.comparison.outcome == Outcome::not_evaluated);
        CHECK(r.comparison.reason.find(what) != std::string::npos);
    };
    {
        auto q = p;
        q.terminal_current_resolution.clear();
        expect(q, s.data, "resolutions not stored");
    }
    {
        auto q = p;
        q.terminal_current.pop_back();
        expect(q, s.data, "currents");
    }
    {
        auto q = p;
        q.thermal_contact_heat.pop_back();  // one of two thermal contacts
        expect(q, s.data, "not one per thermal contact");
    }
    {
        auto q = p;
        q.bands.valence_eV.clear();
        expect(q, s.data, "bands not stored");
    }
    {
        auto q = p;
        q.fields.temperature_K.clear();
        expect(q, s.data, "lattice temperature not stored");
    }
    {
        auto d = s.data;
        d.thermopower_p.clear();
        expect(p, d, "thermopower");
    }
}

TEST_CASE("C5 fails on injected defects", "[verify][consistency]") {
    const Solved& s = diode();
    const auto& p = s.point(3);  // 0.7 V
    const double tol = ruled_tolerance(p, s.data);
    auto fails = [&](const NiTCAD::results::BiasPoint& q, const c::DeviceData& d, const char* what) {
        const auto r = c::check_c5(q, s.mesh(), d);
        CAPTURE(what, r.parts.empty() ? 0.0 : r.parts[0].residual, r.parts.empty() ? 0.0 : r.parts[0].tolerance);
        CHECK(r.comparison.outcome == Outcome::fail);
    };
    {   // one sink's heat moved by twice the tolerance
        auto q = p;
        q.thermal_contact_heat[1] += 2.0 * tol;
        fails(q, s.data, "heat + 2 tol");
    }
    {   // the same, downwards, on the other sink (an absolute value, not a signed bound)
        auto q = p;
        q.thermal_contact_heat[0] -= 2.0 * tol;
        fails(q, s.data, "heat - 2 tol");
    }
    {   // a sink's heat dropped (set to 0)
        auto q = p;
        q.thermal_contact_heat[1] = 0.0;
        fails(q, s.data, "R_th heat dropped");
    }
    {   // a terminal current moved by 2 tol / V on the biased contact
        auto q = p;
        q.terminal_current[0] += 2.0 * tol / q.bias_V[0];
        fails(q, s.data, "current");
    }
    {   // the heat's sign flipped (heat entering instead of leaving)
        auto q = p;
        for (auto& h : q.thermal_contact_heat) h = -h;
        fails(q, s.data, "sign");
    }
    {   // NaN heat
        auto q = p;
        q.thermal_contact_heat[0] = std::nan("");
        fails(q, s.data, "NaN heat");
    }
    {   // without resolutions the bound is the rounding term alone: an imbalance of 1e-12 relative fails
        auto q = p;
        for (auto& r : q.terminal_current_resolution) r = 0.0;
        q.thermal_contact_heat[0] = power(q) - q.thermal_contact_heat[1];
        const auto ok = c::check_c5(q, s.mesh(), s.data);
        CHECK(ok.comparison.outcome == Outcome::pass);
        CHECK(std::abs(ok.parts[0].tolerance - ruled_tolerance(q, s.data)) <= 1e-12 * ok.parts[0].tolerance);
        q.thermal_contact_heat[0] *= 1.0 + 1e-12;
        fails(q, s.data, "rounding term");
    }
}

TEST_CASE("C5 tolerance terms as ruled", "[verify][consistency]") {
    // The tolerance's terms, each moved alone against the restatement: the thermopower exponents (|r + 5/2|, both
    // carriers), the lattice temperature's maximum, the insulator nodes left out of E_max.
    const Solved& s = diode();
    const auto& p = s.point(4);
    auto agrees = [&](const NiTCAD::results::BiasPoint& q, const c::DeviceData& d, const char* what) {
        const auto r = c::check_c5(q, s.mesh(), d);
        REQUIRE(r.parts.size() == 1);
        const double ruled = ruled_tolerance(q, d);
        CAPTURE(what, r.parts[0].tolerance, ruled);
        CHECK(std::abs(r.parts[0].tolerance - ruled) <= 1e-12 * ruled);
    };
    agrees(p, s.data, "as solved");
    {
        auto d = s.data;
        for (auto& r : d.thermopower_n) r = 1.5;
        agrees(p, d, "r_n = 1.5");
    }
    {
        auto d = s.data;
        for (auto& r : d.thermopower_p) r = -6.0;  // |r + 5/2| = 3.5
        agrees(p, d, "r_p = -6");
    }
    {
        auto q = p;
        q.fields.temperature_K[q.fields.temperature_K.size() / 2] = 900.0;
        agrees(q, s.data, "T_max");
    }
    {
        auto q = p;
        auto d = s.data;
        q.bands.conduction_eV[7] = 50.0;  // larger than every other energy, on a node marked insulator
        d.insulator[7] = 1;
        agrees(q, d, "insulator node skipped");
        d.insulator[7] = 0;
        agrees(q, d, "semiconductor node counted");
    }
    {
        auto q = p;
        q.bias_V[1] = -40.0;  // a contact's |V| larger than every band energy
        agrees(q, s.data, "contact voltage");
    }
}

TEST_CASE("Finding C5: the compensated resistor at the default Newton tolerance", "[.finding][verify][consistency]") {
    // Hidden: reports the measured balance on the fully compensated bar (DECISIONS 11.14, C5 finding). The
    // terminal resolution there is ~5e-15 A/cm^2, so E_max R is ~1e-14 W/cm^2. Run at the default Newton tolerance
    // (1e-8) and at 1e-12: a miss that moves with the tolerance is the heat rows' Newton residual (U1: sum K0 f_T =
    // -sum V I at any state), which the ruled bound leaves out. Each point's outcome is checked.
    for (const char* extra : {"", "\n[solver]\ntol_update = 1e-12\n"}) {
        const auto s = solve_bar(resistor_doping, "x_min", extra);
        require_complete(s);
        for (std::size_t k = 0; k < anode_volts.size(); ++k) {
            const auto& p = s->point(k);
            const auto r = c::check_c5(p, s->mesh(), s->data);
            REQUIRE(r.parts.size() == 1);
            const auto& it = p.convergence.iterations;
            std::printf("C5 resistor (%s) %.1f V: power %.6e W/cm^2, R %.3e, iterations %zu, u_last %.3e, residual "
                        "%.3e tolerance %.3e (%.3g)\n",
                        *extra ? "tol 1e-12" : "tol 1e-8", anode_volts[k], power(p),
                        p.terminal_current_resolution[0], it.size(), it.empty() ? 0.0 : it.back().update,
                        r.parts[0].residual, r.parts[0].tolerance, r.parts[0].residual / r.parts[0].tolerance);
            CHECK(r.comparison.outcome == Outcome::pass);
        }
    }
}
