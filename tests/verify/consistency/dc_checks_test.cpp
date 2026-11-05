// DC consistency checks C1, C2, C4, C12 (DECISIONS Unit 55 section 11.6): each fails on its injected defect and
// passes, or is NOT APPLICABLE with the right predicate, on a real solve of examples/pn_diode_1d.deck (at 0 V and
// 0.3 V on the anode; the deck's own sweep is replaced through the API, the deck file is not edited).
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <string>

namespace c = NiTCAD::verify::consistency;
namespace v = NiTCAD::verify;
using v::Outcome;

namespace {

struct Solved {
    NiTCAD::deck::Plan plan;
    NiTCAD::results::Sweep sweep;
    c::DeviceData data;
};

// Solved once per test run (Catch2 runs the cases in one process).
const Solved& pn_diode() {
    static const std::optional<Solved> solved = [] () -> std::optional<Solved> {
        const auto deck = NiTCAD::deck::read(std::filesystem::path(NITCAD_SOURCE_DIR) / "examples/pn_diode_1d.deck");
        if (!deck) return std::nullopt;
        auto plan = NiTCAD::deck::build(*deck);
        if (!plan) return std::nullopt;
        plan->outputs = {};
        plan->points = {{0.0, 0.0}, {0.0, 0.3}};  // cathode, anode
        const auto report = NiTCAD::deck::run(*plan);
        if (!report || !report->sweep) return std::nullopt;
        Solved s{std::move(*plan), *report->sweep, {}};
        s.data = NiTCAD::verify::benchmark::device_data(s.plan.device);
        return s;
    }();
    REQUIRE(solved.has_value());
    REQUIRE(solved->sweep.points.size() == 2);
    return *solved;
}

const NiTCAD::mesh::Mesh& mesh() { return pn_diode().plan.device.mesh(); }

const c::CheckPart* part(const c::CheckResult& r, const char* name) {
    for (const auto& p : r.parts) {
        if (p.name == name) return &p;
    }
    return nullptr;
}

void print(const c::CheckResult& r) {
    std::printf("%s: %s %s\n", r.id.c_str(), std::string(v::to_string(r.comparison.outcome)).c_str(),
                r.comparison.reason.c_str());
    for (const auto& p : r.parts) {
        std::printf("  %s: %s residual %.6g tolerance %.6g ratio %.3g at %s %s\n", p.name.c_str(),
                    std::string(v::to_string(p.outcome)).c_str(), p.residual, p.tolerance,
                    p.tolerance > 0.0 ? p.residual / p.tolerance : 0.0, p.where.c_str(), p.reason.c_str());
    }
}

}  // namespace

TEST_CASE("device data copied for the consistency layer", "[verify][consistency]") {
    const auto& s = pn_diode();
    CHECK(c::validate(s.data, mesh()).empty());
    REQUIRE(s.data.contacts.size() == 2);
    CHECK(s.data.contacts[0].name == "cathode");
    CHECK(s.data.contacts[1].kind == c::ContactKind::ohmic);
}

TEST_CASE("C1 equilibrium on the pn diode", "[verify][consistency]") {
    const auto& s = pn_diode();
    const auto& eq = s.sweep.points[0];
    const auto r = c::check_c1(eq, mesh(), s.data);
    print(r);
    CHECK(r.comparison.outcome == Outcome::pass);
    REQUIRE(r.parts.size() == 3);

    SECTION("the Newton term carries the measured split (owner ruling)") {
        // At 0 V the split is V_T times the last Newton update (measured 1.4695e-14 eV against 64 eps scale
        // 1.4395e-14 eV, DECISIONS 11.13): without the term the part would fail.
        const double vt = NiTCAD::base::k_B_eV_per_K * s.data.temperature_K;
        const double u = eq.convergence.iterations.back().update;
        const auto* f = part(r, "quasi-Fermi levels");
        REQUIRE(f != nullptr);
        CHECK(f->residual > f->tolerance - vt * u);
    }
    SECTION("the equilibrium solve itself, held to rounding") {
        const auto e = NiTCAD::solve::solve_equilibrium(s.plan.device);
        REQUIRE(e.has_value());
        const auto q = c::check_c1(*e, mesh(), s.data);
        print(q);
        CHECK(q.comparison.outcome == Outcome::pass);
        auto bad = *e;
        bad.bands.hole_fermi_eV[50] -= 1e-13;
        const auto f = c::check_c1(bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "quasi-Fermi levels")->where == "node 50");
    }
    SECTION("no Newton record: the quasi-Fermi part NOT EVALUATED") {
        auto bad = eq;
        bad.convergence.iterations.clear();
        CHECK(c::check_c1(bad, mesh(), s.data).comparison.outcome == Outcome::not_evaluated);
    }
    SECTION("forward bias: NOT APPLICABLE, the predicate named") {
        const auto f = c::check_c1(s.sweep.points[1], mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::not_applicable);
        CHECK(f.comparison.reason == "all biases 0 = 0.3 (threshold 0)");
        CHECK(f.parts.empty());
    }
    SECTION("defect: a terminal current above its resolution") {
        auto bad = eq;
        bad.terminal_current[1] = 2.0 * bad.terminal_current_resolution[1];
        const auto f = c::check_c1(bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "terminal currents")->outcome == Outcome::fail);
        CHECK(part(f, "terminal currents")->where == "contact anode");
    }
    SECTION("defect: one edge current") {
        auto bad = eq;
        double tol = bad.terminal_current_resolution[0];
        for (double x : bad.terminal_current_resolution) tol = std::min(tol, x);
        bad.edge_current_p[150] = 1.01 * tol;
        const auto f = c::check_c1(bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "edge currents")->where == "hole edge 150");
    }
    SECTION("defect: F_n shifted on one node") {
        auto bad = eq;
        bad.bands.electron_fermi_eV[200] += 1e-12;
        const auto f = c::check_c1(bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "quasi-Fermi levels")->where == "node 200");
    }
    SECTION("missing bands: NOT EVALUATED") {
        auto bad = eq;
        bad.bands = {};
        const auto f = c::check_c1(bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::not_evaluated);
        CHECK(f.comparison.reason.find("band diagram not stored") != std::string::npos);
    }
    SECTION("thermal contacts at different temperatures: NOT APPLICABLE") {
        auto hot = eq;
        hot.thermal_bias_K = {300.0, 310.0};
        const auto f = c::check_c1(hot, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::not_applicable);
        CHECK(f.comparison.reason == "thermal contacts at one temperature = 10 (threshold 0)");
    }
}

TEST_CASE("C2 rates at equilibrium on the pn diode", "[verify][consistency]") {
    const auto& s = pn_diode();
    const auto& eq = s.sweep.points[0];
    const auto r = c::check_c2(eq, mesh(), s.data);
    print(r);
    CHECK(r.comparison.outcome == Outcome::pass);
    CHECK(part(r, "srh") != nullptr);

    SECTION("defect: one node recombining") {
        auto bad = eq;
        // q R V at one node equal to twice the smallest resolution.
        double tol = 0.0;
        for (double x : bad.terminal_current_resolution) {
            if (x > 0.0 && (tol == 0.0 || x < tol)) tol = x;
        }
        bad.rates.srh[100] = 2.0 * tol / (NiTCAD::base::q_C * mesh().volumes()[100]);
        const auto f = c::check_c2(bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "srh")->outcome == Outcome::fail);
    }
    SECTION("volume-weighted: half the resolution at one node passes") {
        auto ok = eq;
        double tol = 0.0;
        for (double x : ok.terminal_current_resolution) {
            if (x > 0.0 && (tol == 0.0 || x < tol)) tol = x;
        }
        ok.rates.srh[100] = 0.5 * tol / (NiTCAD::base::q_C * mesh().volumes()[100]);
        CHECK(c::check_c2(ok, mesh(), s.data).comparison.outcome == Outcome::pass);
    }
    SECTION("defect: generation and recombination that cancel in a signed sum") {
        auto bad = eq;
        double tol = 0.0;
        for (double x : bad.terminal_current_resolution) {
            if (x > 0.0 && (tol == 0.0 || x < tol)) tol = x;
        }
        REQUIRE(mesh().volumes()[100] == mesh().volumes()[101]);
        bad.rates.srh[100] = 2.0 * tol / (NiTCAD::base::q_C * mesh().volumes()[100]);
        bad.rates.srh[101] = -bad.rates.srh[100];
        CHECK(c::check_c2(bad, mesh(), s.data).comparison.outcome == Outcome::fail);
    }
    SECTION("no rates stored: NOT EVALUATED") {
        auto bad = eq;
        bad.rates = {};
        const auto f = c::check_c2(bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::not_evaluated);
        CHECK(f.comparison.reason.find("no rates stored") != std::string::npos);
    }
    SECTION("forward bias: NOT APPLICABLE") {
        CHECK(c::check_c2(s.sweep.points[1], mesh(), s.data).comparison.outcome == Outcome::not_applicable);
    }
}

TEST_CASE("C4 Kirchhoff on the pn diode", "[verify][consistency]") {
    const auto& s = pn_diode();
    for (const auto& point : s.sweep.points) {
        const auto r = c::check_c4(point, mesh(), s.data);
        print(r);
        CHECK(r.comparison.outcome == Outcome::pass);
    }
    SECTION("defect: one current corrupted") {
        auto bad = s.sweep.points[1];
        bad.terminal_current[0] += 3.0 * (bad.terminal_current_resolution[0] + bad.terminal_current_resolution[1]);
        const auto f = c::check_c4(bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "terminal currents")->outcome == Outcome::fail);
    }
    SECTION("defect: a NaN current fails") {
        auto bad = s.sweep.points[1];
        bad.terminal_current[1] = NAN;
        CHECK(c::check_c4(bad, mesh(), s.data).comparison.outcome == Outcome::fail);
    }
    SECTION("currents missing: NOT EVALUATED") {
        auto bad = s.sweep.points[1];
        bad.terminal_current.clear();
        CHECK(c::check_c4(bad, mesh(), s.data).comparison.outcome == Outcome::not_evaluated);
    }
}

TEST_CASE("C12 positivity on the pn diode", "[verify][consistency]") {
    const auto& s = pn_diode();
    for (const auto& point : s.sweep.points) {
        const auto r = c::check_c12(point.fields, &point, mesh(), s.data);
        print(r);
        CHECK(r.comparison.outcome == Outcome::pass);
    }
    const auto& pt = s.sweep.points[1];
    SECTION("defect: a negative density") {
        auto bad = pt;
        bad.fields.n_cm3[17] = -1.0;
        const auto f = c::check_c12(bad.fields, &bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "densities")->where == "electron node 17 n = -1");
    }
    SECTION("defect: a NaN density") {
        auto bad = pt;
        bad.fields.p_cm3[3] = NAN;
        CHECK(c::check_c12(bad.fields, &bad, mesh(), s.data).comparison.outcome == Outcome::fail);
    }
    SECTION("defect: a zero electron density") {
        auto bad = pt;
        bad.fields.n_cm3[4] = 0.0;
        CHECK(c::check_c12(bad.fields, &bad, mesh(), s.data).comparison.outcome == Outcome::fail);
    }
    SECTION("defect: a zero density") {
        auto bad = pt;
        bad.fields.p_cm3[3] = 0.0;
        CHECK(c::check_c12(bad.fields, &bad, mesh(), s.data).comparison.outcome == Outcome::fail);
    }
    SECTION("defect: a non-positive temperature") {
        auto bad = pt;
        bad.fields.temperature_K.assign(mesh().node_count(), 300.0);
        bad.fields.temperature_K[9] = 0.0;
        const auto f = c::check_c12(bad.fields, &bad, mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "temperatures")->where == "lattice node 9 T = 0");
    }
    SECTION("interface traps: occupancies NOT EVALUATED") {
        auto traps = s.data;
        traps.interface_traps = true;
        const auto f = c::check_c12(pt.fields, &pt, mesh(), traps);
        CHECK(f.comparison.outcome == Outcome::not_evaluated);
        CHECK(f.comparison.reason.find("occupancies not stored") != std::string::npos);
    }
}
