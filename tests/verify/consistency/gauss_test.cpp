// C3, the dual-cell Gauss law (DECISIONS Unit 55 section 11.6 and the two C3 rulings of 11.13): passes on real
// solves of examples/pn_diode_1d.deck (ohmic contacts skipped), examples/moscap_cv.deck (a lumped gate as a group,
// at DC points set through the API), decks/schottky_diode_1d.deck (a Schottky contact skipped) and a 1D
// silicon-oxide device written here (a meshed electrode and an interface with fixed charge and traps, one group);
// each injected defect fails it at the named place, and each missing input makes it NOT EVALUATED.
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/base/constants.hpp"
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
                            std::optional<std::vector<std::vector<double>>> points) {
    if (!deck) return std::nullopt;
    auto plan = NiTCAD::deck::build(*deck);
    if (!plan) return std::nullopt;
    plan->outputs = {};
    if (points) {
        plan->analysis = NiTCAD::deck::Analysis::sweep;
        plan->points = *points;
    }
    auto report = NiTCAD::deck::run(*plan);
    if (!report) return std::nullopt;
    Solved s{std::move(*plan), std::move(*report), {}};
    s.data = NiTCAD::verify::benchmark::device_data(s.plan.device, s.plan.options.models);
    return s;
}

std::filesystem::path source(const char* file) { return std::filesystem::path(NITCAD_SOURCE_DIR) / file; }

const Solved& pn_diode() {
    static const auto s = solve(NiTCAD::deck::read(source("examples/pn_diode_1d.deck")),
                                std::vector<std::vector<double>>{{0.0, 0.0}, {0.0, 0.3}});
    REQUIRE(s.has_value());
    REQUIRE(s->report.sweep->points.size() == 2);
    return *s;
}

const Solved& moscap() {
    // gate, body: accumulation and depletion (the deck's own analysis is small-signal; DC points replace it).
    static const auto s = solve(NiTCAD::deck::read(source("examples/moscap_cv.deck")),
                                std::vector<std::vector<double>>{{0.0, 0.0}, {-0.5, 0.0}, {-1.0, 0.0}, {0.25, 0.0}});
    REQUIRE(s.has_value());
    REQUIRE(s->report.sweep->points.size() == 4);
    return *s;
}

const Solved& schottky() {
    static const auto s = solve(NiTCAD::deck::read(source("decks/schottky_diode_1d.deck")),
                                std::vector<std::vector<double>>{{0.0, 0.0}, {-0.5, 0.0}, {0.2, 0.0}});
    REQUIRE(s.has_value());
    return *s;
}

// A 1D silicon / oxide stack: the interface carries a fixed charge, a trap level and a trap band, so its two edge
// ends form one group with the meshed electrode behind the oxide.
constexpr const char* oxide_deck = R"(
[mesh]
type = tensor
x = linspace(0, 1um, 101) + linspace(1um, 1.02um, 11)
straddle x = 1um 1nm

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 1um 1.02um

[doping]
profile = acceptor constant 1e17

[contact body]
kind = ohmic
patch = x_min

[contact gate]
kind = electrode
patch = x_max
electrode = metal
work_function = 4.6

[interface silicon oxide]
fixed_charge = 5e11
trap_level = acceptor 1e11 0.2
trap_band = donor 1e11 -0.4 -0.1 1e-15 2e-15

[sweep]
contact gate = 0 : 0.5 : 0.25
)";

const Solved& oxide() {
    static const auto s = solve(NiTCAD::deck::parse(oxide_deck),
                                std::vector<std::vector<double>>{{0.0, 0.0}, {0.0, -0.05}, {0.0, 0.05}});  // body, gate
    REQUIRE(s.has_value());
    REQUIRE(s->report.sweep.has_value());
    if (s->report.sweep->stopped) std::printf("oxide sweep stopped: %s\n", s->report.sweep->stopped->message.c_str());
    REQUIRE(s->report.sweep->points.size() == 3);
    return *s;
}

const c::CheckPart* part(const c::CheckResult& r, std::string_view name) {
    for (const auto& p : r.parts) {
        if (p.name.starts_with(name)) return &p;
    }
    return nullptr;
}

void print(const char* what, const c::CheckResult& r) {
    std::printf("%s %s: %s %s\n", what, r.id.c_str(), std::string(v::to_string(r.comparison.outcome)).c_str(),
                r.comparison.reason.c_str());
    for (const auto& p : r.parts) {
        std::printf("  %s: %s residual %.6g tolerance %.6g ratio %.3g at %s %s\n", p.name.c_str(),
                    std::string(v::to_string(p.outcome)).c_str(), p.residual, p.tolerance,
                    p.tolerance > 0.0 ? p.residual / p.tolerance : 0.0, p.where.c_str(), p.reason.c_str());
    }
}

}  // namespace

TEST_CASE("C3 on the pn diode: cells pass, ohmic contacts skipped", "[verify][consistency]") {
    const auto& s = pn_diode();
    for (std::size_t k = 0; k < 2; ++k) {
        const auto r = c::check_c3(s.point(k), s.mesh(), s.data);
        print(k == 0 ? "pn 0 V" : "pn 0.3 V", r);
        CHECK(r.comparison.outcome == Outcome::pass);
        REQUIRE(part(r, "cells") != nullptr);
        CHECK(part(r, "groups") == nullptr);  // no gate, electrode or charged interface
        REQUIRE(part(r, "global") != nullptr);
        CHECK(part(r, "global")->where.find("cathode") != std::string::npos);
    }
    SECTION("the equilibrium solve") {
        const auto eq = NiTCAD::solve::solve_equilibrium(s.plan.device);
        REQUIRE(eq.has_value());
        const auto r = c::check_c3(*eq, s.mesh(), s.data);
        print("pn equilibrium", r);
        CHECK(r.comparison.outcome == Outcome::pass);
    }
    SECTION("a potential shifted at one node fails at that node") {
        auto bad = s.point(1);
        bad.fields.potential_V[200] += 1e-9;
        const auto r = c::check_c3(bad, s.mesh(), s.data);
        CHECK(r.comparison.outcome == Outcome::fail);
        CHECK(part(r, "cells")->where == "node 200");
    }
    SECTION("a density changed at one node fails at that node") {
        auto bad = s.point(1);
        bad.fields.n_cm3[150] *= 1.0 + 1e-6;
        bad.fields.p_cm3[150] *= 1.0 - 1e-6;
        const auto r = c::check_c3(bad, s.mesh(), s.data);
        CHECK(r.comparison.outcome == Outcome::fail);
        CHECK(part(r, "cells")->where == "node 150");
    }
    SECTION("a contact cell's corruption is not judged (its charge is not stored)") {
        auto bad = s.point(1);
        bad.fields.n_cm3[0] *= 2.0;
        const auto r = c::check_c3(bad, s.mesh(), s.data);
        CHECK(r.comparison.outcome == Outcome::pass);
    }
    SECTION("not converged: NOT APPLICABLE, the predicate named") {
        auto bad = s.point(1);
        bad.convergence.converged = false;
        const auto r = c::check_c3(bad, s.mesh(), s.data);
        CHECK(r.comparison.outcome == Outcome::not_applicable);
        CHECK(r.comparison.reason == "converged state = 0 (threshold 1)");
    }
    SECTION("no Newton record, missing fields, incomplete ionization: NOT EVALUATED") {
        auto bad = s.point(1);
        bad.convergence.iterations.clear();
        CHECK(c::check_c3(bad, s.mesh(), s.data).comparison.outcome == Outcome::not_evaluated);
        auto empty = s.point(1);
        empty.fields.n_cm3.clear();
        CHECK(c::check_c3(empty, s.mesh(), s.data).comparison.outcome == Outcome::not_evaluated);
        auto ii = s.data;
        ii.incomplete_ionization = true;
        const auto r = c::check_c3(s.point(1), s.mesh(), ii);
        CHECK(r.comparison.outcome == Outcome::not_evaluated);
        CHECK(r.comparison.reason.find("incomplete ionization") != std::string::npos);
    }
}

TEST_CASE("C3 on the MOS capacitor: the lumped gate as a group", "[verify][consistency]") {
    const auto& s = moscap();
    CHECK(c::validate(s.data, s.mesh()).empty());
    for (std::size_t k = 0; k < 4; ++k) {
        const auto r = c::check_c3(s.point(k), s.mesh(), s.data);
        print(("moscap point " + std::to_string(k)).c_str(), r);
        CHECK(r.comparison.outcome == Outcome::pass);
        REQUIRE(part(r, "groups") != nullptr);
        CHECK(part(r, "groups")->where.find("gate gate") != std::string::npos);
    }
    SECTION("a dropped gate charge fails the group") {
        auto bad = s.point(2);
        bad.gate_charge[0] = 0.0;
        const auto r = c::check_c3(bad, s.mesh(), s.data);
        CHECK(r.comparison.outcome == Outcome::fail);
        CHECK(part(r, "groups")->outcome == Outcome::fail);
    }
    SECTION("a dropped fixed oxide charge fails the group") {
        auto d = s.data;
        d.contacts[0].gate_fixed_charge_cm2 = 1e10;  // the deck has none: an extra one is as wrong
        const auto r = c::check_c3(s.point(2), s.mesh(), d);
        CHECK(r.comparison.outcome == Outcome::fail);
    }
    SECTION("gate charges not stored: NOT EVALUATED") {
        auto bad = s.point(2);
        bad.gate_charge.clear();
        const auto r = c::check_c3(bad, s.mesh(), s.data);
        CHECK(r.comparison.outcome == Outcome::not_evaluated);
        CHECK(r.comparison.reason.find("gate charges not stored") != std::string::npos);
    }
}

TEST_CASE("C3 on the Schottky diode: the Schottky contact skipped", "[verify][consistency]") {
    const auto& s = schottky();
    for (std::size_t k = 0; k < s.report.sweep->points.size(); ++k) {
        const auto r = c::check_c3(s.point(k), s.mesh(), s.data);
        print(("schottky point " + std::to_string(k)).c_str(), r);
        CHECK(r.comparison.outcome == Outcome::pass);
        CHECK(part(r, "global")->where.find("anode") != std::string::npos);
    }
}

TEST_CASE("C3 on a charged silicon-oxide interface: one group with the electrode", "[verify][consistency]") {
    const auto& s = oxide();
    REQUIRE(c::validate(s.data, s.mesh()).empty());
    REQUIRE(s.data.interfaces.size() == 1);
    CHECK(s.data.interfaces[0].charged);
    CHECK(s.data.interfaces[0].trap_density_cm2 > 1e11);
    const auto& points = s.report.sweep->points;
    REQUIRE(!points.empty());
    for (std::size_t k = 0; k < points.size(); ++k) {
        const auto r = c::check_c3(points[k], s.mesh(), s.data);
        print(("oxide point " + std::to_string(k)).c_str(), r);
        CHECK(r.comparison.outcome == Outcome::pass);
        REQUIRE(part(r, "groups") != nullptr);
        CHECK(part(r, "groups")->where.find("interface 0") != std::string::npos);
    }
    SECTION("the equilibrium solve") {
        const auto eq = NiTCAD::solve::solve_equilibrium(s.plan.device);
        REQUIRE(eq.has_value());
        const auto r = c::check_c3(*eq, s.mesh(), s.data);
        print("oxide equilibrium", r);
        CHECK(r.comparison.outcome == Outcome::pass);
    }
    SECTION("a dropped interface trap charge fails the group") {
        auto bad = points.back();
        REQUIRE(bad.interface_trap_charge.size() == 1);
        bad.interface_trap_charge[0] = 0.0;
        const auto r = c::check_c3(bad, s.mesh(), s.data);
        CHECK(r.comparison.outcome == Outcome::fail);
        CHECK(part(r, "groups")->outcome == Outcome::fail);
    }
    SECTION("a dropped fixed interface charge fails the group") {
        auto d = s.data;
        d.interfaces[0].fixed_charge_cm2 = 0.0;
        const auto r = c::check_c3(points.back(), s.mesh(), d);
        CHECK(r.comparison.outcome == Outcome::fail);
        CHECK(part(r, "groups")->outcome == Outcome::fail);
    }
    SECTION("a dropped electrode charge fails the group") {
        auto bad = points.back();
        bad.gate_charge[1] = 0.0;
        CHECK(c::check_c3(bad, s.mesh(), s.data).comparison.outcome == Outcome::fail);
    }
    SECTION("trap charges not stored: NOT EVALUATED") {
        auto bad = points.back();
        bad.interface_trap_charge.clear();
        const auto r = c::check_c3(bad, s.mesh(), s.data);
        CHECK(r.comparison.outcome == Outcome::not_evaluated);
        CHECK(r.comparison.reason.find("interface trap charges not stored") != std::string::npos);
    }
}

namespace {

// The oxide stack without the interface model: its silicon-oxide edge is an ordinary edge, whose permittivity is
// the harmonic mean of its ends' (the discretization's).
std::string uncharged_oxide_deck() {
    std::string s = oxide_deck;
    return s.erase(s.find("[interface")) + "[equilibrium]\n";
}

// The MOS capacitor with a fixed oxide charge under its lumped gate.
constexpr const char* charged_moscap_deck = R"(
[mesh]
type = tensor
x = linspace(0, 20nm, 41) + linspace(20nm, 0.5um, 97) + linspace(0.5um, 2um, 31)
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e17
[contact gate]
kind = gate
patch = x_min
oxide_thickness = 10nm
electrode = n_poly
fixed_charge = 3e11
[contact body]
kind = ohmic
patch = x_max
[equilibrium]
)";

// A 2D silicon / oxide device whose interface has many charged edges (one group across them), at equilibrium.
constexpr const char* oxide_2d_deck = R"(
[mesh]
type = tensor
x = linspace(0, 1um, 11)
y = linspace(0, 0.3um, 7) + 0.35um, 0.4um
straddle y = 0.3um 10nm
[region silicon]
material = silicon
[region oxide]
material = sio2
box = 0 1um   0.3um 0.4um
[doping]
profile = acceptor constant 1e17
profile = donor gaussian 1e19 axis=y center=0 width=50nm x=0:0.3um lateral=20nm
[contact substrate]
kind = ohmic
patch = y_min
box = 0.2um 0.8um  0 0
[contact gate]
kind = electrode
patch = y_max
electrode = metal
work_function = 4.6
[interface silicon oxide]
fixed_charge = 5e10
trap_level = acceptor 1e10 0.2
trap_band = donor 1e11 -0.4 -0.1 1e-15 2e-15
[equilibrium]
)";

const Solved& equilibrium_of(std::string_view text) {
    static std::vector<std::pair<std::string, std::optional<Solved>>> cache;
    for (const auto& [t, s] : cache) {
        if (t == text) return *s;
    }
    cache.emplace_back(std::string(text), solve(NiTCAD::deck::parse(text), std::nullopt));
    REQUIRE(cache.back().second.has_value());
    REQUIRE(cache.back().second->report.equilibrium.has_value());
    return *cache.back().second;
}

}  // namespace

TEST_CASE("C3 on an uncharged oxide edge: the harmonic-mean permittivity", "[verify][consistency]") {
    const std::string text = uncharged_oxide_deck();
    const auto& s = equilibrium_of(text);
    REQUIRE(s.data.interfaces.empty());
    const auto r = c::check_c3(*s.report.equilibrium, s.mesh(), s.data);
    print("uncharged oxide equilibrium", r);
    CHECK(r.comparison.outcome == Outcome::pass);
    REQUIRE(part(r, "groups") != nullptr);  // the electrode alone
    CHECK(part(r, "groups")->where.find("electrode gate") != std::string::npos);
}

TEST_CASE("C3 on a lumped gate with fixed oxide charge", "[verify][consistency]") {
    const auto& s = equilibrium_of(charged_moscap_deck);
    REQUIRE(s.data.contacts[0].gate_fixed_charge_cm2 == 3e11);
    const auto r = c::check_c3(*s.report.equilibrium, s.mesh(), s.data);
    print("charged moscap equilibrium", r);
    CHECK(r.comparison.outcome == Outcome::pass);
    auto d = s.data;
    d.contacts[0].gate_fixed_charge_cm2 = 0.0;  // dropped
    const auto f = c::check_c3(*s.report.equilibrium, s.mesh(), d);
    CHECK(f.comparison.outcome == Outcome::fail);
    CHECK(part(f, "groups")->outcome == Outcome::fail);
}

TEST_CASE("C3 on a 2D charged interface: its edges as one group", "[verify][consistency]") {
    const auto& s = equilibrium_of(oxide_2d_deck);
    REQUIRE(c::validate(s.data, s.mesh()).empty());
    const auto r = c::check_c3(*s.report.equilibrium, s.mesh(), s.data);
    print("2D oxide equilibrium", r);
    CHECK(r.comparison.outcome == Outcome::pass);
    SECTION("a dropped trap charge fails the group") {
        auto bad = *s.report.equilibrium;
        bad.interface_trap_charge[0] = 0.0;
        CHECK(c::check_c3(bad, s.mesh(), s.data).comparison.outcome == Outcome::fail);
    }
    SECTION("a potential shifted at an interface node fails the group") {
        auto bad = *s.report.equilibrium;
        const auto edges = s.mesh().edges();
        std::size_t node = 0;
        for (const auto& e : edges) {
            const auto a = static_cast<std::size_t>(e.first), b = static_cast<std::size_t>(e.second);
            if (s.data.insulator[a] != s.data.insulator[b]) {
                node = s.data.insulator[a] != 0 ? b : a;  // the semiconductor end
                if (s.data.donors[node] < 1e18) break;      // away from the source's doping
            }
        }
        bad.fields.potential_V[node] += 1e-9;
        const auto f = c::check_c3(bad, s.mesh(), s.data);
        CHECK(f.comparison.outcome == Outcome::fail);
        CHECK(part(f, "groups")->outcome == Outcome::fail);
    }
}

namespace {

// The ruled tolerance of 11.13, restated here independently of checks.cpp, for a node set G: 64 eps s + N u with
// s = sum over G's nodes of their faces' C (|psi_i| + |psi_j|) (charged edges also 2 eps_max A/d (|psi_a| + |psi_b|))
// + q V (n + p + N_D + N_A), and N = sum over G's crossing faces of 2 C V_T + q V (n + p); `extra_scale` and
// `extra_newton` carry a group's stored charges and their Newton terms.
struct Tol {
    double scale = 0.0, newton = 0.0;
};

Tol ruled(const NiTCAD::results::NodeFields& f, const NiTCAD::mesh::Mesh& m, const c::DeviceData& d,
          const std::vector<std::size_t>& G, bool charged_edges_in_G) {
    const double q = NiTCAD::base::q_C, e0 = NiTCAD::base::eps0_F_per_cm;
    const double vt = NiTCAD::base::k_B_eV_per_K * d.temperature_K;
    auto in = [&G](std::size_t i) { return std::find(G.begin(), G.end(), i) != G.end(); };
    Tol t;
    for (std::size_t i : G) {
        if (d.insulator[i] != 0) continue;
        const double V = m.volumes()[i];
        t.scale += q * V * (f.n_cm3[i] + f.p_cm3[i] + d.donors[i] + d.acceptors[i]);
        t.newton += q * V * (f.n_cm3[i] + f.p_cm3[i]);
    }
    for (const auto& e : m.edges()) {
        const auto a = static_cast<std::size_t>(e.first), b = static_cast<std::size_t>(e.second);
        if (!in(a) && !in(b)) continue;
        const double ea = d.relative_permittivity[a], eb = d.relative_permittivity[b];
        const double C = e0 * (2.0 * ea * eb / (ea + eb)) * e.coupling_area / e.length;
        double mag = C * (std::abs(f.potential_V[a]) + std::abs(f.potential_V[b]));
        if (charged_edges_in_G && in(a) && in(b) && d.insulator[a] != d.insulator[b]) {
            mag += 2.0 * e0 * std::max(ea, eb) * e.coupling_area / e.length *
                   (std::abs(f.potential_V[a]) + std::abs(f.potential_V[b]));
        }
        t.scale += (in(a) ? mag : 0.0) + (in(b) ? mag : 0.0);
        if (in(a) != in(b)) t.newton += 2.0 * C * vt;
    }
    return t;
}

double tolerance(const Tol& t, double u) { return 64.0 * DBL_EPSILON * t.scale + t.newton * u; }

bool close(double a, double b) { return std::abs(a - b) <= 1e-9 * std::abs(b); }

}  // namespace

TEST_CASE("C3 tolerances are the ruled ones (11.13)", "[verify][consistency]") {
    constexpr double u = 1e-6;  // a last update large enough that the Newton term shows
    SECTION("a cell of the pn diode") {
        const auto& s = pn_diode();
        auto x = s.point(1);
        x.convergence.iterations.back().update = u;
        x.fields.potential_V[150] += 1.0;  // made the worst cell; its tolerance is then reported
        const auto r = c::check_c3(x, s.mesh(), s.data);
        const auto* p = part(r, "cells");
        REQUIRE(p->where == "node 150");
        CHECK(close(p->tolerance, tolerance(ruled(x.fields, s.mesh(), s.data, {150}, false), u)));
        x.convergence.iterations.back().update = 0.0;  // the rounding floor alone
        const auto r0 = c::check_c3(x, s.mesh(), s.data);
        const auto* p0 = part(r0, "cells");
        CHECK(close(p0->tolerance, tolerance(ruled(x.fields, s.mesh(), s.data, {150}, false), 0.0)));
    }
    SECTION("the moscap's lumped gate: its oxide's Newton term C_ox A V_T") {
        const auto& s = moscap();
        auto x = s.point(2);
        x.convergence.iterations.back().update = u;
        x.gate_charge[0] *= 2.0;
        const auto& g = s.data.contacts[0];
        REQUIRE(g.nodes.size() == 1);
        const auto i = static_cast<std::size_t>(g.nodes[0]);
        Tol t = ruled(x.fields, s.mesh(), s.data, {i}, false);
        const double cox = g.oxide_capacitance_F_cm2 * g.gate_areas[0];
        t.scale += std::abs(x.gate_charge[0]) + 2.0 * cox * std::abs(x.fields.potential_V[i]);
        t.newton += cox * NiTCAD::base::k_B_eV_per_K * s.data.temperature_K;
        const auto rg = c::check_c3(x, s.mesh(), s.data);
        const auto* p = part(rg, "groups");
        CHECK(close(p->tolerance, tolerance(t, u)));
    }
    SECTION("the oxide's interface: its traps' Newton term 2 q N_t A") {
        const auto& s = oxide();
        auto x = s.point(0);
        x.convergence.iterations.back().update = u;
        x.interface_trap_charge[0] += 1e-8;  // C/cm^2, a fifth of the fixed charge
        std::vector<std::size_t> G;
        double area = 0.0;
        for (const auto& e : s.mesh().edges()) {
            const auto a = static_cast<std::size_t>(e.first), b = static_cast<std::size_t>(e.second);
            if (s.data.insulator[a] != s.data.insulator[b]) {
                G = {a, b};
                area += e.coupling_area;
            }
        }
        REQUIRE(G.size() == 2);
        Tol t = ruled(x.fields, s.mesh(), s.data, G, true);
        const auto& fi = s.data.interfaces[0];
        const double q = NiTCAD::base::q_C, Nt = q * fi.trap_density_cm2 * area;
        t.scale += std::abs(q * fi.fixed_charge_cm2 * area) + std::abs(x.interface_trap_charge[0]) + Nt;
        t.newton += 2.0 * Nt;
        const auto rr = c::check_c3(x, s.mesh(), s.data);
        print("oxide tolerance", rr);
        const auto* p = part(rr, "groups");
        REQUIRE(p->where.find("interface 0") != std::string::npos);
        CHECK(close(p->tolerance, tolerance(t, u)));
    }
}
