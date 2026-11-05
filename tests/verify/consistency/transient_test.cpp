// C6, transient charge bookkeeping (DECISIONS Unit 55 section 11.6 and the three C6 rulings of 11.13): passes on
// transients of examples/moscap_cv.deck run through the API (the deck is not edited) with BDF2 (a gate ramp and a
// jump, so backward-Euler restarts and variable-step BDF2) and with backward Euler; is NOT APPLICABLE on a pn diode
// transient (no gate); each injected defect fails it; each missing input makes it NOT EVALUATED. The trapezoidal
// diagnostic is reported, never judged.
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace c = NiTCAD::verify::consistency;
namespace v = NiTCAD::verify;
namespace solve = NiTCAD::solve;
using v::Outcome;

namespace {

struct Run {
    NiTCAD::deck::Plan plan;
    NiTCAD::results::Transient result;
    c::DeviceData data;
    const NiTCAD::mesh::Mesh& mesh() const { return plan.device.mesh(); }
};

std::filesystem::path source(const char* file) { return std::filesystem::path(NITCAD_SOURCE_DIR) / file; }

std::optional<Run> run(const char* file, std::vector<solve::Waveform> waveforms, solve::Integrator integrator,
                       double t_end) {
    auto deck = NiTCAD::deck::read(source(file));
    if (!deck) return std::nullopt;
    auto plan = NiTCAD::deck::build(*deck);
    if (!plan) return std::nullopt;
    solve::TransientOptions o;
    o.steady = plan->options;
    o.integrator = integrator;
    o.t_end_s = t_end;
    o.rtol = 1e-2;  // fewer steps (Debug cost, VR10); the check is per step, so the steps' number does not matter
    auto r = solve::solve_transient(plan->device, waveforms, o);
    if (!r) {
        std::printf("transient failed: %s\n", r.error().message.c_str());
        return std::nullopt;
    }
    if (r->stopped) std::printf("transient stopped: %s\n", r->stopped->message.c_str());
    Run out{std::move(*plan), std::move(*r), {}};
    out.data = NiTCAD::verify::benchmark::device_data(out.plan.device, out.plan.options.models);
    return out;
}

// Gate, body: the gate ramps from 0 to -0.5 V over 100 ps, holds, and jumps to 0.1 V at 200 ps (a breakpoint with
// two corners at one time); the body stays at 0.
std::vector<solve::Waveform> moscap_waveforms() {
    auto gate = solve::Waveform::piecewise_linear({{0.0, 0.0}, {1e-10, -0.5}, {2e-10, -0.5}, {2e-10, 0.1}});
    REQUIRE(gate.has_value());
    return {*gate, solve::Waveform::constant(0.0)};
}

const Run& moscap(solve::Integrator integrator) {
    static std::optional<Run> bdf2, be;
    auto& slot = integrator == solve::Integrator::bdf2 ? bdf2 : be;
    if (!slot) slot = run("examples/moscap_cv.deck", moscap_waveforms(), integrator, 4e-10);
    REQUIRE(slot.has_value());
    REQUIRE_FALSE(slot->result.stopped.has_value());
    REQUIRE(slot->result.points.size() > 10);
    return *slot;
}

void print(const char* label, const Run& r, const c::CheckResult& k) {
    std::size_t orders[3] = {0, 0, 0};
    std::size_t restarts = 0;
    const auto& pts = r.result.points;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        ++orders[pts[i].order];
        if (pts[i].order == 1 && pts[i - 1].order == 2) ++restarts;
    }
    std::printf("C6 %s: %zu points (BE %zu, BDF2 %zu, BE after BDF2 %zu), rejected %zu\n", label, pts.size(), orders[1],
                orders[2], restarts, r.result.rejected_steps);
    for (const auto& p : k.parts) {
        std::printf("  %s: residual %.3e tolerance %.3e at %s\n", p.name.c_str(), p.residual, p.tolerance,
                    p.where.c_str());
    }
    for (const auto& d : k.diagnostics) {
        std::printf("  [diagnostic] %s %s: %.3e (scale %.3e) %s\n", d.name.c_str(), d.where.c_str(), d.value, d.scale,
                    d.note.c_str());
    }
}

const c::CheckPart* part(const c::CheckResult& r, std::string_view name) {
    for (const auto& p : r.parts) {
        if (p.name == name) return &p;
    }
    return nullptr;
}

std::size_t first_bdf2(const NiTCAD::results::Transient& t) {
    for (std::size_t i = 2; i < t.points.size(); ++i) {
        if (t.points[i].order == 2) return i;
    }
    return 0;  // none: the callers REQUIRE a nonzero index
}

}  // namespace

TEST_CASE("C6 contract", "[verify][consistency]") {
    const auto info = c::check_info("C6");
    CHECK(v::validate(info).empty());
    CHECK(info.cls == v::ReferenceClass::D);
    CHECK(info.applicability.size() == 1);
}

TEST_CASE("C6 passes on moscap transients", "[verify][consistency]") {
    for (const auto integrator : {solve::Integrator::bdf2, solve::Integrator::backward_euler}) {
        const Run& r = moscap(integrator);
        const auto k = c::check_c6(r.result, r.mesh(), r.data);
        print(integrator == solve::Integrator::bdf2 ? "moscap BDF2" : "moscap BE", r, k);
        CHECK(k.comparison.outcome == Outcome::pass);
        REQUIRE(k.parts.size() == 2);
        CHECK(part(k, "displacement identity") != nullptr);
        CHECK(part(k, "gate and electrode currents") != nullptr);
        // Diagnostics: two per contact (cumulative identity, trapezoidal integral), finite, never judged.
        REQUIRE(k.diagnostics.size() == 2 * r.data.contacts.size());
        for (const auto& d : k.diagnostics) CHECK(std::isfinite(d.value));
    }
    // The BDF2 run exercises both formulas and a restart after the jump.
    const auto& pts = moscap(solve::Integrator::bdf2).result.points;
    bool restart = false;
    std::size_t bdf2 = 0;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        if (pts[i].order == 2) ++bdf2;
        if (pts[i].order == 1 && pts[i - 1].order == 2) restart = true;
    }
    CHECK(bdf2 > 3);
    CHECK(restart);
    for (const auto& p : moscap(solve::Integrator::backward_euler).result.points) CHECK(p.order <= 1);
}

TEST_CASE("C6 is NOT APPLICABLE without a gate", "[verify][consistency]") {
    auto anode = solve::Waveform::ramp(0.0, 0.2, 0.0, 1e-10);
    REQUIRE(anode.has_value());
    // examples/pn_diode_1d.deck: cathode, anode.
    static const auto r = run("examples/pn_diode_1d.deck", {solve::Waveform::constant(0.0), *anode},
                              solve::Integrator::bdf2, 1e-10);
    REQUIRE(r.has_value());
    const auto k = c::check_c6(r->result, r->mesh(), r->data);
    CHECK(k.comparison.outcome == Outcome::not_applicable);
    CHECK(k.comparison.reason.find("transient with a gate or electrode") != std::string::npos);
    CHECK(k.parts.empty());
}

TEST_CASE("C6 reports NOT EVALUATED for missing inputs", "[verify][consistency]") {
    const Run& r = moscap(solve::Integrator::bdf2);
    auto expect = [&](const NiTCAD::results::Transient& t, const char* what) {
        const auto k = c::check_c6(t, r.mesh(), r.data);
        CAPTURE(what, k.comparison.reason);
        CHECK(k.comparison.outcome == Outcome::not_evaluated);
        CHECK(k.comparison.reason.find(what) != std::string::npos);
    };
    {
        auto t = r.result;
        t.points.resize(1);
        expect(t, "fewer than two time points");
    }
    {
        auto t = r.result;
        t.points[3].contact_charge.pop_back();
        expect(t, "not one per contact at point 3");
    }
    {
        auto t = r.result;
        t.points[5].displacement_current.clear();
        expect(t, "not one per contact at point 5");
    }
    {
        auto t = r.result;
        t.points[1].order = 2;  // BDF2 needs two earlier points
        expect(t, "order 2 without its earlier steps at point 1");
    }
    {
        auto t = r.result;
        t.points[4].order = 3;
        expect(t, "order 3");
    }
    {
        auto t = r.result;
        t.points[2].step_s = 0.0;
        expect(t, "step not positive at point 2");
    }
}

TEST_CASE("C6 fails on injected defects", "[verify][consistency]") {
    const Run& r = moscap(solve::Integrator::bdf2);
    const std::size_t gate = 0, body = 1;
    const std::size_t b2 = first_bdf2(r.result);
    REQUIRE(b2 >= 2);
    auto fails = [&](const NiTCAD::results::Transient& t, const c::DeviceData& d, const char* what,
                     const char* part_name) {
        const auto k = c::check_c6(t, r.mesh(), d);
        const auto* p = part(k, part_name);
        CAPTURE(what, k.comparison.reason);
        CHECK(k.comparison.outcome == Outcome::fail);
        REQUIRE(p != nullptr);
        CHECK(p->outcome == Outcome::fail);
    };
    {   // the gate's displacement current at a BDF2 step, moved by twice the ruled bound (beta h dI = 2 x 64 eps x
        // the operands, at most 4 |Q|max of the three charges for a step ratio up to 2, plus beta h |I|; a relative change of I can sit below the
        // rounding of Q when the step's charge change is small against Q)
        auto t = r.result;
        const auto& p = t.points;
        const double q = std::max({std::abs(p[b2].contact_charge[gate]), std::abs(p[b2 - 1].contact_charge[gate]),
                                   std::abs(p[b2 - 2].contact_charge[gate])});
        t.points[b2].displacement_current[gate] += 2.0 * 64.0 * DBL_EPSILON * 4.0 * q / (0.5 * p[b2].step_s);
        t.points[b2].terminal_current[gate] = t.points[b2].displacement_current[gate];
        fails(t, r.data, "gate I_disp", "displacement identity");
    }
    {   // the body's displacement current at the first (backward Euler) step
        auto t = r.result;
        t.points[1].displacement_current[body] *= 1.0 + 1e-9;
        fails(t, r.data, "body I_disp", "displacement identity");
    }
    {   // a stored charge moved (it enters three steps' identities)
        auto t = r.result;
        t.points[b2].contact_charge[gate] *= 1.0 + 1e-9;
        fails(t, r.data, "charge", "displacement identity");
    }
    {   // a BDF2 step relabelled backward Euler
        auto t = r.result;
        t.points[b2].order = 1;
        fails(t, r.data, "order", "displacement identity");
    }
    {   // the previous step's length off by 1e-6 (the BDF2 step ratio)
        auto t = r.result;
        t.points[b2 - 1].step_s *= 1.0 + 1e-6;
        fails(t, r.data, "step ratio", "displacement identity");
    }
    {   // a NaN charge
        auto t = r.result;
        t.points[b2].contact_charge[body] = std::nan("");
        fails(t, r.data, "NaN charge", "displacement identity");
    }
    {   // conduction current on the gate
        auto t = r.result;
        t.points[b2].conduction_current[gate] = 1e-300;
        fails(t, r.data, "gate conduction", "gate and electrode currents");
    }
    {   // the gate's terminal current not its displacement current
        auto t = r.result;
        t.points[0].terminal_current[gate] = std::nextafter(t.points[0].displacement_current[gate], 1.0);
        fails(t, r.data, "gate terminal", "gate and electrode currents");
    }
    {   // a tunnelling gate carries conduction current: not judged by the gate part
        auto t = r.result;
        auto d = r.data;
        d.contacts[gate].tunnelling = true;
        t.points[b2].conduction_current[gate] = 1e-3;
        const auto k = c::check_c6(t, r.mesh(), d);
        CHECK(k.comparison.outcome == Outcome::pass);
        CHECK(part(k, "gate and electrode currents") == nullptr);
    }
}

TEST_CASE("C6 diagnostics are reported, not judged", "[verify][consistency]") {
    // A trapezoidal integral far off the charge change leaves the outcome a pass.
    const Run& r = moscap(solve::Integrator::bdf2);
    auto t = r.result;
    for (auto& s : t.run.settings) {
        if (s.first == "transient.rtol") s.second = 1e-30;  // the scale only
    }
    const auto k = c::check_c6(t, r.mesh(), r.data);
    CHECK(k.comparison.outcome == Outcome::pass);
    bool trapezoid = false;
    for (const auto& d : k.diagnostics) {
        if (d.name == "trapezoidal integral") {
            trapezoid = true;
            CHECK(d.note.find("estimated steps") != std::string::npos);
        }
    }
    CHECK(trapezoid);
    // Without rtol in the run record the scale is NaN and the note says so.
    t.run.settings.clear();
    const auto n = c::check_c6(t, r.mesh(), r.data);
    for (const auto& d : n.diagnostics) {
        if (d.name == "trapezoidal integral") {
            CHECK(std::isnan(d.scale));
            CHECK(d.note.find("rtol not in the run record") != std::string::npos);
        }
    }
}

TEST_CASE("C6 closes the mutation gaps", "[verify][consistency]") {
    const Run& r = moscap(solve::Integrator::bdf2);
    {   // a meshed electrode satisfies the predicate and is judged as a gate
        auto d = r.data;
        d.contacts[0].kind = c::ContactKind::electrode;
        const auto k = c::check_c6(r.result, r.mesh(), d);
        CHECK(k.comparison.outcome == Outcome::pass);
        CHECK(part(k, "gate and electrode currents") != nullptr);
        auto t = r.result;
        t.points[7].conduction_current[0] = 1e-300;
        CHECK(c::check_c6(t, r.mesh(), d).comparison.outcome == Outcome::fail);
    }
    {   // a defect at the last point
        auto t = r.result;
        t.points.back().contact_charge[1] *= 1.0 + 1e-9;
        CHECK(c::check_c6(t, r.mesh(), r.data).comparison.outcome == Outcome::fail);
    }
    {   // the bound includes beta h |I|: two points, Q 0 -> 1, h = 1 s (backward Euler), I = 1 + 100 eps; the
        // residual 100 eps lies between 64 eps |Q_1| and 64 eps (|Q_1| + h |I|)
        NiTCAD::results::Transient t;
        for (int k = 0; k < 2; ++k) {
            NiTCAD::results::TimePoint p{static_cast<double>(k), static_cast<double>(k), k, std::nan(""), {0.0, 0.0}};
            p.contact_charge = {static_cast<double>(k), 0.0};
            p.displacement_current = {k == 0 ? 0.0 : 1.0 + 100.0 * DBL_EPSILON, 0.0};
            p.conduction_current = {0.0, 0.0};
            p.terminal_current = p.displacement_current;
            t.points.push_back(p);
        }
        const auto k = c::check_c6(t, r.mesh(), r.data);
        CHECK(k.comparison.outcome == Outcome::pass);
        const auto* p = part(k, "displacement identity");
        REQUIRE(p != nullptr);
        CHECK(p->tolerance == 64.0 * DBL_EPSILON * (1.0 + 1.0 + 100.0 * DBL_EPSILON));
    }
    {   // the diagnostic scale is proportional to the run's rtol
        auto t = r.result;
        double rtol = 0.0;
        for (auto& s : t.run.settings) {
            if (s.first == "transient.rtol") rtol = s.second;
        }
        REQUIRE(rtol > 0.0);
        const auto base = c::check_c6(t, r.mesh(), r.data);
        for (auto& s : t.run.settings) {
            if (s.first == "transient.rtol") s.second = 2.0 * rtol;
        }
        const auto twice = c::check_c6(t, r.mesh(), r.data);
        REQUIRE(base.diagnostics.size() == twice.diagnostics.size());
        for (std::size_t i = 0; i < base.diagnostics.size(); ++i) {
            if (base.diagnostics[i].name != "trapezoidal integral") continue;
            CHECK(base.diagnostics[i].scale > 0.0);
            CHECK(twice.diagnostics[i].scale == 2.0 * base.diagnostics[i].scale);
        }
    }
}
