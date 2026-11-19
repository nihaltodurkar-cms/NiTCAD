// compare_sweeps (DECISIONS.md Unit 60, XG1 and XG2).
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "NiTCAD/analysis/compare.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/results/run.hpp"

namespace {

using NiTCAD::analysis::CompareOutcome;
using NiTCAD::analysis::CompareTolerances;
using NiTCAD::analysis::compare_sweeps;
using NiTCAD::analysis::QuantityComparison;
using NiTCAD::analysis::SweepComparison;
using NiTCAD::results::BiasPoint;
using NiTCAD::results::Sweep;

// Two contacts and a gate, three points, four nodes, every optional quantity present.
Sweep make_sweep() {
    Sweep s;
    s.run.input_identity = 0x1234;
    s.run.mesh_digest = 0x5678;
    for (std::size_t k = 0; k < 3; ++k) {
        BiasPoint p;
        const double v = 0.1 * static_cast<double>(k);
        p.bias_V = {0.0, v, 1.0};
        p.terminal_current = {-1e-3 * (1.0 + v), 1e-3 * (1.0 + v), 0.0};
        p.terminal_current_resolution = {1e-12, 1e-12, 0.0};
        p.weighted_current = p.terminal_current;
        p.weighted_current_resolution = {1e-14, 1e-14, 0.0};
        p.gate_charge = {0.0, 0.0, 2e-7};
        p.interface_trap_charge = {-3e-8};
        p.fields.potential_V = {0.0, 0.2 + v, 0.4, 0.6};
        p.fields.n_cm3 = {1e18, 1e10, 1e4, 1e2};
        p.fields.p_cm3 = {1e2, 1e10, 1e16, 1e16};
        p.fields.temperature_K = {300.0, 301.0, 302.0, 300.5};
        p.fields.electron_temperature_K = {300.0, 350.0, 400.0, 300.0};
        p.fields.hole_temperature_K = {300.0, 310.0, 320.0, 300.0};
        s.points.push_back(std::move(p));
    }
    return s;
}

const QuantityComparison* find(const SweepComparison& c, const std::string& name) {
    for (const auto& q : c.quantities)
        if (q.name == name) return &q;
    return nullptr;
}

// One quantity's value at point 1, index 2 (or the last index the quantity has), and its allowance
// under the default tolerances.
struct Probe {
    std::string name;
    std::function<double&(BiasPoint&)> value;
    std::function<double(const BiasPoint&)> allowance;
    std::size_t index;
};

std::vector<Probe> probes() {
    const CompareTolerances t;
    return {
        {"terminal_current", [](BiasPoint& p) -> double& { return p.terminal_current[1]; },
         [t](const BiasPoint& p) { return t.relative * p.terminal_current[1] + p.terminal_current_resolution[1]; }, 1},
        {"weighted_current", [](BiasPoint& p) -> double& { return p.weighted_current[0]; },
         [t](const BiasPoint& p) {
             return t.relative * -p.weighted_current[0] + p.weighted_current_resolution[0];
         }, 0},
        {"gate_charge", [](BiasPoint& p) -> double& { return p.gate_charge[2]; },
         [t](const BiasPoint& p) { return t.relative * p.gate_charge[2] + t.charge_floor; }, 2},
        {"interface_trap_charge", [](BiasPoint& p) -> double& { return p.interface_trap_charge[0]; },
         [t](const BiasPoint& p) { return t.relative * -p.interface_trap_charge[0] + t.charge_floor; }, 0},
        {"potential_V", [](BiasPoint& p) -> double& { return p.fields.potential_V[2]; },
         [t](const BiasPoint&) { return t.potential_V; }, 2},
        {"n_cm3", [](BiasPoint& p) -> double& { return p.fields.n_cm3[0]; },
         [t](const BiasPoint& p) { return t.relative * p.fields.n_cm3[0] + t.density_floor * 1e18; }, 0},
        {"p_cm3", [](BiasPoint& p) -> double& { return p.fields.p_cm3[2]; },
         [t](const BiasPoint& p) { return t.relative * p.fields.p_cm3[2] + t.density_floor * 1e18; }, 2},
        {"temperature_K", [](BiasPoint& p) -> double& { return p.fields.temperature_K[2]; },
         [t](const BiasPoint&) { return t.temperature_K; }, 2},
        {"electron_temperature_K", [](BiasPoint& p) -> double& { return p.fields.electron_temperature_K[2]; },
         [t](const BiasPoint&) { return t.temperature_K; }, 2},
        {"hole_temperature_K", [](BiasPoint& p) -> double& { return p.fields.hole_temperature_K[2]; },
         [t](const BiasPoint&) { return t.temperature_K; }, 2},
    };
}

}  // namespace

TEST_CASE("identical sweeps match with zero differences (XG1)", "[analysis][compare]") {
    const Sweep a = make_sweep();
    const SweepComparison c = compare_sweeps(a, a);
    CHECK(c.outcome == CompareOutcome::match);
    CHECK(c.reason.empty());
    CHECK(c.identity_equal);
    CHECK(c.mesh_digest_equal);
    REQUIRE(c.quantities.size() == 10);
    for (const auto& q : c.quantities) {
        INFO(q.name);
        CHECK(q.max_difference == 0.0);
        CHECK(q.max_ratio == 0.0);
        CHECK(q.compared > 0);
    }
}

TEST_CASE("each quantity passes at 0.9x its allowance and fails at 1.1x (XG1)", "[analysis][compare]") {
    for (const Probe& probe : probes()) {
        for (const double factor : {0.9, 1.1}) {
            INFO(probe.name << " at " << factor);
            const Sweep a = make_sweep();
            Sweep b = a;
            BiasPoint& point = b.points[1];
            const double allowance = probe.allowance(a.points[1]);
            probe.value(point) += factor * allowance;
            const SweepComparison c = compare_sweeps(a, b);
            CHECK(c.outcome == (factor < 1.0 ? CompareOutcome::match : CompareOutcome::differ));
            for (const auto& q : c.quantities) {
                INFO(q.name);
                if (q.name != probe.name) {
                    CHECK(q.max_ratio == 0.0);
                    continue;
                }
                CHECK(q.pass() == (factor < 1.0));
                // The allowance uses max(|a|, |b|): b's larger magnitude raises it by at most
                // relative * factor * allowance.
                CHECK(q.max_ratio > factor * 0.999);
                CHECK(q.max_ratio <= factor * 1.000001);
                CHECK(q.point == 1);
                CHECK(q.index == probe.index);
            }
        }
    }
}

TEST_CASE("a current's resolution is its floor (XG1)", "[analysis][compare]") {
    const Sweep a = make_sweep();
    Sweep b = a;
    b.points[0].terminal_current[2] = 0.0;
    b.points[0].terminal_current_resolution[2] = 2e-12;  // the floor is the larger of the two
    b.points[0].terminal_current[2] += 0.5 * 2e-12;
    const SweepComparison c = compare_sweeps(a, b);
    CHECK(c.outcome == CompareOutcome::match);
    const QuantityComparison* q = find(c, "terminal_current");
    REQUIRE(q != nullptr);
    CHECK(q->max_ratio > 0.49);
    CHECK(q->max_ratio < 0.51);

    b.points[0].terminal_current[2] = 3e-12;  // above the floor
    CHECK(compare_sweeps(a, b).outcome == CompareOutcome::differ);
}

TEST_CASE("a NaN on one side differs; quantities absent from both are not listed (XG1)",
          "[analysis][compare]") {
    Sweep a = make_sweep();
    Sweep b = a;
    b.points[2].fields.potential_V[1] = std::numeric_limits<double>::quiet_NaN();
    CHECK(compare_sweeps(a, b).outcome == CompareOutcome::differ);

    for (Sweep* s : {&a, &b})
        for (auto& p : s->points) {
            p.weighted_current.clear();
            p.weighted_current_resolution.clear();
            p.fields.temperature_K.clear();
        }
    const SweepComparison c = compare_sweeps(a, a);
    CHECK(c.outcome == CompareOutcome::match);
    CHECK(find(c, "weighted_current") == nullptr);
    CHECK(find(c, "temperature_K") == nullptr);
    CHECK(find(c, "potential_V") != nullptr);
}

TEST_CASE("every structural mismatch is not comparable, with its reason (XG1)", "[analysis][compare]") {
    const Sweep a = make_sweep();
    const std::vector<std::pair<std::string, std::function<void(Sweep&)>>> cases = {
        {"point count", [](Sweep& s) { s.points.pop_back(); }},
        {"stopped", [](Sweep& s) { s.stopped = NiTCAD::base::Error{NiTCAD::base::ErrorCode::invalid_input, "x", {}}; }},
        {"contact count", [](Sweep& s) { s.points[1].bias_V.pop_back(); }},
        {"the bias differs", [](Sweep& s) { s.points[1].bias_V[1] += 1e-15; }},
        {"node count", [](Sweep& s) { s.points[1].fields.potential_V.pop_back(); }},
        {"n_cm3", [](Sweep& s) { s.points[1].fields.n_cm3.pop_back(); }},
        {"temperature_K", [](Sweep& s) { s.points[2].fields.temperature_K.clear(); }},
        {"weighted_current", [](Sweep& s) { s.points[0].weighted_current.clear(); }},
        {"gate_charge", [](Sweep& s) { s.points[0].gate_charge.pop_back(); }},
    };
    for (const auto& [what, change] : cases) {
        INFO(what);
        Sweep b = a;
        change(b);
        const SweepComparison c = compare_sweeps(a, b);
        CHECK(c.outcome == CompareOutcome::not_comparable);
        CHECK(c.reason.find(what) != std::string::npos);
        CHECK(c.quantities.empty());
        const SweepComparison r = compare_sweeps(b, a);  // symmetric
        CHECK(r.outcome == CompareOutcome::not_comparable);
    }
}

TEST_CASE("identity and mesh digest are reported, not gated (XQ6)", "[analysis][compare]") {
    const Sweep a = make_sweep();
    Sweep b = a;
    b.run.input_identity ^= 1;
    b.run.mesh_digest ^= 1;
    const SweepComparison c = compare_sweeps(a, b);
    CHECK(c.outcome == CompareOutcome::match);
    CHECK_FALSE(c.identity_equal);
    CHECK_FALSE(c.mesh_digest_equal);
}

TEST_CASE("a sweep matches its own saved and loaded copy exactly (XG2)", "[analysis][compare]") {
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "nitcad_analysis_test" / "compare_round_trip";
    fs::create_directories(dir);
    const fs::path file = dir / "sweep.nitsweep";
    const Sweep a = make_sweep();
    REQUIRE(NiTCAD::io::save_sweep(file, a).has_value());
    const auto b = NiTCAD::io::load_sweep(file);
    REQUIRE(b.has_value());
    const SweepComparison c = compare_sweeps(a, *b);
    CHECK(c.outcome == CompareOutcome::match);
    CHECK(c.identity_equal);
    for (const auto& q : c.quantities) {
        INFO(q.name);
        CHECK(q.max_difference == 0.0);
    }
    fs::remove_all(dir);
}
