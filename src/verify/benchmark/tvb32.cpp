// TVB-32 (tvb32.hpp). Moved from tests/verify/consistency/decks_test.cpp (RC-Q1); the test holds these entries to the
// expectations frozen in DECISIONS 11.13 before any run.
#include "NiTCAD/verify/benchmark/tvb32.hpp"

#include <algorithm>
#include <cstdio>
#include <string>
#include <utility>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "verify_decks.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace c = verify::consistency;

bool at_zero(const results::BiasPoint& p) {
    if (std::any_of(p.bias_V.begin(), p.bias_V.end(), [](double b) { return b != 0.0; })) return false;
    return p.thermal_bias_K.empty() || std::all_of(p.thermal_bias_K.begin(), p.thermal_bias_K.end(),
                                                   [&](double t) { return t == p.thermal_bias_K.front(); });
}

std::string bias_label(const results::BiasPoint& p) {
    std::string s = "bias (";
    for (std::size_t i = 0; i < p.bias_V.size(); ++i) s += (i == 0 ? "" : ", ") + shortest(p.bias_V[i]);
    return s + ") V";
}

struct Built {
    deck::Plan plan;
    c::DeviceData data;
};

std::optional<Built> build(std::string_view text, const char* name, std::string& error) {
    auto parsed = deck::parse(text, name);
    if (!parsed) {
        error = std::string(name) + ": " + parsed.error().message;
        return std::nullopt;
    }
    auto plan = deck::build(*parsed);
    if (!plan) {
        error = std::string(name) + ": " + plan.error().message;
        return std::nullopt;
    }
    plan->outputs = {};  // nothing written
    auto data = device_data(plan->device, plan->options.models);
    return Built{std::move(*plan), std::move(data)};
}

struct Collector {
    std::vector<Tvb32Entry>& out;
    std::string deck;
    std::size_t contacts;
    void add(const std::string& at, bool zero, c::CheckResult r) {
        out.push_back({deck, at, zero, contacts, std::move(r)});
    }
    void dc(const results::BiasPoint& p, const mesh::Mesh& m, const c::DeviceData& d, const std::string& at) {
        const bool z = at_zero(p);
        add(at, z, c::check_c1(p, m, d));
        add(at, z, c::check_c2(p, m, d));
        add(at, z, c::check_c3(p, m, d));
        add(at, z, c::check_c4(p, m, d));
        add(at, z, c::check_c5(p, m, d));
        add(at, z, c::check_c12(p.fields, &p, m, d));
    }
    void ac(const results::SmallSignal& run, const mesh::Mesh& m, const c::DeviceData& d, const std::string& prefix) {
        for (std::size_t k = 0; k < run.points.size(); ++k) {
            const auto& p = run.points[k];
            const std::string at = prefix + bias_label(p.dc);
            dc(p.dc, m, d, at);
            add(at, at_zero(p.dc), c::check_c7(run, k, m, d));
            add(at, at_zero(p.dc), c::check_c8(run, k, m, d));
        }
    }
};

bool equilibrium(const Built& b, Collector& col, std::string& error) {
    solve::EquilibriumOptions eo;
    eo.models = b.plan.options.models;
    const auto eq = solve::solve_equilibrium(b.plan.device, eo);
    if (!eq) {
        error = col.deck + " equilibrium: " + eq.error().message;
        return false;
    }
    col.add("equilibrium", true, c::check_c1(*eq, b.plan.device.mesh(), b.data));
    col.add("equilibrium", true, c::check_c3(*eq, b.plan.device.mesh(), b.data));
    return true;
}

}  // namespace

Tvb32Run tvb32_entries() {
    Tvb32Run run;
    std::string error;
    for (const auto& [text, name] : {std::pair{decks::pn_diode_1d, "pn_diode_1d"},
                                     std::pair{decks::schottky_diode_1d, "schottky_diode_1d"}}) {
        const auto b = build(text, name, error);
        if (!b) return {std::move(run.entries), error};
        Collector col{run.entries, name, b->plan.device.contacts().size()};
        const auto r = deck::run(b->plan);
        if (!r || !r->sweep || r->sweep->stopped || r->sweep->points.size() != b->plan.points.size()) {
            return {std::move(run.entries), std::string(name) + ": the sweep did not complete"};
        }
        if (!equilibrium(*b, col, error)) return {std::move(run.entries), error};
        for (const auto& p : r->sweep->points) col.dc(p, b->plan.device.mesh(), b->data, bias_label(p));
    }
    {
        const char* name = "moscap_cv";
        const auto b = build(decks::moscap_cv, name, error);
        if (!b) return {std::move(run.entries), error};
        Collector col{run.entries, name, b->plan.device.contacts().size()};
        const auto r = deck::run(b->plan);
        if (!r || !r->small_signal || r->small_signal->stopped ||
            r->small_signal->points.size() != b->plan.points.size()) {
            return {std::move(run.entries), std::string(name) + ": the C-V run did not complete"};
        }
        if (!equilibrium(*b, col, error)) return {std::move(run.entries), error};
        col.ac(*r->small_signal, b->plan.device.mesh(), b->data, "");
        // The zero-bias AC point added through the API (TVR1-Q5), at the deck's own frequencies; the deck not edited.
        auto o = b->plan.small_signal;
        o.steady = b->plan.options;
        const std::vector<std::vector<double>> zero{std::vector<double>(b->plan.device.contacts().size(), 0.0)};
        const auto added = solve::solve_small_signal(b->plan.device, zero, o);
        if (!added || added->stopped) return {std::move(run.entries), "moscap_cv: the added AC point did not solve"};
        col.ac(*added, b->plan.device.mesh(), b->data, "added AC point, ");
    }
    return run;
}

BenchmarkResult run_tvb32() {
    BenchmarkResult r;
    r.id = "TVB-32";
    r.title = "The consistency suite C1-C12 on the three TVR-1 decks";
    Tvb32Run t = tvb32_entries();
    r.info = consistency::check_info("C1");  // class D, I2; each check judged with its own contract
    r.info.id = "TVB-32";
    r.info.citation = "DECISIONS Unit 55 11.6: checks C1-C12, each judged with its own contract (check_info)";
    for (Tvb32Entry& e : t.entries) {
        double worst = 0.0;
        for (const auto& p : e.result.parts) {
            if (p.tolerance > 0.0) worst = std::max(worst, p.residual / p.tolerance);
        }
        ComparisonRecord rec{e.result.id + " [" + e.deck + "] " + e.at, "residual / tolerance", 0.0, worst, 0.0,
                             0.0, worst, "every part's residual <= its tolerance", std::move(e.result.comparison)};
        r.comparisons.push_back(std::move(rec));
    }
    r.notes.push_back("C6, C9, C10, C11: input not produced by the decks (their real-solve evidence is in their tests)");
    r.run_error = std::move(t.run_error);
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
