// C9, time-frequency consistency at low frequency (DECISIONS Unit 55 section 11.6 and the two C9 rulings of 11.13):
// on examples/moscap_cv.deck (a lumped gate; its row's capacitance is judged, the body's is not stored) and
// examples/pn_diode_1d.deck (conductance only: no gate, so the check is NOT EVALUATED overall), each with an AC
// point and four DC points about it solved through the API (the decks are not edited): the judged parts pass; each
// injected defect fails; each missing input makes it NOT EVALUATED. The Kramers-Kronig diagnostic is shown on a
// synthetic Debye admittance (reconstructed to a small fraction of its dispersion) and reported for the MOS-C.
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace c = NiTCAD::verify::consistency;
namespace v = NiTCAD::verify;
namespace solve = NiTCAD::solve;
namespace results = NiTCAD::results;
using v::Outcome;

namespace {

std::filesystem::path source(const char* file) { return std::filesystem::path(NITCAD_SOURCE_DIR) / file; }

// 0 Hz and 100 Hz to 1e12 Hz at 8 per decade. (From 1 Hz, the gate column's resolution / omega1, whose rounding
// part does not fall with omega, was 3e-3 of the capacitance: measured.)
std::vector<double> band() {
    std::vector<double> f{0.0};
    for (int k = 0; k <= 80; ++k) f.push_back(std::pow(10.0, 2.0 + k / 8.0));
    return f;
}

struct Case {
    NiTCAD::deck::Plan plan;
    results::SmallSignal ac;
    results::Sweep dc;  // V - 2h, V - h, V + h, V + 2h on the stencil's contact
    c::DeviceData data;
    c::DcStencil stencil;
    const NiTCAD::mesh::Mesh& mesh() const { return plan.device.mesh(); }
};

std::optional<Case> solve_case(const char* file, std::vector<double> bias, std::size_t contact, double h,
                               std::vector<double> frequencies) {
    auto deck = NiTCAD::deck::read(source(file));
    if (!deck) return std::nullopt;
    auto plan = NiTCAD::deck::build(*deck);
    if (!plan) return std::nullopt;
    const solve::SmallSignalOptions o{.steady = plan->options, .frequencies_Hz = std::move(frequencies)};
    const std::vector<std::vector<double>> at{bias};
    auto ac = solve::solve_small_signal(plan->device, at, o);
    if (!ac || ac->stopped) return std::nullopt;
    std::vector<std::vector<double>> stencil;
    for (const double m : {-2.0, -1.0, 1.0, 2.0}) {
        auto b = bias;
        b[contact] += m * h;
        stencil.push_back(std::move(b));
    }
    auto dc = solve::sweep_bias(plan->device, stencil, plan->options);
    if (!dc || dc->stopped || dc->points.size() != 4) return std::nullopt;
    Case out{std::move(*plan), std::move(*ac), std::move(*dc), {}, {}};
    out.data = NiTCAD::verify::benchmark::device_data(out.plan.device, out.plan.options.models);
    out.stencil.contact = contact;
    out.stencil.h_V = h;
    return out;
}

// The stencil's pointers into this case's own sweep (set after any copy).
c::DcStencil stencil_of(const Case& s) {
    c::DcStencil st = s.stencil;
    for (std::size_t q = 0; q < 4; ++q) st.points[q] = &s.dc.points[q];
    return st;
}

const Case& moscap() {
    // gate, body: depletion at -0.5 V on the gate (the deck's C-V sweep runs -2 to 0.25 V).
    static const std::optional<Case> s = solve_case("examples/moscap_cv.deck", {-0.5, 0.0}, 0, 2e-3, band());
    REQUIRE(s.has_value());
    return *s;
}

const Case& diode() {
    // cathode, anode: forward 0.5 V.
    static const std::optional<Case> s =
        solve_case("examples/pn_diode_1d.deck", {0.0, 0.5}, 1, 1e-3, {0.0, 1.0, 10.0});
    REQUIRE(s.has_value());
    return *s;
}

const c::CheckPart& part(const c::CheckResult& r, const std::string& prefix) {
    const c::CheckPart* found = nullptr;
    for (const auto& p : r.parts) {
        if (found == nullptr && p.name.rfind(prefix, 0) == 0) found = &p;
    }
    CAPTURE(prefix);
    REQUIRE(found != nullptr);
    return *found;
}

const c::Diagnostic* diagnostic(const c::CheckResult& r, const std::string& name, const std::string& where = "") {
    for (const auto& d : r.diagnostics) {
        if (d.name == name && d.where.rfind(where, 0) == 0) return &d;
    }
    return nullptr;
}

}  // namespace

TEST_CASE("C9 passes on the MOS-C's gate row", "[verify][consistency]") {
    const Case& s = moscap();
    const auto r = c::check_c9(s.ac, 0, stencil_of(s), s.mesh(), s.data);
    const auto& g = part(r, "conductance");
    const auto& q = part(r, "capacitance");
    CAPTURE(r.comparison.reason, g.residual, g.tolerance, g.where, q.residual, q.tolerance, q.where);
    REQUIRE(r.id == "C9");
    REQUIRE(r.comparison.outcome == Outcome::pass);
    REQUIRE(g.outcome == Outcome::pass);
    REQUIRE(q.outcome == Outcome::pass);
    // Not vacuous: the bound is a small fraction of the capacitance it judges.
    const double C = s.ac.points[0].admittance[1][0].imag() / (2.0 * std::numbers::pi * s.ac.frequency_Hz[1]);
    const auto& Qs = s.dc.points;
    const double D4 = (8.0 * (Qs[2].gate_charge[0] - Qs[1].gate_charge[0]) - (Qs[3].gate_charge[0] - Qs[0].gate_charge[0])) / (12.0 * s.stencil.h_V);
    std::printf("C9 MOS-C: C(f1) %.17g, D4 %.17g, Q %.17g %.17g %.17g %.17g\n", C, D4, Qs[0].gate_charge[0], Qs[1].gate_charge[0], Qs[2].gate_charge[0], Qs[3].gate_charge[0]);
    CAPTURE(C);
    REQUIRE(q.tolerance < 1e-3 * std::abs(C));
    REQUIRE(diagnostic(r, "ohmic rows' capacitance") != nullptr);
    for (const auto& d : r.diagnostics) {
        std::printf("C9 MOS-C diagnostic %s %s: %.3e (scale %.3e) %s\n", d.name.c_str(), d.where.c_str(), d.value,
                    d.scale, d.note.c_str());
    }
    REQUIRE(diagnostic(r, "Kramers-Kronig", "contact gate") != nullptr);
}

TEST_CASE("C9 on the pn diode: conductance passes, capacitance NOT EVALUATED", "[verify][consistency]") {
    const Case& s = diode();
    const auto r = c::check_c9(s.ac, 0, stencil_of(s), s.mesh(), s.data);
    const auto& g = part(r, "conductance");
    CAPTURE(g.residual, g.tolerance, g.where);
    REQUIRE(g.outcome == Outcome::pass);
    // The conductance is resolved and the bound is a small fraction of it. (The DC resolution term, ~ 1/h, and
    // |D4 - D2|, ~ h^2, balance near h = 1 mV at about 1e-3 of G here: measured 6e-3 at 1 mV, 1.6e-2 at 0.3 mV.)
    REQUIRE(g.tolerance < 1e-2 * std::abs(s.ac.points[0].admittance[0][1 * 2 + 1].real()));
    REQUIRE(part(r, "capacitance").outcome == Outcome::not_evaluated);
    REQUIRE(r.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(r.comparison.reason.find("no gate or electrode") != std::string::npos);
    const auto* kk = diagnostic(r, "Kramers-Kronig");
    REQUIRE(kk != nullptr);
    REQUIRE(kk->note.find("not estimated") != std::string::npos);
}

// The ruled conductance bound of row i (DECISIONS 11.13, C9), restated here so a defect can be sized against the
// row it is injected in (the part reports only its worst row).
double conductance_bound(const Case& s, std::size_t i) {
    const double h = s.stencil.h_V;
    std::array<double, 4> I{}, r{};
    for (std::size_t q = 0; q < 4; ++q) {
        I[q] = s.dc.points[q].terminal_current[i];
        r[q] = s.dc.points[q].terminal_current_resolution[i];
    }
    const double d4 = (8.0 * (I[2] - I[1]) - (I[3] - I[0])) / (12.0 * h), d2 = (I[2] - I[1]) / (2.0 * h);
    const double terms = 8.0 * (std::abs(I[1]) + std::abs(I[2])) + std::abs(I[0]) + std::abs(I[3]);
    return (8.0 * (r[1] + r[2]) + r[0] + r[3]) / (12.0 * h) + std::abs(d4 - d2) +
           64.0 * DBL_EPSILON * terms / (12.0 * h) + s.ac.points[0].resolution[0][s.stencil.contact];
}

TEST_CASE("C9 fails on injected defects", "[verify][consistency]") {
    const Case& s = moscap();
    const auto base = c::check_c9(s.ac, 0, stencil_of(s), s.mesh(), s.data);
    const double g_tol = conductance_bound(s, 1), q_tol = part(base, "capacitance").tolerance;  // body row; gate row
    CAPTURE(g_tol, q_tol, part(base, "conductance").tolerance);
    REQUIRE(g_tol > 0.0);
    const std::size_t k0 = 0, k1 = 1;  // 0 Hz, 100 Hz
    const double w1 = 2.0 * std::numbers::pi * s.ac.frequency_Hz[k1];

    auto g = s.ac;  // Re Y(body, gate) at 0 Hz moved by twice the body row's bound
    g.points[0].admittance[k0][1 * 2 + 0] += 2.0 * g_tol;
    const auto rg = c::check_c9(g, 0, stencil_of(s), s.mesh(), s.data);
    REQUIRE(part(rg, "conductance").outcome == Outcome::fail);
    REQUIRE(part(rg, "conductance").where == "Y(body, gate)");
    REQUIRE(rg.comparison.outcome == Outcome::fail);

    // Im Y(gate, gate) scaled at every frequency so C(f1) moves by three times the capacitance bound (a systematic
    // error: the omega^2 estimate, formed from C(f1) and C(f2), is unchanged by a common scale to first order). A
    // change at f1 alone is absorbed by that estimate (it grows by f1^2 / (f2^2 - f1^2) = 1.29 of the change at 8
    // per decade): the second half, recorded as a limit of the ruled bound.
    const double C1 = s.ac.points[0].admittance[k1][0].imag() / w1;
    auto q = s.ac;
    for (auto& y : q.points[0].admittance) y[0] = {y[0].real(), y[0].imag() * (1.0 + 3.0 * q_tol / std::abs(C1))};
    const auto rq = c::check_c9(q, 0, stencil_of(s), s.mesh(), s.data);
    CAPTURE(part(rq, "capacitance").residual, part(rq, "capacitance").tolerance);
    REQUIRE(part(rq, "capacitance").outcome == Outcome::fail);
    auto single = s.ac;
    single.points[0].admittance[k1][0] += std::complex<double>(0.0, 3.0 * q_tol * w1);
    const auto rs = c::check_c9(single, 0, stencil_of(s), s.mesh(), s.data);
    CAPTURE(part(rs, "capacitance").residual, part(rs, "capacitance").tolerance);
    REQUIRE(part(rs, "capacitance").outcome == Outcome::pass);  // the limit: absorbed

    // A DC value at V + h (or V - h) moved by dq: D4 moves by 8 dq / 12h and |D4 - D2| by at most 2 dq / 12h, so
    // dq = 6 h tol leaves a residual >= 3 tol against a bound <= 2 tol.
    Case dc = s;  // the gate charge at V + h
    dc.dc.points[2].gate_charge[0] += 6.0 * s.stencil.h_V * q_tol;
    REQUIRE(part(c::check_c9(dc.ac, 0, stencil_of(dc), dc.mesh(), dc.data), "capacitance").outcome == Outcome::fail);

    Case di = s;  // the body current at V - h
    di.dc.points[1].terminal_current[1] -= 6.0 * s.stencil.h_V * g_tol;
    REQUIRE(part(c::check_c9(di.ac, 0, stencil_of(di), di.mesh(), di.data), "conductance").outcome == Outcome::fail);
}

TEST_CASE("C9's bounds at their edge: 0.99 passes, 1.01 fails", "[verify][consistency]") {
    // Conductance: the diode's anode row (its bound restated; every term of it is a visible share there). The
    // residual is set to f times the bound by moving Re Y(anode, anode) at 0 Hz away from D4.
    const Case& s = diode();
    const std::size_t i = 1, nc = 2;
    const double bound = conductance_bound(s, i);
    const double h = s.stencil.h_V;
    std::array<double, 4> I{};
    for (std::size_t q = 0; q < 4; ++q) I[q] = s.dc.points[q].terminal_current[i];
    const double d4 = (8.0 * (I[2] - I[1]) - (I[3] - I[0])) / (12.0 * h);
    const double ac_res = s.ac.points[0].resolution[0][s.stencil.contact];
    CAPTURE(bound, ac_res, d4);  // the AC resolution is 1.2e-5 of this bound: its share is shown below instead
    for (const double f : {0.99, 1.01}) {
        auto y = s.ac;
        const double now = y.points[0].admittance[0][i * nc + i].real();
        y.points[0].admittance[0][i * nc + i] = {d4 + (now >= d4 ? 1.0 : -1.0) * f * bound, 0.0};
        const auto r = c::check_c9(y, 0, stencil_of(s), s.mesh(), s.data);
        CAPTURE(f, part(r, "conductance").residual, part(r, "conductance").tolerance);
        REQUIRE(part(r, "conductance").outcome == (f < 1.0 ? Outcome::pass : Outcome::fail));
    }

    // Conductance on the MOS-C's gate row: a gate carries no DC current and its resolution is 0, so the bound is
    // the column's AC resolution alone (D4 = 0).
    const Case& m = moscap();
    const double gate_bound = conductance_bound(m, 0);
    CAPTURE(gate_bound, m.ac.points[0].resolution[0][0]);
    REQUIRE(gate_bound == m.ac.points[0].resolution[0][0]);
    for (const double f : {0.99, 1.01}) {
        auto y = m.ac;
        y.points[0].admittance[0][0] = {f * gate_bound, 0.0};
        const auto r = c::check_c9(y, 0, stencil_of(m), m.mesh(), m.data);
        CAPTURE(f, part(r, "conductance").residual, part(r, "conductance").tolerance);
        REQUIRE(part(r, "conductance").outcome == (f < 1.0 ? Outcome::pass : Outcome::fail));
    }

    // Capacitance: the MOS-C's gate row; Im Y(gate, gate) scaled at every frequency (the omega^2 term nearly
    // unchanged) so the residual is f times the bound.
    const auto base = c::check_c9(m.ac, 0, stencil_of(m), m.mesh(), m.data);
    const auto& q = part(base, "capacitance");
    const double w1 = 2.0 * std::numbers::pi * m.ac.frequency_Hz[1];
    const double c1 = m.ac.points[0].admittance[1][0].imag() / w1;
    const auto& Q = m.dc.points;
    const double D4 = (8.0 * (Q[2].gate_charge[0] - Q[1].gate_charge[0]) - (Q[3].gate_charge[0] - Q[0].gate_charge[0])) /
                      (12.0 * m.stencil.h_V);
    for (const double f : {0.99, 1.01}) {
        const double target = D4 + (c1 >= D4 ? 1.0 : -1.0) * f * q.tolerance;
        const double scale = target / c1;
        auto y = m.ac;
        for (auto& row : y.points[0].admittance) row[0] = {row[0].real(), row[0].imag() * scale};
        const auto r = c::check_c9(y, 0, stencil_of(m), m.mesh(), m.data);
        CAPTURE(f, part(r, "capacitance").residual, part(r, "capacitance").tolerance);
        REQUIRE(part(r, "capacitance").outcome == (f < 1.0 ? Outcome::pass : Outcome::fail));
    }
}

TEST_CASE("C9 is NOT EVALUATED without its inputs", "[verify][consistency]") {
    const Case& s = moscap();
    const auto& m = s.mesh();
    REQUIRE(c::check_c9(s.ac, 1, stencil_of(s), m, s.data).comparison.outcome == Outcome::not_evaluated);

    auto wrong_h = stencil_of(s);  // the stencil read with another step
    wrong_h.h_V = 2e-2;
    const auto a = c::check_c9(s.ac, 0, wrong_h, m, s.data);
    REQUIRE(a.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(a.comparison.reason.find("stencil not about the AC point") != std::string::npos);

    auto missing = stencil_of(s);
    missing.points[2] = nullptr;
    REQUIRE(c::check_c9(s.ac, 0, missing, m, s.data).comparison.outcome == Outcome::not_evaluated);

    auto other = stencil_of(s);  // the stencil on the body while the points step the gate
    other.contact = 1;
    REQUIRE(c::check_c9(s.ac, 0, other, m, s.data).comparison.outcome == Outcome::not_evaluated);

    auto off = stencil_of(s);  // the step read 1e-9 relative off: the points are not at V + m h to 4 eps
    off.h_V = s.stencil.h_V * (1.0 + 1e-9);
    REQUIRE(c::check_c9(s.ac, 0, off, m, s.data).comparison.outcome == Outcome::not_evaluated);

    auto beyond = stencil_of(s);  // a contact past the device's, every stencil point at the AC point's biases
    beyond.contact = 5;
    for (auto& p : beyond.points) p = &s.ac.points[0].dc;
    const auto b5 = c::check_c9(s.ac, 0, beyond, m, s.data);
    REQUIRE(b5.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(b5.comparison.reason.find("stencil contact") != std::string::npos);

    auto no_zero = s.ac;  // f = 0 relabelled: the conductance part and Kramers-Kronig have no reference
    no_zero.frequency_Hz[0] = 1e-3;
    const auto z = c::check_c9(no_zero, 0, stencil_of(s), m, s.data);
    REQUIRE(part(z, "conductance").outcome == Outcome::not_evaluated);
    REQUIRE(part(z, "conductance").reason.find("f = 0") != std::string::npos);

    Case nq = s;  // no DC gate charges
    for (auto& p : nq.dc.points) p.gate_charge.clear();
    REQUIRE(part(c::check_c9(nq.ac, 0, stencil_of(nq), nq.mesh(), nq.data), "capacitance").outcome ==
            Outcome::not_evaluated);

    auto one = s.ac;  // a single nonzero frequency: no omega^2 term
    one.frequency_Hz.resize(2);
    one.points[0].admittance.resize(2);
    one.points[0].resolution.resize(2);
    const auto o = c::check_c9(one, 0, stencil_of(s), m, s.data);
    REQUIRE(part(o, "capacitance").outcome == Outcome::not_evaluated);
    const auto* kk = diagnostic(o, "Kramers-Kronig");
    REQUIRE(kk != nullptr);
    REQUIRE(kk->note.find("not estimated") != std::string::npos);
}

TEST_CASE("C9's Kramers-Kronig diagnostic on a synthetic Debye admittance", "[verify][consistency]") {
    // Y = G0 + j omega (C_inf + dC / (1 + j omega tau)) on one contact (a causal response): C~ = (Y - Y(0)) / (j
    // omega) is the Debye function, reconstructed over 1e-4 / tau .. 1e4 / tau at 8 per decade to a small fraction
    // of dC. (G0 is subtracted with Y(0); left in, it would put -j G0 / omega into C~.) A non-causal change (a real
    // bump in C~ with no imaginary partner) is seen.
    const double C_inf = 1e-7, dC = 4e-7, tau = 1e-6, G0 = 1e-4;
    results::SmallSignal run;
    run.contacts = 1;
    run.frequency_Hz.push_back(0.0);
    for (int k = 0; k <= 64; ++k) run.frequency_Hz.push_back(std::pow(10.0, -4.0 + k / 8.0) / (2.0 * std::numbers::pi * tau));
    results::SmallSignalPoint p;
    p.dc.bias_V = {0.0};
    for (const double f : run.frequency_Hz) {
        const std::complex<double> jw(0.0, 2.0 * std::numbers::pi * f);
        p.admittance.push_back({G0 + jw * (C_inf + dC / (1.0 + jw * tau))});
        p.resolution.push_back({0.0});
    }
    run.points.push_back(p);
    c::DeviceData d;
    d.contacts = {{.name = "a"}};
    const auto m = NiTCAD::mesh::Mesh::from_parts(1, {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, {0.5, 0.5},
                                                  {{0, 1, 1.0, 1.0}}, {});
    REQUIRE(m.has_value());
    std::array<results::BiasPoint, 4> dc{};
    c::DcStencil st{.contact = 0, .h_V = 1e-3, .points = {&dc[0], &dc[1], &dc[2], &dc[3]}};
    for (std::size_t q = 0; q < 4; ++q) dc[q].bias_V = {std::array{-2.0, -1.0, 1.0, 2.0}[q] * 1e-3};

    const auto r = c::check_c9(run, 0, st, *m, d);
    const auto* kk = diagnostic(r, "Kramers-Kronig", "contact a");
    REQUIRE(kk != nullptr);
    CAPTURE(kk->value, kk->scale, kk->where);
    REQUIRE(std::isfinite(kk->value));
    REQUIRE(kk->value < 1e-2 * dC);

    auto bumped = run;  // Re C~ raised by 0.1 dC at one mid-band sample, Im unchanged
    const std::size_t k = 33;
    const std::complex<double> jw(0.0, 2.0 * std::numbers::pi * run.frequency_Hz[k]);
    bumped.points[0].admittance[k][0] += jw * (0.1 * dC);
    const auto rb = c::check_c9(bumped, 0, st, *m, d);
    const auto* kb = diagnostic(rb, "Kramers-Kronig", "contact a");
    REQUIRE(kb != nullptr);
    CAPTURE(kb->value);
    REQUIRE(kb->value > 0.05 * dC);

    // The sampling rules: 3 decades at 8 per decade, and 8 decades at 4 per decade, are not estimated.
    const auto sampled = [&](double decades, int per_decade) {
        results::SmallSignal s;
        s.contacts = 1;
        s.frequency_Hz.push_back(0.0);
        results::SmallSignalPoint sp;
        sp.dc.bias_V = {0.0};
        sp.admittance.push_back({0.0});
        sp.resolution.push_back({0.0});
        const int n = static_cast<int>(decades * per_decade);
        for (int k = 0; k <= n; ++k) {
            const double f = std::pow(10.0, -4.0 + static_cast<double>(k) / per_decade) / (2.0 * std::numbers::pi * tau);
            const std::complex<double> jw(0.0, 2.0 * std::numbers::pi * f);
            s.frequency_Hz.push_back(f);
            sp.admittance.push_back({jw * (C_inf + dC / (1.0 + jw * tau))});
            sp.resolution.push_back({0.0});
        }
        s.points.push_back(sp);
        return c::check_c9(s, 0, st, *m, d);
    };
    const auto r3 = sampled(3.0, 8);
    const auto* few_decades = diagnostic(r3, "Kramers-Kronig");
    REQUIRE(few_decades != nullptr);
    REQUIRE(few_decades->note.find("fewer than 4 decades") != std::string::npos);
    const auto r4 = sampled(8.0, 4);
    const auto* sparse = diagnostic(r4, "Kramers-Kronig");
    REQUIRE(sparse != nullptr);
    REQUIRE(sparse->note.find("fewer than 8 samples per decade") != std::string::npos);
}
