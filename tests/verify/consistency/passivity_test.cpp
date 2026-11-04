// C8, passivity at equilibrium (DECISIONS Unit 55 section 11.6 and the C8 rulings of 11.13): the Hermitian part of
// the admittance of C7's three-contact 2D MOS (small_signal_test's mos2d) at 0 V is positive semidefinite to the
// ruled bound at every frequency; the device biased, or with its thermal contacts at two temperatures, is NOT
// APPLICABLE; a conductance lowered at one frequency fails there; a synthetic Y with one negative eigenvalue passes
// at 0.99 of the bound and fails at 1.01; a NaN fails; each missing input, and more than 16 contacts, is NOT
// EVALUATED. The three TVR-1 decks (two contacts) pass with a zero-bias AC point added through the API and carry the
// two-contacts diagnostic.
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"
#include "NiTCAD/verify/reference/jacobi_eigen.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace c = NiTCAD::verify::consistency;
namespace v = NiTCAD::verify;
namespace solve = NiTCAD::solve;
namespace device = NiTCAD::device;
namespace mesh = NiTCAD::mesh;
namespace physics = NiTCAD::physics;
using v::Outcome;

namespace {

const std::vector<double> frequencies{0.0, 1e2, 1e5, 1e8, 1e11, 1e13};

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> x(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) x[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return x;
}

// tests/solve/small_signal_test.cpp's mos2d (Unit 21's 2D MOS structure): oxide on [-5 nm, 0) with the electrode
// on part of x_min; p-silicon (1e17) with an n+ pocket (1e19) under the source on y_max; the substrate on x_max; a
// fixed charge, two trap levels, a trap band and surface recombination at the interface.
device::Device mos2d() {
    const auto x = *mesh::straddle_interface(uniform(-5e-7, 1.2e-5, 30), 0.0, 1e-7);
    const auto t = uniform(0.0, 3e-6, 4);
    mesh::Mesh m = *mesh::make_tensor_grid(x, t);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        const bool oxide = p[0] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (oxide) continue;
        acceptors[i] = 1e17;
        if (p[1] > 2e-6 && p[0] < 3e-6) donors[i] = 1e19;
    }
    std::vector<mesh::NodeId> electrode, source;
    for (const mesh::NodeId node : m.find_boundary("x_min")->nodes) {
        if (m.points()[static_cast<std::size_t>(node)][1] < 1.5e-6) electrode.push_back(node);
    }
    for (const mesh::NodeId node : m.find_boundary("y_max")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(node)][0];
        if (px > 0.0 && px < 3e-6) source.push_back(node);
    }
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = 3e11;
    f.traps.levels = {{.type = physics::TrapType::donor, .density_cm2 = 2e11, .energy_eV = -0.2,
                       .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 3e-16},
                      {.type = physics::TrapType::acceptor, .density_cm2 = 1e11, .energy_eV = 0.25,
                       .sigma_n_cm2 = 2e-16, .sigma_p_cm2 = 1e-15}};
    f.traps.bands = {{.type = physics::TrapType::acceptor, .density_cm2_eV = 5e11, .energy_low_eV = -0.1,
                      .energy_high_eV = 0.4}};
    f.recombination_velocity_n_cm_s = 1e3;
    f.recombination_velocity_p_cm_s = 4e2;
    device::Contact gate{"gate", device::ContactKind::electrode, std::move(electrode)};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)},
                      {"source", device::ContactKind::ohmic, std::move(source)}},
         .interfaces = {std::move(f)}});
}

std::optional<NiTCAD::results::SmallSignal> small_signal(const device::Device& d, std::vector<double> bias,
                                                         solve::BiasOptions steady = {}) {
    const solve::SmallSignalOptions o{.steady = std::move(steady), .frequencies_Hz = frequencies};
    const std::vector<std::vector<double>> points{std::move(bias)};
    auto r = solve::solve_small_signal(d, points, o);
    if (!r) {
        std::printf("small signal failed: %s\n", r.error().message.c_str());
        return std::nullopt;
    }
    if (r->stopped) {
        std::printf("small signal stopped: %s\n", r->stopped->message.c_str());
        return std::nullopt;
    }
    return std::move(*r);
}

struct Mos {
    device::Device device;
    NiTCAD::results::SmallSignal zero, biased;
    c::DeviceData data;
};

const Mos& mos() {
    static const std::optional<Mos> s = []() -> std::optional<Mos> {
        device::Device d = mos2d();
        auto zero = small_signal(d, {0.0, 0.0, 0.0});
        auto biased = small_signal(d, {-0.4, 0.0, 0.1});
        if (!zero || !biased) return std::nullopt;
        Mos out{std::move(d), std::move(*zero), std::move(*biased), {}};
        out.data = NiTCAD::verify::benchmark::device_data(out.device);
        return out;
    }();
    REQUIRE(s.has_value());
    return *s;
}

std::filesystem::path source(const char* file) { return std::filesystem::path(NITCAD_SOURCE_DIR) / file; }

// The ruled tolerance's resolution term at frequency k: sqrt(n) ||resolution[k]||_2.
double weyl(const NiTCAD::results::SmallSignal& r, std::size_t k) {
    double s2 = 0.0;
    for (const double x : r.points[0].resolution[k]) s2 += x * x;
    return std::sqrt(static_cast<double>(r.contacts) * s2);
}

}  // namespace

TEST_CASE("C8 passes on a three-contact MOS at equilibrium", "[verify][consistency]") {
    const Mos& s = mos();
    const auto r = c::check_c8(s.zero, 0, s.device.mesh(), s.data);
    REQUIRE(r.id == "C8");
    REQUIRE(r.parts.size() == 1);
    REQUIRE(r.diagnostics.size() == 1);  // the smallest eigenvalue; three contacts: no two-contacts note
    const auto& lo = r.diagnostics.front();
    std::printf("C8 mos2d: residual %.3e tolerance %.3e at %s; smallest eigenvalue %.3e (tolerance %.3e) at %s\n",
                r.parts.front().residual, r.parts.front().tolerance, r.parts.front().where.c_str(), lo.value,
                lo.scale, lo.where.c_str());
    CAPTURE(r.comparison.reason, lo.value, lo.scale, lo.where);
    REQUIRE(r.comparison.outcome == Outcome::pass);
    REQUIRE(lo.name == "smallest eigenvalue");
    // Not vacuous from 1e5 Hz: the largest eigenvalue of the Hermitian part is far above the resolution term, so a
    // lowered conductance can go negative. At 0 and 100 Hz it is below it: jwC is symmetric and drops out of
    // (Y + Y^H)/2, and the dissipative conductance at equilibrium (gate current 0, the source junction at zero bias)
    // is below the resolution, so the check is vacuous there (recorded, DECISIONS 11.14). The ratios are printed.
    for (std::size_t k = 0; k < frequencies.size(); ++k) {
        const std::size_t n = s.zero.contacts;
        std::vector<std::complex<double>> h(n * n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                h[i * n + j] = 0.5 * (s.zero.admittance(0, k, i, j) + std::conj(s.zero.admittance(0, k, j, i)));
            }
            h[i * n + i] = {h[i * n + i].real(), 0.0};
        }
        const auto e = NiTCAD::verify::reference::jacobi_eigen(h, n);
        REQUIRE(e.status == NiTCAD::verify::reference::JacobiStatus::ok);
        std::printf("C8 mos2d: f = %g Hz, lambda_max %.3e, resolution term %.3e (ratio %.3g)\n", frequencies[k],
                    e.eigenvalues.back(), weyl(s.zero, k), e.eigenvalues.back() / weyl(s.zero, k));
        CAPTURE(k, e.eigenvalues.back(), weyl(s.zero, k));
        if (frequencies[k] >= 1e5) REQUIRE(e.eigenvalues.back() > 1e3 * weyl(s.zero, k));
        if (frequencies[k] < 1e5) REQUIRE(e.eigenvalues.back() < weyl(s.zero, k));  // vacuous here: recorded
    }
}

TEST_CASE("C8 is NOT APPLICABLE off equilibrium", "[verify][consistency]") {
    const Mos& s = mos();
    const auto na = c::check_c8(s.biased, 0, s.device.mesh(), s.data);
    REQUIRE(na.comparison.outcome == Outcome::not_applicable);
    REQUIRE(na.comparison.reason.find("all biases 0") != std::string::npos);
    auto hot = s.zero;
    hot.points[0].dc.thermal_bias_K = {300.0, 301.0};
    const auto t = c::check_c8(hot, 0, s.device.mesh(), s.data);
    REQUIRE(t.comparison.outcome == Outcome::not_applicable);
    REQUIRE(t.comparison.reason.find("thermal contacts at one temperature") != std::string::npos);
}

TEST_CASE("C8 fails on a lowered conductance at the named frequency", "[verify][consistency]") {
    const Mos& s = mos();
    const std::size_t k = 3, nc = s.zero.contacts;  // 1e8 Hz
    const auto base = c::check_c8(s.zero, 0, s.device.mesh(), s.data);
    // Re Y(gate, gate) lowered by delta: along Y's null vector u = (1, 1, 1)/sqrt(3), u^H H u falls by delta/3, so
    // lambda_min falls to at most (its value) - delta/3; delta = 6 x the larger tolerance puts it below -2 tol.
    const double delta = 6.0 * std::max(base.parts.front().tolerance, base.diagnostics.front().scale);
    auto bad = s.zero;
    bad.points[0].admittance[k][0 * nc + 0] -= delta;
    const auto r = c::check_c8(bad, 0, s.device.mesh(), s.data);
    CAPTURE(r.parts.front().residual, r.parts.front().tolerance, r.parts.front().where);
    REQUIRE(r.comparison.outcome == Outcome::fail);
    REQUIRE(r.parts.front().where == "f = 1e+08 Hz");
    REQUIRE(r.diagnostics.front().where == "f = 1e+08 Hz");
    REQUIRE(r.diagnostics.front().value == -r.parts.front().residual);  // the same eigenvalue and its tolerance
    REQUIRE(r.diagnostics.front().scale == r.parts.front().tolerance);

    // The same defect at the first (f = 0) and the last frequency, each with its own tolerance.
    for (const std::size_t q : {std::size_t{0}, frequencies.size() - 1}) {
        auto end = s.zero;
        const double rq = weyl(s.zero, q);
        double scale = 0.0;  // n x the largest |Y| entry bounds ||H||_2, for the eigensolver term
        for (const auto& yv : s.zero.points[0].admittance[q]) scale = std::max(scale, std::abs(yv));
        end.points[0].admittance[q][0] -= 6.0 * (rq + 8.0 * 3.0 * 3.0 * DBL_EPSILON * scale);  // 6 x its tolerance
        const auto rq_check = c::check_c8(end, 0, s.device.mesh(), s.data);
        CAPTURE(q, rq_check.parts.front().residual, rq_check.parts.front().tolerance, rq_check.parts.front().where);
        REQUIRE(rq_check.comparison.outcome == Outcome::fail);
        REQUIRE(rq_check.parts.front().where == (q == 0 ? "f = 0 Hz" : "f = 1e+13 Hz"));
    }

    // A skew change (Y_20 += d, Y_02 -= d, real) leaves H unchanged: still a pass.
    auto skew = s.zero;
    skew.points[0].admittance[k][2 * nc + 0] += 1e3 * delta;
    skew.points[0].admittance[k][0 * nc + 2] -= 1e3 * delta;
    REQUIRE(c::check_c8(skew, 0, s.device.mesh(), s.data).comparison.outcome == Outcome::pass);

    // A NaN admittance fails, it does not pass.
    auto nan = s.zero;
    nan.points[0].admittance[1][1 * nc + 2] = {std::numeric_limits<double>::quiet_NaN(), 0.0};
    REQUIRE(c::check_c8(nan, 0, s.device.mesh(), s.data).comparison.outcome == Outcome::fail);
}

TEST_CASE("C8's bound is sqrt(n) times the resolution's 2-norm plus TC9's accuracy", "[verify][consistency]") {
    // At one frequency Y = a P - c u u^T (P = I - u u^T, u = (1, 1, 1)/sqrt(3)): eigenvalues -c, a, a. With a = 10
    // x the resolution term the rounding terms are about 1e-14 of it; c = 0.99 x the bound passes, 1.01 fails.
    const Mos& s = mos();
    const std::size_t k = 2, nc = s.zero.contacts;
    const double rt = weyl(s.zero, k), a = 10.0 * rt;
    REQUIRE(rt > 0.0);
    for (const double factor : {0.99, 1.01}) {
        const double cc = factor * (rt + 8.0 * 3.0 * DBL_EPSILON * a);
        auto copy = s.zero;
        auto& y = copy.points[0].admittance[k];
        for (std::size_t i = 0; i < nc; ++i) {
            for (std::size_t j = 0; j < nc; ++j) y[i * nc + j] = (i == j ? a : 0.0) - (a + cc) / 3.0;
        }
        const auto r = c::check_c8(copy, 0, s.device.mesh(), s.data);
        CAPTURE(factor, r.parts.front().residual, r.parts.front().tolerance, r.parts.front().where);
        REQUIRE(r.comparison.outcome == (factor < 1.0 ? Outcome::pass : Outcome::fail));
        if (factor > 1.0) REQUIRE(r.parts.front().where == "f = 1e+05 Hz");
    }
    // The eigensolver's term: resolution zeroed, Y = a P (eigenvalues 0, a, a; the 0 computed to rounding) passes.
    auto exact = s.zero;
    exact.points[0].resolution[k].assign(nc, 0.0);
    auto& y = exact.points[0].admittance[k];
    for (std::size_t i = 0; i < nc; ++i) {
        for (std::size_t j = 0; j < nc; ++j) y[i * nc + j] = (i == j ? a : 0.0) - a / 3.0;
    }
    const auto z = c::check_c8(exact, 0, s.device.mesh(), s.data);
    CAPTURE(z.parts.front().residual, z.parts.front().tolerance);
    REQUIRE(z.comparison.outcome == Outcome::pass);
}

TEST_CASE("C8 is NOT EVALUATED without its inputs", "[verify][consistency]") {
    const Mos& s = mos();
    const auto& m = s.device.mesh();
    const auto reason = [](const c::CheckResult& r, const char* text) {
        CAPTURE(r.parts.front().reason, text);
        REQUIRE(r.comparison.outcome == Outcome::not_evaluated);
        REQUIRE(r.parts.front().reason.find(text) != std::string::npos);
    };
    reason(c::check_c8(s.zero, 1, m, s.data), "point index past");

    auto no_res = s.zero;
    no_res.points[0].resolution.clear();
    reason(c::check_c8(no_res, 0, m, s.data), "not one per frequency");

    auto short_y = s.zero;
    short_y.points[0].admittance[2].pop_back();
    reason(c::check_c8(short_y, 0, m, s.data), "frequency 2");

    auto short_r = s.zero;
    short_r.points[0].resolution[4].pop_back();
    reason(c::check_c8(short_r, 0, m, s.data), "frequency 4");

    auto fewer = s.data;
    fewer.contacts.pop_back();
    reason(c::check_c8(s.zero, 0, m, fewer), "contacts are not the device's");

    // 17 contacts (the run and the data padded alike): beyond TC9's scope.
    auto many = s.data;
    while (many.contacts.size() < 17) many.contacts.push_back(many.contacts.back());
    auto wide = s.zero;
    wide.contacts = 17;
    reason(c::check_c8(wide, 0, m, many), "more than 16 contacts");
}

TEST_CASE("C8 on the three TVR-1 decks with a zero-bias AC point", "[verify][consistency]") {
    for (const char* file : {"examples/pn_diode_1d.deck", "examples/moscap_cv.deck", "decks/schottky_diode_1d.deck"}) {
        CAPTURE(file);
        auto deck = NiTCAD::deck::read(source(file));
        REQUIRE(deck.has_value());
        auto plan = NiTCAD::deck::build(*deck);
        REQUIRE(plan.has_value());
        const auto run = small_signal(plan->device, std::vector<double>(plan->device.contacts().size(), 0.0),
                                      plan->options);
        REQUIRE(run.has_value());
        const auto data = NiTCAD::verify::benchmark::device_data(plan->device, plan->options.models);
        const auto r = c::check_c8(*run, 0, plan->device.mesh(), data);
        REQUIRE(r.diagnostics.size() == 2);
        const auto& lo = r.diagnostics.front();
        std::printf("C8 %s: residual %.3e tolerance %.3e at %s; smallest eigenvalue %.3e (tolerance %.3e) at %s\n",
                    file, r.parts.front().residual, r.parts.front().tolerance, r.parts.front().where.c_str(),
                    lo.value, lo.scale, lo.where.c_str());
        CAPTURE(r.comparison.reason);
        REQUIRE(r.comparison.outcome == Outcome::pass);
        REQUIRE(r.diagnostics.back().name == "two contacts");
    }
}

TEST_CASE("C8 is NOT EVALUATED with no contacts", "[verify][consistency]") {
    const Mos& s = mos();
    auto none = s.data;
    none.contacts.clear();
    auto empty = s.zero;
    empty.contacts = 0;
    for (auto& y : empty.points[0].admittance) y.clear();
    for (auto& r : empty.points[0].resolution) r.clear();
    const auto r = c::check_c8(empty, 0, s.device.mesh(), none);
    REQUIRE(r.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(r.parts.front().reason == "no contacts");
}
