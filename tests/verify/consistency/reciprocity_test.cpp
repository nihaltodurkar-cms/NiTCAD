// C7, reciprocity of the admittance at equilibrium (DECISIONS Unit 55 section 11.6 and the four C7 rulings of
// 11.13): passes on a three-contact 2D MOS written here (an electrode gate, a source pocket and the substrate, with
// a charged, trapped interface; small_signal_test's mos2d) at 0 V, where it is a real test of reciprocity; the same
// device biased is NOT APPLICABLE, and judged with the predicate bypassed on a copy it fails (a biased device is
// not reciprocal); each injected defect fails it at the named place; each missing input makes it NOT EVALUATED. The
// three TVR-1 decks (two contacts each) are run with a zero-bias AC point added through the API (the decks are not
// edited): they pass, and carry the diagnostic that two contacts show conservation, not reciprocity.
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

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
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

// The largest ratio |Y_ij - Y_ji| / (res_i + res_j) over the run's one point (reported).
double worst_ratio(const NiTCAD::results::SmallSignal& r) {
    double w = 0.0;
    const std::size_t nc = r.contacts;
    for (std::size_t k = 0; k < r.frequency_Hz.size(); ++k) {
        const auto& res = r.points[0].resolution[k];
        for (std::size_t i = 0; i < nc; ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                w = std::max(w, std::abs(r.admittance(0, k, i, j) - r.admittance(0, k, j, i)) / (res[i] + res[j]));
            }
        }
    }
    return w;
}

std::filesystem::path source(const char* file) { return std::filesystem::path(NITCAD_SOURCE_DIR) / file; }

}  // namespace

TEST_CASE("C7 passes on a three-contact MOS at equilibrium", "[verify][consistency]") {
    const Mos& s = mos();
    const auto r = c::check_c7(s.zero, 0, s.device.mesh(), s.data);
    CAPTURE(r.comparison.reason, r.parts.front().residual, r.parts.front().tolerance, r.parts.front().where,
            worst_ratio(s.zero));
    REQUIRE(r.id == "C7");
    REQUIRE(r.comparison.outcome == Outcome::pass);
    REQUIRE(r.parts.size() == 1);
    REQUIRE(r.parts.front().outcome == Outcome::pass);
    REQUIRE(r.parts.front().tolerance > 0.0);
    REQUIRE(r.diagnostics.empty());  // three contacts: no conservation-only note
    // The check is not vacuous here: the off-diagonal admittances are well above their resolution at high f.
    const std::size_t top = frequencies.size() - 1;
    const auto& res = s.zero.points[0].resolution[top];
    double off = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        scale = std::max(scale, res[i]);
        for (std::size_t j = 0; j < 3; ++j) {
            if (i != j) off = std::max(off, std::abs(s.zero.admittance(0, top, i, j)));
        }
    }
    CAPTURE(off, scale);
    REQUIRE(off > 1e3 * scale);
}

TEST_CASE("C7 is NOT APPLICABLE off equilibrium, and a biased device is not reciprocal", "[verify][consistency]") {
    const Mos& s = mos();
    const auto na = c::check_c7(s.biased, 0, s.device.mesh(), s.data);
    REQUIRE(na.comparison.outcome == Outcome::not_applicable);
    REQUIRE(na.comparison.reason.find("all biases 0") != std::string::npos);

    // Defect: the predicate bypassed on a copy (the biases relabelled 0); the biased admittance is judged and fails.
    auto copy = s.biased;
    copy.points[0].dc.bias_V.assign(copy.points[0].dc.bias_V.size(), 0.0);
    const auto r = c::check_c7(copy, 0, s.device.mesh(), s.data);
    CAPTURE(r.parts.front().residual, r.parts.front().tolerance, r.parts.front().where, worst_ratio(s.biased));
    REQUIRE(r.comparison.outcome == Outcome::fail);
    REQUIRE(r.parts.front().residual > 10.0 * r.parts.front().tolerance);

    // Thermal contacts at two temperatures (a copy of the 0 V point).
    auto hot = s.zero;
    hot.points[0].dc.thermal_bias_K = {300.0, 301.0};
    const auto t = c::check_c7(hot, 0, s.device.mesh(), s.data);
    REQUIRE(t.comparison.outcome == Outcome::not_applicable);
    REQUIRE(t.comparison.reason.find("thermal contacts at one temperature") != std::string::npos);
}

TEST_CASE("C7 fails on an injected asymmetry at the named place", "[verify][consistency]") {
    const Mos& s = mos();
    const std::size_t k = 3;  // 1e8 Hz
    const auto& res = s.zero.points[0].resolution[k];
    const std::size_t nc = s.zero.contacts;

    // Y(source, gate) moved by twice the pair's resolution.
    auto bad = s.zero;
    bad.points[0].admittance[k][2 * nc + 0] += std::complex<double>(0.0, 2.0 * (res[0] + res[2]));
    const auto r = c::check_c7(bad, 0, s.device.mesh(), s.data);
    CAPTURE(r.parts.front().where);
    REQUIRE(r.comparison.outcome == Outcome::fail);
    REQUIRE(r.parts.front().where == "f = 1e+08 Hz, Y(source, gate)");

    // The other triangle (Y(gate, substrate)) at f = 0, real part.
    auto low = s.zero;
    const auto& r0 = s.zero.points[0].resolution[0];
    low.points[0].admittance[0][0 * nc + 1] += 2.0 * (r0[0] + r0[1]);
    const auto q = c::check_c7(low, 0, s.device.mesh(), s.data);
    REQUIRE(q.comparison.outcome == Outcome::fail);
    REQUIRE(q.parts.front().where == "f = 0 Hz, Y(substrate, gate)");

    // A NaN admittance fails, it does not pass.
    auto nan = s.zero;
    nan.points[0].admittance[1][1 * nc + 2] = {std::numeric_limits<double>::quiet_NaN(), 0.0};
    REQUIRE(c::check_c7(nan, 0, s.device.mesh(), s.data).comparison.outcome == Outcome::fail);

    // A symmetric change is not an asymmetry: still a pass.
    auto sym = s.zero;
    const std::complex<double> dy(5.0 * (res[0] + res[2]), 0.0);
    sym.points[0].admittance[k][2 * nc + 0] += dy;
    sym.points[0].admittance[k][0 * nc + 2] += dy;
    REQUIRE(c::check_c7(sym, 0, s.device.mesh(), s.data).comparison.outcome == Outcome::pass);
}

TEST_CASE("C7 is NOT EVALUATED without its inputs", "[verify][consistency]") {
    const Mos& s = mos();
    const auto& m = s.device.mesh();
    REQUIRE(c::check_c7(s.zero, 1, m, s.data).comparison.outcome == Outcome::not_evaluated);

    auto no_res = s.zero;
    no_res.points[0].resolution.clear();
    const auto a = c::check_c7(no_res, 0, m, s.data);
    REQUIRE(a.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(a.parts.front().reason.find("not one per frequency") != std::string::npos);

    auto short_y = s.zero;
    short_y.points[0].admittance[2].pop_back();
    const auto b = c::check_c7(short_y, 0, m, s.data);
    REQUIRE(b.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(b.parts.front().reason.find("frequency 2") != std::string::npos);

    auto no_y = s.zero;
    no_y.points[0].admittance.clear();
    const auto e = c::check_c7(no_y, 0, m, s.data);
    REQUIRE(e.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(e.parts.front().reason.find("not one per frequency") != std::string::npos);

    auto fewer = s.data;
    fewer.contacts.pop_back();
    const auto f = c::check_c7(s.zero, 0, m, fewer);
    REQUIRE(f.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(f.parts.front().reason.find("contacts are not the device's") != std::string::npos);
}

TEST_CASE("C7's bound is the sum of the two columns' resolutions", "[verify][consistency]") {
    // On a symmetrized copy, an asymmetry of 0.99 (res_i + res_j) passes and 1.01 fails, at the pair and
    // frequency whose two resolutions are most nearly equal (so neither alone is the bound).
    const Mos& s = mos();
    const std::size_t nc = s.zero.contacts;
    auto sym = s.zero;
    std::size_t kb = 0, ib = 1, jb = 0;
    double share = 0.0;
    for (std::size_t k = 0; k < frequencies.size(); ++k) {
        auto& y = sym.points[0].admittance[k];
        const auto& res = sym.points[0].resolution[k];
        for (std::size_t i = 0; i < nc; ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                y[j * nc + i] = y[i * nc + j];
                const double q = std::min(res[i], res[j]) / (res[i] + res[j]);
                if (q > share) {
                    share = q;
                    kb = k;
                    ib = i;
                    jb = j;
                }
            }
        }
    }
    CAPTURE(share, kb, ib, jb);
    REQUIRE(share > 0.05);
    const auto& res = sym.points[0].resolution[kb];
    const double bound = res[ib] + res[jb];
    for (const double factor : {0.99, 1.01}) {
        auto copy = sym;
        copy.points[0].admittance[kb][ib * nc + jb] += std::complex<double>(0.6 * factor * bound, 0.8 * factor * bound);
        const auto r = c::check_c7(copy, 0, s.device.mesh(), s.data);
        CAPTURE(factor, r.parts.front().residual, r.parts.front().tolerance);
        REQUIRE(r.comparison.outcome == (factor < 1.0 ? Outcome::pass : Outcome::fail));
    }
}

TEST_CASE("C7 on the three TVR-1 decks with a zero-bias AC point", "[verify][consistency]") {
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
        const auto r = c::check_c7(*run, 0, plan->device.mesh(), data);
        CAPTURE(r.comparison.reason, r.parts.front().residual, r.parts.front().tolerance, r.parts.front().where,
                worst_ratio(*run));
        REQUIRE(r.comparison.outcome == Outcome::pass);
        REQUIRE(r.diagnostics.size() == 1);
        REQUIRE(r.diagnostics.front().name == "two contacts");
    }
}
