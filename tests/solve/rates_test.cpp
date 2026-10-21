// Generation-recombination rates in results (Unit 45, RG1, RG2): at a converged steady state the
// net current out of each node's control volume is q U V (electrons) and -q U V (holes), so the
// stored rates and the stored edge currents, two outputs of one solve, check each other; at
// equilibrium every rate vanishes; each mechanism's vector is there exactly when its model is on;
// the rates can be switched off, and change nothing else.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/transient.hpp"

using namespace NiTCAD;

namespace {

std::vector<double> axis(double length, int nodes) {
    std::vector<double> a;
    for (int i = 0; i < nodes; ++i) a.push_back(length * i / (nodes - 1));
    return a;
}

// A 2D pn diode, 2 x 0.5 um: p (1e17) for x < 1 um, n (1e17) beyond; anode on x_min; optionally
// isothermal heat sinks at 300 K on both ends (the electrothermal model).
device::Device diode2d(bool heat_sinks = false) {
    mesh::Mesh m = *mesh::make_tensor_grid(axis(2e-4, 41), axis(0.5e-4, 6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::move(donors),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                 {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
                                    .thermal_contacts = heat_sinks ? std::vector<device::ThermalContact>{
                                        {"left", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
                                        {"right", "x_max", device::ThermalContactKind::isothermal, 300.0, 0.0}}
                                                                   : std::vector<device::ThermalContact>{}});
}

std::vector<double> axis_from(double a, double b, int nodes) {
    std::vector<double> v;
    for (int i = 0; i < nodes; ++i) v.push_back(a + (b - a) * i / (nodes - 1));
    return v;
}

// The small-signal tests' oxide-silicon device (small_signal_test.cpp mos2d): an electrode on a meshed
// oxide, an n+ source pocket and the substrate, with interface traps, fixed charge and surface
// recombination on the oxide-silicon interface.
device::Device mos2d() {
    const auto x = *mesh::straddle_interface(axis_from(-5e-7, 1.2e-5, 30), 0.0, 1e-7);
    const auto t = axis_from(0.0, 3e-6, 4);
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
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        if (m.points()[static_cast<std::size_t>(v)][1] < 1.5e-6) electrode.push_back(v);
    }
    for (const mesh::NodeId v : m.find_boundary("y_max")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        if (px > 0.0 && px < 3e-6) source.push_back(v);
    }
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = 3e11;
    f.traps.levels = {{.type = physics::TrapType::donor, .density_cm2 = 2e11, .energy_eV = -0.2,
                       .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 3e-16},
                      {.type = physics::TrapType::acceptor, .density_cm2 = 1e11,
                       .energy_eV = 0.25, .sigma_n_cm2 = 2e-16, .sigma_p_cm2 = 1e-15}};
    f.traps.bands = {{.type = physics::TrapType::acceptor, .density_cm2_eV = 5e11,
                      .energy_low_eV = -0.1, .energy_high_eV = 0.4}};
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

// Per node, the conventional current out of its control volume along its edges (A / cm^(3-D)).
std::vector<double> outflow(const mesh::Mesh& m, const std::vector<double>& edge_current) {
    std::vector<double> out(m.node_count(), 0.0);
    const auto edges = m.edges();
    for (std::size_t e = 0; e < edges.size(); ++e) {
        out[static_cast<std::size_t>(edges[e].first)] += edge_current[e];
        out[static_cast<std::size_t>(edges[e].second)] -= edge_current[e];
    }
    return out;
}

bool on_contact(const device::Device& d, std::size_t i) {
    for (const auto& c : d.contacts()) {
        if (std::ranges::binary_search(c.nodes, static_cast<mesh::NodeId>(i))) return true;
    }
    return false;
}

// Checks div J = q U against the edge currents at every node off the contacts: each carrier's miss
// within `tolerance` of the larger of the currents through the node and q U V, plus 5% of the run's
// current resolution (BiasPoint::terminal_current_resolution, Unit 15): a node's net current below
// it is the rounding of large opposing fluxes (majority carriers in reverse bias), not a measure of
// U. Reports the worst ratio of miss to that bound (pass: at most 1).
void check_balance(const device::Device& d, const results::BiasPoint& p, double tolerance) {
    const auto& m = d.mesh();
    const auto vol = m.volumes();
    const std::vector<double> out_n = outflow(m, p.edge_current_n);
    const std::vector<double> out_p = outflow(m, p.edge_current_p);
    std::vector<double> scale_n(m.node_count(), 0.0), scale_p(m.node_count(), 0.0);
    for (std::size_t e = 0; e < m.edges().size(); ++e) {
        for (const auto v : {m.edges()[e].first, m.edges()[e].second}) {
            scale_n[static_cast<std::size_t>(v)] += std::abs(p.edge_current_n[e]);
            scale_p[static_cast<std::size_t>(v)] += std::abs(p.edge_current_p[e]);
        }
    }
    const double unresolved = 0.05 * *std::ranges::max_element(p.terminal_current_resolution);
    double worst = 0.0, worst_relative = 0.0;
    for (std::size_t i = 0; i < m.node_count(); ++i) {
        if (on_contact(d, i)) continue;
        const double qUV_n = base::q_C * p.rates.net_n[i] * vol[i];
        const double qUV_p = base::q_C * p.rates.net_p[i] * vol[i];
        const double miss_n = std::abs(out_n[i] - qUV_n), miss_p = std::abs(out_p[i] + qUV_p);
        const double sn = std::max(scale_n[i], std::abs(qUV_n)), sp = std::max(scale_p[i], std::abs(qUV_p));
        worst = std::max({worst, miss_n / (tolerance * sn + unresolved), miss_p / (tolerance * sp + unresolved)});
        if (miss_n > unresolved) worst_relative = std::max(worst_relative, miss_n / sn);
        if (miss_p > unresolved) worst_relative = std::max(worst_relative, miss_p / sp);
    }
    CAPTURE(worst, worst_relative, unresolved);
    CHECK(worst <= 1.0);
}

}  // namespace

TEST_CASE("rates: div J = q U at every node of a forward-biased diode (Unit 45, RG1)") {
    const device::Device d = diode2d();
    const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.3, 0.0}, {0.6, 0.0}};
    const auto run = solve::sweep_bias(d, points);
    REQUIRE(run.has_value());
    REQUIRE(!run->stopped);
    const results::BiasPoint& p = run->points.back();
    const results::RateFields& r = p.rates;
    const std::size_t n = d.mesh().node_count();
    // SRH, Auger and radiative on by default; no interfaces, impact or tunnelling.
    CHECK(r.srh.size() == n);
    CHECK(r.auger.size() == n);
    CHECK(r.radiative.size() == n);
    CHECK(r.impact.empty());
    CHECK(r.interface_n.empty());
    CHECK(r.band_to_band_n.empty());
    REQUIRE(r.net_n.size() == n);
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(r.net_n[i] == r.srh[i] + r.auger[i] + r.radiative[i]);
        CHECK(r.net_p[i] == r.net_n[i]);  // these mechanisms change both carriers alike
    }
    // Forward bias: net recombination in the junction, as large as some 1e20 cm^-3 s^-1.
    CHECK(*std::ranges::max_element(r.srh) > 1e18);
    check_balance(d, p, 1e-6);
}

TEST_CASE("rates: every rate vanishes at equilibrium (Unit 45, RG2)") {
    const device::Device d = diode2d();
    const auto run = solve::sweep_bias(d, std::vector<std::vector<double>>{{0.0, 0.0}});
    REQUIRE(run.has_value());
    const results::RateFields& r = run->points.front().rates;
    // Against the scale n_i / tau (about 1e10 / 1e-6 = 1e16 cm^-3 s^-1): rounding only.
    for (const auto* v : {&r.srh, &r.auger, &r.radiative, &r.net_n, &r.net_p}) {
        for (const double x : *v) CHECK(std::abs(x) <= 1e4);
    }
}

TEST_CASE("rates: impact ionization is a generation, and the balance holds with it (Unit 45)") {
    const device::Device d = diode2d();
    solve::BiasOptions o;
    o.models.impact_ionization = true;
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 8; ++k) points.push_back({-0.5 * k, 0.0});
    const auto run = solve::sweep_bias(d, points, o);
    REQUIRE(run.has_value());
    REQUIRE(!run->stopped);
    const results::BiasPoint& p = run->points.back();
    REQUIRE(p.rates.impact.size() == d.mesh().node_count());
    CHECK(*std::ranges::min_element(p.rates.impact) < 0.0);  // generation is negative R - G
    CHECK(*std::ranges::max_element(p.rates.impact) <= 0.0);
    check_balance(d, p, 1e-6);
}

TEST_CASE("rates: off by option, and nothing else changes (Unit 45, RQ1)") {
    const device::Device d = diode2d();
    const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.4, 0.0}};
    solve::BiasOptions off;
    off.store_rates = false;
    const auto a = solve::sweep_bias(d, points);
    const auto b = solve::sweep_bias(d, points, off);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    CHECK(b->points.back().rates.net_n.empty());
    CHECK(b->points.back().rates.srh.empty());
    for (std::size_t k = 0; k < points.size(); ++k) {
        CHECK(a->points[k].fields.n_cm3 == b->points[k].fields.n_cm3);
        CHECK(a->points[k].terminal_current == b->points[k].terminal_current);
    }
    // A run without rates is another output: its identity differs; with them it is the default's.
    CHECK(a->run.input_identity != b->run.input_identity);
    CHECK(solve::make_run_record(d, solve::BiasOptions{}, points).input_identity == a->run.input_identity);
}

TEST_CASE("rates: transient snapshots carry them (Unit 45)") {
    const device::Device d = diode2d();
    solve::TransientOptions o{.t_end_s = 1e-10, .dt_initial_s = 5e-11, .adaptive = false};
    const std::vector<solve::Waveform> wf{*solve::Waveform::step(0.0, 0.4, 0.0), solve::Waveform::constant(0.0)};
    const auto run = solve::solve_transient(d, wf, o);
    REQUIRE(run.has_value());
    REQUIRE(!run->snapshots.empty());
    CHECK(run->snapshots.back().rates.net_n.size() == d.mesh().node_count());
    o.steady.store_rates = false;
    const auto off = solve::solve_transient(d, wf, o);
    REQUIRE(off.has_value());
    CHECK(off->snapshots.back().rates.net_n.empty());
}

TEST_CASE("rates: interfaces and band-to-band tunnelling, with the balance (Unit 45, RG1)") {
    // Interface traps and surface recombination: per-carrier interface rates on the interface
    // nodes only, and div J = q U everywhere.
    const device::Device mos = mos2d();
    const auto run = solve::sweep_bias(mos, std::vector<std::vector<double>>{{0.0, 0.0, 0.0}, {-0.6, 0.0, 0.3}});
    REQUIRE(run.has_value());
    REQUIRE(!run->stopped);
    const results::BiasPoint& p = run->points.back();
    REQUIRE(p.rates.interface_n.size() == mos.mesh().node_count());
    REQUIRE(p.rates.interface_p.size() == mos.mesh().node_count());
    bool some = false;
    for (std::size_t i = 0; i < mos.mesh().node_count(); ++i) {
        const double x = mos.mesh().points()[i][0];
        if (std::abs(x) > 2e-7) {
            CHECK(p.rates.interface_n[i] == 0.0);  // away from the interface
        } else if (p.rates.interface_n[i] != 0.0) {
            some = true;
        }
    }
    CHECK(some);
    check_balance(mos, p, 1e-6);

    // Local band-to-band tunnelling in a reverse-biased 1e19 / 1e19 junction.
    mesh::Mesh m = *mesh::make_tensor_grid(axis(0.2e-4, 81), axis(0.05e-4, 3));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 0.1e-4 ? acceptors : donors)[i] = 1e19;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    solve::BiasOptions o;
    o.models.btbt_local = true;
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 4; ++k) points.push_back({-0.25 * k, 0.0});
    const auto btbt = solve::sweep_bias(d, points, o);
    REQUIRE(btbt.has_value());
    REQUIRE(!btbt->stopped);
    const results::BiasPoint& b = btbt->points.back();
    REQUIRE(b.rates.band_to_band_n.size() == n);
    CHECK(*std::ranges::min_element(b.rates.band_to_band_n) < 0.0);  // a generation
    CHECK(b.rates.band_to_band_n == b.rates.band_to_band_p);       // local: pairs at one node
    check_balance(d, b, 1e-6);
}

TEST_CASE("rates: the electrothermal terms, with the balance (Unit 45, RG1)") {
    const device::Device d = diode2d(true);
    solve::BiasOptions o;
    o.models.electrothermal = true;
    const auto run = solve::sweep_bias(d, std::vector<std::vector<double>>{{0.0, 0.0}, {0.3, 0.0}, {0.6, 0.0}}, o);
    REQUIRE(run.has_value());
    REQUIRE(!run->stopped);
    const results::BiasPoint& p = run->points.back();
    CHECK(p.rates.srh.size() == d.mesh().node_count());
    CHECK(*std::ranges::max_element(p.rates.srh) > 1e18);
    check_balance(d, p, 1e-6);
}

TEST_CASE("rates: what storing them costs (Unit 45, RG5)", "[.benchmark]") {
    // Sweep time with and without the rates (the best of three runs each) on a 2D diode of 241 x 61
    // nodes and a 3D diode of 41 x 11 x 11, 0 to 0.6 V in 0.1 V steps.
    using clock = std::chrono::steady_clock;
    const auto pn = [](mesh::Mesh m) {
        const std::size_t n = m.node_count();
        std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
        auto anode = m.find_boundary("x_min")->nodes;
        auto cathode = m.find_boundary("x_max")->nodes;
        return *device::Device::create({.mesh = std::move(m),
                                        .temperature_K = 300.0,
                                        .regions = {{"silicon", physics::silicon()}},
                                        .node_region = std::vector<device::RegionId>(n, 0),
                                        .donors = std::move(donors),
                                        .acceptors = std::move(acceptors),
                                        .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                     {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    };
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 6; ++k) points.push_back({0.1 * k, 0.0});
    for (const bool three_d : {false, true}) {
        const device::Device d = three_d ? pn(*mesh::make_tensor_grid(axis(2e-4, 41), axis(0.5e-4, 11), axis(0.5e-4, 11)))
                                         : pn(*mesh::make_tensor_grid(axis(2e-4, 241), axis(0.5e-4, 61)));
        double best[2] = {1e300, 1e300};
        for (int run = 0; run < 3; ++run) {
            for (const bool store : {true, false}) {
                solve::BiasOptions o;
                o.store_rates = store;
                const auto t0 = clock::now();
                const auto s = solve::sweep_bias(d, points, o);
                const double t = std::chrono::duration<double>(clock::now() - t0).count();
                REQUIRE(s.has_value());
                best[store ? 0 : 1] = std::min(best[store ? 0 : 1], t);
            }
        }
        std::printf("rates cost, %s diode (%zu nodes, 7 points): with %.3f s, without %.3f s (%+.1f%%)\n",
                    three_d ? "3D" : "2D", d.mesh().node_count(), best[0], best[1], 100.0 * (best[0] / best[1] - 1.0));
    }
}

TEST_CASE("rates: nonlocal tunnelling puts electrons and holes apart; the balance holds (Unit 45)") {
    // The BTBT tests' 1D Zener junction (btbt_test.cpp zener: 100 nm, 5e19 each side, its cells),
    // reverse biased to 1.5 V with the nonlocal Kane model: the electron and hole rates differ (each
    // path's holes at its start, its electrons where it ends), and div J = q U holds per carrier.
    const auto x = axis(1e-5, 201);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    auto cells = *mesh::TensorCells::create(m, x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 0.5e-5 ? acceptors : donors)[i] = 5e19;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"semiconductor", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .cells = std::move(cells)});
    solve::BiasOptions o;
    o.models.btbt_nonlocal = assemble::NonlocalTunnelling::kane;
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 6; ++k) points.push_back({-0.25 * k, 0.0});
    const auto run = solve::sweep_bias(d, points, o);
    REQUIRE(run.has_value());
    REQUIRE(!run->stopped);
    const results::BiasPoint& p = run->points.back();
    REQUIRE(p.rates.band_to_band_n.size() == n);
    CHECK(p.rates.band_to_band_n != p.rates.band_to_band_p);
    CHECK(p.rates.net_n != p.rates.net_p);
    check_balance(d, p, 1e-6);
}
