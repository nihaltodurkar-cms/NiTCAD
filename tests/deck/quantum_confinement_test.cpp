// MLDA quantum confinement on the example devices (Unit 47; DECISIONS.md QC1, QC5, QC6, QC10,
// QC12): the run record names the model only when on; and, Release only (hidden tag [.mlda]): the
// 2D MOSFET example at equilibrium with the model on, Newton's order in strong inversion, its Id-Vg
// reflected (x and y exchanged), the small-signal conductances against DC differences, and the
// reported times, currents and the FinFET deck's distance search.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/quantum_confinement.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

using deck_test::must;

namespace {

const auto mlda = assemble::QuantumConfinement::mlda;

// A plan's sweep with quantum confinement set, no files written, and its wall time.
std::pair<results::Sweep, double> sweep(deck::Plan plan, assemble::QuantumConfinement q) {
    plan.options.models.quantum_confinement = q;
    plan.outputs = {};
    const auto t0 = std::chrono::steady_clock::now();
    const auto report = deck::run(plan);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    REQUIRE_FALSE(report->sweep->stopped.has_value());
    return {std::move(*report->sweep), seconds};
}

deck::Plan example(const fs::path& file) {
    const auto d = deck::read(file);
    if (!d) FAIL(d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

// The 2D MOSFET example with x and y exchanged (as surface_mobility_test.cpp's).
const std::string reflected_mosfet = R"(
[mesh]
type = tensor
x = linspace(0, 0.9um, 19) + linspace(0.9um, 0.95um, 6) + linspace(0.95um, 1um, 3)
y = linspace(0, 2um, 41)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0.95um 1um  0 2um

[doping]
profile = acceptor constant 1e16
profile = donor gaussian 2e19 axis=x center=0.95um width=0.15um y=0:0.5um lateral=20nm
profile = donor gaussian 2e19 axis=x center=0.95um width=0.15um y=1.5um:2um lateral=20nm

[contact source]
kind = ohmic
patch = y_min
box = 0.6um 0.94um  0 0

[contact drain]
kind = ohmic
patch = y_max
box = 0.6um 0.94um  2um 2um

[contact gate]
kind = electrode
patch = x_max
box = 1um 1um  0.5um 1.5um
electrode = metal
work_function = 4.1

[contact body]
kind = ohmic
patch = x_min

[sweep]
contact gate = 0 : 1.5 : 0.1
contact drain = 0.1
)";

}  // namespace

TEST_CASE("quantum confinement: the run record names the model and its masses only when on "
          "(QC1)",
          "[deck]") {
    deck::Plan p = example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck");
    const auto off = solve::make_run_record(p.device, p.options, p.points);
    p.options.models.quantum_confinement = mlda;
    const auto on = solve::make_run_record(p.device, p.options, p.points);
    p.options.models.mlda_mass_n = 0.916;
    const auto other_mass = solve::make_run_record(p.device, p.options, p.points);
    const auto named = [](const results::RunRecord& r, const char* key) {
        return std::ranges::any_of(r.settings, [&](const auto& s) { return s.first == key; });
    };
    for (const char* key :
         {"models.quantum_confinement", "models.mlda_mass_n", "models.mlda_mass_p"}) {
        CHECK(named(on, key));
        CHECK_FALSE(named(off, key));
    }
    CHECK(on.input_identity != off.input_identity);
    CHECK(other_mass.input_identity != on.input_identity);
    // A mass changed with the model off is not read, so the identity is unchanged.
    deck::Plan q = example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck");
    q.options.models.mlda_mass_p = 0.3;
    CHECK(solve::make_run_record(q.device, q.options, q.points).input_identity ==
          off.input_identity);
}

TEST_CASE("quantum confinement: the deck keys", "[deck]") {
    const auto with = [](const std::string& models) {
        return deck_test::diode_device + "[models]\n" + models + "\n[sweep]\ncontact anode = 0\n";
    };
    const deck::Plan p = must(with("quantum_confinement = mlda\nmlda_mass_n = 0.916\nmlda_mass_p = 0.3"));
    CHECK(p.options.models.quantum_confinement == mlda);
    CHECK(p.options.models.mlda_mass_n == 0.916);
    CHECK(p.options.models.mlda_mass_p == 0.3);
    const deck::Plan q = must(with("quantum_confinement = off"));
    CHECK(q.options.models.quantum_confinement == assemble::QuantumConfinement::off);
    CHECK(q.options.models.mlda_mass_n == 0.98);
    // Unit 62: density_gradient is a value now, with its own keys (QP1).
    const deck::Plan d = must(with("quantum_confinement = density_gradient\ndg_mass_n = 0.916\ndg_mass_p = 0.3\n"
                                   "dg_gamma_n = 0.5\ndg_gamma_p = 2"));
    CHECK(d.options.models.quantum_confinement == assemble::QuantumConfinement::density_gradient);
    CHECK(d.options.models.dg_mass_n == 0.916);
    CHECK(d.options.models.dg_mass_p == 0.3);
    CHECK(d.options.models.dg_gamma_n == 0.5);
    CHECK(d.options.models.dg_gamma_p == 2.0);
    CHECK(q.options.models.dg_mass_n == 0.98);
    CHECK(q.options.models.dg_mass_p == 0.49);
    CHECK(q.options.models.dg_gamma_n == 1.0);
    const auto bad = deck::parse(with("quantum_confinement = schrodinger"), "bad.deck");
    REQUIRE(bad.has_value());
    const auto built = deck::build(*bad);
    REQUIRE_FALSE(built.has_value());
    CHECK(built.error().message.find("quantum_confinement is off, mlda or density_gradient") != std::string::npos);
}

TEST_CASE("quantum confinement: the MOSFET examples, reflected, small signal, reports "
          "(QC5, QC6, QC10, QC12)",
          "[.mlda]") {
    const deck::Plan p = example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck");
    {
        const auto c = assemble::mlda_factors(p.device, [&] {
            auto m = p.options.models;
            m.quantum_confinement = mlda;
            return m;
        }());
        REQUIRE(c.has_value());
        double nearest = 1.0;
        for (std::size_t i = 0; i < c->distance_cm.size(); ++i) {
            if (!p.device.is_insulator(static_cast<mesh::NodeId>(i))) {
                nearest = std::min(nearest, c->distance_cm[i]);
            }
        }
        std::printf("QC12 2D MOSFET example: nearest silicon node %.2f nm from the oxide\n",
                    nearest * 1e7);
    }
    const auto [off, t_off] = sweep(p, assemble::QuantumConfinement::off);
    const auto [on, t_on] = sweep(p, mlda);
    const auto id = [](const results::Sweep& s, std::size_t k) {
        return s.points[k].terminal_current[1];
    };
    const std::size_t top = 15;  // Vg 1.5 V
    REQUIRE(on.points.size() == 16);
    std::printf("QC12 2D MOSFET example: %.2f s off, %.2f s with MLDA; Id(1.5 V) %.6e A/cm with, "
                "%.6e without (ratio %.6f)\n",
                t_off, t_on, id(on, top), id(off, top), id(on, top) / id(off, top));

    // QC6: Newton at Vg = 1.5 V converges quadratically over its last three corrections.
    {
        const auto& steps = on.points[top].convergence.iterations;
        for (const auto& r : steps) {
            std::printf("QC6 iteration %d: correction %.3e\n", r.iteration, r.update);
        }
        REQUIRE(steps.size() >= 3);
        const double u1 = steps[steps.size() - 3].update, u2 = steps[steps.size() - 2].update,
                     u3 = steps.back().update;
        const double order = std::log(u3 / u2) / std::log(u2 / u1);
        std::printf("QC6 convergence order %.2f\n", order);
        CHECK(order >= 1.8);
    }

    // QC5: equilibrium at Vg = 1 V (source, drain, body at 0) carries no current.
    {
        auto options = p.options;
        options.models.quantum_confinement = mlda;
        const auto eq = solve::solve_bias(p.device, std::vector<double>{0.0, 0.0, 1.0, 0.0}, options);
        if (!eq) FAIL(eq.error().message);
        for (std::size_t c = 0; c < 4; ++c) {
            CAPTURE(c, eq->terminal_current[c], eq->terminal_current_resolution[c]);
            CHECK(std::abs(eq->terminal_current[c]) <= eq->terminal_current_resolution[c]);
        }
    }

    // QC10: the reflected device gives the same currents.
    const auto [reflected, t_reflected] = sweep(must(reflected_mosfet), mlda);
    for (const std::size_t k : {std::size_t{5}, top}) {
        const double a = id(on, k), b = id(reflected, k);
        CAPTURE(k, a, b);
        CHECK(std::abs(b / a - 1.0) <= 1e-9);
    }

    // QC10: the small-signal drain conductance and transconductance against DC differences.
    solve::SmallSignalOptions ac;
    ac.steady = p.options;
    ac.steady.models.quantum_confinement = mlda;
    ac.frequencies_Hz = {0.0};
    const auto y = solve::solve_small_signal(p.device, p.points, ac);
    REQUIRE(y.has_value());
    REQUIRE_FALSE(y->stopped.has_value());
    const auto& Y = y->points[top].admittance[0];
    const std::size_t contacts = 4;
    const double h = 1e-4;
    for (const std::size_t j : {std::size_t{1}, std::size_t{2}}) {  // drain, gate
        deck::Plan q = p;
        auto up = p.points[top], down = p.points[top];
        up[j] += h;
        down[j] -= h;
        q.points.push_back(up);
        q.points.push_back(down);
        const auto [dc, t_dc] = sweep(q, mlda);
        const double g = (id(dc, top + 1) - id(dc, top + 2)) / (2.0 * h);
        const double ac_g = Y[1 * contacts + j].real();
        std::printf("QC10 dI_drain/dV_%zu: small signal %.8e, DC %.8e\n", j, ac_g, g);
        CHECK(std::abs(ac_g / g - 1.0) <= 1e-5);
    }

    // QC12: the 3D example.
    const deck::Plan p3 = example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_3d.deck");
    const auto [off3, t_off3] = sweep(p3, assemble::QuantumConfinement::off);
    const auto [on3, t_on3] = sweep(p3, mlda);
    std::printf("QC12 3D MOSFET example: %.2f s off, %.2f s with MLDA; Id(1.2 V) ratio %.6f\n",
                t_off3, t_on3,
                on3.points.back().terminal_current[1] / off3.points.back().terminal_current[1]);

    // QC12: the distance search on the FinFET deck's device (45,375 nodes).
    const deck::Plan fin = example(fs::path(NITCAD_DECKS_DIR) / "finfet_3d.deck");
    auto m = fin.options.models;
    m.quantum_confinement = mlda;
    const auto t0 = std::chrono::steady_clock::now();
    const auto c = assemble::mlda_factors(fin.device, m);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    REQUIRE(c.has_value());
    double nearest = 1.0;
    for (std::size_t i = 0; i < c->distance_cm.size(); ++i) {
        if (!fin.device.is_insulator(static_cast<mesh::NodeId>(i))) {
            nearest = std::min(nearest, c->distance_cm[i]);
        }
    }
    std::printf("QC12 FinFET deck: %zu nodes, distance search %.3f s, nearest silicon node %.3f nm\n",
                fin.device.mesh().node_count(), seconds, nearest * 1e7);
}

TEST_CASE("density gradient: the run record names the model, its masses and gammas only when on (Unit 62, QP1)",
          "[deck]") {
    deck::Plan p = example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck");
    const auto off = solve::make_run_record(p.device, p.options, p.points);
    p.options.models.quantum_confinement = mlda;
    const auto with_mlda = solve::make_run_record(p.device, p.options, p.points);
    p.options.models.quantum_confinement = assemble::QuantumConfinement::density_gradient;
    const auto on = solve::make_run_record(p.device, p.options, p.points);
    p.options.models.dg_gamma_n = 0.5;
    const auto other_gamma = solve::make_run_record(p.device, p.options, p.points);
    const auto named = [](const results::RunRecord& r, const char* key) {
        return std::ranges::any_of(r.settings, [&](const auto& s) { return s.first == key; });
    };
    for (const char* key : {"models.dg_mass_n", "models.dg_mass_p", "models.dg_gamma_n", "models.dg_gamma_p"}) {
        CHECK(named(on, key));
        CHECK_FALSE(named(off, key));
        CHECK_FALSE(named(with_mlda, key));
    }
    CHECK(named(on, "models.quantum_confinement"));
    CHECK_FALSE(named(on, "models.mlda_mass_n"));
    CHECK(on.input_identity != off.input_identity);
    CHECK(on.input_identity != with_mlda.input_identity);
    CHECK(other_gamma.input_identity != on.input_identity);
    // DG's parameters changed with MLDA on or the model off are not read: the identities are unchanged.
    deck::Plan q = example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck");
    q.options.models.dg_mass_n = 0.3;
    CHECK(solve::make_run_record(q.device, q.options, q.points).input_identity == off.input_identity);
    q.options.models.quantum_confinement = mlda;
    CHECK(solve::make_run_record(q.device, q.options, q.points).input_identity == with_mlda.input_identity);
}

// Unit 62 (QV4, QV9, QV13 reported): density gradient on the 2D MOSFET example. Hidden as the MLDA
// case above (the sweeps take minutes); run by name.
TEST_CASE("density gradient: the 2D MOSFET example, equilibrium and reflected (QV4, QV9)", "[.dg]") {
    const auto dg = assemble::QuantumConfinement::density_gradient;
    deck::Plan p = example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck");
    // QV4: equilibrium at Vg = 1 V (source, drain, body at 0): no current, flat quasi-Fermi levels.
    {
        auto options = p.options;
        options.models.quantum_confinement = dg;
        const auto t0 = std::chrono::steady_clock::now();
        const auto eq = solve::solve_bias(p.device, std::vector<double>{0.0, 0.0, 1.0, 0.0}, options);
        if (!eq) FAIL(eq.error().message);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        for (std::size_t c = 0; c < 4; ++c) {
            CAPTURE(c, eq->terminal_current[c], eq->terminal_current_resolution[c]);
            CHECK(std::abs(eq->terminal_current[c]) <= eq->terminal_current_resolution[c]);
        }
        double flat = 0.0, peak_n = 0.0, far = 0.0;
        for (std::size_t i = 0; i < p.device.mesh().node_count(); ++i) {
            if (p.device.is_insulator(static_cast<mesh::NodeId>(i))) continue;
            flat = std::max({flat, std::abs(eq->bands.electron_fermi_eV[i]), std::abs(eq->bands.hole_fermi_eV[i])});
            peak_n = std::max(peak_n, eq->fields.quantum_potential_n_eV[i]);
            far = std::max({far, std::abs(eq->fields.quantum_potential_n_eV[i]),
                            std::abs(eq->fields.quantum_potential_p_eV[i])});
        }
        std::printf("QV4 2D MOSFET example at Vg 1 V: %.1f s, %zu iterations, E_F spread %.2e eV, "
                    "peak Lambda_n %.4f eV, largest |Lambda| %.4f eV\n",
                    seconds, eq->convergence.iterations.size(), flat, peak_n, far);
        CHECK(flat <= 1e-9);
    }
    // QV9: the reflected device gives the same drain current. Three points of the sweep.
    p.points = {{0.0, 0.1, 0.5, 0.0}, {0.0, 0.1, 1.0, 0.0}, {0.0, 0.1, 1.5, 0.0}};
    auto r = must(reflected_mosfet);
    r.points = p.points;
    const auto [a, ta] = sweep(p, dg);
    const auto [b, tb] = sweep(std::move(r), dg);
    REQUIRE(a.points.size() == 3);
    REQUIRE(b.points.size() == 3);
    for (std::size_t k = 0; k < 3; ++k) {
        const double ia = a.points[k].terminal_current[1], ib = b.points[k].terminal_current[1];
        std::printf("QV9 Vg %.1f V: Id %.9e and reflected %.9e A/cm (%.2e)\n", p.points[k][2], ia, ib, ib / ia - 1.0);
        CHECK(std::abs(ib / ia - 1.0) <= 1e-9);
    }
    std::printf("QV13 2D MOSFET example, 3 points: %.1f s and %.1f s reflected\n", ta, tb);
}

// The lumped-gate MOS-C of tests/decks/lumped_gate_1d.deck (0.5 nm surface spacing), no analysis.
const std::string mos_device = R"(
[device]
temperature = 300
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
[contact body]
kind = ohmic
patch = x_max
)";

TEST_CASE("density gradient: results files, plot fields and a restart (Unit 62, QV12)", "[deck]") {
    const fs::path dir = deck_test::scratch("dg_restart");
    const std::string models = "[models]\nquantum_confinement = density_gradient\n";
    const std::string device = mos_device + models;
    const auto run = [&](const std::string& text) { return deck::run(must(text, dir)); };
    // The uninterrupted sweep, its file read back with the quantum potentials.
    const deck::Plan plan = must(device + "[sweep]\ncontact gate = 0 : -1.2 : -0.1\n[output]\nsweep = full.nitsweep\n", dir);
    const auto full = deck::run(plan);
    if (!full) FAIL(full.error().message);
    REQUIRE(full->sweep.has_value());
    const auto loaded = io::load_sweep(dir / "full.nitsweep");
    if (!loaded) FAIL(loaded.error().message);
    if (full->sweep->stopped) FAIL(full->sweep->stopped->message);
    REQUIRE(loaded->points.size() == 13);
    for (std::size_t k = 0; k < 13; ++k) {
        CHECK(loaded->points[k].fields.quantum_potential_n_eV == full->sweep->points[k].fields.quantum_potential_n_eV);
        CHECK(loaded->points[k].fields.quantum_potential_p_eV == full->sweep->points[k].fields.quantum_potential_p_eV);
    }
    const auto& last = full->sweep->points.back();
    REQUIRE(last.fields.quantum_potential_n_eV.size() == plan.device.mesh().node_count());
    double peak = 0.0;
    for (const double l : last.fields.quantum_potential_p_eV) peak = std::max(peak, std::abs(l));
    CAPTURE(peak);
    CHECK(peak > 1e-3);  // accumulation: Lambda_p resolved
    // The plot fields: Lambda and the effective band edges.
    const auto lambda = render::field_values(plan.device, last, render::Field::quantum_potential_n);
    REQUIRE(lambda.has_value());
    CHECK(*lambda == last.fields.quantum_potential_n_eV);
    const auto ec = render::field_values(plan.device, last, render::Field::effective_conduction_band);
    const auto ev = render::field_values(plan.device, last, render::Field::effective_valence_band);
    REQUIRE(ec.has_value());
    REQUIRE(ev.has_value());
    for (std::size_t i = 0; i < plan.device.mesh().node_count(); ++i) {
        CHECK((*ec)[i] == last.bands.conduction_eV[i] + last.fields.quantum_potential_n_eV[i]);
        CHECK((*ev)[i] == last.bands.valence_eV[i] - last.fields.quantum_potential_p_eV[i]);
    }
    CHECK(render::find_field("Ln") == render::Field::quantum_potential_n);
    CHECK(render::find_field("Ev_eff") == render::Field::effective_valence_band);

    // A checkpoint at -0.6 V; the restart's points equal the uninterrupted sweep's to 1e-9 relative
    // (as Unit 61's steady restart, DK9; a steady checkpoint is fields in physical units, not the
    // solver state, so a sweep restart is not bit for bit: the transient's is, transient_restart_test).
    const auto first = run(device + "[sweep]\ncontact gate = 0 : -0.6 : -0.1\n[output]\ncheckpoint = at.nitckpt\n");
    if (!first) FAIL(first.error().message);
    const auto resumed = run(device + "[sweep]\ncontact gate = 0 : -1.2 : -0.1\nrestart = at.nitckpt\n");
    if (!resumed) FAIL(resumed.error().message);
    REQUIRE(resumed->sweep.has_value());
    REQUIRE(resumed->sweep->points.size() == 6);
    const auto worst = [](const std::vector<double>& u, const std::vector<double>& v) {
        double d = 0.0, s = 1e-300;
        for (std::size_t i = 0; i < u.size(); ++i) {
            d = std::max(d, std::abs(u[i] - v[i]));
            s = std::max(s, std::abs(v[i]));
        }
        return d / s;
    };
    for (std::size_t k = 0; k < 6; ++k) {
        const auto& a = resumed->sweep->points[k];
        const auto& b = full->sweep->points[k + 7];
        const double current = worst(a.terminal_current, b.terminal_current);  // reported
        const double psi = worst(a.fields.potential_V, b.fields.potential_V);
        double quantum = 0.0;  // absolute, eV: as psi in V
        for (std::size_t i = 0; i < plan.device.mesh().node_count(); ++i) {
            quantum = std::max({quantum, std::abs(a.fields.quantum_potential_n_eV[i] - b.fields.quantum_potential_n_eV[i]),
                               std::abs(a.fields.quantum_potential_p_eV[i] - b.fields.quantum_potential_p_eV[i])});
        }
        CAPTURE(k, current, psi, quantum);
        // A MOS-C at DC: the currents are at their resolution, so within it (not relative).
        for (std::size_t c = 0; c < a.terminal_current.size(); ++c) {
            CHECK(std::abs(a.terminal_current[c] - b.terminal_current[c]) <=
                  std::max(a.terminal_current_resolution[c], b.terminal_current_resolution[c]));
        }
        CHECK(psi <= 1e-9);
        CHECK(quantum <= 1e-9);
    }
    // The model changed: refused both ways.
    const auto off = run(mos_device + "[sweep]\ncontact gate = 0 : -1.2 : -0.1\nrestart = at.nitckpt\n");
    REQUIRE_FALSE(off.has_value());
    CHECK(off.error().message.find("was solved with quantum_confinement density_gradient and this run does not "
                                   "solve it") != std::string::npos);
    const auto classical = run(mos_device +
                               "[sweep]\ncontact gate = 0 : -0.6 : -0.1\n[output]\ncheckpoint = plain.nitckpt\n");
    if (!classical) FAIL(classical.error().message);
    const auto on = run(device + "[sweep]\ncontact gate = 0 : -1.2 : -0.1\nrestart = plain.nitckpt\n");
    REQUIRE_FALSE(on.has_value());
    CHECK(on.error().message.find("was solved without quantum_confinement density_gradient and this run solves it") !=
          std::string::npos);
}

namespace {

// Unit 62 QV13 / QV15: a sweep's report with the model at q. Converged or stopped (the stop and the run's
// warnings are printed, not failed: QV15 gates convergence below).
struct ExampleRun {
    results::Sweep sweep;
    double seconds = 0.0;
    std::size_t iterations = 0;
};

ExampleRun run_example(deck::Plan plan, assemble::QuantumConfinement q) {
    plan.options.models.quantum_confinement = q;
    plan.outputs = {};
    const auto t0 = std::chrono::steady_clock::now();
    auto report = deck::run(plan);
    ExampleRun r;
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    r.sweep = std::move(*report->sweep);
    for (const auto& p : r.sweep.points) r.iterations += p.convergence.iterations.size();
    return r;
}

// V_th by the largest-g_m linear extrapolation (minus V_D / 2), the smallest subthreshold swing [mV/dec], over the
// gate column `g` and the drain current column `d`.
std::pair<double, double> vth_swing(const results::Sweep& s, std::size_t g, std::size_t d) {
    const auto& P = s.points;
    double best = 0.0, vth = std::nan(""), swing = std::numeric_limits<double>::infinity();
    for (std::size_t k = 1; k + 1 < P.size(); ++k) {
        const double gm = (P[k + 1].terminal_current[d] - P[k - 1].terminal_current[d]) /
                          (P[k + 1].bias_V[g] - P[k - 1].bias_V[g]);
        if (gm > best) {
            best = gm;
            vth = P[k].bias_V[g] - P[k].terminal_current[d] / gm - 0.5 * P[k].bias_V[d];
        }
    }
    for (std::size_t k = 1; k < P.size(); ++k) {
        const double a = P[k - 1].terminal_current[d], b = P[k].terminal_current[d];
        if (a <= 0.0 || b <= a) continue;
        swing = std::min(swing, 1e3 * (P[k].bias_V[g] - P[k - 1].bias_V[g]) / std::log10(b / a));
    }
    return {vth, swing};
}

// The largest |Lambda| / kT over a sweep's points, and where.
void report_lambda(const char* name, const deck::Plan& p, const results::Sweep& s) {
    const double kT = 0.025851999786435;
    double largest = 0.0;
    std::size_t node = 0, point = 0;
    for (std::size_t k = 0; k < s.points.size(); ++k) {
        const auto& f = s.points[k].fields;
        for (std::size_t i = 0; i < f.quantum_potential_n_eV.size(); ++i) {
            const double l = std::max(std::abs(f.quantum_potential_n_eV[i]), std::abs(f.quantum_potential_p_eV[i])) / kT;
            if (l > largest) {
                largest = l;
                node = i;
                point = k;
            }
        }
    }
    const auto& x = p.device.mesh().points()[node];
    std::printf("QV15 %s: largest |Lambda| %.3f kT at node %zu (", name, largest, node);
    for (int a = 0; a < p.device.mesh().dimension(); ++a) std::printf("%s%.4f", a ? ", " : "", x[a] * 1e4);
    std::printf(" um), point %zu\n", point);
    for (const auto& w : s.run.warnings) std::printf("QV15 %s: warning: %s\n", name, w.c_str());
    if (s.stopped) std::printf("QV15 %s: stopped after %zu points: %s\n", name, s.points.size(), s.stopped->message.c_str());
}

}  // namespace

// Unit 62 QV13 (reported) and QV15 (gated: every point converges; reported: Lambda, the warnings, the I_d change),
// the examples on their own meshes, and the 2D MOSFET refined to 0.1 nm at the surface. Minutes; run by name.
TEST_CASE("density gradient: the examples on their meshes, and refined (QV13, QV15)", "[.dg_examples]") {
    const auto off = assemble::QuantumConfinement::off, dg = assemble::QuantumConfinement::density_gradient;
    // The 1D diode: no insulator, so Lambda acts only at the junction and contacts.
    {
        const deck::Plan p = example(fs::path(NITCAD_EXAMPLES_DIR) / "pn_diode_1d.deck");
        const auto c = run_example(p, off);
        const auto q = run_example(p, dg);
        report_lambda("1D diode", p, q.sweep);
        CHECK_FALSE(q.sweep.stopped.has_value());
        if (!q.sweep.stopped) {
            std::printf("QV15 1D diode: I(0.7 V) change %.3e; %.2f s, %zu iterations (classical %.2f s, %zu)\n",
                        q.sweep.points.back().terminal_current[0] / c.sweep.points.back().terminal_current[0] - 1.0,
                        q.seconds, q.iterations, c.seconds, c.iterations);
        }
    }
    // The MOSFETs: contacts source, drain, gate, body.
    const auto mosfet = [&](const char* name, const deck::Plan& p, bool gated) {
        const auto c = run_example(p, off);
        const auto m = run_example(p, mlda);
        const auto q = run_example(p, dg);
        report_lambda(name, p, q.sweep);
        if (gated) CHECK_FALSE(q.sweep.stopped.has_value());  // QV15 gates the examples' own meshes
        REQUIRE_FALSE(c.sweep.stopped.has_value());
        REQUIRE_FALSE(m.sweep.stopped.has_value());
        for (const auto* r : {&c, &m, &q}) {
            const char* model = r == &c ? "classical" : r == &m ? "MLDA" : "DG";
            if (r->sweep.stopped) continue;
            const auto [vth, swing] = vth_swing(r->sweep, 2, 1);
            std::printf("QV13 %s %s: V_th %.4f V, swing %.1f mV/dec, I_on %.6e; %.2f s, %zu iterations "
                        "(%.4f s per iteration)\n",
                        name, model, vth, swing, r->sweep.points.back().terminal_current[1], r->seconds, r->iterations,
                        r->seconds / static_cast<double>(r->iterations));
        }
        if (!q.sweep.stopped) {
            std::printf("QV15 %s: I_on change against classical %.4f\n", name,
                        q.sweep.points.back().terminal_current[1] / c.sweep.points.back().terminal_current[1] - 1.0);
        }
    };
    mosfet("2D MOSFET", example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck"), true);
    {
        // Refined: 10 nm to 0.94 um, 1 nm to 0.948 um, 0.1 nm to the interface at 0.95 um.
        std::ifstream in(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck", std::ios::binary);
        std::string text{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
        const std::string coarse = "y = linspace(0, 0.9um, 19) + linspace(0.9um, 0.95um, 6) + linspace(0.95um, 1um, 3)";
        const auto at = text.find(coarse);
        REQUIRE(at != std::string::npos);
        text.replace(at, coarse.size(),
                     "y = linspace(0, 0.9um, 19) + linspace(0.9um, 0.94um, 5) + linspace(0.94um, 0.948um, 9) + "
                     "linspace(0.948um, 0.95um, 21) + linspace(0.95um, 1um, 3)");
        mosfet("2D MOSFET refined", must(text), false);  // reported
        deck::Plan cut = must(text);  // with bias-step cutting (a deck option), every model
        cut.options.step_control.max_cuts = 8;
        mosfet("2D MOSFET refined, max_cuts 8", cut, false);
    }
    mosfet("3D MOSFET", example(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_3d.deck"), true);
}
