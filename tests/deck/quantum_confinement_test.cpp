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
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/quantum_confinement.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
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
    const auto bad = deck::parse(with("quantum_confinement = density_gradient"), "bad.deck");
    REQUIRE(bad.has_value());
    const auto built = deck::build(*bad);
    REQUIRE_FALSE(built.has_value());
    CHECK(built.error().message.find("quantum_confinement is off or mlda") != std::string::npos);
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
