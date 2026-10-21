// The batch runner (Unit 31, deck/run.hpp and the nitcad executable): a deck's run equals the
// hand-built solve bit for bit, its outputs read back, a restart continues a stopped sweep, and the
// executable's exit status says how the run ended.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_device;
using deck_test::hand_diode;
using deck_test::must;
using deck_test::scratch;
namespace fs = std::filesystem;

namespace {

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
    }
    return true;
}

double relative(const std::vector<double>& u, const std::vector<double>& v) {
    REQUIRE(u.size() == v.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        worst = std::max(worst, std::abs(u[i] - v[i]) / std::max({std::abs(u[i]), std::abs(v[i]), 1e-300}));
    }
    return worst;
}

std::string text_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

struct Cli {
    int status;
    std::string output;
};

// Runs the nitcad executable with these arguments, stdout and stderr together.
Cli cli(const std::string& arguments) {
    const std::string command = "\"\"" NITCAD_CLI_PATH "\" " + arguments + " 2>&1\"";
    FILE* pipe = _popen(command.c_str(), "r");
    REQUIRE(pipe != nullptr);
    std::string out;
    char buffer[512];
    while (std::fgets(buffer, sizeof buffer, pipe)) out += buffer;
    return {_pclose(pipe), out};
}

}  // namespace

TEST_CASE("deck run: the diode sweep equals the hand-built sweep; outputs read back", "[deck]") {
    const fs::path dir = scratch("run_diode");
    const deck::Plan plan = must(diode_device + R"(
[sweep]
contact anode = 0 : 0.3 : 0.1

[output]
sweep = out.nitsweep
device = out.nitdev
mesh = out.nitmesh
checkpoint = out.nitckpt
table = iv.csv
)", dir);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    const results::Sweep& sweep = *report->sweep;
    REQUIRE_FALSE(sweep.stopped.has_value());
    REQUIRE(sweep.points.size() == 4);

    const auto hand = solve::sweep_bias(hand_diode(), plan.points);
    REQUIRE(hand.has_value());
    CHECK(sweep.run.input_identity == hand->run.input_identity);
    for (std::size_t k = 0; k < 4; ++k) {
        CHECK(same_bits(sweep.points[k].terminal_current, hand->points[k].terminal_current));
        CHECK(same_bits(sweep.points[k].fields.potential_V, hand->points[k].fields.potential_V));
    }

    CHECK(report->written == std::vector<fs::path>{dir / "out.nitdev", dir / "out.nitmesh",
                                                   dir / "out.nitsweep", dir / "iv.csv",
                                                   dir / "out.nitckpt"});
    const auto loaded = io::load_sweep(dir / "out.nitsweep");
    REQUIRE(loaded.has_value());
    CHECK(loaded->run.input_identity == sweep.run.input_identity);
    CHECK(same_bits(loaded->points[3].terminal_current, sweep.points[3].terminal_current));
    const auto device = io::load_device(dir / "out.nitdev");
    REQUIRE(device.has_value());
    CHECK(solve::make_run_record(*device, plan.options, plan.points).input_identity == sweep.run.input_identity);
    const auto m = io::load_mesh(dir / "out.nitmesh");
    REQUIRE(m.has_value());
    CHECK(io::mesh_digest(*m) == io::mesh_digest(plan.device.mesh()));
    const auto c = io::load_checkpoint(dir / "out.nitckpt", plan.device.mesh());
    REQUIRE(c.has_value());
    CHECK(c->bias_V == std::vector<double>{0.3, 0.0});
    CHECK(c->run_identity == sweep.run.input_identity);

    const std::string table = text_of(dir / "iv.csv");
    CHECK(table.starts_with("point,anode V,cathode V,anode A,cathode A,substeps\n0,0,0,"));
    CHECK(std::count(table.begin(), table.end(), '\n') == 5);
    CHECK(table.find("\n3,0.29999999999999999,0,") != std::string::npos);  // 17 digits
    CHECK_FALSE(fs::exists(dir / "iv.csv.tmp"));
}

TEST_CASE("deck run: the weighted current in the sweep and the table (Unit 32)", "[deck]") {
    const fs::path dir = scratch("run_weighted");
    const deck::Plan plan = must(diode_device + R"(
[solver]
weighted_current = on
[sweep]
contact anode = 0 : 0.2 : 0.1
[output]
table = iv.csv
)", dir);
    CHECK(plan.options.weighted_current);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    const auto& last = report->sweep->points.back();
    REQUIRE(last.weighted_current.size() == 2);
    const std::string table = text_of(dir / "iv.csv");
    CHECK(table.starts_with("point,anode V,cathode V,anode A,cathode A,anode weighted A,cathode "
                            "weighted A,anode weighted resolution A,cathode weighted resolution "
                            "A,substeps\n"));
    char row[64];
    std::snprintf(row, sizeof row, ",%.17g,", last.weighted_current[0]);
    CHECK(table.find(row) != std::string::npos);
}

TEST_CASE("deck run: a restart continues after the checkpoint's point", "[deck]") {
    const fs::path dir = scratch("run_restart");
    const deck::Plan first = must(diode_device + R"(
[sweep]
contact anode = 0 : 0.3 : 0.1
[output]
checkpoint = at.nitckpt
)", dir);
    REQUIRE(deck::run(first).has_value());

    const std::string rest = diode_device + R"(
[sweep]
contact anode = 0 : 0.6 : 0.1
restart = at.nitckpt
)";
    const deck::Plan second = must(rest, dir);
    const auto resumed = deck::run(second);
    if (!resumed) FAIL(resumed.error().message);
    CHECK(resumed->skipped_points == 4);
    REQUIRE(resumed->sweep.has_value());
    REQUIRE(resumed->sweep->points.size() == 3);
    CHECK(resumed->sweep->points[0].bias_V == std::vector<double>{0.4, 0.0});

    const auto full = solve::sweep_bias(second.device, second.points);
    REQUIRE(full.has_value());
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(resumed->sweep->points[k].bias_V == full->points[k + 4].bias_V);
        CHECK(relative(resumed->sweep->points[k].terminal_current, full->points[k + 4].terminal_current) < 1e-9);
        CHECK(relative(resumed->sweep->points[k].fields.n_cm3, full->points[k + 4].fields.n_cm3) < 1e-9);
    }

    // The checkpoint's bias is not a point of this sweep; or it is the last point.
    const auto off = deck::run(must(diode_device + "[sweep]\ncontact anode = 0 : 0.6 : 0.2\nrestart = at.nitckpt\n", dir));
    REQUIRE_FALSE(off.has_value());
    CHECK(off.error().message.find("not a point of the sweep") != std::string::npos);
    const auto done = deck::run(must(diode_device + "[sweep]\ncontact anode = 0 : 0.3 : 0.1\nrestart = at.nitckpt\n", dir));
    REQUIRE_FALSE(done.has_value());
    CHECK(done.error().message.find("nothing is left") != std::string::npos);
}

TEST_CASE("deck run: the equilibrium analysis", "[deck]") {
    const fs::path dir = scratch("run_equilibrium");
    const deck::Plan plan = must(diode_device + "[equilibrium]\n[output]\ncheckpoint = eq.nitckpt\n", dir);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->equilibrium.has_value());
    CHECK_FALSE(report->sweep.has_value());
    const auto hand = solve::solve_equilibrium(hand_diode());
    REQUIRE(hand.has_value());
    CHECK(same_bits(report->equilibrium->fields.potential_V, hand->fields.potential_V));
    const auto c = io::load_checkpoint(dir / "eq.nitckpt", plan.device.mesh());
    REQUIRE(c.has_value());
    CHECK(c->bias_V == std::vector<double>{0.0, 0.0});
    CHECK(same_bits(c->fields.potential_V, hand->fields.potential_V));
}

TEST_CASE("deck run: a stopped sweep writes its completed points", "[deck]") {
    const fs::path dir = scratch("run_stopped");
    // 1 V in one step with a few Newton iterations: the first point converges, the jump does not.
    const deck::Plan plan = must(diode_device + R"(
[solver]
max_iterations = 8
[sweep]
contact anode = 0, 1.0
[output]
sweep = s.nitsweep
checkpoint = s.nitckpt
)", dir);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    REQUIRE(report->sweep->stopped.has_value());
    CHECK(report->sweep->points.size() == 1);
    const auto loaded = io::load_sweep(dir / "s.nitsweep");
    REQUIRE(loaded.has_value());
    REQUIRE(loaded->stopped.has_value());
    CHECK(loaded->stopped->code == report->sweep->stopped->code);
    const auto c = io::load_checkpoint(dir / "s.nitckpt", plan.device.mesh());
    REQUIRE(c.has_value());
    CHECK(c->bias_V == std::vector<double>{0.0, 0.0});
}

TEST_CASE("nitcad: check, run and the exit status", "[deck]") {
    const fs::path dir = scratch("cli");
    const auto write = [&](const char* name, const std::string& text) {
        std::ofstream(dir / name, std::ios::binary) << text;
        return "\"" + (dir / name).string() + "\"";
    };
    const std::string ok = write("ok.deck", deck_test::diode_deck + "[output]\nsweep = ok.nitsweep\n");
    const Cli check = cli("check " + ok);
    CHECK(check.status == 0);
    CHECK(check.output.find("mesh: 1D, 201 nodes") != std::string::npos);
    CHECK(check.output.find("contact anode: ohmic, 1 nodes") != std::string::npos);
    CHECK(check.output.find("analysis: sweep of 3 point(s)") != std::string::npos);
    CHECK_FALSE(fs::exists(dir / "ok.nitsweep"));

    const Cli run = cli("run " + ok);
    CHECK(run.status == 0);
    CHECK(run.output.find("point 3 of 3: converged") != std::string::npos);
    CHECK(run.output.find("wrote ") != std::string::npos);
    CHECK(io::load_sweep(dir / "ok.nitsweep").has_value());

    const Cli bad = cli("check " + write("bad.deck", "[mesh]\ntype = tensor\nx = 0, 1um\n[region]\n"));
    CHECK(bad.status == 1);
    CHECK(bad.output.find("nitcad: bad.deck:4: [region] needs a name") != std::string::npos);

    const Cli stopped = cli("run " + write("stop.deck", diode_device + "[solver]\nmax_iterations = 8\n"
                                                                       "[sweep]\ncontact anode = 0, 1.0\n"));
    CHECK(stopped.status == 2);
    CHECK(stopped.output.find("the sweep stopped after 1 of 2 point(s)") != std::string::npos);

    CHECK(cli("").status == 64);
    CHECK(cli("simulate " + ok).status == 64);
}

TEST_CASE("nitcad: the transient, small-signal and trace analyses (Unit 43)", "[deck]") {
    const fs::path dir = scratch("cli_analyses");
    const auto write = [&](const char* name, const std::string& text) {
        std::ofstream(dir / name, std::ios::binary) << text;
        return "\"" + (dir / name).string() + "\"";
    };
    const std::string transient = write("tr.deck", diode_device + R"(
[transient]
contact anode = ramp(0, 0.3, 0, 2e-10)
t_end = 4e-10
dt_initial = 1e-10
adaptive = off
[output]
transient = tr.nittran
)");
    const Cli tc = cli("check " + transient);
    CHECK(tc.status == 0);
    CHECK(tc.output.find("analysis: transient to 4e-10 s") != std::string::npos);
    CHECK(tc.output.find("  contact anode: piecewise linear") != std::string::npos);
    CHECK(tc.output.find("  contact cathode: constant") != std::string::npos);
    const Cli tr = cli("run " + transient);
    CHECK(tr.status == 0);
    CHECK(tr.output.find("t = 4e-10 s: converged") != std::string::npos);
    CHECK(tr.output.find("build: nitcad ") != std::string::npos);
    CHECK(io::load_transient(dir / "tr.nittran").has_value());

    const std::string ac = write("ac.deck", diode_device + "[small_signal]\ncontact anode = 0\n"
                                                          "frequencies = logspace(1e3, 1e6, 2)\n");
    const Cli ac_check = cli("check " + ac);
    CHECK(ac_check.output.find("analysis: small signal at 1 point(s), 2 frequencies from 1000 to 1e+06 Hz") !=
          std::string::npos);
    const Cli ac_run = cli("run " + ac);
    CHECK(ac_run.status == 0);
    CHECK(ac_run.output.find("point 1 of 1, 1e+06 Hz: solved") != std::string::npos);

    const std::string trace = write("trace.deck", diode_device + "[trace]\ncontact = anode\nend = 0.2\n"
                                                                 "start cathode = 0\nmax_step = 0.1\n");
    CHECK(cli("check " + trace).output.find("analysis: trace of contact anode from 0 to 0.2 V") != std::string::npos);
    const Cli trace_run = cli("run " + trace);
    CHECK(trace_run.status == 0);
    CHECK(trace_run.output.find("traced point 2: converged") != std::string::npos);

    // Stopped early: exit 2, the reason on stderr. A 1 V jump with few Newton iterations and no
    // room to cut the step (transient), or a 1 V operating point (small signal).
    const Cli tstop = cli("run " + write("tstop.deck", diode_device + R"(
[solver]
max_iterations = 8
[transient]
contact anode = step(0, 1.0, 1e-10)
t_end = 1e-9
dt_initial = 1e-10
dt_min = 1e-10
adaptive = off
)"));
    CHECK(tstop.status == 2);
    CHECK(tstop.output.find("nitcad: the transient stopped at t = ") != std::string::npos);
    const Cli astop = cli("run " + write("astop.deck", diode_device + "[solver]\nmax_iterations = 8\n"
                                                                      "[small_signal]\ncontact anode = 0, 1.0\n"
                                                                      "frequencies = 1e6\n"));
    CHECK(astop.status == 2);
    CHECK(astop.output.find("nitcad: the small-signal run stopped after 1 of 2 point(s)") != std::string::npos);
}

TEST_CASE("deck run: a stopped transient writes its accepted steps (Unit 43)", "[deck]") {
    const fs::path dir = scratch("run_transient_stopped");
    const deck::Plan plan = must(diode_device + R"(
[solver]
max_iterations = 8
[transient]
contact anode = step(0, 1.0, 2e-10)
t_end = 1e-9
dt_initial = 1e-10
dt_min = 1e-10
adaptive = off
[output]
transient = t.nittran
table = t.csv
)", dir);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->transient.has_value());
    REQUIRE(report->transient->stopped.has_value());
    // The start, then 1e-10 and 2e-10 (the step that lands on the breakpoint sees the value before
    // the jump); the jump to 1 V fails.
    CHECK(report->transient->points.size() == 3);
    const auto loaded = io::load_transient(dir / "t.nittran");
    REQUIRE(loaded.has_value());
    REQUIRE(loaded->stopped.has_value());
    CHECK(loaded->stopped->code == report->transient->stopped->code);
    CHECK(loaded->points.size() == 3);
    CHECK(fs::exists(dir / "t.csv"));
}

TEST_CASE("deck run: a 3D diode deck with backend = gmres_ilut (Unit 44, LG5)", "[deck]") {
    // The deck keys reach the solver: the 3D sweep runs on GMRES-ILUT in strict mode (a fallback
    // cannot finish it) and its currents agree with the Eigen SparseLU run of the same deck.
    const std::string device = R"(
[mesh]
type = tensor
x = linspace(0, 2um, 21)
y = linspace(0, 1um, 6)
z = linspace(0, 1um, 6)

[region silicon]
material = silicon

[doping]
profile = acceptor constant 1e17 x=-1:1um
profile = donor constant 1e17 x=1um:3um

[contact anode]
kind = ohmic
patch = x_min

[contact cathode]
kind = ohmic
patch = x_max

[sweep]
contact anode = 0 : 0.4 : 0.2
)";
    const deck::Plan gmres = must(device + R"(
[solver]
backend = gmres_ilut
linear_fallback = none
)");
    REQUIRE(gmres.device.mesh().dimension() == 3);
    const auto g = deck::run(gmres);
    if (!g) FAIL(g.error().message);
    const auto e = deck::run(must(device));
    if (!e) FAIL(e.error().message);
    REQUIRE_FALSE(g->sweep->stopped.has_value());
    REQUIRE(g->sweep->points.size() == 3);
    double worst = 0.0;
    for (std::size_t k = 1; k < 3; ++k) {
        worst = std::max(worst, relative(g->sweep->points[k].terminal_current,
                                         e->sweep->points[k].terminal_current));
    }
    std::printf("gmres 3D deck: worst relative current difference %.2e\n", worst);
    CHECK(worst < 1e-8);
}
