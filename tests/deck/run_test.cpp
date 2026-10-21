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
#include "NiTCAD/render/image.hpp"
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
profile = acceptor constant 1e17 axis=y x=-1:1um
profile = donor constant 1e17 axis=y x=1um:3um

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

TEST_CASE("nitcad plot: images from decks and result files; a foreign mesh refused (Unit 25a, VG5)", "[deck]") {
    const fs::path dir = scratch("cli_plot");
    const auto write = [&](const char* name, const std::string& text) {
        std::ofstream(dir / name, std::ios::binary) << text;
        return "\"" + (dir / name).string() + "\"";
    };
    const auto out = [&](const char* name) { return "\"" + (dir / name).string() + "\""; };
    const auto image_ok = [&](const char* name, int w, int h) {
        const auto img = render::load_png(dir / name);
        REQUIRE(img.has_value());
        CHECK(img->width == w);
        CHECK(img->height == h);
    };
    // A 2D pn diode: the deck's sweep and device files.
    const std::string pn2d = R"(
[mesh]
type = tensor
x = linspace(0, 2um, 21)
y = linspace(0, 1um, 6)
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e17 axis=y x=-1:1um
profile = donor constant 1e17 axis=y x=1um:3um
[contact anode]
kind = ohmic
patch = x_min
[contact cathode]
kind = ohmic
patch = x_max
[sweep]
contact anode = 0 : 0.4 : 0.2
[output]
sweep = pn.nitsweep
device = pn.nitdev
)";
    const std::string deck = write("pn.deck", pn2d);
    // Views needing results before the run: refused with the reason; the regions need none.
    const Cli early = cli("plot " + deck + " --field n --adapter warp --out " + out("early.png"));
    CHECK(early.status == 1);
    CHECK(early.output.find("(run the deck first)") != std::string::npos);
    const Cli solid = cli("plot " + deck + " --adapter warp --size 320x240 --out " + out("solid.png"));
    CHECK(solid.status == 0);
    image_ok("solid.png", 320, 240);
    REQUIRE(cli("run " + deck).status == 0);

    for (const auto& [args, name] :
         {std::pair{std::string("--field n --scale log --contours"), "n.png"},
          std::pair{std::string("--view mesh"), "mesh.png"},
          std::pair{std::string("--field doping"), "doping.png"},
          std::pair{std::string("--view iv --log"), "iv.png"}}) {
        const Cli c = cli("plot " + deck + " " + args + " --adapter warp --size 400x300 --out " + out(name));
        CAPTURE(args, c.output);
        CHECK(c.status == 0);
        CHECK(c.output.find("wrote ") != std::string::npos);
        image_ok(name, 400, 300);
    }
    // From the result file with its device file, at a chosen point.
    const Cli from_file = cli("plot " + out("pn.nitsweep") + " --device " + out("pn.nitdev") +
                              " --field potential --point 1 --adapter warp --size 400x300 --out " + out("psi1.png"));
    CHECK(from_file.status == 0);
    image_ok("psi1.png", 400, 300);
    const Cli beyond = cli("plot " + out("pn.nitsweep") + " --device " + out("pn.nitdev") +
                           " --point 9 --field n --out " + out("x.png"));
    CHECK(beyond.status == 1);
    CHECK(beyond.output.find("--point 9 is out of range (3 points)") != std::string::npos);

    // Another device's mesh: refused (by the result's mesh digest since Unit 45).
    std::string other = pn2d;
    other.replace(other.find("linspace(0, 1um, 6)"), 19, "linspace(0, 1um, 4)");
    other.replace(other.find("pn.nitsweep"), 11, "o.nitsweep");
    other.replace(other.find("pn.nitdev"), 9, "o.nitdev");
    REQUIRE(cli("run " + write("o.deck", other)).status == 0);
    const Cli foreign = cli("plot " + out("pn.nitsweep") + " --device " + out("o.nitdev") + " --field n --out " +
                            out("foreign.png"));
    CHECK(foreign.status == 1);
    CHECK(foreign.output.find("its mesh digest is not the device's") != std::string::npos);  // Unit 45
    CHECK_FALSE(fs::exists(dir / "foreign.png"));


    // A generated mesh: from the result file with its device and --mesh (the simplices), and
    // refused with another deck's mesh file.
    const auto generated = [&](const char* tag, const char* size) {
        return std::string(R"(
[mesh]
type = generate2d
patch left = 0 0, 0 1um
patch right = 2um 0, 2um 1um
max_size = )") + size + R"(
interface_resolution = 20nm
[region p]
material = silicon
polygon = 0 0, 1um 0, 1um 1um, 0 1um
[region n]
material = silicon
polygon = 1um 0, 2um 0, 2um 1um, 1um 1um
[doping]
profile = acceptor constant 1e16 x=0:1um axis=y
profile = donor constant 1e16 x=1um:2um axis=y
[contact left]
kind = ohmic
patch = left
[contact right]
kind = ohmic
patch = right
[sweep]
contact left = 0, 0.1
[output]
mesh = )" + tag + R"(.nitmesh
device = )" + tag + R"(.nitdev
sweep = )" + tag + R"(.nitsweep
)";
    };
    REQUIRE(cli("run " + write("g.deck", generated("g", "0.25um"))).status == 0);
    REQUIRE(cli("run " + write("h.deck", generated("h", "0.2um"))).status == 0);
    const Cli gen = cli("plot " + out("g.nitsweep") + " --device " + out("g.nitdev") + " --mesh " + out("g.nitmesh") +
                        " --field potential --adapter warp --size 400x300 --out " + out("g.png"));
    CAPTURE(gen.output);
    CHECK(gen.status == 0);
    image_ok("g.png", 400, 300);
    const Cli wrong = cli("plot " + out("g.nitsweep") + " --device " + out("g.nitdev") + " --mesh " +
                          out("h.nitmesh") + " --field potential --out " + out("wrong.png"));
    CHECK(wrong.status == 1);
    CHECK(wrong.output.find("is not the device's mesh") != std::string::npos);

    // Unit 45: a mesh with the same node count is refused by its digest; the rates draw; a file
    // written before Unit 45 falls back to the node count and holds no rates.
    std::string same_count = pn2d;
    same_count.replace(same_count.find("linspace(0, 2um, 21)"), 20, "linspace(0, 3um, 21)");
    same_count.replace(same_count.find("pn.nitsweep"), 11, "s.nitsweep");
    same_count.replace(same_count.find("pn.nitdev"), 9, "s.nitdev");
    REQUIRE(cli("run " + write("s.deck", same_count)).status == 0);
    const Cli digest = cli("plot " + out("pn.nitsweep") + " --device " + out("s.nitdev") + " --field n --out " +
                           out("digest.png"));
    CHECK(digest.status == 1);
    CHECK(digest.output.find("its mesh digest is not the device's") != std::string::npos);
    const Cli rate = cli("plot " + deck + " --field R --adapter warp --size 400x300 --out " + out("r.png"));
    CAPTURE(rate.output);
    CHECK(rate.status == 0);
    image_ok("r.png", 400, 300);
    const fs::path old = fs::path(NITCAD_IO_FIXTURES) / "pre_unit42";
    const std::string old_sweep = "\"" + (old / "bar.nitsweep").string() + "\"";
    const std::string old_device = "\"" + (old / "bar.nitdev").string() + "\"";
    const Cli fallback = cli("plot " + old_sweep + " --device " + old_device +
                             " --field potential --adapter warp --size 400x300 --out " + out("old.png"));
    CAPTURE(fallback.output);
    CHECK(fallback.status == 0);
    const Cli no_rates = cli("plot " + old_sweep + " --device " + old_device + " --field R --out " + out("old_r.png"));
    CHECK(no_rates.status == 1);
    CHECK(no_rates.output.find("holds no R rates") != std::string::npos);
    // An old file with another device: the node counts refuse it, named in the message.
    const Cli counts = cli("plot " + old_sweep + " --device " + out("pn.nitdev") + " --field n --out " + out("c.png"));
    CHECK(counts.status == 1);
    CHECK(counts.output.find("they do not belong together") != std::string::npos);
    // Usage: exit 64.
    CHECK(cli("plot").status == 64);
    CHECK(cli("plot " + deck).status == 64);  // no --out
    CHECK(cli("plot " + deck + " --view nothing --out " + out("u.png")).status == 64);
    CHECK(cli("plot " + deck + " --field Jz --out " + out("u.png")).status == 64);
    CHECK(cli("plot " + out("pn.nitsweep") + " --out " + out("u.png")).status == 64);  // no --device
}

TEST_CASE("nitcad plot: |E|, |J|, |Jn| and |Jp| (Unit 25b, BG6)", "[deck]") {
    const fs::path dir = scratch("cli_plot_vectors");
    const auto write = [&](const char* name, const std::string& text) {
        std::ofstream(dir / name, std::ios::binary) << text;
        return "\"" + (dir / name).string() + "\"";
    };
    const auto out = [&](const char* name) { return "\"" + (dir / name).string() + "\""; };
    const auto pn = [](const char* mesh, const char* tag, const char* extra, int D = 2) {
        return std::string("[mesh]\ntype = tensor\n") + mesh + R"(
[region silicon]
material = silicon
[doping]
)" + std::string(D == 1 ? "profile = acceptor constant 1e17\nprofile = donor erfc 2e17 center=1um width=0.02um\n"
                                   : "profile = acceptor constant 1e17 axis=y x=-1:1um\n"
                                     "profile = donor constant 1e17 axis=y x=1um:3um\n") + R"([contact anode]
kind = ohmic
patch = x_min
[contact cathode]
kind = ohmic
patch = x_max
[output]
sweep = )" + tag + ".nitsweep\ndevice = " + tag + ".nitdev\n" + extra;
    };
    // 2D: each magnitude drawn from the deck and from the result file.
    const std::string deck2 = write("pn2.deck", pn("x = linspace(0, 2um, 21)\ny = linspace(0, 1um, 6)\n", "pn2",
                                                   "[sweep]\ncontact anode = 0 : 0.4 : 0.2\n"));
    REQUIRE(cli("run " + deck2).status == 0);
    for (const char* f : {"E", "J", "Jn", "Jp"}) {
        const std::string name = std::string(f) + ".png";
        const Cli c = cli("plot " + deck2 + " --field " + f + " --adapter warp --size 400x300 --out " + out(name.c_str()));
        CAPTURE(f, c.output);
        CHECK(c.status == 0);
        CHECK(render::load_png(dir / name).has_value());
    }
    const Cli from_file = cli("plot " + out("pn2.nitsweep") + " --device " + out("pn2.nitdev") +
                              " --field J --point 0 --adapter warp --size 400x300 --out " + out("j0.png"));
    CAPTURE(from_file.output);
    CHECK(from_file.status == 0);
    // 1D: a profile of |E|.
    const std::string deck1 = write("pn1.deck", pn("x = linspace(0, 2um, 41)\n", "pn1",
                                                   "[sweep]\ncontact anode = 0, 0.3\n", 1));
    REQUIRE(cli("run " + deck1).status == 0);
    const Cli profile = cli("plot " + deck1 + " --field E --adapter warp --size 400x300 --out " + out("e1.png"));
    CAPTURE(profile.output);
    CHECK(profile.status == 0);
    CHECK(render::load_png(dir / "e1.png").has_value());
    // A quasi-static result's edge currents are zero: |J| drawn as zero (a linear 0..1 scale), |E| as
    // usual.
    const std::string qs = write("qs.deck", pn("x = linspace(0, 2um, 21)\ny = linspace(0, 1um, 6)\n", "qs",
                                               "[sweep]\ncontact anode = 0\n[solver]\nequations = equilibrium_poisson\n"));
    REQUIRE(cli("run " + qs).status == 0);
    const Cli zero_j = cli("plot " + qs + " --field J --adapter warp --size 400x300 --out " + out("qs_j.png"));
    CAPTURE(zero_j.output);
    CHECK(zero_j.status == 0);
    CHECK(cli("plot " + qs + " --field E --adapter warp --size 400x300 --out " + out("qs_e.png")).status == 0);
}

TEST_CASE("nitcad plot: isosurfaces and contours at levels (Unit 25b, BG6)", "[deck]") {
    const fs::path dir = scratch("cli_plot_iso");
    const auto write = [&](const char* name, const std::string& text) {
        std::ofstream(dir / name, std::ios::binary) << text;
        return "\"" + (dir / name).string() + "\"";
    };
    const auto out = [&](const char* name) { return "\"" + (dir / name).string() + "\""; };
    const auto pn = [](const char* mesh, const char* tag) {
        return std::string("[mesh]\ntype = tensor\n") + mesh + R"(
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e17 axis=y x=-1:1um
profile = donor constant 1e17 axis=y x=1um:3um
[contact anode]
kind = ohmic
patch = x_min
[contact cathode]
kind = ohmic
patch = x_max
[sweep]
contact anode = 0, 0.3
[output]
sweep = )" + tag + ".nitsweep\ndevice = " + tag + ".nitdev\n";
    };
    const std::string d3 = write("pn3.deck", pn("x = linspace(0, 2um, 11)\ny = linspace(0, 1um, 5)\nz = linspace(0, 0.5um, 3)\n", "pn3"));
    REQUIRE(cli("run " + d3).status == 0);
    for (const auto& [args, name] : {std::pair{std::string("--field n --levels 1e10,1e15"), "iso_n.png"},
                                     std::pair{std::string("--field potential --levels 0.1"), "iso_psi.png"}}) {
        const Cli c = cli("plot " + d3 + " --view iso " + args + " --adapter warp --size 400x300 --out " + out(name));
        CAPTURE(args, c.output);
        CHECK(c.status == 0);
        CHECK(render::load_png(dir / name).has_value());
    }
    const std::string d2 = write("pn2.deck", pn("x = linspace(0, 2um, 21)\ny = linspace(0, 1um, 6)\n", "pn2"));
    REQUIRE(cli("run " + d2).status == 0);
    const Cli contours = cli("plot " + d2 + " --view iso --field n --levels 1e12 --adapter warp --size 400x300 --out " +
                             out("iso2.png"));
    CAPTURE(contours.output);
    CHECK(contours.status == 0);
    // Refused: no levels (usage), a level off the log scale, a 1D device.
    CHECK(cli("plot " + d3 + " --view iso --field n --out " + out("x.png")).status == 64);
    CHECK(cli("plot " + d3 + " --view iso --field n --levels 1e10,x --out " + out("x.png")).status == 64);
    const Cli negative = cli("plot " + d3 + " --view iso --field n --levels -5 --out " + out("x.png"));
    CHECK(negative.status == 1);
    CHECK(negative.output.find("not finite on the field's scale") != std::string::npos);
    const std::string d1 = write("pn1.deck", "[mesh]\ntype = tensor\nx = linspace(0, 1um, 11)\n[region silicon]\n"
                                             "material = silicon\n[doping]\nprofile = donor constant 1e16\n"
                                             "[contact a]\nkind = ohmic\npatch = x_min\n[contact b]\nkind = ohmic\n"
                                             "patch = x_max\n[sweep]\ncontact a = 0\n[output]\nsweep = pn1.nitsweep\n");
    REQUIRE(cli("run " + d1).status == 0);
    const Cli one = cli("plot " + d1 + " --view iso --field n --levels 1e16 --out " + out("x.png"));
    CHECK(one.status == 1);
    CHECK(one.output.find("isosurfaces need a 2D or 3D device") != std::string::npos);
}
