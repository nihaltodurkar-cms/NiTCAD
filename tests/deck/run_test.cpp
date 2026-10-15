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
