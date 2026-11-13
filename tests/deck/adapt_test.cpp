// [adapt] in the deck, the runner and nitcad (Unit 17, slice 6; DECISIONS.md AD7, AD8, AD13, AR8;
// AC9, AC10): the section's options and refusals with their lines, the deck rebuilt on a refined
// mesh, the outputs of an adaptive run (the final mesh's, the report, the passes' meshes, the
// certificates in the table) and their round trip, the run identity, determinism and the exit status.
#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/io/adapt_report.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/render/image.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;
using deck_test::plan_of;
using deck_test::scratch;
namespace fs = std::filesystem;

namespace {

// An abrupt 1e17 diode on a coarse 2D tensor mesh (21 x 2 nodes), forward to 0.4 V. The doping is a
// window across y (axis=y) so x=... bounds it; a profile ignores its own depth axis's window.
const std::string coarse_diode = R"(
[mesh]
type = tensor
x = linspace(0, 2um, 21)
y = 0, 0.1um

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

std::string bytes_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

// The error of a deck that must not build; its message must hold `what`.
std::string refused(const std::string& text) {
    auto p = plan_of(text);
    if (p) FAIL("the deck built");
    return p.error().message;
}

struct Cli {
    int status;
    std::string output;
};

Cli cli(const std::string& arguments) {
    const std::string command = "\"\"" NITCAD_CLI_PATH "\" " + arguments + " 2>&1\"";
    FILE* pipe = _popen(command.c_str(), "r");
    REQUIRE(pipe != nullptr);
    std::string out;
    char buffer[512];
    while (std::fgets(buffer, sizeof buffer, pipe)) out += buffer;
    return {_pclose(pipe), out};
}

void write_text(const fs::path& p, const std::string& text) {
    std::ofstream f(p, std::ios::binary);
    f << text;
}

// A report with the measured seconds, what is fitted to them and the provenance cleared: what two
// runs must agree on bit for bit (AD9; AR7's exclusion).
results::AdaptReport deterministic(results::AdaptReport a) {
    a.run.provenance = {};
    a.coarse_seconds = a.verification_seconds = 0.0;
    for (results::AdaptPassRow& p : a.passes) {
        p.seconds_mesh = p.seconds_transfer = p.seconds_solve = p.seconds_estimate = 0.0;
        p.seconds_estimate_factorize = p.seconds_transposed = 0.0;
        p.predicted_seconds = p.predicted_seconds_nodes = 0.0;
        p.cost_a = p.cost_b = p.cost_c = p.cost_d = 0.0;
        p.cost_mesh_per_node = p.cost_transfer_per_node = p.cost_residual = 0.0;
        p.cost_alpha = p.cost_beta = p.cost_residual_nodes = 0.0;
        for (results::AdaptCandidateRow& c : p.candidates) c.cost = 0.0;
    }
    return a;
}

}  // namespace

TEST_CASE("deck adapt: [adapt] gives the loop's options; the deck rebuilds on a refined mesh (AD7)", "[deck]") {
    const deck::Plan plain = must(coarse_diode);
    CHECK_FALSE(plain.adapt.has_value());

    const deck::Plan p = must(coarse_diode + R"(
[adapt]
goal = current anode
tolerance = 1e-2
at = 0, 2
indicators = potential, doping
weighting = adjoint
verify = richardson
certify = required
max_passes = 4
max_nodes = 500
min_size = 1nm
grading = 1.4
theta = 0.4
)");
    REQUIRE(p.adapt.has_value());
    const solve::AdaptOptions& o = p.adapt->options;
    REQUIRE(o.goals.size() == 1);
    CHECK(o.goals[0].contact == 0);
    CHECK(o.tolerance == 1e-2);
    CHECK(o.at == std::vector<std::size_t>{0, 2});
    CHECK(o.indicators.enabled == std::array<bool, solve::indicator_count>{true, false, false, true, false});
    CHECK(o.weighting == solve::Weighting::adjoint);
    CHECK(o.policy == solve::Policy::basic);  // ID9: not the default until AC14-AC17 pass
    CHECK(o.verify);
    CHECK(p.adapt->certify_required);
    CHECK(o.max_passes == 4);
    CHECK(o.max_nodes == 500);
    CHECK(o.min_size == 1e-7);
    CHECK(o.grading == 1.4);
    CHECK(o.theta == 0.4);
    CHECK(p.adapt->kind == solve::AdaptMeshKind::tensor);
    REQUIRE(p.adapt->deck != nullptr);

    // Defaults; at = auto with the intelligent policy chooses up to 3 points (ID3, AR8).
    const deck::Plan d = must(coarse_diode + "[adapt]\ngoal = current anode\npolicy = intelligent\n");
    CHECK(d.adapt->options.tolerance == 1e-3);
    CHECK(d.adapt->options.at.empty());
    CHECK(d.adapt->options.intelligent.max_points == 3);
    CHECK(d.adapt->options.weighting == solve::Weighting::none);
    CHECK_FALSE(d.adapt->certify_required);
    CHECK(d.adapt->options.max_passes == 6);

    // The deck rebuilt: an empty refinement is the deck's mesh (AD10), a refined axis replaces it.
    auto again = deck::build(*p.adapt->deck, solve::MeshRefinement{});
    REQUIRE(again.has_value());
    CHECK(mesh::graph_digest(again->device.mesh()) == mesh::graph_digest(p.device.mesh()));
    CHECK(mesh::graph_digest(p.device.mesh()) == mesh::graph_digest(plain.device.mesh()));
    solve::MeshRefinement r;
    r.axes[0] = deck_test::linspace(0.0, 2e-4, 41);
    auto finer = deck::build(*p.adapt->deck, r);
    REQUIRE(finer.has_value());
    CHECK(finer->device.mesh().node_count() == 82);
    const solve::AdaptableMesh m = deck::adaptable(std::move(*finer));
    CHECK(m.kind == solve::AdaptMeshKind::tensor);
    auto built = deck::device_builder(p)(r);
    REQUIRE(built.has_value());
    CHECK(mesh::graph_digest(built->device.mesh()) == mesh::graph_digest(m.device.mesh()));
}

TEST_CASE("deck adapt: a tensor straddle is protected and kept when the axis is refined", "[deck]") {
    const deck::Plan p = must(R"(
[mesh]
type = tensor
x = linspace(0, 2um, 21)
straddle x = 1um 2nm
[region silicon]
material = silicon
box = 0 1um
[region silicon2]
material = silicon
[contact anode]
kind = ohmic
patch = x_min
[contact cathode]
kind = ohmic
patch = x_max
[sweep]
contact anode = 0
[adapt]
goal = current anode
)");
    REQUIRE(p.adapt.has_value());
    REQUIRE(p.adapt->axis_protected[0].size() == 1);
    CHECK(std::abs(p.adapt->axis_protected[0][0].lo - (1e-4 - 1e-7)) < 1e-18);
    CHECK(std::abs(p.adapt->axis_protected[0][0].hi - (1e-4 + 1e-7)) < 1e-18);
    // The refined axis holds the straddle's lines already: taken whole, not straddled again.
    const auto x = p.device.cells()->axis(0);
    solve::MeshRefinement r;
    r.axes[0].assign(x.begin(), x.end());
    auto again = deck::build(*p.adapt->deck, r);
    REQUIRE(again.has_value());
    CHECK(mesh::graph_digest(again->device.mesh()) == mesh::graph_digest(p.device.mesh()));
}

TEST_CASE("deck adapt: refusals name their line (AC9, AC10)", "[deck]") {
    const std::string base = coarse_diode + "[adapt]\n";
    const int line = static_cast<int>(std::count(base.begin(), base.end(), '\n')) + 1;  // the first [adapt] key
    const auto at = [&](int k) { return "t.deck:" + std::to_string(line + k) + ":"; };
    const auto holds = [](const std::string& message, const std::string& what) {
        INFO(message);
        CHECK(message.find(what) != std::string::npos);
    };
    {
        const std::string m = refused(base + "goal = vt\n");
        holds(m, at(0));
        holds(m, "not adapted in this unit");
    }
    holds(refused(base + "goal = current nowhere\n"), "no contact named 'nowhere'");
    holds(refused(base + "goal = charge anode\n"), "unknown goal");
    holds(refused(base + "goal = current anode\nfoo = 1\n"), at(1) + " unknown key 'foo' in [adapt]");
    holds(refused(base + "goal = current anode\nindicators = curvature\n"), "unknown indicator 'curvature'");
    holds(refused(base + "goal = current anode\nat = 7\n"), "not a point of the sweep");
    holds(refused(base + "goal = current anode\ntolerance = 2\n"), "in (0, 1)");
    holds(refused(base + "goal = current anode\nmode = budget\nbudget_dofs = 100\n"), "needs policy = intelligent");
    holds(refused(base + "goal = current anode\npolicy = intelligent\nmode = budget\n"), "budget_dofs or budget_seconds");
    holds(refused(base + "goal = current anode\nbudget_dofs = 100\n"), "is for mode = budget");
    holds(refused(base + "goal = current anode\npredictor = m.nitlearn\n"), "Unit 17b");
    holds(refused(base + "goal = current anode\nmax_points = 2\n"), "max_points is for at = auto");
    holds(refused(base + "tolerance = 1e-2\n"), "at least one 'goal = current <contact>'");
    holds(refused(base + "goal = current anode\ngoal = current anode\n"), "a second goal on 'anode'");
    // Goal weighting with a backend that cannot transpose and no fallback (AC10).
    {
        const std::string m = refused(coarse_diode + "[solver]\nbackend = gmres_ilut\nlinear_fallback = none\n[adapt]\n"
                                                     "goal = current anode\nweighting = adjoint\n");
        holds(m, "weighting = adjoint:");
    }
    // Another analysis, a restart, a general 3D mesh (AQ2), the adaptation's outputs without it.
    {
        std::string eq = coarse_diode;
        eq.replace(eq.find("[sweep]\ncontact anode = 0 : 0.4 : 0.2"), 37, "[equilibrium]");
        holds(refused(eq + "[adapt]\ngoal = current anode\n"), "not a sweep");
    }
    holds(refused(coarse_diode + "restart = c.nitckpt\n[adapt]\ngoal = current anode\n"), "with a restart");
    holds(refused(R"(
[mesh]
type = generate3d
cuboid = 0 1um 0 1um 0 1um
max_size = 0.5um
[region silicon]
material = silicon
[contact anode]
kind = ohmic
patch = x_min
[contact cathode]
kind = ohmic
patch = x_max
[sweep]
contact anode = 0
[adapt]
goal = current anode
)"),
          "general 3D or file mesh is not adapted in this unit");
    holds(refused(coarse_diode + "[output]\nadapt = r.nitadapt\n"), "the deck has no [adapt]");
}

TEST_CASE("deck adapt: identity stamp covers the settings and the final mesh (AD8)", "[deck]") {
    results::RunRecord r0;
    r0.input_identity = 0x1234;
    r0.mesh_digest = 77;
    solve::AdaptOptions a;
    a.goals = {{0}};
    a.tolerance = 1e-2;
    results::RunRecord r1 = r0, r2 = r0, r3 = r0, r4 = r0;
    solve::stamp_adapt(r1, a, false);
    solve::stamp_adapt(r2, a, false);
    CHECK(r1.input_identity == r2.input_identity);
    CHECK(r1.input_identity != r0.input_identity);
    solve::AdaptOptions b = a;
    b.tolerance = 2e-2;
    solve::stamp_adapt(r3, b, false);
    CHECK(r3.input_identity != r1.input_identity);
    solve::stamp_adapt(r4, a, true);
    CHECK(r4.input_identity != r1.input_identity);
    results::RunRecord r5 = r0;
    r5.mesh_digest = 78;
    solve::stamp_adapt(r5, a, false);
    CHECK(r5.input_identity != r1.input_identity);
    const auto has = [&](const char* name, double v) {
        return std::ranges::any_of(r1.settings, [&](const auto& s) { return s.first == name && s.second == v; });
    };
    CHECK(has("adapt.tolerance", 1e-2));
    CHECK(has("adapt.goal", 0.0));
    CHECK(has("adapt.certify_required", 0.0));
}

TEST_CASE("deck adapt: the run writes the final mesh's outputs, the report and the passes' meshes (AC9)", "[deck]") {
    const fs::path dir = scratch("adapt_run");
    const std::string text = coarse_diode + R"(
[adapt]
goal = current anode
tolerance = 1e-2
at = 2

[output]
sweep = out.nitsweep
device = out.nitdev
mesh = out.nitmesh
checkpoint = out.nitckpt
table = iv.csv
adapt = out.nitadapt
adapt_meshes = passes
)";
    const deck::Plan plan = must(text, dir);
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->adapt.has_value());
    REQUIRE(report->adapted.has_value());
    REQUIRE(report->sweep.has_value());
    const results::AdaptReport& a = *report->adapt;
    INFO("passes " << a.passes.size() << ", stop " << a.stop << ", status " << a.status);
    CHECK(a.stop == "converged");
    CHECK(a.status == "certified");
    REQUIRE(a.passes.size() >= 2);
    REQUIRE(a.certificates.size() == 1);
    CHECK(a.certificates[0].point == 2);
    CHECK(a.goal_contacts == std::vector<std::string>{"anode"});

    // The final mesh is the result's, everywhere.
    const std::uint64_t fin = mesh::graph_digest(report->adapted->device.mesh());
    CHECK(a.final_mesh_digest == fin);
    CHECK(report->sweep->run.mesh_digest == fin);
    CHECK(fin != mesh::graph_digest(plan.device.mesh()));
    const auto m = io::load_mesh(dir / "out.nitmesh");
    REQUIRE(m.has_value());
    CHECK(mesh::graph_digest(*m) == fin);
    const auto dev = io::load_device(dir / "out.nitdev");
    REQUIRE(dev.has_value());
    CHECK(mesh::graph_digest(dev->mesh()) == fin);
    const auto c = io::load_checkpoint(dir / "out.nitckpt", report->adapted->device.mesh());
    REQUIRE(c.has_value());
    CHECK(c->run_identity == report->sweep->run.input_identity);
    const auto s = io::load_sweep(dir / "out.nitsweep");
    REQUIRE(s.has_value());
    CHECK(s->run.input_identity == report->sweep->run.input_identity);
    CHECK(std::ranges::any_of(s->run.settings, [](const auto& x) { return x.first == "adapt.tolerance" && x.second == 1e-2; }));

    // The passes' meshes, one device and mesh per pass, as the report names them.
    for (std::size_t k = 0; k < a.passes.size(); ++k) {
        char name[32];
        std::snprintf(name, sizeof name, "pass_%02zu", k);
        const auto pd = io::load_device(dir / "passes" / (std::string(name) + ".nitdev"));
        REQUIRE(pd.has_value());
        CHECK(mesh::graph_digest(pd->mesh()) == a.passes[k].mesh_digest);
        CHECK(fs::exists(dir / "passes" / (std::string(name) + ".nitmesh")));
    }

    // The report round-trips exactly: saved again from the loaded copy, the same bytes.
    const auto loaded = io::load_adapt_report(dir / "out.nitadapt");
    REQUIRE(loaded.has_value());
    CHECK(loaded->passes.size() == a.passes.size());
    CHECK(loaded->status == a.status);
    CHECK(loaded->run.input_identity == report->sweep->run.input_identity);
    REQUIRE(io::save_adapt_report(dir / "again.nitadapt", *loaded).has_value());
    CHECK(bytes_of(dir / "again.nitadapt") == bytes_of(dir / "out.nitadapt"));
    // Another content kind is refused.
    CHECK_FALSE(io::load_adapt_report(dir / "out.nitsweep").has_value());

    // The certificates in the table (AD13), at the adaptation point only.
    const std::string table = bytes_of(dir / "iv.csv");
    CHECK(table.find(",substeps,anode certificate,anode estimate A,anode verified error A,anode allowed A,anode order\n") !=
          std::string::npos);
    CHECK(table.find("\n0,") != std::string::npos);
    const auto row2 = table.find("\n2,");
    REQUIRE(row2 != std::string::npos);
    CHECK(table.find(",certified,", row2) != std::string::npos);
    const auto row0 = table.find("\n0,");
    CHECK(table.substr(row0, table.find('\n', row0 + 1) - row0).ends_with(",,,,,"));
}

TEST_CASE("deck adapt: two adaptive runs are bit-identical (AC10)", "[deck]") {
    const fs::path dir = scratch("adapt_twice");
    const deck::Plan plan = must(coarse_diode + "[adapt]\ngoal = current anode\ntolerance = 1e-2\n", dir);
    const auto r1 = deck::run(plan);
    const auto r2 = deck::run(plan);
    REQUIRE(r1.has_value());
    REQUIRE(r2.has_value());
    REQUIRE(io::save_adapt_report(dir / "1.nitadapt", deterministic(*r1->adapt)).has_value());
    REQUIRE(io::save_adapt_report(dir / "2.nitadapt", deterministic(*r2->adapt)).has_value());
    CHECK(bytes_of(dir / "1.nitadapt") == bytes_of(dir / "2.nitadapt"));
    CHECK(r1->sweep->run.input_identity == r2->sweep->run.input_identity);
    CHECK(mesh::graph_digest(r1->adapted->device.mesh()) == mesh::graph_digest(r2->adapted->device.mesh()));
}

TEST_CASE("nitcad: an uncertified result exits 3 when certification is required (AD13)", "[deck]") {
    const fs::path dir = scratch("adapt_cli");
    // One pass: a budget stop, unverified (AD13).
    const std::string adapt = coarse_diode + "[adapt]\ngoal = current anode\ntolerance = 1e-2\nmax_passes = 1\n";
    write_text(dir / "required.deck", adapt + "certify = required\n");
    write_text(dir / "report.deck", adapt + "certify = report\n");
    const Cli check = cli("check \"" + (dir / "required.deck").string() + "\"");
    CHECK(check.status == 0);
    CHECK(check.output.find("adapt: current of anode to 0.01 (relative), policy basic, weighting none, verify richardson, "
                            "certify required") != std::string::npos);
    const Cli required = cli("run \"" + (dir / "required.deck").string() + "\"");
    INFO(required.output);
    CHECK(required.status == 3);
    CHECK(required.output.find("certificate: current of anode at point 2") != std::string::npos);
    CHECK(required.output.find("adapt: unverified") != std::string::npos);
    CHECK(required.output.find("requires it certified") != std::string::npos);
    const Cli report = cli("run \"" + (dir / "report.deck").string() + "\"");
    CHECK(report.status == 0);
    // The deck's errors name their line.
    write_text(dir / "bad.deck", coarse_diode + "[adapt]\ngoal = gm\n");
    const Cli bad = cli("check \"" + (dir / "bad.deck").string() + "\"");
    CHECK(bad.status == 1);
    CHECK(bad.output.find("bad.deck:") != std::string::npos);
    CHECK(bad.output.find("not adapted in this unit") != std::string::npos);
}

TEST_CASE("nitcad plot: an adaptive run's estimate, its passes' meshes and its final mesh (AD8)", "[deck]") {
    const fs::path dir = scratch("adapt_plot");
    write_text(dir / "a.deck", coarse_diode + R"(
[adapt]
goal = current anode
tolerance = 1e-2
at = 2

[output]
sweep = a.nitsweep
device = a.nitdev
mesh = a.nitmesh
adapt = a.nitadapt
adapt_meshes = passes
)");
    const std::string deck = "\"" + (dir / "a.deck").string() + "\"";
    const auto out = [&](const char* name) { return "\"" + (dir / name).string() + "\""; };
    const auto image = [&](const char* name) {
        auto img = render::load_png(dir / name);
        REQUIRE(img.has_value());
        return std::move(*img);
    };
    // Before the run: no report.
    const Cli early = cli("plot " + deck + " --view adapt --adapter warp --out " + out("early.png"));
    CHECK(early.status == 1);
    const Cli run = cli("run " + deck);
    INFO(run.output);
    REQUIRE(run.status == 0);
    const auto report = io::load_adapt_report(dir / "a.nitadapt");
    REQUIRE(report.has_value());
    const std::size_t passes = report->passes.size();
    REQUIRE(passes >= 2);

    const Cli estimate = cli("plot " + deck + " --view adapt --adapter warp --size 400x300 --out " + out("estimate.png"));
    INFO(estimate.output);
    CHECK(estimate.status == 0);
    CHECK(estimate.output.find(std::to_string(passes) + " passes") != std::string::npos);
    const render::Image e = image("estimate.png");
    CHECK(e.width == 400);
    CHECK(e.height == 300);
    // The report file alone gives the same picture.
    const Cli from_report = cli("plot " + out("a.nitadapt") + " --view adapt --adapter warp --size 400x300 --out " +
                                out("estimate2.png"));
    CHECK(from_report.status == 0);
    CHECK(image("estimate2.png").rgba == e.rgba);

    // Each pass's mesh, tiled in rows of --columns.
    const Cli tiles = cli("plot " + deck + " --view passes --columns 2 --adapter warp --size 200x150 --out " +
                          out("passes.png"));
    INFO(tiles.output);
    CHECK(tiles.status == 0);
    const render::Image t = image("passes.png");
    CHECK(t.width == 400);
    CHECK(t.height == static_cast<int>(150 * ((passes + 1) / 2)));

    // The other views draw the final mesh with the results computed on it.
    const Cli field = cli("plot " + deck + " --field n --adapter warp --size 320x240 --out " + out("n.png"));
    INFO(field.output);
    CHECK(field.status == 0);

    // Refusals: the passes need the deck; a report plots only its estimate; a pass file changed.
    CHECK(cli("plot " + out("a.nitadapt") + " --view passes --adapter warp --out " + out("x.png")).status == 64);
    CHECK(cli("plot " + out("a.nitadapt") + " --view mesh --adapter warp --out " + out("x.png")).status == 64);
    fs::copy_file(deck::pass_device_file(dir / "passes", passes - 1), deck::pass_device_file(dir / "passes", 0),
                  fs::copy_options::overwrite_existing);
    const Cli changed = cli("plot " + deck + " --view passes --adapter warp --out " + out("x.png"));
    CHECK(changed.status == 1);
    CHECK(changed.output.find("digest differs") != std::string::npos);
}
