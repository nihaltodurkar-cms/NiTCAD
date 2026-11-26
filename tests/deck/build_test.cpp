// The deck's meaning (Unit 31, deck/build.hpp): decks give the devices, options and points built by
// hand (equal run identities, so every input the solver reads is equal), and every mistake is an
// error naming its line.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/materials.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using base::ErrorCode;
namespace fs = std::filesystem;

using deck_test::diode_deck;
using deck_test::hand_diode;
using deck_test::must;
using deck_test::plan_of;
using deck_test::scratch;

namespace {

// The message of building `text`, which must fail.
std::string failure(const std::string& text, const fs::path& dir = {}) {
    const auto p = plan_of(text, dir);
    REQUIRE_FALSE(p.has_value());
    return p.error().message;
}

}  // namespace

TEST_CASE("deck build: the diode deck is the hand-built diode", "[deck]") {
    const deck::Plan p = must(diode_deck);
    const device::Device hand = hand_diode();
    CHECK(p.analysis == deck::Analysis::sweep);
    REQUIRE(p.points.size() == 3);
    CHECK(p.points[1] == std::vector<double>{0.1, 0.0});
    CHECK(p.points[2] == std::vector<double>{0.2, 0.0});
    CHECK(p.restart.empty());
    CHECK(p.outputs.sweep.empty());
    CHECK_FALSE(p.generated.has_value());
    // The run identity digests the mesh, the device, the options and the points.
    CHECK(solve::make_run_record(p.device, p.options, p.points).input_identity ==
          solve::make_run_record(hand, solve::BiasOptions{}, p.points).input_identity);
    CHECK(io::mesh_digest(p.device.mesh()) == io::mesh_digest(hand.mesh()));
    REQUIRE(p.device.cells() != nullptr);
}

TEST_CASE("deck build: models, solver options, series and outputs", "[deck]") {
    const fs::path dir = scratch("options");
    const deck::Plan p = must(diode_deck + R"(
[models]
field_mobility = on
auger = off
btbt_nonlocal = kane
impact_current_resolution = 1e-12
hole_mobility = masetti

[solver]
equations = drift_diffusion
max_iterations = 40
tol_update = 1e-9
max_update = 2
line_search = natural
line_search.max_backtracks = 12
natural_norm = weighted
max_cuts = 8
growth = 1.5
cut_max_iterations = 30
equilibrate = off
min_pivot_ratio = 0
Ns = 1e17
weighted_current = on

[output]
sweep = out/diode.nitsweep
table = "iv table.csv"
)", dir);
    const auto& m = p.options.models;
    CHECK(m.field_mobility);
    CHECK_FALSE(m.auger);
    CHECK(m.srh);
    CHECK(m.btbt_nonlocal == assemble::NonlocalTunnelling::kane);
    CHECK(m.impact_current_resolution == 1e-12);
    CHECK(m.hole_mobility == assemble::HoleMobility::masetti);
    const auto& o = p.options;
    CHECK(o.newton.max_iterations == 40);
    CHECK(o.newton.tol_update == 1e-9);
    CHECK(o.newton.max_update == 2.0);
    CHECK(o.newton.line_search.enabled);
    CHECK(o.newton.line_search.test == solve::LineSearchTest::natural);
    CHECK(o.newton.line_search.max_backtracks == 12);
    CHECK(o.newton.line_search.norm == solve::NaturalNorm::weighted);
    CHECK(o.step_control.max_cuts == 8);
    CHECK(o.step_control.growth == 1.5);
    CHECK(o.step_control.max_iterations == 30);
    CHECK_FALSE(o.linear.equilibrate);
    CHECK(o.linear.min_pivot_ratio == 0.0);
    REQUIRE(o.Ns_override.has_value());
    CHECK(*o.Ns_override == 1e17);
    CHECK(o.weighted_current);
    CHECK(p.outputs.sweep == dir / "out/diode.nitsweep");
    CHECK(p.outputs.table == dir / "iv table.csv");

    const deck::Plan q = must(R"(
[mesh]
type = tensor
x = 0, 1um, 2um
[region si]
material = silicon
[contact a]
kind = ohmic
patch = x_min
[contact b]
kind = ohmic
patch = x_max
[sweep]
contact b = 0.5
contact a = 0.1, 0.2, 0.4
)");
    // A single value is held; a contact not named would be 0 V.
    CHECK(q.points == std::vector<std::vector<double>>{{0.1, 0.5}, {0.2, 0.5}, {0.4, 0.5}});

    // Axis pieces joined by '+' keep a repeated join point once; '+' in an exponent is a number's.
    const deck::Plan joined = must(R"(
[mesh]
type = tensor
x = linspace(0, 1um, 3) + linspace(1um, 2um, 3) + 3e+0um
[region si]
material = silicon
[contact a]
kind = ohmic
patch = x_min
[equilibrium]
)");
    std::vector<double> x;
    for (const auto& point : joined.device.mesh().points()) x.push_back(point[0]);
    CHECK(x == std::vector<double>{0, 0.5e-4, 1e-4, 1e-4 + (2e-4 - 1e-4) * 0.5, 2e-4, 3.0 * 1e-4});
}

TEST_CASE("deck build: a 2D tensor device with an oxide, an electrode and interface traps", "[deck]") {
    const deck::Plan p = must(R"(
[mesh]
type = tensor
x = linspace(0, 1um, 11)
y = linspace(0, 0.3um, 7) + 0.35um, 0.4um
straddle y = 0.3um 10nm

[region silicon]
material = silicon
eps_r = 11.9
electron_mobility.mu_max = 1400

[region oxide]
material = sio2
box = 0 1um   0.3um 0.4um
thermal.conductivity_W_cmK = 0.015

[doping]
profile = acceptor constant 1e17
profile = donor gaussian 1e19 axis=y center=0 width=50nm x=0:0.3um lateral=20nm

[contact substrate]
kind = ohmic
patch = y_min
box = 0.2um 0.8um  0 0

[contact drain]
kind = ohmic
patch = x_max

[contact gate]
kind = electrode
patch = y_max
electrode = metal
work_function = 4.6

[contact source]
kind = ohmic
patch = x_min
box = 0 0  0.1um 0.2um

[interface silicon oxide]
fixed_charge = 5e10
trap_level = acceptor 1e10 0.2
trap_band = donor 1e11 -0.4 -0.1 1e-15 2e-15
recombination_velocity_n = 100

[equilibrium]
)");
    const device::Device& d = p.device;
    CHECK(p.analysis == deck::Analysis::equilibrium);
    REQUIRE(d.regions().size() == 2);
    const auto& si = std::get<physics::Semiconductor>(d.regions()[0].material).parameters();
    CHECK(si.eps_r == 11.9);
    CHECK(si.electron_mobility.mu_max == 1400.0);
    CHECK(si.hole_mobility == physics::silicon_parameters.hole_mobility);
    const auto& ox = std::get<physics::Insulator>(d.regions()[1].material).parameters();
    CHECK(ox.eps_r == 3.9);
    CHECK(ox.thermal.conductivity_W_cmK == 0.015);
    // Every node at y > 0.3 um is oxide; the straddle put nodes at 0.3 um -+ 5 nm.
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        const double y = d.mesh().points()[i][1];
        CHECK(d.node_region()[i] == (y > 0.3e-4 ? 1 : 0));
        CHECK(std::abs(y - 0.3e-4) > 4e-7);
    }
    const auto* gate = d.find_contact("gate");
    REQUIRE(gate != nullptr);
    CHECK(gate->kind == device::ContactKind::electrode);
    CHECK(gate->electrode.kind == device::GateElectrode::metal);
    CHECK(gate->electrode.work_function_eV == 4.6);
    CHECK(gate->nodes.size() == 11);
    const auto* source = d.find_contact("source");
    REQUIRE(source != nullptr);
    // The drain's patch reaches into the oxide: an ohmic contact takes its semiconductor nodes only.
    const auto* drain = d.find_contact("drain");
    REQUIRE(drain != nullptr);
    CHECK(drain->nodes.size() == 7);  // y = 0, 0.05, ..., 0.25, 0.295 um
    for (const auto k : drain->nodes) CHECK(d.node_region()[static_cast<std::size_t>(k)] == 0);
    // The x_min nodes with 0.1 um <= y <= 0.2 um (lengths read as 0.1 * 1e-4 and 0.2 * 1e-4).
    REQUIRE_FALSE(source->nodes.empty());
    for (const auto k : source->nodes) {
        const double y = d.mesh().points()[static_cast<std::size_t>(k)][1];
        CHECK(y >= 0.1 * 1e-4);
        CHECK(y <= 0.2 * 1e-4);
    }
    REQUIRE(d.interfaces().size() == 1);
    const auto& f = d.interfaces()[0];
    CHECK(f.fixed_charge_cm2 == 5e10);
    CHECK(f.recombination_velocity_n_cm_s == 100.0);
    REQUIRE(f.traps.levels.size() == 1);
    CHECK(f.traps.levels[0] == physics::TrapLevel{physics::TrapType::acceptor, 1e10, 0.2, 1e-15, 1e-15});
    REQUIRE(f.traps.bands.size() == 1);
    CHECK(f.traps.bands[0] == physics::TrapBand{physics::TrapType::donor, 1e11, -0.4, -0.1, 1e-15, 2e-15});
    // The Gaussian sits at the surface under x <= 0.3 um only.
    const std::size_t corner = 0;  // (0, 0)
    CHECK(d.donors()[corner] > 1e18);
    CHECK(d.acceptors()[corner] == 1e17);
}

// Unit 26k (GD1, GE2): a layered generate2d deck (a background, boxes and polygons in deck order, each replacing
// what it covers) builds to the same mesh bits as its hand-tiled equivalent; the errors.
TEST_CASE("deck build: a layered generate2d deck equals its hand tiling (Unit 26k)", "[deck]") {
    const std::string tail = R"(
[doping]
profile = acceptor constant 1e16

[contact source]
kind = ohmic
patch = source

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";
    const std::string layered = R"(
[mesh]
type = generate2d
layered = on
domain = 0 2um 0 1um
patch source = 0 0, 0 0.6um
patch body = 0.2um 0, 0.8um 0
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon

[region n]
material = silicon
polygon = 1um 0, 2um 0, 2um 0.5um, 2um 0.8um, 1um 0.8um

[region oxide]
material = sio2
box = 0 2um 0.8um 1um
)" + tail;
    // By hand: the junction under the oxide; the oxide carries the junction's end on its lower edge.
    const std::string tiled = R"(
[mesh]
type = generate2d
patch source = 0 0, 0 0.6um
patch body = 0.2um 0, 0.8um 0
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon
polygon = 0 0, 1um 0, 1um 0.8um, 0 0.8um

[region n]
material = silicon
polygon = 1um 0, 2um 0, 2um 0.8um, 1um 0.8um

[region oxide]
material = sio2
polygon = 0 0.8um, 1um 0.8um, 2um 0.8um, 2um 1um, 0 1um
)" + tail;
    const deck::Plan a = must(layered), b = must(tiled);
    CHECK(io::mesh_digest(a.device.mesh()) == io::mesh_digest(b.device.mesh()));
    CHECK(std::ranges::equal(a.device.node_region(), b.device.node_region()));

    const auto fails = [](const std::string& text, std::string_view message) {
        const auto p = plan_of(text);
        REQUIRE_FALSE(p.has_value());
        CAPTURE(p.error().message);
        CHECK(p.error().message.find(message) != std::string::npos);
    };
    std::string t = layered;
    fails(t.replace(t.find("layered = on\n"), 13, ""), "domain is for a layered generate2d mesh");
    t = tiled;
    fails(t.replace(t.find("polygon = 0 0.8um, 1um"), 22, "box = 0.5um 1.5um 0.9um 1um\n#"), "box with layered = on");
    t = layered;
    fails(t.replace(t.find("box = 0 2um 0.8um 1um"), 21, ""), "two regions without a shape");
    t = layered;
    fails(t.replace(t.find("polygon = 1um 0, 2um 0"), 22, "polygon = 1um 0, 1um 0, 2um 0"), "repeats a vertex");
    // A slanted well builds since slice 4 (GD10): its left edge from (1.1, 0) up to (1, 0.8).
    t = layered;
    const deck::Plan s = must(t.replace(t.find("polygon = 1um 0, 2um 0"), 22, "polygon = 1.1um 0, 2um 0"));
    CHECK(s.device.mesh().node_count() > 0);
}

TEST_CASE("deck build: a layered hole under a full-width contact meshes (Unit 26k slice 1)", "[deck]") {
    // A buried box in the background: the ring around it is cut by two rays, kept off the body patch along the
    // whole bottom (a straight cut down through the box would put a vertex inside it).
    const std::string deck = R"(
[mesh]
type = generate2d
layered = on
domain = 0 2um 0 1um
patch body = 0 0, 2um 0
patch top = 0.8um 1um, 1.2um 1um
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon

[region box]
material = sio2
box = 0.8um 1.2um 0.3um 0.6um

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[contact top]
kind = ohmic
patch = top

[equilibrium]
)";
    const auto p = plan_of(deck);
    INFO((p ? std::string() : p.error().message));
    REQUIRE(p);
    // Every node of the body contact on y = 0, from x = 0 to 2 um.
    CHECK(p->device.mesh().node_count() > 0);
}

TEST_CASE("deck build: generated 2D and 3D meshes and mesh files", "[deck]") {
    const fs::path dir = scratch("generated");
    const std::string two_d = R"(
[mesh]
type = generate2d
patch left = 0 0, 0 1um
patch right = 2um 0, 2um 1um
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon
polygon = 0 0, 1um 0, 1um 1um, 0 1um

[region n]
material = silicon
polygon = 1um 0, 2um 0, 2um 1um, 1um 1um

[doping]
profile = acceptor constant 1e16 x=0:1um
profile = donor constant 1e16 x=1um:2um

[contact left]
kind = ohmic
patch = left

[contact right]
kind = ohmic
patch = right

[output]
mesh = pn.nitmesh

[equilibrium]
)";
    const deck::Plan p = must(two_d, dir);
    REQUIRE(p.generated.has_value());
    REQUIRE(p.device.interface_certificate() != nullptr);
    mesh::Geometry2D g;
    g.polygons = {{{{0, 0, 0}, {1e-4, 0, 0}, {1e-4, 1e-4, 0}, {0, 1e-4, 0}}, 0},
                  {{{1e-4, 0, 0}, {2e-4, 0, 0}, {2e-4, 1e-4, 0}, {1e-4, 1e-4, 0}}, 1}};
    g.patches = {{"left", {0, 0, 0}, {0, 1e-4, 0}}, {"right", {2e-4, 0, 0}, {2e-4, 1e-4, 0}}};
    mesh::MeshOptions2D o;
    o.max_size = 0.25e-4;
    o.interface_resolution = 20 * 1e-7;
    const auto hand = mesh::generate_mesh(g, o);
    REQUIRE(hand.has_value());
    CHECK(io::mesh_digest(p.device.mesh()) == io::mesh_digest(hand->mesh));
    CHECK(std::vector<std::int32_t>(p.device.node_region().begin(), p.device.node_region().end()) ==
          hand->node_region);

    // The same mesh from a file: the regions and the certificate come with it.
    REQUIRE(io::save_generated_mesh(dir / "pn.nitmesh", *p.generated).has_value());
    std::string from_file = two_d;
    from_file.replace(from_file.find("type = generate2d"), 17, "type = file\npath = pn.nitmesh");
    for (const char* line : {"patch left = 0 0, 0 1um\n", "patch right = 2um 0, 2um 1um\n",
                             "max_size = 0.25um\n", "interface_resolution = 20nm\n",
                             "polygon = 0 0, 1um 0, 1um 1um, 0 1um\n",
                             "polygon = 1um 0, 2um 0, 2um 1um, 1um 1um\n"}) {
        from_file.erase(from_file.find(line), std::string(line).size());
    }
    const deck::Plan f = must(from_file, dir);
    CHECK(io::mesh_digest(f.device.mesh()) == io::mesh_digest(p.device.mesh()));
    REQUIRE(f.device.interface_certificate() != nullptr);
    const std::vector<std::vector<double>> zero{{0.0, 0.0}};
    CHECK(solve::make_run_record(f.device, f.options, zero).input_identity ==
          solve::make_run_record(p.device, p.options, zero).input_identity);
    // A plain mesh file is a one-region mesh.
    REQUIRE(io::save_mesh(dir / "plain.nitmesh", p.device.mesh()).has_value());
    CHECK(failure(from_file.replace(from_file.find("pn.nitmesh\n"), 10, "plain.nitmesh"), dir)
              .find("a plain mesh file has one region") != std::string::npos);

    const deck::Plan cube = must(R"(
[mesh]
type = generate3d
cuboid = 0 1um  0 1um  0 1um
max_size = 0.5um
[region si]
material = silicon
[doping]
profile = donor constant 1e16
[contact bottom]
kind = ohmic
patch = z_min
[contact top]
kind = ohmic
patch = z_max
[equilibrium]
)");
    const auto& m = cube.device.mesh();
    CHECK(m.dimension() == 3);
    const double volume = std::accumulate(m.volumes().begin(), m.volumes().end(), 0.0);
    CHECK(std::abs(volume - 1e-12) < 1e-24);
    for (const char* face : {"x_min", "x_max", "y_min", "y_max", "z_min", "z_max"}) {
        const auto* b = m.find_boundary(face);
        REQUIRE(b != nullptr);
        CHECK(std::abs(std::accumulate(b->areas.begin(), b->areas.end(), 0.0) - 1e-8) < 1e-20);
    }

    // Two regions split by a plane, seeds naming them, and the top face tiled into two patches.
    const deck::Plan split = must(R"(
[mesh]
type = generate3d
cuboid = 0 1um  0 1um  0 1um
plane = 0 0 0.5um, 1um 0 0.5um, 0 1um 0.5um
patch z_max left = 0 0 1um, 0.5um 0 1um, 0.5um 1um 1um, 0 1um 1um
patch z_max right = 0.5um 0 1um, 1um 0 1um, 1um 1um 1um, 0.5um 1um 1um
max_size = 0.5um
interface_resolution = 50nm
[region lower]
material = silicon
seed = 0.5um 0.5um 0.25um
[region upper]
material = silicon
seed = 0.5um 0.5um 0.75um
[doping]
profile = acceptor constant 1e16 z=0:0.5um
profile = donor constant 1e16 z=0.5um:1um
[contact bottom]
kind = ohmic
patch = z_min
[contact top]
kind = ohmic
patch = left
[equilibrium]
)");
    const auto& sm = split.device.mesh();
    REQUIRE(split.device.interface_certificate() != nullptr);
    CHECK(sm.find_boundary("z_max") == nullptr);
    for (const char* name : {"left", "right"}) {
        const auto* b = sm.find_boundary(name);
        REQUIRE(b != nullptr);
        CHECK(std::abs(std::accumulate(b->areas.begin(), b->areas.end(), 0.0) - 0.5e-8) < 1e-20);
    }
    for (std::size_t i = 0; i < sm.node_count(); ++i) {
        CHECK(split.device.node_region()[i] == (sm.points()[i][2] > 0.5e-4 ? 1 : 0));
    }
}

TEST_CASE("deck build: every mistake names its line", "[deck]") {
    const std::string head = "[mesh]\ntype = tensor\nx = 0, 1um, 2um\n[region si]\nmaterial = silicon\n"
                             "[contact a]\nkind = ohmic\npatch = x_min\n[contact b]\nkind = ohmic\npatch = x_max\n";
    // head is 11 lines; the next line is 12.
    const auto expect = [](const std::string& message, const std::string& what) {
        INFO(message);
        CHECK(message.find(what) != std::string::npos);
    };
    expect(failure(head + "[sweep]\ncontact a = 1\n[bogus]\n"), "t.deck:14: unknown section [bogus]");
    expect(failure(head + "[sweep]\ncontact c = 1\n"), "t.deck:13: no contact named 'c'");
    expect(failure(head + "[sweep]\ncontact a = 1\ncontact a = 2\n"), "t.deck:14: a second series for 'a'");
    expect(failure(head + "[sweep]\ncontact a = 0,1\ncontact b = 0,1,2\n"),
           "t.deck:14: this series has 3 values, the one on line 13 has 2");
    expect(failure(head + "[sweep]\nvoltage a = 1\n"), "t.deck:13: unknown key 'voltage a' in [sweep]");
    expect(failure(head), "the deck needs exactly one analysis: [equilibrium], [sweep], [transient], [small_signal] or [trace]");
    expect(failure(head + "[equilibrium]\n[sweep]\n"), "the deck needs exactly one analysis: [equilibrium], [sweep], [transient], [small_signal] or [trace]");
    expect(failure(head + "[equilibrium]\n[mesh]\n"), "t.deck:13: a second [mesh] section");
    expect(failure(head + "[equilibrium]\n[models]\nsrh = maybe\n"), "t.deck:14: 'maybe' is not on or off");
    expect(failure(head + "[equilibrium]\n[models]\nquantum = on\n"), "t.deck:14: unknown key 'quantum' in [models]");
    expect(failure(head + "[equilibrium]\n[models]\nhole_mobility = arora\n"),
           "t.deck:14: hole_mobility is automatic, caughey_thomas, masetti, klaassen or klaassen_equilibrium");
    expect(failure(head + "[equilibrium]\n[solver]\nbackend = pardiso\n"), "t.deck:13: backend = pardiso needs 'mkl'");
    expect(failure(head + "[equilibrium]\n[solver]\nmax_cuts = 99\n"), "t.deck:14: 'max_cuts' must lie in 0..52");
    expect(failure(head + "[equilibrium]\n[solver]\nnatural_norm = rms\n"),
           "t.deck:14: natural_norm is convergence or weighted");
    expect(failure(head + "[equilibrium]\n[doping]\nprofile = donor flat 1e16\n"),
           "t.deck:14: the shape is constant, gaussian or erfc");
    expect(failure(head + "[equilibrium]\n[doping]\nprofile = donor gaussian 1e16 axis=x\n"),
           "t.deck:14: [doping]: ");  // the profile library's error, located
    expect(failure(head + "[equilibrium]\n[output]\nsweep = a.nitsweep\n"), "an equilibrium run writes no sweep");
    expect(failure("[mesh]\ntype = tensor\nx = 0, 1um\n[region si]\nmaterial = unobtainium\n[equilibrium]\n"),
           "t.deck:5: unknown material 'unobtainium'");
    expect(failure("[mesh]\ntype = tensor\nx = 0, 1um\n[region ox]\nmaterial = sio2\nmu_max = 3\n[equilibrium]\n"),
           "t.deck:6: unknown key 'mu_max' in [region ox] (not a parameter of its material)");
    expect(failure("[mesh]\ntype = tensor\nx = 0, 1um\n[region si]\nmaterial = silicon\neps_r = -1\n[equilibrium]\n"),
           "t.deck:5: [region si]: ");  // Semiconductor::create's error, located at the material
    expect(failure("[mesh]\ntype = tensor\nx = 0, 1um\n[region si]\n[equilibrium]\n"),
           "t.deck:4: [region si] has no 'material'");
    expect(failure("[mesh]\ntype = tensor\nx = 0, 1um, 2um\n[region a]\nmaterial = silicon\nbox = 0 0.5um\n"
                   "[region b]\nmaterial = silicon\nbox = 1.5um 2um\n[equilibrium]\n"),
           "t.deck:1: the node at (0.000100");
    expect(failure("[mesh]\ntype = tensor\nx = 0, 1um, 2um\nstraddle x = 5um 1nm\n[region a]\nmaterial = silicon\n[equilibrium]\n"),
           "t.deck:4: straddle: ");
    expect(failure("[mesh]\ntype = tensor\nx = 0, 1um, 2um\nstraddle y = 1um 1nm\n[region a]\nmaterial = silicon\n[equilibrium]\n"),
           "t.deck:4: straddle on an axis the mesh does not have");
    expect(failure("[mesh]\ntype = tensor\nx = 0, linspace(1um)\n[region a]\nmaterial = silicon\n[equilibrium]\n"),
           "t.deck:3: ");
    expect(failure("[mesh]\ntype = hexagonal\n[region a]\nmaterial = silicon\n[equilibrium]\n"),
           "t.deck:2: unknown mesh type 'hexagonal'");
    expect(failure("[region a]\nmaterial = silicon\n[equilibrium]\n"), "the deck has no [mesh] section");
    expect(failure(head + "[contact c]\nkind = ohmic\npatch = y_min\n[equilibrium]\n"),
           "t.deck:14: the mesh has no boundary patch 'y_min'");
    expect(failure(head + "[contact c]\nkind = electrode\npatch = x_min\n[equilibrium]\n"),
           "t.deck:14: patch 'x_min' has no insulator node for this contact");
    expect(failure(head + "[contact c]\nkind = ohmic\npatch = x_min\noxide_thickness = 2nm\n[equilibrium]\n"),
           "t.deck:12: oxide_thickness, oxide_permittivity and fixed_charge are for a gate");
    // Device::create's own errors come located at line 1 with their message.
    expect(failure(head + "[contact a]\nkind = ohmic\npatch = x_max\n[equilibrium]\n"), "t.deck:1: device: ");
}

TEST_CASE("deck build: silicon, silicon_legacy and the hole mobility default (Unit 35)", "[deck]") {
    // silicon is the published band set, silicon_legacy the legacy's; hole_mobility defaults to
    // automatic.
    const deck::Plan p = must(diode_deck);
    const auto& si = std::get<physics::Semiconductor>(p.device.regions()[0].material).parameters();
    CHECK(si == physics::silicon_parameters);
    CHECK(p.options.models.hole_mobility == assemble::HoleMobility::automatic);
    std::string legacy = diode_deck;
    legacy.replace(legacy.find("material = silicon"), 18, "material = silicon_legacy");
    const deck::Plan q = must(legacy + "[models]\nhole_mobility = caughey_thomas\n");
    const auto& old = std::get<physics::Semiconductor>(q.device.regions()[0].material).parameters();
    CHECK(old == physics::legacy_silicon_parameters);
    CHECK(q.options.models.hole_mobility == assemble::HoleMobility::caughey_thomas);
    // A Passler parameter is overridable by its member path.
    std::string theta = diode_deck;
    theta.replace(theta.find("material = silicon"), 18, "material = silicon\npassler_gap.Theta_K = 400");
    const deck::Plan r = must(theta);
    CHECK(std::get<physics::Semiconductor>(r.device.regions()[0].material).parameters().passler_gap.Theta_K ==
          400.0);
}

TEST_CASE("deck build: ionization_levels and its run record (Unit 37, I5)", "[deck]") {
    const deck::Plan fixed = must(diode_deck);
    CHECK(fixed.options.models.ionization_levels == assemble::IonizationLevels::fixed);
    const deck::Plan a = must(diode_deck + "[models]\nincomplete_ionization = on\nionization_levels = altermatt\n");
    CHECK(a.options.models.ionization_levels == assemble::IonizationLevels::altermatt);
    const deck::Plan b = must(diode_deck + "[models]\nincomplete_ionization = on\n");
    const deck::Plan c = must(diode_deck + "[models]\nionization_levels = altermatt\n");
    const auto record = [](const deck::Plan& p) {
        return solve::make_run_record(p.device, p.options, p.points);
    };
    const auto named = [](const results::RunRecord& r) {
        return std::ranges::any_of(r.settings,
                                   [](const auto& s) { return s.first == "models.ionization_levels"; });
    };
    // Named, and the identity changed, only where it is read and not the default.
    CHECK(named(record(a)));
    CHECK_FALSE(named(record(b)));
    CHECK_FALSE(named(record(c)));
    CHECK(record(a).input_identity != record(b).input_identity);
    CHECK(record(c).input_identity == record(fixed).input_identity);
    CHECK(failure(diode_deck + "[models]\nionization_levels = schenk\n").find(
              "ionization_levels is fixed or altermatt") != std::string::npos);
}

TEST_CASE("deck build: electron_mobility (Unit 38)", "[deck]") {
    CHECK(must(diode_deck).options.models.electron_mobility ==
          assemble::ElectronMobility::caughey_thomas);
    CHECK(must(diode_deck + "[models]\nelectron_mobility = masetti\n").options.models.electron_mobility ==
          assemble::ElectronMobility::masetti);
    CHECK(failure(diode_deck + "[models]\nelectron_mobility = arora\n").find(
              "electron_mobility is caughey_thomas, masetti, klaassen or klaassen_equilibrium") != std::string::npos);
}

TEST_CASE("deck build: klaassen_equilibrium mobility (Unit 66, KC8, KC9)", "[deck]") {
    const auto built =
        must(diode_deck + "[models]\nelectron_mobility = klaassen_equilibrium\n"
                          "hole_mobility = klaassen_equilibrium\n");
    CHECK(built.options.models.electron_mobility == assemble::ElectronMobility::klaassen_equilibrium);
    CHECK(built.options.models.hole_mobility == assemble::HoleMobility::klaassen_equilibrium);
    // The run record carries the values KD7 fixes: electrons 3, holes 4.
    const auto r = solve::make_run_record(built.device, built.options, built.points);
    const auto has = [&](const char* key, double value) {
        for (const auto& [k, v] : r.settings) {
            if (k == key) return v == value;
        }
        return false;
    };
    CHECK(has("models.electron_mobility", 3.0));
    CHECK(has("models.hole_mobility", 4.0));
    // `klaassen` (carrier-dependent): electrons 2, holes 3.
    const auto k = must(diode_deck + "[models]\nelectron_mobility = klaassen\nhole_mobility = klaassen\n");
    CHECK(k.options.models.electron_mobility == assemble::ElectronMobility::klaassen);
    CHECK(k.options.models.hole_mobility == assemble::HoleMobility::klaassen);
    const auto rk = solve::make_run_record(k.device, k.options, k.points);
    const auto value = [&](const char* key) {
        for (const auto& [name, v] : rk.settings) {
            if (name == key) return v;
        }
        return -1.0;
    };
    CHECK(value("models.electron_mobility") == 2.0);
    CHECK(value("models.hole_mobility") == 3.0);
}

TEST_CASE("deck build: surface_mobility (Unit 46, LD1)", "[deck]") {
    CHECK(must(diode_deck).options.models.surface_mobility == assemble::SurfaceMobility::off);
    CHECK(must(diode_deck + "[models]\nsurface_mobility = lombardi\n").options.models.surface_mobility ==
          assemble::SurfaceMobility::lombardi);
    CHECK(must(diode_deck + "[models]\nsurface_mobility = off\n").options.models.surface_mobility ==
          assemble::SurfaceMobility::off);
    CHECK(failure(diode_deck + "[models]\nsurface_mobility = darwish\n").find(
              "surface_mobility is off or lombardi") != std::string::npos);
}

TEST_CASE("deck build: energy_balance, the relaxation times and the run record (Unit 52, EC16)", "[deck]") {
    CHECK_FALSE(must(diode_deck).options.models.energy_balance);
    const deck::Plan on = must(diode_deck + "[models]\nenergy_balance = on\n");
    CHECK(on.options.models.energy_balance);
    // Each region's relaxation times by their member paths (ED5: deck overridable per region).
    std::string tau = diode_deck;
    tau.replace(tau.find("material = silicon"), 18,
                "material = silicon\nenergy_relaxation.tau_n_s = 0.3e-12\nenergy_relaxation.tau_p_s = 0.45e-12");
    const deck::Plan given = must(tau + "[models]\nenergy_balance = on\n");
    const auto& e =
        std::get<physics::Semiconductor>(given.device.regions()[0].material).parameters().energy_relaxation;
    CHECK(e.tau_n_s == 0.3e-12);
    CHECK(e.tau_p_s == 0.45e-12);
    CHECK(failure(tau + "[models]\nenergy_balance = maybe\n").find("'maybe' is not on or off") != std::string::npos);
    std::string negative = diode_deck;
    negative.replace(negative.find("material = silicon"), 18,
                     "material = silicon\nenergy_relaxation.tau_n_s = -1e-12");
    CHECK(failure(negative).find("tau_n_s") != std::string::npos);
    // The run record names the model, and its identity reads the relaxation times only when on.
    const auto record = [](const deck::Plan& p) {
        return solve::make_run_record(p.device, p.options, p.points);
    };
    const auto named = [](const results::RunRecord& r) {
        return std::ranges::any_of(r.settings,
                                   [](const auto& s) { return s.first == "models.energy_balance"; });
    };
    const deck::Plan off = must(diode_deck);
    const deck::Plan given_off = must(tau);
    CHECK(named(record(on)));
    CHECK_FALSE(named(record(off)));
    CHECK(record(on).input_identity != record(off).input_identity);
    CHECK(record(given).input_identity != record(on).input_identity);
    CHECK(record(given_off).input_identity == record(off).input_identity);
    // Each carrier's time on its own (EC17: a digest of tau_p twice survived the checks above).
    for (const char* other : {"energy_relaxation.tau_n_s = 0.31e-12\nenergy_relaxation.tau_p_s = 0.45e-12",
                              "energy_relaxation.tau_n_s = 0.3e-12\nenergy_relaxation.tau_p_s = 0.46e-12"}) {
        std::string changed = diode_deck;
        changed.replace(changed.find("material = silicon"), 18, std::string("material = silicon\n") + other);
        CAPTURE(other);
        CHECK(record(must(changed + "[models]\nenergy_balance = on\n")).input_identity != record(given).input_identity);
    }
}

TEST_CASE("deck build: the GMRES-ILUT keys and their run record (Unit 44, LG5)", "[deck]") {
    const std::string solver = R"(
[solver]
backend = gmres_ilut
restart = 40
max_linear_iterations = 500
ilut_drop = 1e-5
ilut_fill = 20
linear_fallback = none
)";
    const deck::Plan p = must(diode_deck + solver);
    const linalg::SolverConfig& l = p.options.linear;
    CHECK(l.backend == linalg::SolverBackend::gmres_ilut);
    CHECK(l.restart == 40);
    CHECK(l.max_iterations == 500);
    CHECK(l.ilut_drop_tolerance == 1e-5);
    CHECK(l.ilut_fill_factor == 20);
    CHECK(l.fallback == linalg::LinearFallback::none);
    const deck::Plan d = must(diode_deck + "[solver]\nbackend = gmres_ilut\n");
    CHECK(d.options.linear.fallback == linalg::LinearFallback::sparse_lu);

    // Each key is in the run's identity with the iterative backend; with a direct backend the keys
    // do not apply and the identity is the plain deck's.
    const auto identity = [](const deck::Plan& q) {
        return solve::make_run_record(q.device, q.options, q.points).input_identity;
    };
    const deck::Plan plain = must(diode_deck);
    CHECK(identity(d) != identity(plain));
    for (const char* key : {"restart = 41", "max_linear_iterations = 501", "ilut_drop = 2e-5",
                            "ilut_fill = 21", "linear_fallback = sparse_lu"}) {
        CAPTURE(key);
        CHECK(identity(must(diode_deck + solver + key + "\n")) != identity(p));
    }
    const deck::Plan direct =
        must(diode_deck + "[solver]\nrestart = 40\nilut_fill = 20\nlinear_fallback = none\n");
    CHECK(identity(direct) == identity(plain));

    const std::string head = "[mesh]\ntype = tensor\nx = 0, 1um, 2um\n[region si]\nmaterial = silicon\n"
                             "[contact a]\nkind = ohmic\npatch = x_min\n[contact b]\nkind = ohmic\npatch = x_max\n"
                             "[equilibrium]\n[solver]\n";
    const auto expect = [&](const std::string& entry, const std::string& what) {
        const std::string message = failure(head + entry + "\n");
        INFO(message);
        CHECK(message.find(what) != std::string::npos);
    };
    expect("backend = cuda", "t.deck:14: backend is eigen, pardiso, gmres or gmres_ilut");
    expect("restart = 0", "t.deck:14: 'restart' must lie in 1..1000");
    expect("max_linear_iterations = 0", "t.deck:14: 'max_linear_iterations' must lie in 1..");
    expect("ilut_drop = 1", "t.deck:14: ilut_drop is in [0, 1)");
    expect("ilut_fill = 1001", "t.deck:14: 'ilut_fill' must lie in 1..1000");
    expect("linear_fallback = pardiso", "t.deck:14: linear_fallback is sparse_lu or none");
}

TEST_CASE("deck build: backend = gmres, the preconditioner and the amg keys (Unit 58 SD6, SQ7)", "[deck]") {
    const auto identity = [](const deck::Plan& q) {
        return solve::make_run_record(q.device, q.options, q.points).input_identity;
    };
    // gmres is gmres_ilut, with an unchanged identity; ilut is the default preconditioner.
    const deck::Plan old = must(diode_deck + "[solver]\nbackend = gmres_ilut\n");
    const deck::Plan gmres = must(diode_deck + "[solver]\nbackend = gmres\n");
    CHECK(gmres.options.linear.backend == linalg::SolverBackend::gmres_ilut);
    CHECK(gmres.options.linear.preconditioner == linalg::Preconditioner::ilut);
    CHECK(identity(gmres) == identity(old));
    CHECK(identity(must(diode_deck + "[solver]\nbackend = gmres\npreconditioner = ilut\n")) == identity(old));
    const std::string cpr = "[solver]\nbackend = gmres\npreconditioner = cpr\n";
    const deck::Plan c = must(diode_deck + cpr);
    CHECK(c.options.linear.preconditioner == linalg::Preconditioner::cpr);
    CHECK(must(diode_deck + "[solver]\nbackend = gmres\npreconditioner = cpr_amg\n").options.linear.preconditioner ==
          linalg::Preconditioner::cpr_amg);  // Unit 63 MD4
    CHECK(must(diode_deck + "[solver]\nbackend = gmres\npreconditioner = amg\n").options.linear.preconditioner ==
          linalg::Preconditioner::amg);
    CHECK(identity(c) != identity(old));
    const deck::Plan k = must(diode_deck + cpr +
                              "amg_strength_threshold = 0.1\namg_smoother_degree = 3\namg_coarse_size = 500\n"
                              "amg_max_levels = 6\nthreads = 4\n");
    const linalg::SolverConfig& l = k.options.linear;
    CHECK(l.amg_strength_threshold == 0.1);
    CHECK(l.amg_smoother_degree == 3);
    CHECK(l.amg_coarse_size == 500);
    CHECK(l.amg_max_levels == 6);
    CHECK(l.threads == 4);
    for (const char* key : {"amg_strength_threshold = 0.1", "amg_smoother_degree = 3", "amg_coarse_size = 500",
                            "amg_max_levels = 6"}) {
        CAPTURE(key);
        CHECK(identity(must(diode_deck + cpr + key + "\n")) != identity(c));
    }

    const std::string head = "[mesh]\ntype = tensor\nx = 0, 1um, 2um\n[region si]\nmaterial = silicon\n"
                             "[contact a]\nkind = ohmic\npatch = x_min\n[contact b]\nkind = ohmic\npatch = x_max\n"
                             "[equilibrium]\n[solver]\n";
    const auto expect = [&](const std::string& entry, const std::string& what) {
        const std::string message = failure(head + entry + "\n");
        INFO(message);
        CHECK(message.find(what) != std::string::npos);
    };
    expect("preconditioner = jacobi", "t.deck:14: preconditioner is ilut, amg, cpr or cpr_amg");
    expect("preconditioner = cpr", "preconditioner = amg, cpr or cpr_amg needs backend = gmres");
    expect("backend = gmres\namg_strength_threshold = 1", "t.deck:15: amg_strength_threshold is in [0, 1)");
    expect("amg_smoother_degree = 11", "t.deck:14: 'amg_smoother_degree' must lie in 1..10");
    expect("amg_coarse_size = 0", "t.deck:14: 'amg_coarse_size' must lie in 1..100000");
    expect("amg_max_levels = 1", "t.deck:14: 'amg_max_levels' must lie in 2..50");
}

TEST_CASE("deck build: rates = off (Unit 45, RQ1)", "[deck]") {
    CHECK(must(diode_deck).options.store_rates);
    CHECK_FALSE(must(diode_deck + "[solver]\nrates = off\n").options.store_rates);
    CHECK(must(diode_deck + "[solver]\nrates = on\n").options.store_rates);
    const std::string head = "[mesh]\ntype = tensor\nx = 0, 1um, 2um\n[region si]\nmaterial = silicon\n"
                             "[contact a]\nkind = ohmic\npatch = x_min\n[contact b]\nkind = ohmic\npatch = x_max\n"
                             "[equilibrium]\n[solver]\nrates = maybe\n";
    CHECK(failure(head).find("t.deck:14: 'maybe' is not on or off") != std::string::npos);
}

TEST_CASE("deck build: minority_floor (Unit 50, RD2)", "[deck]") {
    CHECK_FALSE(must(diode_deck).options.minority_floor);
    CHECK(must(diode_deck + "[solver]\nminority_floor = on\n").options.minority_floor);
    CHECK_FALSE(must(diode_deck + "[solver]\nminority_floor = off\n").options.minority_floor);
}

// The decks of decks/ (heavy devices, not in the app's examples): each still reads and builds as
// the deck format evolves. Built, not run (the FinFET sweep takes about 18 min with PARDISO).
TEST_CASE("deck build: the FinFET deck of decks/ builds") {
    const auto d = deck::read(fs::path(NITCAD_DECKS_DIR) / "finfet_3d.deck");
    if (!d) FAIL(d.error().message);
    const auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    CHECK(p->device.mesh().node_count() == 45375);
    CHECK(p->device.contacts().size() == 6);
    CHECK(p->points.size() == 17);
    CHECK(p->options.models.fermi_dirac);
    CHECK(p->options.models.surface_mobility == assemble::SurfaceMobility::off);
}

// Unit 53 (DECISIONS.md FD13, FC11): the sectioned extrusion.
namespace {

// Two stacked slabs along x, 1 um wide (y) and 1 um high (z), the interface at z = 0.5 um.
const std::string extruded_head = R"([mesh]
type = extrude3d
axis = x
stations = 0, 1um
cell a = 0 0, 1um 0, 1um 0.5um, 0 0.5um
cell b = 0 0.5um, 1um 0.5um, 1um 1um, 0 1um
patch bottom = -1um -1um -1um 2um 2um 0
patch top = -1um -1um 1um 2um 2um 2um
max_size = 0.5um
interface_resolution = 0.1um
[region lower]
material = silicon
extrude = a 0:1um
[region upper]
material = silicon
extrude = b 0:1um
[doping]
profile = acceptor constant 1e16 z=0:0.5um
profile = donor constant 1e16 z=0.5um:1um
[contact bottom]
kind = ohmic
patch = bottom
[contact top]
kind = ohmic
patch = top
[equilibrium]
)";

// The deck with `from` replaced by `to` (which must occur once).
std::string changed(std::string text, const std::string& from, const std::string& to) {
    const auto at = text.find(from);
    REQUIRE(at != std::string::npos);
    text.replace(at, from.size(), to);
    return text;
}

}  // namespace

TEST_CASE("deck build: an extrude3d mesh is a certified two-region mesh with named patches", "[deck][extrude]") {
    const deck::Plan p = must(extruded_head);
    REQUIRE(p.generated.has_value());
    const auto& m = p.device.mesh();
    CHECK(m.dimension() == 3);
    REQUIRE(p.device.interface_certificate() != nullptr);
    CHECK(p.device.regions().size() == 2);
    for (std::size_t i = 0; i < m.node_count(); ++i) {
        const double z = m.points()[i][2];
        if (std::abs(z - 0.5e-4) > 1e-9) CHECK(p.device.node_region()[i] == (z > 0.5e-4 ? 1 : 0));
    }
    for (const char* name : {"bottom", "top"}) {
        const auto* b = m.find_boundary(name);
        REQUIRE(b != nullptr);
        CHECK(std::abs(std::accumulate(b->areas.begin(), b->areas.end(), 0.0) - 1e-8) < 1e-20);
    }
    // The unnamed side faces keep the default patch name.
    CHECK(m.find_boundary("boundary") != nullptr);
}

TEST_CASE("deck build: an extrude3d position spelled in other units is the same station or box face", "[deck][extrude]") {
    // 0.3um, 300nm and 3e-5 (cm) need not round to one double. The station takes the lowest
    // spelling, the box face the highest (so an exact test would miss the facets at the station),
    // the ranges the others.
    std::vector<std::pair<double, std::string>> spellings{{0.3 * 1e-4, "0.3um"}, {300.0 * 1e-7, "300nm"}, {3e-5, "3e-5"}};
    std::ranges::sort(spellings);
    REQUIRE(spellings.front().first < spellings.back().first);
    std::string text = changed(extruded_head, "stations = 0, 1um", "stations = 0, " + spellings[0].second);
    text = changed(text, "extrude = a 0:1um", "extrude = a 0:" + spellings[1].second);
    text = changed(text, "extrude = b 0:1um", "extrude = b 0:" + spellings[2].second);
    text = changed(text, "interface_resolution = 0.1um", "interface_resolution = 0.04um");  // the 6 s rule on 0.3 um
    // A box whose lower x face is the end station: it holds the two end facets.
    text = changed(text, "patch top =", "patch far = " + spellings[2].second + " -1um -1um 1um 2um 2um\npatch top =");
    const deck::Plan p = must(text);
    const auto* far = p.device.mesh().find_boundary("far");
    REQUIRE(far != nullptr);
    CHECK(std::abs(std::accumulate(far->areas.begin(), far->areas.end(), 0.0) - 1e-8) < 1e-20);
}

TEST_CASE("deck build: every extrude3d mistake names its line", "[deck][extrude]") {
    const auto expect = [](const std::string& message, const std::string& what) {
        INFO(message);
        CHECK(message.find(what) != std::string::npos);
    };
    expect(failure(changed(extruded_head, "axis = x\n", "")), "t.deck:1: an extrude3d mesh needs 'axis'");
    expect(failure(changed(extruded_head, "stations = 0, 1um\n", "")), "needs 'stations'");
    expect(failure(changed(extruded_head, "axis = x", "axis = w")), "t.deck:3: axis is x, y or z");
    expect(failure(changed(extruded_head, "extrude = a 0:1um", "extrude = c 0:1um")), "no cell named 'c'");
    expect(failure(changed(extruded_head, "extrude = a 0:1um", "extrude = a 0:0.5um")), "the range's ends must be stations");
    expect(failure(changed(extruded_head, "extrude = a 0:1um", "extrude = a 1um:0")), "needs t_lo < t_hi");
    expect(failure(changed(extruded_head, "extrude = a 0:1um", "extrude = a 1um:1um")), "needs t_lo < t_hi");  // FC13
    expect(failure(changed(extruded_head, "cell a = 0 0, 1um 0, 1um 0.5um, 0 0.5um\ncell b = 0 0.5um, 1um 0.5um, 1um 1um, 0 1um\n", "")),
           "t.deck:1: an extrude3d mesh needs at least one 'cell'");
    expect(failure(changed(extruded_head, "patch top = -1um -1um 1um 2um 2um 2um", "patch top = -1um -1um 2um 2um 2um 1um")),
           "a box needs lo <= hi on each axis");
    expect(failure(changed(extruded_head, "extrude = b 0:1um", "extrude = b 0:1um\nextrude = a 0:1um")),
           "cell 'a' is already assigned in section 0");
    expect(failure(changed(extruded_head, "patch top = -1um -1um 1um 2um 2um 2um", "patch top = -1um -1um 3um 2um 2um 4um")),
           "t.deck:1: [mesh] extrude3d: a patch box selects no boundary facet");
    expect(failure(changed(extruded_head, "material = silicon\nextrude = a 0:1um", "material = silicon\nbox = 0 1um\nextrude = a 0:1um")),
           "an extrude3d mesh takes extrude in its regions");
    expect(failure(changed(extruded_head, "type = extrude3d\naxis = x\nstations = 0, 1um\ncell a = 0 0, 1um 0, 1um 0.5um, 0 0.5um\n"
                                          "cell b = 0 0.5um, 1um 0.5um, 1um 1um, 0 1um\n",
                           "type = generate3d\ncuboid = 0 1um 0 1um 0 1um\n")),
           "extrude is for an extrude3d mesh");
}

TEST_CASE("deck build: the run identity covers the extruded complex (FC11)", "[deck][extrude]") {
    const deck::Plan a = must(extruded_head);
    const deck::Plan again = must(extruded_head);
    const auto id = [](const deck::Plan& p) { return solve::make_run_record(p.device, p.options, p.points).input_identity; };
    CHECK(id(a) == id(again));  // deterministic: the same complex and options, the same mesh
    // A moved station and a renamed patch each change it.
    const std::string longer = changed(changed(changed(extruded_head, "stations = 0, 1um", "stations = 0, 1.02um"), "extrude = a 0:1um", "extrude = a 0:1.02um"),
                                         "extrude = b 0:1um", "extrude = b 0:1.02um");
    CHECK(id(a) != id(must(longer)));
    CHECK(id(a) != id(must(changed(changed(extruded_head, "patch top =", "patch roof ="), "patch = top", "patch = roof"))));
}

TEST_CASE("deck run: the run record names the extrude3d mesh and its generation digest (FD14)", "[deck][extrude]") {
    const fs::path dir = deck_test::scratch("extrude_generation");
    const std::string swept = changed(extruded_head, "[equilibrium]\n",
                                      "[sweep]\ncontact top = 0 : 0.05 : 0.05\n[output]\nsweep = s.nitsweep\nmesh = m.nitmesh\n");
    const deck::Plan p = must(swept, dir);
    REQUIRE(p.generated.has_value());
    CHECK(p.mesh_type == "extrude3d");
    CHECK(p.generated->generation_digest != 0);
    // max_nodes changes the generation digest, not the mesh, so not the run identity.
    const deck::Plan capped = must(changed(swept, "max_size = 0.5um\n", "max_size = 0.5um\nmax_nodes = 100000\n"), dir);
    CHECK(capped.generated->generation_digest != p.generated->generation_digest);
    CHECK(solve::make_run_record(capped.device, capped.options, capped.points).input_identity ==
          solve::make_run_record(p.device, p.options, p.points).input_identity);

    const auto report = deck::run(p);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    CHECK(report->sweep->run.mesh_type == "extrude3d");
    CHECK(report->sweep->run.mesh_generation_digest == p.generated->generation_digest);
    const auto sweep = io::load_sweep(dir / "s.nitsweep");
    if (!sweep) FAIL(sweep.error().message);
    CHECK(sweep->run.mesh_type == "extrude3d");
    CHECK(sweep->run.mesh_generation_digest == p.generated->generation_digest);
    const auto mesh = io::load_generated_mesh(dir / "m.nitmesh");
    if (!mesh) FAIL(mesh.error().message);
    CHECK(mesh->generation_digest == p.generated->generation_digest);
    // The saved mesh read back as a file mesh: the record names the generator its digest belongs to.
    const std::string tail = extruded_head.substr(extruded_head.find("[doping]"));
    const deck::Plan reread = must("[mesh]\ntype = file\npath = m.nitmesh\n[region lower]\nmaterial = silicon\n[region upper]\nmaterial = silicon\n" +
                                       changed(tail, "[equilibrium]\n", "[sweep]\ncontact top = 0 : 0.05 : 0.05\n[output]\nsweep = r.nitsweep\n"),
                                   dir);
    CHECK(reread.mesh_type == "extrude3d");
    const auto rerun = deck::run(reread);
    if (!rerun) FAIL(rerun.error().message);
    CHECK(rerun->sweep->run.mesh_type == "extrude3d");
    CHECK(rerun->sweep->run.mesh_generation_digest == p.generated->generation_digest);
    // A run on any other mesh records neither.
    const auto plain = deck::run(must(deck_test::diode_device + "[sweep]\ncontact anode = 0 : 0.1 : 0.1\n"));
    REQUIRE(plain.has_value());
    CHECK(plain->sweep->run.mesh_type.empty());
    CHECK(plain->sweep->run.mesh_generation_digest == 0);
}


// Unit 54 (DECISIONS.md UD9-UD11, UC10): the swept mesh from a deck.
namespace {

// The two slabs swept: five lines, a layer straddling the interface z = 0.5 um (u, v = y, z).
// (Built without changed(): its REQUIRE may not run before main.)
const std::string swept_head = [] {
    std::string t = extruded_head;
    const std::string after = "interface_resolution = 0.1um\n";
    t.insert(t.find(after) + after.size(),
             "sweep = on\nlines = linspace(0, 1um, 5)\n"
             "layer = 0 0.5um, 1um 0.5um along=linspace(0, 1um, 5) lines=-0.3um, -0.2um + 0.2um, 0.3um straddle=0.05um\n");
    return t;
}();

// What the swept deck says, built by hand.
std::pair<mesh::ExtrudedGeometry, mesh::SweptOptions> swept_by_hand() {
    const double um = 1e-4;
    mesh::ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1 * um};
    x.cells = {{"a", {{0, 0, 0}, {1 * um, 0, 0}, {1 * um, 0.5 * um, 0}, {0, 0.5 * um, 0}}},
               {"b", {{0, 0.5 * um, 0}, {1 * um, 0.5 * um, 0}, {1 * um, 1 * um, 0}, {0, 1 * um, 0}}}};
    x.regions = {{0}, {1}};
    x.patches = {{"bottom", {-1 * um, -1 * um, -1 * um}, {2 * um, 2 * um, 0}}, {"top", {-1 * um, -1 * um, 1 * um}, {2 * um, 2 * um, 2 * um}}};
    mesh::SweptOptions o;
    o.section.max_size = 0.5 * um;
    o.section.interface_resolution = 0.1 * um;
    const auto linspace = [](double a, double b, int n) {
        std::vector<double> v;
        for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(n - 1));
        return v;
    };
    o.lines = linspace(0, 1 * um, 5);
    mesh::BoundaryLayer l;
    l.a = {0, 0.5 * um, 0};
    l.b = {1 * um, 0.5 * um, 0};
    l.along = linspace(0, 1 * um, 5);
    l.lines = {-0.3 * um, -0.2 * um, 0.2 * um, 0.3 * um};
    l.straddle = 0.05 * um;
    o.section.layers = {l};
    return {x, o};
}

}  // namespace

TEST_CASE("deck build: sweep = on gives the swept mesh of the deck's extrusion and options (UC10)", "[deck][extrude][swept]") {
    const deck::Plan p = must(swept_head);
    REQUIRE(p.generated.has_value());
    CHECK(p.mesh_type == "extrude3d");
    const auto [x, o] = swept_by_hand();
    const auto m = mesh::generate_swept(x, o);
    if (!m) FAIL(m.error().message);
    // Every value read as written: the same digest, so the same inputs, and the same mesh.
    CHECK(p.generated->generation_digest == mesh::generation_digest(x, o));
    CHECK(mesh::graph_digest(p.device.mesh()) == mesh::graph_digest(m->mesh));
    CHECK(p.generated->node_region == m->node_region);
    REQUIRE(p.device.interface_certificate() != nullptr);
    const auto& r = p.generated->report;
    CHECK(r.axis_lines == 5);
    CHECK(r.prisms * 3 == r.tetrahedra);
    CHECK(r.prisms % 4 == 0);  // the cross-section's triangles times four intervals
    CHECK(r.layers == 1);
    CHECK(r.corner_blocks == 0);
    CHECK(p.device.regions().size() == 2);
    // Each layer line, with the straddles at +-0.025 um, is a row of nodes at its z.
    for (const double z : {0.2, 0.3, 0.475, 0.525, 0.7, 0.8}) {
        const auto n = std::ranges::count_if(p.device.mesh().points(), [&](const mesh::Point& q) { return std::abs(q[2] - z * 1e-4) < 1e-16; });
        CAPTURE(z);
        CHECK(n == 25);  // 5 along the base times 5 axis lines
    }
    // A box in the cross-section refines it; sweep = off is Unit 53's generator, unchanged (UC1).
    const deck::Plan boxed = must(changed(swept_head, "max_size = 0.5um\n", "max_size = 0.5um\nbox = 0 0 1um 0.2um 0.1um\n"));
    CHECK(boxed.device.mesh().node_count() > p.device.mesh().node_count());
    const deck::Plan off = must(changed(extruded_head, "max_size = 0.5um\n", "max_size = 0.5um\nsweep = off\n"));
    CHECK(off.generated->report.prisms == 0);
    CHECK(off.generated->generation_digest == must(extruded_head).generated->generation_digest);
}

TEST_CASE("deck build: every swept-mesh mistake names its line (UC2, UD10)", "[deck][extrude][swept]") {
    const auto expect = [](const std::string& message, const std::string& what) {
        INFO(message);
        CHECK(message.find(what) != std::string::npos);
    };
    // The swept keys without sweep = on, each on line 10.
    for (const char* entry : {"lines = 0, 1um", "straddle = 0.5um 0.1um", "layer = 0 0.5um, 1um 0.5um along=0, 1um lines=-0.1um, 0.1um",
                              "min_angle = 25"}) {
        const std::string key = std::string(entry).substr(0, std::string(entry).find(' '));
        expect(failure(changed(extruded_head, "max_size = 0.5um\n", "max_size = 0.5um\n" + std::string(entry) + "\n")),
               "t.deck:10: '" + key + "' is for a swept mesh (sweep = on)");
    }
    expect(failure(changed(swept_head, "sweep = on\n", "sweep = on\nradius_edge = 2\n")), "t.deck:12: radius_edge is for sweep = off");
    expect(failure(changed(swept_head, "lines = linspace(0, 1um, 5)\n", "")), "t.deck:1: a swept extrude3d mesh (sweep = on) needs 'lines'");
    expect(failure(changed(swept_head, "lines = linspace(0, 1um, 5)\n", "lines = 0, 1um\nlines = 0, 1um\n")), "t.deck:13: a second 'lines'");
    expect(failure(changed(swept_head, "sweep = on", "sweep = maybe")), "t.deck:11: 'maybe' is not on or off");
    expect(failure(changed(swept_head, "sweep = on\n", "sweep = on\nbox = 0 0 0 1um 1um 1um 0.1um\n")), "t.deck:12: 'box' needs 5 lengths, has 7");
    expect(failure(changed(swept_head, "sweep = on\n", "sweep = on\nstraddle = 0.5um\n")), "t.deck:12: 'straddle' needs 2 lengths, has 1");
    expect(failure(changed(swept_head, "lines = linspace(0, 1um, 5)", "lines = linspace(0, 1um, 1)")), "t.deck:12: linspace needs a count");
    // The layer's form, on line 13.
    const std::string layer = swept_head.substr(swept_head.find("layer = "), swept_head.find('\n', swept_head.find("layer = ")) - swept_head.find("layer = "));
    const auto with_layer = [&](const std::string& l) { return failure(changed(swept_head, layer, l)); };
    expect(with_layer("layer = 0 0.5um, 1um 0.5um"), "t.deck:13: a layer is '<u v>, <u v> along=");
    expect(with_layer("layer = 0 0.5um along=0, 1um lines=-0.1um, 0.1um"), "t.deck:13: a layer's base is two points");
    expect(with_layer("layer = 0 0.5um, 1um along=0, 1um lines=-0.1um, 0.1um"), "t.deck:13: a point needs 2 coordinates");
    expect(with_layer("layer = 0 0.5um, 1um 0.5um along=0, 1um lines=-0.1um, 0.1um width=1um"), "t.deck:13: unknown layer option 'width='");
    expect(with_layer("layer = 0 0.5um, 1um 0.5um along=0, 1um along=0, 1um lines=-0.1um, 0.1um"), "t.deck:13: a second 'along=' in a layer");
    expect(with_layer("layer = 0 0.5um, 1um 0.5um along=0, 1um"), "t.deck:13: a layer needs along= and lines=");
    expect(with_layer("layer = 0 0.5um, 1um 0.5um along= lines=-0.1um, 0.1um"), "t.deck:13: the layer option 'along=' needs a value");
    expect(with_layer("layer = 0 0.5um, 1um 0.5um along=0, 1um lines=-0.1qq, 0.1um"), "t.deck:13: '-0.1qq' is not a number");
    expect(with_layer("layer = 0 0.5um, 1um 0.5um along=0, 1um lines=-0.1um, 0.1um straddle=wide"), "t.deck:13: 'wide' is not a number");
    // The generators' refusals, located at [mesh]: the sweep's own and the cross-section's.
    const auto at_mesh = failure(changed(swept_head, "sweep = on\n", "sweep = on\nstraddle = 0 0.1um\n"));
    expect(at_mesh, "t.deck:1: [mesh] extrude3d (sweep): ");
    const auto section = with_layer("layer = 0 0.5um, 1um 0.5um along=0, 1um lines=-0.6um, 0.3um");
    expect(section, "t.deck:1: [mesh] extrude3d (sweep): the cross-section: ");
}

TEST_CASE("deck run: a swept mesh's run record, files and identity (UC10, UD11)", "[deck][extrude][swept]") {
    const fs::path dir = deck_test::scratch("swept_generation");
    const std::string run_text = changed(swept_head, "[equilibrium]\n",
                                         "[sweep]\ncontact top = 0 : 0.05 : 0.05\n[output]\nsweep = s.nitsweep\nmesh = m.nitmesh\n");
    const deck::Plan p = must(run_text, dir);
    REQUIRE(p.generated.has_value());
    const auto id = [](const deck::Plan& q) { return solve::make_run_record(q.device, q.options, q.points).input_identity; };
    // The identity covers the lines, the straddles and the layers through the mesh they make; max_nodes
    // changes the generation digest only.
    CHECK(id(p) == id(must(run_text, dir)));
    CHECK(id(p) != id(must(changed(run_text, "lines = linspace(0, 1um, 5)", "lines = linspace(0, 1um, 6)"), dir)));
    CHECK(id(p) != id(must(changed(run_text, "straddle=0.05um", "straddle=0.04um"), dir)));
    CHECK(id(p) != id(must(changed(run_text, "lines=-0.3um, -0.2um", "lines=-0.35um, -0.2um"), dir)));
    const deck::Plan capped = must(changed(run_text, "sweep = on\n", "sweep = on\nmax_nodes = 100000\n"), dir);
    CHECK(capped.generated->generation_digest != p.generated->generation_digest);
    CHECK(id(capped) == id(p));

    const auto report = deck::run(p);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    CHECK(report->sweep->run.mesh_type == "extrude3d");
    CHECK(report->sweep->run.mesh_generation_digest == p.generated->generation_digest);
    const auto sweep = io::load_sweep(dir / "s.nitsweep");
    if (!sweep) FAIL(sweep.error().message);
    CHECK(sweep->run.mesh_generation_digest == p.generated->generation_digest);
    // The mesh file: the swept report, regions and certificate come back.
    const auto mesh = io::load_generated_mesh(dir / "m.nitmesh");
    if (!mesh) FAIL(mesh.error().message);
    CHECK(mesh->generation_digest == p.generated->generation_digest);
    CHECK(mesh->report.prisms == p.generated->report.prisms);
    CHECK(mesh->report.axis_lines == 5);
    CHECK(mesh->report.layers == 1);
    CHECK(mesh->node_region == p.generated->node_region);
    REQUIRE(mesh->certificate.has_value());
    // Read back as a file mesh, the run names the generator its digest belongs to.
    const std::string tail = swept_head.substr(swept_head.find("[doping]"));
    const deck::Plan reread = must("[mesh]\ntype = file\npath = m.nitmesh\n[region lower]\nmaterial = silicon\n[region upper]\nmaterial = silicon\n" +
                                       changed(tail, "[equilibrium]\n", "[sweep]\ncontact top = 0 : 0.05 : 0.05\n[output]\nsweep = r.nitsweep\n"),
                                   dir);
    CHECK(reread.mesh_type == "extrude3d");
    CHECK(id(reread) == id(p));
    const auto rerun = deck::run(reread);
    if (!rerun) FAIL(rerun.error().message);
    CHECK(rerun->sweep->run.mesh_generation_digest == p.generated->generation_digest);
}

// Unit 54 slice 4: the two FinFET decks as swept decks (UC8's and UC9's meshes) build and are certified,
// with P6 where the gate patch ends at the straddled stations; a few seconds each in Release but minutes in Debug
// (UC13), so hidden; run it by name.
TEST_CASE("deck build: the swept FinFET decks mesh and are certified", "[.benchmark][deck][extrude][swept]") {
    for (const auto& [name, axis_lines] :
         {std::pair{"finfet_3d_swept.deck", std::size_t{25}}, std::pair{"finfet_3d_shaped_swept.deck", std::size_t{37}}}) {
        CAPTURE(name);
        const auto d = deck::read(fs::path(NITCAD_DECKS_DIR) / name);
        if (!d) FAIL(d.error().message);
        const auto p = deck::build(*d);
        if (!p) FAIL(p.error().message);
        REQUIRE(p->generated.has_value());
        const auto& r = p->generated->report;
        CHECK(r.layers == 4);
        CHECK(r.corner_blocks == 4);
        CHECK(r.axis_lines == axis_lines);
        const auto* c = p->device.interface_certificate();
        REQUIRE(c != nullptr);
        CHECK(c->patches.has_value());
        CHECK(p->device.contacts().size() == 4);
        CHECK(p->device.regions().size() == 5);
        UNSCOPED_INFO(name << ": nodes " << p->device.mesh().node_count() << ", cross-section " << r.nodes / r.axis_lines
                           << ", P4 flat " << c->report.max_flat_deviation << ", corner " << c->report.max_corner_deviation);
    }
}
// Unit 53 slice 5: the extruded FinFET deck of decks/ meshes and is certified (P4 and P6). About twenty
// minutes in Release, so hidden; run it by name.
TEST_CASE("deck build: the reduced extruded FinFET deck meshes and is certified", "[.benchmark][deck][extrude]") {
    const auto d = deck::read(fs::path(NITCAD_DECKS_DIR) / "finfet_3d_reduced.deck");
    if (!d) FAIL(d.error().message);
    const auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    REQUIRE(p->generated.has_value());
    const auto* c = p->device.interface_certificate();
    REQUIRE(c != nullptr);
    CHECK(c->patches.has_value());
    CHECK(p->device.contacts().size() == 4);
    CHECK(p->device.mesh().find_boundary("gate") != nullptr);
    UNSCOPED_INFO("nodes " << p->device.mesh().node_count() << ", slivers " << p->generated->report.slivers);
}

// The tapered deck (FC10's mesh): certified since the owner's ruling on slice 5 (DECISIONS.md: the P4b bound at the
// interfacial ring's oblique corners). About twenty minutes in Release, so hidden; run it by name.
TEST_CASE("deck build: the shaped extruded FinFET deck meshes and is certified", "[.benchmark][deck][extrude]") {
    const auto d = deck::read(fs::path(NITCAD_DECKS_DIR) / "finfet_3d_shaped.deck");
    if (!d) FAIL(d.error().message);
    const auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    CHECK(p->device.interface_certificate() != nullptr);
}

// Unit 26l slice 0 (SD1): a layered generate3d deck.
TEST_CASE("deck build: a layered generate3d deck builds; a prism and a box of one solid mesh alike (Unit 26l)", "[deck]") {
    const std::string mesh = R"([mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
patch gate = 0.3um 0 0.6um  0.7um 0.6um 0.6um
patch body = 0 0 0  1um 0.6um 0
max_size = 0.15um
interface_resolution = 0.02um
)";
    const std::string tail = R"(
[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[contact gate]
kind = electrode
patch = gate
electrode = metal
work_function = 4.4

[equilibrium]
)";
    const std::string boxed = mesh + R"(
[region si]
material = silicon

[region ox]
material = sio2
box = 0 1um 0 0.6um 0.4um 0.6um
)" + tail;
    const deck::Plan a = must(boxed);
    CHECK(a.device.regions().size() == 2);
    CHECK(a.device.mesh().find_boundary("gate") != nullptr);
    CHECK(a.device.mesh().find_boundary("x_min") != nullptr);  // the faces not in a patch box: by name
    std::string prism = boxed;
    prism.replace(prism.find("box = 0 1um 0 0.6um 0.4um 0.6um"), 31, "prism = z 0.4um 0.6um, 0 0, 1um 0, 1um 0.6um, 0 0.6um");
    const deck::Plan b = must(prism);
    CHECK(io::mesh_digest(a.device.mesh()) == io::mesh_digest(b.device.mesh()));
    std::string along_x = boxed;  // the same oxide as a prism along x, its polygon in (y, z)
    along_x.replace(along_x.find("box = 0 1um 0 0.6um 0.4um 0.6um"), 31, "prism = x 0 1um, 0 0.4um, 0.6um 0.4um, 0.6um 0.6um, 0 0.6um");
    CHECK(io::mesh_digest(a.device.mesh()) == io::mesh_digest(must(along_x).device.mesh()));

    const auto fails = [](const std::string& text, std::string_view message) {
        const std::string m = failure(text);
        CAPTURE(m);
        CHECK(m.find(message) != std::string::npos);
    };
    std::string t = boxed;
    fails(t.replace(t.find("layered = on\n"), 13, ""), "t.deck:3: domain is for a layered generate3d mesh (layered = on)");
    t = boxed;
    fails(t.replace(t.find("domain = 0 1um 0 0.6um 0 0.6um\n"), 31, "cuboid = 0 1um 0 0.6um 0 0.6um\n"),
          "t.deck:4: 'cuboid' is for generate3d without layers: a layered one takes domain");
    t = boxed;
    fails(t.replace(t.find("domain = 0 1um 0 0.6um 0 0.6um\n"), 31, ""), "a layered generate3d mesh needs a domain");
    t = boxed;  // a slanted oxide builds since slice 4 (SD6): a trapezoid along y, its polygon in (z, x)
    t.replace(t.find("box = 0 1um 0 0.6um 0.4um 0.6um"), 31, "prism = y 0 0.6um, 0.4um 0.1um, 0.4um 0.9um, 0.6um 0.8um, 0.6um 0.2um");
    CHECK(must(t).device.regions().size() == 2);
    t = prism;
    fails(t.replace(t.find("prism = z"), 9, "prism = w"), "t.deck:15: a prism is '<axis x|y|z> <from> <to>, <u v>, <u v>, ...'");
    t = boxed;
    fails(t.replace(t.find("patch gate = 0.3um 0 0.6um  0.7um 0.6um 0.6um"), 45, "patch gate = 0.3um 0.1um 0.3um  0.7um 0.5um 0.5um"),
          "t.deck:5: [mesh] layered: a patch box selects no boundary facet");
    t = boxed;
    fails(t.replace(t.find("[region si]\nmaterial = silicon\n"), 30, "[region si]\nmaterial = silicon\nseed = 0 0 0\n"),
          "a layered generate3d mesh takes box and prism in its regions");
    fails(R"([mesh]
type = tensor
x = 0 1um
[region si]
material = silicon
prism = z 0 1um, 0 0, 1um 0, 1um 1um
[equilibrium]
)", "t.deck:6: prism is for a layered generate3d mesh (layered = on)");
}

// Unit 26l (SE2): the reduced FinFET written as layers (decks/finfet_3d_reduced.deck's [mesh] and regions redrawn as
// boxes); the complex's volumes and patch areas equal the extrusion's (the mesh test). Meshing takes about twenty
// minutes in Release, so hidden; run it by name.
TEST_CASE("deck build: the reduced FinFET as layered generate3d meshes and is certified (Unit 26l)", "[.benchmark][deck][layered3d]") {
    std::ifstream in(fs::path(NITCAD_DECKS_DIR) / "finfet_3d_reduced.deck", std::ios::binary);
    REQUIRE(in);
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const auto mesh_at = text.find("[mesh]"), mesh_end = text.find("[region buried_oxide]");
    text.replace(mesh_at, mesh_end - mesh_at, R"([mesh]
type = generate3d
layered = on
domain = 0 85nm 0 52nm 0 14nm
patch source = 0 20nm 2nm  0 50nm 12nm
patch drain = 85nm 20nm 2nm  85nm 50nm 12nm
patch gate = 30nm 20nm 0  55nm 52nm 14nm
patch substrate = 0 0 0  85nm 0 14nm
max_size = 8nm
interface_resolution = 0.08nm

)");
    for (const auto& [from, to] : std::vector<std::pair<std::string, std::string>>{
             {"extrude = box 0:85nm\n", ""},
             {"extrude = il hk 0:30nm\nextrude = il hk 55nm:85nm\n", "box = 0 30nm 20nm 52nm 0 14nm\nbox = 55nm 85nm 20nm 52nm 0 14nm\n"},
             {"extrude = hk 30nm:55nm\n", "box = 30nm 55nm 20nm 52nm 0 14nm\n"},
             {"extrude = il 30nm:55nm\n", "box = 30nm 55nm 20nm 50.5nm 1.5nm 12.5nm\n"},
             {"extrude = fin 0:85nm\n", "box = 0 85nm 20nm 50nm 2nm 12nm\n"}}) {
        const auto at = text.find(from);
        REQUIRE(at != std::string::npos);
        text.replace(at, from.size(), to);
    }
    const auto p = plan_of(text);
    if (!p) FAIL(p.error().message);
    CHECK(p->device.interface_certificate() != nullptr);
    CHECK(p->device.contacts().size() == 4);
}

TEST_CASE("deck build: overlap = fill takes only what no earlier shape covers (Unit 26m MT4, MC5)", "[deck][overlap]") {
    const std::string head = R"(
[mesh]
type = generate2d
layered = on
domain = 0 2um 0 1um
patch source = 0 0, 0 0.6um
patch body = 0.2um 0, 0.8um 0
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon

[region n]
material = silicon
polygon = 1um 0, 2um 0, 2um 0.5um, 2um 0.8um, 1um 0.8um
)";
    const std::string tail = R"(
[doping]
profile = acceptor constant 1e16

[contact source]
kind = ohmic
patch = source

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";
    const auto same = [](const deck::Plan& a, const deck::Plan& b) {
        return io::mesh_digest(a.device.mesh()) == io::mesh_digest(b.device.mesh()) && std::ranges::equal(a.device.node_region(), b.device.node_region());
    };
    // Over the n well too, filling: only the p part is taken, as the box cut to it.
    const deck::Plan fill = must(head + "\n[region cap]\nmaterial = sio2\noverlap = fill\nbox = 0.5um 1.5um 0.6um 0.8um\n" + tail);
    const deck::Plan cut = must(head + "\n[region cap]\nmaterial = sio2\nbox = 0.5um 1um 0.6um 0.8um\n" + tail);
    const deck::Plan over = must(head + "\n[region cap]\nmaterial = sio2\nbox = 0.5um 1.5um 0.6um 0.8um\n" + tail);
    CHECK(same(fill, cut));
    CHECK_FALSE(same(fill, over));
    // Over the background only: the same bits as replace; `overlap = replace` the default.
    const deck::Plan a = must(head + "\n[region cap]\nmaterial = sio2\noverlap = fill\nbox = 0 1um 0.8um 1um\n" + tail);
    const deck::Plan b = must(head + "\n[region cap]\nmaterial = sio2\nbox = 0 1um 0.8um 1um\n" + tail);
    const deck::Plan c = must(head + "\n[region cap]\nmaterial = sio2\noverlap = replace\nbox = 0 1um 0.8um 1um\n" + tail);
    CHECK(same(a, b));
    CHECK(same(b, c));
    // The errors, on their lines.
    const auto fails = [](const std::string& text, std::string_view message) {
        const auto p = plan_of(text);
        REQUIRE_FALSE(p.has_value());
        CAPTURE(p.error().message);
        CHECK(p.error().message.find(message) != std::string::npos);
    };
    fails(head + "\n[region cap]\nmaterial = sio2\noverlap = under\nbox = 0 1um 0.8um 1um\n" + tail, ":20: overlap is 'replace' or 'fill'");
    fails(head + "\n[region cap]\nmaterial = sio2\noverlap = fill\noverlap = fill\nbox = 0 1um 0.8um 1um\n" + tail, ":21: a second 'overlap'");
    fails("[mesh]\ntype = tensor\nx = 0, 1um, 2um\n[region si]\nmaterial = silicon\noverlap = fill\n[contact a]\nkind = ohmic\npatch = x_min\n"
          "[contact b]\nkind = ohmic\npatch = x_max\n[equilibrium]\n",
          ":6: overlap is for a layered mesh");
}

// Unit 68 (MAD4, MAC5): user materials, inline and from a .nitmat file.
TEST_CASE("deck build: [material] user materials (MAC5)", "[deck][unit68]") {
    const std::string head = "[mesh]\ntype = tensor\nx = 0, 1um\n";
    const auto region_parameters = [](const deck::Plan& p) {
        return std::get<physics::Semiconductor>(p.device.regions()[0].material).parameters();
    };
    // A deck that must build: an ohmic contact added before its analysis.
    const auto built = [](std::string text, const fs::path& dir = {}) {
        text.insert(text.find("[equilibrium]"), "[contact a]\nkind = ohmic\npatch = x_min\n");
        return must(text, dir);
    };
    SECTION("a copy of silicon with no overrides is silicon bit for bit") {
        const deck::Plan p = built(head + "[material mysi]\nbase = silicon\n[region si]\nmaterial = mysi\n[equilibrium]\n");
        const deck::Plan q = built(head + "[region si]\nmaterial = silicon\n[equilibrium]\n");
        CHECK(region_parameters(p) == region_parameters(q));
        CHECK(region_parameters(p) == physics::silicon_parameters);
    }
    SECTION("an override changes exactly that member; a region may override again") {
        const deck::Plan p = built(head + "[material mysi]\nbase = silicon\neps_r = 12.0\nhole_mobility.mu_max = 400\n"
                                         "[region si]\nmaterial = mysi\nlifetime.tau_n0 = 1e-7\n[equilibrium]\n");
        physics::SemiconductorParameters want = physics::silicon_parameters;
        want.eps_r = 12.0;
        want.hole_mobility.mu_max = 400;
        want.lifetime.tau_n0 = 1e-7;
        CHECK(region_parameters(p) == want);
    }
    SECTION("a base with composition, and an insulator base") {
        const deck::Plan p = built(head + "[material barrier]\nbase = algaas 0.3\n[region a]\nmaterial = barrier\n[equilibrium]\n");
        CHECK(region_parameters(p) == *physics::algaas_parameters(0.3));
        const deck::Plan q = built("[mesh]\ntype = tensor\nx = 0, 1um, 2um\n[material hk]\nbase = insulator\neps_r = 22\n"
                                  "[region si]\nmaterial = silicon\nbox = 0 1um\n[region ox]\nmaterial = hk\nbox = 1um 2um\n"
                                  "[equilibrium]\n");
        CHECK(std::get<physics::Insulator>(q.device.regions()[1].material).parameters().eps_r == 22.0);
    }
    SECTION("an include file equals the inline section") {
        const fs::path dir = scratch("user_material");
        {
            std::ofstream f(dir / "mysi.nitmat", std::ios::binary);
            f << "\xEF\xBB\xBF# a user material\nbase = silicon\neps_r = 12.0\nhole_mobility.mu_max = 400\n";
        }
        const deck::Plan p = built(head + "[material mysi]\ninclude = mysi.nitmat\n[region si]\nmaterial = mysi\n[equilibrium]\n", dir);

        const deck::Plan q = built(head + "[material mysi]\nbase = silicon\neps_r = 12.0\nhole_mobility.mu_max = 400\n"
                                         "[region si]\nmaterial = mysi\n[equilibrium]\n");
        CHECK(region_parameters(p) == region_parameters(q));
        {
            std::ofstream f(dir / "bad.nitmat", std::ios::binary);
            f << "base = silicon\n\nbogus = 1\n";
        }
        CHECK(failure(head + "[material m]\ninclude = bad.nitmat\n[region si]\nmaterial = m\n[equilibrium]\n", dir)
                  .find("bad.nitmat:3: unknown key 'bogus' in [material m]") != std::string::npos);
        {
            std::ofstream f(dir / "nested.nitmat", std::ios::binary);
            f << "base = silicon\ninclude = mysi.nitmat\n";
        }
        CHECK(failure(head + "[material m]\ninclude = nested.nitmat\n[region si]\nmaterial = m\n[equilibrium]\n", dir)
                  .find("nested.nitmat:2: an included material file cannot include another") != std::string::npos);
        {
            std::ofstream f(dir / "sections.nitmat", std::ios::binary);
            f << "base = silicon\n[region x]\n";
        }
        CHECK(failure(head + "[material m]\ninclude = sections.nitmat\n[region si]\nmaterial = m\n[equilibrium]\n", dir)
                  .find("t.deck:5: include: a material file holds key lines only") != std::string::npos);
        CHECK(failure(head + "[material m]\ninclude = none.nitmat\n[region si]\nmaterial = m\n[equilibrium]\n", dir)
                  .find("t.deck:5: cannot open the material file") != std::string::npos);
    }
    SECTION("refusals name the key or section") {
        const auto expect = [](const std::string& message, const std::string& what) {
            INFO(message);
            CHECK(message.find(what) != std::string::npos);
        };
        const std::string tail = "[region si]\nmaterial = m\n[equilibrium]\n";
        expect(failure(head + "[material silicon]\nbase = silicon\n[region si]\nmaterial = silicon\n[equilibrium]\n"),
               "t.deck:4: [material silicon]: 'silicon' is a registry material");
        expect(failure(head + "[material m]\nbase = silicon\n[material m]\nbase = gaas\n" + tail), "t.deck:6: a second [material m]");
        expect(failure(head + "[material m]\neps_r = 3\n" + tail), "t.deck:4: [material m] has no 'base'");
        expect(failure(head + "[material m]\nbase = silicon\nbase = gaas\n" + tail), "t.deck:6: a second 'base'");
        expect(failure(head + "[material m]\nbase = m2\n" + tail), "t.deck:5: unknown material 'm2'");
        expect(failure(head + "[material m]\nbase = silicon\nbox = 0 1um\n" + tail),
               "t.deck:6: unknown key 'box' in [material m] (not a parameter of its base)");
        expect(failure(head + "[material m]\nbase = silicon\neps_r = -1\n" + tail), "t.deck:5: [material m]: ");
        expect(failure(head + "[material m]\nbase = silicon\neps_r = -1\n" + tail), "eps_r");
        expect(failure(head + "[material m]\nbase = silicon\ninclude = a.nitmat\n" + tail),
               "t.deck:6: a [material] with 'include' takes no other key");
        expect(failure(head + "[material m]\nbase = silicon\n[region si]\nmaterial = m 0.3\n[equilibrium]\n"),
               "t.deck:7: 'm' takes no arguments");
        expect(failure(head + "[material]\nbase = silicon\n" + tail), "t.deck:4: [material] needs a name");
    }
}

// Unit 68 (MAD3, MAC6): the plan records each region's user material; the pre-Unit-68 names are
// not recorded (their result files stay byte-identical, MAC1).
TEST_CASE("deck build: user materials are recorded with their composition and hash (MAC6)", "[deck][unit68]") {
    const std::string head = "[mesh]\ntype = tensor\nx = 0, 1um\n";
    const std::string tail = "[contact a]\nkind = ohmic\npatch = x_min\n[equilibrium]\n";
    const deck::Plan p = must(head + "[material barrier]\nbase = algaas 0.3\n[region a]\nmaterial = barrier\n" + tail);
    REQUIRE(p.materials.size() == 1);
    CHECK(p.materials[0].region == "a");
    CHECK(p.materials[0].material == "barrier");
    CHECK(p.materials[0].composition == std::vector<double>{0.3});
    const auto& used = std::get<physics::Semiconductor>(p.device.regions()[0].material).parameters();
    CHECK(p.materials[0].parameter_hash == physics::parameter_hash(physics::MaterialParameters{used}));
    // Two decks differing only in x differ in hash; a region override changes it too.
    const deck::Plan q = must(head + "[material barrier]\nbase = algaas 0.31\n[region a]\nmaterial = barrier\n" + tail);
    CHECK(q.materials[0].parameter_hash != p.materials[0].parameter_hash);
    const deck::Plan r = must(head + "[material barrier]\nbase = algaas 0.3\n[region a]\nmaterial = barrier\neps_r = 12\n" + tail);
    CHECK(r.materials[0].parameter_hash != p.materials[0].parameter_hash);
    // The registry names before Unit 68 are not recorded, with or without composition.
    CHECK(must(head + "[region a]\nmaterial = algaas 0.3\n" + tail).materials.empty());
    CHECK(must(head + "[region a]\nmaterial = silicon\n" + tail).materials.empty());
}

TEST_CASE("deck build: an unknown material suggests the nearest name (MAC13)", "[deck][unit68]") {
    const std::string head = "[mesh]\ntype = tensor\nx = 0, 1um\n";
    CHECK(failure(head + "[region a]\nmaterial = silicn\n[equilibrium]\n")
              .find("t.deck:5: unknown material 'silicn' (did you mean 'silicon'?)") != std::string::npos);
    CHECK(failure(head + "[material barrier]\nbase = algaas 0.3\n[region a]\nmaterial = barier\n[equilibrium]\n")
              .find("unknown material 'barier' (did you mean 'barrier'?)") != std::string::npos);
    CHECK(failure(head + "[region a]\nmaterial = unobtainium\n[equilibrium]\n").find("did you mean") == std::string::npos);
}

// Unit 68 (MAD14, slice 5a): `metal = <preset>` on a contact is its preset's work function or
// polysilicon electrode, exactly as written out; the refusals name the reason.
TEST_CASE("deck build: metal presets on contacts (MAD14)", "[deck][unit68]") {
    const auto deck = [](const std::string& kind, const std::string& keys) {
        return "[mesh]\ntype = tensor\nx = linspace(0, 1um, 11)\n[region si]\nmaterial = silicon\n"
               "[doping]\nprofile = donor constant 1e16\n[contact a]\nkind = " + kind + "\npatch = x_min\n" + keys +
               "[contact b]\nkind = ohmic\npatch = x_max\n[equilibrium]\n";
    };
    const auto contact = [](const deck::Plan& p) { return p.device.contacts().front(); };
    const auto same = [&](const std::string& kind, const std::string& preset, const std::string& written) {
        const device::Contact a = contact(must(deck(kind, "metal = " + preset + "\n")));
        const device::Contact b = contact(must(deck(kind, written)));
        CHECK(a.gate.electrode == b.gate.electrode);
        CHECK(a.gate.work_function_eV == b.gate.work_function_eV);
        CHECK(a.schottky.by_work_function == b.schottky.by_work_function);
        CHECK(a.schottky.work_function_eV == b.schottky.work_function_eV);
    };
    same("schottky", "pt", "work_function = 5.65\n");
    same("schottky", "tin", "work_function = 4.70\n");
    same("gate", "al\noxide_thickness = 2nm", "oxide_thickness = 2nm\nelectrode = metal\nwork_function = 4.28\n");
    same("gate", "p_poly\noxide_thickness = 2nm", "oxide_thickness = 2nm\nelectrode = p_poly\n");
    CHECK(contact(must(deck("gate", "oxide_thickness = 2nm\nmetal = n_poly\n"))).gate.electrode == device::GateElectrode::n_poly);
    CHECK(failure(deck("schottky", "metal = pt\nwork_function = 5\n")).find("'metal' replaces") != std::string::npos);
    CHECK(failure(deck("gate", "oxide_thickness = 2nm\nmetal = al\nelectrode = metal\n")).find("'metal' replaces") != std::string::npos);
    CHECK(failure(deck("ohmic", "metal = al\n")).find("'metal' is for a gate, electrode or Schottky contact") != std::string::npos);
    CHECK(failure(deck("schottky", "metal = n_poly\n")).find("polysilicon electrode, not a Schottky metal") != std::string::npos);
    CHECK(failure(deck("schottky", "metal = tln\n")).find("unknown metal 'tln' (did you mean 'tin'?)") != std::string::npos);
}

// Unit 68 (MAD6, MAC7): a gate written as a metal region equals the same gate as an electrode patch:
// the potential at every node they share agrees to 1e-9 V, at equilibrium and under gate bias.
TEST_CASE("deck build: a metal-region gate equals its electrode patch (MAC7)", "[deck][unit68]") {
    const std::string body = "[region si]\nmaterial = silicon\n[region ox]\nmaterial = sio2\nbox = 1um 1.01um\n"
                             "[doping]\nprofile = acceptor constant 1e17\n";
    const std::string tail = "[contact b]\nkind = ohmic\npatch = x_min\n[sweep]\ncontact g = 0 : 0.5 : 1\n";
    const std::string patch = "[mesh]\ntype = tensor\nx = linspace(0, 1um, 101) + linspace(1um, 1.01um, 11)\n" + body +
                              "[contact g]\nkind = electrode\npatch = x_max\nelectrode = metal\nwork_function = 4.6\n" + tail;
    const std::string metal = "[mesh]\ntype = tensor\nx = linspace(0, 1um, 101) + linspace(1um, 1.01um, 11) + "
                              "linspace(1.01um, 1.02um, 6)\n" + body +
                              "[region gm]\nmaterial = metal\nbox = 1.01um 1.02um\nwork_function = 4.6\n"
                              "[contact g]\nkind = electrode\nregion = gm\n" + tail;
    const deck::Plan pp = must(patch);
    const deck::Plan pm = must(metal);
    REQUIRE(pm.device.contacts().front().nodes.size() == 6);
    const auto a = deck::run(pp);
    const auto b = deck::run(pm);
    if (!a) FAIL(a.error().message);
    if (!b) FAIL(b.error().message);
    REQUIRE(a->sweep->points.size() == b->sweep->points.size());
    const std::size_t shared = pp.device.mesh().node_count();  // the patch mesh's nodes come first in x
    for (std::size_t k = 0; k < a->sweep->points.size(); ++k) {
        const auto& va = a->sweep->points[k].fields.potential_V;
        const auto& vb = b->sweep->points[k].fields.potential_V;
        double worst = 0.0;
        for (std::size_t i = 0; i < shared; ++i) worst = std::max(worst, std::abs(va[i] - vb[i]));
        INFO("point " << k);
        CHECK(worst <= 1e-9);
    }
    // A preset is its work function; the refusals.
    CHECK(must(changed(metal, "work_function = 4.6\n[contact", "metal = mo\n[contact")).device.contacts().front().electrode.work_function_eV == 4.6);
    CHECK(failure(changed(metal, "[contact g]\nkind = electrode\nregion = gm\n", "")).find("held by no electrode contact") != std::string::npos);
    CHECK(failure(changed(metal, "region = gm", "region = ox")).find("is not a metal region") != std::string::npos);
    CHECK(failure(changed(metal, "region = gm", "region = gm\npatch = x_max")).find("takes no 'patch'") != std::string::npos);
    CHECK(failure(changed(metal, "region = gm", "region = gm\nwork_function = 4")).find("gives the electrode's work function") != std::string::npos);
    CHECK(failure(changed(metal, "work_function = 4.6\n[contact", "[contact")).find("needs 'work_function' or 'metal") != std::string::npos);
    CHECK(failure(changed(metal, "work_function = 4.6\n[contact", "work_function = 4.6\neps_r = 2\n[contact")).find("takes 'work_function' or 'metal', not 'eps_r'") != std::string::npos);
    CHECK(failure(changed(metal, "material = silicon\n", "material = silicon\nwork_function = 4\n")).find("is for a metal region") != std::string::npos);
    CHECK(failure(changed(metal, "box = 1um 1.01um", "box = 0.5um 0.51um")).find("touches a semiconductor") != std::string::npos);
    CHECK(failure(metal + "[material m2]\nbase = metal\n").find("not a [material] base") != std::string::npos);
}
