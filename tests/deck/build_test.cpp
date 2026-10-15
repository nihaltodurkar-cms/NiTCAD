// The deck's meaning (Unit 31, deck/build.hpp): decks give the devices, options and points built by
// hand (equal run identities, so every input the solver reads is equal), and every mistake is an
// error naming its line.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
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

[solver]
equations = drift_diffusion
max_iterations = 40
tol_update = 1e-9
max_update = 2
line_search = natural
line_search.max_backtracks = 12
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
    const auto& o = p.options;
    CHECK(o.newton.max_iterations == 40);
    CHECK(o.newton.tol_update == 1e-9);
    CHECK(o.newton.max_update == 2.0);
    CHECK(o.newton.line_search.enabled);
    CHECK(o.newton.line_search.test == solve::LineSearchTest::natural);
    CHECK(o.newton.line_search.max_backtracks == 12);
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
    expect(failure(head), "the deck needs exactly one of [equilibrium] and [sweep]");
    expect(failure(head + "[equilibrium]\n[sweep]\n"), "the deck needs exactly one of [equilibrium] and [sweep]");
    expect(failure(head + "[equilibrium]\n[mesh]\n"), "t.deck:13: a second [mesh] section");
    expect(failure(head + "[equilibrium]\n[models]\nsrh = maybe\n"), "t.deck:14: 'maybe' is not on or off");
    expect(failure(head + "[equilibrium]\n[models]\nquantum = on\n"), "t.deck:14: unknown key 'quantum' in [models]");
    expect(failure(head + "[equilibrium]\n[solver]\nbackend = pardiso\n"), "t.deck:13: backend = pardiso needs 'mkl'");
    expect(failure(head + "[equilibrium]\n[solver]\nmax_cuts = 99\n"), "t.deck:14: 'max_cuts' must lie in 0..52");
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
