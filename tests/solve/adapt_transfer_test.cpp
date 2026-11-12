// Gates for state transfer (Unit 17, slice 3; DECISIONS.md Unit 17 AD5; AC3): linear potentials and
// quasi-Fermi levels transfer exactly on tensor and generated meshes; a diode's and a MOSFET's state
// transferred to a refined mesh converge to the solution from scratch in fewer Newton iterations;
// the point locator reproduces linear fields.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/locate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/solve/adapt_transfer.hpp"
#include "NiTCAD/solve/bias.hpp"

using namespace NiTCAD;

namespace {

deck::Plan plan_of(const std::string& text) {
    auto d = deck::parse(text, "t.deck", {});
    if (!d) FAIL(d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}



// A synthetic state: psi, E_Fn, E_Fp linear in (x, y) (um), densities consistent with them.
struct Linear {
    double c_n = std::log(1e10), c_p = std::log(1e10);
    double psi(const mesh::Point& p) const { return 0.1 * p[0] / 1e-4 + 0.05 * p[1] / 1e-4; }
    double fn(const mesh::Point& p) const { return 0.02 * p[0] / 1e-4 - 0.03 * p[1] / 1e-4; }
    double fp(const mesh::Point& p) const { return -0.04 * p[0] / 1e-4; }
};

void synthetic(const device::Device& d, results::NodeFields& s, results::BandDiagram& b) {
    const Linear L;
    const double vt = base::thermal_voltage(300.0);
    const auto pts = d.mesh().points();
    const std::size_t n = pts.size();
    s.potential_V.resize(n);
    s.n_cm3.resize(n);
    s.p_cm3.resize(n);
    b.electron_fermi_eV.resize(n);
    b.hole_fermi_eV.resize(n);
    b.conduction_eV.assign(n, 0.0);
    b.valence_eV.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        s.potential_V[i] = L.psi(pts[i]);
        b.electron_fermi_eV[i] = L.fn(pts[i]);
        b.hole_fermi_eV[i] = L.fp(pts[i]);
        s.n_cm3[i] = std::exp(L.c_n + (L.psi(pts[i]) + L.fn(pts[i])) / vt);
        s.p_cm3[i] = std::exp(L.c_p - (L.psi(pts[i]) + L.fp(pts[i])) / vt);
    }
}

void check_linear(const device::Device& to, const results::NodeFields& t) {
    const Linear L;
    const double vt = base::thermal_voltage(300.0);
    const auto pts = to.mesh().points();
    for (std::size_t i = 0; i < pts.size(); ++i) {
        CHECK(std::abs(t.potential_V[i] - L.psi(pts[i])) <= 1e-12);
        const double fn = vt * (std::log(t.n_cm3[i]) - L.c_n) - t.potential_V[i];
        const double fp = -vt * (std::log(t.p_cm3[i]) - L.c_p) - t.potential_V[i];
        CHECK(std::abs(fn - L.fn(pts[i])) <= 1e-12);
        CHECK(std::abs(fp - L.fp(pts[i])) <= 1e-12);
    }
}

}  // namespace

TEST_CASE("transfer: linear fields are reproduced by the locator on tensor and simplex meshes", "[adapt_transfer]") {
    // Tensor (multilinear: exact for linear).
    const auto m = mesh::make_tensor_grid(std::vector<double>{0, 1, 3}, std::vector<double>{0, 2, 3});
    REQUIRE(m.has_value());
    const std::vector<double> x{0, 1, 3}, y{0, 2, 3};
    auto tc = mesh::TensorCells::create(*m, x, y);
    REQUIRE(tc.has_value());
    const auto L = mesh::PointLocator::create(*tc);
    CHECK(L.node_count() == 9);
    const auto f = [](const mesh::Point& p) { return 2.0 + 3.0 * p[0] - 1.5 * p[1]; };
    for (const mesh::Point q : {mesh::Point{0.5, 0.5, 0}, mesh::Point{2.9, 2.1, 0}, mesh::Point{3, 3, 0}, mesh::Point{0, 0, 0}}) {
        const auto it = L.locate(q);
        double v = 0.0;
        for (int k = 0; k < it.count; ++k) v += it.weights[static_cast<std::size_t>(k)] * f(m->points()[static_cast<std::size_t>(it.nodes[static_cast<std::size_t>(k)])]);
        CHECK(std::abs(v - f(q)) <= 1e-12 * 10);
        CHECK(it.inside);
    }
    CHECK(!L.locate({-1, 0, 0}).inside);
    // Simplex 2D: a generated mesh.
    mesh::Geometry2D g;
    g.polygons = {{{{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}}, 0}};
    mesh::MeshOptions2D o;
    o.max_size = 0.3;
    auto gm = mesh::generate_mesh(g, o);
    REQUIRE(gm.has_value());
    auto L2 = mesh::PointLocator::create(gm->simplices);
    REQUIRE(L2.has_value());
    for (int i = 0; i <= 10; ++i) {
        for (int j = 0; j <= 10; ++j) {
            const mesh::Point q{0.2 * i, 0.1 * j, 0};
            const auto it = L2->locate(q);
            double v = 0.0, w = 0.0;
            for (int k = 0; k < it.count; ++k) {
                v += it.weights[static_cast<std::size_t>(k)] * f(gm->simplices.points[static_cast<std::size_t>(it.nodes[static_cast<std::size_t>(k)])]);
                w += it.weights[static_cast<std::size_t>(k)];
            }
            CHECK(std::abs(w - 1.0) <= 1e-12);
            CHECK(std::abs(v - f(q)) <= 1e-12 * 10);
        }
    }
    // Simplex 3D: a swept mesh's tetrahedra.
    mesh::ExtrudedGeometry xg;
    xg.axis = 2;
    xg.stations = {0.0, 2.0};
    xg.cells = {{"r", {{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}}}};
    xg.regions = {{0}};
    mesh::SweptOptions so;
    so.section.max_size = 0.4;
    so.lines = {0.0, 0.7, 2.0};
    auto sw = mesh::generate_swept(xg, so);
    REQUIRE(sw.has_value());
    REQUIRE(sw->section.has_value());
    CHECK(sw->section->points.size() * 3 == sw->simplices.points.size());
    auto L3 = mesh::PointLocator::create(sw->simplices);
    REQUIRE(L3.has_value());
    const auto f3 = [](const mesh::Point& p) { return 1.0 + p[0] - 2.0 * p[1] + 0.5 * p[2]; };
    for (const mesh::Point q : {mesh::Point{0.3, 0.3, 0.3}, mesh::Point{1.9, 0.9, 1.9}, mesh::Point{1.0, 0.5, 0.7}}) {
        const auto it = L3->locate(q);
        CHECK(it.inside);
        double v = 0.0;
        for (int k = 0; k < it.count; ++k) v += it.weights[static_cast<std::size_t>(k)] * f3(sw->simplices.points[static_cast<std::size_t>(it.nodes[static_cast<std::size_t>(k)])]);
        CHECK(std::abs(v - f3(q)) <= 1e-12 * 10);
    }
}

TEST_CASE("transfer: a linear potential and linear quasi-Fermi levels transfer exactly (AC3)", "[adapt_transfer]") {
    // Tensor to tensor, not nested.
    const auto tensor = [](int nx, int ny) {
        return "type = tensor\nx = linspace(0, 1um, " + std::to_string(nx) + ")\ny = linspace(0, 1um, " +
               std::to_string(ny) + ")\n";
    };
    const std::string body = R"(
[region silicon]
material = silicon
[doping]
profile = donor constant 1e16
[contact a]
kind = ohmic
patch = x_min
[equilibrium]
)";
    const deck::Plan a = plan_of("[mesh]\n" + tensor(5, 4) + body);
    const deck::Plan b = plan_of("[mesh]\n" + tensor(9, 7) + body);
    results::NodeFields s;
    results::BandDiagram bd;
    synthetic(a.device, s, bd);
    const auto loc = mesh::PointLocator::create(*a.device.cells());
    solve::TransferReport rep;
    auto t = solve::transfer_state(a.device, loc, s, bd, b.device, &rep);
    REQUIRE(t.has_value());
    CHECK(rep.nodes == b.device.mesh().node_count());
    CHECK(rep.outside == 0);
    CHECK(rep.other_region == 0);
    check_linear(b.device, *t);

    // Generated 2D to generated 2D (a coarser and a finer h).
    const auto gen = [](const std::string& h) {
        return "type = generate2d\nmax_size = " + h + "\npatch left = 0 0, 0 1um\n";
    };
    const std::string gbody = R"(
[region silicon]
material = silicon
polygon = 0 0, 1um 0, 1um 1um, 0 1um
[doping]
profile = donor constant 1e16
[contact a]
kind = ohmic
patch = left
[equilibrium]
)";
    const deck::Plan ga = plan_of("[mesh]\n" + gen("0.3um") + gbody);
    const deck::Plan gb = plan_of("[mesh]\n" + gen("0.11um") + gbody);
    synthetic(ga.device, s, bd);
    REQUIRE(ga.generated.has_value());
    auto gl = mesh::PointLocator::create(ga.generated->simplices);
    REQUIRE(gl.has_value());
    auto gt = solve::transfer_state(ga.device, *gl, s, bd, gb.device, &rep);
    REQUIRE(gt.has_value());
    check_linear(gb.device, *gt);

    // Refusals.
    results::NodeFields bad = s;
    bad.p_cm3.pop_back();
    CHECK(!solve::transfer_state(ga.device, *gl, bad, bd, gb.device).has_value());
    CHECK(!solve::transfer_state(ga.device, loc, s, bd, gb.device).has_value());  // locator of another mesh
}

namespace {

std::string diode2d(int nx) {
    return R"(
[device]
temperature = 300
[mesh]
type = tensor
x = linspace(0, 2um, )" + std::to_string(nx) + R"()
y = 0, 0.1um
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e17 axis=y x=-1um:1um
profile = donor constant 1e17 axis=y x=1um:3um
[contact anode]
kind = ohmic
patch = x_min
[contact cathode]
kind = ohmic
patch = x_max
[sweep]
contact anode = 0.5
)";
}

// A lateral n+ p n+ MOSFET on a tensor mesh: source and drain full-depth columns, the gate on top.
std::string mosfet(int nx, int ny) {
    return R"(
[device]
temperature = 300
[mesh]
type = tensor
x = linspace(0, 1um, )" + std::to_string(nx) + R"()
y = linspace(0, 0.3um, )" + std::to_string(ny) + R"()
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e17
profile = donor constant 1e20 axis=y x=-1um:0.25um
profile = donor constant 1e20 axis=y x=0.75um:2um
[contact source]
kind = ohmic
patch = x_min
[contact drain]
kind = ohmic
patch = x_max
[contact body]
kind = ohmic
patch = y_min
box = 0.35um 0.65um 0 0
[contact gate]
kind = gate
patch = y_max
box = 0.25um 0.75um 0.3um 0.3um
oxide_thickness = 5nm
electrode = n_poly
[sweep]
contact gate = 1.0
contact drain = 0.5
[solver]
max_cuts = 20
)";
}

int iterations(const results::Sweep& s) {
    int n = 0;
    for (const auto& p : s.points) n += static_cast<int>(p.convergence.iterations.size());
    return n;
}

void transfer_converges(const std::string& coarse_text, const std::string& fine_text) {
    deck::Plan coarse = plan_of(coarse_text);
    deck::Plan fine = plan_of(fine_text);
    auto c = deck::run(coarse);
    REQUIRE(c.has_value());
    REQUIRE(c->sweep.has_value());
    if (c->sweep->stopped) FAIL("coarse sweep stopped: " + c->sweep->stopped->message);
    // From scratch on the fine mesh (its own starting equilibrium and bias steps).
    auto scratch = solve::sweep_bias(fine.device, fine.points, fine.options);
    REQUIRE(scratch.has_value());
    REQUIRE(!scratch->stopped);
    // From the transferred state, at the same bias.
    const results::BiasPoint& p = c->sweep->points.back();
    const auto loc = mesh::PointLocator::create(*coarse.device.cells());
    auto init = solve::transfer_state(coarse.device, loc, p.fields, p.bands, fine.device);
    REQUIRE(init.has_value());
    solve::BiasOptions o = fine.options;
    o.step_control.initial_bias_V = p.bias_V;
    auto warm = solve::sweep_bias(fine.device, fine.points, o, &*init);
    REQUIRE(warm.has_value());
    REQUIRE(!warm->stopped);
    const auto& a = scratch->points.back().terminal_current;
    const auto& b = warm->points.back().terminal_current;
    REQUIRE(a.size() == b.size());
    double scale = 0.0;
    for (const double x : a) scale = std::max(scale, std::abs(x));
    for (std::size_t k = 0; k < a.size(); ++k) CHECK(std::abs(a[k] - b[k]) <= 1e-9 * scale);
    const int from_scratch = iterations(*scratch);
    const int transferred = iterations(*warm);
    UNSCOPED_INFO("Newton iterations: from scratch " << from_scratch << " (sweep points only), transferred " << transferred);
    CHECK(transferred < from_scratch);
    WARN("AC3 Newton iterations: from scratch " << from_scratch << ", transferred " << transferred);
}

}  // namespace

TEST_CASE("transfer: a diode's state on a refined mesh converges to the solution from scratch, faster (AC3)",
          "[adapt_transfer]") {
    transfer_converges(diode2d(41), diode2d(97));
}

TEST_CASE("transfer: a MOSFET's state on a refined mesh converges to the solution from scratch, faster (AC3)",
          "[adapt_transfer]") {
    transfer_converges(mosfet(41, 16), mosfet(73, 29));
}
