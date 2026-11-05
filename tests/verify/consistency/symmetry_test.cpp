// C11, symmetry and invariance (DECISIONS Unit 55 section 11.6 and the C11 ruling of 11.13): a 2D MOS structure
// symmetric under y -> W - y (a meshed gate over the whole oxide face, source and drain pockets that swap, the
// substrate), solved with the drain at 0.1 V and, mirrored, the source at 0.1 V: the currents and the gate charge
// agree to the ruled bounds; a generated unstructured pn diode rotated by 30 degrees, rebuilt by make_dual_mesh,
// gives the same terminal currents (the ruled translation is refused by P2 on this mesh: see Rotated below). A
// broken symmetry (doping, a bias, the contact map) is NOT APPLICABLE with the item named; a corrupted current or
// charge fails; a missing Newton record leaves the charges NOT EVALUATED.
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace c = NiTCAD::verify::consistency;
namespace v = NiTCAD::verify;
namespace solve = NiTCAD::solve;
namespace device = NiTCAD::device;
namespace mesh = NiTCAD::mesh;
namespace physics = NiTCAD::physics;
namespace results = NiTCAD::results;
using v::Outcome;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> x(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) x[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return x;
}

results::BiasPoint solved(const device::Device& d, std::vector<double> bias) {
    solve::BiasOptions o;
    o.step_control.max_cuts = 20;
    auto p = solve::solve_bias(d, bias, o);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

constexpr double W = 3e-6;  // the MOS structure's width (y)

// Oxide on [-5 nm, 0), p-silicon (1e17) to 120 nm with n+ pockets (1e19) under the source (y < 0.75 um) and the
// drain (y > 2.25 um) to 30 nm; the gate electrode on all of x_min; source and drain contacts on the pockets' y_min
// and y_max faces; the substrate on x_max; a fixed interface charge. Contacts: gate, substrate, source, drain.
device::Device symmetric_mos() {
    const auto x = *mesh::straddle_interface(uniform(-5e-7, 1.2e-5, 30), 0.0, 1e-7);
    const auto t = uniform(0.0, W, 7);
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
        if ((p[1] < 0.75e-6 || p[1] > 2.25e-6) && p[0] < 3e-6) donors[i] = 1e19;
    }
    std::vector<mesh::NodeId> source, drain;
    for (const auto node : m.find_boundary("y_min")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(node)][0];
        if (px > 0.0 && px < 3e-6) source.push_back(node);
    }
    for (const auto node : m.find_boundary("y_max")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(node)][0];
        if (px > 0.0 && px < 3e-6) drain.push_back(node);
    }
    auto gate_nodes = m.find_boundary("x_min")->nodes;
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = 3e11;
    device::Contact gate{"gate", device::ContactKind::electrode, std::move(gate_nodes)};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)},
                      {"source", device::ContactKind::ohmic, std::move(source)},
                      {"drain", device::ContactKind::ohmic, std::move(drain)}},
         .interfaces = {std::move(f)}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

// The node of (x, W - y) for each node (x, y).
std::vector<mesh::NodeId> mirror_map(const mesh::Mesh& m) {
    const auto pts = m.points();
    std::vector<mesh::NodeId> map(pts.size(), -1);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        for (std::size_t j = 0; j < pts.size(); ++j) {
            if (pts[j][0] == pts[i][0] && std::abs(pts[j][1] - (W - pts[i][1])) <= 1e-9 * W) {
                map[i] = static_cast<mesh::NodeId>(j);
            }
        }
    }
    return map;
}

struct Mirror {
    device::Device device;
    results::BiasPoint a, b;
    c::DeviceData data;
    c::SymmetryMap map;
};

const Mirror& mirror() {
    static const std::optional<Mirror> s = []() -> std::optional<Mirror> {
        device::Device d = symmetric_mos();
        auto a = solved(d, {1.0, 0.0, 0.0, 0.1});
        auto b = solved(d, {1.0, 0.0, 0.1, 0.0});
        Mirror out{std::move(d), std::move(a), std::move(b), {}, {}};
        out.data = NiTCAD::verify::benchmark::device_data(out.device);
        out.map.node = mirror_map(out.device.mesh());
        out.map.contact = {0, 1, 3, 2};
        return out;
    }();
    REQUIRE(s.has_value());
    return *s;
}

// A generated unstructured pn diode on [0, L] x [0, H] (acceptors for x < L/2, donors above), ohmic contacts on
// x = 0 and x = L; and the same simplices rotated by 30 degrees about the origin, rebuilt by make_dual_mesh,
// solved forward at 0.4 V. (A translation, alone or with the rotation, is refused by P2 on this mesh: the corner
// cell at (L, H) has four exactly cocircular nodes and the rounding of the moved coordinates breaks the tie
// inward; DECISIONS 11.14, C11.)
struct Rotated {
    std::optional<device::Device> a, b;
    results::BiasPoint pa, pb;
    c::DeviceData da, db;
    c::SymmetryMap map;
    std::string refused;  // make_dual_mesh's or Device::create's refusal of the moved mesh, if any
};

constexpr double L = 2e-4, H = 1e-4;

std::optional<device::Device> diode_on(mesh::Mesh m, const std::vector<double>& donors,
                                       const std::vector<double>& acceptors, std::string& error) {
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("left");
    auto right = m.find_boundary("right");
    if (left == nullptr || right == nullptr) {
        error = "a contact patch is missing";
        return std::nullopt;
    }
    auto ln = left->nodes;
    auto rn = right->nodes;
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = donors,
                                     .acceptors = acceptors,
                                     .contacts = {{"left", device::ContactKind::ohmic, std::move(ln)},
                                                  {"right", device::ContactKind::ohmic, std::move(rn)}}});
    if (!d) {
        error = d.error().message;
        return std::nullopt;
    }
    return std::move(*d);
}

Rotated make_rotated() {
    Rotated out;
    mesh::Geometry2D g;
    g.polygons = {{.vertices = {{0.0, 0.0, 0.0}, {L, 0.0, 0.0}, {L, H, 0.0}, {0.0, H, 0.0}}, .region = 0}};
    g.patches = {{"left", {0.0, 0.0, 0.0}, {0.0, H, 0.0}}, {"right", {L, 0.0, 0.0}, {L, H, 0.0}}};
    auto gen = mesh::generate_mesh(g, {.max_size = L / 12.0});
    if (!gen) {
        out.refused = "generate_mesh: " + gen.error().message;
        return out;
    }
    const std::size_t n = gen->mesh.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (gen->mesh.points()[i][0] < 0.5 * L ? acceptors : donors)[i] = 1e17;
    mesh::SimplexMesh moved = gen->simplices;
    const double angle = std::numbers::pi / 6.0, cs = std::cos(angle), sn = std::sin(angle);
    for (auto& p : moved.points) {
        const double x = p[0], y = p[1];
        p = {cs * x - sn * y, sn * x + cs * y, 0.0};
    }
    auto dual = mesh::make_dual_mesh(moved);
    if (!dual) {
        out.refused = "make_dual_mesh: " + dual.error().message;
        return out;
    }
    std::string error;
    out.a = diode_on(gen->mesh, donors, acceptors, error);
    out.b = diode_on(std::move(*dual), donors, acceptors, error);
    if (!out.a || !out.b) {
        out.refused = "Device::create: " + error;
        return out;
    }
    out.pa = solved(*out.a, {0.4, 0.0});  // forward: the p side (x < L/2) on contact left
    out.pb = solved(*out.b, {0.4, 0.0});
    out.da = NiTCAD::verify::benchmark::device_data(*out.a);
    out.db = NiTCAD::verify::benchmark::device_data(*out.b);
    for (std::size_t i = 0; i < n; ++i) out.map.node.push_back(static_cast<mesh::NodeId>(i));
    out.map.contact = {0, 1};
    return out;
}

const Rotated& rotated() {
    static const Rotated s = make_rotated();
    return s;
}

c::CheckResult check(const Mirror& s, const results::BiasPoint& b, const c::DeviceData& db,
                     const c::SymmetryMap& map) {
    return c::check_c11(s.a, s.device.mesh(), s.data, b, s.device.mesh(), db, map);
}

}  // namespace

TEST_CASE("C11 passes on a mirrored MOS structure", "[verify][consistency]") {
    const Mirror& s = mirror();
    for (const auto j : s.map.node) REQUIRE(j >= 0);
    const auto r = check(s, s.b, s.data, s.map);
    CAPTURE(r.comparison.reason);
    for (const auto& p : r.parts) {
        std::printf("C11 mirror: %s residual %.3e tolerance %.3e at %s\n", p.name.c_str(), p.residual, p.tolerance,
                    p.where.c_str());
    }
    REQUIRE(r.id == "C11");
    REQUIRE(r.comparison.outcome == Outcome::pass);
    REQUIRE(r.parts.size() == 2);  // currents and the gate charge
    // Not vacuous: the drain current is far above its bound, and the two drains carry opposite-signed currents
    // in the unmapped order.
    REQUIRE(std::abs(s.a.terminal_current[3]) > 1e3 * r.parts.front().tolerance);

    // The gate charge's bound restated from the ruling: 64 eps max |Q| + (sum over the gate's edges of C_f 2 V_T)
    // (u_A + u_B), C_f = eps_0 eps_h A_f / d_f (the meshed electrode has no lumped term).
    const auto& gate = s.data.contacts[0];
    REQUIRE(gate.kind == c::ContactKind::electrode);
    const double vt = NiTCAD::base::k_B_eV_per_K * s.data.temperature_K;
    std::vector<std::uint8_t> member(s.device.mesh().node_count(), 0);
    for (const auto v : gate.nodes) member[static_cast<std::size_t>(v)] = 1;
    double term = 0.0;
    for (const auto& e : s.device.mesh().edges()) {
        const auto i = static_cast<std::size_t>(e.first), j = static_cast<std::size_t>(e.second);
        if (member[i] == 0 && member[j] == 0) continue;
        const double ei = s.data.relative_permittivity[i], ej = s.data.relative_permittivity[j];
        term += NiTCAD::base::eps0_F_per_cm * (2.0 * ei * ej / (ei + ej)) * e.coupling_area / e.length * 2.0 * vt;
    }
    const double u = s.a.convergence.iterations.back().update + s.b.convergence.iterations.back().update;
    const double q = std::max(std::abs(s.a.gate_charge[0]), std::abs(s.b.gate_charge[0]));
    const double bound = 64.0 * DBL_EPSILON * q + term * u;
    CAPTURE(bound, term * u, q);
    REQUIRE(term * u > 0.0);
    REQUIRE(std::abs(r.parts[1].tolerance - bound) <= 1e-12 * bound);
}

TEST_CASE("C11 passes on a rotated unstructured mesh", "[verify][consistency]") {
    const Rotated& s = rotated();
    CAPTURE(s.refused);
    REQUIRE(s.refused.empty());
    const auto r = c::check_c11(s.pa, s.a->mesh(), s.da, s.pb, s.b->mesh(), s.db, s.map);
    CAPTURE(r.comparison.reason);
    for (const auto& p : r.parts) {
        std::printf("C11 rotation: %s residual %.3e tolerance %.3e at %s (I = %.6e)\n", p.name.c_str(), p.residual,
                    p.tolerance, p.where.c_str(), s.pa.terminal_current[0]);
    }
    REQUIRE(r.comparison.outcome == Outcome::pass);
    REQUIRE(std::abs(s.pa.terminal_current[0]) > 1e3 * r.parts.front().tolerance);
}

TEST_CASE("C11 is NOT APPLICABLE when the map is not a symmetry", "[verify][consistency]") {
    const Mirror& s = mirror();
    const auto names = [&](const c::CheckResult& r, const char* item) {
        CAPTURE(r.comparison.reason, item);
        REQUIRE(r.comparison.outcome == Outcome::not_applicable);
        REQUIRE(r.comparison.reason.find("the map is a symmetry of the device") != std::string::npos);
        REQUIRE(r.comparison.reason.find(item) != std::string::npos);
    };
    auto doping = s.data;  // one node's donors changed in B's data
    doping.donors[s.map.node[40]] *= 1.0 + 1e-12;
    names(check(s, s.b, doping, s.map), "doping differs");

    auto bias = s.b;  // B's drain bias not A's source bias
    bias.bias_V[2] = 0.1 + 1e-12;
    names(check(s, bias, s.data, s.map), "bias differs");

    auto contacts = s.map;  // the contact map not a bijection
    contacts.contact = {0, 1, 2, 2};
    names(check(s, s.b, s.data, contacts), "contact map is not a bijection");

    auto identity = s.map;  // the identity node map: the contacts' nodes do not correspond under {0, 1, 3, 2}
    for (std::size_t i = 0; i < identity.node.size(); ++i) identity.node[i] = static_cast<mesh::NodeId>(i);
    names(check(s, s.b, s.data, identity), "contact nodes do not correspond");

    auto nodes = s.map;  // two nodes mapped to one
    nodes.node[1] = nodes.node[0];
    names(check(s, s.b, s.data, nodes), "node map is not a bijection");

    auto eps = s.data;  // one node's permittivity
    eps.relative_permittivity[s.map.node[3]] *= 1.0 + 1e-12;
    names(check(s, s.b, eps, s.map), "permittivity differs");
}

TEST_CASE("C11 is NOT APPLICABLE when the meshes or contacts do not correspond", "[verify][consistency]") {
    const Mirror& s = mirror();
    const mesh::Mesh& m = s.device.mesh();
    const auto names = [&](const c::CheckResult& r, const char* item) {
        CAPTURE(r.comparison.reason, item);
        REQUIRE(r.comparison.outcome == Outcome::not_applicable);
        REQUIRE(r.comparison.reason.find(item) != std::string::npos);
    };
    // B's mesh: A's with one entry changed by more than 64 eps relative (and within from_parts' 1e-12 on lengths).
    const auto rebuilt = [&](auto change) {
        std::vector<mesh::Point> points(m.points().begin(), m.points().end());
        std::vector<double> volumes(m.volumes().begin(), m.volumes().end());
        std::vector<mesh::Edge> edges(m.edges().begin(), m.edges().end());
        std::vector<mesh::BoundaryPatch> boundary(m.boundary().begin(), m.boundary().end());
        change(volumes, edges);
        auto out = mesh::Mesh::from_parts(m.dimension(), std::move(points), std::move(volumes), std::move(edges),
                                          std::move(boundary), m.origin());
        REQUIRE(out.has_value());
        return std::move(*out);
    };
    const auto with_mesh = [&](const mesh::Mesh& mb) {
        return c::check_c11(s.a, m, s.data, s.b, mb, s.data, s.map);
    };
    const std::size_t e = m.edges().size() / 2;
    const auto node = static_cast<std::size_t>(s.map.node[40]);
    names(with_mesh(rebuilt([&](auto& v, auto&) { v[node] *= 1.0 + 1e-12; })), "volume differs");
    names(with_mesh(rebuilt([&](auto&, auto& ed) { ed[e].length *= 1.0 + 1e-13; })), "length differs");
    names(with_mesh(rebuilt([&](auto&, auto& ed) { ed[e].coupling_area *= 1.0 + 1e-12; })), "coupling differs");
    // The unchanged mesh rebuilt: still a symmetry (the comparisons are not identity tests).
    REQUIRE(with_mesh(rebuilt([](auto&, auto&) {})).comparison.outcome == Outcome::pass);

    auto region = s.data;  // one node's region in B's data
    region.region[node] = 7;
    names(check(s, s.b, region, s.map), "region differs");
    auto insulator = s.data;  // one node's insulator flag in B's data
    insulator.insulator[node] = static_cast<std::uint8_t>(1 - insulator.insulator[node]);
    names(check(s, s.b, insulator, s.map), "region differs");
    auto kind = s.data;  // the substrate's kind in B's data
    kind.contacts[1].kind = c::ContactKind::schottky;
    names(check(s, s.b, kind, s.map), "contact kind differs (substrate)");
}

TEST_CASE("C11 fails on a corrupted current or charge", "[verify][consistency]") {
    const Mirror& s = mirror();
    const auto base = check(s, s.b, s.data, s.map);
    const double i_tol = base.parts[0].tolerance;  // the worst current's bound
    const double res = s.a.terminal_current_resolution[3] + s.b.terminal_current_resolution[2];
    auto current = s.b;  // B's source current (A's drain image) moved by twice its pair's bound
    current.terminal_current[2] += 2.0 * (res + 64.0 * DBL_EPSILON * std::abs(current.terminal_current[2]) * 1.01);
    const auto rc = check(s, current, s.data, s.map);
    CAPTURE(i_tol, res, rc.parts[0].where);
    REQUIRE(rc.comparison.outcome == Outcome::fail);
    REQUIRE(rc.parts[0].where == "contact drain -> source");

    auto charge = s.b;  // B's gate charge moved by twice the charge part's bound
    charge.gate_charge[0] += 2.0 * base.parts[1].tolerance;
    const auto rq = check(s, charge, s.data, s.map);
    REQUIRE(rq.parts[1].outcome == Outcome::fail);

    auto no_newton = s.b;  // no Newton record: the charges NOT EVALUATED, the check with them
    no_newton.convergence.iterations.clear();
    const auto rn = check(s, no_newton, s.data, s.map);
    REQUIRE(rn.parts[1].outcome == Outcome::not_evaluated);
    REQUIRE(rn.comparison.outcome == Outcome::not_evaluated);
}
