// Device and solver gates on generated 2D meshes (Unit 16; DECISIONS.md "Gates"). Thresholds are
// those predefined in DECISIONS.md; nothing here was tuned after a run.
// - Linear-field consistency (ctest): Laplace in an oxide between two electrodes is solved
//   exactly by the box method on a Voronoi dual, whatever the triangles.
// - [.unstructured], Release only (CMakeLists.txt), minutes each:
//   - convergence of a smooth-doped 2D diode against a 4x finer tensor reference: discrete H1
//     order >= 0.9, L2 >= 1.5 (equilibrium potential), terminal current order >= 0.9;
//   - an abrupt, laterally uniform diode against the 1D reference: current error no larger than a
//     2D tensor grid's with the same node count;
//   - the MOSFET: |Q_u - Q_ref| <= 1.5 |Q_t - Q_ref| + |Q_ref - Q_t| / 3 for I_D (resolved points),
//     V_T, maximum g_m and the subthreshold swing; every point converged.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "legacy_graded_mesh.hpp"
#include "legacy_moscap.hpp"

using namespace NiTCAD;

namespace {

mesh::RegionPolygon rect(double x0, double y0, double x1, double y1, std::int32_t region) {
    return {{{x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}}, region};
}

std::vector<mesh::NodeId> patch_nodes(const mesh::Mesh& m, const char* name) {
    const mesh::BoundaryPatch* p = m.find_boundary(name);
    REQUIRE(p != nullptr);
    return p->nodes;
}

std::vector<double> uniform(double length, std::size_t cells) {
    std::vector<double> v(cells + 1);
    for (std::size_t i = 0; i <= cells; ++i) v[i] = length * static_cast<double>(i) / static_cast<double>(cells);
    return v;
}

// Bilinear interpolation of a node field on the uniform tensor grid (x fastest).
double bilinear(const std::vector<double>& f, const std::vector<double>& x, const std::vector<double>& y, double px,
                double py) {
    const double hx = x[1] - x[0], hy = y[1] - y[0];
    const auto i = std::min(static_cast<std::size_t>(std::max(0.0, std::floor(px / hx))), x.size() - 2);
    const auto j = std::min(static_cast<std::size_t>(std::max(0.0, std::floor(py / hy))), y.size() - 2);
    const double s = (px - x[i]) / hx, t = (py - y[j]) / hy;
    const auto at = [&](std::size_t a, std::size_t b) { return f[a + x.size() * b]; };
    return (1 - s) * (1 - t) * at(i, j) + s * (1 - t) * at(i + 1, j) + (1 - s) * t * at(i, j + 1) + s * t * at(i + 1, j + 1);
}

}  // namespace

TEST_CASE("unstructured: Laplace in an oxide is exact on a Voronoi dual", "[solve][unstructured_fast]") {
    // An oxide between two electrodes (1 V apart) and, separately, a silicon block with an ohmic
    // contact (a device needs a semiconductor). The oxide carries no charge, so the potential is
    // linear in x; the box method reproduces it exactly when every dual face is orthogonal to its
    // edge, which the circumcentric dual of a Delaunay mesh guarantees for skew triangles too.
    const double L = 1e-4, H = 5e-5;
    mesh::Geometry2D g{{rect(0, 0, L, H, 0), rect(2 * L, 0, 3 * L, H, 1)},
                       {{"e1", {0, 0, 0}, {0, H, 0}}, {"e2", {L, 0, 0}, {L, H, 0}}, {"ground", {2 * L, 0, 0}, {2 * L, H, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = 1.3e-5;  // not commensurate with L: genuinely unstructured
    const auto gm = mesh::generate_mesh(g, o);
    REQUIRE(gm.has_value());
    REQUIRE(gm->certificate.has_value());  // two regions on a general mesh
    const std::vector<device::Region> regions{{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}};
    const std::vector<device::DopingProfile> profiles{{.dopant = device::Dopant::donor, .peak_cm3 = 1e16}};
    const auto doping = device::sample(profiles, gm->mesh, gm->node_region, regions);
    REQUIRE(doping.has_value());
    device::DeviceDescription d{.mesh = gm->mesh,
                                .temperature_K = 300.0,
                                .regions = regions,
                                .node_region = gm->node_region,
                                .donors = doping->donors,
                                .acceptors = doping->acceptors,
                                .contacts = {{"e1", device::ContactKind::electrode, patch_nodes(gm->mesh, "e1")},
                                             {"e2", device::ContactKind::electrode, patch_nodes(gm->mesh, "e2")},
                                             {"ground", device::ContactKind::ohmic, patch_nodes(gm->mesh, "ground")}}};
    d.interface_certificate = gm->certificate;
    const auto dev = device::Device::create(std::move(d));
    REQUIRE(dev.has_value());
    const std::vector<double> bias{0.0, 1.0, 0.0};
    const auto p = solve::solve_bias(*dev, bias);
    REQUIRE(p.has_value());
    const auto& psi = p->fields.potential_V;
    double psi0 = 0.0, psi1 = 0.0;
    for (const mesh::NodeId v : patch_nodes(gm->mesh, "e1")) psi0 = psi[static_cast<std::size_t>(v)];
    for (const mesh::NodeId v : patch_nodes(gm->mesh, "e2")) psi1 = psi[static_cast<std::size_t>(v)];
    REQUIRE(std::abs(psi1 - psi0 - 1.0) < 1e-12);
    double worst = 0.0;
    std::size_t oxide_nodes = 0;
    for (std::size_t i = 0; i < gm->mesh.node_count(); ++i) {
        if (gm->node_region[i] != 0) continue;
        ++oxide_nodes;
        const double exact = psi0 + (psi1 - psi0) * gm->mesh.points()[i][0] / L;
        worst = std::max(worst, std::abs(psi[i] - exact));
    }
    CHECK(oxide_nodes > 20);
    INFO("worst deviation from the linear potential: " << worst << " V");
    CHECK(worst < 1e-10);
}

TEST_CASE("unstructured: a smooth-doped diode converges at the predefined orders", "[.unstructured]") {
    const double W = 2e-4, H = 1e-4;
    // Acceptors 1e16 everywhere; a donor Gaussian from the right edge (sigma 30 um... 3e-5 cm),
    // limited in y by an erfc window: a 2D junction with smooth doping.
    std::vector<device::DopingProfile> profiles{{.dopant = device::Dopant::acceptor, .peak_cm3 = 1e16}};
    device::DopingProfile n{.dopant = device::Dopant::donor, .peak_cm3 = 1e18, .shape = device::DepthShape::gaussian,
                            .axis = 0, .center = W, .width = 3e-5, .lateral = 1e-5};
    n.hi[1] = 6e-5;
    profiles.push_back(n);
    const std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    const auto build = [&](const mesh::Mesh& m, std::vector<mesh::NodeId> anode, std::vector<mesh::NodeId> cathode) {
        const std::vector<device::RegionId> region(m.node_count(), 0);
        const auto doping = device::sample(profiles, m, region, regions);
        REQUIRE(doping.has_value());
        return *device::Device::create({.mesh = m,
                                        .temperature_K = 300.0,
                                        .regions = regions,
                                        .node_region = region,
                                        .donors = doping->donors,
                                        .acceptors = doping->acceptors,
                                        .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                     {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    };
    // The reference: a uniform tensor grid 4x finer than the finest unstructured mesh.
    const std::vector<double> rx = uniform(W, 400), ry = uniform(H, 200);
    const mesh::Mesh ref_mesh = *mesh::make_tensor_grid(rx, ry);
    std::vector<mesh::NodeId> ref_cathode;
    for (const mesh::NodeId v : patch_nodes(ref_mesh, "x_max")) {
        if (ref_mesh.points()[static_cast<std::size_t>(v)][1] <= 4e-5) ref_cathode.push_back(v);
    }
    const device::Device ref = build(ref_mesh, patch_nodes(ref_mesh, "x_min"), ref_cathode);
    const auto ref_eq = solve::solve_equilibrium(ref);
    REQUIRE(ref_eq.has_value());
    const std::vector<double> forward{0.3, 0.0};
    const auto ref_on = solve::solve_bias(ref, forward);
    REQUIRE(ref_on.has_value());
    const double I_ref = ref_on->terminal_current[0];

    std::vector<double> h1, l2, cur;
    for (const double h : {8e-6, 4e-6, 2e-6}) {
        mesh::Geometry2D g{{rect(0, 0, W, H, 0)}, {{"anode", {0, 0, 0}, {0, H, 0}}, {"cathode", {W, 0, 0}, {W, 4e-5, 0}}}};
        mesh::MeshOptions2D o;
        o.max_size = h;
        const auto gm = mesh::generate_mesh(g, o);
        REQUIRE(gm.has_value());
        const device::Device d = build(gm->mesh, patch_nodes(gm->mesh, "anode"), patch_nodes(gm->mesh, "cathode"));
        const auto eq = solve::solve_equilibrium(d);
        REQUIRE(eq.has_value());
        std::vector<double> e(gm->mesh.node_count());
        double l2_err = 0.0, l2_ref = 0.0;
        for (std::size_t i = 0; i < e.size(); ++i) {
            const auto& p = gm->mesh.points()[i];
            const double r = bilinear(ref_eq->fields.potential_V, rx, ry, p[0], p[1]);
            e[i] = eq->fields.potential_V[i] - r;
            l2_err += gm->mesh.volumes()[i] * e[i] * e[i];
            l2_ref += gm->mesh.volumes()[i] * r * r;
        }
        double h1_err = 0.0;
        for (const mesh::Edge& ed : gm->mesh.edges()) {
            const double de = e[static_cast<std::size_t>(ed.second)] - e[static_cast<std::size_t>(ed.first)];
            h1_err += ed.coupling_area * de * de / ed.length;
        }
        const auto on = solve::solve_bias(d, forward);
        REQUIRE(on.has_value());
        h1.push_back(std::sqrt(h1_err));
        l2.push_back(std::sqrt(l2_err / l2_ref));
        cur.push_back(std::abs(on->terminal_current[0] - I_ref) / std::abs(I_ref));
        std::printf("unstructured diode h %.1e: %zu nodes, H1 %.4e, L2 %.4e, current %.4e\n", h,
                    gm->mesh.node_count(), h1.back(), l2.back(), cur.back());
    }
    for (std::size_t k = 1; k < 3; ++k) {
        const double oh1 = std::log2(h1[k - 1] / h1[k]);
        const double ol2 = std::log2(l2[k - 1] / l2[k]);
        const double oc = std::log2(cur[k - 1] / cur[k]);
        std::printf("orders %zu: H1 %.3f, L2 %.3f, current %.3f\n", k, oh1, ol2, oc);
        CHECK(oh1 >= 0.9);
        CHECK(ol2 >= 1.5);
        CHECK(oc >= 0.9);
    }
}

TEST_CASE("unstructured: an abrupt diode against the 1D reference, versus a tensor grid", "[.unstructured]") {
    // N_A = N_D = 1e17, junction at x = 1e-4 of a 2e-4 long, 2e-5 wide bar; forward 0.5 V.
    const double L = 2e-4, Wd = 2e-5;
    const std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    const auto build = [&](mesh::Mesh m) {
        const std::size_t n = m.node_count();
        std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
        const bool grid = m.origin() == mesh::MeshOrigin::tensor_grid;
        auto anode = patch_nodes(m, grid ? "x_min" : "anode");
        auto cathode = patch_nodes(m, grid ? "x_max" : "cathode");
        return *device::Device::create({.mesh = std::move(m),
                                        .temperature_K = 300.0,
                                        .regions = regions,
                                        .node_region = std::vector<device::RegionId>(n, 0),
                                        .donors = donors,
                                        .acceptors = acceptors,
                                        .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                     {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    };
    const std::vector<double> forward{0.5, 0.0};
    // The 1D reference: the legacy junction mesh (1e-8 at the junction).
    const auto x1 = legacy_graded_mesh(L, {1e-4}, 1e-8, 1e-6, 1.2);
    const auto J_ref = solve::solve_bias(build(*mesh::make_tensor_grid(x1)), forward)->terminal_current[0];
    // Unstructured, with a size box at the junction.
    mesh::Geometry2D g{{rect(0, 0, L, Wd, 0)}, {{"anode", {0, 0, 0}, {0, Wd, 0}}, {"cathode", {L, 0, 0}, {L, Wd, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = 2e-6;
    o.boxes = {{{0.9e-4, 0, 0}, {1.1e-4, Wd, 0}, 3e-7}};
    const auto gm = mesh::generate_mesh(g, o);
    REQUIRE(gm.has_value());
    const std::size_t nodes = gm->mesh.node_count();
    const double J_u = solve::solve_bias(build(gm->mesh), forward)->terminal_current[0] / Wd;
    // A 2D tensor grid with the same node count: the legacy-style graded x axis, ny rows.
    // The same graded x spacing as the size box asks for, and as many rows as give the same count.
    const auto xt = legacy_graded_mesh(L, {1e-4}, 3e-7, 2e-6, 1.2);
    const std::size_t ny = std::max<std::size_t>(2, (nodes + xt.size() / 2) / xt.size());
    const std::vector<double> yt = uniform(Wd, ny - 1);
    const double J_t = solve::solve_bias(build(*mesh::make_tensor_grid(xt, yt)), forward)->terminal_current[0] / Wd;
    std::printf("abrupt diode: J_ref %.6e, unstructured %.6e (%zu nodes, error %.3e), tensor %.6e (%zu nodes, "
                "error %.3e)\n",
                J_ref, J_u, nodes, std::abs(J_u / J_ref - 1), J_t, xt.size() * ny, std::abs(J_t / J_ref - 1));
    CHECK(std::abs(J_u - J_ref) <= std::abs(J_t - J_ref));
}

namespace {

// The MOSFET of mosfet_test.cpp (Unit 12; legacy mosfet.build_mosfet): Lg 600 nm, Lsd 300 nm,
// depth 200 nm, N_A 1e17, n+ source and drain (1e19 peak, Gaussian sigma 50 nm in depth, erfc
// lateral roll-off sigma 10 nm), t_ox 5 nm lumped gate, contacts source, drain, body, gate.
constexpr double Lg = 6e-5, Lsd = 3e-5, depth = 2e-5, N_A = 1e17, t_ox = 5e-7;
constexpr double N_peak = 1e19, sigma_y = 5e-6, sigma_lat = 1e-6;

std::vector<device::DopingProfile> mosfet_doping() {
    const double lateral = std::sqrt(2.0) * sigma_lat;
    device::DopingProfile source{.dopant = device::Dopant::donor, .peak_cm3 = N_peak,
                                 .shape = device::DepthShape::gaussian, .axis = 1, .center = 0.0,
                                 .width = sigma_y, .lateral = lateral};
    source.hi[0] = Lsd;
    device::DopingProfile drain = source;
    drain.hi[0] = INFINITY;
    drain.lo[0] = Lsd + Lg;
    return {source, drain, {.dopant = device::Dopant::acceptor, .peak_cm3 = N_A}};
}

device::Device mosfet_on(mesh::Mesh m, std::vector<mesh::NodeId> source, std::vector<mesh::NodeId> drain,
                         std::vector<mesh::NodeId> gate, std::vector<mesh::NodeId> body, const char* gate_patch) {
    const std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    const std::vector<device::RegionId> region(m.node_count(), 0);
    const auto doping = device::sample(mosfet_doping(), m, region, regions);
    REQUIRE(doping.has_value());
    device::Contact g{"gate", device::ContactKind::gate, std::move(gate)};
    g.gate = {.boundary = gate_patch, .oxide_thickness_cm = t_ox};
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = regions,
                                    .node_region = region,
                                    .donors = doping->donors,
                                    .acceptors = doping->acceptors,
                                    .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                                                 {"drain", device::ContactKind::ohmic, std::move(drain)},
                                                 {"body", device::ContactKind::ohmic, std::move(body)},
                                                 std::move(g)}});
}

device::Device tensor_mosfet(int nx, int ny) {
    const double L = 2 * Lsd + Lg;
    const auto x = legacy_graded_mesh(L, {Lsd, Lsd + Lg}, L / (nx * 20), L / nx, 1.15);
    const auto y = legacy_graded_mesh(depth, 0.0, depth / (ny * 20), depth / ny, 1.15);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    std::vector<mesh::NodeId> source, drain, gate;
    for (const mesh::NodeId v : patch_nodes(m, "y_min")) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        (px <= Lsd ? source : (px >= Lsd + Lg ? drain : gate)).push_back(v);
    }
    auto body = patch_nodes(m, "y_max");
    return mosfet_on(std::move(m), source, drain, gate, body, "y_min");
}

struct Transfer {
    std::vector<double> Vg, Id, resolution;
    double vt = 0.0, gm_max = 0.0, ss = 0.0;
    std::size_t nodes = 0;
};

Transfer transfer(const device::Device& d) {
    const std::vector<double> Vg = [] {
        std::vector<double> v;
        for (int k = 0; k <= 20; ++k) v.push_back(-1.0 + 0.1 * k);
        return v;
    }();
    std::vector<std::vector<double>> points;
    for (const double v : Vg) points.push_back({0.0, 0.05, 0.0, v});
    const auto sweep = solve::sweep_bias(d, points);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());  // every point converged
    const auto curve = analysis::current_curve(*sweep, 3, 1);
    REQUIRE(curve.has_value());
    Transfer t{Vg, curve->y, curve->resolution};
    t.vt = analysis::threshold_max_gm(*curve, 0.05)->value;
    const auto gm = numpy_gradient(curve->y, Vg);
    t.gm_max = *std::max_element(gm.begin(), gm.end());
    t.ss = analysis::subthreshold_swing(*curve)->value;
    t.nodes = d.mesh().node_count();
    return t;
}

// The predefined acceptance: the unstructured error within k = 1.5 of the tensor error at no more
// nodes, plus the reference's Richardson estimate |Q_ref - Q_t| / 3.
bool accepted(double q_u, double q_t, double q_ref) {
    return std::abs(q_u - q_ref) <= 1.5 * std::abs(q_t - q_ref) + std::abs(q_ref - q_t) / 3.0;
}

}  // namespace

// UNMET GATE (DECISIONS.md U14), marked [!shouldfail]: measured 2026-10-14 at 8031 nodes against
// the tensor 9490, the on current is 4 % below the reference while the tensor grid is within
// 0.05 % (V_T 4.7 mV off against 0.3 mV; SS within 0.04 mV/decade for both). The tensor grid has
// 2e-8 cm cells normal to the surface; the generator makes isotropic triangles, and resolving the
// inversion layer isotropically needs more nodes than the budget. If this test starts passing,
// Catch2 reports it as a failure, so the marker must then be removed.
TEST_CASE("unstructured: the MOSFET against the tensor MOSFET (predefined k = 1.5)",
          "[.unstructured][!shouldfail]") {
    const double L = 2 * Lsd + Lg;
    mesh::Geometry2D g{{rect(0, 0, L, depth, 0)},
                       {{"source", {0, 0, 0}, {Lsd, 0, 0}}, {"channel", {Lsd, 0, 0}, {Lsd + Lg, 0, 0}},
                        {"drain", {Lsd + Lg, 0, 0}, {L, 0, 0}}, {"body", {0, depth, 0}, {L, depth, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = 3e-6;
    o.boxes = {{{Lsd - 2e-6, 0, 0}, {Lsd + Lg + 2e-6, 3e-7, 0}, 1e-7}, {{0, 0, 0}, {L, 3e-6, 0}, 6e-7}};
    const auto gm = mesh::generate_mesh(g, o);
    REQUIRE(gm.has_value());
    const auto& m = gm->mesh;
    std::vector<mesh::NodeId> source = patch_nodes(m, "source"), drain = patch_nodes(m, "drain"), gate;
    for (const mesh::NodeId v : patch_nodes(m, "channel")) {
        if (!std::ranges::binary_search(source, v) && !std::ranges::binary_search(drain, v)) gate.push_back(v);
    }
    const device::Device unstructured = mosfet_on(m, source, drain, gate, patch_nodes(m, "body"), "channel");

    const Transfer u = transfer(unstructured);
    const Transfer t = transfer(tensor_mosfet(90, 50));
    const Transfer ref = transfer(tensor_mosfet(180, 100));
    std::printf("mosfet nodes: unstructured %zu, tensor %zu, reference %zu\n", u.nodes, t.nodes, ref.nodes);
    std::printf("mosfet V_T: u %.5f t %.5f ref %.5f; gm_max: u %.5e t %.5e ref %.5e; SS: u %.3f t %.3f ref %.3f\n",
                u.vt, t.vt, ref.vt, u.gm_max, t.gm_max, ref.gm_max, u.ss, t.ss, ref.ss);
    REQUIRE(u.nodes <= t.nodes);  // the comparison is at no more nodes
    CHECK(accepted(u.vt, t.vt, ref.vt));
    CHECK(accepted(u.gm_max, t.gm_max, ref.gm_max));
    CHECK(accepted(u.ss, t.ss, ref.ss));
    std::size_t resolved = 0, passed = 0;
    for (std::size_t k = 0; k < ref.Id.size(); ++k) {
        if (!(ref.Id[k] > ref.resolution[k] && t.Id[k] > t.resolution[k] && u.Id[k] > u.resolution[k])) continue;
        ++resolved;
        const bool ok = accepted(u.Id[k], t.Id[k], ref.Id[k]);
        passed += ok ? 1 : 0;
        std::printf("  Vg %.1f: Id u %.6e t %.6e ref %.6e %s\n", u.Vg[k], u.Id[k], t.Id[k], ref.Id[k], ok ? "" : "FAIL");
    }
    CHECK(resolved > 5);
    CHECK(passed == resolved);
}
