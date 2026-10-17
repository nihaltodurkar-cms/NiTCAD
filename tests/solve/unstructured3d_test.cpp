// Device and solver gates on 3D meshes (Unit 16, commit 6e; DECISIONS.md "Gates", M10, M13). The
// thresholds are those predefined in DECISIONS.md; nothing here was tuned after a run.
// - [unstructured_fast] (ctest): Kuhn-split boxes give the tensor grid's diode solve to 1e-12.
// - [.unstructured3d] (run explicitly; beyond the ctest node budget of M13), the owner's choice of
//   2026-10-15: a diode with doping smooth in x and y and uniform in z, on generated meshes at h,
//   h/2, h/4, against a 2D tensor solve 4x finer than the finest: discrete H1 order >= 0.9, L2 >=
//   1.5, terminal current order >= 0.9; and the same gates on perturbed lattices (the sliver
//   stress test of M10), with their sliver counts.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "assembly_timing.hpp"
#include "triangulation3d.hpp"

using namespace NiTCAD;

namespace {

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

// Bilinear interpolation of a node field on a uniform 2D tensor grid (x fastest).
double bilinear(const std::vector<double>& f, const std::vector<double>& x, const std::vector<double>& y, double px,
                double py) {
    const double hx = x[1] - x[0], hy = y[1] - y[0];
    const auto i = std::min(static_cast<std::size_t>(std::max(0.0, std::floor(px / hx))), x.size() - 2);
    const auto j = std::min(static_cast<std::size_t>(std::max(0.0, std::floor(py / hy))), y.size() - 2);
    const double s = (px - x[i]) / hx, t = (py - y[j]) / hy;
    const auto at = [&](std::size_t a, std::size_t b) { return f[a + x.size() * b]; };
    return (1 - s) * (1 - t) * at(i, j) + s * (1 - t) * at(i + 1, j) + (1 - s) * t * at(i, j + 1) + s * t * at(i + 1, j + 1);
}

// The smallest dihedral angle of a tetrahedron [degrees].
double min_dihedral(const std::array<mesh::Point, 4>& p) {
    double best = 180.0;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = i + 1; j < 4; ++j) {
            std::size_t k = 0;
            while (k == i || k == j) ++k;
            std::size_t l = k + 1;
            while (l == i || l == j) ++l;
            std::array<double, 3> e{}, u{}, v{};
            for (std::size_t c = 0; c < 3; ++c) e[c] = p[j][c] - p[i][c];
            const double ee = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
            const auto perp = [&](const mesh::Point& q, std::array<double, 3>& out) {
                std::array<double, 3> r{};
                for (std::size_t c = 0; c < 3; ++c) r[c] = q[c] - p[i][c];
                const double t = (r[0] * e[0] + r[1] * e[1] + r[2] * e[2]) / ee;
                for (std::size_t c = 0; c < 3; ++c) out[c] = r[c] - t * e[c];
            };
            perp(p[k], u);
            perp(p[l], v);
            const double c = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) /
                             std::sqrt((u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
            best = std::min(best, std::acos(std::clamp(c, -1.0, 1.0)) * 180.0 / 3.141592653589793);
        }
    }
    return best;
}

// The diode of the convergence gates: W x H x D, acceptors 1e16, a donor Gaussian from x = W
// limited in y; the anode on x = 0, the cathode on x = W for y <= 2e-5. Smooth in x and y,
// uniform in z, so the exact solution is the 2D one.
constexpr double W = 1e-4, H = 5e-5, D = 1e-5, cathode_top = 2e-5;

std::vector<device::DopingProfile> diode_profiles() {
    std::vector<device::DopingProfile> profiles{{.dopant = device::Dopant::acceptor, .peak_cm3 = 1e16}};
    device::DopingProfile n{.dopant = device::Dopant::donor, .peak_cm3 = 1e18, .shape = device::DepthShape::gaussian,
                            .axis = 0, .center = W, .width = 1.5e-5, .lateral = 5e-6};
    n.hi[1] = 3e-5;
    profiles.push_back(n);
    return profiles;
}

device::Device diode(const mesh::Mesh& m, std::vector<mesh::NodeId> anode, std::vector<mesh::NodeId> cathode) {
    const std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    const std::vector<device::RegionId> region(m.node_count(), 0);
    const auto doping = device::sample(diode_profiles(), m, region, regions);
    REQUIRE(doping.has_value());
    auto d = device::Device::create({.mesh = m,
                                     .temperature_K = 300.0,
                                     .regions = regions,
                                     .node_region = region,
                                     .donors = doping->donors,
                                     .acceptors = doping->acceptors,
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    REQUIRE(d.has_value());
    return std::move(*d);
}

// The 2D reference, 4x finer than the finest 3D mesh (h / 4 = 6.25e-7): equilibrium potential
// and the forward current per unit depth.
struct Reference {
    std::vector<double> x, y, psi;
    double current = 0.0;
};

const Reference& reference() {
    static const Reference r = [] {
        Reference out;
        out.x = uniform(W, 160);
        out.y = uniform(H, 80);
        const mesh::Mesh m = *mesh::make_tensor_grid(out.x, out.y);
        std::vector<mesh::NodeId> cathode;
        for (const mesh::NodeId v : patch_nodes(m, "x_max")) {
            if (m.points()[static_cast<std::size_t>(v)][1] <= cathode_top) cathode.push_back(v);
        }
        const device::Device dev = diode(m, patch_nodes(m, "x_min"), cathode);
        const auto eq = solve::solve_equilibrium(dev);
        REQUIRE(eq.has_value());
        out.psi = eq->fields.potential_V;
        const auto on = solve::solve_bias(dev, std::vector<double>{0.3, 0.0});
        REQUIRE(on.has_value());
        out.current = on->terminal_current[0];
        return out;
    }();
    return r;
}

// H1 and relative L2 errors of the equilibrium potential against the reference, and the relative
// error of the forward current (per unit depth).
std::array<double, 3> errors(const mesh::Mesh& m, const device::Device& d) {
    const Reference& ref = reference();
    const auto eq = solve::solve_equilibrium(d);
    REQUIRE(eq.has_value());
    std::vector<double> e(m.node_count());
    double l2_err = 0.0, l2_ref = 0.0;
    for (std::size_t i = 0; i < e.size(); ++i) {
        const auto& p = m.points()[i];
        const double r = bilinear(ref.psi, ref.x, ref.y, p[0], p[1]);
        e[i] = eq->fields.potential_V[i] - r;
        l2_err += m.volumes()[i] * e[i] * e[i];
        l2_ref += m.volumes()[i] * r * r;
    }
    double h1 = 0.0;
    for (const mesh::Edge& ed : m.edges()) {
        const double de = e[static_cast<std::size_t>(ed.second)] - e[static_cast<std::size_t>(ed.first)];
        h1 += ed.coupling_area * de * de / ed.length;
    }
    const auto on = solve::solve_bias(d, std::vector<double>{0.3, 0.0});
    REQUIRE(on.has_value());
    return {std::sqrt(h1 / D), std::sqrt(l2_err / l2_ref), std::abs(on->terminal_current[0] / D - ref.current) / std::abs(ref.current)};
}

void check_orders(const std::vector<std::array<double, 3>>& err, const char* label) {
    for (std::size_t k = 1; k < err.size(); ++k) {
        const double oh1 = std::log2(err[k - 1][0] / err[k][0]);
        const double ol2 = std::log2(err[k - 1][1] / err[k][1]);
        const double oc = std::log2(err[k - 1][2] / err[k][2]);
        std::printf("%s orders %zu: H1 %.3f, L2 %.3f, current %.3f\n", label, k, oh1, ol2, oc);
        CHECK(oh1 >= 0.9);
        CHECK(ol2 >= 1.5);
        CHECK(oc >= 0.9);
    }
}

// A lattice of spacing h on the diode's box with every interior point moved by up to 0.25 h in
// each coordinate (boundary points stay, so every boundary triangle keeps an empty diametral
// ball), tetrahedralized by the generator's Delaunay; patches anode and cathode.
mesh::SimplexMesh perturbed_lattice(double h, std::uint64_t seed) {
    const auto nx = static_cast<std::size_t>(std::lround(W / h));
    const auto ny = static_cast<std::size_t>(std::lround(H / h));
    const auto nz = static_cast<std::size_t>(std::lround(D / h));
    const std::vector<double> xs = uniform(W, nx), ys = uniform(H, ny), zs = uniform(D, nz);
    std::uint64_t state = seed;
    const auto next = [&]() {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return static_cast<double>((z ^ (z >> 31)) >> 11) * 0x1p-53;
    };
    mesh::SimplexMesh s;
    s.dimension = 3;
    for (std::size_t k = 0; k <= nz; ++k) {
        for (std::size_t j = 0; j <= ny; ++j) {
            for (std::size_t i = 0; i <= nx; ++i) {
                mesh::Point p{xs[i], ys[j], zs[k]};
                if (i > 0 && i < nx && j > 0 && j < ny && k > 0 && k < nz) {
                    for (double& c : p) c += (next() - 0.5) * 0.5 * h;
                }
                s.points.push_back(p);
            }
        }
    }
    mesh::detail::Triangulation3D t({0, 0, 0}, {W, H, D});
    for (const mesh::Point& p : s.points) t.insert(p);
    for (std::size_t c = 0; c < t.tet_slots(); ++c) {
        const auto& tt = t.tet(static_cast<int>(c));
        if (!tt.alive || std::ranges::any_of(tt.v, [&](int v) { return t.is_super(v); })) continue;
        for (const int v : tt.v) s.cells.push_back(static_cast<mesh::NodeId>(v - 4));
        for (std::size_t f = 0; f < 4; ++f) {
            std::array<mesh::NodeId, 3> fv{};
            for (std::size_t a = 0, m = 0; a < 4; ++a) {
                if (a != f) fv[m++] = static_cast<mesh::NodeId>(tt.v[a] - 4);
            }
            const auto P = [&](mesh::NodeId v) -> const mesh::Point& { return s.points[static_cast<std::size_t>(v)]; };
            const bool anode = std::ranges::all_of(fv, [&](mesh::NodeId v) { return P(v)[0] == 0.0; });
            const bool cathode =
                std::ranges::all_of(fv, [&](mesh::NodeId v) { return P(v)[0] == W && P(v)[1] <= cathode_top; });
            if (!anode && !cathode) continue;
            s.facets.insert(s.facets.end(), fv.begin(), fv.end());
            s.facet_patch.push_back(anode ? 0 : 1);
        }
    }
    s.patch_names = {"anode", "cathode"};
    return s;
}

}  // namespace

TEST_CASE("unstructured 3D: Kuhn-split boxes give the tensor grid's diode solve", "[solve][unstructured_fast]") {
    // Kuhn split (six tetrahedra per box along the main diagonal): every box is cospherical, so the
    // circumcentric dual is the tensor grid's box method exactly (diagonals have zero coupling).
    const std::vector<double> xs = uniform(W, 10), ys = uniform(H, 5), zs = uniform(D, 2);
    const mesh::Mesh tensor = *mesh::make_tensor_grid(xs, ys, zs);
    mesh::SimplexMesh s;
    s.dimension = 3;
    s.points.assign(tensor.points().begin(), tensor.points().end());
    const auto v = [&](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<mesh::NodeId>(i + xs.size() * (j + ys.size() * k));
    };
    const std::array<std::array<int, 3>, 6> perms{{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}};
    for (std::size_t k = 0; k + 1 < zs.size(); ++k) {
        for (std::size_t j = 0; j + 1 < ys.size(); ++j) {
            for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
                for (const auto& perm : perms) {
                    std::array<std::size_t, 3> at{i, j, k};
                    s.cells.push_back(v(at[0], at[1], at[2]));
                    for (const int axis : perm) {
                        ++at[static_cast<std::size_t>(axis)];
                        s.cells.push_back(v(at[0], at[1], at[2]));
                    }
                }
            }
        }
    }
    // The contacts: the boundary triangles on x = 0 and on x = W with y <= cathode_top.
    for (std::size_t c = 0; c < s.cells.size(); c += 4) {
        for (std::size_t f = 0; f < 4; ++f) {
            std::array<mesh::NodeId, 3> fv{};
            for (std::size_t a = 0, m = 0; a < 4; ++a) {
                if (a != f) fv[m++] = s.cells[c + a];
            }
            const auto P = [&](mesh::NodeId n) -> const mesh::Point& { return s.points[static_cast<std::size_t>(n)]; };
            const bool anode = std::ranges::all_of(fv, [&](mesh::NodeId n) { return P(n)[0] == 0.0; });
            const bool cathode = std::ranges::all_of(fv, [&](mesh::NodeId n) { return P(n)[0] == W && P(n)[1] <= cathode_top; });
            if (!anode && !cathode) continue;
            s.facets.insert(s.facets.end(), fv.begin(), fv.end());
            s.facet_patch.push_back(anode ? 0 : 1);
        }
    }
    s.patch_names = {"anode", "cathode"};
    const auto kuhn = mesh::make_dual_mesh(s);
    REQUIRE(kuhn.has_value());
    std::vector<mesh::NodeId> tensor_cathode;
    for (const mesh::NodeId n : patch_nodes(tensor, "x_max")) {
        if (tensor.points()[static_cast<std::size_t>(n)][1] <= cathode_top) tensor_cathode.push_back(n);
    }
    const device::Device a = diode(tensor, patch_nodes(tensor, "x_min"), tensor_cathode);
    const device::Device b = diode(*kuhn, patch_nodes(*kuhn, "anode"), patch_nodes(*kuhn, "cathode"));
    const auto ra = solve::solve_bias(a, std::vector<double>{0.3, 0.0});
    const auto rb = solve::solve_bias(b, std::vector<double>{0.3, 0.0});
    REQUIRE(ra.has_value());
    REQUIRE(rb.has_value());
    double worst = 0.0;
    for (std::size_t i = 0; i < tensor.node_count(); ++i) {
        worst = std::max(worst, std::abs(ra->fields.potential_V[i] - rb->fields.potential_V[i]));
    }
    INFO("potential difference " << worst << " V; currents " << ra->terminal_current[0] << " and " << rb->terminal_current[0]);
    CHECK(worst <= 1e-12);
    CHECK(std::abs(ra->terminal_current[0] - rb->terminal_current[0]) <= 1e-12 * std::abs(ra->terminal_current[0]));
}

TEST_CASE("unstructured 3D: a z-uniform diode on generated meshes converges at the predefined orders",
          "[.unstructured3d]") {
    std::vector<std::array<double, 3>> err;
    for (const double h : {1e-5, 5e-6, 2.5e-6}) {
        mesh::Geometry3D g;
        g.domain.vertices = {{0, 0, 0}, {W, 0, 0}, {W, H, 0}, {0, H, 0}, {0, 0, D}, {W, 0, D}, {W, H, D}, {0, H, D}};
        g.domain.faces = {{{0, 3, 7, 4}, "anode"}, {{1, 2, 6, 5}, "x_max"}, {{0, 1, 5, 4}, "y_min"},
                          {{3, 2, 6, 7}, "y_max"}, {{0, 1, 2, 3}, "z_min"}, {{4, 5, 6, 7}, "z_max"}};
        g.domain.faces[1].patches = {{"cathode", {{W, 0, 0}, {W, cathode_top, 0}, {W, cathode_top, D}, {W, 0, D}}},
                                     {"x_max_rest", {{W, cathode_top, 0}, {W, H, 0}, {W, H, D}, {W, cathode_top, D}}}};
        mesh::MeshOptions3D o;
        o.max_size = h;
        const auto gm = mesh::generate_mesh(g, o);
        REQUIRE(gm.has_value());
        const device::Device d = diode(gm->mesh, patch_nodes(gm->mesh, "anode"), patch_nodes(gm->mesh, "cathode"));
        err.push_back(errors(gm->mesh, d));
        std::printf("generated 3D diode h %.2e: %zu nodes, %zu slivers (min dihedral %.2f), H1 %.4e, L2 %.4e, current %.4e\n",
                    h, gm->mesh.node_count(), gm->report.slivers, gm->report.min_dihedral_deg, err.back()[0],
                    err.back()[1], err.back()[2]);
    }
    check_orders(err, "generated 3D diode");
}

TEST_CASE("unstructured 3D: perturbed lattices (sliver stress, M10) pass the same gates", "[.unstructured3d]") {
    std::vector<std::array<double, 3>> err;
    for (const double h : {1e-5, 5e-6, 2.5e-6}) {
        const mesh::SimplexMesh s = perturbed_lattice(h, 0x5EEDull);
        std::size_t slivers = 0;
        double least = 180.0;
        for (std::size_t c = 0; c < s.cells.size(); c += 4) {
            std::array<mesh::Point, 4> p{};
            for (std::size_t a = 0; a < 4; ++a) p[a] = s.points[static_cast<std::size_t>(s.cells[c + a])];
            const double dih = min_dihedral(p);
            least = std::min(least, dih);
            slivers += dih < 5.0 ? 1 : 0;
        }
        const auto m = mesh::make_dual_mesh(s);
        REQUIRE(m.has_value());
        double volume = 0.0;
        for (const double x : m->volumes()) volume += x;
        CHECK(std::abs(volume - W * H * D) <= 1e-12 * W * H * D);
        const device::Device d = diode(*m, patch_nodes(*m, "anode"), patch_nodes(*m, "cathode"));
        err.push_back(errors(*m, d));
        std::printf("perturbed lattice h %.2e: %zu nodes, %zu tetrahedra, %zu slivers (min dihedral %.2f), H1 %.4e, "
                    "L2 %.4e, current %.4e\n",
                    h, m->node_count(), s.cells.size() / 4, slivers, least, err.back()[0], err.back()[1], err.back()[2]);
        // A stress test only if it has slivers. The coarsest lattice has one cell through z, so no
        // interior point moves and it has none (found on the first run; this check was narrowed
        // to the meshes with interior points then, the gates unchanged).
        if (std::lround(D / h) >= 2) CHECK(slivers > 0);
    }
    check_orders(err, "perturbed lattice");
}

TEST_CASE("unstructured 3D: assembly cost on a generated mesh", "[.assembly_benchmark]") {
    // Unit 40 (DECISIONS.md FG4): the z-uniform diode on the generated mesh of size 5e-6 at 0.3 V.
    mesh::Geometry3D g;
    g.domain.vertices = {{0, 0, 0}, {W, 0, 0}, {W, H, 0}, {0, H, 0}, {0, 0, D}, {W, 0, D}, {W, H, D}, {0, H, D}};
    g.domain.faces = {{{0, 3, 7, 4}, "anode"}, {{1, 2, 6, 5}, "x_max"}, {{0, 1, 5, 4}, "y_min"},
                      {{3, 2, 6, 7}, "y_max"}, {{0, 1, 2, 3}, "z_min"}, {{4, 5, 6, 7}, "z_max"}};
    g.domain.faces[1].patches = {{"cathode", {{W, 0, 0}, {W, cathode_top, 0}, {W, cathode_top, D}, {W, 0, D}}},
                                 {"x_max_rest", {{W, cathode_top, 0}, {W, H, 0}, {W, H, D}, {W, cathode_top, D}}}};
    mesh::MeshOptions3D o;
    o.max_size = 5e-6;
    const auto gm = mesh::generate_mesh(g, o);
    REQUIRE(gm.has_value());
    const device::Device d = diode(gm->mesh, patch_nodes(gm->mesh, "anode"), patch_nodes(gm->mesh, "cathode"));
    const std::vector<double> bias{0.3, 0.0};
    const auto p = solve::solve_bias(d, bias);
    REQUIRE(p.has_value());
    const auto scaling = *assemble::make_scaling(d);
    auto system = *assemble::DriftDiffusion::create(d, scaling);
    REQUIRE(system.set_bias(bias).has_value());
    time_assembly("3D generated diode", system, scaled_state(system, p->fields, scaling));
}
