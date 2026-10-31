// FC8 (DECISIONS.md Unit 53; owner rulings SQ3, SQ5, SQ7, 2026-10-30): cross-section convergence of
// the complex generator, the gate for FD4 and FD7 together. A bulk-style tapered fin, uniform along
// z: silicon (a substrate under a fin, vertical beside the STI and tapered 3.6 degrees above it: a
// side straight through the STI top would make two touching pieces parallel only to rounding,
// which the 4 s rule refuses), STI oxide beside the fin's foot, a 1 nm interfacial oxide and 2 nm
// of HfO2 wrapped over the fin and the STI, an ohmic body contact on the bottom face (parallel to
// the axis, so the exact 3D solution is the 2D one; the end faces are Neumann). Gate at 0.5 V, the
// body at 0: equilibrium Poisson.
//
// The gate (SQ7): a flat top face over an oxide fill above the HfO2. A gate wrapped on the HfO2's
// outer surface bends 266 degrees at the fin's foot, a corner where the solution goes as r^0.68,
// so no quasi-uniform mesh reaches the orders below (generate2d measured H1 0.86, L2 1.15 and 1.36
// there); that case is run and reported, ungated.
//
// 3D sectioned extrusions at (h, s), (h/2, s/2), (h/4, s/4) against a 2D generate2d solve of the
// cross-section at (h/16, s/16): discrete H1 order of the potential >= 0.9, relative L2 >= 1.5, the
// gate charge per unit length >= 0.9 (Unit 16's thresholds, as tests/solve/unstructured3d_test.cpp);
// every 3D mesh certified (P4, P6). The silicon is one cell (two same-region cells meeting on a
// boundary face would put a vertex on a trace, a slice-5 limitation), and it reaches the whole
// bottom face so the 2D body patch has no end at an interface junction (refused by generate2d).
// Hidden: beyond the ctest node budget (M13), run explicitly in Release.
//
// FC8 fails on the flat gate (H1 orders 0.43, -0.04; L2 1.02, 0.47; gate charge 0.91, 2.17): an
// error floor of about 5 mV in the interfacial oxide at the fin's top corners, present only with a
// permittivity contrast. Owner ruling 2026-10-31 ("skip fc8"): an accepted limitation, the case
// [!shouldfail] with its thresholds unchanged (DECISIONS.md Unit 53, slice 6).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"

using namespace NiTCAD;

namespace {

constexpr double nm = 1e-7;
constexpr double length = 4 * nm;  // along z
constexpr double gate_V = 0.5;
constexpr double top = 25 * nm;    // the flat gate's face

enum class Gate { flat, wrapped };

// The cross-section in (x, y) [nm]: one polygon per region, shared vertices bit-identical.
struct Section {
    std::vector<mesh::Point> silicon, sti_left, sti_right, interfacial, hfo2, fill;  // fill: flat gate only
    std::vector<std::array<mesh::Point, 2>> gate;  // the gate's boundary edges
};

Section section(Gate gate) {
    const auto P = [](double x, double y) { return mesh::Point{x * nm, y * nm, 0.0}; };
    Section s;
    s.silicon = {P(0, 0),  P(20, 0),   P(20, 4),   P(13.75, 4), P(13.75, 8),
                 P(13, 20), P(7, 20),  P(6.25, 8), P(6.25, 4),  P(0, 4)};
    s.sti_left = {P(0, 4), P(6.25, 4), P(6.25, 8), P(5.25, 8), P(0, 8)};
    s.sti_right = {P(13.75, 4), P(20, 4), P(20, 8), P(14.75, 8), P(13.75, 8)};
    s.interfacial = {P(5.25, 8),  P(6.25, 8),  P(7, 20),        P(13, 20),
                     P(13.75, 8), P(14.75, 8), P(13.9375, 21), P(6.0625, 21)};
    s.hfo2 = {P(0, 8),   P(5.25, 8),   P(6.0625, 21),  P(13.9375, 21), P(14.75, 8),  P(20, 8),
              P(20, 10), P(16.625, 10), P(15.8125, 23), P(4.1875, 23), P(3.375, 10), P(0, 10)};
    // The HfO2's outer surface, from the domain's left side to its right.
    const std::vector<mesh::Point> outer{P(0, 10), P(3.375, 10), P(4.1875, 23), P(15.8125, 23), P(16.625, 10),
                                         P(20, 10)};
    if (gate == Gate::wrapped) {
        for (std::size_t k = 0; k + 1 < outer.size(); ++k) s.gate.push_back({outer[k], outer[k + 1]});
    } else {
        s.fill = outer;
        s.fill.push_back(P(20, 25));
        s.fill.push_back(P(0, 25));
        s.gate.push_back({P(0, 25), P(20, 25)});
    }
    return s;
}

// Regions: 0 silicon, 1 STI, 2 interfacial oxide, 3 HfO2, 4 the oxide fill (flat gate only).
std::vector<device::Region> regions(Gate gate) {
    physics::InsulatorParameters hfo2 = physics::silicon_dioxide_parameters;
    hfo2.eps_r = 22.0;
    std::vector<device::Region> r{{"silicon", physics::silicon()},
                                  {"sti", physics::silicon_dioxide()},
                                  {"interfacial", physics::silicon_dioxide()},
                                  {"hfo2", *physics::Insulator::create(hfo2)}};
    if (gate == Gate::flat) r.push_back({"fill", physics::silicon_dioxide()});
    return r;
}

std::vector<mesh::NodeId> sorted_unique(std::vector<mesh::NodeId> v) {
    std::ranges::sort(v);
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

device::Device make_device(Gate kind, const mesh::GeneratedMesh& g, std::vector<mesh::NodeId> body,
                           std::vector<mesh::NodeId> gate) {
    const std::size_t n = g.mesh.node_count();
    std::vector<device::RegionId> region(g.node_region.begin(), g.node_region.end());
    std::vector<double> acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (region[i] == 0) acceptors[i] = 1e18;
    }
    device::Contact gc{"gate", device::ContactKind::electrode, std::move(gate)};
    gc.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto d = device::Device::create({.mesh = g.mesh,
                                     .temperature_K = 300.0,
                                     .regions = regions(kind),
                                     .node_region = std::move(region),
                                     .donors = std::vector<double>(n, 0.0),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"body", device::ContactKind::ohmic, std::move(body)},
                                                  std::move(gc)},
                                     .interface_certificate = g.certificate});
    INFO((d.has_value() ? std::string() : d.error().message));
    REQUIRE(d.has_value());
    return std::move(*d);
}

struct Solution {
    std::vector<double> psi;
    double gate_charge = 0.0;  // per unit length along the axis [C/cm]
};

Solution solve_gated(const device::Device& d, double per_length) {
    solve::BiasOptions o;
    o.equations = solve::Equations::equilibrium_poisson;
    const auto r = solve::solve_bias(d, std::vector<double>{0.0, gate_V}, o);
    INFO((r.has_value() ? std::string() : r.error().message));
    REQUIRE(r.has_value());
    return {r->fields.potential_V, r->gate_charge[1] / per_length};
}

// The 2D reference and a point locator over its triangles (a uniform bucket grid).
struct Reference {
    std::optional<mesh::GeneratedMesh> mesh;
    Solution solution;
    double cell = 0;
    std::size_t nx = 0, ny = 0;
    std::vector<std::vector<std::size_t>> buckets;

    // The linear interpolant at (px, py): the triangle holding it, or the nearest by barycentric
    // coordinates among the neighbouring buckets' (a point on the domain to rounding).
    double at(double px, double py) const {
        const auto& pts = mesh->simplices.points;
        const auto& tri = mesh->simplices.cells;
        const auto clamp = [](double v, std::size_t n) {
            return std::min(static_cast<std::size_t>(std::max(0.0, std::floor(v))), n - 1);
        };
        const std::size_t bi = clamp(px / cell, nx), bj = clamp(py / cell, ny);
        double best = -std::numeric_limits<double>::infinity(), value = 0.0;
        for (std::size_t dj = (bj > 0 ? bj - 1 : 0); dj <= std::min(bj + 1, ny - 1); ++dj) {
            for (std::size_t di = (bi > 0 ? bi - 1 : 0); di <= std::min(bi + 1, nx - 1); ++di) {
                for (const std::size_t t : buckets[di + nx * dj]) {
                    const auto& a = pts[static_cast<std::size_t>(tri[3 * t])];
                    const auto& b = pts[static_cast<std::size_t>(tri[3 * t + 1])];
                    const auto& c = pts[static_cast<std::size_t>(tri[3 * t + 2])];
                    const double det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
                    const double l1 = ((px - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (py - a[1])) / det;
                    const double l2 = ((b[0] - a[0]) * (py - a[1]) - (px - a[0]) * (b[1] - a[1])) / det;
                    const double l0 = 1.0 - l1 - l2;
                    const double worst = std::min({l0, l1, l2});
                    if (worst > best) {
                        best = worst;
                        const auto f = [&](std::size_t k) {
                            return solution.psi[static_cast<std::size_t>(tri[3 * t + k])];
                        };
                        value = l0 * f(0) + l1 * f(1) + l2 * f(2);
                    }
                    if (best >= 0.0) return value;
                }
            }
        }
        REQUIRE(best > -1e-6);
        return value;
    }
};

Reference build_2d(Gate kind, double h, double s) {
    const Section sec = section(kind);
    mesh::Geometry2D g;
    g.polygons = {{sec.silicon, 0}, {sec.sti_left, 1}, {sec.sti_right, 1}, {sec.interfacial, 2}, {sec.hfo2, 3}};
    if (kind == Gate::flat) g.polygons.push_back({sec.fill, 4});
    g.patches.push_back({"body", sec.silicon[0], sec.silicon[1]});
    for (std::size_t k = 0; k < sec.gate.size(); ++k) {
        g.patches.push_back({"gate" + std::to_string(k), sec.gate[k][0], sec.gate[k][1]});
    }
    mesh::MeshOptions2D o;
    o.max_size = h;
    o.interface_resolution = s;
    const auto t0 = std::chrono::steady_clock::now();
    auto gm = mesh::generate_mesh(g, o);
    INFO((gm.has_value() ? std::string() : gm.error().message));
    REQUIRE(gm.has_value());
    Reference out;
    out.mesh = std::move(*gm);
    std::vector<mesh::NodeId> gate;
    for (std::size_t k = 0; k < sec.gate.size(); ++k) {
        const auto* p = out.mesh->mesh.find_boundary("gate" + std::to_string(k));
        REQUIRE(p != nullptr);
        gate.insert(gate.end(), p->nodes.begin(), p->nodes.end());
    }
    const device::Device d =
        make_device(kind, *out.mesh, out.mesh->mesh.find_boundary("body")->nodes, sorted_unique(std::move(gate)));
    out.solution = solve_gated(d, 1.0);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("2D h %.4g nm, s %.4g nm: %zu nodes, gate charge %.6e C/cm, %.1f s\n", h / nm, s / nm,
                out.mesh->mesh.node_count(), out.solution.gate_charge, secs);
    // The bucket grid: cells 4 h, each triangle in every bucket its bounding box meets.
    out.cell = 4.0 * h;
    out.nx = static_cast<std::size_t>(std::ceil(20 * nm / out.cell)) + 1;
    out.ny = static_cast<std::size_t>(std::ceil(top / out.cell)) + 1;
    out.buckets.resize(out.nx * out.ny);
    const auto& pts = out.mesh->simplices.points;
    const auto& tri = out.mesh->simplices.cells;
    for (std::size_t t = 0; 3 * t < tri.size(); ++t) {
        double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
        for (std::size_t k = 0; k < 3; ++k) {
            const auto& p = pts[static_cast<std::size_t>(tri[3 * t + k])];
            for (std::size_t c = 0; c < 2; ++c) {
                lo[c] = std::min(lo[c], p[c]);
                hi[c] = std::max(hi[c], p[c]);
            }
        }
        const auto i0 = static_cast<std::size_t>(lo[0] / out.cell), i1 = static_cast<std::size_t>(hi[0] / out.cell);
        const auto j0 = static_cast<std::size_t>(lo[1] / out.cell), j1 = static_cast<std::size_t>(hi[1] / out.cell);
        for (std::size_t j = j0; j <= std::min(j1, out.ny - 1); ++j) {
            for (std::size_t i = i0; i <= std::min(i1, out.nx - 1); ++i) out.buckets[i + out.nx * j].push_back(t);
        }
    }
    return out;
}

// H1 (per unit length) and relative L2 errors of the potential, and the relative gate charge error.
std::array<double, 3> errors(const mesh::Mesh& m, const Solution& u, const Reference& ref, double per_length) {
    std::vector<double> e(m.node_count());
    double l2_err = 0.0, l2_ref = 0.0;
    for (std::size_t i = 0; i < e.size(); ++i) {
        const auto& p = m.points()[i];
        const double r = ref.at(p[0], p[1]);
        e[i] = u.psi[i] - r;
        l2_err += m.volumes()[i] * e[i] * e[i];
        l2_ref += m.volumes()[i] * r * r;
    }
    double h1 = 0.0;
    for (const mesh::Edge& ed : m.edges()) {
        const double de = e[static_cast<std::size_t>(ed.second)] - e[static_cast<std::size_t>(ed.first)];
        h1 += ed.coupling_area * de * de / ed.length;
    }
    return {std::sqrt(h1 / per_length), std::sqrt(l2_err / l2_ref),
            std::abs(u.gate_charge - ref.solution.gate_charge) / std::abs(ref.solution.gate_charge)};
}

constexpr double h0 = 2 * nm, s0 = 0.16 * nm;

// The three 3D meshes against the reference; their errors in order.
std::vector<std::array<double, 3>> ladder(Gate kind, const char* label) {
    const Reference ref = build_2d(kind, h0 / 16, s0 / 16);
    const Section sec = section(kind);
    mesh::ExtrudedGeometry x;
    x.axis = 2;
    x.stations = {0.0, length};
    x.cells = {{"silicon", sec.silicon}, {"sti_left", sec.sti_left}, {"sti_right", sec.sti_right},
               {"interfacial", sec.interfacial}, {"hfo2", sec.hfo2}};
    x.regions = {{0}, {1}, {1}, {2}, {3}};
    if (kind == Gate::flat) {
        x.cells.push_back({"fill", sec.fill});
        x.regions.push_back({4});
        x.patches = {{"body", {0, 0, 0}, {20 * nm, 0, length}}, {"gate", {0, top, 0}, {20 * nm, top, length}}};
    } else {
        x.patches = {{"body", {0, 0, 0}, {20 * nm, 0, length}}, {"gate", {0, 10 * nm, 0}, {20 * nm, 23 * nm, length}}};
    }
    const auto complex = mesh::extrude(x);
    INFO((complex.has_value() ? std::string() : complex.error().message));
    REQUIRE(complex.has_value());
    std::vector<std::array<double, 3>> err;
    for (const int k : {1, 2, 4}) {
        mesh::MeshOptions3D o;
        o.max_size = h0 / k;
        o.interface_resolution = s0 / k;
        const auto t0 = std::chrono::steady_clock::now();
        const auto gm = mesh::generate_mesh(*complex, o);
        INFO((gm.has_value() ? std::string() : gm.error().message));
        REQUIRE(gm.has_value());
        REQUIRE(gm->certificate.has_value());
        const double mesh_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const device::Device d =
            make_device(kind, *gm, gm->mesh.find_boundary("body")->nodes, gm->mesh.find_boundary("gate")->nodes);
        const auto t1 = std::chrono::steady_clock::now();
        const Solution u = solve_gated(d, length);
        const double solve_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
        err.push_back(errors(gm->mesh, u, ref, length));
        std::printf("%s 3D h %.3g nm, s %.3g nm: %zu nodes, %zu tetrahedra, %zu slivers, mesh %.1f s, solve %.1f s; "
                    "H1 %.4e, L2 %.4e, gate charge %.6e C/cm (error %.4e)\n",
                    label, o.max_size / nm, o.interface_resolution / nm, gm->mesh.node_count(), gm->report.tetrahedra,
                    gm->report.slivers, mesh_s, solve_s, err.back()[0], err.back()[1], u.gate_charge, err.back()[2]);
        std::fflush(stdout);
    }
    return err;
}

std::array<double, 3> orders(const std::vector<std::array<double, 3>>& err, std::size_t k, const char* label) {
    const std::array<double, 3> o{std::log2(err[k - 1][0] / err[k][0]), std::log2(err[k - 1][1] / err[k][1]),
                                  std::log2(err[k - 1][2] / err[k][2])};
    std::printf("%s orders %zu: H1 %.3f, L2 %.3f, gate charge %.3f\n", label, k, o[0], o[1], o[2]);
    return o;
}

}  // namespace

TEST_CASE("complex convergence: a gated tapered fin uniform along its axis against its 2D cross-section (FC8)",
          "[.complex_convergence][!shouldfail]") {
    const auto err = ladder(Gate::flat, "FC8");
    for (std::size_t k = 1; k < err.size(); ++k) {
        const auto o = orders(err, k, "FC8");
        CHECK(o[0] >= 0.9);
        CHECK(o[1] >= 1.5);
        CHECK(o[2] >= 0.9);
    }
}

// SQ7: the wrapped gate's 266-degree corner caps the orders; measured and reported, not gated.
TEST_CASE("complex convergence: the same fin under a wrapped gate, reported (FC8, SQ7)", "[.complex_convergence]") {
    const auto err = ladder(Gate::wrapped, "FC8 wrapped");
    for (std::size_t k = 1; k < err.size(); ++k) (void)orders(err, k, "FC8 wrapped");
}

// UC5 (DECISIONS.md Unit 54): the flat-gate fin above, swept along z on three lines (Unit 54's swept
// mesh), against an equilibrium solve on the very 2D mesh it was swept from. The product dual
// (UD1) makes the uniform extension of the 2D discrete solution solve the 3D discrete equations, so
// the potential on every line equals the 2D one, and the gate charge per unit length the 2D charge,
// within the solver's tolerance.
TEST_CASE("swept mesh: a fin uniform along the axis solves as its 2D cross-section (UC5)", "[solve][swept]") {
    const Section sec = section(Gate::flat);
    mesh::ExtrudedGeometry x;
    x.axis = 2;
    x.stations = {0.0, length};
    x.cells = {{"silicon", sec.silicon}, {"sti_left", sec.sti_left}, {"sti_right", sec.sti_right},
               {"interfacial", sec.interfacial}, {"hfo2", sec.hfo2}, {"fill", sec.fill}};
    x.regions = {{0}, {1}, {1}, {2}, {3}, {4}};
    x.patches = {{"body", {0, 0, 0}, {20 * nm, 0, length}}, {"gate", {0, top, 0}, {20 * nm, top, length}}};
    mesh::SweptOptions o;
    o.section.max_size = 2 * nm;
    o.section.interface_resolution = 0.16 * nm;
    o.lines = {0.0, 0.5 * length, length};
    const auto gm = mesh::generate_swept(x, o);
    INFO((gm.has_value() ? std::string() : gm.error().message));
    REQUIRE(gm.has_value());
    REQUIRE(gm->certificate.has_value());
    const device::Device d3 = make_device(Gate::flat, *gm, gm->mesh.find_boundary("body")->nodes, gm->mesh.find_boundary("gate")->nodes);
    const Solution u3 = solve_gated(d3, length);

    // The same cross-section mesh, from the same polygons, 2D regions and options (UD2).
    mesh::Geometry2D g2;
    g2.polygons = {{sec.silicon, 0}, {sec.sti_left, 1}, {sec.sti_right, 1}, {sec.interfacial, 2}, {sec.hfo2, 3}, {sec.fill, 4}};
    g2.patches = {{"body", sec.silicon[0], sec.silicon[1]}, {"gate", sec.gate[0][0], sec.gate[0][1]}};
    const auto m2 = mesh::generate_mesh(g2, o.section);
    REQUIRE(m2.has_value());
    const std::size_t n2 = m2->mesh.node_count();
    REQUIRE(gm->mesh.node_count() == 3 * n2);
    std::vector<mesh::NodeId> body, gate;
    for (std::size_t i = 0; i < n2; ++i) {
        const auto& p = m2->mesh.points()[i];
        REQUIRE(p[0] == gm->mesh.points()[i][0]);
        REQUIRE(p[1] == gm->mesh.points()[i][1]);
        if (p[1] == 0.0) body.push_back(static_cast<mesh::NodeId>(i));
        if (p[1] == top) gate.push_back(static_cast<mesh::NodeId>(i));
    }
    const device::Device d2 = make_device(Gate::flat, *m2, body, gate);
    const Solution u2 = solve_gated(d2, 1.0);
    double worst = 0.0, scale = 0.0;
    for (std::size_t j = 0; j < 3; ++j) {
        for (std::size_t i = 0; i < n2; ++i) {
            worst = std::max(worst, std::abs(u3.psi[j * n2 + i] - u2.psi[i]));
            scale = std::max(scale, std::abs(u2.psi[i]));
        }
    }
    const double charge = std::abs(u3.gate_charge - u2.gate_charge) / std::abs(u2.gate_charge);
    std::printf("UC5: %zu 2D nodes x 3 lines; potential deviation %.3e V of %.3e; gate charge %.6e against %.6e C/cm (%.3e)\n",
                n2, worst, scale, u3.gate_charge, u2.gate_charge, charge);
    CHECK(worst <= 1e-9 * scale);
    CHECK(charge <= 1e-9);
}
