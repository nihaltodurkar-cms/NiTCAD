// Unit 64 IC3: a gate-all-around wire whose gate faces are holed interface pieces (slice 1, IP2), meshed by the
// polyhedral generator, against the same device meshed by the swept generator (Unit 54; IQ3 a). A p-type silicon
// wire (1e17) along y through an HfO2 block (0.2 to 0.8 µm in x and z) in oxide, an ohmic body contact on the wire's y = 0 end, a metal gate
// (4.4 eV) on the four faces along the wire at 0.5 V: equilibrium Poisson. The potential and carrier densities are compared
// at 200 seeded points (the densities inside the silicon), as linear interpolants on each mesh. The reference's own
// discretization is measured by the swept mesh at half its sizes, for the record.
//
// IC3 is missed (DECISIONS.md Unit 64, slice 1): the two generators converge to different potentials (80 mV in the
// HfO2 between wire and gate, the swept mesh still moving at h/4), and do so without the block, so without holes,
// on the polyhedral path slice 1 leaves bit-identical. The case is [!shouldfail], its 1e-6 thresholds unchanged, as
// Unit 53's FC8; the open question is which mesher, or this swept set-up, is off.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <utility>
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

constexpr double um = 1e-4;
constexpr double gate_V = 0.5;

// Regions: 0 oxide, 1 silicon, 2 HfO2.
std::vector<device::Region> regions() {
    physics::InsulatorParameters hfo2 = physics::silicon_dioxide_parameters;
    hfo2.eps_r = 22.0;
    return {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}, {"hfo2", *physics::Insulator::create(hfo2)}};
}

mesh::PatchBox body_box() { return {"body", {0.4 * um, 0, 0.4 * um}, {0.6 * um, 0, 0.6 * um}}; }
mesh::PatchBox far_box() { return {"far", {0.4 * um, 2 * um, 0.4 * um}, {0.6 * um, 2 * um, 0.6 * um}}; }  // the wire's other end, no contact
// The gate: the four faces along the wire (x = 0, x = 1, z = 0, z = 1), four boxes of one name; the end faces carry
// no contact but the body, so no contact edge lies near the wire.
std::vector<mesh::PatchBox> patches() {
    return {body_box(), far_box(), {"gate", {0, 0, 0}, {0, 2 * um, 1 * um}}, {"gate", {1 * um, 0, 0}, {1 * um, 2 * um, 1 * um}},
            {"gate", {0, 0, 0}, {1 * um, 2 * um, 0}}, {"gate", {0, 0, 1 * um}, {1 * um, 2 * um, 1 * um}}};
}

mesh::LayeredSolid box3(double x0, double x1, double y0, double y1, double z0, double z1, std::int32_t region) {
    mesh::LayeredSolid s;
    s.kind = mesh::LayeredSolid::Kind::box;
    s.lo = {x0 * um, y0 * um, z0 * um};
    s.hi = {x1 * um, y1 * um, z1 * um};
    s.region = region;
    return s;
}

mesh::GeneratedMesh polyhedral(double h, double s) {
    mesh::LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {1 * um, 2 * um, 1 * um}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(box3(0.2, 0.8, 0.5, 1.5, 0.2, 0.8, 2));
    g.solids.push_back(box3(0.4, 0.6, 0, 2, 0.4, 0.6, 1));
    g.patches = patches();
    const auto c = mesh::layered_partition3d(g);
    REQUIRE(c.has_value());
    mesh::MeshOptions3D o;
    o.max_size = h * um;
    o.interface_resolution = s * um;
    auto m = mesh::generate_mesh(*c, o);
    INFO((m.has_value() ? std::string() : m.error().message));
    REQUIRE(m.has_value());
    return std::move(*m);
}

// The same device as an extrusion along y, sections at the block's faces. The cross-section in (u, v) = (z, x): the
// wire's square, and two rings of four trapezoids (the block's, the oxide's) cut along the diagonals, so every
// internal edge joins two corners (the 2D generator's P4 failed where a grid's internal edges ended on an interface).
mesh::GeneratedMesh swept(double h, double s) {
    mesh::ExtrudedGeometry x;
    x.axis = 1;
    x.stations = {0, 0.5 * um, 1.5 * um, 2 * um};
    const auto P = [](double z, double xx) { return mesh::Point{z * um, xx * um, 0.0}; };
    x.cells.push_back({"wire", {P(0.4, 0.4), P(0.6, 0.4), P(0.6, 0.6), P(0.4, 0.6)}});
    x.regions.push_back({1, 1, 1});
    for (const auto& [a, b, block] : {std::tuple{0.2, 0.4, true}, std::tuple{0.0, 0.2, false}}) {
        const std::array<mesh::Point, 4> o{P(a, a), P(1 - a, a), P(1 - a, 1 - a), P(a, 1 - a)};
        const std::array<mesh::Point, 4> i{P(b, b), P(1 - b, b), P(1 - b, 1 - b), P(b, 1 - b)};
        for (std::size_t k = 0; k < 4; ++k) {
            x.cells.push_back({std::format("{}{}", block ? "block" : "oxide", k), {o[k], o[(k + 1) % 4], i[(k + 1) % 4], i[k]}});
            x.regions.push_back({0, block ? 2 : 0, 0});
        }
    }
    x.patches = patches();
    mesh::SweptOptions o;
    o.section.max_size = h * um;
    o.section.interface_resolution = s * um;
    // Lines at most h apart, those within s of the gate's faces left to the straddles there (+- s / 2).
    const int n = static_cast<int>(std::ceil(2.0 / h));
    for (int k = 0; k <= n; ++k) {
        const double y = 2.0 * k / n;
        if (std::abs(y - 0.5) > s && std::abs(y - 1.5) > s) o.lines.push_back(y * um);
    }
    o.straddles = {{0.5 * um, s * um}, {1.5 * um, s * um}};
    auto m = mesh::generate_swept(x, o);
    INFO((m.has_value() ? std::string() : m.error().message));
    REQUIRE(m.has_value());
    return std::move(*m);
}

struct Solved {
    mesh::GeneratedMesh mesh;
    std::vector<double> psi, n, p;
};

// The contact nodes per face, for the record (x0, x1, z0, z1, y0, y1, inside).
std::string contact_faces(const mesh::GeneratedMesh& m, const char* name) {
    const auto* b = m.mesh.find_boundary(name);
    std::array<int, 7> c{};
    for (const auto v : b->nodes) {
        const mesh::Point& p = m.mesh.points()[static_cast<std::size_t>(v)];
        const int k = p[0] == 0 ? 0 : p[0] == 1 * um ? 1 : p[2] == 0 ? 2 : p[2] == 1 * um ? 3 : p[1] == 0 ? 4 : p[1] == 2 * um ? 5 : 6;
        ++c[static_cast<std::size_t>(k)];
    }
    return std::format("{} {}: {} {} {} {} {} {} {}", name, b->nodes.size(), c[0], c[1], c[2], c[3], c[4], c[5], c[6]);
}

Solved solve_device(mesh::GeneratedMesh m) {
    const std::size_t nn = m.mesh.node_count();
    std::vector<device::RegionId> region(m.node_region.begin(), m.node_region.end());
    std::vector<double> acceptors(nn, 0.0);
    for (std::size_t i = 0; i < nn; ++i) {
        if (region[i] == 1) acceptors[i] = 1e17;
    }
    const auto* body = m.mesh.find_boundary("body");
    const auto* gate = m.mesh.find_boundary("gate");
    REQUIRE(body != nullptr);
    REQUIRE(gate != nullptr);
    device::Contact gc{"gate", device::ContactKind::electrode, gate->nodes};
    gc.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto d = device::Device::create({.mesh = m.mesh,
                                     .temperature_K = 300.0,
                                     .regions = regions(),
                                     .node_region = std::move(region),
                                     .donors = std::vector<double>(nn, 0.0),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"body", device::ContactKind::ohmic, body->nodes}, std::move(gc)},
                                     .interface_certificate = m.certificate});
    INFO((d.has_value() ? std::string() : d.error().message));
    REQUIRE(d.has_value());
    solve::BiasOptions o;
    o.equations = solve::Equations::equilibrium_poisson;
    const auto r = solve::solve_bias(*d, std::vector<double>{0.0, gate_V}, o);
    INFO((r.has_value() ? std::string() : r.error().message));
    REQUIRE(r.has_value());
    return {std::move(m), r->fields.potential_V, r->fields.n_cm3, r->fields.p_cm3};
}

// The tetrahedron holding x and its barycentric weights (the best by the smallest weight: on a face to rounding).
struct Hit {
    std::array<std::size_t, 4> v{};
    std::array<double, 4> w{};
};

Hit locate(const mesh::SimplexMesh& sm, const mesh::Point& x) {
    Hit best;
    double score = -std::numeric_limits<double>::infinity();
    const auto& P = sm.points;
    for (std::size_t t = 0; 4 * t < sm.cells.size(); ++t) {
        std::array<std::size_t, 4> v{};
        for (std::size_t k = 0; k < 4; ++k) v[k] = static_cast<std::size_t>(sm.cells[4 * t + k]);
        bool far = false;
        for (std::size_t c = 0; c < 3 && !far; ++c) {
            const double lo = std::min({P[v[0]][c], P[v[1]][c], P[v[2]][c], P[v[3]][c]});
            const double hi = std::max({P[v[0]][c], P[v[1]][c], P[v[2]][c], P[v[3]][c]});
            far = x[c] < lo - 1e-12 || x[c] > hi + 1e-12;
        }
        if (far) continue;
        const auto det3 = [](const mesh::Point& a, const mesh::Point& b, const mesh::Point& c) {
            return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
        };
        const auto sub = [](const mesh::Point& a, const mesh::Point& b) { return mesh::Point{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
        const mesh::Point e1 = sub(P[v[1]], P[v[0]]), e2 = sub(P[v[2]], P[v[0]]), e3 = sub(P[v[3]], P[v[0]]), r = sub(x, P[v[0]]);
        const double D = det3(e1, e2, e3);
        const double w1 = det3(r, e2, e3) / D, w2 = det3(e1, r, e3) / D, w3 = det3(e1, e2, r) / D;
        const std::array<double, 4> w{1.0 - w1 - w2 - w3, w1, w2, w3};
        const double m = std::min({w[0], w[1], w[2], w[3]});
        if (m > score) {
            score = m;
            best = {v, w};
            if (m >= 0.0) break;
        }
    }
    return best;
}

double at(const Hit& h, const std::vector<double>& f) {
    double s = 0.0;
    for (std::size_t k = 0; k < 4; ++k) s += h.w[k] * f[h.v[k]];
    return s;
}

bool all_silicon(const Hit& h, const mesh::GeneratedMesh& m) {
    return std::ranges::all_of(h.v, [&](std::size_t v) { return m.node_region[v] == 1; });
}

struct Differences {
    double psi = 0.0;  // max |d psi| / (|psi_ref| + 1e-3 V): <= 1e-6 iff |d psi| <= 1e-6 |psi_ref| + 1e-9 V (IC3)
    double psi_abs = 0.0, psi_range = 0.0;  // max |d psi| [V]; max |psi_ref| [V]
    double psi_rms = 0.0;                   // root mean square |d psi| [V]
    mesh::Point worst_at{};                 // where |d psi| is largest [cm]
    double worst_min_weight = 1.0;  // the smallest barycentric weight found (< 0: a point off both meshes to rounding)
    double n = 0.0, p = 0.0;  // max relative
    int density_points = 0;
};

// Points with y in [y0, y1] µm.
Differences compare(const Solved& a, const Solved& ref, int count, double y0 = 0.0, double y1 = 2.0) {
    std::mt19937 rng(64);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    Differences d;
    for (int k = 0; k < count; ++k) {
        // Half the points in the wire's box (where the densities live), half anywhere.
        const mesh::Point x = k % 2 == 0 ? mesh::Point{(0.4 + 0.2 * U(rng)) * um, (y0 + (y1 - y0) * U(rng)) * um, (0.4 + 0.2 * U(rng)) * um}
                                         : mesh::Point{U(rng) * um, (y0 + (y1 - y0) * U(rng)) * um, U(rng) * um};
        const Hit ha = locate(a.mesh.simplices, x), hr = locate(ref.mesh.simplices, x);
        const double pa = at(ha, a.psi), pr = at(hr, ref.psi);
        d.psi = std::max(d.psi, std::abs(pa - pr) / (std::abs(pr) + 1e-3));
        if (std::abs(pa - pr) > d.psi_abs) d.worst_at = x;
        d.psi_abs = std::max(d.psi_abs, std::abs(pa - pr));
        d.psi_rms += (pa - pr) * (pa - pr) / count;
        d.psi_range = std::max(d.psi_range, std::abs(pr));
        d.worst_min_weight = std::min({d.worst_min_weight, *std::ranges::min_element(ha.w), *std::ranges::min_element(hr.w)});
        if (all_silicon(ha, a.mesh) && all_silicon(hr, ref.mesh)) {
            ++d.density_points;
            d.n = std::max(d.n, std::abs(at(ha, a.n) - at(hr, ref.n)) / std::abs(at(hr, ref.n)));
            d.p = std::max(d.p, std::abs(at(ha, a.p) - at(hr, ref.p)) / std::abs(at(hr, ref.p)));
        }
    }
    d.psi_rms = std::sqrt(d.psi_rms);
    return d;
}

std::string line(const char* what, const Differences& d) {
    return std::format("IC3 {}: psi {:.3e} (max {:.3e} V at ({:.3f}, {:.3f}, {:.3f}) um, rms {:.3e} V), n {:.3e}, p {:.3e} ({} density points)",
                       what, d.psi, d.psi_abs, d.worst_at[0] / um, d.worst_at[1] / um, d.worst_at[2] / um, d.psi_rms, d.n, d.p, d.density_points);
}

}  // namespace

TEST_CASE("interface pieces accuracy: a holed gate face's wire against its swept mesh (Unit 64 IC3, recorded)",
          "[.interface_pieces_accuracy][solve][!shouldfail]") {
    const double h = 0.15, s = 0.02;
    const Solved holed = solve_device(polyhedral(h, s));
    const Solved ref = solve_device(swept(h, s));
    const Solved fine = solve_device(swept(h / 2, s / 2));
    const Solved holed2 = solve_device(polyhedral(h / 2, s / 2));
    const Differences d = compare(holed, ref, 200);
    for (const Solved* x : {&holed, &ref}) WARN(std::format("IC3 contacts {} / {}", contact_faces(x->mesh, "gate"), contact_faces(x->mesh, "body")));
    WARN(std::format("IC3: nodes holed {}, holed at h/2 {}, swept {}, swept at h/2 {}; psi up to {:.3e} V", holed.mesh.mesh.node_count(),
                     holed2.mesh.mesh.node_count(), ref.mesh.mesh.node_count(), fine.mesh.mesh.node_count(), d.psi_range));
    WARN(line("holed vs swept", d));
    WARN(line("swept vs swept at h/2", compare(ref, fine, 200)));
    WARN(line("holed vs swept at h/2", compare(holed, fine, 200)));
    WARN(line("holed at h/2 vs swept at h/2", compare(holed2, fine, 200)));
    WARN(line("holed vs holed at h/2", compare(holed, holed2, 200)));
    // Away from the wire's ends (the body contact's edge is singular): y in [0.3, 1.7] µm, the holed faces at 0.5 and 1.5.
    for (const mesh::Point& q : {mesh::Point{0.5 * um, 1.0 * um, 0.5 * um}, mesh::Point{0.5 * um, 0.3 * um, 0.5 * um}, mesh::Point{0.681 * um, 0.364 * um, 0.647 * um}, mesh::Point{0.3 * um, 1.0 * um, 0.5 * um}, mesh::Point{0.1 * um, 1.0 * um, 0.5 * um}}) {
        WARN(std::format("IC3 probe ({:.3f}, {:.3f}, {:.3f}): holed {:.6f} V, h/2 {:.6f}; swept {:.6f}, h/2 {:.6f}", q[0] / um, q[1] / um, q[2] / um, at(locate(holed.mesh.simplices, q), holed.psi), at(locate(holed2.mesh.simplices, q), holed2.psi), at(locate(ref.mesh.simplices, q), ref.psi), at(locate(fine.mesh.simplices, q), fine.psi)));
    }
    WARN(line("mid: holed vs swept", compare(holed, ref, 200, 0.3, 1.7)));
    WARN(line("mid: swept vs swept at h/2", compare(ref, fine, 200, 0.3, 1.7)));
    WARN(line("mid: holed vs holed at h/2", compare(holed, holed2, 200, 0.3, 1.7)));
    WARN(line("mid: holed at h/2 vs swept at h/2", compare(holed2, fine, 200, 0.3, 1.7)));
    CHECK(d.psi <= 1e-6);
    CHECK(d.n <= 1e-6);
    CHECK(d.p <= 1e-6);
}

namespace {

// Unit 64 slice 3 (IP4): SD6's FinFET on a dyadic grid (units of 2^-14 cm times k; tests/mesh/interface_pieces_test.cpp
// dyadic_fin): at k = 1 the gate's bottom cuts the slanted walls exactly coplanar, at k = 1 - 1e-9 it folds them.
// Body contact on z = 0, gate contact on z = 1. IC3 is missed (DECISIONS.md Unit 64, slice 3): the generator places nodes
// differently under any ulp-level change of the layout, folded or not (the uncrossed FinFET under the same scaling
// moves as much), so the two meshes differ and their solutions agree only to discretization. [!shouldfail], 1e-6 kept.
mesh::GeneratedMesh dyadic_finfet(double k, bool crossed) {
    const double u = std::ldexp(1.0, -14) * k;
    const auto B = [&](double x0, double x1, double y0, double y1, double z0, double z1, std::int32_t r) {
        mesh::LayeredSolid s;
        s.kind = mesh::LayeredSolid::Kind::box;
        s.lo = {x0 * u, y0 * u, z0 * u};
        s.hi = {x1 * u, y1 * u, z1 * u};
        s.region = r;
        return s;
    };
    mesh::LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {1 * u, 2 * u, 1 * u}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(B(0, 1, 0, 2, 0, 0.3125, 1));
    g.solids.push_back(B(0, 1, 0.5, 1.5, crossed ? 0.5 : 0.3125, 1, 2));
    mesh::LayeredSolid fin;
    fin.kind = mesh::LayeredSolid::Kind::prism;
    fin.axis = 1;
    fin.from = 0;
    fin.to = 2 * u;
    for (const auto& [z, x] : std::vector<std::pair<double, double>>{{0.3125, 0.3125}, {0.3125, 0.6875}, {0.8125, 0.5625}, {0.8125, 0.4375}}) {
        fin.polygon.push_back({z * u, x * u, 0.0});
    }
    fin.region = 1;
    g.solids.push_back(fin);
    g.patches = {{"body", {0, 0, 0}, {1 * u, 2 * u, 0}}, {"gate", {0, 0, 1 * u}, {1 * u, 2 * u, 1 * u}}};
    const auto c = mesh::layered_partition3d(g);
    REQUIRE(c.has_value());
    mesh::MeshOptions3D o;
    o.max_size = 0.15 * u;
    o.interface_resolution = u / 64;
    auto m = mesh::generate_mesh(*c, o);
    INFO((m.has_value() ? std::string() : m.error().message));
    REQUIRE(m.has_value());
    return std::move(*m);
}

// As compare(), on the dyadic FinFET: a point x of the reference is k x on `a`.
Differences compare_scaled(const Solved& a, double k, const Solved& ref, int count) {
    const double u = std::ldexp(1.0, -14);
    std::mt19937 rng(64);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    Differences d;
    for (int i = 0; i < count; ++i) {
        // Half the points in the fin's box, half anywhere.
        const mesh::Point x = i % 2 == 0 ? mesh::Point{(0.3125 + 0.375 * U(rng)) * u, 2 * U(rng) * u, (0.3125 + 0.5 * U(rng)) * u}
                                         : mesh::Point{U(rng) * u, 2 * U(rng) * u, U(rng) * u};
        const Hit ha = locate(a.mesh.simplices, {k * x[0], k * x[1], k * x[2]}), hr = locate(ref.mesh.simplices, x);
        const double pa = at(ha, a.psi), pr = at(hr, ref.psi);
        d.psi = std::max(d.psi, std::abs(pa - pr) / (std::abs(pr) + 1e-3));
        if (std::abs(pa - pr) > d.psi_abs) d.worst_at = {x[0] / u * um, x[1] / u * um, x[2] / u * um};
        d.psi_abs = std::max(d.psi_abs, std::abs(pa - pr));
        d.psi_rms += (pa - pr) * (pa - pr) / count;
        d.psi_range = std::max(d.psi_range, std::abs(pr));
        d.worst_min_weight = std::min({d.worst_min_weight, *std::ranges::min_element(ha.w), *std::ranges::min_element(hr.w)});
        if (all_silicon(ha, a.mesh) && all_silicon(hr, ref.mesh)) {
            ++d.density_points;
            d.n = std::max(d.n, std::abs(at(ha, a.n) - at(hr, ref.n)) / std::abs(at(hr, ref.n)));
            d.p = std::max(d.p, std::abs(at(ha, a.p) - at(hr, ref.p)) / std::abs(at(hr, ref.p)));
        }
    }
    d.psi_rms = std::sqrt(d.psi_rms);
    return d;
}

}  // namespace

TEST_CASE("interface pieces accuracy: the folded FinFET against the same FinFET cut exactly coplanar (Unit 64 IC3, recorded)",
          "[.interface_pieces_accuracy][solve][!shouldfail]") {
    const double k = 1.0 - 1e-9;
    const Solved exact = solve_device(dyadic_finfet(1.0, true));
    const Solved bent = solve_device(dyadic_finfet(k, true));
    // The floor: the uncrossed FinFET (no fold at either k) under the same scaling.
    const Solved plain = solve_device(dyadic_finfet(1.0, false));
    const Solved plain_k = solve_device(dyadic_finfet(k, false));
    const Differences d = compare_scaled(bent, k, exact, 200);
    WARN(std::format("IC3 fold: nodes exact {}, folded {}; uncrossed {}, scaled {}; psi up to {:.3e} V (coordinates in units of 0.61 um)", exact.mesh.mesh.node_count(),
                     bent.mesh.mesh.node_count(), plain.mesh.mesh.node_count(), plain_k.mesh.mesh.node_count(), d.psi_range));
    WARN(line("folded vs exact", d));
    WARN(line("uncrossed scaled vs uncrossed (the floor)", compare_scaled(plain_k, k, plain, 200)));
    CHECK(d.psi <= 1e-6);
    CHECK(d.n <= 1e-6);
    CHECK(d.p <= 1e-6);
}
