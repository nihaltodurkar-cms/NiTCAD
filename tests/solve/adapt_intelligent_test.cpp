// AC14, AC15 and AC16 (Unit 17, slice 5; DECISIONS.md Unit 17 ID1-ID6, AD12, AR5-AR7): the a priori
// prediction against the region the basic loop finally refines (recall, precision, the false-positive
// DOFs and their cost, the mesh overhead, the passes) on the legacy diode, a MOS-C and the 2D MOSFET
// on a tensor mesh (AR6); the benefit and cost estimates against the realized ones; the allocation's
// DOFs against basic's in tolerance and budget modes. Release only (hidden tag; CMakeLists.txt).
// AC16 on AC7's well MOSFET is in adapt_goal_test.cpp.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/locate.hpp"
#include "NiTCAD/mesh/size_field.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/adapt.hpp"
#include "NiTCAD/solve/adapt_policy.hpp"
#include "NiTCAD/solve/bias.hpp"

using namespace NiTCAD;

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();

double seconds_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

// --- devices -----------------------------------------------------------------------------------

device::Device legacy_diode(const std::vector<double>& x) {
    auto m = mesh::make_tensor_grid(x);
    REQUIRE(m.has_value());
    auto cells = mesh::TensorCells::create(*m, x);
    REQUIRE(cells.has_value());
    const std::size_t n = m->node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m->points()[i][0] < 1e-4 ? acceptors[i] : donors[i]) = 1e17;
    auto anode = m->find_boundary("x_min")->nodes;
    auto cathode = m->find_boundary("x_max")->nodes;
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::legacy_silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
                                     .cells = std::move(*cells)});
    REQUIRE(d.has_value());
    return std::move(*d);
}

std::vector<double> uniform(double L, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = i == nodes - 1 ? L : L * i / (nodes - 1);
    return a;
}

std::string axis_text(const std::vector<double>& a) {
    std::string s;
    for (std::size_t i = 0; i < a.size(); ++i) s += (i ? ", " : "") + std::format("{:.17g}cm", a[i]);
    return s;
}

// A coarse MOS-C (adapt_indicators_test.cpp's device): p 1e17 under 1.5 nm of oxide, the gate at 1 V.
deck::Plan mosc_plan(const std::vector<double>& x, const std::vector<double>& y) {
    const std::string text = "[device]\ntemperature = 300\n[mesh]\ntype = tensor\nx = " + axis_text(x) +
                             "\ny = " + axis_text(y) + R"(
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e17
[contact gate]
kind = gate
patch = x_min
oxide_thickness = 1.5nm
electrode = n_poly
[contact substrate]
kind = ohmic
patch = x_max
[sweep]
contact gate = 0.5, 1.0
[solver]
equations = equilibrium_poisson
)";
    auto d = deck::parse(text, "mosc.deck", {});
    if (!d) FAIL(d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

// The MOSFET of unstructured_test.cpp (Unit 12's; Unit 16's deck shape): Lg 600 nm, Lsd 300 nm,
// depth 200 nm, N_A 1e17, n+ source and drain (1e19 Gaussian, sigma 50 nm), t_ox 5 nm.
constexpr double Lg = 6e-5, Lsd = 3e-5, depth = 2e-5, t_ox = 5e-7;
constexpr double h0 = 2e-6;  // the start's size, 20 nm

std::vector<device::DopingProfile> mosfet_doping() {
    const double lateral = std::sqrt(2.0) * 1e-6;
    device::DopingProfile source{.dopant = device::Dopant::donor, .peak_cm3 = 1e19,
                                 .shape = device::DepthShape::gaussian, .axis = 1, .center = 0.0,
                                 .width = 5e-6, .lateral = lateral};
    source.hi[0] = Lsd;
    device::DopingProfile drain = source;
    drain.hi[0] = INFINITY;
    drain.lo[0] = Lsd + Lg;
    return {source, drain, {.dopant = device::Dopant::acceptor, .peak_cm3 = 1e17}};
}

std::expected<solve::AdaptableMesh, base::Error> mosfet_mesh(const solve::MeshRefinement& r) {
    const double L = 2 * Lsd + Lg;
    const auto rect = [](double x0, double y0, double x1, double y1) {
        return mesh::RegionPolygon{{{x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}}, 0};
    };
    mesh::Geometry2D g{{rect(0, 0, L, depth)},
                       {{"source", {0, 0, 0}, {Lsd, 0, 0}}, {"channel", {Lsd, 0, 0}, {Lsd + Lg, 0, 0}},
                        {"drain", {Lsd + Lg, 0, 0}, {L, 0, 0}}, {"body", {0, depth, 0}, {L, depth, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = h0;
    o.size_field = r.size_field;
    auto gm = mesh::generate_mesh(g, o);
    if (!gm) return std::unexpected(gm.error());
    const mesh::Mesh& m = gm->mesh;
    const auto nodes = [&](const char* name) { return m.find_boundary(name)->nodes; };
    std::vector<mesh::NodeId> source = nodes("source"), drain = nodes("drain"), gate;
    for (const mesh::NodeId v : nodes("channel")) {
        if (!std::ranges::binary_search(source, v) && !std::ranges::binary_search(drain, v)) gate.push_back(v);
    }
    const std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    const std::vector<device::RegionId> region(m.node_count(), 0);
    auto doping = device::sample(mosfet_doping(), m, region, regions);
    if (!doping) return std::unexpected(doping.error());
    device::Contact gc{"gate", device::ContactKind::gate, std::move(gate)};
    gc.gate = {.boundary = "channel", .oxide_thickness_cm = t_ox};
    auto d = device::Device::create({.mesh = m,
                                     .temperature_K = 300.0,
                                     .regions = regions,
                                     .node_region = region,
                                     .donors = doping->donors,
                                     .acceptors = doping->acceptors,
                                     .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                                                  {"drain", device::ContactKind::ohmic, std::move(drain)},
                                                  {"body", device::ContactKind::ohmic, nodes("body")},
                                                  std::move(gc)}});
    if (!d) return std::unexpected(d.error());
    return solve::AdaptableMesh{.device = std::move(*d), .kind = solve::AdaptMeshKind::generated2d,
                                .generated = std::move(*gm), .max_size = h0};
}

// The same MOSFET on a tensor mesh (AR6: AC14's MOSFET; the generated mesh's isotropic refinement
// does not converge within a gate's budget): from 25 x 11 uniform nodes (50 nm by 20 nm), the
// contact edges on lines.
std::expected<solve::AdaptableMesh, base::Error> tensor_mosfet_mesh(const solve::MeshRefinement& r) {
    const std::vector<double> x = r.axes[0].empty() ? uniform(2 * Lsd + Lg, 25) : r.axes[0];
    const std::vector<double> y = r.axes[1].empty() ? uniform(depth, 11) : r.axes[1];
    auto m = mesh::make_tensor_grid(x, y);
    if (!m) return std::unexpected(m.error());
    auto cells = mesh::TensorCells::create(*m, x, y);
    if (!cells) return std::unexpected(cells.error());
    std::vector<mesh::NodeId> source, drain, gate;
    for (const mesh::NodeId v : m->find_boundary("y_min")->nodes) {
        const double px = m->points()[static_cast<std::size_t>(v)][0];
        (px <= Lsd ? source : (px >= Lsd + Lg ? drain : gate)).push_back(v);
    }
    auto body = m->find_boundary("y_max")->nodes;
    const std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    const std::vector<device::RegionId> region(m->node_count(), 0);
    auto doping = device::sample(mosfet_doping(), *m, region, regions);
    if (!doping) return std::unexpected(doping.error());
    device::Contact gc{"gate", device::ContactKind::gate, std::move(gate)};
    gc.gate = {.boundary = "y_min", .oxide_thickness_cm = t_ox};
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = 300.0,
                                     .regions = regions,
                                     .node_region = region,
                                     .donors = doping->donors,
                                     .acceptors = doping->acceptors,
                                     .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                                                  {"drain", device::ContactKind::ohmic, std::move(drain)},
                                                  {"body", device::ContactKind::ohmic, std::move(body)},
                                                  std::move(gc)},
                                     .cells = std::move(*cells)});
    if (!d) return std::unexpected(d.error());
    return solve::AdaptableMesh{.device = std::move(*d), .kind = solve::AdaptMeshKind::tensor};
}

// --- sizes and metrics -------------------------------------------------------------------------

// The size of a tensor axis at x: the interval containing it; at a line, the larger neighbour.
double size_at(std::span<const double> lines, double x) {
    const auto it = std::ranges::lower_bound(lines, x);
    const auto k = static_cast<std::size_t>(it - lines.begin());
    if (k < lines.size() && std::abs(lines[k] - x) <= 1e-12 * (std::abs(x) + 1e-12)) {
        double s = 0.0;
        if (k > 0) s = std::max(s, lines[k] - lines[k - 1]);
        if (k + 1 < lines.size()) s = std::max(s, lines[k + 1] - lines[k]);
        return s;
    }
    if (k == 0) return lines[1] - lines[0];
    if (k >= lines.size()) return lines.back() - lines[lines.size() - 2];
    return lines[k] - lines[k - 1];
}

// The local size of a mesh at a point: per axis for tensor meshes, isotropic for generated ones.
class Sizes {
public:
    explicit Sizes(const solve::AdaptableMesh& m) : tensor_(m.kind == solve::AdaptMeshKind::tensor) {
        if (tensor_) {
            const mesh::TensorCells& c = *m.device.cells();
            for (int a = 0; a < c.dimension(); ++a) {
                const auto ax = c.axis(a);
                axes_.emplace_back(ax.begin(), ax.end());
            }
        } else {
            local_ = mesh::local_sizes(m.generated->simplices);
            locator_ = std::make_shared<mesh::PointLocator>(*mesh::PointLocator::create(m.generated->simplices));
        }
    }
    [[nodiscard]] bool tensor() const { return tensor_; }
    [[nodiscard]] std::size_t dimension() const { return tensor_ ? axes_.size() : 1; }  // generated: one, isotropic
    // Along axis a (tensor), or the isotropic size (generated; a ignored).
    [[nodiscard]] double at(const mesh::Point& p, std::size_t a) const {
        if (tensor_) return size_at(axes_[a], p[a]);
        const mesh::Interpolant w = locator_->locate(p);
        double s = 0.0;
        for (int j = 0; j < w.count; ++j) s += w.weights[static_cast<std::size_t>(j)] * local_[static_cast<std::size_t>(w.nodes[static_cast<std::size_t>(j)])];
        return s;
    }

private:
    bool tensor_;
    std::vector<std::vector<double>> axes_;
    std::vector<double> local_;
    std::shared_ptr<mesh::PointLocator> locator_;
};

// Refined at p along axis a: smaller than the start (tensor), or below 0.9 of it (generated, one
// isotropic axis: the generator's sizes scatter around the field).
bool refined(const Sizes& start, const Sizes& fin, const mesh::Point& p, std::size_t a) {
    const double s = start.at(p, a), f = fin.at(p, a);
    return start.tensor() ? f < s * (1 - 1e-9) : f < 0.9 * s;
}

// Predicted to be refined at p along axis a: the target below the start's size there, a being the
// normal's axis (any axis for an isotropic target; the one axis of a generated mesh).
bool predicted(const Sizes& start, double target, const std::array<double, 3>& n, const mesh::Point& p, std::size_t a) {
    if (!std::isfinite(target)) return false;
    if (!start.tensor()) return target < 0.9 * start.at(p, 0);
    std::size_t axis = 0;
    double best = 0.0;
    for (std::size_t b = 0; b < start.dimension(); ++b) {
        if (std::abs(n[b]) > best) {
            best = std::abs(n[b]);
            axis = b;
        }
    }
    if (best > 0.0 && axis != a) return false;
    return target < start.at(p, a);
}

struct Gate14 {
    double recall = 0.0, precision = 0.0;
    std::size_t basic_refined = 0, predicted_nodes = 0, false_positive_nodes = 0;
    std::size_t false_positive_dofs = 0, false_positive_nnz = 0, intelligent_dofs = 0, basic_dofs = 0;
    int passes_basic = 0, passes_intelligent = 0;
    double seconds_basic = 0.0, seconds_intelligent = 0.0;
    double final_pass_basic = 0.0, final_pass_intelligent = 0.0;
};

std::size_t dofs(const device::Device& d) {
    std::size_t n = 0;
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) n += d.is_insulator(static_cast<mesh::NodeId>(i)) ? 1 : 3;
    return n;
}

Gate14 gate14(const char* name, const solve::DeviceBuilder& build, const std::vector<std::vector<double>>& pts,
              const solve::BiasOptions& bias, solve::AdaptOptions o) {
    Gate14 g;
    auto start = build({});
    REQUIRE(start.has_value());
    o.policy = solve::Policy::basic;
    auto t = std::chrono::steady_clock::now();
    auto b = solve::adapt(build, pts, bias, o);
    g.seconds_basic = seconds_since(t);
    if (!b) FAIL(name << " basic: " << b.error().message);
    o.policy = solve::Policy::intelligent;
    t = std::chrono::steady_clock::now();
    auto r = solve::adapt(build, pts, bias, o);
    g.seconds_intelligent = seconds_since(t);
    if (!r) FAIL(name << " intelligent: " << r.error().message);
    const Sizes s0(*start), sb(*b->mesh), si(*r->mesh);
    // Recall, per node and axis: of the basic loop's refinements, those the a priori field (on that
    // mesh) predicts along the same axis.
    const device::Device& db = b->mesh->device;
    auto pb = solve::predict_sizes(db, pts, bias);
    REQUIRE(pb.has_value());
    std::size_t covered = 0;
    for (std::size_t i = 0; i < db.mesh().node_count(); ++i) {
        const auto& p = db.mesh().points()[i];
        for (std::size_t a = 0; a < s0.dimension(); ++a) {
            if (!refined(s0, sb, p, a)) continue;
            ++g.basic_refined;
            covered += predicted(s0, pb->targets.size[i], pb->targets.normal[i], p, a) ? 1 : 0;
        }
    }
    g.recall = g.basic_refined > 0 ? static_cast<double>(covered) / static_cast<double>(g.basic_refined) : 1.0;
    // Precision and false positives on the intelligent final mesh, per node and axis; a node with a
    // predicted refinement the basic loop did not make there is a false positive.
    const device::Device& di = r->mesh->device;
    auto pi = solve::predict_sizes(di, pts, bias);
    REQUIRE(pi.has_value());
    std::vector<char> fp(di.mesh().node_count(), 0);
    std::size_t agreed = 0;
    for (std::size_t i = 0; i < di.mesh().node_count(); ++i) {
        const auto& p = di.mesh().points()[i];
        for (std::size_t a = 0; a < s0.dimension(); ++a) {
            if (!predicted(s0, pi->targets.size[i], pi->targets.normal[i], p, a)) continue;
            ++g.predicted_nodes;
            if (refined(s0, sb, p, a)) {
                ++agreed;
            } else if (fp[i] == 0) {
                fp[i] = 1;
                ++g.false_positive_nodes;
                g.false_positive_dofs += di.is_insulator(static_cast<mesh::NodeId>(i)) ? 1 : 3;
            }
        }
    }
    for (const mesh::Edge& e : di.mesh().edges()) {
        if (fp[static_cast<std::size_t>(e.first)] != 0 || fp[static_cast<std::size_t>(e.second)] != 0) g.false_positive_nnz += 2 * 9;
    }
    g.false_positive_nnz += 9 * g.false_positive_nodes;
    g.precision = g.predicted_nodes > 0 ? static_cast<double>(agreed) / static_cast<double>(g.predicted_nodes) : 1.0;
    g.intelligent_dofs = dofs(di);
    g.basic_dofs = dofs(db);
    g.passes_basic = static_cast<int>(b->passes.size());
    g.passes_intelligent = static_cast<int>(r->passes.size());
    g.final_pass_basic = b->passes.back().seconds_solve;
    g.final_pass_intelligent = r->passes.back().seconds_solve;
    for (std::size_t k = 0; k < r->passes.size(); ++k) {
        const auto& p = r->passes[k];
        WARN("AC14 " << name << " intelligent pass " << k << ": " << solve::to_string(p.source) << ", nodes " << p.nodes
                     << ", relative estimate " << p.relative_estimate << (p.miss ? " (miss)" : ""));
    }
    for (std::size_t k = 0; k < b->passes.size(); ++k) {
        WARN("AC14 " << name << " basic pass " << k << ": nodes " << b->passes[k].nodes << ", relative estimate "
                     << b->passes[k].relative_estimate);
    }
    WARN("AC14 " << name << ": recall " << g.recall << " (" << covered << " of " << g.basic_refined << "), precision "
                 << g.precision << " (" << agreed << " of " << g.predicted_nodes << "), false-positive DOFs "
                 << g.false_positive_dofs << " of " << g.intelligent_dofs << " (" << g.false_positive_nnz
                 << " nonzeros), DOFs intelligent " << g.intelligent_dofs << " basic " << g.basic_dofs << ", passes "
                 << g.passes_intelligent << " vs " << g.passes_basic << ", stop " << solve::to_string(r->stop) << " vs "
                 << solve::to_string(b->stop) << ", certificates " << solve::to_string(r->certificates[0].status) << " vs "
                 << solve::to_string(b->certificates[0].status) << ", final pass " << g.final_pass_intelligent << " s vs "
                 << g.final_pass_basic << " s, total " << g.seconds_intelligent << " s vs " << g.seconds_basic << " s;"
                 << " prediction: " << r->junction_points << " junction and " << r->interface_points
                 << " interface points, " << r->prediction_rounds << " rounds; fallback " << solve::to_string(r->fallback)
                 << ", misses " << r->misses);
    CHECK(b->stop == solve::AdaptStop::converged);
    CHECK(r->stop == solve::AdaptStop::converged);
    for (const auto& c : r->certificates) CHECK(c.status == solve::CertificateStatus::certified);
    return g;
}

solve::BiasOptions legacy_options() {
    solve::BiasOptions o{.models = {.hole_mobility = assemble::HoleMobility::caughey_thomas}};
    o.step_control.max_cuts = 20;
    o.weighted_current = true;
    return o;
}

}  // namespace

// Fails as frozen (AR6, owner ruling pending; re-measured after slice 5b's ID4-ID6): measured
// 2026-11-12, recall per node and axis 0.87 / 0.83 / 0.51 (the tensor MOSFET's prediction is
// anisotropic, the basic loop refines both axes), DOFs 1116 vs 540 and 504 vs 216 (diode, MOS-C),
// passes 11 vs 10, 2 vs 2, 5 vs 6 (saved 0).
TEST_CASE("intelligent: the a priori prediction covers the basic loop's refinement (AC14)", "[.adapt_intelligent][!shouldfail]") {
    // The diode (AC6's run).
    const std::vector<std::vector<double>> diode_points{{-1.0, 0.0}, {-0.5, 0.0}, {0.0, 0.0}, {0.3, 0.0}, {0.6, 0.0}};
    solve::AdaptOptions od;
    od.goals = {{0}};
    od.tolerance = 1e-3;
    od.at = {0, 4};
    od.max_passes = 14;
    const solve::DeviceBuilder diode = [](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        const std::vector<double> x = r.axes[0].empty() ? uniform(2e-4, 21) : r.axes[0];
        return solve::AdaptableMesh{.device = legacy_diode(x), .kind = solve::AdaptMeshKind::tensor};
    };
    const Gate14 a = gate14("diode", diode, diode_points, legacy_options(), od);

    // The MOS-C: the substrate's current (zero) is the only goal; the basic loop refines once.
    const deck::Plan mp = mosc_plan(uniform(1e-4, 21), {0.0, 0.5e-4, 1e-4});
    solve::AdaptOptions om;
    om.goals = {{1}};
    om.tolerance = 1e-2;
    om.max_passes = 8;
    const solve::DeviceBuilder mosc = [](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        const std::vector<double> x = r.axes[0].empty() ? uniform(1e-4, 21) : r.axes[0];
        const std::vector<double> y = r.axes[1].empty() ? std::vector<double>{0.0, 0.5e-4, 1e-4} : r.axes[1];
        deck::Plan p = mosc_plan(x, y);
        return solve::AdaptableMesh{.device = std::move(p.device), .kind = solve::AdaptMeshKind::tensor};
    };
    const Gate14 b = gate14("mos-c", mosc, mp.points, mp.options, om);

    // The 2D MOSFET on a tensor mesh (AR6), from 25 x 11 nodes: I_D at Vg 0.5 and 1 V, Vd 0.05 V.
    const std::vector<std::vector<double>> mosfet_points{{0.0, 0.05, 0.0, 0.0}, {0.0, 0.05, 0.0, 0.5}, {0.0, 0.05, 0.0, 1.0}};
    solve::BiasOptions mb;
    mb.step_control.max_cuts = 20;
    solve::AdaptOptions of;
    of.goals = {{1}};
    of.tolerance = 1e-2;
    of.at = {1, 2};
    of.max_passes = 10;
    const Gate14 c = gate14("mosfet", tensor_mosfet_mesh, mosfet_points, mb, of);

    for (const Gate14* g : {&a, &b, &c}) {
        CHECK(g->recall >= 0.8);
        CHECK(static_cast<double>(g->false_positive_dofs) <= 0.3 * static_cast<double>(g->intelligent_dofs));
        CHECK(static_cast<double>(g->intelligent_dofs) <= 1.1 * static_cast<double>(g->basic_dofs));
    }
    const int saved = (a.passes_basic - a.passes_intelligent) + (b.passes_basic - b.passes_intelligent) +
                      (c.passes_basic - c.passes_intelligent);
    WARN("AC14 passes saved, summed over the three: " << saved);
    CHECK(saved >= 1);
}

TEST_CASE("intelligent: probe the 2D MOSFET's passes", "[.adapt_probe]") {
    const std::vector<std::vector<double>> pts{{0.0, 0.05, 0.0, 0.0}, {0.0, 0.05, 0.0, 0.5}, {0.0, 0.05, 0.0, 1.0}};
    solve::BiasOptions mb;
    mb.step_control.max_cuts = 20;
    solve::AdaptOptions o;
    o.goals = {{1}};
    o.tolerance = 1e-2;
    o.at = {1, 2};
    o.max_passes = 4;  // measured to 12 passes: 84715 nodes, I_D still moving, 47 minutes (AR6)
    o.max_nodes = 60000;
    for (const solve::Policy policy : {solve::Policy::intelligent, solve::Policy::basic}) {
        o.policy = policy;
        const auto t = std::chrono::steady_clock::now();
        auto r = solve::adapt(mosfet_mesh, pts, mb, o);
        if (!r) FAIL(r.error().message);
        for (std::size_t k = 0; k < r->passes.size(); ++k) {
            const auto& p = r->passes[k];
            WARN("probe " << solve::to_string(policy) << " pass " << k << ": " << solve::to_string(p.source) << ", nodes "
                          << p.nodes << ", I_D " << p.goal[0] << " " << p.goal[1] << ", relative " << p.relative_estimate
                          << ", mesh " << p.seconds_mesh << " s, solve " << p.seconds_solve << " s, marked "
                          << p.marking.count);
        }
        WARN("probe " << solve::to_string(policy) << ": stop " << solve::to_string(r->stop) << ", rounds "
                      << r->prediction_rounds << ", " << seconds_since(t) << " s");
    }
}

namespace {

// --- AC15, AC16 ------------------------------------------------------------------------------------

std::vector<double> ranks(const std::vector<double>& v) {
    std::vector<std::size_t> o(v.size());
    for (std::size_t i = 0; i < o.size(); ++i) o[i] = i;
    std::ranges::stable_sort(o, [&](std::size_t a, std::size_t b) { return v[a] < v[b]; });
    std::vector<double> r(v.size());
    for (std::size_t i = 0; i < o.size();) {
        std::size_t j = i;
        while (j + 1 < o.size() && v[o[j + 1]] == v[o[i]]) ++j;
        for (std::size_t k = i; k <= j; ++k) r[o[k]] = 0.5 * static_cast<double>(i + j);  // ties: the mean rank
        i = j + 1;
    }
    return r;
}

double spearman(const std::vector<double>& x, const std::vector<double>& y) {
    const std::vector<double> rx = ranks(x), ry = ranks(y);
    const double n = static_cast<double>(x.size());
    double mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        mx += rx[i] / n;
        my += ry[i] / n;
    }
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        sxy += (rx[i] - mx) * (ry[i] - my);
        sxx += (rx[i] - mx) * (rx[i] - mx);
        syy += (ry[i] - my) * (ry[i] - my);
    }
    return sxy / std::sqrt(sxx * syy);
}

// AC15's measures of a run: the predicted against the realized reductions of its listed candidates
// (rank correlation), and from the third pass on the passes' times against the cost model's and the
// node-only model's predictions (each fitted on the earlier passes).
struct Estimates {
    double rank = std::numeric_limits<double>::quiet_NaN();
    std::size_t pairs = 0;
    std::size_t timed = 0;
    double worst = 0.0;        // the cost model's largest relative error
    double mare = 0.0;         // its mean absolute relative error
    double mare_nodes = 0.0;   // the node-only model's
};

Estimates estimates_of(const solve::AdaptResult& r, const char* name) {
    Estimates e;
    std::vector<double> predicted, realized;
    for (const auto& p : r.passes) {
        for (const auto& c : p.candidates) {
            if (!std::isfinite(c.realized)) continue;
            predicted.push_back(c.benefit);
            realized.push_back(c.realized);
        }
    }
    e.pairs = predicted.size();
    if (e.pairs >= 3) e.rank = spearman(predicted, realized);
    for (std::size_t k = 2; k < r.passes.size(); ++k) {
        const auto& p = r.passes[k];
        const double t = p.seconds_mesh + p.seconds_transfer + p.seconds_solve;
        const double a = std::abs(p.predicted_seconds / t - 1.0), b = std::abs(p.predicted_seconds_nodes / t - 1.0);
        WARN("AC15 " << name << " pass " << k << ": " << t << " s, cost model " << p.predicted_seconds << " ("
                     << a << "), node-only " << p.predicted_seconds_nodes << " (" << b << ")");
        ++e.timed;
        e.worst = std::max(e.worst, a);
        e.mare += a;
        e.mare_nodes += b;
    }
    if (e.timed > 0) {
        e.mare /= static_cast<double>(e.timed);
        e.mare_nodes /= static_cast<double>(e.timed);
    }
    WARN("AC15 " << name << ": rank correlation " << e.rank << " over " << e.pairs << " candidates; cost model from the "
                 << "third pass worst " << e.worst << ", mean " << e.mare << " vs node-only " << e.mare_nodes << " over "
                 << e.timed << " passes");
    return e;
}

const std::vector<std::vector<double>> ac16_points{{0.0, 0.05, 0.0, 0.0}, {0.0, 0.05, 0.0, 0.5}, {0.0, 0.05, 0.0, 1.0}};

solve::AdaptOptions ac16_options() {
    solve::AdaptOptions o;
    o.goals = {{1}};
    o.tolerance = 1e-2;
    o.at = {1, 2};
    o.max_passes = 10;
    o.weighting = solve::Weighting::adjoint;
    return o;
}

solve::BiasOptions ac16_bias() {
    solve::BiasOptions b;
    b.step_control.max_cuts = 20;
    return b;
}

// I_D on the tensor MOSFET's start refined uniformly (every interval split k times).
std::vector<double> mosfet_reference(int k) {
    solve::MeshRefinement r;
    r.axes[0] = uniform(2 * Lsd + Lg, 24 * k + 1);
    r.axes[1] = uniform(depth, 10 * k + 1);
    auto m = tensor_mosfet_mesh(r);
    REQUIRE(m.has_value());
    auto s = solve::sweep_bias(m->device, ac16_points, ac16_bias());
    REQUIRE(s.has_value());
    REQUIRE(!s->stopped);
    return {s->points[1].terminal_current[1], s->points[2].terminal_current[1]};
}

double mosfet_error(const solve::AdaptResult& r, const std::vector<double>& ref) {
    return std::max(std::abs(r.sweep.points[1].terminal_current[1] - ref[0]) / std::abs(ref[0]),
                    std::abs(r.sweep.points[2].terminal_current[1] - ref[1]) / std::abs(ref[1]));
}

}  // namespace

TEST_CASE("intelligent: the allocation needs fewer DOFs than basic, and does better on a budget (AC16)",
          "[.adapt_intelligent]") {
    // The 2D MOSFET on the tensor mesh (AR6), adjoint-weighted (AR7), against its start refined 8
    // times uniformly (4 times reported as the reference's spread).
    const std::vector<double> ref = mosfet_reference(8), ref4 = mosfet_reference(4);
    WARN("AC16 reference I_D (8x): " << ref[0] << ", " << ref[1] << "; 4x " << ref4[0] << ", " << ref4[1]);
    solve::AdaptOptions o = ac16_options();
    o.policy = solve::Policy::basic;
    auto b = solve::adapt(tensor_mosfet_mesh, ac16_points, ac16_bias(), o);
    if (!b) FAIL(b.error().message);
    o.policy = solve::Policy::intelligent;
    auto i = solve::adapt(tensor_mosfet_mesh, ac16_points, ac16_bias(), o);
    if (!i) FAIL(i.error().message);
    const std::size_t db = dofs(b->mesh->device), di = dofs(i->mesh->device);
    const double eb = mosfet_error(*b, ref), ei = mosfet_error(*i, ref);
    WARN("AC16 tolerance mode: DOFs intelligent " << di << " basic " << db << " (ratio "
                                                  << static_cast<double>(di) / static_cast<double>(db) << "), errors "
                                                  << ei << " vs " << eb << ", passes " << i->passes.size() << " vs "
                                                  << b->passes.size() << ", misses " << i->misses << ", fallback "
                                                  << solve::to_string(i->fallback));
    // The same verified accuracy: both certified at the tolerance. Their errors against the 8x
    // reference are reported, not gated here: 4x and 8x differ by 2.8% at Vg 0.5 V, so the reference
    // is not converged to the tolerance (AR7; certification against references is AC20's, slice 7).
    CHECK(b->stop == solve::AdaptStop::converged);
    CHECK(i->stop == solve::AdaptStop::converged);
    for (const auto& c : b->certificates) CHECK(c.status == solve::CertificateStatus::certified);
    for (const auto& c : i->certificates) CHECK(c.status == solve::CertificateStatus::certified);
    CHECK(static_cast<double>(di) <= 0.8 * static_cast<double>(db));
    // Budget mode at half basic's final DOFs; basic at the same budget (max_nodes, all nodes silicon).
    o.intelligent.dof_budget = db / 2;
    auto ib = solve::adapt(tensor_mosfet_mesh, ac16_points, ac16_bias(), o);
    if (!ib) FAIL(ib.error().message);
    solve::AdaptOptions ob = ac16_options();
    ob.max_nodes = db / 2 / 3;
    auto bb = solve::adapt(tensor_mosfet_mesh, ac16_points, ac16_bias(), ob);
    if (!bb) FAIL(bb.error().message);
    const double eib = mosfet_error(*ib, ref), ebb = mosfet_error(*bb, ref);
    WARN("AC16 budget " << db / 2 << " DOFs: intelligent " << dofs(ib->mesh->device) << " DOFs, error " << eib << " ("
                        << solve::to_string(ib->stop) << "); basic " << dofs(bb->mesh->device) << " DOFs, error " << ebb
                        << " (" << solve::to_string(bb->stop) << ")");
    CHECK(dofs(ib->mesh->device) <= db / 2);
    CHECK(eib <= ebb);
}

// Fails as frozen (AR7, owner ruling pending), measured 2026-11-13: rank correlation 0.97 (MOSFET,
// 14 candidates) and 0.997 (diode, 37); the cost model within 30% from the third pass (worst 5% and
// 15%) but not better than node-only (mean 3.3% vs 2.7%, 8.0% vs 7.6%: both small, within timing
// noise); the diode's orders 1.10, 1.09 and 0.81, and hole density's 1.64 from 2 samples.
TEST_CASE("intelligent: the benefit and cost estimates track the realized ones (AC15)",
          "[.adapt_intelligent][!shouldfail]") {
    // AC5's run is the generated MOSFET's, which does not converge in a gate's time (AR6): the tensor
    // MOSFET stands in; AC6's diode, adjoint-weighted; AC8 (the FinFET) is slice 7's.
    solve::AdaptOptions o = ac16_options();
    o.policy = solve::Policy::intelligent;
    auto m = solve::adapt(tensor_mosfet_mesh, ac16_points, ac16_bias(), o);
    if (!m) FAIL(m.error().message);
    const Estimates em = estimates_of(*m, "mosfet");
    solve::AdaptOptions od;
    od.goals = {{0}};
    od.tolerance = 1e-3;
    od.at = {0, 4};
    od.max_passes = 14;
    od.weighting = solve::Weighting::adjoint;
    od.policy = solve::Policy::intelligent;
    const std::vector<std::vector<double>> diode_points{{-1.0, 0.0}, {-0.5, 0.0}, {0.0, 0.0}, {0.3, 0.0}, {0.6, 0.0}};
    const solve::DeviceBuilder diode = [](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        const std::vector<double> x = r.axes[0].empty() ? uniform(2e-4, 21) : r.axes[0];
        return solve::AdaptableMesh{.device = legacy_diode(x), .kind = solve::AdaptMeshKind::tensor};
    };
    auto d = solve::adapt(diode, diode_points, legacy_options(), od);
    if (!d) FAIL(d.error().message);
    const Estimates ed = estimates_of(*d, "diode");
    // Rank correlation over the MOSFET's refined candidates.
    REQUIRE(em.pairs >= 5);
    CHECK(em.rank >= 0.6);
    // The diode's fitted orders (abrupt: first order expected).
    for (const auto& f : d->passes.back().orders) {
        WARN("AC15 diode order: kind " << static_cast<int>(f.kind) << ", region " << f.region << ": " << f.order << " from "
                                       << f.samples << " samples");
        CHECK(f.order >= 0.8);
        CHECK(f.order <= 1.3);
    }
    // The cost model within 30% from the third pass on, and better than node-only.
    for (const Estimates* e : {&em, &ed}) {
        if (e->timed == 0) continue;
        CHECK(e->worst <= 0.3);
        CHECK(e->mare < e->mare_nodes);
    }
}

TEST_CASE("intelligent: probe the allocation", "[.adapt_probe5b]") {
    const std::vector<std::vector<double>> diode_points{{-1.0, 0.0}, {-0.5, 0.0}, {0.0, 0.0}, {0.3, 0.0}, {0.6, 0.0}};
    const solve::DeviceBuilder diode = [](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        const std::vector<double> x = r.axes[0].empty() ? uniform(2e-4, 21) : r.axes[0];
        return solve::AdaptableMesh{.device = legacy_diode(x), .kind = solve::AdaptMeshKind::tensor};
    };
    const std::vector<std::vector<double>> mosfet_points{{0.0, 0.05, 0.0, 0.0}, {0.0, 0.05, 0.0, 0.5}, {0.0, 0.05, 0.0, 1.0}};
    solve::BiasOptions mb;
    mb.step_control.max_cuts = 20;
    struct Case {
        const char* name;
        solve::DeviceBuilder build;
        std::vector<std::vector<double>> pts;
        solve::BiasOptions bias;
        solve::AdaptOptions o;
    };
    std::vector<Case> cases;
    {
        solve::AdaptOptions o;
        o.goals = {{0}};
        o.tolerance = 1e-3;
        o.at = {0, 4};
        o.max_passes = 14;
        cases.push_back({"diode", diode, diode_points, legacy_options(), o});
        o.weighting = solve::Weighting::adjoint;
        cases.push_back({"diode-adjoint", diode, diode_points, legacy_options(), o});
    }
    {
        solve::AdaptOptions o;
        o.goals = {{1}};
        o.tolerance = 1e-2;
        o.at = {1, 2};
        o.max_passes = 10;
        o.weighting = solve::Weighting::adjoint;
        cases.push_back({"mosfet-adjoint", tensor_mosfet_mesh, mosfet_points, mb, o});
    }
    for (Case& c : cases) {
        for (const solve::Policy policy : {solve::Policy::basic, solve::Policy::intelligent}) {
            c.o.policy = policy;
            const auto t = std::chrono::steady_clock::now();
            auto r = solve::adapt(c.build, c.pts, c.bias, c.o);
            if (!r) FAIL(r.error().message);
            for (std::size_t k = 0; k < r->passes.size(); ++k) {
                const auto& p = r->passes[k];
                std::string orders;
                for (const auto& f : p.orders) {
                    orders += std::format(" [{} r{} p{:.3f} n{}]", static_cast<int>(f.kind), f.region, f.order, f.samples);
                }
                std::string cands;
                for (std::size_t q = 0; q < std::min<std::size_t>(p.candidates.size(), 4); ++q) {
                    const auto& x = p.candidates[q];
                    cands += std::format(" (ax{} i{} b{:.3g} re{:.3g} w{:.3g} c{:.3g})", x.axis, x.index, x.benefit, x.realized,
                                         x.work, x.cost);
                }
                const double measured = p.seconds_mesh + p.seconds_transfer + p.seconds_solve;
                WARN(std::format("{} {} pass {}: {} nodes {} dofs {} nnz {} fill {} newton {}, rel {:.4g} pred_rel {:.4g} "
                                 "pred_dofs {}, cand {} taken {}, t {:.4g} pred {:.4g} pred_nodes {:.4g}{}{}",
                                 c.name, solve::to_string(policy), k, solve::to_string(p.source), p.nodes, p.dofs,
                                 p.matrix_nonzeros, p.factor_fill, p.newton_iterations, p.relative_estimate,
                                 p.predicted_relative, p.predicted_dofs, p.candidate_count, p.candidates_taken, measured,
                                 p.predicted_seconds, p.predicted_seconds_nodes, orders, cands));
            }
            WARN(std::format("{} {}: stop {}, {} passes, cert {}, misses {}, fallback {} at {}, {:.3g} s", c.name,
                             solve::to_string(policy), solve::to_string(r->stop), r->passes.size(),
                             solve::to_string(r->certificates[0].status), r->misses, solve::to_string(r->fallback),
                             r->fallback_pass, seconds_since(t)));
        }
    }
}
