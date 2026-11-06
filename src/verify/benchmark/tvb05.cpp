// TVB-05, the coaxial capacitor on 2D generated meshes (DECISIONS Unit 55 §6, 16.4; criteria 16.8, frozen).
//
// Code under test: Poisson's equation on unstructured (generated) 2D meshes, through the small-signal path, on an
// insulator-only device (16.22, DT-Q19 A): SiO2 between two coaxial metal electrodes, radii a = 1 um and b = 2 um,
// per unit length. The annulus is tiled by M quadrilateral sectors of one region (the 2D generator takes polygons
// that meet along whole edges; their union has the hole), so both circles are M-gons, and M doubles with each
// level: the faceting error and the discretization error are both O(h^2), which is the frozen p_th = 2 (§6: the
// faceting is included in p_th).
// Reference: C = 2 pi eps0 eps_r / ln(b / a) of the circles (verify/reference/capacitor; PD, second method at L0).
// Levels: M = 32, 64, 128, 256 with h = 2 pi a / M (full tier only).
// Criteria (E, frozen 16.8): p-hat of the relative error of C in [1.8, 2.3] with the asymptotic check (4 levels);
// on the finest level, relative error <= 1e-3. C = -Im Y_ab / omega at 1 MHz. The AC solve is linear and direct,
// so the iterative error is the linear solve's backward error (relative).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <numbers>
#include <string>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/capacitor.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double a_cm = 1e-4, b_cm = 2e-4, eps_ox = 3.9, frequency_Hz = 1e6;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-05";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I2;
    r.provenance = Provenance::PD;
    r.citation = "Gauss's law for a coaxial capacitor: D = Q / (2 pi r)";
    r.version = "DECISIONS Unit 55 16.3, 16.4, 16.8";
    r.equations = {1};
    r.form = "C = 2 pi eps0 eps_r / ln(b / a) per unit length";
    r.assumptions = {"circular electrodes (the mesh's are M-gons; faceting O(h^2), in p_th)", "per unit length"};
    r.derivation = "div D = 0 with D radial: D = Q / (2 pi r), V = Q ln(b/a) / (2 pi eps)";
    r.verification = "against a Gauss-Kronrod integral of dr / (2 pi eps r) (capacitor_test)";
    r.second_method = "Gauss-Kronrod quadrature of the elastance density dr / (2 pi eps r) in double-double";
    return r;
}

struct Level {
    double h = 0.0, error = 0.0, u_iter = 0.0, C = 0.0;
};

// The annulus as M sectors and its 2M boundary patches; the electrodes' nodes from those patches.
std::expected<Level, std::string> run_level(int M, double C_ref) {
    std::vector<mesh::Point> inner(static_cast<std::size_t>(M)), outer(static_cast<std::size_t>(M));
    for (int l = 0; l < M; ++l) {
        const double t = 2.0 * std::numbers::pi * l / M;
        inner[static_cast<std::size_t>(l)] = {a_cm * std::cos(t), a_cm * std::sin(t), 0.0};
        outer[static_cast<std::size_t>(l)] = {b_cm * std::cos(t), b_cm * std::sin(t), 0.0};
    }
    mesh::Geometry2D g;
    for (int l = 0; l < M; ++l) {
        const auto i0 = static_cast<std::size_t>(l), i1 = static_cast<std::size_t>((l + 1) % M);
        g.polygons.push_back({{inner[i0], outer[i0], outer[i1], inner[i1]}, 0});
        g.patches.push_back({"in_" + std::to_string(l), inner[i0], inner[i1]});
        g.patches.push_back({"out_" + std::to_string(l), outer[i0], outer[i1]});
    }
    const double h = 2.0 * std::numbers::pi * a_cm / M;
    auto gm = mesh::generate_mesh(g, {.max_size = h});
    if (!gm) return std::unexpected("M = " + std::to_string(M) + ": generate_mesh: " + gm.error().message);
    const double spacing = study::mean_spacing(gm->mesh);
    const auto nodes_of = [&](const char* prefix) {
        std::vector<mesh::NodeId> v;
        for (int l = 0; l < M; ++l) {
            const auto* p = gm->mesh.find_boundary(prefix + std::to_string(l));
            if (p != nullptr) v.insert(v.end(), p->nodes.begin(), p->nodes.end());
        }
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        return v;
    };
    device::Contact ca{"inner", device::ContactKind::electrode, nodes_of("in_")};
    device::Contact cb{"outer", device::ContactKind::electrode, nodes_of("out_")};
    ca.electrode = cb.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.5};
    const std::size_t n = gm->mesh.node_count();
    device::DeviceDescription dd{.mesh = std::move(gm->mesh),
                                 .temperature_K = 300.0,
                                 .regions = {{"oxide", physics::silicon_dioxide()}},
                                 .node_region = std::vector<device::RegionId>(n, 0),
                                 .donors = std::vector<double>(n, 0.0),
                                 .acceptors = std::vector<double>(n, 0.0),
                                 .contacts = {std::move(ca), std::move(cb)}};
    dd.reference_material = physics::silicon();
    const auto d = device::Device::create(std::move(dd));
    if (!d) return std::unexpected("M = " + std::to_string(M) + ": " + d.error().message);
    const std::vector<std::vector<double>> points{{0.0, 1.0}};
    const auto s = solve::solve_small_signal(*d, points, {.frequencies_Hz = {frequency_Hz}});
    if (!s || s->stopped) {
        return std::unexpected("M = " + std::to_string(M) + ": " + (s ? s->stopped->message : s.error().message));
    }
    Level lv;
    lv.h = spacing;
    lv.C = -s->admittance(0, 0, 0, 1).imag() / (2.0 * std::numbers::pi * frequency_Hz);
    lv.error = std::abs(lv.C - C_ref) / C_ref;
    lv.u_iter = s->points[0].backward_error[0];
    return lv;
}

}  // namespace

BenchmarkResult run_tvb05(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-05";
    r.title = "Coaxial capacitor on generated 2D meshes";
    r.info = info();
    const std::string band = "p-hat in [1.8, 2.3] (p_th 2, faceting included)", accuracy = "relative <= 1e-3";
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.4)");
        return r;
    }
    const double C_ref = reference::coaxial_capacitance(a_cm, b_cm, eps_ox).hi;
    std::vector<Level> levels;  // finest first
    for (const int M : {256, 128, 64, 32}) {
        const auto lv = run_level(M, C_ref);
        if (!lv) {
            r.run_error = lv.error();
            return r;
        }
        r.notes.push_back("M = " + std::to_string(M) + ": h = " + shortest(lv->h) + " cm, C = " + shortest(lv->C) +
                          " F/cm, relative error " + shortest(lv->error) + ", backward error " +
                          shortest(lv->u_iter));
        levels.push_back(*lv);
    }
    std::vector<double> hs, e, u;
    for (const Level& l : levels) {
        hs.push_back(l.h);
        e.push_back(l.error);
        u.push_back(l.u_iter);
    }
    const study::OrderStudy st = study::order_study(hs, e, u);
    const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
    std::string orders;
    for (const double p : st.orders) orders += " " + shortest(p);
    r.notes.push_back("orders (finest pair first):" + orders + "; spread " + shortest(st.order_spread));
    ComparisonRecord co{"C order", "1", 4.0, st.observed, 2.0, 0.0, st.observed, band, {}};
    co.comparison = v.judged ? judge(r.info, {}, v.pass, band) : Comparison{Outcome::not_evaluated, v.reason, {}};
    r.comparisons.push_back(std::move(co));
    ComparisonRecord ca{"C finest (M = 256)", "F/cm", 256.0, levels.front().C, C_ref, 0.0, levels.front().error,
                        accuracy, {}};
    ca.comparison = judge(r.info, {}, levels.front().error <= 1e-3, accuracy);
    r.comparisons.push_back(std::move(ca));
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
