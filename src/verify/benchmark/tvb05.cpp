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
// The (h, tol) factorial (16.2 (p), DT-Q23 A) at M = 256: tol = 1e-12, 1e-14 and M = 256, 512; eps_int enters U, and
// above 10 % of U it blocks both E passes (DECISIONS 4.1 (c)). Nightly tier only (M = 512 is about 4 minutes a solve);
// in full it is not evaluated, so the finest pass is held at not_evaluated.
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
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double a_cm = 1e-4, b_cm = 2e-4, eps_ox = 3.9, frequency_Hz = 1e6;
constexpr double tol_newton = 1e-12;  // the operating point's Newton tolerance (the factorial's tol)

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

// The annulus as M sectors and its 2M boundary patches; the electrodes' nodes from those patches. One level per Newton
// tolerance in tols, on one mesh.
std::expected<std::vector<Level>, std::string> run_levels(int M, double C_ref, const std::vector<double>& tols) {
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
    std::vector<Level> out;
    for (const double tol : tols) {
        solve::SmallSignalOptions o{.frequencies_Hz = {frequency_Hz}};
        o.steady.newton.tol_update = tol;
        const auto s = solve::solve_small_signal(*d, points, o);
        if (!s || s->stopped) {
            return std::unexpected("M = " + std::to_string(M) + ": " + (s ? s->stopped->message : s.error().message));
        }
        Level lv;
        lv.h = spacing;
        lv.C = -s->admittance(0, 0, 0, 1).imag() / (2.0 * std::numbers::pi * frequency_Hz);
        lv.error = std::abs(lv.C - C_ref) / C_ref;
        lv.u_iter = s->points[0].backward_error[0];
        out.push_back(lv);
    }
    return out;
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
        const auto lvs = run_levels(M, C_ref, {tol_newton});
        if (!lvs) {
            r.run_error = lvs.error();
            return r;
        }
        const Level* lv = &lvs->front();
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
    // The (h, tol) factorial at M = 256 (16.2 (p), DT-Q23 A): the operating point's Newton tol and tol / 100, M = 256
    // and 512. U at (h, tol): U_disc (GCI of the finest pair with p-hat, F_s = 1.25) and U_iter (U_ref is rounding).
    // M = 512 takes about 4 minutes per solve (131k nodes; 16.30), beyond 16.7's 120 s for tvr4_full on CI. So it runs
    // in the nightly tier (the local verify_full.tvb05); in full the factorial is not evaluated and the finest pass is
    // held at not_evaluated (4.3: without its ladder U_num is not estimated).
    const double u_disc = reference::gci_safety_factor_three_grids * std::abs(levels[1].error - levels[0].error) /
                          (std::pow(levels[1].h / levels[0].h, st.observed) - 1.0);
    Factorial fac;
    if (tier < Tier::nightly) {
        fac.note = "factorial (h, tol): not evaluated in this tier (M = 512 runs in nightly; 16.7, 16.30)";
    } else {
        const auto at_h = run_levels(256, C_ref, {0.01 * tol_newton});
        const auto at_h2 = run_levels(512, C_ref, {tol_newton, 0.01 * tol_newton});
        if (!at_h || !at_h2) {
            r.run_error = (at_h ? at_h2.error() : at_h.error()) + " (factorial)";
            return r;
        }
        const std::vector<double> parts{u_disc, levels.front().u_iter};
        fac = factorial(levels.front().error, at_h->front().error, (*at_h2)[0].error, (*at_h2)[1].error, parts);
        fac.note += " (U_disc " + shortest(u_disc) + ")";
    }
    r.notes.push_back(fac.note);
    for (ComparisonRecord& c : r.comparisons) block_if_interacting(c, fac);
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
