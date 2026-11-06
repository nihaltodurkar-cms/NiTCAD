// TVB-06, the dielectric stack (DECISIONS Unit 55 §6, 16.4; criteria 16.8, frozen).
//
// Code under test: Poisson's equation through the small-signal path on an insulator-only device (16.22, DT-Q19 A):
// SiO2 (eps_r 3.9) and a high-k layer (eps_r 22, the spec's value) in series between two metal electrodes,
// scaled by the device's reference_material. The potential is linear in each layer and D is continuous, so the
// box method with the interface on a dual face (the midpoint of the straddling edge) reproduces, to rounding,
//     C = eps0 / (t1 / eps1 + t2 / eps2)     and     E1 / E2 = eps2 / eps1,
// t1 the distance from the first electrode to that midpoint.
// Levels: fast, 1D N = 32 and 128; full, + 2D 32 x 32 (per unit depth, width 10 nm).
// Criterion (D, frozen 16.8): C and the field ratio, each relative <= 64 eps N^2, eps = 2^-52, N the nodes across
// the stack. C is -Im Y_ab / omega at 1 MHz (the off-diagonal admittance: Y_ab = -j omega C for an ideal
// capacitor; the diagonal and the conductance are checked as well, against the same bound).
// Reference: verify/reference/capacitor (PD with a second method at L0); the inputs (thicknesses, eps_r) are the
// spec's, the same numbers NiTCAD is given (§2.2).
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/capacitor.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double eps = 0x1p-52;
constexpr double stack_cm = 15.5e-7, width_cm = 10e-7, eps_ox = 3.9, eps_hk = 22.0, frequency_Hz = 1e6;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-06";
    r.cls = ReferenceClass::D;
    r.independence = Independence::I2;
    r.provenance = Provenance::PD;
    r.citation = "Gauss's law in linear dielectrics: D uniform across a planar stack";
    r.version = "DECISIONS Unit 55 16.3, 16.4, 16.8";
    r.equations = {1};
    r.form = "C = eps0 / sum t_i / eps_i; E1 / E2 = eps2 / eps1";
    r.assumptions = {"planar layers", "electrodes equipotential", "interface on a dual face"};
    r.derivation = "div D = 0 with D normal: D uniform, E_i = D / eps_i, V = sum E_i t_i";
    r.verification = "series form against a Gauss-Kronrod integral of dx / eps(x) (capacitor_test)";
    r.second_method = "Gauss-Kronrod quadrature of the elastance density dx / eps(x) in double-double";
    return r;
}

ComparisonRecord record(const std::string& quantity, const std::string& unit, int nodes, double sim, double ref,
                        double criterion, bool extra = true) {
    const double e = std::abs(sim - ref) / std::abs(ref);
    ComparisonRecord c{quantity, unit, static_cast<double>(nodes), sim, ref, 0.0, e,
                       "relative <= 64 eps N^2 = " + shortest(criterion), {}};
    c.comparison = judge(info(), {}, e <= criterion && extra, c.criterion);
    return c;
}

}  // namespace

BenchmarkResult run_tvb06(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-06";
    r.title = "Dielectric stack between electrodes";
    r.info = info();
    physics::InsulatorParameters hk_p = physics::silicon_dioxide_parameters;
    hk_p.eps_r = eps_hk;
    const auto hk = physics::Insulator::create(hk_p);
    if (!hk) {
        r.run_error = "high-k insulator: " + hk.error().message;
        return r;
    }
    struct Level {
        int dimension, nodes;
    };
    std::vector<Level> levels{{1, 32}, {1, 128}};
    if (tier != Tier::fast) levels.push_back({2, 32});
    for (const Level lv : levels) {
        const int N = lv.nodes;
        const double h = stack_cm / (N - 1);
        const int k = static_cast<int>(std::lround(0.4 * (N - 1) - 0.5));  // last SiO2 node
        std::vector<double> x(static_cast<std::size_t>(N));
        for (int i = 0; i < N; ++i) x[static_cast<std::size_t>(i)] = h * i;
        const double t1 = 0.5 * (x[static_cast<std::size_t>(k)] + x[static_cast<std::size_t>(k + 1)]);
        const double t2 = stack_cm - t1;
        std::vector<double> y;
        for (int j = 0; j < N; ++j) y.push_back(width_cm * j / (N - 1));
        mesh::Mesh m = lv.dimension == 1 ? *mesh::make_tensor_grid(x) : *mesh::make_tensor_grid(x, y);
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n, 0);
        for (std::size_t i = 0; i < n; ++i) region[i] = m.points()[i][0] > t1 ? 1 : 0;
        device::Contact a{"a", device::ContactKind::electrode, m.find_boundary("x_min")->nodes};
        device::Contact b{"b", device::ContactKind::electrode, m.find_boundary("x_max")->nodes};
        a.electrode = b.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.5};
        device::DeviceDescription dd{.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"oxide", physics::silicon_dioxide()}, {"highk", *hk}},
                                     .node_region = std::move(region),
                                     .donors = std::vector<double>(n, 0.0),
                                     .acceptors = std::vector<double>(n, 0.0),
                                     .contacts = {std::move(a), std::move(b)},
                                     .interfaces = {{"oxide", "highk"}}};
        dd.reference_material = physics::silicon();
        const std::string label = std::to_string(lv.dimension) + "D N = " + std::to_string(N);
        const auto created = device::Device::create(std::move(dd));
        if (!created) {
            r.run_error = label + ": " + created.error().message;
            return r;
        }
        const std::vector<std::vector<double>> points{{0.0, 1.0}};
        const auto s = solve::solve_small_signal(*created, points, {.frequencies_Hz = {frequency_Hz}});
        if (!s || s->stopped) {
            r.run_error = label + ": " + (s ? s->stopped->message : s.error().message);
            return r;
        }
        const double criterion = 64.0 * eps * N * N;
        const double area = lv.dimension == 1 ? 1.0 : width_cm;
        const double C_ref =
            reference::series_capacitance(std::vector<double>{t1, t2}, std::vector<double>{eps_ox, eps_hk}).hi * area;
        const double w = 2.0 * std::numbers::pi * frequency_Hz;
        const auto y_ab = s->admittance(0, 0, 0, 1), y_aa = s->admittance(0, 0, 0, 0);
        const std::string unit = lv.dimension == 1 ? "F/cm^2" : "F/cm";
        r.comparisons.push_back(record("C (-Im Y_ab / omega), " + label, unit, N, -y_ab.imag() / w, C_ref, criterion,
                                       std::abs(y_ab.real()) <= criterion * w * C_ref));
        r.comparisons.push_back(record("C (Im Y_aa / omega), " + label, unit, N, y_aa.imag() / w, C_ref, criterion,
                                       std::abs(y_aa.real()) <= criterion * w * C_ref));
        // The field ratio from the operating point's potential: the mean field of each layer, the straddling edge
        // excluded (it holds the interface).
        const auto& psi = s->points[0].dc.fields.potential_V;
        const auto& pts = created->mesh().points();
        const auto at = [&](int xi) {  // the node at x[xi] on the y = 0 row, by its coordinates
            for (std::size_t i = 0; i < pts.size(); ++i) {
                if (pts[i][0] == x[static_cast<std::size_t>(xi)] && pts[i][1] == 0.0) return i;
            }
            return pts.size();
        };
        const std::size_t first = at(0), last_ox = at(k), first_hk = at(k + 1), last = at(N - 1);
        if (last == pts.size() || last_ox == pts.size() || first_hk == pts.size() || first == pts.size()) {
            r.run_error = label + ": a row node was not found";
            return r;
        }
        const double E1 = (psi[last_ox] - psi[first]) / (pts[last_ox][0] - pts[first][0]);
        const double E2 = (psi[last] - psi[first_hk]) / (pts[last][0] - pts[first_hk][0]);
        r.comparisons.push_back(record("field ratio E_ox / E_hk, " + label, "1", N, E1 / E2, eps_hk / eps_ox, criterion));
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
