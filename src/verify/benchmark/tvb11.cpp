// TVB-11, the uniform resistor (DECISIONS Unit 55 §6, 14.4; criteria 14.7).
//
// Code under test: the drift-diffusion system and solve_bias on a uniformly doped bar with ohmic contacts at both
// ends: psi is linear, n = n0 and p = p0 everywhere, and the current is pure drift,
//     J = q (mu_n n0 + mu_p p0) V / L,
// which Scharfetter-Gummel reproduces exactly at constant field and density. Silicon, N_D = 1e16 cm^-3, L = 1 um,
// V = 0.1 V; constant mobilities (doping_mobility off, Caughey-Thomas holes: the material's mu_max at 300 K).
// Levels: 1D with N = 32 and 128 nodes, 2D 32 x 32 (the current per unit width; width 1 um).
// Criterion (D, frozen): |I - I_ref| / |I_ref| and max_i |psi_i - psi_lin,i| / max|psi| each <= 64 eps N^2,
// eps = 2^-52, N the nodes along the current path.
// Reference independence I1: n0, p0 and the contact potentials are NiTCAD's Dirichlet values (input data), the
// mobilities the material's parameters; the rest is exact arithmetic.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double length_cm = 1e-4, width_cm = 1e-4, doping = 1e16, bias_V = 0.1;
constexpr double eps = 0x1p-52;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-11";
    r.cls = ReferenceClass::D;
    r.independence = Independence::I1;
    r.provenance = Provenance::PA;
    r.citation = "Ohm's law for a uniform semiconductor bar: J = q (mu_n n + mu_p p) E";
    r.version = "DECISIONS Unit 55 14.4, 14.7";
    r.equations = {1, 2, 3};
    r.form = "I = q (mu_n n0 + mu_p p0) V A / L; psi linear";
    r.assumptions = {"uniform doping", "constant mobilities", "ohmic contacts at both ends"};
    r.derivation = "constant field and density: the diffusion current vanishes";
    r.verification = "D: exact to rounding";
    return r;
}

std::vector<double> uniform(double length, int nodes) {
    std::vector<double> v;
    for (int i = 0; i < nodes; ++i) v.push_back(length * i / (nodes - 1));
    return v;
}

}  // namespace

BenchmarkResult run_tvb11() {
    BenchmarkResult r;
    r.id = "TVB-11";
    r.title = "Uniform resistor with ohmic contacts";
    r.info = info();
    const physics::Semiconductor si = physics::silicon();
    const double mu_n = si.parameters().electron_mobility.mu_max, mu_p = si.parameters().hole_mobility.mu_max;
    struct Level {
        int dimension, nodes;
    };
    for (const Level lv : {Level{1, 32}, Level{1, 128}, Level{2, 32}}) {
        const auto x = uniform(length_cm, lv.nodes);
        mesh::Mesh m = lv.dimension == 1 ? *mesh::make_tensor_grid(x)
                                         : *mesh::make_tensor_grid(x, uniform(width_cm, lv.nodes));
        const std::size_t n = m.node_count();
        auto left = m.find_boundary("x_min")->nodes;
        auto right = m.find_boundary("x_max")->nodes;
        const auto created = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", si}},
             .node_region = std::vector<device::RegionId>(n, 0),
             .donors = std::vector<double>(n, doping),
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                          {"right", device::ContactKind::ohmic, std::move(right)}}});
        const std::string label = std::to_string(lv.dimension) + "D N = " + std::to_string(lv.nodes);
        if (!created) {
            r.run_error = label + ": " + created.error().message;
            return r;
        }
        const device::Device& d = *created;
        solve::BiasOptions o;
        o.models.doping_mobility = false;
        o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
        o.newton.tol_update = 1e-12;
        const auto s = solve::solve_bias(d, std::vector<double>{0.0, bias_V}, o);
        const double criterion = 64.0 * eps * lv.nodes * lv.nodes;
        const std::string crit = "relative <= 64 eps N^2 = " + shortest(criterion);
        if (!s) {
            r.notes.push_back(label + ": the solve failed: " + s.error().message);
            for (const char* q : {"current", "psi"}) {
                ComparisonRecord c{std::string(q) + ", " + label, "1", static_cast<double>(lv.nodes), 0.0, 0.0,
                                   0.0, 0.0, crit, {}};
                c.comparison = judge(r.info, {}, false, c.criterion);
                r.comparisons.push_back(std::move(c));
            }
            continue;
        }
        const auto& f = s->fields;
        const std::size_t a = static_cast<std::size_t>(d.contacts()[0].nodes.front());
        const std::size_t b = static_cast<std::size_t>(d.contacts()[1].nodes.front());
        const double n0 = f.n_cm3[a], p0 = f.p_cm3[a];
        const double J = base::q_C * (mu_n * n0 + mu_p * p0) * (f.potential_V[b] - f.potential_V[a]) / length_cm;
        const double I_ref = lv.dimension == 1 ? J : J * width_cm;
        const double I = s->terminal_current[1];
        const double e_I = std::abs(std::abs(I) - std::abs(I_ref)) / std::abs(I_ref);
        ComparisonRecord ci{"current, " + label, lv.dimension == 1 ? "A/cm^2" : "A/cm", static_cast<double>(lv.nodes),
                            std::abs(I), std::abs(I_ref), 0.0, e_I, crit, {}};
        ci.comparison = judge(r.info, {}, e_I <= criterion && I * s->terminal_current[0] <= 0.0, ci.criterion);
        r.comparisons.push_back(std::move(ci));
        double worst = 0.0, scale = 0.0;
        const auto& pts = d.mesh().points();
        for (std::size_t i = 0; i < n; ++i) {
            const double lin = f.potential_V[a] + (f.potential_V[b] - f.potential_V[a]) * (pts[i][0] / length_cm);
            worst = std::max(worst, std::abs(f.potential_V[i] - lin));
            scale = std::max(scale, std::abs(f.potential_V[i]));
        }
        ComparisonRecord cp{"psi linear, " + label, "1", static_cast<double>(lv.nodes), worst, 0.0, 0.0,
                            worst / scale, crit, {}};
        cp.comparison = judge(r.info, {}, worst / scale <= criterion, cp.criterion);
        r.comparisons.push_back(std::move(cp));
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
