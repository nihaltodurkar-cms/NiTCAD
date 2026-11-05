// TVB-12, Scharfetter-Gummel constant-field exactness (DECISIONS Unit 55 §6, 14.4; criteria 14.7).
//
// Code under test: the drift-diffusion system's SG carrier fluxes and box-method Poisson rows (assemble), solved by
// solve_bias. Construction (14.4): 1D, N = 64 uniform nodes on [0, 1 um]; constant mobilities (doping_mobility
// off, Caughey-Thomas holes: mu_max), no generation or recombination, no band-gap narrowing, Boltzmann. The ohmic
// contacts fix psi, n, p at both ends (left N_D = 1e18, right N_D = 1e15 cm^-3; interior 1e16). The exact solution
// with constant field and constant currents is
//     psi linear;  n = n_L + (n_R - n_L) (e^{(psi - psi_L)/V_T} - 1) / (e^{(psi_R - psi_L)/V_T} - 1),
//     and p likewise with -psi,
// the two homogeneous solutions (e^{+psi/V_T} for electrons, e^{-psi/V_T} for holes) that SG reproduces exactly at
// constant field, written without overflow. The manufactured charge n - p - (N_D - N_A) at each node (the TVQ3 hook)
// makes psi linear the discrete Poisson solution. The bias sets the Peclet number per edge, |d psi / V_T| in
// {0.01, 1, 10, 40}. The solve starts from the exact state perturbed by 1e-3 and must return to it.
// Criterion (D, frozen): max over nodes of |n - n_exact| / n_exact and the same for p <= 64 eps N^2, eps = 2^-52.
// Reference independence I1: the boundary values come from NiTCAD's own equilibrium contact values (input data);
// the interior is exact arithmetic on them.
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
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

constexpr int nodes = 64;
constexpr double length_cm = 1e-4;
constexpr double eps = 0x1p-52;
constexpr double criterion = 64.0 * eps * nodes * nodes;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-12";
    r.cls = ReferenceClass::D;
    r.independence = Independence::I1;
    r.provenance = Provenance::PA;
    r.citation = "Scharfetter and Gummel, IEEE Trans. Electron Devices 16 (1969): exact for constant field and current";
    r.version = "DECISIONS Unit 55 14.4, 14.7";
    r.equations = {1, 2, 3};
    r.form = "psi linear; n = A e^{psi/V_T} + K, p = A' e^{-psi/V_T} + K' through the contact values";
    r.assumptions = {"1D uniform mesh", "constant mobilities", "no generation or recombination", "Boltzmann",
                     "manufactured charge cancels the space charge (TVQ3)"};
    r.derivation = "the homogeneous solutions of the constant-field continuity equations";
    r.verification = "D: exact to rounding";
    return r;
}

// (e^{u} - 1) / (e^{w} - 1) for u between 0 and w, without overflow.
double fraction(double u, double w) {
    if (w == 0.0) return 0.0;
    if (w > 0.0) return std::exp(u - w) * std::expm1(-u) / std::expm1(-w);
    return std::expm1(u) / std::expm1(w);
}

device::Device make_device(std::vector<double>& x) {
    x.clear();
    for (int i = 0; i < nodes; ++i) x.push_back(length_cm * i / (nodes - 1));
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    std::vector<double> donors(nodes, 1e16), acceptors(nodes, 0.0);
    donors.front() = 1e18;
    donors.back() = 1e15;
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(nodes, 0),
                                    .donors = std::move(donors),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                 {"right", device::ContactKind::ohmic, std::move(right)}}});
}

solve::BiasOptions options() {
    solve::BiasOptions o;
    o.models.doping_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.newton.tol_update = 1e-12;
    return o;
}

}  // namespace

BenchmarkResult run_tvb12() {
    BenchmarkResult r;
    r.id = "TVB-12";
    r.title = "Scharfetter-Gummel constant-field exactness";
    r.info = info();
    std::vector<double> x;
    const device::Device d = make_device(x);
    const double VT = base::thermal_voltage(300.0);
    const auto eq = solve::solve_bias(d, std::vector<double>{0.0, 0.0}, options());
    if (!eq) {
        r.run_error = "equilibrium: " + eq.error().message;
        return r;
    }
    const std::size_t L = 0, R = nodes - 1;
    const double psiL = eq->fields.potential_V[L], psiR0 = eq->fields.potential_V[R];
    const double nL = eq->fields.n_cm3[L], nR = eq->fields.n_cm3[R];
    const double pL = eq->fields.p_cm3[L], pR = eq->fields.p_cm3[R];
    for (const double peclet : {0.01, 1.0, 10.0, 40.0}) {
        const double V = peclet * (nodes - 1) * VT - (psiR0 - psiL);
        const double psiR = psiR0 + V;
        const double w = (psiR - psiL) / VT;
        results::NodeFields exact;
        auto sources = std::make_shared<assemble::ManufacturedSources>();
        sources->charge_cm3.assign(nodes, 0.0);
        for (std::size_t i = 0; i < nodes; ++i) {
            const double psi = psiL + (psiR - psiL) * (x[i] / length_cm);
            const double u = (psi - psiL) / VT;
            const double n = nL + (nR - nL) * fraction(u, w);
            const double p = pL + (pR - pL) * fraction(-u, -w);
            exact.potential_V.push_back(psi);
            exact.n_cm3.push_back(n);
            exact.p_cm3.push_back(p);
            sources->charge_cm3[i] = n - p - (d.donors()[i] - d.acceptors()[i]);
        }
        solve::BiasOptions o = options();
        o.models.manufactured = sources.get();  // kept alive by sources through the solve
        results::NodeFields start = exact;
        for (std::size_t i = 0; i < nodes; ++i) {
            start.potential_V[i] += 1e-3 * VT * ((i % 2 == 0) ? 1.0 : -1.0);
            start.n_cm3[i] *= 1.0 + 1e-3;
            start.p_cm3[i] *= 1.0 - 1e-3;
        }
        const std::string label = "Peclet " + shortest(peclet);
        const auto s = solve::solve_bias(d, std::vector<double>{0.0, V}, o, &start);
        if (!s) {
            r.notes.push_back(label + ": the solve failed: " + s.error().message);
            for (const char* q : {"n", "p"}) {
                ComparisonRecord c{std::string("nodal ") + q + ", " + label, "1", peclet, 0.0, 0.0, 0.0, 0.0,
                                   "relative <= 64 eps N^2", {}};
                c.comparison = judge(r.info, {}, false, c.criterion);
                r.comparisons.push_back(std::move(c));
            }
            continue;
        }
        for (const bool electrons : {true, false}) {
            const auto& sim = electrons ? s->fields.n_cm3 : s->fields.p_cm3;
            const auto& ref = electrons ? exact.n_cm3 : exact.p_cm3;
            double worst = 0.0, at_sim = 0.0, at_ref = 0.0;
            for (std::size_t i = 0; i < nodes; ++i) {
                const double e = std::abs(sim[i] - ref[i]) / ref[i];
                if (e >= worst) {
                    worst = e;
                    at_sim = sim[i];
                    at_ref = ref[i];
                }
            }
            ComparisonRecord c{std::string("nodal ") + (electrons ? "n" : "p") + ", " + label, "cm^-3", peclet,
                               at_sim, at_ref, 0.0, worst, "relative <= 64 eps N^2 = " + shortest(criterion), {}};
            c.comparison = judge(r.info, {}, worst <= criterion, c.criterion);
            r.comparisons.push_back(std::move(c));
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
