// TVB-04, the Poisson manufactured solution with a dielectric jump (DECISIONS Unit 55 §6, 14.4, CV-Q2 A; criteria
// 14.7).
//
// Code under test: the drift-diffusion system's box-method Poisson rows (edge fluxes, the harmonic-mean coupling of
// an insulator-semiconductor edge, the lumped space charge) and solve_bias, at thermal equilibrium. Construction:
// SiO2 on x in [-L, 0] with a metal electrode at x = -L, p-silicon (N_A = 1e17 cm^-3) on [0, L] with an ohmic contact
// at x = L, L = 50 nm; in 2D and 3D the cross-section is [0, L]^(D-1) with natural (zero-flux) sides. The interface
// lies at the midpoint of the edges that join the two regions (mesh::straddle_interface, the decks' axis). The exact
// potential is
//     psi = psi_c + A g(x / L) H,   H = h(y / L) h(z / L),   h(v) = 1 + cos(pi v) / 2   (h' = 0 on the sides),
//     g(u) = (1 - u) e^{alpha u} in silicon,   g(u) = (1 + u) e^{beta u} in the oxide,
// with g = 0 at both contacts, g continuous at u = 0 and eps_ox (1 + beta) = eps_si (alpha - 1) so the normal
// displacement is continuous (no interface charge); alpha = 2, A = 0.2 V. At equilibrium with Boltzmann statistics the
// carriers are n = n_c e^{(psi - psi_c) / V_T}, p = p_c e^{-(psi - psi_c) / V_T}, which make every Scharfetter-Gummel
// flux vanish, so the continuity rows hold and the problem is Poisson-Boltzmann. The manufactured charge (the TVQ3
// hook, lumped at the nodes, CV-Q5) is
//     S = -eps Lap(psi) / q + n - p - (N_D - N_A)   in silicon,   S = -eps_ox Lap(psi) / q   in the oxide,
// so psi solves -div(eps grad psi) = q (p - n + N_D - N_A + S) exactly.
// Criteria (E, frozen 14.7): on tensor ladders the L2 and L-infinity orders p-hat in [1.8, 2.3] (p_th = 2), and the
// finest level's L2 error relative to the potential scale A + V_T at most 1e-4. Errors are against the exact psi at
// the nodes; L2 = sqrt(sum V_i e_i^2 / sum V_i). U_iter per level is the last Newton correction (scaled potential,
// times V_T) over the scale, and an order is computed only when U_iter <= e / 10 (CV-Q8).
// Tiers (CV-Q9): fast, 1D N = 32 ... 256; full, + 2D tensor 16^2 ... 128^2 and 3D tensor 8^3 ... 64^3 (the 14.7
// ladder column, owner ruling CV-Q13 in 14.9), the 3D levels solved by MKL PARDISO (ruling 2026-11-05: Eigen SparseLU
// does not fit 64^3); without the MKL runtime loaded they are not_evaluated. The finest 2D and 3D L2 misses are known
// finding (iv) (CV-Q12).
// Reference independence I1: psi_c, n_c, p_c and the electrode's potential at zero bias are NiTCAD's own equilibrium
// contact values (input data); the rest is analytic. Second method: the analytic first and second derivatives of g
// and h are checked against complex-step derivatives at run time (verify/reference/complex_step.hpp).
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <expected>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double L = 5e-6;  // cm: oxide and silicon thickness, and the cross-section's side
constexpr double acceptors_cm3 = 1e17;
constexpr double amplitude_V = 0.2;
constexpr double alpha = 2.0;
constexpr double l2_tolerance = 1e-4;
constexpr double pi = std::numbers::pi;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-04";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PA;
    r.citation = "Method of manufactured solutions (Roache, Verification and Validation in Computational Science and "
                 "Engineering, 1998; Salari and Knupp, SAND2000-1444); box-method orders per DECISIONS Unit 55 5.3 "
                 "(Eymard, Gallouet and Herbin, Finite Volume Methods, 2000 [unverified])";
    r.version = "DECISIONS Unit 55 14.4, 14.7";
    r.equations = {1};
    r.form = "psi = psi_c + A g(x/L) h(y/L) h(z/L); Poisson-Boltzmann with the manufactured charge S";
    r.assumptions = {"thermal equilibrium", "Boltzmann statistics", "planar SiO2-Si interface at an edge midpoint",
                     "manufactured charge lumped at the nodes (TVQ3, CV-Q5)"};
    r.derivation = "DECISIONS Unit 55 14.4 and this file's header: S from the exact Laplacian";
    r.verification = "E: observed order and finest-level accuracy on tensor ladders";
    r.second_method = "complex-step derivatives of the profile against the analytic ones, checked at run time";
    return r;
}

struct Materials {
    double eps_si, eps_ox, beta;
};

Materials materials() {
    const double si = physics::silicon().parameters().eps_r;
    const double ox = physics::silicon_dioxide().parameters().eps_r;
    return {si, ox, si / ox * (alpha - 1.0) - 1.0};
}

// g on the side of u (silicon for u >= 0), with its analytic derivatives.
template <class T>
T g(T u, bool silicon, double beta) {
    using std::exp;
    return silicon ? (1.0 - u) * exp(alpha * u) : (1.0 + u) * exp(beta * u);
}
template <class T>
T g1(T u, bool silicon, double beta) {
    using std::exp;
    return silicon ? exp(alpha * u) * (alpha * (1.0 - u) - 1.0) : exp(beta * u) * (beta * (1.0 + u) + 1.0);
}
double g2(double u, bool silicon, double beta) {
    return silicon ? std::exp(alpha * u) * (alpha * alpha * (1.0 - u) - 2.0 * alpha)
                   : std::exp(beta * u) * (beta * beta * (1.0 + u) + 2.0 * beta);
}
template <class T>
T h(T v) {
    using std::cos;
    return 1.0 + 0.5 * cos(pi * v);
}
template <class T>
T h1(T v) {
    using std::sin;
    return -0.5 * pi * sin(pi * v);
}
double h2(double v) { return -0.5 * pi * pi * std::cos(pi * v); }

// The analytic derivatives against complex-step ones; empty when they agree.
std::string check_derivatives(double beta) {
    using reference::complex_d;
    const auto close = [](double a, double b) { return std::abs(a - b) <= 1e-13 * std::max(1.0, std::abs(b)); };
    for (const bool silicon : {true, false}) {
        for (int k = 0; k <= 8; ++k) {
            const double u = (silicon ? 1.0 : -1.0) * k / 8.0;
            const auto d1 = reference::complex_step([&](complex_d z) { return g(z, silicon, beta); }, u);
            const auto d2 = reference::complex_step([&](complex_d z) { return g1(z, silicon, beta); }, u);
            if (!close(d1.derivative, g1(u, silicon, beta)) || !close(d2.derivative, g2(u, silicon, beta))) {
                return "g's derivatives disagree with complex step at u = " + shortest(u);
            }
        }
    }
    for (int k = 0; k <= 8; ++k) {
        const double v = k / 8.0;
        const auto d1 = reference::complex_step([](complex_d z) { return h(z); }, v);
        const auto d2 = reference::complex_step([](complex_d z) { return h1(z); }, v);
        if (!close(d1.derivative, h1(v)) || !close(d2.derivative, h2(v))) {
            return "h's derivatives disagree with complex step at v = " + shortest(v);
        }
    }
    return {};
}

// The x axis: uniform spacing L / n_side on [-L, L], then the interface straddled at 0 (2 n_side + 2 nodes).
std::vector<double> x_axis(int n_side) {
    std::vector<double> x;
    for (int k = -n_side; k <= n_side; ++k) x.push_back(L * k / n_side);
    return *mesh::straddle_interface(std::move(x), 0.0, L / n_side);
}

std::vector<double> uniform(int nodes) {
    std::vector<double> v;
    for (int k = 0; k < nodes; ++k) v.push_back(L * k / (nodes - 1));
    return v;
}

device::Device make_device(mesh::Mesh m) {
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const bool oxide = m.points()[i][0] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (!oxide) acceptors[i] = acceptors_cm3;
    }
    device::Contact gate{"gate", device::ContactKind::electrode, m.find_boundary("x_min")->nodes};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate), {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
}

solve::BiasOptions options() {
    solve::BiasOptions o;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.newton.tol_update = 1e-12;
    return o;
}

struct Contacts {
    double psi_c, n_c, p_c, gate_bias_V;
};

struct LevelError {
    double h, l2, linf, u_iter;
    std::size_t nodes;
};

std::expected<LevelError, std::string> run_level(const device::Device& d, double h_nominal, const Contacts& c,
                                                 const Materials& mat, bool pardiso) {
    const std::size_t n = d.mesh().node_count();
    const int D = d.mesh().dimension();
    const double VT = base::thermal_voltage(300.0);
    results::NodeFields exact;
    auto sources = std::make_shared<assemble::ManufacturedSources>();
    sources->charge_cm3.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& pt = d.mesh().points()[i];
        const double u = pt[0] / L;
        const bool silicon = !d.is_insulator(static_cast<mesh::NodeId>(i));
        double H = 1.0, lapH = 0.0;  // H and its Laplacian over the cross-section, in 1/L^2 units
        if (D >= 2) {
            const double vy = pt[1] / L;
            if (D == 2) {
                H = h(vy);
                lapH = h2(vy);
            } else {
                const double vz = pt[2] / L;
                H = h(vy) * h(vz);
                lapH = h2(vy) * h(vz) + h(vy) * h2(vz);
            }
        }
        const double gu = g(u, silicon, mat.beta);
        const double psi = c.psi_c + amplitude_V * gu * H;
        const double lap = amplitude_V * (g2(u, silicon, mat.beta) * H + gu * lapH) / (L * L);  // V/cm^2
        const double eps = (silicon ? mat.eps_si : mat.eps_ox) * base::eps0_F_per_cm;
        double S = -eps * lap / base::q_C;
        const double nn = c.n_c * std::exp((psi - c.psi_c) / VT);
        const double pp = c.p_c * std::exp(-(psi - c.psi_c) / VT);
        if (silicon) S += nn - pp - (d.donors()[i] - d.acceptors()[i]);
        exact.potential_V.push_back(psi);
        exact.n_cm3.push_back(nn);
        exact.p_cm3.push_back(pp);
        sources->charge_cm3[i] = S;
    }
    solve::BiasOptions o = options();
    o.models.manufactured = sources.get();  // kept alive by sources through the solve
    if (pardiso) o.linear = {.backend = linalg::SolverBackend::mkl_pardiso, .threads = pardiso_threads()};
    const auto s = solve::solve_bias(d, std::vector<double>{c.gate_bias_V, 0.0}, o, &exact);
    if (!s) return std::unexpected(s.error().message);
    const double scale = amplitude_V + VT;
    double sum = 0.0, volume = 0.0, worst = 0.0;
    const auto V = d.mesh().volumes();
    for (std::size_t i = 0; i < n; ++i) {
        const double e = s->fields.potential_V[i] - exact.potential_V[i];
        sum += V[i] * e * e;
        volume += V[i];
        worst = std::max(worst, std::abs(e));
    }
    const auto& it = s->convergence.iterations;
    const double u_iter = it.empty() ? 0.0 : it.back().update * VT / scale;
    return LevelError{h_nominal, std::sqrt(sum / volume) / scale, worst / scale, u_iter, n};
}

struct Family {
    std::string name;
    int dimension;
    std::vector<int> n_side;  // finest first: n_side per x half; the cross-section has 2 n_side + 2 nodes per axis
    bool pardiso = false;     // the 3D levels: MKL PARDISO (CV-Q13 ruling; Eigen SparseLU does not fit 64^3)
};

void judge_family(BenchmarkResult& r, const Family& f, const Contacts& c, const Materials& mat) {
    if (f.pardiso && !linalg::loaded_mkl_runtime()) {  // not run: its input, the MKL runtime, is missing
        const std::string why = "the MKL runtime for the " + f.name + " levels' PARDISO solves (nitcad verify --mkl)";
        for (const char* q : {"psi L2 order, ", "psi Linf order, ", "psi L2 finest, "}) {
            ComparisonRecord rec{q + f.name, "1", 0.0, 0.0, 0.0, 0.0, 0.0, "needs PARDISO (CV-Q13)", {}};
            rec.comparison = missing_input(why);
            r.comparisons.push_back(std::move(rec));
        }
        r.notes.push_back(f.name + ": not evaluated without " + why);
        return;
    }
    std::vector<LevelError> levels;
    for (const int k : f.n_side) {
        const std::vector<double> x = x_axis(k);
        const int nodes = static_cast<int>(x.size());
        std::expected<mesh::Mesh, base::Error> m = f.dimension == 1   ? mesh::make_tensor_grid(x)
                                                   : f.dimension == 2 ? mesh::make_tensor_grid(x, uniform(nodes))
                                                                      : mesh::make_tensor_grid(x, uniform(nodes),
                                                                                               uniform(nodes));
        if (!m) {
            r.run_error = f.name + ": " + m.error().message;
            return;
        }
        const auto e = run_level(make_device(std::move(*m)), L / k, c, mat, f.pardiso);
        if (!e) {
            r.notes.push_back(f.name + " N = " + std::to_string(nodes) + ": the solve failed: " + e.error());
            levels.clear();
            break;
        }
        levels.push_back(*e);
        r.notes.push_back(f.name + " N = " + std::to_string(nodes) + " (" + std::to_string(e->nodes) +
                          " nodes): L2 " + shortest(e->l2) + ", Linf " + shortest(e->linf) + ", U_iter " +
                          shortest(e->u_iter));
    }
    const auto record = [&](const std::string& quantity, double sim, double ref, double metric, const std::string& crit,
                            Comparison cmp) {
        ComparisonRecord rec{quantity, "1", levels.empty() ? 0.0 : levels.front().h, sim, ref, 0.0, metric, crit, {}};
        rec.comparison = std::move(cmp);
        r.comparisons.push_back(std::move(rec));
    };
    const std::string band = "p-hat in [1.8, 2.3]";
    if (levels.size() != f.n_side.size()) {
        for (const char* norm : {"L2", "Linf"}) {
            record(std::string("psi ") + norm + " order, " + f.name, 0.0, 2.0, 0.0, band,
                   judge(r.info, {}, false, band));
        }
        record("psi L2 finest, " + f.name, 0.0, 0.0, 0.0, "relative <= 1e-4", judge(r.info, {}, false, "1e-4"));
        return;
    }
    std::vector<double> hs, u;
    for (const auto& l : levels) {
        hs.push_back(l.h);
        u.push_back(l.u_iter);
    }
    for (const bool l2 : {true, false}) {
        std::vector<double> e;
        for (const auto& l : levels) e.push_back(l2 ? l.l2 : l.linf);
        const study::OrderStudy st = study::order_study(hs, e, u);
        const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
        std::string spread = st.asymptotic_evaluated ? ", order spread " + shortest(st.order_spread) : "";
        r.notes.push_back(f.name + (l2 ? " L2" : " Linf") + " orders:" + [&] {
            std::string s;
            for (const double p : st.orders) s += " " + shortest(p);
            return s;
        }() + spread);
        record(std::string("psi ") + (l2 ? "L2" : "Linf") + " order, " + f.name, st.observed, 2.0, st.observed, band,
               v.judged ? judge(r.info, {}, v.pass, band) : Comparison{Outcome::not_evaluated, v.reason, {}});
    }
    const double finest = levels.front().l2;
    record("psi L2 finest, " + f.name, finest, 0.0, finest, "relative <= 1e-4",
           judge(r.info, {}, finest <= l2_tolerance, "relative L2 <= 1e-4"));
}

}  // namespace

BenchmarkResult run_tvb04(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-04";
    r.title = "Poisson manufactured solution with a dielectric jump";
    r.info = info();
    const Materials mat = materials();
    if (const std::string bad = check_derivatives(mat.beta); !bad.empty()) {
        r.run_error = "reference self-check: " + bad;
        return r;
    }
    // The contact values at zero bias, from the coarsest 1D device (they do not depend on the mesh).
    const device::Device d0 = make_device(*mesh::make_tensor_grid(x_axis(15)));
    const auto eq = solve::solve_bias(d0, std::vector<double>{0.0, 0.0}, options());
    if (!eq) {
        r.run_error = "equilibrium: " + eq.error().message;
        return r;
    }
    const std::size_t gate = static_cast<std::size_t>(d0.contacts()[0].nodes.front());
    const std::size_t sub = static_cast<std::size_t>(d0.contacts()[1].nodes.front());
    Contacts c{eq->fields.potential_V[sub], eq->fields.n_cm3[sub], eq->fields.p_cm3[sub], 0.0};
    c.gate_bias_V = c.psi_c - eq->fields.potential_V[gate];  // the electrode's potential then equals psi_c
    std::vector<Family> families{{"1D tensor", 1, {127, 63, 31, 15}}};
    if (tier >= Tier::full) {
        families.push_back({"2D tensor", 2, {63, 31, 15, 7}});
        families.push_back({"3D tensor", 3, {31, 15, 7, 3}, true});  // CV-Q13: in full, with PARDISO
    }
    for (const Family& f : families) {
        judge_family(r, f, c, mat);
        if (r.run_error) return r;
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
