// TVB-07, the equilibrium abrupt pn junction (DECISIONS Unit 55 6, 14.4; criteria 14.7).
//
// Code under test: the drift-diffusion system at thermal equilibrium (the box-method Poisson rows with nodal doping,
// the Boltzmann densities, the ohmic contacts' neutral potentials) and solve_bias. Construction: silicon on
// [-L, L], L = 1 um, N_A = 1e16 cm^-3 for x < 0 and N_D = 2e16 cm^-3 for x > 0, ohmic contacts at both ends, no
// bias; Boltzmann statistics, no band-gap narrowing. Nodes are uniform, N of them with N even, so x = 0 is the
// midpoint of the middle edge: the doping step lies on a dual-cell face, where the box method puts it exactly.
// Reference (verify/reference/pn_junction.hpp): the infinite-domain solution by the first integral (Brent for the
// junction potential, tanh-sinh for x(psi), in double-double), cross-checked at run time against Chebyshev
// collocation of the finite-domain problem (I3).
// Criteria (frozen 14.7):
// - D: V_bi = psi(x = L) - psi(x = -L), the contacts' potential difference, against V_T [asinh(N_D / (2 n_i)) +
//   asinh(N_A / (2 n_i))] (the exact neutral levels, in double-double) relative <= 64 eps, at every level;
// - E (a): psi L-infinity (relative to V_bi) and E_max orders p-hat in [1.8, 2.3]; E_max is the middle edge's
//   (psi_{N/2} - psi_{N/2-1}) / h against |psi'(0)|;
// - E (b): on the finest level psi L-infinity relative <= 1e-4 and E_max relative <= 1e-3.
// Ladder: N = 64, 128, 256, 512 (full); fast runs the coarsest two with V_bi only, the orders not judged (14.7).
// Reference uncertainty: the domain truncation V_T t at the contacts (t the reduced tail of the infinite-domain
// solution there), which must be at most 1e-3 of the psi criterion (TB-Q4), else the E comparisons are
// not_evaluated; the collocation agrees with the first integral to within that plus its own tail estimate.
// Reference independence I1: q, eps_si, V_T and n_i are NiTCAD's material data (input); the rest is analytic.
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
#include "NiTCAD/verify/reference/pn_junction.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double L = 1e-4;  // cm, each side
constexpr double acceptors_cm3 = 1e16, donors_cm3 = 2e16;
constexpr double eps = 0x1p-52;
constexpr double psi_tolerance = 1e-4, field_tolerance = 1e-3;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-07";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PA;
    r.citation = "First integral of the 1D Poisson-Boltzmann equation for an abrupt junction (e.g. Sze and Ng, "
                 "Physics of Semiconductor Devices, 3rd ed., 2007, 2.2; exact, not the depletion approximation)";
    r.version = "DECISIONS Unit 55 14.4, 14.7";
    r.equations = {1, 4};
    r.form = "(1/2) u'^2 = (q / (eps V_T)) [n_b G(u - u_b) + p_b G(u_b - u)], G(t) = e^t - 1 - t";
    r.assumptions = {"thermal equilibrium", "Boltzmann statistics", "abrupt step doping on a dual-cell face",
                     "constant permittivity", "infinite-domain reference (truncation measured)"};
    r.derivation = "verify/reference/pn_junction.hpp";
    r.verification = "D: V_bi to rounding; E: psi and E_max orders and finest-level accuracy";
    r.second_method = "Chebyshev collocation of the finite-domain problem (I3), at run time";
    return r;
}

std::vector<double> nodes(int N) {  // x_k = L (2k - (N - 1)) / (N - 1): antisymmetric, 0 at the middle edge's midpoint
    std::vector<double> x;
    for (int k = 0; k < N; ++k) x.push_back(L * (2.0 * k - (N - 1)) / (N - 1));
    return x;
}

solve::BiasOptions options() {
    solve::BiasOptions o;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.newton.tol_update = 1e-12;
    return o;
}

struct Level {
    int N;
    double h, vbi, psi_linf, field, u_iter, tail;
};

}  // namespace

BenchmarkResult run_tvb07(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-07";
    r.title = "Equilibrium abrupt pn junction";
    r.info = info();
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(300.0);
    reference::PnJunctionInput in{.q_C = base::q_C,
                                  .eps_F_per_cm = si.parameters().eps_r * base::eps0_F_per_cm,
                                  .thermal_voltage_V = VT,
                                  .ni_cm3 = physics::intrinsic_density(si, 300.0),
                                  .acceptors_cm3 = acceptors_cm3,
                                  .donors_cm3 = donors_cm3};
    const auto ref = reference::pn_junction(in);
    if (!ref) {
        r.run_error = "reference: " + ref.error();
        return r;
    }
    const double vbi_ref = ref->built_in_V.to_double();
    const double field_ref = ref->field_max_V_per_cm.to_double();
    r.notes.push_back("reference: V_bi " + shortest(vbi_ref) + " V, psi_j " + shortest(ref->psi_j.to_double()) +
                      " V, E_max " + shortest(field_ref) + " V/cm, depletion edges x_p " +
                      shortest(ref->depletion_p_cm.to_double()) + " cm and x_n " +
                      shortest(ref->depletion_n_cm.to_double()) + " cm (eps E_max / (q N)), Debye lengths " +
                      shortest(ref->debye_p_cm.to_double()) + " and " + shortest(ref->debye_n_cm.to_double()) +
                      " cm, Brent " + std::to_string(ref->brent_evaluations) + " evaluations");
    // The I3 cross-check on the coarsest grid's nodes.
    {
        const std::vector<double> x = nodes(64);
        const auto fi = reference::pn_potential(*ref, x);
        const auto col = reference::pn_collocation(*ref, -L, L, 4, 32, x);
        if (!fi || !col) {
            r.run_error = "reference cross-check: " + (!fi ? fi.error() : col.error());
            return r;
        }
        double worst = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            worst = std::max(worst, std::abs((*fi)[i].psi_V.to_double() - col->psi_V[i]));
        }
        const double contact_tail = VT * std::max((*fi)[0].tail.to_double(), (*fi)[x.size() - 1].tail.to_double());
        const double allowed = contact_tail + col->tail_estimate_V + 1e-10 * vbi_ref;
        r.notes.push_back("I3: collocation (4 + 4 pieces, degree 32, " + std::to_string(col->newton_iterations) +
                          " Newton steps, tail " + shortest(col->tail_estimate_V) + " V) against the first integral: " +
                          "max |difference| " + shortest(worst) + " V; truncation V_T t at the contacts " +
                          shortest(contact_tail) + " V");
        if (!(worst <= allowed)) {
            r.run_error = "reference cross-check: collocation and first integral differ by " + shortest(worst) +
                          " V, more than " + shortest(allowed) + " V";
            return r;
        }
    }
    std::vector<int> ladder{64, 128};
    if (tier >= Tier::full) ladder = {64, 128, 256, 512};
    std::vector<Level> levels;
    double u_ref = 0.0;  // the reference's uncertainty, relative to V_bi (truncation)
    for (const int N : ladder) {
        const std::vector<double> x = nodes(N);
        auto m = mesh::make_tensor_grid(x);
        if (!m) {
            r.run_error = m.error().message;
            return r;
        }
        const std::size_t n = m->node_count();
        std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (m->points()[i][0] < 0.0) {
                acceptors[i] = acceptors_cm3;
            } else {
                donors[i] = donors_cm3;
            }
        }
        auto left = m->find_boundary("x_min")->nodes;
        auto right = m->find_boundary("x_max")->nodes;
        const auto d = device::Device::create({.mesh = std::move(*m),
                                               .temperature_K = 300.0,
                                               .regions = {{"silicon", si}},
                                               .node_region = std::vector<device::RegionId>(n, 0),
                                               .donors = std::move(donors),
                                               .acceptors = std::move(acceptors),
                                               .contacts = {{"anode", device::ContactKind::ohmic, std::move(left)},
                                                            {"cathode", device::ContactKind::ohmic, std::move(right)}}});
        if (!d) {
            r.run_error = "N = " + std::to_string(N) + ": " + d.error().message;
            return r;
        }
        const auto s = solve::solve_bias(*d, std::vector<double>{0.0, 0.0}, options());
        if (!s) {
            r.run_error = "N = " + std::to_string(N) + ": the solve failed: " + s.error().message;
            return r;
        }
        const auto exact = reference::pn_potential(*ref, x);
        if (!exact) {
            r.run_error = "reference: " + exact.error();
            return r;
        }
        const auto& psi = s->fields.potential_V;
        double worst = 0.0;
        for (std::size_t i = 0; i < n; ++i) worst = std::max(worst, std::abs(psi[i] - (*exact)[i].psi_V.to_double()));
        const double h = x[1] - x[0];
        const std::size_t a = static_cast<std::size_t>(N / 2 - 1), b = a + 1;
        const double field = (psi[b] - psi[a]) / (x[b] - x[a]);
        const double tail = VT * std::max((*exact)[0].tail.to_double(), (*exact)[n - 1].tail.to_double()) / vbi_ref;
        u_ref = std::max(u_ref, tail);
        const auto& it = s->convergence.iterations;
        const double u_iter = it.empty() ? 0.0 : it.back().update * VT / vbi_ref;
        levels.push_back({N, h, psi[n - 1] - psi[0], worst / vbi_ref, std::abs(field - field_ref) / field_ref, u_iter,
                          tail});
        r.notes.push_back("N = " + std::to_string(N) + ": V_bi " + shortest(psi[n - 1] - psi[0]) + " V, psi Linf " +
                          shortest(worst / vbi_ref) + ", E_max " + shortest(field) + " V/cm (relative error " +
                          shortest(levels.back().field) + "), U_iter " + shortest(u_iter) + ", truncation " +
                          shortest(tail));
    }
    const auto record = [&](const std::string& quantity, double hh, double sim, double refv, double metric,
                            const std::string& crit, Comparison cmp) {
        ComparisonRecord rec{quantity, "1", hh, sim, refv, 0.0, metric, crit, {}};
        rec.comparison = std::move(cmp);
        r.comparisons.push_back(std::move(rec));
    };
    // D: V_bi at every level.
    for (const Level& l : levels) {
        const double e = std::abs(l.vbi - vbi_ref) / vbi_ref;
        const std::string crit = "relative <= 64 eps = " + shortest(64.0 * eps);
        record("V_bi, N = " + std::to_string(l.N), l.h, l.vbi, vbi_ref, e, crit, judge(r.info, {}, e <= 64.0 * eps, crit));
    }
    if (tier < Tier::full) {
        apply_known_findings(r, known_findings());
        return r;
    }
    // E: orders and finest accuracy, unless the reference's uncertainty exceeds 1e-3 of the psi criterion (TB-Q4).
    const bool ref_ok = u_ref <= 1e-3 * psi_tolerance;
    r.notes.push_back("reference uncertainty (truncation, relative to V_bi) " + shortest(u_ref) + " against " +
                      shortest(1e-3 * psi_tolerance));
    std::vector<double> hs, u;
    for (auto l = levels.rbegin(); l != levels.rend(); ++l) {  // finest first, as order_study takes them
        hs.push_back(l->h);
        u.push_back(l->u_iter);
    }
    const std::string band = "p-hat in [1.8, 2.3]";
    for (const bool field : {false, true}) {
        std::vector<double> e;
        for (auto l = levels.rbegin(); l != levels.rend(); ++l) e.push_back(field ? l->field : l->psi_linf);
        const study::OrderStudy st = study::order_study(hs, e, u);
        const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
        std::string orders;
        for (const double p : st.orders) orders += " " + shortest(p);
        if (st.asymptotic_evaluated) orders += ", order spread " + shortest(st.order_spread);
        r.notes.push_back(std::string(field ? "E_max" : "psi Linf") + " orders:" + orders);
        Comparison c = !ref_ok      ? missing_input("a reference within 1e-3 of the criterion (truncation)")
                       : v.judged   ? judge(r.info, {}, v.pass, band)
                                    : Comparison{Outcome::not_evaluated, v.reason, {}};
        record(field ? "E_max order" : "psi Linf order", hs.front(), st.observed, 2.0, st.observed, band, std::move(c));
    }
    const Level& fin = levels.back();
    const auto finest = [&](const std::string& q, double metric, double tol) {
        const std::string crit = "relative <= " + shortest(tol);
        record(q, fin.h, metric, 0.0, metric, crit,
               ref_ok ? judge(r.info, {}, metric <= tol, crit)
                      : missing_input("a reference within 1e-3 of the criterion (truncation)"));
    };
    finest("psi Linf finest", fin.psi_linf, psi_tolerance);
    finest("E_max finest", fin.field, field_tolerance);
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
