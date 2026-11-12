// TVB-24, band-to-band tunnelling: the generated pairs of the local and nonlocal Kane rates and of the direct-gap
// WKB rate in a reverse-biased p-i-n (DECISIONS Unit 55 6, 17.2 (m)-(o), 17.4, 17.13 (e), (h); criteria 17.8, frozen
// before code; the construction, the QE predicate and the pre-registered finding S2 in 17.24, before the first run).
//
// Code under test: the drift-diffusion system with band-to-band generation (models.btbt_local; models.btbt_nonlocal
// kane and direct_wkb with the traced paths, tunnel_paths.hpp), solve_bias / sweep_bias, the stored rates (Unit 45)
// and paths (results::TunnelPath).
// Construction, fixed before any run (17.24):
//   - 1D, one material, an abrupt p+ i n+ diode: N_A = 1e19 cm^-3 on x < 0, undoped on (0, d), N_D = 1e19 on x > d,
//     d = 10 nm; ohmic contacts at both ends (the anode, contact 0, on the p side, at the bias; the cathode at 0 V);
//     300 K; Boltzmann, no band-gap narrowing, no recombination (SRH, Auger, radiative off), constant mobility (the
//     doping dependence off); the weighted terminal current (Unit 32), the minority floor (Unit 50) and the stored
//     rates on.
//   - mesh: uniform, spacing h = d / m, nodes at x_j = (j + 1/2) h for j = -J - 1 ... m + J, so that x = 0 and x = d
//     are edge midpoints (the boxes' boundaries) at every level; m = 18 2^k, J = 54 2^k (each side 3 d = 30 nm, about
//     13 nm beyond the widest depletion layer); k = 3, 2, 1, 0 (N = 1010, 506, 254, 128 nodes), and k = 4 (N = 2018)
//     in the factorial. Tensor cells (the nonlocal paths need them).
//   - Three cases, two reverse biases each, reached by 0.5 V steps from 0 V:
//     local Kane on silicon (NiTCAD's) at -3 and -6 V; nonlocal Kane on silicon at -3 and -6 V; direct WKB on GaAs
//     (NiTCAD's set) with the spec's tunnelling masses m_c = 0.067, m_v = 0.5 (17.2 (o); the masses' provenance is not
//     at issue at this level) at -2 and -4 V.
// Reference (verify/reference/pin_junction.hpp, tunnel_generation.hpp): the QE p-i-n potential with NiTCAD's q, eps,
// V_T, n_i, gap and (A, B) as inputs (I1), and on it the local pairs int G(E) dx, the nonlocal Kane pairs and paths,
// and the WKB pairs and paths with the flat field 1e3 V/cm (tunnel_paths.cpp's `flat_field_V_cm`).
// The simulated pairs P are -sum_i band_to_band_p[i] V_i (the hole vector; a generation is negative).
// E (17.8), per case and bias:
//   - (a) p-hat in [1.8, 2.3] of the relative error of P on k = 3 ... 0;
//   - (b) at k = 3, the relative error of P <= 1e-3;
//   - nonlocal cases, at k = 3: per path, l and F relative <= 1e-3 for every path that reached its crossing with a
//     rate >= 1e-6 of the largest (judged on the largest of the errors), against the reference's path from the level
//     at the start node's position;
//   - the current, at k = 3: |J - q P_ref| / (q P_ref) <= 1e-3, J the anode's weighted current, judged where |J| >= 1e3
//     times its resolution (else NOT EVALUATED).
//   Applicability (17.24): P_QE = q sum_i |Delta rho_i| V_i / (eps E_max) <= 1e-6 on every level, Delta rho the
//   difference of p - n from its value with the electrons at the cathode's level and the holes at the anode's (the
//   reference's QE model), E_max the reference's largest field; else NOT APPLICABLE. A comparison whose U_ref or U_iter
//   exceeds a tenth of its tolerance is NOT EVALUATED.
//   U_ref: the reference's own error estimate and its difference from the second method. U_iter: |P(k = 3, tol / 100)
//   - P(k = 3, tol)| / P_ref, the Newton tolerance's share, measured at the finest level and used at every level.
//   U_disc: the three-grid GCI (F_s 1.25) on k = 3, 2, 1 at the observed order. The (h, tol) factorial (16.2 (p)) on
//   k = 3, 4 and tol = 1e-9, 1e-11, per pairs comparison.
// S2 (17.24, pre-registered): NiTCAD's WKB slope |dE_v/dx| is the forward difference over the path's first segment,
// first order where E' != 0, so the WKB pairs carry a leading relative error c_1 h (c_1 from the reference:
// wkb_slope_coefficient). The criteria stand (IG-Q25 A); the predicted c_1 h is printed beside each level's error.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/pin_junction.hpp"
#include "NiTCAD/verify/reference/tunnel_generation.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace ref = reference;

constexpr double T = 300.0;
constexpr double d_cm = 1e-6, doping_cm3 = 1e19;
constexpr int cells_i = 18, cells_side = 54;  // at k = 0; each doubles with k
constexpr double tolerance = 1e-3, tol_newton = 1e-9;
constexpr double limit_qe = 1e-6, resolution_factor = 1e3, path_share = 1e-6;
constexpr double step_V = 0.5, flat_field = 1e3;
constexpr std::array<int, 4> ladder{3, 2, 1, 0};  // finest first
constexpr int half_k = 4;                          // the factorial's h / 2
constexpr double gaas_mc = 0.067, gaas_mv = 0.5;

enum class Model { local, kane, wkb };
struct Case {
    Model model;
    const char* name;
    bool gaas;
    std::array<double, 2> biases;
};
constexpr std::array<Case, 3> cases{{{Model::local, "local Kane, Si", false, {-3.0, -6.0}},
                                     {Model::kane, "nonlocal Kane, Si", false, {-3.0, -6.0}},
                                     {Model::wkb, "direct WKB, GaAs", true, {-2.0, -4.0}}}};

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-24";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = "The QE p-i-n first integral (DECISIONS Unit 55 17.13 (e)) and on it Kane's rate A F^2 e^(-B / F), "
                 "local and along paths, and the direct-gap WKB path rate of Esseni et al. (2017) eqs (9), (11), (12) "
                 "as TVR-3's model_generation";
    r.version = "DECISIONS Unit 55 17.3, 17.8, 17.13 (e), (h), 17.24";
    r.equations = {1, 2, 3};
    r.form = "P = int G(E(x)) dx (local); P = int G(E_g / l(x_a)) dx_a with l the start's distance to E_c = E_v(x_a) "
             "(nonlocal Kane); P = int G_WKB(I_k, I_ik, |dE_v/dx|, k_m^2) dx_a (direct WKB)";
    r.assumptions = {"the carriers at the contacts' levels (QE; P_QE measured)", "one material, no band-gap narrowing",
                     "300 K", "the spec's GaAs tunnelling masses 0.067 and 0.5"};
    r.derivation = "verify/reference/pin_junction.hpp, verify/reference/tunnel_generation.hpp";
    r.verification = "E: orders and finest-level pairs per case and bias; per-path l and F (nonlocal); the current "
                     "where resolved";
    r.second_method = "three-piece Chebyshev collocation (the potential, L0); tanh-sinh for the pairs and lengths, "
                      "Kronrod-15 in a smoothing substitution for the WKB integrals (at run time)";
    return r;
}

double spacing(int k) { return d_cm / static_cast<double>(cells_i << k); }

std::vector<double> axis(int k) {
    const int m = cells_i << k, J = cells_side << k;
    const double h = spacing(k);
    std::vector<double> x;
    for (int j = -J - 1; j <= m + J; ++j) x.push_back((j + 0.5) * h);
    return x;
}

std::expected<device::Device, std::string> diode(int k, const physics::Semiconductor& material) {
    const std::vector<double> x = axis(k);
    auto m = mesh::make_tensor_grid(x);
    if (!m) return std::unexpected(m.error().message);
    auto cells = mesh::TensorCells::create(*m, x);
    if (!cells) return std::unexpected(cells.error().message);
    const std::size_t n = m->node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double xi = m->points()[i][0];
        if (xi < 0.0) acceptors[i] = doping_cm3;
        if (xi > d_cm) donors[i] = doping_cm3;
    }
    auto anode = m->find_boundary("x_min")->nodes;
    auto cathode = m->find_boundary("x_max")->nodes;
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = T,
                                     .regions = {{"semiconductor", material}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
                                     .cells = std::move(*cells)});
    if (!d) return std::unexpected(d.error().message);
    return std::move(*d);
}

solve::BiasOptions options(Model model, double tol) {
    solve::BiasOptions o;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.doping_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.btbt_local = model == Model::local;
    o.models.btbt_nonlocal = model == Model::kane  ? assemble::NonlocalTunnelling::kane
                             : model == Model::wkb ? assemble::NonlocalTunnelling::direct_wkb
                                                   : assemble::NonlocalTunnelling::off;
    o.weighted_current = true;
    o.minority_floor = true;
    o.store_rates = true;
    o.newton.tol_update = tol;
    return o;
}

// The reference at one bias.
struct RefPoint {
    double V = 0.0;
    ref::PinJunction j;
    double pairs = 0.0, u_ref = 0.0, E_max = 0.0;
    double c1 = 0.0;  // WKB: the slope coefficient [1/cm] (S2)
};

// One bias's values on one level.
struct Sample {
    double pairs = 0.0, path_pairs = 0.0, J = 0.0, J_res = 0.0, qe = 0.0;
    std::vector<results::TunnelPath> paths;
};

struct Level {
    int k = 0;
    double h = 0.0;
    std::vector<double> x;  // the node positions
    std::array<Sample, 2> s{};
};

// P_QE (17.24): q sum |Delta rho| V / (eps E_max), Delta rho = (p - n) - (p_QE - n_QE) with the electrons at the
// cathode's Fermi level (0 eV) and the holes at the anode's (-V_anode eV; the band diagram's convention).
double qe_share(const results::BiasPoint& p, const std::vector<double>& volume, double V_anode, double VT, double eps,
                double E_max) {
    double sum = 0.0;
    for (std::size_t i = 0; i < volume.size(); ++i) {
        const double n = p.fields.n_cm3[i], h = p.fields.p_cm3[i];
        const double nq = n * std::exp(-p.bands.electron_fermi_eV[i] / VT);
        const double hq = h * std::exp((p.bands.hole_fermi_eV[i] + V_anode) / VT);
        sum += std::abs((h - n) - (hq - nq)) * volume[i];
    }
    return base::q_C * sum / (eps * E_max);
}

}  // namespace

BenchmarkResult run_tvb24(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-24";
    r.title = "Band-to-band tunnelling: local and nonlocal Kane, direct WKB";
    r.info = info();
    if (tier < Tier::full) return r;
    const double VT = base::thermal_voltage(T);
    physics::SemiconductorParameters gp = physics::gallium_arsenide_parameters;
    gp.band_to_band.electron_mass = gaas_mc;
    gp.band_to_band.hole_mass = gaas_mv;
    const auto gaas = physics::Semiconductor::create(gp);
    if (!gaas) {
        r.run_error = "GaAs with the spec's masses: " + gaas.error().message;
        return r;
    }
    const physics::Semiconductor si = physics::silicon();

    const auto record = [&](const std::string& quantity, const std::string& unit, double parameter, double sim,
                            double refv, double uref, double metric, const std::string& crit) -> ComparisonRecord& {
        r.comparisons.push_back(ComparisonRecord{quantity, unit, parameter, sim, refv, uref, metric, crit, {}});
        return r.comparisons.back();
    };

    for (const Case& c : cases) {
        const physics::Semiconductor& material = c.gaas ? *gaas : si;
        const physics::SemiconductorParameters& mp = material.parameters();
        const double eps = mp.eps_r * base::eps0_F_per_cm;
        const double ni = physics::intrinsic_density(material, T);
        const double gap = physics::band_gap_eV(material, T);
        const ref::KaneInput kane{mp.band_to_band.A_per_cm3_s, mp.band_to_band.B_V_per_cm};
        const ref::WkbInput wkb{gap, mp.band_to_band.electron_mass, mp.band_to_band.hole_mass, flat_field};
        const ref::PinInput in{base::q_C, eps, VT, ni, doping_cm3, 0.0, doping_cm3, d_cm};
        r.notes.push_back(std::string(c.name) + ": spec (NiTCAD's material data, I1): eps_r " + shortest(mp.eps_r) +
                          ", n_i " + shortest(ni) + " cm^-3, E_g " + shortest(gap) + " eV" +
                          (c.model == Model::wkb ? ", m_c " + shortest(wkb.electron_mass) + ", m_v " +
                                                       shortest(wkb.hole_mass)
                                                 : ", A " + shortest(kane.A_per_cm3_s) + " cm^-3 s^-1, B " +
                                                       shortest(kane.B_V_per_cm) + " V/cm"));

        // The reference at both biases.
        std::array<RefPoint, 2> refs{};
        for (std::size_t b = 0; b < 2; ++b) {
            const double V = c.biases[b];
            const std::string at = std::string(c.name) + " at " + shortest(V) + " V";
            const auto j = ref::pin_junction(in, V);
            if (!j) {
                r.run_error = "reference: " + at + ": " + j.error();
                return r;
            }
            const auto P = c.model == Model::local ? ref::local_kane_pairs(*j, kane)
                           : c.model == Model::kane ? ref::nonlocal_kane_pairs(*j, gap, kane)
                                                    : ref::wkb_pairs(*j, wkb);
            if (!P) {
                r.run_error = "reference pairs: " + at + ": " + P.error();
                return r;
            }
            RefPoint& rp = refs[b];
            rp.V = V;
            rp.j = *j;
            rp.pairs = P->pairs_cm2_s.to_double();
            rp.u_ref = (std::abs((P->second_cm2_s - P->pairs_cm2_s).to_double()) + P->abs_error) / rp.pairs;
            rp.E_max = std::max(j->field_p_edge_V_per_cm.to_double(), j->field_n_edge_V_per_cm.to_double());
            if (c.model == Model::wkb) {
                const auto c1 = ref::wkb_slope_coefficient(*j, wkb);
                if (!c1) {
                    r.run_error = "reference slope coefficient: " + at + ": " + c1.error();
                    return r;
                }
                rp.c1 = c1->to_double();
            }
            r.notes.push_back("reference " + at + ": pairs " + shortest(rp.pairs) + " cm^-2 s^-1 (q P " +
                              shortest(base::q_C * rp.pairs) + " A/cm^2), U_ref " + shortest(rp.u_ref) + ", E_max " +
                              shortest(rp.E_max) + " V/cm, W_p " + shortest(j->depletion_p_cm.to_double()) +
                              " cm, matching residual " + shortest(j->matching_residual) +
                              (c.model == Model::wkb ? ", c_1 " + shortest(rp.c1) + " /cm" : ""));
        }

        // The ramp: 0.5 V steps to each bias in turn.
        std::vector<std::vector<double>> ramp;
        std::array<std::size_t, 2> at{};
        {
            double v = 0.0;
            for (std::size_t b = 0; b < 2; ++b) {
                while (v > c.biases[b] + 1e-12) {
                    v = std::max(v - step_V, c.biases[b]);
                    ramp.push_back({v, 0.0});
                }
                at[b] = ramp.size() - 1;
            }
        }

        const auto run_level = [&](int k, double tol, bool keep_paths) -> std::expected<Level, std::string> {
            const std::string where = std::string(c.name) + ", k = " + std::to_string(k) + ", tol " + shortest(tol);
            const auto d = diode(k, material);
            if (!d) return std::unexpected(where + ": " + d.error());
            const auto sw = solve::sweep_bias(*d, ramp, options(c.model, tol));
            if (!sw) return std::unexpected(where + ": " + sw.error().message);
            if (sw->stopped) return std::unexpected(where + ": " + sw->stopped->message);
            Level lv;
            lv.k = k;
            lv.h = spacing(k);
            for (const auto& p : d->mesh().points()) lv.x.push_back(p[0]);
            const std::vector<double> volume(d->mesh().volumes().begin(), d->mesh().volumes().end());
            for (std::size_t b = 0; b < 2; ++b) {
                const results::BiasPoint& p = sw->points[at[b]];
                if (p.rates.band_to_band_p.size() != volume.size()) {
                    return std::unexpected(where + ": the band-to-band rates were not stored");
                }
                Sample& s = lv.s[b];
                for (std::size_t i = 0; i < volume.size(); ++i) s.pairs -= p.rates.band_to_band_p[i] * volume[i];
                for (const auto& t : p.tunnel_paths) s.path_pairs += t.generation_cm3_s * volume[t.start_node];
                s.J = std::abs(p.weighted_current[0]);
                s.J_res = p.weighted_current_resolution[0];
                s.qe = qe_share(p, volume, c.biases[b], VT, eps, refs[b].E_max);
                if (keep_paths) s.paths = p.tunnel_paths;
            }
            return lv;
        };

        std::vector<Level> levels;  // finest first
        for (const int k : ladder) {
            auto lv = run_level(k, tol_newton, k == ladder.front());
            if (!lv) {
                r.run_error = lv.error();
                return r;
            }
            levels.push_back(std::move(*lv));
        }
        std::array<Level, 3> fac;  // (h, tol / 100), (h / 2, tol), (h / 2, tol / 100)
        {
            std::size_t n = 0;
            for (const auto& [k, tol] : {std::pair{ladder.front(), 0.01 * tol_newton}, std::pair{half_k, tol_newton},
                                         std::pair{half_k, 0.01 * tol_newton}}) {
                auto lv = run_level(k, tol, false);
                if (!lv) {
                    r.run_error = lv.error() + " (factorial)";
                    return r;
                }
                fac[n++] = std::move(*lv);
            }
        }

        for (std::size_t b = 0; b < 2; ++b) {
            const RefPoint& rp = refs[b];
            const std::string at_b = std::string(c.name) + ", " + shortest(rp.V) + " V";
            const auto signed_error = [&](const Level& l) { return (l.s[b].pairs - rp.pairs) / rp.pairs; };
            const double u_iter = std::abs(fac[0].s[b].pairs - levels.front().s[b].pairs) / rp.pairs;
            double qe = 0.0;
            std::vector<double> hs, es, us;
            for (const Level& l : levels) {
                hs.push_back(l.h);
                es.push_back(std::abs(signed_error(l)));
                us.push_back(u_iter);
                qe = std::max(qe, l.s[b].qe);
                std::string line = at_b + ", k = " + std::to_string(l.k) + " (h " + shortest(l.h) + " cm): pairs " +
                                   shortest(l.s[b].pairs) + ", error " + shortest(signed_error(l));
                if (c.model != Model::local) line += ", paths' pairs " + shortest(l.s[b].path_pairs);
                if (c.model == Model::wkb) line += ", S2's c_1 h " + shortest(rp.c1 * l.h);
                line += ", J " + shortest(l.s[b].J) + " (res " + shortest(l.s[b].J_res) + "), P_QE " +
                        shortest(l.s[b].qe);
                r.notes.push_back(line);
            }
            for (const Level& l : fac) qe = std::max(qe, l.s[b].qe);
            const study::OrderStudy st = study::order_study(hs, es, us);
            const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
            std::string orders;
            for (const double p : st.orders) orders += " " + shortest(p);
            const double u_disc = ref::gci_safety_factor_three_grids * std::abs(es[1] - es[0]) /
                                  (std::pow(hs[1] / hs[0], st.observed) - 1.0);
            const std::vector<double> parts{u_disc, u_iter, rp.u_ref};
            const Factorial fx = factorial(signed_error(levels.front()), signed_error(fac[0]), signed_error(fac[1]),
                                           signed_error(fac[2]), parts);
            r.notes.push_back(at_b + ": orders (finest pair first):" + orders + "; spread " +
                              shortest(st.order_spread) + "; " + fx.note + " (U_disc " + shortest(u_disc) + ")");
            std::optional<Comparison> blocked;
            if (qe > limit_qe) {
                blocked = Comparison{Outcome::not_applicable, "P_QE " + shortest(qe) + " above 1e-6", {}};
            } else if (rp.u_ref > 0.1 * tolerance) {
                blocked = Comparison{Outcome::not_evaluated,
                                     "U_ref " + shortest(rp.u_ref) + " above a tenth of the tolerance", {}};
            } else if (u_iter > 0.1 * tolerance) {
                blocked = Comparison{Outcome::not_evaluated,
                                     "U_iter " + shortest(u_iter) + " above a tenth of the tolerance", {}};
            }
            const std::string band = "p-hat in [1.8, 2.3] (p_th 2)";
            {
                ComparisonRecord& rec =
                    record("pairs order, " + at_b, "1", rp.V, st.observed, 2.0, rp.u_ref, st.observed, band);
                rec.comparison = blocked    ? *blocked
                                 : v.judged ? judge(r.info, {}, v.pass, band)
                                            : Comparison{Outcome::not_evaluated, v.reason, {}};
                block_if_interacting(rec, fx);
            }
            const Level& fine = levels.front();
            {
                const std::string crit = "pairs relative <= 1e-3";
                ComparisonRecord& rec = record("pairs finest (N = " + std::to_string(fine.x.size()) + "), " + at_b,
                                               "cm^-2 s^-1", rp.V, fine.s[b].pairs, rp.pairs, rp.u_ref, es[0], crit);
                rec.comparison = blocked ? *blocked : judge(r.info, {}, es[0] <= tolerance, crit);
                block_if_interacting(rec, fx);
            }
            if (c.model != Model::local) {
                // Per path on the finest level: the reference's path from the level at the start node's position.
                double G_max = 0.0;
                for (const auto& t : fine.s[b].paths) G_max = std::max(G_max, t.generation_cm3_s);
                double worst_l = 0.0, worst_F = 0.0, worst_u = 0.0;
                std::size_t judged = 0;
                std::string missing;
                for (const auto& t : fine.s[b].paths) {
                    if (!t.reached || !(t.generation_cm3_s >= path_share * G_max)) continue;
                    const double x = fine.x[t.start_node];
                    const auto s = ref::pin_level(rp.j, x);
                    if (!s) {
                        r.run_error = "reference level at " + shortest(x) + " cm: " + s.error();
                        return r;
                    }
                    const auto p = c.model == Model::kane ? ref::kane_path(rp.j, gap, kane, *s)
                                                          : ref::wkb_path(rp.j, wkb, *s);
                    if (!p) {
                        r.run_error = "reference path at " + shortest(x) + " cm: " + p.error();
                        return r;
                    }
                    ++judged;
                    if (!p->exists) {
                        worst_l = worst_F = INFINITY;
                        missing = "; the reference has no crossing from " + shortest(x) + " cm";
                        continue;
                    }
                    const double l = p->length_cm.to_double(), F = p->field_V_per_cm.to_double();
                    worst_l = std::max(worst_l, std::abs(t.length_cm - l) / l);
                    worst_F = std::max(worst_F, std::abs(t.field_V_per_cm - F) / F);
                    worst_u = std::max(worst_u, std::abs((p->second_length_cm - p->length_cm).to_double()) / l);
                }
                r.notes.push_back(at_b + ": " + std::to_string(judged) + " of " +
                                  std::to_string(fine.s[b].paths.size()) + " paths judged (rate >= 1e-6 of " +
                                  shortest(G_max) + " cm^-3 s^-1); worst l " + shortest(worst_l) + ", worst F " +
                                  shortest(worst_F) + ", reference length uncertainty " + shortest(worst_u) + missing);
                std::optional<Comparison> path_blocked;
                if (fine.s[b].qe > limit_qe) {
                    path_blocked = Comparison{Outcome::not_applicable, "P_QE " + shortest(fine.s[b].qe) + " above 1e-6", {}};
                } else if (judged == 0) {
                    path_blocked = Comparison{Outcome::not_evaluated, "no path to judge", {}};
                } else if (worst_u > 0.1 * tolerance) {
                    path_blocked = Comparison{Outcome::not_evaluated,
                                              "the reference's length uncertainty above a tenth of the tolerance", {}};
                }
                for (const bool length : {true, false}) {
                    const std::string crit = std::string(length ? "l" : "F") + " relative <= 1e-3 on every path";
                    const double worst = length ? worst_l : worst_F;
                    ComparisonRecord& rec = record(std::string("path ") + (length ? "l" : "F") + " finest, worst of " +
                                                       std::to_string(judged) + ", " + at_b,
                                                   "1", rp.V, worst, 0.0, worst_u, worst, crit);
                    rec.comparison = path_blocked ? *path_blocked : judge(r.info, {}, worst <= tolerance, crit);
                }
            }
            {
                const double Jref = base::q_C * rp.pairs;
                const double e = std::abs(fine.s[b].J - Jref) / Jref;
                const std::string crit = "|J - q P_ref| / (q P_ref) <= 1e-3";
                ComparisonRecord& rec = record("current finest, " + at_b, "A/cm^2", rp.V, fine.s[b].J, Jref,
                                               rp.u_ref, e, crit);
                if (fine.s[b].qe > limit_qe) {
                    rec.comparison = Comparison{Outcome::not_applicable, "P_QE above 1e-6", {}};
                } else if (fine.s[b].J < resolution_factor * fine.s[b].J_res) {
                    rec.comparison = Comparison{Outcome::not_evaluated,
                                                "J below 1e3 times its resolution (resolution / J " +
                                                    shortest(fine.s[b].J_res / fine.s[b].J) + ")",
                                                {}};
                } else if (rp.u_ref > 0.1 * tolerance) {
                    rec.comparison = Comparison{Outcome::not_evaluated, "U_ref above a tenth of the tolerance", {}};
                } else {
                    rec.comparison = judge(r.info, {}, e <= tolerance, crit);
                }
            }
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
