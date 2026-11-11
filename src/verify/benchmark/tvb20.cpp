// TVB-20 (the long-channel MOSFET against Pao-Sah) and TVB-21 (its subthreshold swing) (DECISIONS Unit 55 16.4, 16.18
// (d), 16.19; construction 16.46).
//
// Code under test: the 2D drift-diffusion system with a lumped gate on a long-channel n-MOSFET: a uniform p-substrate
// (N_A = 1e17 cm^-3), abrupt n+ source and drain boxes (1e19 cm^-3, 0.5 um wide, 0.1 um deep), t_ox = 10 nm, a metal
// gate (4.5 eV) over the channel, the body contact at the bottom (1 um deep); 300 K; constant mobility (doping, field
// and surface mobility off), Boltzmann, no BGN, no generation-recombination (16.2 (i)).
// Reference: verify/reference/pao_sah (Pao and Sah's double integral, form 1, cutoff psi_B; PB) with the same
// material values NiTCAD is given (I1).
// The channel-length limit (16.18 (d), DT-Q13): NiTCAD at L = L0 {4, 8, 16, 32} (DT-Q42 A:
// L0 {1, 2, 4, 8} left U_L at 1.5 to 2 % in weak inversion), L0 = 2.5 um (above 20 W_dmax, W_dmax about 0.11 um); Q(L) = I_D L per unit width; the order of Q in 1/L from the ladder's three differences, in [0.8, 1.3]
// with spread <= 0.1; Q_inf by Richardson at the finest pair, U_L its three-grid GCI.
// Accuracy (A, 16.19): |Q_inf - Q_PS| <= U_num + U_ref + eps_approx, eps_approx = U_L + eps_cutoff, eps_cutoff =
// |Q_PS(c = 1) - Q_PS(c = 1/2)| + |Q_PS(form 1) - Q_PS(form 2)|; NOT EVALUATED when eps_approx > 1e-2 Q_PS or the
// order is not judged a pass; NOT APPLICABLE where psi_s(V_S) < psi_B + 3 V_T.
// U_num at the longest L: the (h, tol) factorial (16.2 (p)); U_disc = the two-grid GCI (F_s 3, p 2) of the h value
// against h/2 (every spacing halved and the grading ratio square-rooted), U_iter = |Q(tol) - Q(tol / 100)|.
// TVB-21: S = 2 h_G / log10(Q(V_G + h_G) / Q(V_G - h_G)), h_G = 20 mV, at V_D = 50 mV (<= 3 V_T), the same stencil on
// NiTCAD's ladder and on the reference; the same rule on S (16.19). The depletion formula ln 10 V_T (1 + C_dep / C_ox)
// is printed beside it (DT-Q6), never judged.
// Nightly tier (16.7): the factorial's h/2 mesh is about 30 000 nodes.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/benchmark/tvb32.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/mos_capacitor.hpp"
#include "NiTCAD/verify/reference/pao_sah.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double NA = 1e17, N_sd = 1e19, t_ox_cm = 1e-6, eps_ox_r = 3.9, phi_m_eV = 4.5, T = 300.0;
constexpr double Lsd = 5e-5, xj = 1e-5, depth = 1e-4, L0 = 2.5e-4;
// tol 1e-6, so the factorial's tol / 100 is 1e-8: below about 1e-9 the holes of the inverted surface, 1e-20 of the
// largest density, hold Newton's update at their rounding (16.47, DT-Q40).
constexpr double tol_newton = 1e-6, h_G = 0.02, approx_limit = 1e-2;
constexpr int levels = 4, first_level = 2;  // L = L0 2^(first_level + k), k < levels: 4 L0 ... 32 L0 (16.47, DT-Q42 A)
constexpr double drains[2] = {0.05, 0.5};
const std::vector<double> judged_gates[2] = {{0.45, 0.6, 0.8, 1.2, 1.6, 2.0}, {0.6, 0.8, 1.2, 1.6, 2.0}};
constexpr double swing_gates[3] = {0.4, 0.5, 0.6};

ReferenceInfo info(const char* id) {
    ReferenceInfo r;
    r.id = id;
    r.cls = ReferenceClass::A;
    r.independence = Independence::I1;
    r.provenance = Provenance::PB;
    r.citation = "Pao and Sah, Solid-State Electronics 9, 927 (1966), through H.-S. P. Wong (Stanford EE 316) and "
                 "Galup-Montoro, Schneider and Cunha (NSTI-Nanotech 2007) (PB, 16.17)";
    r.version = "DECISIONS Unit 55 16.2 (i), 16.17, 16.18 (d), 16.19, 16.46; DT-Q6, DT-Q13";
    r.equations = {1};
    r.form = "I_D L / W = mu int_0^V_D (-Q_I(V)) dV, Q_I = -q int_{psi_B}^{psi_s(V)} n / E dpsi";
    r.assumptions = {"gradual channel (the L -> infinity limit of NiTCAD's ladder)", "constant mobility",
                     "uniform substrate", "Boltzmann statistics", "no generation-recombination",
                     "the bulk term's ambiguity bounded (eps_cutoff)"};
    r.derivation = "verify/reference/pao_sah.hpp";
    r.verification = "pao_sah_test: V = 0 against Kingston-Neustadter, the limits, the bulk ambiguity";
    r.second_method = "g_d = mu (-Q_I(V_D)) against a fourth-order difference of Q in V_D, at run time";
    r.approximation_error = "eps_approx = U_L (the ladder's GCI) + eps_cutoff (16.18 (d))";
    return r;
}

// Spacings growing geometrically from h0 at both ends of [a, b] up to hm, scaled to fit exactly.
void segment(std::vector<double>& x, double a, double b, double h0a, double h0b, double hm, double r) {
    std::vector<double> left, right;
    double sl = 0.0, sr = 0.0;
    const double half = 0.5 * (b - a);
    for (double h = h0a; sl < half; h = std::min(h * r, hm)) {
        left.push_back(h);
        sl += h;
    }
    for (double h = h0b; sr < half; h = std::min(h * r, hm)) {
        right.push_back(h);
        sr += h;
    }
    std::vector<double> steps = left;
    steps.insert(steps.end(), right.rbegin(), right.rend());
    const double scale = (b - a) / (sl + sr);
    double p = a;
    for (std::size_t k = 0; k + 1 < steps.size(); ++k) {
        p += steps[k] * scale;
        x.push_back(p);
    }
    x.push_back(b);
}

// The x and y axes at refinement m (1: the ladder's mesh, 2: every spacing halved, the ratio square-rooted).
std::vector<double> x_axis(double L, int m) {
    const double r = std::pow(1.2, 1.0 / m), hj = 2e-7 / m;
    std::vector<double> x{0.0};
    segment(x, 0.0, Lsd, 5e-6 / m, hj, 5e-6 / m, r);
    segment(x, Lsd, Lsd + L, hj, hj, L / 40.0 / m, r);
    segment(x, Lsd + L, 2.0 * Lsd + L, hj, 5e-6 / m, 5e-6 / m, r);
    return x;
}

std::vector<double> y_axis(int m) {
    std::vector<double> y{0.0};
    segment(y, 0.0, depth, 2e-8 / m, 4e-6 / m, 4e-6 / m, std::pow(1.15, 1.0 / m));
    return y;
}

// Contacts in the order source, drain, body, gate.
device::Device mosfet(double L, int m) {
    mesh::Mesh mesh = *mesh::make_tensor_grid(x_axis(L, m), y_axis(m));
    const std::size_t n = mesh.node_count();
    std::vector<double> donors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = mesh.points()[i];
        if (p[1] <= xj && (p[0] <= Lsd || p[0] >= Lsd + L)) donors[i] = N_sd;
    }
    std::vector<mesh::NodeId> source, drain, gate;
    for (const mesh::NodeId v : mesh.find_boundary("y_min")->nodes) {
        const double px = mesh.points()[static_cast<std::size_t>(v)][0];
        if (px <= Lsd) {
            source.push_back(v);
        } else if (px >= Lsd + L) {
            drain.push_back(v);
        } else {
            gate.push_back(v);
        }
    }
    auto body = mesh.find_boundary("y_max")->nodes;
    device::Contact g{"gate", device::ContactKind::gate, std::move(gate)};
    g.gate = {.boundary = "y_min",
              .oxide_thickness_cm = t_ox_cm,
              .oxide_relative_permittivity = eps_ox_r,
              .electrode = device::GateElectrode::metal,
              .work_function_eV = phi_m_eV};
    return *device::Device::create({.mesh = std::move(mesh),
                                    .temperature_K = T,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::move(donors),
                                    .acceptors = std::vector<double>(n, NA),
                                    .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                                                 {"drain", device::ContactKind::ohmic, std::move(drain)},
                                                 {"body", device::ContactKind::ohmic, std::move(body)},
                                                 std::move(g)}});
}

solve::BiasOptions options(double tol) {
    solve::BiasOptions o;
    o.models.doping_mobility = false;
    o.models.field_mobility = false;
    o.models.electron_mobility = assemble::ElectronMobility::caughey_thomas;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.newton.tol_update = tol;
    o.step_control.max_cuts = 30;
    o.weighted_current = true;
    return o;
}

// The gate biases solved at drain bias index k: the judged ones and, at V_D = 50 mV, the swing stencils.
std::vector<double> sweep_gates(std::size_t k) {
    std::vector<double> g = judged_gates[k];
    if (k == 0) {
        for (const double v : swing_gates) {
            g.push_back(v - h_G);
            g.push_back(v + h_G);
        }
    }
    std::sort(g.begin(), g.end());
    g.erase(std::unique(g.begin(), g.end()), g.end());
    return g;
}

using Values = std::map<std::pair<std::size_t, double>, double>;  // (drain index, V_G) -> Q [A]

std::expected<Values, std::string> solve_level(double L, int m, double tol) {
    const device::Device d = mosfet(L, m);
    Values out;
    for (std::size_t k = 0; k < std::size(drains); ++k) {
        const std::vector<double> gates = sweep_gates(k);
        std::vector<std::vector<double>> points;
        for (const double vg : gates) points.push_back({0.0, drains[k], 0.0, vg});
        const auto s = solve::sweep_bias(d, points, options(tol));
        const std::string at = "L = " + shortest(L) + " cm, refinement " + std::to_string(m) + ", V_D " +
                               shortest(drains[k]) + " V";
        if (!s || s->stopped || s->points.size() != points.size()) {
            std::string history;  // the stopped point's last Newton updates: a rounding floor or a divergence
            if (s && s->unfinished) {
                const auto& it = s->unfinished->iterations;
                history = "; last updates";
                for (std::size_t i = it.size() > 6 ? it.size() - 6 : 0; i < it.size(); ++i) {
                    history += " " + shortest(it[i].update);
                }
            }
            return std::unexpected(at + ": " +
                                   (!s ? s.error().message : s->stopped ? s->stopped->message : "short sweep") +
                                   history);
        }
        for (std::size_t q = 0; q < gates.size(); ++q) {
            out[{k, gates[q]}] = std::abs(s->points[q].weighted_current[1]) * L;
        }
    }
    return out;
}

reference::PaoSahInput reference_input(int form, double cutoff) {
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T), ni = physics::intrinsic_density(si, T);
    const double cox = eps_ox_r * base::eps0_F_per_cm / t_ox_cm;
    const auto vfb = reference::flat_band_voltage({.work_function_eV = phi_m_eV,
                                                   .electron_affinity_eV = si.parameters().electron_affinity_eV,
                                                   .band_gap_eV = physics::band_gap_eV(si, T),
                                                   .conduction_dos_cm3 = physics::conduction_band_dos(si, T),
                                                   .valence_dos_cm3 = physics::valence_band_dos(si, T),
                                                   .thermal_voltage_V = VT,
                                                   .ni_cm3 = ni,
                                                   .net_doping_cm3 = -NA,
                                                   .fixed_charge_cm2 = 0.0,
                                                   .q_C = base::q_C,
                                                   .oxide_capacitance_F_per_cm2 = cox});
    return {.q_C = base::q_C,
            .eps_F_per_cm = si.parameters().eps_r * base::eps0_F_per_cm,
            .thermal_voltage_V = VT,
            .ni_cm3 = ni,
            .acceptors_cm3 = NA,
            .oxide_capacitance_F_per_cm2 = cox,
            .flat_band_V = vfb ? vfb->to_double() : std::nan(""),
            .mobility_cm2_Vs = si.parameters().electron_mobility.mu_max,
            .form = form,
            .cutoff_fraction = cutoff};
}

// One quantity's values on the ladder (4 L0 first), its factorial at the longest L, and the reference's three forms.
struct Quantity {
    std::string name;
    double ladder[levels] = {};
    double h_tol = 0.0, h_tolfine = 0.0, h2_tol = 0.0, h2_tolfine = 0.0;  // the factorial, at 32 L0
    double ref = 0.0, ref_cut = 0.0, ref_form = 0.0, u_ref = 0.0;
    bool applicable = true;
    double applicability_value = 0.0, applicability_threshold = 0.0;
};

void judge_limit(BenchmarkResult& r, const Quantity& q) {
    const double d1 = q.ladder[1] - q.ladder[0], d2 = q.ladder[2] - q.ladder[1], d3 = q.ladder[3] - q.ladder[2];
    const double p1 = std::log2(std::abs(d1 / d2)), p2 = std::log2(std::abs(d2 / d3));
    const double spread = std::abs(p2 - p1);
    const bool finite = std::isfinite(p1) && std::isfinite(p2);
    const double denom = std::pow(2.0, p2) - 1.0;
    const double q_inf = q.ladder[3] + d3 / denom;
    const double U_L = reference::gci_safety_factor_three_grids * std::abs(d3) / denom;
    // Roache's two-grid estimate of the h value (16.36's form; DT-Q41 A): F_s |f(h/2) - f(h)| r^p / (r^p - 1),
    // r = p = 2.
    const double u_disc = reference::gci_safety_factor_two_grids * std::abs(q.h2_tol - q.h_tol) * 4.0 / 3.0;
    const double u_iter = std::abs(q.h_tol - q.h_tolfine);
    const double eps_cut = std::abs(q.ref - q.ref_cut) + std::abs(q.ref - q.ref_form);
    const double eps_approx = U_L + eps_cut;
    const double bound = u_disc + u_iter + q.u_ref + eps_approx;
    const double s = std::abs(q.ref);
    const Factorial f = factorial((q.h_tol - q.ref) / s, (q.h_tolfine - q.ref) / s, (q.h2_tol - q.ref) / s,
                                  (q.h2_tolfine - q.ref) / s, std::vector<double>{u_disc / s, u_iter / s});
    std::string ladder;
    for (const double v : q.ladder) ladder += " " + shortest(v);
    r.notes.push_back(q.name + ": ladder (4 L0 ... 32 L0)" + ladder + "; orders " + shortest(p1) + " " + shortest(p2) +
                      "; Q_inf " + shortest(q_inf) + ", reference " + shortest(q.ref) + " (cutoff psi_B / 2 " +
                      shortest(q.ref_cut) + ", form 2 " + shortest(q.ref_form) + "); U_L " + shortest(U_L) +
                      ", eps_cutoff " + shortest(eps_cut) + ", U_disc " + shortest(u_disc) + ", U_iter " +
                      shortest(u_iter) + ", U_ref " + shortest(q.u_ref) + "; " + f.note);
    const std::vector<Predicate> pr{{"psi_s(V_S) >= psi_B + 3 V_T", q.applicability_value,
                                     q.applicability_threshold, q.applicable}};
    const std::string band = "order in 1/L in [0.8, 1.3], spread <= 0.1 (16.18 (d))";
    ComparisonRecord co{q.name + " order in 1/L", "1", 0.0, p2, 1.0, 0.0, p2, band, {}};
    if (!q.applicable) {
        co.comparison = judge(r.info, pr, false, band);
    } else if (!finite || spread > 0.1) {
        co.comparison = Comparison{Outcome::not_evaluated, "not in the asymptotic range (spread " + shortest(spread) +
                                                               ")", {}};
    } else {
        co.comparison = judge(r.info, pr, p2 >= 0.8 && p2 <= 1.3, band);
    }
    const std::string acc = "|Q_inf - Q_PS| <= U_num + U_ref + eps_approx (16.19)";
    ComparisonRecord ca{q.name + " limit against Pao-Sah", "A", 0.0, q_inf, q.ref, q.u_ref,
                        std::abs(q_inf - q.ref) / bound, acc, {}};
    if (!q.applicable) {
        ca.comparison = judge(r.info, pr, false, acc);
    } else if (co.comparison.outcome != Outcome::pass) {
        ca.comparison = Comparison{Outcome::not_evaluated, "Q_inf not established: the order is not a pass", {}};
    } else if (eps_approx > approx_limit * s) {
        ca.comparison = Comparison{Outcome::not_evaluated, "eps_approx " + shortest(eps_approx / s) +
                                                               " of Q_PS, above 1e-2", {}};
    } else {
        ca.comparison = judge(r.info, pr, std::abs(q_inf - q.ref) <= bound, acc);
    }
    for (ComparisonRecord* c : {&co, &ca}) {
        if (c->comparison.outcome != Outcome::not_applicable) block_if_interacting(*c, f);
        r.comparisons.push_back(std::move(*c));
    }
}

// Everything TVB-20 and TVB-21 share: NiTCAD's ladder and factorial, and the reference's three forms at each point.
struct Shared {
    std::vector<Values> ladder;          // 4 L0 ... 32 L0
    Values h_tol, h_tolfine, h2_tol, h2_tolfine;
    std::vector<std::string> notes;
    std::optional<std::string> error;
};

const Shared& shared_run() {
    static const Shared s = [] {
        Shared out;
        for (int k = 0; k < levels; ++k) {
            auto v = solve_level(L0 * std::ldexp(1.0, first_level + k), 1, tol_newton);
            if (!v) {
                out.error = v.error();
                return out;
            }
            out.ladder.push_back(std::move(*v));
        }
        const double Lmax = L0 * std::ldexp(1.0, first_level + levels - 1);
        auto a = solve_level(Lmax, 1, 0.01 * tol_newton), b = solve_level(Lmax, 2, tol_newton),
             c = solve_level(Lmax, 2, 0.01 * tol_newton);
        if (!a || !b || !c) {
            out.error = (!a ? a.error() : !b ? b.error() : c.error()) + " (factorial)";
            return out;
        }
        out.h_tol = out.ladder.back();
        out.h_tolfine = std::move(*a);
        out.h2_tol = std::move(*b);
        out.h2_tolfine = std::move(*c);
        return out;
    }();
    return s;
}

struct ReferenceValue {
    double Q = 0.0, cut = 0.0, form = 0.0, u = 0.0, psi_s = 0.0, psi_B = 0.0;
};

std::expected<ReferenceValue, std::string> reference_at(double VG, double VD) {
    const auto a = reference::pao_sah(reference_input(1, 1.0), VG, VD);
    const auto b = reference::pao_sah(reference_input(1, 0.5), VG, VD);
    const auto c = reference::pao_sah(reference_input(2, 1.0), VG, VD);
    if (!a || !b || !c) return std::unexpected(!a ? a.error() : !b ? b.error() : c.error());
    return ReferenceValue{a->current.to_double(), b->current.to_double(), c->current.to_double(),
                          a->quadrature_error * std::abs(a->current.to_double()),
                          a->surface_source_V.to_double(), a->bulk_potential_V.to_double()};
}

}  // namespace

BenchmarkResult run_tvb20(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-20";
    r.title = "Long-channel MOSFET against Pao-Sah: the channel-length limit";
    r.info = info("TVB-20");
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const reference::PaoSahInput in = reference_input(1, 1.0);
    if (!std::isfinite(in.flat_band_V)) {
        r.run_error = "reference: the flat-band voltage";
        return r;
    }
    const double VT = in.thermal_voltage_V;
    // The reference's second method at run time (16.18 (d)), at the judged corners.
    for (const auto [VG, VD] : {std::pair{0.6, 0.05}, std::pair{2.0, 0.5}}) {
        const auto e = reference::pao_sah_conductance_identity(in, VG, VD, 0x1p-17);
        if (!e || *e > 1e-10) {
            r.run_error = "reference: the conductance identity at V_G " + shortest(VG) + " V, V_D " + shortest(VD) +
                          " V: " + (e ? shortest(*e) : e.error());
            return r;
        }
        r.notes.push_back("conductance identity at V_G " + shortest(VG) + ", V_D " + shortest(VD) + " V: " +
                          shortest(*e));
    }
    const Shared& s = shared_run();
    if (s.error) {
        r.run_error = *s.error;
        return r;
    }
    r.notes.push_back("V_FB " + shortest(in.flat_band_V) + " V, mu_n " + shortest(in.mobility_cm2_Vs) + " cm^2/Vs");
    for (std::size_t k = 0; k < std::size(drains); ++k) {
        for (const double VG : judged_gates[k]) {
            const auto ref = reference_at(VG, drains[k]);
            if (!ref) {
                r.run_error = "reference at V_G " + shortest(VG) + ": " + ref.error();
                return r;
            }
            Quantity q;
            q.name = "Q, V_G " + shortest(VG) + " V, V_D " + shortest(drains[k]) + " V";
            for (int l = 0; l < levels; ++l) q.ladder[l] = s.ladder[static_cast<std::size_t>(l)].at({k, VG});
            q.h_tol = s.h_tol.at({k, VG});
            q.h_tolfine = s.h_tolfine.at({k, VG});
            q.h2_tol = s.h2_tol.at({k, VG});
            q.h2_tolfine = s.h2_tolfine.at({k, VG});
            q.ref = ref->Q;
            q.ref_cut = ref->cut;
            q.ref_form = ref->form;
            q.u_ref = ref->u;
            q.applicability_value = ref->psi_s;
            q.applicability_threshold = ref->psi_B + 3.0 * VT;
            q.applicable = ref->psi_s >= q.applicability_threshold;
            judge_limit(r, q);
        }
    }
    const auto cs = reference::charge_sheet(in, 1.2, 0.05);
    if (cs) {
        const auto ps = reference_at(1.2, 0.05);
        r.notes.push_back("Brews' charge sheet at V_G 1.2, V_D 0.05 V (diagnostic, 16.18 (e)): " +
                          shortest(cs->current.to_double()) + " A against Pao-Sah " + (ps ? shortest(ps->Q) : "?") +
                          "; its closed form against Gauss-Kronrod " + shortest(cs->second_method));
    }
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb21(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-21";
    r.title = "Subthreshold swing of the long-channel MOSFET";
    r.info = info("TVB-21");
    r.info.form = "S = 2 h_G / log10(Q(V_G + h_G) / Q(V_G - h_G)) on Pao-Sah, h_G = 20 mV";
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const reference::PaoSahInput in = reference_input(1, 1.0);
    if (!std::isfinite(in.flat_band_V)) {
        r.run_error = "reference: the flat-band voltage";
        return r;
    }
    const Shared& s = shared_run();
    if (s.error) {
        r.run_error = *s.error;
        return r;
    }
    const double VT = in.thermal_voltage_V, VD = drains[0];
    const auto S_of = [&](double hi, double lo) { return 2.0 * h_G / std::log10(hi / lo); };
    for (const double VG : swing_gates) {
        const auto a = reference_at(VG + h_G, VD), b = reference_at(VG - h_G, VD), mid = reference_at(VG, VD);
        if (!a || !b || !mid) {
            r.run_error = "reference at V_G " + shortest(VG) + ": " + (!a ? a.error() : !b ? b.error() : mid.error());
            return r;
        }
        const auto at = [&](const Values& v) { return S_of(v.at({0, VG + h_G}), v.at({0, VG - h_G})); };
        Quantity q;
        q.name = "S, V_G " + shortest(VG) + " V, V_D " + shortest(VD) + " V";
        for (int l = 0; l < levels; ++l) q.ladder[l] = at(s.ladder[static_cast<std::size_t>(l)]);
        q.h_tol = at(s.h_tol);
        q.h_tolfine = at(s.h_tolfine);
        q.h2_tol = at(s.h2_tol);
        q.h2_tolfine = at(s.h2_tolfine);
        q.ref = S_of(a->Q, b->Q);
        q.ref_cut = S_of(a->cut, b->cut);
        q.ref_form = S_of(a->form, b->form);
        // The reference's quadrature estimate carried to S: dS = S^2 / (2 h_G ln 10) (u_a / Q_a + u_b / Q_b).
        q.u_ref = q.ref * q.ref / (2.0 * h_G * std::log(10.0)) * (a->u / a->Q + b->u / b->Q);
        q.applicability_value = mid->psi_s;
        q.applicability_threshold = mid->psi_B + 3.0 * VT;
        q.applicable = mid->psi_s >= q.applicability_threshold;
        judge_limit(r, q);
        const double Cdep = std::sqrt(in.q_C * in.eps_F_per_cm * NA / (2.0 * mid->psi_s));
        r.notes.push_back(q.name + ": the depletion formula ln 10 V_T (1 + C_dep / C_ox) " +
                          shortest(std::log(10.0) * VT * (1.0 + Cdep / in.oxide_capacitance_F_per_cm2)) +
                          " V/dec (diagnostic, DT-Q6); psi_s " + shortest(mid->psi_s) + " V, psi_B " +
                          shortest(mid->psi_B) + " V");
    }
    apply_known_findings(r, known_findings());
    return r;
}

// TVB-32 on this benchmark's device (16.44): L = L0, the V_D = 50 mV sweep from zero bias; the equilibrium.
Tvb32Run tvb20_consistency() {
    Tvb32Run run;
    const device::Device d = mosfet(L0, 1);
    const std::string name = "TVB-20 MOSFET L = " + shortest(L0) + " cm";
    const solve::BiasOptions o = options(tol_newton);
    std::vector<std::vector<double>> points{{0.0, 0.0, 0.0, 0.0}};
    for (const double vg : judged_gates[0]) points.push_back({0.0, drains[0], 0.0, vg});
    const auto s = solve::sweep_bias(d, points, o);
    if (!s || s->stopped || s->points.size() != points.size()) {
        return {{}, name + ": " + (!s ? s.error().message : s->stopped ? s->stopped->message : "short sweep")};
    }
    Tvb32Collector col{run.entries, name, d.contacts().size()};
    const consistency::DeviceData data = device_data(d, o.models);
    for (const auto& p : s->points) col.dc(p, d.mesh(), data, bias_label(p));
    std::string error;
    if (!col.equilibrium(d, o.models, error)) return {std::move(run.entries), error};
    return run;
}

}  // namespace NiTCAD::verify::benchmark
