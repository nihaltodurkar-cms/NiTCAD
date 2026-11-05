// TVB-14, the full drift-diffusion manufactured solution (DECISIONS Unit 55 §6, 14.4; criteria 14.7: "as TVB-04, for
// psi, phi_n and phi_p").
//
// Code under test: the coupled drift-diffusion system (box-method Poisson, Scharfetter-Gummel electron and hole
// fluxes, SRH recombination with Scharfetter lifetimes, Boltzmann statistics, ohmic contacts) and solve_bias, out of
// equilibrium. Construction: n-silicon, N_D = 1e16 cm^-3, on [0, L]^D, L = 1 um; ohmic contacts at x = 0 (0 V) and
// x = L (V_R = 0.1 V), natural (zero-flux) sides. With u = x / L, v = y / L, w = z / L, the exact state is
//     psi - psi_c = V_R u + A b(u) H,   phi_n = V_R u + a_n b(u) H,   phi_p = V_R u + a_p b(u) H,
//     b(u) = 4 u (1 - u),   H = h(v) h(w),   h(v) = 1 + (3 v^2 - 2 v^3) / 2   (h' = 0 on the sides),
// A = 0.1 V, a_n = 0.08 V, a_p = 0.15 V; every field takes its contact values (b = 0) and has zero normal derivative on
// the sides. The densities are n = n_c e^{(psi - psi_c - phi_n) / V_T} and p = p_c e^{(phi_p - psi + psi_c) / V_T},
// which equal the contacts' n_c and p_c there. The manufactured sources (the TVQ3 hook, lumped at the nodes, CV-Q5):
//     S   = -eps Lap(psi) / q + n - p - N_D                        (charge; -div(eps grad psi) = q (p - n + N_D + S))
//     G_n = R + mu_n (grad n . grad phi_n + n Lap(phi_n))         (div J_n = q (R - G_n), J_n = -q mu_n n grad phi_n)
//     G_p = R - mu_p (grad p . grad phi_p + p Lap(phi_p))         (div J_p = -q (R - G_p), J_p = -q mu_p p grad phi_p)
// with R the SRH rate of the exact state, (n p - n_ie^2) / (tau_p (n + n_ie) + tau_n (p + n_ie)), n_ie^2 = n_c p_c,
// tau the material's Scharfetter lifetimes at N_D, and constant mobilities (doping_mobility off, Caughey-Thomas holes:
// the material's mu_max). The sources are evaluated in double-double (verify/reference/double_double.hpp) and rounded
// once.
// Errors at the nodes: psi against the exact psi; phi_n = psi - psi_c - V_T ln(n / n_c) and phi_p = psi - psi_c +
// V_T ln(p / p_c) from NiTCAD's psi, n and p against the exact phi_n and phi_p. Normalisers (14.7): for psi the potential
// scale max|psi - psi_eq| + V_T (psi_eq = psi_c, uniform doping), for phi_n and phi_p the same norm of the exact field
// at the level's nodes; the H1 seminorm of each is relative to the exact field's (CV-Q17).
// Criteria (E, frozen 14.7, "as TVB-04", for each of psi, phi_n, phi_p): on tensor ladders the L2 and L-infinity
// orders p-hat in [1.8, 2.3] and the finest L2 at most 1e-4; on the unstructured families (2D generated, 3D generated,
// 3D swept) the H1 order p-hat in [0.8, 1.3] and the finest H1 at most 1e-2, their L2 and L-infinity orders reported.
// U_iter per level is Newton's last update (the potential's in V_T, the densities' relative, so times V_T in volts)
// over the normaliser, and an order is computed only when U_iter <= e / 10 (CV-Q8).
// Ladders (14.7, "as TVB-04"): fast, 1D N = 32 ... 256 nodes; full, + 2D tensor 16^2 ... 128^2, 3D tensor 8^3 ... 64^3
// (PARDISO), 2D generated h_max = L/16 ... L/128, 3D generated and swept L/8 ... L/32 (PARDISO; CV-Q15). 14.4's
// families column names 1D and 2D tensor for fast; 14.7's frozen ladder column ("as TVB-04") governs, and the
// difference is recorded in DECISIONS 14.8.
// Reference independence I1: psi_c, n_c and p_c are NiTCAD's equilibrium contact values, the mobilities, lifetimes and
// permittivity the material's parameters (input data); the rest is analytic. Second method: the analytic gradients and
// Laplacians of the three fields, and the gradient of n, against complex-step derivatives at run time.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

using reference::DD;

constexpr double L = 1e-4;  // cm
constexpr double donors_cm3 = 1e16;
constexpr double bias_V = 0.1;  // V_R
constexpr double amp_psi = 0.1, amp_n = 0.08, amp_p = 0.15;
constexpr double l2_tolerance = 1e-4;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-14";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PA;
    r.citation = "Method of manufactured solutions (Roache, Verification and Validation in Computational Science and "
                 "Engineering, 1998; Salari and Knupp, SAND2000-1444); box-method orders per DECISIONS Unit 55 5.3 "
                 "(Eymard, Gallouet and Herbin, Finite Volume Methods, 2000 [unverified])";
    r.version = "DECISIONS Unit 55 14.4, 14.7";
    r.equations = {1, 2, 3};
    r.form = "psi, phi_n, phi_p = V_R x/L + a b(x/L) h(y/L) h(z/L); drift-diffusion with manufactured S, G_n, G_p";
    r.assumptions = {"steady state", "Boltzmann statistics", "SRH with Scharfetter lifetimes, mid-gap trap",
                     "constant mobilities", "manufactured sources lumped at the nodes (TVQ3, CV-Q5)"};
    r.derivation = "DECISIONS Unit 55 14.4 and this file's header: the sources from the exact state, in double-double";
    r.verification = "E: observed order and finest-level accuracy, tensor and unstructured families";
    r.second_method = "complex-step derivatives of the fields and of n against the analytic ones, checked at run time";
    return r;
}

// The profiles and their analytic derivatives, generic over double, DD and complex<double>.
template <class T>
T b0(T u) { return 4.0 * u * (1.0 - u); }
template <class T>
T b1(T u) { return 4.0 - 8.0 * u; }
template <class T>
T b2(T) { return T(-8.0); }
template <class T>
T h0(T v) { return 1.0 + 0.5 * (3.0 * v * v - 2.0 * v * v * v); }
template <class T>
T h1(T v) { return 3.0 * v * (1.0 - v); }
template <class T>
T h2(T v) { return 3.0 - 6.0 * v; }

// A field V_R u + a b(u) H(v, w) at reduced coordinates, with its gradient [V/cm] and Laplacian [V/cm^2].
template <class T>
struct Field {
    T value;
    T grad[3];
    T lap;
};

template <class T>
Field<T> field(T u, T v, T w, int D, double a) {
    T H(1.0), Hy(0.0), Hz(0.0), lapH(0.0);
    if (D == 2) {
        H = h0(v);
        Hy = h1(v);
        lapH = h2(v);
    } else if (D == 3) {
        H = h0(v) * h0(w);
        Hy = h1(v) * h0(w);
        Hz = h0(v) * h1(w);
        lapH = h2(v) * h0(w) + h0(v) * h2(w);
    }
    const T b = b0(u);
    Field<T> f;
    f.value = bias_V * u + a * b * H;
    f.grad[0] = (bias_V + a * b1(u) * H) / L;
    f.grad[1] = a * b * Hy / L;
    f.grad[2] = a * b * Hz / L;
    f.lap = a * (b2(u) * H + b * lapH) / (L * L);
    return f;
}

struct Contacts {
    double psi_c, n_c, p_c;
};

struct Physics {
    double VT, eps, mu_n, mu_p, tau_n, tau_p;
};

// The analytic derivatives against complex-step ones; empty when they agree.
std::string check_derivatives(const Physics& ph, const Contacts& c) {
    using reference::complex_d;
    const auto close = [](double x, double y, double scale) {
        return std::abs(x - y) <= 1e-12 * std::max(scale, std::abs(y));
    };
    const auto cs = [](auto&& f, double at) { return reference::complex_step(f, at).derivative; };
    for (const int D : {1, 2, 3}) {
        for (const double a : {amp_psi, amp_n, amp_p}) {
            for (int i = 0; i <= 4; ++i) {
                for (int j = 0; j <= 4; ++j) {
                    const double u = i / 4.0, v = j / 4.0, w = (4 - j) / 4.0;
                    const Field<double> f = field(u, v, w, D, a);
                    const double gx = cs([&](complex_d z) { return field(z, complex_d(v), complex_d(w), D, a).value; }, u);
                    const double gy = cs([&](complex_d z) { return field(complex_d(u), z, complex_d(w), D, a).value; }, v);
                    const double gz = cs([&](complex_d z) { return field(complex_d(u), complex_d(v), z, D, a).value; }, w);
                    const double xx = cs([&](complex_d z) { return field(z, complex_d(v), complex_d(w), D, a).grad[0]; }, u);
                    const double yy = cs([&](complex_d z) { return field(complex_d(u), z, complex_d(w), D, a).grad[1]; }, v);
                    const double zz = cs([&](complex_d z) { return field(complex_d(u), complex_d(v), z, D, a).grad[2]; }, w);
                    const double gs = std::abs(f.grad[0]) + 1.0 / L, ls = std::abs(f.lap) + 1.0 / (L * L);
                    if (!close(gx / L, f.grad[0], gs) || !close(gy / L, f.grad[1], gs) ||
                        !close(gz / L, f.grad[2], gs) || !close((xx + yy + zz) / L, f.lap, ls)) {
                        return "the field derivatives disagree with complex step at D = " + std::to_string(D) +
                               ", u = " + shortest(u) + ", v = " + shortest(v);
                    }
                    // grad n = n (grad psi - grad phi_n) / V_T, along x.
                    const auto n_of = [&](complex_d z) {
                        const auto ps = field(z, complex_d(v), complex_d(w), D, amp_psi).value;
                        const auto fn = field(z, complex_d(v), complex_d(w), D, amp_n).value;
                        return c.n_c * std::exp((ps - fn) / ph.VT);
                    };
                    const Field<double> ps = field(u, v, w, D, amp_psi), fn = field(u, v, w, D, amp_n);
                    const double n = c.n_c * std::exp((ps.value - fn.value) / ph.VT);
                    const double dn = n * (ps.grad[0] - fn.grad[0]) / ph.VT;
                    if (!close(cs(n_of, u) / L, dn, std::abs(dn) + n / L)) {
                        return "grad n disagrees with complex step at D = " + std::to_string(D) + ", u = " + shortest(u);
                    }
                }
            }
        }
    }
    return {};
}

struct NodeSource {
    double psi, n, p, S, Gn, Gp, R;
};

// The exact state and the three sources at one node, in double-double.
NodeSource source_at(const mesh::Point& pt, int D, const Physics& ph, const Contacts& c) {
    const DD u(pt[0] / L), v(D >= 2 ? pt[1] / L : 0.0), w(D >= 3 ? pt[2] / L : 0.0);
    const Field<DD> ps = field(u, v, w, D, amp_psi), fn = field(u, v, w, D, amp_n), fp = field(u, v, w, D, amp_p);
    const DD VT(ph.VT);
    const DD n = DD(c.n_c) * exp((ps.value - fn.value) / VT);
    const DD p = DD(c.p_c) * exp((fp.value - ps.value) / VT);
    DD gn_dot(0.0), gp_dot(0.0);
    for (int k = 0; k < 3; ++k) {
        gn_dot += n * (ps.grad[k] - fn.grad[k]) / VT * fn.grad[k];  // grad n . grad phi_n
        gp_dot += p * (fp.grad[k] - ps.grad[k]) / VT * fp.grad[k];  // grad p . grad phi_p
    }
    const DD nie = sqrt(DD(c.n_c) * DD(c.p_c));
    const DD R = (n * p - nie * nie) / (DD(ph.tau_p) * (n + nie) + DD(ph.tau_n) * (p + nie));
    const DD Gn = R + DD(ph.mu_n) * (gn_dot + n * fn.lap);
    const DD Gp = R - DD(ph.mu_p) * (gp_dot + p * fp.lap);
    const DD S = -DD(ph.eps) * ps.lap / DD(base::q_C) + n - p - DD(donors_cm3);
    return {(DD(c.psi_c) + ps.value).to_double(), n.to_double(), p.to_double(), S.to_double(), Gn.to_double(),
            Gp.to_double(), R.to_double()};
}

// The potential scale max|psi - psi_eq| + V_T over the domain (b >= 0 and H >= 1, so the maximum is at the largest H,
// v = w = 1, over u), by dense sampling of the polynomial.
double potential_scale(int D, double VT) {
    const double Hmax = std::pow(h0(1.0), D - 1);
    double m = 0.0;
    for (int k = 0; k <= 100000; ++k) {
        const double u = k / 100000.0;
        m = std::max(m, bias_V * u + amp_psi * b0(u) * Hmax);
    }
    return m + VT;
}

mesh::Point xyz(double x, double y, double z = 0.0) { return {x, y, z}; }

std::vector<double> uniform(int nodes) {
    std::vector<double> v;
    for (int k = 0; k < nodes; ++k) v.push_back(k == nodes - 1 ? L : L * k / (nodes - 1));
    return v;
}

std::expected<device::Device, std::string> make_device(mesh::Mesh m,
                                                       std::optional<mesh::InterfaceCertificate> certificate = {},
                                                       const std::vector<std::int32_t>* node_region = nullptr) {
    const std::size_t n = m.node_count();
    if (node_region) {
        for (std::size_t i = 0; i < n; ++i) {
            if ((*node_region)[i] != 0) return std::unexpected("node " + std::to_string(i) + " is not in region 0");
        }
    }
    const auto* xmin = m.find_boundary("x_min");
    const auto* xmax = m.find_boundary("x_max");
    if (!xmin || !xmax) return std::unexpected("the mesh has no x_min or x_max patch");
    auto left = xmin->nodes;
    auto right = xmax->nodes;
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::vector<double>(n, donors_cm3),
                                     .acceptors = std::vector<double>(n, 0.0),
                                     .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                  {"right", device::ContactKind::ohmic, std::move(right)}},
                                     .interface_certificate = std::move(certificate)});
    if (!d) return std::unexpected(d.error().message);
    return std::move(*d);
}

std::expected<mesh::GeneratedMesh, base::Error> generated_2d(double hmax) {
    mesh::Geometry2D g;
    g.polygons = {{{xyz(0, 0), xyz(L, 0), xyz(L, L), xyz(0, L)}, 0}};
    g.patches = {{"x_min", xyz(0, 0), xyz(0, L)}, {"x_max", xyz(L, 0), xyz(L, L)}};
    mesh::MeshOptions2D o;
    o.max_size = hmax;
    o.interface_resolution = hmax;
    return mesh::generate_mesh(g, o);
}

std::expected<mesh::GeneratedMesh, base::Error> generated_3d(double hmax) {
    mesh::Geometry3D g;
    g.domain.vertices = {xyz(0, 0, 0), xyz(L, 0, 0), xyz(L, L, 0), xyz(0, L, 0),
                         xyz(0, 0, L), xyz(L, 0, L), xyz(L, L, L), xyz(0, L, L)};
    g.domain.faces = {{{0, 3, 7, 4}, "x_min", {}}, {{1, 2, 6, 5}, "x_max", {}}, {{0, 1, 5, 4}, "y_min", {}},
                      {{3, 2, 6, 7}, "y_max", {}}, {{0, 1, 2, 3}, "z_min", {}}, {{4, 5, 6, 7}, "z_max", {}}};
    mesh::MeshOptions3D o;
    o.max_size = hmax;
    o.interface_resolution = hmax;
    return mesh::generate_mesh(g, o);
}

std::expected<mesh::GeneratedMesh, base::Error> swept(double hmax) {
    mesh::ExtrudedGeometry x;
    x.axis = 2;
    x.stations = {0.0, L};
    x.cells = {{"silicon", {xyz(0, 0), xyz(L, 0), xyz(L, L), xyz(0, L)}}};
    x.regions = {{0}};
    x.patches = {{"x_min", xyz(0, 0, 0), xyz(0, L, L)}, {"x_max", xyz(L, 0, 0), xyz(L, L, L)}};
    mesh::SweptOptions o;
    o.section.max_size = hmax;
    o.section.interface_resolution = hmax;
    const int k = static_cast<int>(std::lround(L / hmax));
    for (int j = 0; j <= k; ++j) o.lines.push_back(j == k ? L : L * j / k);
    return mesh::generate_swept(x, o);
}

solve::BiasOptions options() {
    solve::BiasOptions o;
    o.models.srh = true;
    o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.doping_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.newton.tol_update = 1e-12;
    return o;
}

constexpr int field_count = 3;  // psi, phi_n, phi_p
constexpr const char* field_names[field_count] = {"psi", "phi_n", "phi_p"};

struct LevelError {
    double h;
    double l2[field_count], linf[field_count], h1[field_count];
    double u_iter[field_count];
    double srh_share;  // max |R| / max |G_p|: the weight of SRH in the hole balance (reported)
    std::size_t nodes;
};

std::expected<LevelError, std::string> run_level(const device::Device& d, double h_level, const Contacts& c,
                                                 const Physics& ph, bool pardiso) {
    const std::size_t n = d.mesh().node_count();
    const int D = d.mesh().dimension();
    results::NodeFields exact;
    std::vector<double> phi_ex[field_count];  // psi - psi_c, phi_n, phi_p
    auto sources = std::make_shared<assemble::ManufacturedSources>();
    double R_max = 0.0, Gp_max = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& pt = d.mesh().points()[i];
        const NodeSource s = source_at(pt, D, ph, c);
        exact.potential_V.push_back(s.psi);
        exact.n_cm3.push_back(s.n);
        exact.p_cm3.push_back(s.p);
        sources->charge_cm3.push_back(s.S);
        sources->electron_generation.push_back(s.Gn);
        sources->hole_generation.push_back(s.Gp);
        R_max = std::max(R_max, std::abs(s.R));
        Gp_max = std::max(Gp_max, std::abs(s.Gp));
        const double u = pt[0] / L, v = D >= 2 ? pt[1] / L : 0.0, w = D >= 3 ? pt[2] / L : 0.0;
        phi_ex[0].push_back(s.psi - c.psi_c);
        phi_ex[1].push_back(field(u, v, w, D, amp_n).value);
        phi_ex[2].push_back(field(u, v, w, D, amp_p).value);
    }
    solve::BiasOptions o = options();
    o.models.manufactured = sources.get();  // kept alive by sources through the solve
    if (pardiso) o.linear = {.backend = linalg::SolverBackend::mkl_pardiso, .threads = pardiso_threads()};
    const auto s = solve::solve_bias(d, std::vector<double>{0.0, bias_V}, o, &exact);
    if (!s) return std::unexpected(s.error().message);
    const auto& f = s->fields;
    std::vector<double> sim[field_count];
    for (std::size_t i = 0; i < n; ++i) {
        const double dpsi = f.potential_V[i] - c.psi_c;
        sim[0].push_back(dpsi);
        sim[1].push_back(dpsi - ph.VT * std::log(f.n_cm3[i] / c.n_c));
        sim[2].push_back(dpsi + ph.VT * std::log(f.p_cm3[i] / c.p_c));
    }
    const auto V = d.mesh().volumes();
    const auto& it = s->convergence.iterations;
    const double update_V = it.empty() ? 0.0 : it.back().update * ph.VT;
    LevelError le{};
    le.h = h_level;
    le.nodes = n;
    le.srh_share = Gp_max > 0.0 ? R_max / Gp_max : 0.0;
    const double psi_scale = potential_scale(D, ph.VT);
    for (int q = 0; q < field_count; ++q) {
        double sum = 0.0, sum_ex = 0.0, volume = 0.0, worst = 0.0, worst_ex = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double e = sim[q][i] - phi_ex[q][i];
            sum += V[i] * e * e;
            sum_ex += V[i] * phi_ex[q][i] * phi_ex[q][i];
            volume += V[i];
            worst = std::max(worst, std::abs(e));
            worst_ex = std::max(worst_ex, std::abs(phi_ex[q][i]));
        }
        double h1_sum = 0.0, h1_ex = 0.0;
        for (const mesh::Edge& ed : d.mesh().edges()) {
            const auto a = static_cast<std::size_t>(ed.first), b = static_cast<std::size_t>(ed.second);
            const double de = (sim[q][b] - phi_ex[q][b]) - (sim[q][a] - phi_ex[q][a]);
            const double dx = phi_ex[q][b] - phi_ex[q][a];
            h1_sum += ed.coupling_area * de * de / ed.length;
            h1_ex += ed.coupling_area * dx * dx / ed.length;
        }
        const double l2_norm = q == 0 ? psi_scale : std::sqrt(sum_ex / volume);
        const double linf_norm = q == 0 ? psi_scale : worst_ex;
        le.l2[q] = std::sqrt(sum / volume) / l2_norm;
        le.linf[q] = worst / linf_norm;
        le.h1[q] = std::sqrt(h1_sum / h1_ex);
        le.u_iter[q] = update_V / std::min(l2_norm, linf_norm);
    }
    return le;
}

enum class Kind { tensor, generated_2d, generated_3d, swept };

struct Family {
    std::string name;
    int dimension;
    Kind kind;
    std::vector<int> levels;  // finest first. Tensor: nodes per axis; generated and swept: k, h_max = L / k.
    bool pardiso = false;
};

bool unstructured(const Family& f) { return f.kind != Kind::tensor; }

std::expected<std::pair<device::Device, double>, std::string> level_device(const Family& f, int k) {
    if (f.kind == Kind::tensor) {
        const std::vector<double> x = uniform(k);
        std::expected<mesh::Mesh, base::Error> m = f.dimension == 1   ? mesh::make_tensor_grid(x)
                                                   : f.dimension == 2 ? mesh::make_tensor_grid(x, x)
                                                                      : mesh::make_tensor_grid(x, x, x);
        if (!m) return std::unexpected(m.error().message);
        auto d = make_device(std::move(*m));
        if (!d) return std::unexpected(d.error());
        return std::pair{std::move(*d), L / (k - 1)};
    }
    const double hmax = L / k;
    auto gm = f.kind == Kind::generated_2d   ? generated_2d(hmax)
              : f.kind == Kind::generated_3d ? generated_3d(hmax)
                                             : swept(hmax);
    if (!gm) return std::unexpected("the generator: " + gm.error().message);
    const double n = static_cast<double>(gm->mesh.node_count());
    const double h = std::pow(std::pow(L, f.dimension) / n, 1.0 / f.dimension);
    auto d = make_device(std::move(gm->mesh), std::move(gm->certificate), &gm->node_region);
    if (!d) return std::unexpected(d.error());
    return std::pair{std::move(*d), h};
}

// The comparisons a family records, in order: per field the gated orders, then per field the finest accuracy.
std::vector<std::string> quantities(const Family& f) {
    std::vector<std::string> q;
    for (const char* name : field_names) {
        const std::string fn(name);
        if (unstructured(f)) {
            q.push_back(fn + " H1 order, " + f.name);
        } else {
            q.push_back(fn + " L2 order, " + f.name);
            q.push_back(fn + " Linf order, " + f.name);
        }
    }
    for (const char* name : field_names) {
        q.push_back(std::string(name) + (unstructured(f) ? " H1 finest, " : " L2 finest, ") + f.name);
    }
    return q;
}

std::string orders_text(const study::OrderStudy& st) {
    std::string s;
    for (const double p : st.orders) s += " " + shortest(p);
    if (st.asymptotic_evaluated) s += ", order spread " + shortest(st.order_spread);
    return s;
}

void judge_family(BenchmarkResult& r, const Family& f, const Contacts& c, const Physics& ph) {
    const std::vector<std::string> q = quantities(f);
    if (f.pardiso && !linalg::loaded_mkl_runtime()) {
        const std::string why = "the MKL runtime for the " + f.name + " levels' PARDISO solves (nitcad verify --mkl)";
        for (const std::string& name : q) {
            ComparisonRecord rec{name, "1", 0.0, 0.0, 0.0, 0.0, 0.0, "needs PARDISO (CV-Q15)", {}};
            rec.comparison = missing_input(why);
            r.comparisons.push_back(std::move(rec));
        }
        r.notes.push_back(f.name + ": not evaluated without " + why);
        return;
    }
    std::vector<LevelError> levels;
    for (const int k : f.levels) {
        const std::string label = f.name + (unstructured(f) ? " h_max = L/" : " N = ") + std::to_string(k);
        auto dev = level_device(f, k);
        if (!dev) {
            r.run_error = label + ": " + dev.error();
            return;
        }
        const auto e = run_level(dev->first, dev->second, c, ph, f.pardiso);
        if (!e) {
            r.notes.push_back(label + ": the solve failed: " + e.error());
            levels.clear();
            break;
        }
        levels.push_back(*e);
        std::string line = label + " (" + std::to_string(e->nodes) + " nodes, h " + shortest(e->h) + " cm):";
        for (int j = 0; j < field_count; ++j) {
            line += std::string(" ") + field_names[j] + " L2 " + shortest(e->l2[j]) + " Linf " + shortest(e->linf[j]) +
                    " H1 " + shortest(e->h1[j]) + ";";
        }
        line += " U_iter " + shortest(e->u_iter[0]) + ", max|R| / max|G_p| " + shortest(e->srh_share);
        r.notes.push_back(line);
    }
    const auto record = [&](const std::string& quantity, double sim, double ref, double metric, const std::string& crit,
                            Comparison cmp) {
        ComparisonRecord rec{quantity, "1", levels.empty() ? 0.0 : levels.front().h, sim, ref, 0.0, metric, crit, {}};
        rec.comparison = std::move(cmp);
        r.comparisons.push_back(std::move(rec));
    };
    const bool un = unstructured(f);
    const double p_th = un ? 1.0 : 2.0;
    const std::string band = un ? "p-hat in [0.8, 1.3]" : "p-hat in [1.8, 2.3]";
    const double tolerance = un ? 1e-2 : l2_tolerance;
    const std::string accuracy = un ? "relative H1 <= 1e-2" : "relative L2 <= 1e-4";
    const std::size_t orders_per_field = un ? 1 : 2;
    if (levels.size() != f.levels.size()) {
        for (std::size_t j = 0; j < q.size(); ++j) {
            const bool order = j < field_count * orders_per_field;
            record(q[j], 0.0, order ? p_th : 0.0, 0.0, order ? band : accuracy,
                   judge(r.info, {}, false, order ? band : accuracy));
        }
        return;
    }
    std::vector<double> hs;
    for (const auto& l : levels) hs.push_back(l.h);
    enum class Norm { l2, linf, h1 };
    const auto name_of = [](Norm n) { return n == Norm::l2 ? " L2" : n == Norm::linf ? " Linf" : " H1"; };
    const auto study_of = [&](int field, Norm norm) {
        std::vector<double> e, u;
        for (const auto& l : levels) {
            e.push_back(norm == Norm::l2 ? l.l2[field] : norm == Norm::linf ? l.linf[field] : l.h1[field]);
            u.push_back(l.u_iter[field]);
        }
        return study::order_study(hs, e, u);
    };
    const std::vector<Norm> gated = un ? std::vector<Norm>{Norm::h1} : std::vector<Norm>{Norm::l2, Norm::linf};
    std::size_t slot = 0;
    for (int field = 0; field < field_count; ++field) {
        for (const Norm norm : gated) {
            const study::OrderStudy st = study_of(field, norm);
            const study::OrderVerdict v = study::judge_order(st, {.p_th = p_th});
            r.notes.push_back(f.name + " " + field_names[field] + name_of(norm) + " orders:" + orders_text(st));
            record(q[slot++], st.observed, p_th, st.observed, band,
                   v.judged ? judge(r.info, {}, v.pass, band) : Comparison{Outcome::not_evaluated, v.reason, {}});
        }
        if (un) {  // reported, not gated (14.7, 5.3)
            for (const Norm norm : {Norm::l2, Norm::linf}) {
                const study::OrderStudy st = study_of(field, norm);
                r.notes.push_back(f.name + " " + field_names[field] + name_of(norm) +
                                  " orders (reported, not gated):" + orders_text(st) + "; observed " +
                                  shortest(st.observed));
            }
        }
    }
    for (int field = 0; field < field_count; ++field) {
        const double finest = un ? levels.front().h1[field] : levels.front().l2[field];
        record(q[slot++], finest, 0.0, finest, accuracy, judge(r.info, {}, finest <= tolerance, accuracy));
    }
}

}  // namespace

BenchmarkResult run_tvb14(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-14";
    r.title = "Full drift-diffusion manufactured solution";
    r.info = info();
    const physics::Semiconductor si = physics::silicon();
    const Physics ph{base::thermal_voltage(300.0),
                     si.parameters().eps_r * base::eps0_F_per_cm,
                     si.parameters().electron_mobility.mu_max,
                     si.parameters().hole_mobility.mu_max,
                     physics::scharfetter_lifetime(si, physics::Carrier::electron, donors_cm3),
                     physics::scharfetter_lifetime(si, physics::Carrier::hole, donors_cm3)};
    // The contact values at zero bias, from a coarse 1D device (they do not depend on the mesh).
    const auto made = make_device(*mesh::make_tensor_grid(uniform(16)));
    if (!made) {
        r.run_error = "the equilibrium device: " + made.error();
        return r;
    }
    const auto eq = solve::solve_bias(*made, std::vector<double>{0.0, 0.0}, options());
    if (!eq) {
        r.run_error = "equilibrium: " + eq.error().message;
        return r;
    }
    const std::size_t left = static_cast<std::size_t>(made->contacts()[0].nodes.front());
    const Contacts c{eq->fields.potential_V[left], eq->fields.n_cm3[left], eq->fields.p_cm3[left]};
    if (const std::string bad = check_derivatives(ph, c); !bad.empty()) {
        r.run_error = "reference self-check: " + bad;
        return r;
    }
    std::vector<Family> families{{"1D tensor", 1, Kind::tensor, {256, 128, 64, 32}}};
    if (tier >= Tier::full) {
        families.push_back({"2D tensor", 2, Kind::tensor, {128, 64, 32, 16}});
        families.push_back({"3D tensor", 3, Kind::tensor, {64, 32, 16, 8}, true});
        families.push_back({"2D generated", 2, Kind::generated_2d, {128, 64, 32, 16}});
        families.push_back({"3D generated", 3, Kind::generated_3d, {32, 16, 8}, true});
        families.push_back({"3D swept", 3, Kind::swept, {32, 16, 8}, true});
    }
    for (const Family& f : families) {
        judge_family(r, f, c, ph);
        if (r.run_error) return r;
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
