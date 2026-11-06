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
// On the unstructured families (2D generated, 3D generated, 3D swept; CV-Q7, h = (V / N)^(1/D)) the gate is the
// discrete H1 seminorm, sum over edges of coupling_area (e_i - e_j)^2 / length, relative to the exact psi's: order
// p-hat in [0.8, 1.3] (p_th = 1) and the finest level at most 1e-2. Their L2 and L-infinity orders are reported in
// the notes, not gated (14.7, 5.3). Their meshes are the generators' (mesh/generate.hpp) on [-L, L] x [0, L]^(D-1),
// size bound h_max = L / k and interface resolution s = h_max; the swept cross-section is (x, y), swept along z on
// lines h_max apart.
// Tiers (CV-Q9): fast, 1D N = 32 ... 256; full, + 2D tensor 16^2 ... 128^2 and 3D tensor 8^3 ... 64^3 (the 14.7
// ladder column, owner ruling CV-Q13 in 14.9), 2D generated k = 16 ... 128, 3D generated and swept k = 8 ... 32. The
// 3D levels are solved by MKL PARDISO (ruling CV-Q15: Eigen SparseLU does not fit 64^3); without the MKL runtime loaded
// they are not_evaluated. The finest 2D and 3D tensor L2 misses are known finding (iv) (CV-Q12).
// Reference independence I1: psi_c, n_c, p_c and the electrode's potential at zero bias are NiTCAD's own equilibrium
// contact values (input data); the rest is analytic. Second method: the analytic first and second derivatives of g
// and h are checked against complex-step derivatives at run time (verify/reference/complex_step.hpp).
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <memory>
#include <numbers>
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

// Region 0 is the oxide (x < 0), region 1 the silicon. A generated mesh brings its P4 certificate; its node regions
// are the generator's, and must agree with the side of x = 0 (no node lies on the interface: it is straddled).
std::expected<device::Device, std::string> make_device(mesh::Mesh m,
                                                       std::optional<mesh::InterfaceCertificate> certificate = {},
                                                       const std::vector<std::int32_t>* node_region = nullptr) {
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = m.points()[i][0];
        if (x == 0.0) return std::unexpected("node " + std::to_string(i) + " lies on the interface");
        const bool oxide = x < 0.0;
        region[i] = oxide ? 0 : 1;
        if (node_region && (*node_region)[i] != region[i]) {
            return std::unexpected("node " + std::to_string(i) + "'s generated region is not its side of x = 0");
        }
        if (!oxide) acceptors[i] = acceptors_cm3;
    }
    const auto* xmin = m.find_boundary("x_min");
    const auto* xmax = m.find_boundary("x_max");
    if (!xmin || !xmax) return std::unexpected("the mesh has no x_min or x_max patch");
    device::Contact gate{"gate", device::ContactKind::electrode, xmin->nodes};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto substrate = xmax->nodes;
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate), {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .interface_certificate = std::move(certificate)});
    if (!d) return std::unexpected(d.error().message);
    return std::move(*d);
}

// The generated families (CV-Q7): the domain [-L, L] x [0, L]^(D-1), h the generator's size bound and the interface
// resolution s = h (the straddle's width, as the tensor axis's L / n_side); patches x_min and x_max at x = -L and L.
mesh::Point xyz(double x, double y, double z = 0.0) { return {x, y, z}; }

std::expected<mesh::GeneratedMesh, base::Error> generated_2d(double hmax) {
    mesh::Geometry2D g;
    g.polygons = {{{xyz(-L, 0), xyz(0, 0), xyz(0, L), xyz(-L, L)}, 0}, {{xyz(0, 0), xyz(L, 0), xyz(L, L), xyz(0, L)}, 1}};
    g.patches = {{"x_min", xyz(-L, 0), xyz(-L, L)}, {"x_max", xyz(L, 0), xyz(L, L)}};
    mesh::MeshOptions2D o;
    o.max_size = hmax;
    o.interface_resolution = hmax;
    return mesh::generate_mesh(g, o);
}

std::expected<mesh::GeneratedMesh, base::Error> generated_3d(double hmax) {
    mesh::Geometry3D g;
    g.domain.vertices = {xyz(-L, 0, 0), xyz(L, 0, 0), xyz(L, L, 0), xyz(-L, L, 0),
                         xyz(-L, 0, L), xyz(L, 0, L), xyz(L, L, L), xyz(-L, L, L)};
    g.domain.faces = {{{0, 3, 7, 4}, "x_min", {}}, {{1, 2, 6, 5}, "x_max", {}}, {{0, 1, 5, 4}, "y_min", {}},
                      {{3, 2, 6, 7}, "y_max", {}}, {{0, 1, 2, 3}, "z_min", {}}, {{4, 5, 6, 7}, "z_max", {}}};
    g.interfaces = {{{xyz(0, 0, 0), xyz(0, L, 0), xyz(0, 0, L)}}};
    g.regions = {{xyz(-L / 2, L / 2, L / 2), 0}, {xyz(L / 2, L / 2, L / 2), 1}};
    mesh::MeshOptions3D o;
    o.max_size = hmax;
    o.interface_resolution = hmax;
    return mesh::generate_mesh(g, o);
}

// Swept: the (x, y) cross-section meshed by the 2D generator (the interface in it), swept along z on uniform lines
// h apart.
std::expected<mesh::GeneratedMesh, base::Error> swept(double hmax) {
    mesh::ExtrudedGeometry x;
    x.axis = 2;
    x.stations = {0.0, L};
    x.cells = {{"oxide", {xyz(-L, 0), xyz(0, 0), xyz(0, L), xyz(-L, L)}}, {"silicon", {xyz(0, 0), xyz(L, 0), xyz(L, L), xyz(0, L)}}};
    x.regions = {{0}, {1}};
    x.patches = {{"x_min", xyz(-L, 0, 0), xyz(-L, L, L)}, {"x_max", xyz(L, 0, 0), xyz(L, L, L)}};
    mesh::SweptOptions o;
    o.section.max_size = hmax;
    o.section.interface_resolution = hmax;
    const int k = static_cast<int>(std::lround(L / hmax));
    for (int j = 0; j <= k; ++j) o.lines.push_back(j == k ? L : L * j / k);
    return mesh::generate_swept(x, o);
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
    double h, l2, linf, h1, u_iter;
    std::size_t nodes;
};

std::expected<LevelError, std::string> run_level(const device::Device& d, double h_level, const Contacts& c,
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
    std::vector<double> err(n);
    const auto V = d.mesh().volumes();
    for (std::size_t i = 0; i < n; ++i) {
        const double e = s->fields.potential_V[i] - exact.potential_V[i];
        err[i] = e;
        sum += V[i] * e * e;
        volume += V[i];
        worst = std::max(worst, std::abs(e));
    }
    // The discrete H1 seminorm, the sum over edges of coupling_area (e_i - e_j)^2 / length, relative to the exact
    // potential's (the seminorm does not see psi_c).
    double h1 = 0.0, h1_exact = 0.0;
    for (const mesh::Edge& ed : d.mesh().edges()) {
        const auto a = static_cast<std::size_t>(ed.first), b = static_cast<std::size_t>(ed.second);
        const double de = err[b] - err[a], dp = exact.potential_V[b] - exact.potential_V[a];
        h1 += ed.coupling_area * de * de / ed.length;
        h1_exact += ed.coupling_area * dp * dp / ed.length;
    }
    const auto& it = s->convergence.iterations;
    const double u_iter = it.empty() ? 0.0 : it.back().update * VT / scale;
    return LevelError{h_level, std::sqrt(sum / volume) / scale, worst / scale, std::sqrt(h1 / h1_exact), u_iter, n};
}

enum class Kind { tensor, generated_2d, generated_3d, swept };

struct Family {
    std::string name;
    int dimension;
    Kind kind;
    std::vector<int> levels;  // finest first. Tensor: n_side per x half (the cross-section has 2 n_side + 2 nodes per
                              // axis). Generated and swept: k, the size bound h = L / k.
    bool pardiso = false;     // the 3D levels: MKL PARDISO (CV-Q15 ruling; Eigen SparseLU does not fit 64^3)
};

bool unstructured(const Family& f) { return f.kind != Kind::tensor; }

// One level's device and its h: L / n_side on tensor grids, (V / N)^(1/D) on the others (CV-Q7).
std::expected<std::pair<device::Device, double>, std::string> level_device(const Family& f, int k) {
    if (f.kind == Kind::tensor) {
        const std::vector<double> x = x_axis(k);
        const int nodes = static_cast<int>(x.size());
        std::expected<mesh::Mesh, base::Error> m = f.dimension == 1   ? mesh::make_tensor_grid(x)
                                                   : f.dimension == 2 ? mesh::make_tensor_grid(x, uniform(nodes))
                                                                      : mesh::make_tensor_grid(x, uniform(nodes),
                                                                                               uniform(nodes));
        if (!m) return std::unexpected(m.error().message);
        auto d = make_device(std::move(*m));
        if (!d) return std::unexpected(d.error());
        return std::pair{std::move(*d), L / k};
    }
    const double hmax = L / k;
    auto gm = f.kind == Kind::generated_2d   ? generated_2d(hmax)
              : f.kind == Kind::generated_3d ? generated_3d(hmax)
                                             : swept(hmax);
    if (!gm) return std::unexpected("the generator: " + gm.error().message);
    const double n = static_cast<double>(gm->mesh.node_count());
    const double h = std::pow(2.0 * std::pow(L, f.dimension) / n, 1.0 / f.dimension);
    auto d = make_device(std::move(gm->mesh), std::move(gm->certificate), &gm->node_region);
    if (!d) return std::unexpected(d.error());
    return std::pair{std::move(*d), h};
}

// The comparisons a family records, in order: the gated orders, then the finest level's accuracy (14.7).
std::vector<std::string> quantities(const Family& f) {
    if (unstructured(f)) return {"psi H1 order, " + f.name, "psi H1 finest, " + f.name};
    return {"psi L2 order, " + f.name, "psi Linf order, " + f.name, "psi L2 finest, " + f.name};
}

std::string orders_text(const study::OrderStudy& st) {
    std::string s;
    for (const double p : st.orders) s += " " + shortest(p);
    if (st.asymptotic_evaluated) s += ", order spread " + shortest(st.order_spread);
    return s;
}

void judge_family(BenchmarkResult& r, const Family& f, const Contacts& c, const Materials& mat) {
    if (f.pardiso && !linalg::loaded_mkl_runtime()) {  // not run: its input, the MKL runtime, is missing
        const std::string why = "the MKL runtime for the " + f.name + " levels' PARDISO solves (nitcad verify --mkl)";
        for (const std::string& q : quantities(f)) {
            ComparisonRecord rec{q, "1", 0.0, 0.0, 0.0, 0.0, 0.0, "needs PARDISO (CV-Q15)", {}};
            rec.comparison = missing_input(why);
            r.comparisons.push_back(std::move(rec));
        }
        r.notes.push_back(f.name + ": not evaluated without " + why);
        return;
    }
    std::vector<LevelError> levels;
    for (const int k : f.levels) {
        const std::string label = f.name + (unstructured(f) ? " h_max = L/" : " N_side = ") + std::to_string(k);
        auto dev = level_device(f, k);
        if (!dev) {
            r.run_error = label + ": " + dev.error();
            return;
        }
        const auto e = run_level(dev->first, dev->second, c, mat, f.pardiso);
        if (!e) {
            r.notes.push_back(label + ": the solve failed: " + e.error());
            levels.clear();
            break;
        }
        levels.push_back(*e);
        r.notes.push_back(label + " (" + std::to_string(e->nodes) + " nodes, h " + shortest(e->h) + " cm): L2 " +
                          shortest(e->l2) + ", Linf " + shortest(e->linf) + ", H1 " + shortest(e->h1) + ", U_iter " +
                          shortest(e->u_iter));
    }
    const auto record = [&](const std::string& quantity, double sim, double ref, double metric, const std::string& crit,
                            Comparison cmp) {
        ComparisonRecord rec{quantity, "1", levels.empty() ? 0.0 : levels.front().h, sim, ref, 0.0, metric, crit, {}};
        rec.comparison = std::move(cmp);
        r.comparisons.push_back(std::move(rec));
    };
    // Tensor: L2 and Linf orders at p_th = 2, finest L2 <= 1e-4. Unstructured: the H1 order at p_th = 1, finest H1
    // <= 1e-2; the L2 and Linf orders are reported in the notes only (14.7).
    const bool un = unstructured(f);
    const double p_th = un ? 1.0 : 2.0;
    const std::string band = un ? "p-hat in [0.8, 1.3]" : "p-hat in [1.8, 2.3]";
    const double tolerance = un ? 1e-2 : l2_tolerance;
    const std::string accuracy = un ? "relative H1 <= 1e-2" : "relative L2 <= 1e-4";
    const std::vector<std::string> q = quantities(f);
    if (levels.size() != f.levels.size()) {
        for (std::size_t j = 0; j + 1 < q.size(); ++j) {
            record(q[j], 0.0, p_th, 0.0, band, judge(r.info, {}, false, band));
        }
        record(q.back(), 0.0, 0.0, 0.0, accuracy, judge(r.info, {}, false, accuracy));
        return;
    }
    std::vector<double> hs, u;
    for (const auto& l : levels) {
        hs.push_back(l.h);
        u.push_back(l.u_iter);
    }
    enum class Norm { l2, linf, h1 };
    const auto name_of = [](Norm n) { return n == Norm::l2 ? " L2" : n == Norm::linf ? " Linf" : " H1"; };
    const auto study_of = [&](Norm norm) {
        std::vector<double> e;
        for (const auto& l : levels) e.push_back(norm == Norm::l2 ? l.l2 : norm == Norm::linf ? l.linf : l.h1);
        return study::order_study(hs, e, u);
    };
    const std::vector<Norm> gated = un ? std::vector<Norm>{Norm::h1} : std::vector<Norm>{Norm::l2, Norm::linf};
    for (std::size_t j = 0; j < gated.size(); ++j) {
        const study::OrderStudy st = study_of(gated[j]);
        const study::OrderVerdict v = study::judge_order(st, {.p_th = p_th});
        r.notes.push_back(f.name + name_of(gated[j]) + " orders:" + orders_text(st));
        record(q[j], st.observed, p_th, st.observed, band,
               v.judged ? judge(r.info, {}, v.pass, band) : Comparison{Outcome::not_evaluated, v.reason, {}});
    }
    if (un) {  // reported, not gated (14.7; 5.3: L2 orders are not guaranteed on unstructured meshes)
        for (const Norm norm : {Norm::l2, Norm::linf}) {
            const study::OrderStudy st = study_of(norm);
            r.notes.push_back(f.name + name_of(norm) + " orders (reported, not gated):" + orders_text(st) +
                              "; observed " + shortest(st.observed));
        }
    }
    const double finest = un ? levels.front().h1 : levels.front().l2;
    record(q.back(), finest, 0.0, finest, accuracy, judge(r.info, {}, finest <= tolerance, accuracy));
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
    const auto made = make_device(*mesh::make_tensor_grid(x_axis(15)));
    if (!made) {
        r.run_error = "the equilibrium device: " + made.error();
        return r;
    }
    const device::Device& d0 = *made;
    const auto eq = solve::solve_bias(d0, std::vector<double>{0.0, 0.0}, options());
    if (!eq) {
        r.run_error = "equilibrium: " + eq.error().message;
        return r;
    }
    const std::size_t gate = static_cast<std::size_t>(d0.contacts()[0].nodes.front());
    const std::size_t sub = static_cast<std::size_t>(d0.contacts()[1].nodes.front());
    Contacts c{eq->fields.potential_V[sub], eq->fields.n_cm3[sub], eq->fields.p_cm3[sub], 0.0};
    c.gate_bias_V = c.psi_c - eq->fields.potential_V[gate];  // the electrode's potential then equals psi_c
    std::vector<Family> families{{"1D tensor", 1, Kind::tensor, {127, 63, 31, 15}}};
    if (tier >= Tier::full) {
        families.push_back({"2D tensor", 2, Kind::tensor, {63, 31, 15, 7}});
        families.push_back({"3D tensor", 3, Kind::tensor, {31, 15, 7, 3}, true});  // CV-Q13: in full, with PARDISO
        // CV-Q7: 4 levels for 2D generated, 3 for 3D generated and swept, r = 2 in h_max.
        families.push_back({"2D generated", 2, Kind::generated_2d, {128, 64, 32, 16}});
        families.push_back({"3D generated", 3, Kind::generated_3d, {32, 16, 8}, true});
        families.push_back({"3D swept", 3, Kind::swept, {32, 16, 8}, true});
    }
    for (const Family& f : families) {
        judge_family(r, f, c, mat);
        if (r.run_error) return r;
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
