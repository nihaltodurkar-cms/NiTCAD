// Unit 62 slice 4 (DECISIONS.md QV6-QV8): density gradient against the in-house 1D Schrodinger-Poisson
// reference (verify/reference/schrodinger_poisson.hpp, QP14) on a 1D MOS-C at equilibrium (the quasi-static
// sweep: E_F = 0 everywhere, the reference's setting), its mesh convergence, and its limits. Release only, hidden
// ([.dg_validation]); run by name. QV6 is reported, not gated, until QP0 fixes its bounds (the placeholders 10%,
// 0.3 nm and 25% are printed beside the measured values, not used for a verdict).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/quantum_confinement.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/reference/schrodinger_poisson.hpp"

using namespace NiTCAD;
namespace ref = NiTCAD::verify::reference;

namespace {

constexpr double T = 300.0;
constexpr double depth = 3e-5;  // [cm] silicon depth, past the depletion edge at 1e17

enum class Model { classical, mlda, dg };

struct Mos {
    double t_ox, N_A;
    double h = 1e-8;               // surface spacing [cm] (0.1 nm)
    double oxide_mass = 0.0;       // with > 0: an insulator of SiO2's permittivity and band data with this mass
    double gamma = 1.0;
};

// x from -t_ox: 20 oxide cells, the interface between nodes at -+h/2, silicon at spacing h to 3 nm, then growing
// by 1.08 up to 2 nm, to `depth`.
std::vector<double> axis(const Mos& m) {
    std::vector<double> x;
    for (int k = 0; k <= 20; ++k) x.push_back(-m.t_ox + k * (m.t_ox / 20.0));
    double v = m.h, step = m.h;
    while (v < depth) {
        x.push_back(v);
        if (v >= 3e-7) step = std::min(step * 1.08, 2e-7);
        v += step;
    }
    x.push_back(depth);
    return *mesh::straddle_interface(x, 0.0, m.h);
}

device::Device device_of(const Mos& m) {
    mesh::Mesh mesh = *mesh::make_tensor_grid(axis(m));
    const std::size_t n = mesh.node_count();
    std::vector<device::RegionId> region(n, 1);
    std::vector<double> acceptors(n, m.N_A);
    for (std::size_t i = 0; i < n; ++i) {
        if (mesh.points()[i][0] < 0.0) {
            region[i] = 0;
            acceptors[i] = 0.0;
        }
    }
    auto substrate = mesh.find_boundary("x_max")->nodes;
    device::Contact gate{"gate", device::ContactKind::electrode, mesh.find_boundary("x_min")->nodes};
    gate.electrode = {.kind = device::GateElectrode::n_poly};
    physics::InsulatorParameters ox = physics::silicon_dioxide_tunnelling_parameters;
    if (m.oxide_mass > 0.0) ox.tunnel_mass_n = ox.tunnel_mass_p = m.oxide_mass;
    std::vector<device::Region> regions;
    if (m.oxide_mass > 0.0) {
        regions.push_back({"oxide", *physics::Insulator::create(ox)});
    } else {
        regions.push_back({"oxide", physics::silicon_dioxide()});
    }
    regions.push_back({"silicon", physics::silicon()});
    return *device::Device::create({.mesh = std::move(mesh),
                                    .temperature_K = T,
                                    .regions = std::move(regions),
                                    .node_region = std::move(region),
                                    .donors = std::vector<double>(n, 0.0),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {std::move(gate),
                                                 {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
}

solve::BiasOptions options(Model model, double gamma) {
    solve::BiasOptions o;
    o.equations = solve::Equations::equilibrium_poisson;
    o.models.quantum_confinement = model == Model::dg     ? assemble::QuantumConfinement::density_gradient
                                   : model == Model::mlda ? assemble::QuantumConfinement::mlda
                                                          : assemble::QuantumConfinement::off;
    o.models.bgn = false;  // the reference has no bandgap narrowing
    o.models.dg_gamma_n = gamma;
    o.models.dg_gamma_p = gamma;
    return o;
}

struct Layer {
    double N_inv = 0.0;     // [cm^-2]
    double centroid = 0.0;  // [cm]
    double surface_n = 0.0; // the first silicon node's n [cm^-3]
    double peak_at = 0.0;   // where n peaks [cm]
    std::vector<double> x, n;
};

// The electrons over the silicon, box widths from the node coordinates (the box of the first silicon node starts
// at the interface, x = 0, the midpoint of the straddling edge).
Layer layer(const device::Device& d, const results::BiasPoint& p) {
    Layer L;
    const auto& pts = d.mesh().points();
    const std::size_t n = d.mesh().node_count();
    double peak = -1.0;
    bool first = true;
    for (std::size_t i = 0; i < n; ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        const double lo = i > 0 ? 0.5 * (pts[i - 1][0] + pts[i][0]) : pts[i][0];
        const double hi = i + 1 < n ? 0.5 * (pts[i][0] + pts[i + 1][0]) : pts[i][0];
        const double w = hi - std::max(lo, 0.0);
        const double c = p.fields.n_cm3[i];
        L.N_inv += c * w;
        L.centroid += pts[i][0] * c * w;
        if (first) L.surface_n = c;
        first = false;
        if (c > peak) {
            peak = c;
            L.peak_at = pts[i][0];
        }
        L.x.push_back(pts[i][0]);
        L.n.push_back(c);
    }
    L.centroid /= L.N_inv;
    return L;
}

std::optional<results::BiasPoint> solve_at(const device::Device& d, Model model, double vg, double gamma = 1.0) {
    const auto p = solve::solve_bias(d, std::vector<double>{vg, 0.0}, options(model, gamma));
    if (!p) {
        UNSCOPED_INFO("V_G = " << vg << ": " << p.error().message);
        return std::nullopt;
    }
    return *p;
}

// V_FB: the classical gate charge's zero, by bisection in [-2, 1] V.
double flat_band(const device::Device& d) {
    double a = -2.0, b = 1.0;
    const auto q = [&](double v) { return solve_at(d, Model::classical, v)->gate_charge[0]; };
    double qa = q(a);
    REQUIRE((qa <= 0.0) != (q(b) <= 0.0));
    for (int k = 0; k < 50; ++k) {
        const double c = 0.5 * (a + b), qc = q(c);
        if ((qc <= 0.0) == (qa <= 0.0)) {
            a = c;
            qa = qc;
        } else {
            b = c;
        }
    }
    return 0.5 * (a + b);
}

// The gate bias (relative to `offset`) at which f(V) = target, f increasing: bisection in [lo, hi].
template <class F>
double bisect(F f, double target, double lo, double hi) {
    for (int k = 0; k < 40; ++k) {
        const double c = 0.5 * (lo + hi);
        if (f(c) < target) {
            lo = c;
        } else {
            hi = c;
        }
    }
    return 0.5 * (lo + hi);
}

ref::MosSpInput sp_input(const Mos& m, ref::SpElectrons electrons) {
    const physics::Semiconductor si = physics::silicon();
    ref::MosSpInput in;
    in.q_C = base::q_C;
    in.kT_eV = base::thermal_voltage(T);
    in.eps0_F_per_cm = base::eps0_F_per_cm;
    in.a_eV_cm2 = base::hbar_J_s * base::hbar_J_s / (2.0 * base::m0_kg) / base::q_C * 1e4;
    in.eps_si_r = si.parameters().eps_r;
    in.eps_ox_r = physics::silicon_dioxide().parameters().eps_r;
    in.band_gap_eV = physics::band_gap_eV(si, T);
    in.conduction_dos_cm3 = physics::conduction_band_dos(si, T);
    in.valence_dos_cm3 = physics::valence_band_dos(si, T);
    const auto& ox = physics::silicon_dioxide_tunnelling_parameters;
    in.oxide_barrier_eV = si.parameters().electron_affinity_eV - ox.electron_affinity_eV;
    in.oxide_mass = ox.tunnel_mass_n;
    in.valleys = {{.mass_z = 0.916, .mass_dos = 0.19, .degeneracy = 2.0},
                  {.mass_z = 0.19, .mass_dos = std::sqrt(0.916 * 0.19), .degeneracy = 4.0}};
    in.acceptors_cm3 = m.N_A;
    in.oxide_thickness_cm = m.t_ox;
    in.oxide_cells = 100;
    in.h0_cm = 2e-9;
    in.ratio = 1.03;
    in.h_max_cm = 1e-7;
    in.depth_cm = depth;
    in.quantum_depth_cm = 40e-7;
    in.electrons = electrons;
    in.statistics = ref::SpStatistics::boltzmann;
    in.subbands_per_valley = 8;
    return in;
}

ref::MosSpResult sp(const ref::MosSpInput& in, double overdrive) {
    auto r = ref::mos_schrodinger_poisson(in, overdrive);
    if (!r) FAIL(r.error());
    return *r;
}

}  // namespace

TEST_CASE("density gradient against Schrodinger-Poisson (QV6, reported until QP0)", "[.dg_validation]") {
    const double target = 1e12;  // N_inv defining the threshold [cm^-2]
    std::printf("QV6 (placeholders, no verdict until QP0: |dN/N| 10%%, |d centroid| 0.3 nm, |d dVth / dVth| 25%%)\n");
    for (const double t_ox : {1.5e-7, 5e-7}) {
        for (const double N_A : {1e17, 1e18}) {
            const Mos m{t_ox, N_A};
            const device::Device d = device_of(m);
            const double vfb = flat_band(d);
            const auto spq = sp_input(m, ref::SpElectrons::quantum);
            const auto spc = sp_input(m, ref::SpElectrons::classical);
            // Thresholds (V_G - V_FB) at N_inv = 1e12.
            const auto n_of = [&](Model model) {
                return [&, model](double ov) {
                    const auto p = solve_at(d, model, vfb + ov);
                    REQUIRE(p.has_value());
                    return layer(d, *p).N_inv;
                };
            };
            const double vt_cl = bisect(n_of(Model::classical), target, 0.0, 3.0);
            const double vt_dg = bisect(n_of(Model::dg), target, 0.0, 3.0);
            const double vt_ml = bisect(n_of(Model::mlda), target, 0.0, 3.0);
            const double vt_spc = bisect([&](double ov) { return sp(spc, ov).inversion_cm2; }, target, 0.0, 3.0);
            const double vt_spq = bisect([&](double ov) { return sp(spq, ov).inversion_cm2; }, target, 0.0, 3.0);
            const double shift_sp = vt_spq - vt_spc;
            std::printf("t_ox %.1f nm, N_A %.0e: V_FB %.4f V; V_th - V_FB classical %.4f (SP classical %.4f); "
                        "shift DG %.4f MLDA %.4f SP %.4f V (DG %+.1f%%, MLDA %+.1f%%)\n",
                        t_ox * 1e7, N_A, vfb, vt_cl, vt_spc, vt_dg - vt_cl, vt_ml - vt_cl, shift_sp,
                        100.0 * ((vt_dg - vt_cl) / shift_sp - 1.0), 100.0 * ((vt_ml - vt_cl) / shift_sp - 1.0));
            CHECK(std::abs(vt_cl - vt_spc) < 5e-3);  // the classical cross-check (same equations, two codes)
            for (const double over : {0.5, 1.0, 1.5}) {
                const double ov = vt_spc + over;
                const auto q = sp(spq, ov);
                const auto c = sp(spc, ov);
                const auto pd = solve_at(d, Model::dg, vfb + ov);
                const auto pm = solve_at(d, Model::mlda, vfb + ov);
                const auto pc = solve_at(d, Model::classical, vfb + ov);
                REQUIRE(pd.has_value());
                REQUIRE(pm.has_value());
                REQUIRE(pc.has_value());
                const Layer ld = layer(d, *pd), lm = layer(d, *pm), lc = layer(d, *pc);
                std::printf("  V_G - V_FB %.3f: N_inv SP %.4e DG %+.2f%% MLDA %+.2f%% (classical SP %.4e, NiTCAD "
                            "%+.3f%%); centroid SP %.3f nm DG %+.3f MLDA %+.3f nm (classical %.3f)\n",
                            ov, q.inversion_cm2, 100.0 * (ld.N_inv / q.inversion_cm2 - 1.0),
                            100.0 * (lm.N_inv / q.inversion_cm2 - 1.0), c.inversion_cm2,
                            100.0 * (lc.N_inv / c.inversion_cm2 - 1.0), q.centroid_cm * 1e7,
                            (ld.centroid - q.centroid_cm) * 1e7, (lm.centroid - q.centroid_cm) * 1e7,
                            lc.centroid * 1e7);
                CHECK(std::abs(lc.N_inv / c.inversion_cm2 - 1.0) < 1e-2);  // classical cross-check
            }
        }
    }
}

TEST_CASE("density gradient: mesh convergence of the inversion charge (QV7)", "[.dg_validation]") {
    // QV6's 5 nm device (1e17), strong inversion (V_G - V_FB = 2 V): surface spacings 0.2, 0.1, 0.05 nm.
    std::vector<double> q;
    double vfb = 0.0;
    for (const double h : {2e-8, 1e-8, 5e-9}) {
        Mos m{5e-7, 1e17};
        m.h = h;
        const device::Device d = device_of(m);
        if (q.empty()) vfb = flat_band(d);
        const auto p = solve_at(d, Model::dg, vfb + 2.0);
        REQUIRE(p.has_value());
        q.push_back(layer(d, *p).N_inv);
        std::printf("QV7: h %.2f nm, N_inv %.6e cm^-2\n", h * 1e7, q.back());
    }
    const double d1 = std::abs(q[1] / q[0] - 1.0), d2 = std::abs(q[2] / q[1] - 1.0);
    std::printf("QV7: |0.1 / 0.2 - 1| %.3e, |0.05 / 0.1 - 1| %.3e\n", d1, d2);
    CHECK(d2 <= 1e-2);
    CHECK(d2 < d1);
}

TEST_CASE("density gradient: the limits x_p -> 0 and gamma -> 0 (QV8)", "[.dg_validation]") {
    // x_p -> 0: the oxide mass m_ox x 4 halves x_p each time; the interface density falls monotonically and the
    // peak moves inward.
    const double ov = 2.0;
    double vfb = 0.0;
    {
        Mos m{5e-7, 1e17};
        vfb = flat_band(device_of(m));
    }
    double last_n = 0.0, last_peak = -1.0;
    bool first = true;
    for (const double mass : {0.5, 2.0, 8.0, 32.0, 128.0}) {
        Mos m{5e-7, 1e17};
        m.oxide_mass = mass;
        const device::Device d = device_of(m);
        const auto p = solve_at(d, Model::dg, vfb + ov);
        REQUIRE(p.has_value());
        const Layer L = layer(d, *p);
        const double xp = physics::barrier_decay_length_cm(mass, 3.0);
        std::printf("QV8: m_ox %.1f, x_p %.4f nm: n(surface) %.4e cm^-3, peak at %.3f nm, N_inv %.5e\n", mass, xp * 1e7,
                    L.surface_n, L.peak_at * 1e7, L.N_inv);
        if (!first) {
            CHECK(L.surface_n < last_n);
            CHECK(L.peak_at >= last_peak);
        }
        first = false;
        last_n = L.surface_n;
        last_peak = L.peak_at;
    }
    // gamma -> 0: the difference from classical, in N_inv and in max |n - n_cl| / max n_cl, proportional to gamma.
    const Mos m{5e-7, 1e17};
    const device::Device d = device_of(m);
    const auto pc = solve_at(d, Model::classical, vfb + ov);
    REQUIRE(pc.has_value());
    const Layer lc = layer(d, *pc);
    double peak_cl = 0.0;
    for (const double c : lc.n) peak_cl = std::max(peak_cl, c);
    std::vector<double> dq, dn;
    // 1e-2 is reported, not gated: there the surface layer still moves by 20%, outside the limit's range.
    const std::vector<double> gammas{1e-2, 1e-3, 1e-4, 1e-5};
    for (const double g : gammas) {
        const auto p = solve_at(d, Model::dg, vfb + ov, g);
        REQUIRE(p.has_value());
        const Layer L = layer(d, *p);
        double worst = 0.0;
        for (std::size_t i = 0; i < L.n.size(); ++i) worst = std::max(worst, std::abs(L.n[i] - lc.n[i]));
        dq.push_back(std::abs(L.N_inv / lc.N_inv - 1.0));
        dn.push_back(worst / peak_cl);
        std::printf("QV8: gamma %.0e: |dN_inv / N_inv| %.4e, max |dn| / max n %.4e\n", g, dq.back(), dn.back());
    }
    for (std::size_t k = 1; k < gammas.size(); ++k) {
        const double oq = std::log(dq[k - 1] / dq[k]) / std::log(gammas[k - 1] / gammas[k]);
        const double on = std::log(dn[k - 1] / dn[k]) / std::log(gammas[k - 1] / gammas[k]);
        std::printf("QV8: order %.3f (N_inv), %.3f (n)\n", oq, on);
        if (k == 1) continue;  // the 1e-2 pair: reported
        CHECK(oq >= 0.9);
        CHECK(oq <= 1.1);
        CHECK(on >= 0.9);
        CHECK(on <= 1.1);
    }
}
