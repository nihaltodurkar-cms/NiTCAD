// Impact ionization in the drift-diffusion assembler (ARCHITECTURE.md section 11, Unit 19): the
// finite-difference Jacobian gate with the generation in 1D, 2D and 3D, both statistics and field
// mobility (the generation's part on its own too); the reconstruction: a uniform field and current
// at 0, 30 and 45 degrees to the grid ionize at alpha(|E|) |j| at every node, a field across the
// current does not ionize, and a current-free state does not; the temperature factor; and with the
// model off, or on for a material without coefficients, nothing changes.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

// A uniformly doped block, length_cm along each axis, ohmic contacts on x_min and x_max.
device::Device block(int D, int nodes, double donors, double acceptors,
                     physics::SemiconductorParameters material = physics::silicon_parameters,
                     double temperature_K = 300.0, double length_cm = 1e-4) {
    const auto a = uniform(0.0, length_cm, nodes);
    mesh::Mesh m = D == 1   ? *mesh::make_tensor_grid(a)
                   : D == 2 ? *mesh::make_tensor_grid(a, a)
                            : *mesh::make_tensor_grid(a, a, a);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = temperature_K,
         .regions = {{"semiconductor", *physics::Semiconductor::create(material)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, donors),
         .acceptors = std::vector<double>(n, acceptors),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}}});
}

struct Noise {
    std::uint64_t state;
    double next() {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        z ^= z >> 31;
        return static_cast<double>(z >> 11) * 0x1p-53 * 2.0 - 1.0;
    }
};

std::vector<std::vector<double>> dense(const linalg::SparseMatrix& a) {
    const auto n = static_cast<std::size_t>(a.rows());
    std::vector<std::vector<double>> d(n, std::vector<double>(n, 0.0));
    for (std::size_t r = 0; r < n; ++r) {
        for (auto k = static_cast<std::size_t>(a.row_offsets()[r]);
             k < static_cast<std::size_t>(a.row_offsets()[r + 1]); ++k) {
            d[r][static_cast<std::size_t>(a.col_indices()[k])] = a.values()[k];
        }
    }
    return d;
}

// A state with a field E0 [V/cm] at angle theta to x (in the x-y plane), n and p uniform, with
// relative noise on every unknown (0 for none). The ohmic nodes are not stamped: their rows are
// Dirichlet rows, and stamped values would put a jump of the potential beside them.
std::vector<double> sloped(const DriftDiffusion& dd, const device::Device& d,
                           const assemble::Scaling& s, double E0, double theta, double n_cm3,
                           double p_cm3, double noise, std::uint64_t seed) {
    std::vector<double> x(dd.unknowns());
    Noise z{seed};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        const auto& p = d.mesh().points()[i];
        const double along = p[0] * std::cos(theta) + (d.mesh().dimension() > 1 ? p[1] : 0.0) *
                                                          std::sin(theta);
        x[3 * i] = -E0 * along / s.V_T * (1.0 + noise * z.next());
        x[3 * i + 1] = n_cm3 / s.Ns * (1.0 + noise * z.next());
        x[3 * i + 2] = p_cm3 / s.Ns * (1.0 + noise * z.next());
    }
    return x;
}

// The finite-difference gate, every column, central differences, the worst error over the
// column's largest entry, each column's error the smallest at steps of 1e-5, 1e-6 and 1e-7
// max(|x|, 1) (a wrong entry is wrong at every step). The probe states carry 3e5 V/cm (potentials
// to 1160 V_T, edge drops of 116 V_T) and the generation puts large values in the continuity rows:
// no single step resolves every column. Measured with the house step 1e-7 alone: 2D Boltzmann
// 3e-6 (3.4e-8 at 1e-5 and 1e-6 alike: rounding); 3D Fermi-Dirac 6.5e-5 (6e-7 at 1e-5; the same
// with the model off); with field mobility 4e-6 at 1e-5 (4e-8 at 1e-6: truncation).
double fd_jacobian_error(const DriftDiffusion& dd, std::vector<double> x) {
    const std::size_t n = dd.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    auto jacobian = dd.make_jacobian();
    dd.evaluate(x, f, jacobian);
    const auto j = dense(jacobian);
    double worst = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        double scale = 1e-300;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][c]));
        double best = std::numeric_limits<double>::infinity();
        for (const double relative : {1e-5, 1e-6, 1e-7}) {
            const double base = x[c];
            const double step = relative * std::max(std::abs(base), 1.0);
            x[c] = base + step;
            const double up = x[c];
            dd.residual(x, fp);
            x[c] = base - step;
            const double down = x[c];
            dd.residual(x, fm);
            x[c] = base;
            double error = 0.0;
            for (std::size_t r = 0; r < n; ++r) {
                error = std::max(error, std::abs((fp[r] - fm[r]) / (up - down) - j[r][c]));
            }
            best = std::min(best, error / scale);
        }
        worst = std::max(worst, best);
    }
    return worst;
}

// The generation the residual carries at node i: (f_on - f_off) of its electron row, in
// cm^-3 s^-1 (the row takes V_i G / R0).
double generation(const std::vector<double>& on, const std::vector<double>& off, std::size_t i,
                  const device::Device& d, const assemble::Scaling& s) {
    const double V = d.mesh().volumes()[i] / std::pow(s.L_D, d.mesh().dimension());
    return (on[3 * i + 1] - off[3 * i + 1]) / V * s.R0;
}

}  // namespace

TEST_CASE("impact ionization assemble: the finite-difference Jacobian gate") {
    // A field of 3e5 V/cm (alpha_n 1.2e4 /cm, holes on their low branch, 25% from the switch)
    // with 5% noise on every unknown, at 0 and 30 degrees, n 1e16 and p 1e15 cm^-3: drift and
    // diffusion currents of both carriers, generation at every node.
    for (const int D : {1, 2, 3}) {
        const device::Device d = block(D, D == 3 ? 6 : 11, 1e16, 0.0);
        const auto s = *assemble::make_scaling(d);
        for (const bool fd : {false, true}) {
            for (const bool field : {false, true}) {
                const assemble::PhysicsModels models{
                    .field_mobility = field, .fermi_dirac = fd, .impact_ionization = true};
                const auto dd = *DriftDiffusion::create(d, s, models);
                for (const double theta : {0.0, std::numbers::pi / 6.0}) {
                    if (D == 1 && theta != 0.0) continue;
                    const auto x = sloped(dd, d, s, 3e5, theta, 1e16, 1e15, 0.05, 7);
                    const double error = fd_jacobian_error(dd, x);
                    CAPTURE(D, fd, field, theta, error);
                    REQUIRE(error <= 1e-6);  // measured at most 6.6e-7
                }
            }
        }
    }
}

TEST_CASE("impact ionization assemble: the generation part of the Jacobian on its own") {
    // With and without the model, so the generation entries cannot hide under the flux entries.
    for (const int D : {1, 2}) {
        const device::Device d = block(D, 11, 1e16, 0.0);
        const auto s = *assemble::make_scaling(d);
        const auto on = *DriftDiffusion::create(d, s, {.impact_ionization = true});
        const auto off = *DriftDiffusion::create(d, s);
        auto x = sloped(on, d, s, 3e5, D == 2 ? 0.4 : 0.0, 1e16, 1e15, 0.05, 11);
        const std::size_t n = on.unknowns();
        auto ja = on.make_jacobian(), jb = off.make_jacobian();
        std::vector<double> fa(n), fb(n), ua(n), ub(n), da(n), db(n);
        on.evaluate(x, fa, ja);
        off.evaluate(x, fb, jb);
        const auto A = dense(ja), B = dense(jb);
        double worst = 0.0, largest = 0.0;
        for (std::size_t c = 0; c < n; ++c) {
            const double base = x[c];
            const double step = 1e-5 * std::max(std::abs(base), 1.0);
            x[c] = base + step;
            const double hi = x[c];
            on.residual(x, ua);
            off.residual(x, ub);
            x[c] = base - step;
            const double lo = x[c];
            on.residual(x, da);
            off.residual(x, db);
            x[c] = base;
            for (std::size_t r = 0; r < n; ++r) {
                const double analytic = A[r][c] - B[r][c];
                const double numeric = ((ua[r] - ub[r]) - (da[r] - db[r])) / (hi - lo);
                largest = std::max(largest, std::abs(analytic));
                worst = std::max(worst, std::abs(numeric - analytic));
            }
        }
        CAPTURE(D, worst, largest);
        REQUIRE(largest > 0.0);
        REQUIRE(worst <= 1e-6 * largest);
    }
}

TEST_CASE("impact ionization assemble: a uniform field ionizes alike at any angle") {
    // n-silicon (1e16), a uniform field of 3e5 V/cm at 0, 30 and 45 degrees, n uniform: every edge
    // carries the drift current q mu n E . t exactly (Scharfetter-Gummel with a uniform density),
    // so the reconstruction gives E and j exactly and G = alpha_n(|E|) mu n |E| at every node off
    // the contacts. A per-edge model would give alpha_n(|E| cos) at 45 degrees, 3.6 times less.
    const device::Device d = block(2, 11, 1e16, 0.0);
    const auto s = *assemble::make_scaling(d);
    const auto on = *DriftDiffusion::create(d, s, {.impact_ionization = true});
    const auto off = *DriftDiffusion::create(d, s);
    const physics::Semiconductor si = physics::silicon();
    const double mu = physics::caughey_thomas_mobility(si, physics::Carrier::electron, 1e16, 300.0);
    const double E0 = 3e5, n0 = 1e16;
    const double expected =
        physics::impact_ionization_coefficient(
            physics::impact_ionization(si, physics::Carrier::electron), 1.0, E0)
            .alpha *
        mu * n0 * E0;
    for (const double theta : {0.0, std::numbers::pi / 6.0, std::numbers::pi / 4.0}) {
        // p negligible (1 cm^-3): the holes do not ionize measurably.
        const auto x = sloped(on, d, s, E0, theta, n0, 1.0, 0.0, 0);
        std::vector<double> fa(on.unknowns()), fb(on.unknowns());
        on.residual(x, fa);
        off.residual(x, fb);
        double worst = 0.0;
        std::size_t checked = 0;
        for (std::size_t i = 0; i < on.node_count(); ++i) {
            const double px = d.mesh().points()[i][0];
            if (px == 0.0 || px == 1e-4) {  // the contacts: their Dirichlet rows do not generate
                for (int c = 0; c < 3; ++c) REQUIRE(fa[3 * i + c] == fb[3 * i + c]);
                continue;
            }
            worst = std::max(worst, std::abs(generation(fa, fb, i, d, s) / expected - 1.0));
            ++checked;
        }
        CAPTURE(theta, worst, checked);
        REQUIRE(checked == 99);
        REQUIRE(worst < 1e-9);
    }
}

TEST_CASE("impact ionization assemble: a field across the current does not ionize") {
    // A field of 3e5 V/cm along x (0.1 um block), the electrons in equilibrium along x (n
    // proportional to e^psi, so no x current) and diffusing along y (n rising with y), the holes in
    // equilibrium: the current is perpendicular to the field and the generation is exactly zero.
    const device::Device d = block(2, 11, 1e16, 0.0, physics::silicon_parameters, 300.0, 1e-5);
    const auto s = *assemble::make_scaling(d);
    const auto on = *DriftDiffusion::create(d, s, {.impact_ionization = true});
    const auto off = *DriftDiffusion::create(d, s);
    std::vector<double> x(on.unknowns());
    for (std::size_t i = 0; i < on.node_count(); ++i) {
        const auto& p = d.mesh().points()[i];
        const double psi = -3e5 * p[0] / s.V_T;
        x[3 * i] = psi;
        x[3 * i + 1] = std::exp(psi) * (1.0 + 0.5 * p[1] / 1e-5);
        x[3 * i + 2] = 1e-40 * std::exp(-psi);
    }
    std::vector<double> fa(on.unknowns()), fb(on.unknowns());
    on.residual(x, fa);
    off.residual(x, fb);
    double largest = 0.0, current = 0.0;
    for (std::size_t i = 0; i < on.node_count(); ++i) {
        largest = std::max(largest, std::abs(fa[3 * i + 1] - fb[3 * i + 1]));
    }
    for (const auto& [jn, jp] : on.edge_currents(x)) current = std::max(current, std::abs(jn));
    // What the same current would generate along the field: V alpha_n(3e5) |j| / (q R0) in the
    // rows' units, |j| the largest y current density (J0 times the scaled flux over the scaled
    // coupling area).
    double j_max = 0.0;
    const auto currents = on.edge_currents(x);
    for (std::size_t k = 0; k < currents.size(); ++k) {
        const double area = d.mesh().edges()[k].coupling_area / s.L_D;
        j_max = std::max(j_max, s.J0 * std::abs(currents[k].first) / area);
    }
    const double alpha =
        physics::impact_ionization_coefficient(
            physics::impact_ionization(physics::silicon(), physics::Carrier::electron), 1.0, 3e5)
            .alpha;
    const double V = d.mesh().volumes()[60] / (s.L_D * s.L_D);
    const double aligned = V * alpha * j_max / base::q_C / s.R0;
    CAPTURE(largest, current, aligned, largest / aligned);
    REQUIRE(current > 0.0);  // the y diffusion current flows
    // The x fluxes are rounding (n e^-psi constant along x), so E . j is rounding against
    // |E| |j|: measured 4e-10 of the aligned generation.
    REQUIRE(largest <= 1e-8 * aligned);
}

TEST_CASE("impact ionization assemble: no current, no generation, and finite partials") {
    // Thermal equilibrium of a 1e19 / 1e19 junction (a built-in field near 1e6 V/cm, where
    // alpha_n is 2e5 /cm): the currents are rounding, so the generation is negligible against the
    // recombination terms' scale and every Jacobian entry is finite.
    const auto a = uniform(0.0, 2e-5, 81);
    mesh::Mesh m = *mesh::make_tensor_grid(a);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-5 ? acceptors : donors)[i] = 1e19;
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"p", device::ContactKind::ohmic, std::move(left)},
                      {"n", device::ContactKind::ohmic, std::move(right)}}});
    const auto s = *assemble::make_scaling(d);
    const auto on = *DriftDiffusion::create(d, s, {.impact_ionization = true});
    const auto off = *DriftDiffusion::create(d, s);
    // The equilibrium carriers of a smooth potential through the junction (tanh profile).
    std::vector<double> psi(n);
    for (std::size_t i = 0; i < n; ++i) {
        psi[i] = 42.0 * std::tanh((d.mesh().points()[i][0] - 1e-5) / 2e-7);
    }
    const auto x = on.state_from_potential(psi);
    auto ja = on.make_jacobian();
    std::vector<double> fa(on.unknowns()), fb(on.unknowns());
    on.evaluate(x, fa, ja);
    off.residual(x, fb);
    for (const double v : ja.values()) REQUIRE(std::isfinite(v));
    double G = 0.0, peak_E = 0.0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        G = std::max(G, std::abs(generation(fa, fb, i, d, s)));
        peak_E = std::max(peak_E, s.V_T * std::abs(psi[i + 1] - psi[i]) /
                                      (d.mesh().points()[i + 1][0] - d.mesh().points()[i][0]));
    }
    // The equilibrium currents are the rounding of the opposing fluxes (no current in the model):
    // the generation is bounded by alpha at the peak field times them, and is negligible against
    // what a majority drift current q mu N E would generate there (measured 6e-10 of it).
    double j_max = 0.0;
    const auto currents = on.edge_currents(x);
    for (const auto& [jn, jp] : currents) {
        j_max = std::max(j_max, s.J0 * (std::abs(jn) + std::abs(jp)));  // 1D: area 1
    }
    const physics::Semiconductor si = physics::silicon();
    const auto rate = [&](physics::Carrier c) {
        return physics::impact_ionization_coefficient(physics::impact_ionization(si, c), 1.0,
                                                      peak_E)
            .alpha;
    };
    const double alpha_max =
        std::max(rate(physics::Carrier::electron), rate(physics::Carrier::hole));
    const double drift = 1e19 * 100.0 * peak_E;  // mu N E with mu 100 cm^2/Vs, per q
    CAPTURE(G, peak_E, j_max, alpha_max, drift, G / (alpha_max * drift));
    REQUIRE(peak_E > 5e5);
    REQUIRE(G <= alpha_max * j_max / base::q_C * (1.0 + 1e-9));
    REQUIRE(G <= 1e-8 * alpha_max * drift);
}

TEST_CASE("impact ionization assemble: the temperature factor") {
    // At 400 K the generation is alpha_n(E; gamma(400 K)) mu(400 K) n |E|.
    const physics::Semiconductor si = physics::silicon();
    const device::Device d = block(1, 21, 1e16, 0.0, physics::silicon_parameters, 400.0);
    const auto s = *assemble::make_scaling(d);
    const auto on = *DriftDiffusion::create(d, s, {.impact_ionization = true});
    const auto off = *DriftDiffusion::create(d, s);
    const auto x = sloped(on, d, s, 3e5, 0.0, 1e16, 1.0, 0.0, 0);
    std::vector<double> fa(on.unknowns()), fb(on.unknowns());
    on.residual(x, fa);
    off.residual(x, fb);
    const double gamma = physics::impact_ionization_temperature_factor(0.063, 400.0);
    const double mu = physics::caughey_thomas_mobility(si, physics::Carrier::electron, 1e16, 400.0);
    const double expected =
        physics::impact_ionization_coefficient(
            physics::impact_ionization(si, physics::Carrier::electron), gamma, 3e5)
            .alpha *
        mu * 1e16 * 3e5;
    const double G = generation(fa, fb, 10, d, s);
    CAPTURE(G, expected, gamma);
    REQUIRE(gamma > 1.0);
    REQUIRE(std::abs(G / expected - 1.0) < 1e-9);
}

TEST_CASE("impact ionization assemble: off, or without coefficients, nothing changes") {
    const device::Device d = block(2, 9, 1e16, 0.0);
    const auto s = *assemble::make_scaling(d);
    const auto off = *DriftDiffusion::create(d, s);
    const auto x = sloped(off, d, s, 3e5, 0.3, 1e16, 1e15, 0.05, 3);
    // The default is off.
    REQUIRE_FALSE(assemble::PhysicsModels{}.impact_ionization);
    // A material without coefficients: the model on changes no value and no pattern.
    physics::SemiconductorParameters none = physics::silicon_parameters;
    none.impact_ionization = {};
    const device::Device dn = block(2, 9, 1e16, 0.0, none);
    const auto a = *DriftDiffusion::create(dn, s);
    const auto b = *DriftDiffusion::create(dn, s, {.impact_ionization = true});
    auto ja = a.make_jacobian(), jb = b.make_jacobian();
    std::vector<double> fa(a.unknowns()), fb(a.unknowns());
    a.evaluate(x, fa, ja);
    b.evaluate(x, fb, jb);
    REQUIRE(fa == fb);
    // The pattern gains the cross entries (n row, neighbour p; p row, neighbour n), all zero.
    REQUIRE(jb.nonzeros() > ja.nonzeros());
    const auto A = dense(ja), B = dense(jb);
    REQUIRE(A == B);
}


TEST_CASE("impact ionization assemble: the current resolution's partials") {
    // eps_c^2 = floor^2 + (rho R_c)^2 depends on the densities and potentials of the node's edges.
    // At the default rho (0) only the floor remains, negligible at these probe states; with
    // rho = 0.5 it is comparable to the currents, and the finite-difference gate checks its
    // partials too.
    for (const int D : {1, 2}) {
        const device::Device d = block(D, 11, 1e16, 0.0);
        const auto s = *assemble::make_scaling(d);
        for (const double rho : {0.5, 0.05}) {
            const auto dd = *DriftDiffusion::create(
                d, s, {.impact_ionization = true, .impact_current_resolution = rho});
            const auto x = sloped(dd, d, s, 3e5, D == 2 ? 0.4 : 0.0, 1e16, 1e15, 0.05, 29);
            const double error = fd_jacobian_error(dd, x);
            CAPTURE(D, rho, error);
            REQUIRE(error <= 1e-6);
        }
    }
    // The resolution changes the generation: it is smaller with a coarser one.
    const device::Device d = block(1, 11, 1e16, 0.0);
    const auto s = *assemble::make_scaling(d);
    const auto fine = *DriftDiffusion::create(d, s, {.impact_ionization = true});
    const auto coarse = *DriftDiffusion::create(
        d, s, {.impact_ionization = true, .impact_current_resolution = 0.5});
    const auto off = *DriftDiffusion::create(d, s);
    const auto x = sloped(fine, d, s, 3e5, 0.0, 1e16, 1e15, 0.05, 29);
    std::vector<double> ff(fine.unknowns()), fc(fine.unknowns()), fo(fine.unknowns());
    fine.residual(x, ff);
    coarse.residual(x, fc);
    off.residual(x, fo);
    REQUIRE(generation(fc, fo, 5, d, s) < generation(ff, fo, 5, d, s));
    REQUIRE(generation(fc, fo, 5, d, s) > 0.0);
    // Its value: a uniform field E0 and uniform n (p negligible) give every edge the drift current
    // j = q mu n E0 and R = 2 j sqrt(1 + Delta^2) / Delta, Delta = E0 h / V_T, so at a node off the
    // contacts G = alpha(E0 j / m) j^2 / m / q with m = sqrt(j^2 + (rho R)^2).
    {
        const double E0 = 3e5, n0 = 1e16, rho = 0.5, h = 1e-5;
        const auto u = sloped(coarse, d, s, E0, 0.0, n0, 1.0, 0.0, 0);
        coarse.residual(u, fc);
        off.residual(u, fo);
        const physics::Semiconductor si = physics::silicon();
        const double mu =
            physics::caughey_thomas_mobility(si, physics::Carrier::electron, 1e16, 300.0);
        const double j = base::q_C * mu * n0 * E0;
        const double delta = E0 * h / s.V_T;
        const double R = 2.0 * j * std::sqrt(1.0 + delta * delta) / delta;
        const double m = std::sqrt(j * j + rho * rho * R * R);
        const double expected =
            physics::impact_ionization_coefficient(
                physics::impact_ionization(si, physics::Carrier::electron), 1.0, E0 * j / m)
                .alpha *
            j * j / m / base::q_C;
        for (const std::size_t i : {std::size_t{3}, std::size_t{5}, std::size_t{8}}) {
            CAPTURE(i, generation(fc, fo, i, d, s), expected);
            REQUIRE(std::abs(generation(fc, fo, i, d, s) / expected - 1.0) < 1e-9);
        }
    }
    // Invalid with the model on; not read with it off.
    REQUIRE(DriftDiffusion::create(d, s, {.impact_ionization = true,
                                          .impact_current_resolution = -1.0})
                .error()
                .code == base::ErrorCode::invalid_input);
    REQUIRE(DriftDiffusion::create(d, s, {.impact_current_resolution = -1.0}).has_value());
}
