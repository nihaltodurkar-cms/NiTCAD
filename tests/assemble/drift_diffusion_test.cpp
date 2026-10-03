// Coupled drift-diffusion residual and Jacobian (ARCHITECTURE.md section 11, Unit 9; FD-Jacobian
// gate of section 10 with the legacy probe), D-generic reduction, and the Newton update hooks.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/bernoulli.hpp"
#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using base::ErrorCode;

namespace {

std::vector<double> graded_axis(int per_side) {
    std::vector<double> half{0.0};
    const double h = 1e-4 * (1.15 - 1.0) / (std::pow(1.15, per_side) - 1.0);
    for (int k = 0; k < per_side; ++k) {
        half.push_back(half.back() + h * std::pow(1.15, per_side - 1 - k));
    }
    half.back() = 1e-4;
    std::vector<double> x = half;
    for (int k = per_side - 1; k >= 0; --k) x.push_back(2e-4 - half[static_cast<std::size_t>(k)]);
    return x;
}

std::vector<double> uniform_axis(double length, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = length * i / (nodes - 1);
    return a;
}

device::Device diode(mesh::Mesh m, double NA = 1e17, double ND = 1e17) {
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 1e-4) {
            acceptors[i] = NA;
        } else {
            donors[i] = ND;
        }
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
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

// A biased system at the legacy probe's state: Boltzmann carriers at the charge-neutral potential
// with the anode at +0.3 V, then psi + 0.02 noise and n, p times (1 + 0.01 noise).
std::pair<DriftDiffusion, std::vector<double>> probe_state(
    const device::Device& d, std::uint64_t seed, const assemble::PhysicsModels& models = {}) {
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling, models);
    const std::vector<double> bias{0.3, 0.0};
    REQUIRE(system.set_bias(bias).has_value());
    const auto poisson = *assemble::EquilibriumPoisson::create(d, scaling);
    auto x = system.state_from_potential(poisson.charge_neutral_potential());
    Noise noise{seed};
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[3 * i] += 0.02 * noise.next();
        x[3 * i + 1] *= 1.0 + 0.01 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.01 * noise.next();
    }
    return {std::move(system), std::move(x)};
}

// Legacy _jacobian_probe: central differences, step 1e-7 max(|u|, 1), every column; per column
// max |fd - J| over the column's largest |J| (+1e-30).
double fd_jacobian_error(const DriftDiffusion& system, std::vector<double> x) {
    const std::size_t n = system.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    auto jacobian = system.make_jacobian();
    system.evaluate(x, f, jacobian);
    const auto j = dense(jacobian);
    double worst = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        const double base = x[c];
        const double step = 1e-7 * std::max(std::abs(base), 1.0);
        x[c] = base + step;
        const double up = x[c];
        system.residual(x, fp);
        x[c] = base - step;
        const double down = x[c];
        system.residual(x, fm);
        x[c] = base;
        double scale = 1e-30, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((fp[r] - fm[r]) / (up - down) - j[r][c]));
        }
        worst = std::max(worst, error / scale);
    }
    return worst;
}

}  // namespace

TEST_CASE("drift-diffusion: FD-Jacobian gate in 1D, 2D and 3D (section 10, <= 5e-5)") {
    const auto x = graded_axis(30);  // 61 nodes, 183 columns
    const auto xs = graded_axis(10);
    const auto y = uniform_axis(1e-4, 3);
    const auto [s1, x1] = probe_state(diode(*mesh::make_tensor_grid(x)), 42);
    const auto [s2, x2] = probe_state(diode(*mesh::make_tensor_grid(xs, y)), 43);
    const auto [s3, x3] = probe_state(diode(*mesh::make_tensor_grid(xs, y, y)), 44);
    const double e1 = fd_jacobian_error(s1, x1);
    const double e2 = fd_jacobian_error(s2, x2);
    const double e3 = fd_jacobian_error(s3, x3);
    UNSCOPED_INFO("worst column errors: 1D " << e1 << ", 2D " << e2 << ", 3D " << e3);
    REQUIRE(s1.unknowns() >= 80);
    REQUIRE(e1 <= 5e-5);
    REQUIRE(e2 <= 5e-5);
    REQUIRE(e3 <= 5e-5);
}

TEST_CASE("drift-diffusion: a y-uniform 2D residual is the 1D residual times the scaled width") {
    const auto x = graded_axis(15);
    const auto y = std::vector<double>{0.0, 2e-5, 3e-5, 7e-5};
    const auto d1 = diode(*mesh::make_tensor_grid(x));
    const auto d2 = diode(*mesh::make_tensor_grid(x, y));
    const auto [s1, x1] = probe_state(d1, 7);
    const auto scaling = *assemble::make_scaling(d2);
    auto s2 = *DriftDiffusion::create(d2, scaling);
    const std::vector<double> bias{0.3, 0.0};
    REQUIRE(s2.set_bias(bias).has_value());
    const std::size_t nx = x.size();
    std::vector<double> x2(s2.unknowns());
    for (std::size_t k = 0; k < s2.node_count(); ++k) {
        for (std::size_t c = 0; c < 3; ++c) x2[3 * k + c] = x1[3 * (k % nx) + c];
    }
    std::vector<double> f1(s1.unknowns()), f2(s2.unknowns());
    s1.residual(x1, f1);
    s2.residual(x2, f2);
    double worst = 0.0;
    for (std::size_t k = 0; k < s2.node_count(); ++k) {
        const std::size_t i = k % nx;
        if (i == 0 || i + 1 == nx) continue;  // contact rows are Dirichlet, not scaled
        const double width = d2.mesh().volumes()[k] / d1.mesh().volumes()[i] / scaling.L_D;
        for (std::size_t c = 0; c < 3; ++c) {
            const double expected = width * f1[3 * i + c];
            worst = std::max(worst, std::abs(f2[3 * k + c] - expected) / std::abs(expected));
        }
    }
    REQUIRE(worst <= 1e-12);
}

TEST_CASE("drift-diffusion: update hooks measure the full correction and keep densities positive") {
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(5)));
    auto system = *DriftDiffusion::create(d, *assemble::make_scaling(d));
    std::vector<double> x(system.unknowns(), 1.0);
    std::vector<double> dx(system.unknowns(), 0.0);
    dx[0] = 7.0;   // psi: clipped to 5
    dx[1] = -3.0;  // n: 1 - 3 < 0, clamped to 0.1
    dx[5] = 20.0;  // p of node 1: clamped to 10
    REQUIRE(system.update_size(x, dx) == 20.0);  // the full relative correction, 20 / 1
    system.apply_update(x, dx, 5.0);
    REQUIRE(x[0] == 6.0);
    REQUIRE(x[1] == 0.1);
    REQUIRE(x[5] == 10.0);
    // A non-positive density cannot converge.
    x[4] = 0.0;
    REQUIRE(system.update_size(x, dx) == std::numeric_limits<double>::infinity());
    // Unit 15: a density below 1e-20 of the largest is measured against that floor (the linear
    // solve does not resolve it relative to itself); above it, against itself.
    std::vector<double> y(system.unknowns(), 1.0);
    std::vector<double> dy(system.unknowns(), 0.0);
    y[7] = 1e-30;
    dy[7] = 1e-29;
    REQUIRE(system.update_size(y, dy) == 1e-29 / 1e-20);
    y[7] = 1e-15;
    dy[7] = 1e-16;
    REQUIRE(system.update_size(y, dy) == 1e-16 / 1e-15);
}

TEST_CASE("drift-diffusion: biases are validated and stamped on the contacts") {
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(5)));
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling);
    const std::vector<double> one{0.1};
    REQUIRE(system.set_bias(one).error().code == ErrorCode::invalid_input);
    const std::vector<double> infinite{INFINITY, 0.0};
    REQUIRE(system.set_bias(infinite).error().code == ErrorCode::invalid_input);
    const std::vector<double> bias{0.25, -0.5};
    REQUIRE(system.set_bias(bias).has_value());
    const std::vector<double> psi(system.node_count(), 0.0);
    const auto x = system.state_from_potential(psi);
    const auto poisson = *assemble::EquilibriumPoisson::create(d, scaling);
    const auto psi0 = poisson.contact_potential();
    const std::size_t last = system.node_count() - 1;
    REQUIRE(std::abs(x[0] - (psi0[0] + 0.25 / scaling.V_T)) <= 1e-14 * std::abs(x[0]));
    const double cathode = psi0[last] - 0.5 / scaling.V_T;
    REQUIRE(std::abs(x[3 * last] - cathode) <= 1e-14 * std::abs(cathode));
    // Densities on the contacts stay at the neutral values whatever the bias.
    REQUIRE(std::abs(x[2] - 1.0) <= 1e-13);  // p = N_A / Ns + n_ie^2 / N_A on the anode
}

TEST_CASE("drift-diffusion: the recombination part of the Jacobian matches finite differences") {
    // The column-normalized gate cannot see it: in a density column the Scharfetter-Gummel
    // coefficients (about D / h^2) dwarf dR/dn (measured: dropping dR/dn from the Jacobian passes
    // that gate). So the SRH contribution is checked on its own: J(srh) - J(no srh) against finite
    // differences of F(srh) - F(no srh). The state is uniform but out of equilibrium
    // (n p != n_ie^2), so every flux is exactly zero there and the subtraction does not cancel
    // large flux sums.
    // What is left grows with an / (k dR/dn V), the flux coefficient over the recombination term:
    // 0.1 um cells keep it small (on the 1e-8 cm junction cells of a graded mesh it reached 7e-5).
    const auto d = diode(*mesh::make_tensor_grid(uniform_axis(2e-4, 21)));
    const auto scaling = *assemble::make_scaling(d);
    auto with = *DriftDiffusion::create(d, scaling);
    auto without = *DriftDiffusion::create(d, scaling, {.srh = false});
    std::vector<double> x(with.unknowns());
    for (std::size_t i = 0; i < with.node_count(); ++i) {
        x[3 * i] = 0.1;
        x[3 * i + 1] = 0.5;
        x[3 * i + 2] = 2e-3;
    }
    const std::size_t n = with.unknowns();
    std::vector<double> f(n), g(n);
    auto ja = with.make_jacobian(), jb = without.make_jacobian();
    with.evaluate(x, f, ja);
    without.evaluate(x, g, jb);
    const auto a = dense(ja), b = dense(jb);
    const auto difference = [&](const std::vector<double>& u, std::vector<double>& out) {
        with.residual(u, f);
        without.residual(u, g);
        for (std::size_t r = 0; r < n; ++r) out[r] = f[r] - g[r];
    };
    std::vector<double> up(n), down(n);
    double worst = 0.0, largest = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        const double base = x[c];
        const double step = 1e-6 * std::max(std::abs(base), 1e-12);
        x[c] = base + step;
        const double hi = x[c];
        difference(x, up);
        x[c] = base - step;
        const double lo = x[c];
        difference(x, down);
        x[c] = base;
        double scale = 1e-300, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(a[r][c] - b[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((up[r] - down[r]) / (hi - lo) - (a[r][c] - b[r][c])));
        }
        largest = std::max(largest, scale);
        if (scale > 1e-300) worst = std::max(worst, error / scale);
    }
    UNSCOPED_INFO("worst recombination column error " << worst);
    REQUIRE(largest > 0.0);
    REQUIRE(worst <= 1e-5);  // measured 7.6e-7; dropping dR/dn gives about 1
}

TEST_CASE("drift-diffusion: the 1D residual is the legacy rows, written out") {
    // Legacy device1d.cpp residual (baseline models: doping mobility, SRH, Auger, band-gap
    // narrowing) in scaled units, on an asymmetric diode so the edge diffusivity's harmonic mean
    // matters at the junction and the 1e18 side has band-gap narrowing:
    //   an_k = hmean(mu_n) V_T / D0 / h_k, Jn_k = an_k (n[k+1] B(d_k) - n[k] B(-d_k)),
    //   Jp_k = -ap_k (p[k+1] B(-dp_k) - p[k] B(dp_k)),
    //   d_k = psi[k+1] - psi[k] + ln(nie[k+1] / nie[k]), dp_k = psi[k+1] - psi[k] - ln(...),
    //   R = recombination_boltzmann(n Ns, p Ns, nie, ..., auger) / R0,
    //   F_n = Jn[i] - Jn[i-1] - R dV,  F_p = Jp[i] - Jp[i-1] + R dV.
    const auto x = graded_axis(15);
    const auto d = diode(*mesh::make_tensor_grid(x), 1e18, 1e16);
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling);
    const std::vector<double> bias{0.3, 0.0};
    REQUIRE(system.set_bias(bias).has_value());
    const auto poisson = *assemble::EquilibriumPoisson::create(d, scaling);
    auto u = system.state_from_potential(poisson.charge_neutral_potential());
    Noise noise{3};
    for (std::size_t i = 1; i + 1 < x.size(); ++i) {
        u[3 * i] += 0.5 * noise.next();
        u[3 * i + 1] *= 1.0 + 0.3 * noise.next();
        u[3 * i + 2] *= 1.0 + 0.3 * noise.next();
    }
    std::vector<double> f(system.unknowns());
    system.residual(u, f);

    const physics::Semiconductor si = physics::silicon();
    const double T = 300.0, VT = scaling.V_T, Ns = scaling.Ns;
    const std::size_t N = x.size();
    std::vector<double> mun(N), mup(N), taun(N), taup(N), nie(N);
    for (std::size_t i = 0; i < N; ++i) {
        const double Nt = x[i] < 1e-4 ? 1e18 : 1e16;
        nie[i] = physics::effective_intrinsic_density(si, Nt, T);
        mun[i] = physics::caughey_thomas_mobility(si, physics::Carrier::electron, Nt, T);
        mup[i] = physics::caughey_thomas_mobility(si, physics::Carrier::hole, Nt, T);
        taun[i] = physics::scharfetter_lifetime(si, physics::Carrier::electron, Nt);
        taup[i] = physics::scharfetter_lifetime(si, physics::Carrier::hole, Nt);
    }
    std::vector<double> Jn(N - 1), Jp(N - 1), h(N - 1);
    for (std::size_t k = 0; k + 1 < N; ++k) {
        h[k] = (x[k + 1] - x[k]) / scaling.L_D;
        const double dlnnie = std::log(nie[k + 1] / nie[k]);
        const double dk = u[3 * (k + 1)] - u[3 * k] + dlnnie;
        const double dpk = u[3 * (k + 1)] - u[3 * k] - dlnnie;
        const double an = 2 * mun[k] * mun[k + 1] / (mun[k] + mun[k + 1]) * VT / h[k];
        const double ap = 2 * mup[k] * mup[k + 1] / (mup[k] + mup[k + 1]) * VT / h[k];
        Jn[k] = an * (u[3 * (k + 1) + 1] * assemble::bernoulli(dk) -
                      u[3 * k + 1] * assemble::bernoulli(-dk));
        Jp[k] = -ap * (u[3 * (k + 1) + 2] * assemble::bernoulli(-dpk) -
                       u[3 * k + 2] * assemble::bernoulli(dpk));
    }
    for (std::size_t i = 1; i + 1 < N; ++i) {
        const double dV = 0.5 * (h[i - 1] + h[i]);
        const double n = u[3 * i + 1] * Ns, p = u[3 * i + 2] * Ns;
        const double excess = n * p - nie[i] * nie[i];
        const double srh = excess / (taup[i] * (n + nie[i]) + taun[i] * (p + nie[i]));
        const double auger = (2.8e-31 * n + 9.9e-32 * p) * excess;
        const double R = (srh + auger) / scaling.R0;
        const double expected_n = Jn[i] - Jn[i - 1] - R * dV;
        const double expected_p = Jp[i] - Jp[i - 1] + R * dV;
        const double scale_n = std::max({std::abs(Jn[i]), std::abs(Jn[i - 1]), std::abs(R * dV)});
        const double scale_p = std::max({std::abs(Jp[i]), std::abs(Jp[i - 1]), std::abs(R * dV)});
        CAPTURE(i);
        REQUIRE(std::abs(f[3 * i + 1] - expected_n) <= 1e-12 * scale_n);
        REQUIRE(std::abs(f[3 * i + 2] - expected_p) <= 1e-12 * scale_p);
    }
}

TEST_CASE("drift-diffusion: the Auger part of the Jacobian matches finite differences") {
    // As the recombination test above, for Auger: J(auger) - J(no auger) against finite
    // differences of the residual difference, at a uniform state of 1e19 cm^-3 carriers where
    // Auger dominates (n = p = 100 Ns).
    const auto d = diode(*mesh::make_tensor_grid(uniform_axis(2e-4, 21)));
    const auto scaling = *assemble::make_scaling(d);
    auto with = *DriftDiffusion::create(d, scaling);
    auto without = *DriftDiffusion::create(d, scaling, {.auger = false});
    std::vector<double> x(with.unknowns());
    for (std::size_t i = 0; i < with.node_count(); ++i) {
        x[3 * i] = 0.1;
        x[3 * i + 1] = 100.0;
        x[3 * i + 2] = 100.0;
    }
    const std::size_t n = with.unknowns();
    std::vector<double> f(n), g(n), up(n), down(n);
    auto ja = with.make_jacobian(), jb = without.make_jacobian();
    with.evaluate(x, f, ja);
    without.evaluate(x, g, jb);
    const auto a = dense(ja), b = dense(jb);
    const auto difference = [&](const std::vector<double>& u, std::vector<double>& out) {
        with.residual(u, f);
        without.residual(u, g);
        for (std::size_t r = 0; r < n; ++r) out[r] = f[r] - g[r];
    };
    double worst = 0.0, largest = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        const double base = x[c];
        const double step = 1e-6 * std::max(std::abs(base), 1e-12);
        x[c] = base + step;
        const double hi = x[c];
        difference(x, up);
        x[c] = base - step;
        const double lo = x[c];
        difference(x, down);
        x[c] = base;
        double scale = 0.0, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(a[r][c] - b[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((up[r] - down[r]) / (hi - lo) - (a[r][c] - b[r][c])));
        }
        largest = std::max(largest, scale);
        if (scale > 0.0) worst = std::max(worst, error / scale);
    }
    UNSCOPED_INFO("worst Auger column error " << worst);
    REQUIRE(largest > 0.0);
    REQUIRE(worst <= 1e-5);
}

TEST_CASE("drift-diffusion: FD-Jacobian gate on a heavily doped diode, all models on") {
    // 1e19 / 1e18: band-gap narrowing on both sides (so ln n_ie varies), Auger significant.
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(30)), 1e19, 1e18);
    const auto [system, x] = probe_state(d, 77);
    const double e = fd_jacobian_error(system, x);
    UNSCOPED_INFO("worst column error " << e);
    REQUIRE(e <= 5e-5);
}

// Unit 13: Canali field-dependent mobility of each edge.

TEST_CASE("drift-diffusion: FD-Jacobian gate with field mobility in 1D, 2D and 3D") {
    // The probe's psi noise (0.02) gives edge fields up to about 1e5 V/cm on the junction cells,
    // where the Canali factor is far from 1.
    const assemble::PhysicsModels fm{.field_mobility = true};
    const auto x = graded_axis(30);
    const auto xs = graded_axis(10);
    const auto y = uniform_axis(1e-4, 3);
    const auto [s1, x1] = probe_state(diode(*mesh::make_tensor_grid(x)), 42, fm);
    const auto [s2, x2] = probe_state(diode(*mesh::make_tensor_grid(xs, y)), 43, fm);
    const auto [s3, x3] = probe_state(diode(*mesh::make_tensor_grid(xs, y, y)), 44, fm);
    const double e1 = fd_jacobian_error(s1, x1);
    const double e2 = fd_jacobian_error(s2, x2);
    const double e3 = fd_jacobian_error(s3, x3);
    UNSCOPED_INFO("worst column errors: 1D " << e1 << ", 2D " << e2 << ", 3D " << e3);
    REQUIRE(e1 <= 5e-5);
    REQUIRE(e2 <= 5e-5);
    REQUIRE(e3 <= 5e-5);
}

TEST_CASE("drift-diffusion: the field-mobility part of the Jacobian matches finite differences") {
    // As for recombination, the column-normalized gate may hide a small term, so the part the
    // field dependence adds is checked on its own: J(fm) - J(no fm) against finite differences of
    // F(fm) - F(no fm), per column relative to that column's largest difference. Only the potential
    // columns: the field factor depends on psi alone, and the fluxes are linear in the densities
    // (their columns are the low-field partials times the factor, covered by the full gate; a
    // minority density near 1e-14 is below what a difference of residuals resolves).
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(20)));
    const auto [with, x0] = probe_state(d, 5, {.field_mobility = true});
    const auto [without, unused] = probe_state(d, 5);
    std::vector<double> x = x0;
    const std::size_t n = with.unknowns();
    std::vector<double> f(n), g(n);
    auto ja = with.make_jacobian(), jb = without.make_jacobian();
    with.evaluate(x, f, ja);
    without.evaluate(x, g, jb);
    const auto a = dense(ja), b = dense(jb);
    const auto difference = [&](const std::vector<double>& u, std::vector<double>& out) {
        with.residual(u, f);
        without.residual(u, g);
        for (std::size_t r = 0; r < n; ++r) out[r] = f[r] - g[r];
    };
    std::vector<double> up(n), down(n);
    double worst = 0.0, largest = 0.0;
    for (std::size_t c = 0; c < n; c += 3) {
        const double base = x[c];
        const double step = 1e-7 * std::max(std::abs(base), 1.0);
        x[c] = base + step;
        const double hi = x[c];
        difference(x, up);
        x[c] = base - step;
        const double lo = x[c];
        difference(x, down);
        x[c] = base;
        double scale = 1e-300, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(a[r][c] - b[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((up[r] - down[r]) / (hi - lo) - (a[r][c] - b[r][c])));
        }
        largest = std::max(largest, scale);
        if (scale > 1e-300) worst = std::max(worst, error / scale);
    }
    UNSCOPED_INFO("worst field-mobility column error " << worst);
    REQUIRE(largest > 0.0);
    REQUIRE(worst <= 1e-5);
}

TEST_CASE("drift-diffusion: each edge's current is scaled by its Canali factor") {
    // The edge factor is the low-field one times mu_C(mu0, E) / mu0, with mu0 the harmonic mean of
    // the end nodes' Caughey-Thomas mobilities and E = V_T |psi_b - psi_a| / length.
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(15)), 1e18, 1e16);
    const auto scaling = *assemble::make_scaling(d);
    const auto [with, x] = probe_state(d, 9, {.field_mobility = true});
    const auto [without, unused] = probe_state(d, 9);
    const auto a = with.edge_currents(x), b = without.edge_currents(x);
    const physics::Semiconductor si = physics::silicon();
    const auto edges = d.mesh().edges();
    double worst = 0.0, smallest_factor = 1.0;
    for (std::size_t k = 0; k < edges.size(); ++k) {
        const auto i = edges[k].first, j = edges[k].second;
        const double E = scaling.V_T *
                         std::abs(x[3 * static_cast<std::size_t>(j)] -
                                  x[3 * static_cast<std::size_t>(i)]) /
                         edges[k].length;
        for (const auto carrier : {physics::Carrier::electron, physics::Carrier::hole}) {
            const double mi =
                physics::caughey_thomas_mobility(si, carrier, d.total_impurity(i), 300.0);
            const double mj =
                physics::caughey_thomas_mobility(si, carrier, d.total_impurity(j), 300.0);
            const double mu0 = 2.0 * mi * mj / (mi + mj);
            const double factor =
                physics::canali_mobility(mu0, E, physics::saturation(si, carrier)).mobility / mu0;
            smallest_factor = std::min(smallest_factor, factor);
            const double ratio = carrier == physics::Carrier::electron ? a[k].first / b[k].first
                                                                       : a[k].second / b[k].second;
            worst = std::max(worst, std::abs(ratio / factor - 1.0));
        }
    }
    CAPTURE(worst, smallest_factor);
    REQUIRE(smallest_factor < 0.5);  // the probe reaches well into saturation
    REQUIRE(worst <= 1e-13);
}

// Unit 14: Fermi-Dirac statistics.

TEST_CASE("drift-diffusion: FD-Jacobian gate under Fermi-Dirac statistics in 1D, 2D and 3D") {
    // 1e20 / 1e17 (legacy G5 fixture), every model on: band-gap narrowing, SRH, Auger, field
    // mobility. The probe noise puts the 1e20 side's holes 2 kT inside the band.
    const assemble::PhysicsModels all{.field_mobility = true, .fermi_dirac = true};
    const auto x = graded_axis(30);
    const auto xs = graded_axis(10);
    const auto y = uniform_axis(1e-4, 3);
    const auto [s1, x1] = probe_state(diode(*mesh::make_tensor_grid(x), 1e20, 1e17), 42, all);
    const auto [s2, x2] = probe_state(diode(*mesh::make_tensor_grid(xs, y), 1e20, 1e17), 43, all);
    const auto [s3, x3] =
        probe_state(diode(*mesh::make_tensor_grid(xs, y, y), 1e20, 1e17), 44, all);
    const double e1 = fd_jacobian_error(s1, x1);
    const double e2 = fd_jacobian_error(s2, x2);
    const double e3 = fd_jacobian_error(s3, x3);
    UNSCOPED_INFO("worst column errors: 1D " << e1 << ", 2D " << e2 << ", 3D " << e3);
    REQUIRE(s1.unknowns() >= 80);
    REQUIRE(e1 <= 5e-5);
    REQUIRE(e2 <= 5e-5);
    REQUIRE(e3 <= 5e-5);
    // And with Fermi-Dirac alone.
    const auto [s4, x4] = probe_state(diode(*mesh::make_tensor_grid(x), 1e20, 1e17), 45,
                                      {.fermi_dirac = true});
    REQUIRE(fd_jacobian_error(s4, x4) <= 5e-5);
}

TEST_CASE("drift-diffusion: the Fermi-Dirac part of the Jacobian matches finite differences") {
    // As for recombination and field mobility, checked on its own so that the column-normalized
    // gate cannot hide it: J(FD) - J(Boltzmann) against finite differences of F(FD) - F(Boltzmann),
    // per density column relative to that column's largest difference. Only the density columns
    // carry the new terms (the degeneracy factors in the driving force and the equilibrium
    // product), and only the majority ones are resolved: a minority density (1e-19 to 1e-17 here,
    // in units of Ns = 1e20) changes rows that also hold majority fluxes of order 1 by less than
    // their rounding (the Unit 13 finding). The full gate above covers those columns with its
    // absolute step. The step here is relative to the density, 1e-6.
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(20)), 1e20, 1e17);
    const auto [with, x0] = probe_state(d, 7, {.fermi_dirac = true});
    const auto [without, unused] = probe_state(d, 7);
    std::vector<double> x = x0;
    const std::size_t n = with.unknowns();
    std::vector<double> f(n), g(n);
    auto ja = with.make_jacobian(), jb = without.make_jacobian();
    with.evaluate(x, f, ja);
    without.evaluate(x, g, jb);
    const auto a = dense(ja), b = dense(jb);
    const auto difference = [&](const std::vector<double>& u, std::vector<double>& out) {
        with.residual(u, f);
        without.residual(u, g);
        for (std::size_t r = 0; r < n; ++r) out[r] = f[r] - g[r];
    };
    std::vector<double> up(n), down(n);
    double worst = 0.0, largest = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        if (c % 3 == 0 || x[c] < 1e-4) continue;  // potentials (unchanged), minority densities
        const double base = x[c];
        const double step = 1e-6 * base;
        x[c] = base + step;
        const double hi = x[c];
        difference(x, up);
        x[c] = base - step;
        const double lo = x[c];
        difference(x, down);
        x[c] = base;
        double scale = 1e-300, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(a[r][c] - b[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((up[r] - down[r]) / (hi - lo) - (a[r][c] - b[r][c])));
        }
        largest = std::max(largest, scale);
        if (scale > 1e-300) worst = std::max(worst, error / scale);
    }
    UNSCOPED_INFO("worst Fermi-Dirac column error " << worst);
    REQUIRE(largest > 0.0);
    REQUIRE(worst <= 1e-5);
}

TEST_CASE("drift-diffusion: Fermi-Dirac contacts, and the equilibrium state from a potential") {
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(15)), 1e20, 1e17);
    const auto scaling = *assemble::make_scaling(d);
    const assemble::PhysicsModels fd{.fermi_dirac = true};
    auto system = *DriftDiffusion::create(d, scaling, fd);
    const auto poisson = *assemble::EquilibriumPoisson::create(d, scaling, fd);
    // At zero bias the contacts are the Poisson system's, and the state from a potential has the
    // Poisson system's carriers. On a contact node the majority is |C| + minority (neutral to
    // rounding), which the density at the root potential reproduces to a few ulp.
    const std::vector<double> psi = poisson.charge_neutral_potential();
    const std::vector<double> x = system.state_from_potential(psi);
    std::vector<double> n(psi.size()), p(psi.size());
    poisson.carriers(psi, n, p);
    for (std::size_t i = 0; i < psi.size(); ++i) {
        CAPTURE(i);
        REQUIRE(x[3 * i] == psi[i]);
        if (poisson.is_contact()[i] != 0) {
            REQUIRE(std::abs(x[3 * i + 1] - n[i]) <= 1e-14 * std::max(n[i], p[i]));
            REQUIRE(std::abs(x[3 * i + 2] - p[i]) <= 1e-14 * std::max(n[i], p[i]));
            continue;
        }
        REQUIRE(x[3 * i + 1] == n[i]);
        REQUIRE(x[3 * i + 2] == p[i]);
    }
    // Biased: the potential moves by the bias, the densities stay.
    const std::vector<double> bias{0.4, 0.0};
    REQUIRE(system.set_bias(bias).has_value());
    std::vector<double> y = x;
    system.stamp_contacts(y);
    REQUIRE(std::abs(y[0] - (psi[0] + 0.4 / scaling.V_T)) <= 1e-14 * std::abs(y[0]));
    REQUIRE(y[1] == x[1]);
    REQUIRE(y[2] == x[2]);
}

TEST_CASE("drift-diffusion: the recombination part of the Jacobian under Fermi-Dirac statistics") {
    // As the recombination test above, with Fermi-Dirac statistics on both sides of the difference,
    // so it isolates SRH and Auger with the Fermi-Dirac equilibrium product and its partials
    // (dE/dn = E d ln gamma_n / dn). A uniform 1e20 n-type state (electrons 2.4 kT inside the band,
    // gamma_n = 0.34) at its neutral potential, with p raised 10% above equilibrium: every flux is
    // zero, n p - E is small, and E's partials are comparable to p, so they are not lost in dR/dn.
    mesh::Mesh m = *mesh::make_tensor_grid(uniform_axis(2e-4, 21));
    const std::size_t nodes = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(nodes, 0),
         .donors = std::vector<double>(nodes, 1e20),
         .acceptors = std::vector<double>(nodes, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)}}});
    const auto scaling = *assemble::make_scaling(d);
    const assemble::PhysicsModels fd{.fermi_dirac = true};
    auto with = *DriftDiffusion::create(d, scaling, fd);
    auto without = *DriftDiffusion::create(d, scaling, {.srh = false, .auger = false,
                                                        .fermi_dirac = true});
    const auto poisson = *assemble::EquilibriumPoisson::create(d, scaling, fd);
    std::vector<double> x = with.state_from_potential(poisson.charge_neutral_potential());
    const double n0 = x[3 * 2 + 1], p0 = x[3 * 2 + 2];  // a bulk node's equilibrium densities
    for (std::size_t i = 0; i < nodes; ++i) {
        x[3 * i + 1] = n0;
        x[3 * i + 2] = 1.1 * p0;
    }
    const std::size_t n = with.unknowns();
    std::vector<double> f(n), g(n), up(n), down(n);
    auto ja = with.make_jacobian(), jb = without.make_jacobian();
    with.evaluate(x, f, ja);
    without.evaluate(x, g, jb);
    const auto a = dense(ja), b = dense(jb);
    const auto difference = [&](const std::vector<double>& u, std::vector<double>& out) {
        with.residual(u, f);
        without.residual(u, g);
        for (std::size_t r = 0; r < n; ++r) out[r] = f[r] - g[r];
    };
    double worst = 0.0, largest = 0.0;
    for (std::size_t c = 3; c < n; ++c) {  // node 0 is the contact
        const double base = x[c];
        const double step = 1e-6 * std::max(std::abs(base), 1e-30);
        x[c] = base + step;
        const double hi = x[c];
        difference(x, up);
        x[c] = base - step;
        const double lo = x[c];
        difference(x, down);
        x[c] = base;
        double scale = 0.0, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(a[r][c] - b[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((up[r] - down[r]) / (hi - lo) - (a[r][c] - b[r][c])));
        }
        largest = std::max(largest, scale);
        if (scale > 0.0) worst = std::max(worst, error / scale);
    }
    UNSCOPED_INFO("worst Fermi-Dirac recombination column error " << worst);
    REQUIRE(largest > 0.0);
    REQUIRE(worst <= 1e-5);
}
