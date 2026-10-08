// Scaling, ohmic contact values and the equilibrium Poisson residual and Jacobian (ARCHITECTURE.md
// section 11, Unit 7 gate: FD-Jacobian, section 10). Scaling references computed here in
// double-double arithmetic (../physics/references.hpp).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/ohmic.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/base/constants.hpp"
#include "../physics/references.hpp"

using namespace NiTCAD;
using assemble::EquilibriumPoisson;
using assemble::Scaling;
using base::ErrorCode;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// Graded axis on [0, 2 um]: spacing shrinks geometrically (ratio 1.15) towards the junction at
// 1 um.
std::vector<double> graded_axis(int per_side) {
    std::vector<double> half{0.0};
    double h = 1e-4 * (1.15 - 1.0) / (std::pow(1.15, per_side) - 1.0);  // sums to 1 um
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

// A 1e17 / 1e17 silicon pn diode on the given mesh, anode on x_min, cathode on x_max.
device::Device diode(mesh::Mesh m, double NA = 1e17, double ND = 1e17) {
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors[i] : donors[i]) = 1.0;
    for (std::size_t i = 0; i < n; ++i) {
        acceptors[i] *= NA;
        donors[i] *= ND;
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

// Deterministic values in [-1, 1] (splitmix64), so the probe does not depend on a library RNG.
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

// The house FD-Jacobian gate (legacy tests/test_m13_solver.py::_jacobian_probe): at a perturbed
// state, central differences with step 1e-7 max(|u|, 1) for every column; per-column error is
// max |fd - J| over the column divided by the column's largest |J| (+1e-30). Returns the worst.
double fd_jacobian_error(const EquilibriumPoisson& system, std::uint64_t seed) {
    std::vector<double> psi = system.charge_neutral_potential();
    Noise noise{seed};
    for (double& v : psi) v += 0.02 * noise.next();
    const std::size_t n = system.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    linalg::SparseMatrix jacobian = system.make_jacobian();
    system.evaluate(psi, f, jacobian);
    const auto j = dense(jacobian);
    double worst = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        const double base = psi[c];
        const double step = 1e-7 * std::max(std::abs(base), 1.0);
        psi[c] = base + step;
        const double up = psi[c];
        system.residual(psi, fp);
        psi[c] = base - step;
        const double down = psi[c];
        system.residual(psi, fm);
        psi[c] = base;
        double scale = 1e-30, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((fp[r] - fm[r]) / (up - down) - j[r][c]));
        }
        worst = std::max(worst, error / scale);
    }
    return worst;
}

EquilibriumPoisson system_for(const device::Device& d) {
    return *EquilibriumPoisson::create(d, *assemble::make_scaling(d));
}

}  // namespace

TEST_CASE("scaling: the legacy definitions on a 1e17 silicon diode") {
    // V_T = kT/q, L_D = sqrt(eps V_T / (q Ns)), J0 = q D0 Ns / L_D, R0 = D0 Ns / L_D^2 and n_i,
    // in double-double arithmetic from the CODATA constants (../physics/references.hpp).
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(30)));
    const Scaling s = *assemble::make_scaling(d);
    using reference::DD;
    const DD VT = reference::thermal_voltage(300.0);
    const DD eps = DD(physics::silicon_parameters.eps_r) * DD(base::eps0_F_per_cm);
    const DD Ns = 1e17, q = base::q_C;
    const DD LD = reference::sqrt(eps * VT / (q * Ns));
    REQUIRE(s.Ns == 1e17);  // max |N_D - N_A|
    REQUIRE(close(s.V_T, VT.value(), 1e-15));
    REQUIRE(close(s.L_D, LD.value(), 1e-14));
    REQUIRE(close(s.J0, (q * Ns / LD).value(), 1e-14));
    REQUIRE(close(s.R0, (Ns / (LD * LD)).value(), 1e-14));
    REQUIRE(s.D0 == 1.0);
    REQUIRE(close(s.n_i, reference::intrinsic_density(physics::silicon_parameters, 300.0).value(),
                  1e-13));
}

TEST_CASE("scaling: Ns override, and n_i for an undoped device") {
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(10)));
    REQUIRE(assemble::make_scaling(d, 5e16)->Ns == 5e16);
    for (const double bad : {0.0, -1.0, std::nan("")}) {
        REQUIRE(assemble::make_scaling(d, bad).error().code == ErrorCode::invalid_input);
    }
    const auto intrinsic = diode(*mesh::make_tensor_grid(graded_axis(10)), 0.0, 0.0);
    const Scaling s = *assemble::make_scaling(intrinsic);
    REQUIRE(s.Ns == s.n_i);
}

TEST_CASE("ohmic: contact values are the neutral equilibrium shifted by the bias") {
    const double C = 1.0, ni = physics::intrinsic_density(physics::silicon(), 300.0) / 1e17;
    const auto zero = assemble::ohmic_contact_value(C, ni, 0.0);
    const auto e = physics::boltzmann_neutral_equilibrium(C, ni);
    REQUIRE(zero.psi == e.eta);
    REQUIRE(zero.n == e.n);
    REQUIRE(zero.p == e.p);
    REQUIRE(assemble::ohmic_contact_value(C, ni, 19.34).psi == 19.34 + e.eta);
    // p side: psi0 = -asinh(N_A / 2 n_ie).
    REQUIRE(assemble::ohmic_contact_value(-C, ni, 0.0).psi == -e.eta);
}

TEST_CASE("poisson: FD-Jacobian gate in 1D, 2D and 3D (section 10, <= 5e-5)") {
    const auto x = graded_axis(45);  // 91 nodes in 1D
    const auto y = uniform_axis(1e-4, 4);
    const auto z = uniform_axis(5e-5, 3);
    const auto xs = graded_axis(14);  // 29 x 4 x 3 = 348 nodes in 3D
    const double e1 = fd_jacobian_error(system_for(diode(*mesh::make_tensor_grid(x))), 1);
    const double e2 = fd_jacobian_error(system_for(diode(*mesh::make_tensor_grid(x, y))), 2);
    const double e3 = fd_jacobian_error(system_for(diode(*mesh::make_tensor_grid(xs, y, z))), 3);
    // Measured: 1D 1.9e-9, 2D 2.7e-9, 3D 7.4e-10.
    UNSCOPED_INFO("worst column errors: 1D " << e1 << ", 2D " << e2 << ", 3D " << e3);
    REQUIRE(e1 <= 5e-5);
    REQUIRE(e2 <= 5e-5);
    REQUIRE(e3 <= 5e-5);
    // Every column is probed; each case has at least the gate's 80 (91, 364 and 348).
    REQUIRE(x.size() == 91);
}

TEST_CASE("poisson: contact rows are Dirichlet and the neutral guess satisfies the bulk") {
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(30)));
    const EquilibriumPoisson system = system_for(d);
    const auto psi = system.charge_neutral_potential();
    std::vector<double> f(system.unknowns());
    linalg::SparseMatrix j = system.make_jacobian();
    system.evaluate(psi, f, j);
    const auto jd = dense(j);
    const std::size_t last = system.unknowns() - 1;
    for (const std::size_t c : {std::size_t{0}, last}) {
        REQUIRE(system.is_contact()[c] == 1);
        REQUIRE(f[c] == 0.0);  // the guess holds the contact value
        for (std::size_t k = 0; k < jd.size(); ++k) REQUIRE(jd[c][k] == (k == c ? 1.0 : 0.0));
    }
    // Anode (p side, N_A = Ns): psi0 = -asinh(1 / (2 n_ie / Ns)).
    const double nie = physics::intrinsic_density(physics::silicon(), 300.0) / 1e17;
    REQUIRE(close(psi[0], -std::asinh(1.0 / (2.0 * nie)), 1e-14));
    // Next to the contacts the guess is neutral and flat: residual at rounding level of the
    // charge term.
    REQUIRE(std::abs(f[1]) <= 1e-12);
    REQUIRE(std::abs(f[last - 1]) <= 1e-12);
    // At the junction it is not: the guess jumps by the built-in potential, ~32 in V_T units.
    REQUIRE(psi[last] - psi[0] > 32.0);
}

TEST_CASE("poisson: the 1D residual is the legacy row, written out") {
    // Legacy device1d.cpp:597-599 in scaled units: xs = x / L_D, h = xs[k+1] - xs[k],
    // dV = (h[i-1] + h[i]) / 2, F = (psi[i+1] - psi[i]) / h[i] - (psi[i] - psi[i-1]) / h[i-1]
    //   - dV (n - p - C), n = n_ie e^psi, p = n_ie e^-psi, et = 1.
    const auto x = graded_axis(25);
    const auto d = diode(*mesh::make_tensor_grid(x));
    const Scaling s = *assemble::make_scaling(d);
    const EquilibriumPoisson system = system_for(d);
    std::vector<double> psi = system.charge_neutral_potential();
    Noise noise{11};
    for (double& v : psi) v += 0.5 * noise.next();
    std::vector<double> f(psi.size());
    system.residual(psi, f);
    const double nie = physics::intrinsic_density(physics::silicon(), 300.0) / s.Ns;
    for (std::size_t i = 1; i + 1 < x.size(); ++i) {
        const double hl = (x[i] - x[i - 1]) / s.L_D;
        const double hr = (x[i + 1] - x[i]) / s.L_D;
        const double C = x[i] < 1e-4 ? -1.0 : 1.0;
        const double left = (psi[i] - psi[i - 1]) / hl;
        const double right = (psi[i + 1] - psi[i]) / hr;
        const double charge =
            0.5 * (hl + hr) * (nie * std::exp(psi[i]) - nie * std::exp(-psi[i]) - C);
        const double expected = right - left - charge;
        const double scale = std::max({std::abs(left), std::abs(right), std::abs(charge)});
        CAPTURE(i);
        REQUIRE(std::abs(f[i] - expected) <= 1e-13 * scale);
    }
}

TEST_CASE("poisson: a y-uniform 2D residual is the 1D residual times the scaled width") {
    const auto x = graded_axis(20);
    const auto y = std::vector<double>{0.0, 2e-5, 3e-5, 7e-5};
    const auto d1 = diode(*mesh::make_tensor_grid(x));
    const auto d2 = diode(*mesh::make_tensor_grid(x, y));
    const EquilibriumPoisson s1 = system_for(d1);
    const EquilibriumPoisson s2 = system_for(d2);
    const double L_D = assemble::make_scaling(d1)->L_D;
    std::vector<double> psi1 = s1.charge_neutral_potential();
    Noise noise{7};
    for (double& v : psi1) v += 0.3 * noise.next();
    const std::size_t nx = x.size();
    std::vector<double> psi2(nx * y.size());
    for (std::size_t k = 0; k < psi2.size(); ++k) psi2[k] = psi1[k % nx];  // x varies fastest
    std::vector<double> f1(nx), f2(psi2.size());
    s1.residual(psi1, f1);
    s2.residual(psi2, f2);
    double worst = 0.0;
    for (std::size_t k = 0; k < psi2.size(); ++k) {
        const std::size_t i = k % nx;
        if (s1.is_contact()[i] != 0) continue;
        // Transverse control width of node k, from the two meshes' volumes.
        const double width = d2.mesh().volumes()[k] / d1.mesh().volumes()[i] / L_D;
        worst = std::max(worst, std::abs(f2[k] - width * f1[i]) / (width * std::abs(f1[i])));
    }
    REQUIRE(worst <= 1e-12);
}

TEST_CASE("poisson: the Jacobian pattern is fixed, so the solver analyzes once") {
    const auto d = diode(*mesh::make_tensor_grid(graded_axis(30), uniform_axis(1e-4, 5)));
    const EquilibriumPoisson system = system_for(d);
    auto solver = *linalg::LinearSolver::create({});
    linalg::SparseMatrix j = system.make_jacobian();
    std::vector<double> psi = system.charge_neutral_potential();
    std::vector<double> f(system.unknowns()), rhs(system.unknowns()), step(system.unknowns());
    for (int iteration = 0; iteration < 3; ++iteration) {
        system.evaluate(psi, f, j);
        REQUIRE(j.has_same_pattern(system.make_jacobian()));
        REQUIRE(solver.factorize(j).has_value());
        for (std::size_t i = 0; i < f.size(); ++i) rhs[i] = -f[i];
        REQUIRE(solver.solve(rhs, step).has_value());
        for (std::size_t i = 0; i < f.size(); ++i) psi[i] += std::clamp(step[i], -5.0, 5.0);
    }
    REQUIRE(solver.analyses() == 1);
    REQUIRE(solver.factorizations() == 3);
}

TEST_CASE("poisson: a heterojunction is accepted (Unit 15), a mismatched scaling is rejected") {
    auto m = *mesh::make_tensor_grid(graded_axis(10));
    const std::size_t n = m.node_count();
    physics::SemiconductorParameters other = physics::silicon_parameters;
    other.Eg0_eV = 1.2;
    std::vector<device::RegionId> regions(n, 0);
    for (std::size_t i = n / 2; i < n; ++i) regions[i] = 1;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const auto hetero = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"a", physics::silicon()}, {"b", *physics::Semiconductor::create(other)}},
         .node_region = std::move(regions),
         .donors = std::vector<double>(n, 1e16),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    REQUIRE(EquilibriumPoisson::create(hetero, *assemble::make_scaling(hetero)).has_value());

    const auto d = diode(*mesh::make_tensor_grid(graded_axis(10)));
    Scaling s = *assemble::make_scaling(d);
    s.temperature_K = 310.0;
    REQUIRE(EquilibriumPoisson::create(d, s).error().code == ErrorCode::invalid_input);
}

// Unit 14: Fermi-Dirac statistics.

TEST_CASE("poisson: FD-Jacobian gate under Fermi-Dirac statistics on a degenerate diode") {
    // 1e20 / 1e17 (legacy G5 fixture): the 1e20 side sits 2.3 kT inside the valence band.
    const assemble::PhysicsModels fd{.fermi_dirac = true};
    const auto x = graded_axis(45);
    const auto y = uniform_axis(1e-4, 4);
    const auto z = uniform_axis(5e-5, 3);
    const auto xs = graded_axis(14);
    const auto make = [&](mesh::Mesh m) {
        const device::Device d = diode(std::move(m), 1e20, 1e17);
        return *EquilibriumPoisson::create(d, *assemble::make_scaling(d), fd);
    };
    const double e1 = fd_jacobian_error(make(*mesh::make_tensor_grid(x)), 11);
    const double e2 = fd_jacobian_error(make(*mesh::make_tensor_grid(x, y)), 12);
    const double e3 = fd_jacobian_error(make(*mesh::make_tensor_grid(xs, y, z)), 13);
    UNSCOPED_INFO("worst column errors: 1D " << e1 << ", 2D " << e2 << ", 3D " << e3);
    REQUIRE(e1 <= 5e-5);
    REQUIRE(e2 <= 5e-5);
    REQUIRE(e3 <= 5e-5);
}

TEST_CASE("poisson: Fermi-Dirac contacts and neutral guess are the neutral Fermi-Dirac root") {
    const device::Device d = diode(*mesh::make_tensor_grid(graded_axis(30)), 1e20, 1e17);
    const Scaling s = *assemble::make_scaling(d);
    const auto system = *EquilibriumPoisson::create(d, s, {.fermi_dirac = true});
    const auto boltzmann = *EquilibriumPoisson::create(d, s, {.fermi_dirac = false});
    const physics::Semiconductor si = physics::silicon();
    const double T = 300.0;
    // Band-gap narrowing is on by default: each side's n_ie is the effective one.
    for (const auto& [node, N] : {std::pair<std::size_t, double>{0, -1e20},
                                  std::pair<std::size_t, double>{60, 1e17}}) {
        CAPTURE(node);
        const double nie = physics::effective_intrinsic_density(si, std::abs(N), T);
        const auto e = physics::fermi_dirac_neutral_equilibrium(
            N / s.Ns, nie / s.Ns, std::log(physics::conduction_band_dos(si, T) / nie),
            std::log(physics::valence_band_dos(si, T) / nie));
        REQUIRE(system.contact_potential()[node] == e.eta);
        REQUIRE(system.charge_neutral_potential()[node] == e.eta);
        // The degenerate contact needs a deeper potential than Boltzmann statistics give.
        REQUIRE(std::abs(system.contact_potential()[node]) >=
                std::abs(boltzmann.contact_potential()[node]));
    }
    REQUIRE(system.contact_potential()[0] - boltzmann.contact_potential()[0] < -0.5);
    // In the neutral bulk the guess leaves no charge: the residual is the (zero) flux balance.
    const std::vector<double> psi = system.charge_neutral_potential();
    std::vector<double> f(system.unknowns());
    system.residual(psi, f);
    std::vector<double> n(system.unknowns()), p(system.unknowns());
    system.carriers(psi, n, p);
    REQUIRE(close(n[3] - p[3], 0.0 - 1e20 / s.Ns, 1e-14));
    REQUIRE(close(n[58] - p[58], 1e17 / s.Ns, 1e-13));
}
