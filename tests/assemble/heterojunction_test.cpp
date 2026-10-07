// Heterointerfaces in the assemblers (ARCHITECTURE.md section 11, Unit 15): the band shift, the
// permittivity step, thermionic-emission fluxes and their Jacobians.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/thermionic_flux.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using assemble::EquilibriumPoisson;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

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

constexpr auto drift_diffusion = device::InterfaceTransport::drift_diffusion;
constexpr auto thermionic = device::InterfaceTransport::thermionic_emission;

// Material `left` on x < 1 um with net doping C_left, `right` beyond with C_right (positive:
// donors); ohmic contacts on x_min and x_max; the interface declared with `transport`.
device::Device junction(mesh::Mesh m, const physics::SemiconductorParameters& left,
                        const physics::SemiconductorParameters& right, double C_left,
                        double C_right, device::InterfaceTransport transport = drift_diffusion) {
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const bool is_left = m.points()[i][0] < 1e-4;
        region[i] = is_left ? 0 : 1;
        const double C = is_left ? C_left : C_right;
        (C > 0.0 ? donors[i] : acceptors[i]) = std::abs(C);
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"left", *physics::Semiconductor::create(left)},
                     {"right", *physics::Semiconductor::create(right)}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .interfaces = {{"left", "right", transport}}});
}

physics::SemiconductorParameters silicon_with(double d_chi, double d_eps = 0.0) {
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.electron_affinity_eV += d_chi;
    p.eps_r += d_eps;
    return p;
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

// The legacy probe state (as drift_diffusion_test.cpp): the system's own equilibrium guess with the
// anode at +0.3 V, psi + 0.02 noise, n and p times (1 + 0.01 noise).
std::pair<DriftDiffusion, std::vector<double>> probe_state(
    const device::Device& d, std::uint64_t seed, const assemble::PhysicsModels& models = {}) {
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling, models);
    const std::vector<double> bias{0.3, 0.0};
    REQUIRE(system.set_bias(bias).has_value());
    const auto poisson = *EquilibriumPoisson::create(d, scaling, models);
    auto x = system.state_from_potential(poisson.charge_neutral_potential());
    Noise noise{seed};
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[3 * i] += 0.02 * noise.next();
        x[3 * i + 1] *= 1.0 + 0.01 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.01 * noise.next();
    }
    return {std::move(system), std::move(x)};
}

// The house gate (legacy _jacobian_probe), every column.
template <class System>
double fd_jacobian_error(const System& system, std::vector<double> x) {
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

TEST_CASE("thermionic flux: zero at equilibrium for either step, exact partials elsewhere") {
    // n2 / n1 = e^delta is equilibrium whatever the step: n2 g2 = n1 g1.
    for (const double log_ratio : {-2.3, 0.0, 1.7}) {
        for (const double delta : {-12.0, -1.0, 0.4, 9.0}) {
            CAPTURE(log_ratio, delta);
            const double ratio = std::exp(-log_ratio);
            const double n1 = 3.0, n2 = n1 * std::exp(delta);
            const auto fn =
                assemble::thermionic_electron_flux(2.5, delta, log_ratio, ratio, n1, n2);
            REQUIRE(std::abs(fn.flux) <= 1e-14 * 2.5 * (n1 + n2));
            const double p1 = 3.0, p2 = p1 * std::exp(-delta);
            const auto fp = assemble::thermionic_hole_flux(2.5, delta, log_ratio, ratio, p1, p2);
            REQUIRE(std::abs(fp.flux) <= 1e-14 * 2.5 * (p1 + p2));
        }
    }
    // Partials against central differences, away from the kink.
    const double K = 1.7, lr = 0.8, r = std::exp(-0.3);
    for (const double delta : {-3.0, -0.2, 1.5, 4.0}) {
        CAPTURE(delta);
        const double n1 = 2.0, n2 = 0.7, h = 1e-6;
        const auto f = assemble::thermionic_electron_flux(K, delta, lr, r, n1, n2);
        const auto at = [&](double d, double a, double b) {
            return assemble::thermionic_electron_flux(K, d, lr, r, a, b).flux;
        };
        REQUIRE(close(f.d_psi2, (at(delta + h, n1, n2) - at(delta - h, n1, n2)) / (2 * h), 1e-8));
        REQUIRE(f.d_psi1 == -f.d_psi2);
        REQUIRE(close(f.d_c1, (at(delta, n1 + h, n2) - at(delta, n1 - h, n2)) / (2 * h), 1e-8));
        REQUIRE(close(f.d_c2, (at(delta, n1, n2 + h) - at(delta, n1, n2 - h)) / (2 * h), 1e-8));
        const auto g = assemble::thermionic_hole_flux(K, delta, lr, r, n1, n2);
        const auto hat = [&](double d, double a, double b) {
            return assemble::thermionic_hole_flux(K, d, lr, r, a, b).flux;
        };
        REQUIRE(close(g.d_psi2, (hat(delta + h, n1, n2) - hat(delta - h, n1, n2)) / (2 * h), 1e-8));
        REQUIRE(close(g.d_c1, (hat(delta, n1 + h, n2) - hat(delta, n1 - h, n2)) / (2 * h), 1e-8));
        REQUIRE(close(g.d_c2, (hat(delta, n1, n2 + h) - hat(delta, n1, n2 - h)) / (2 * h), 1e-8));
    }
    // Emission over a barrier is cut by e^-barrier; down a step it is not.
    const auto up = assemble::thermionic_electron_flux(1.0, -5.0, 0.0, 1.0, 1.0, 0.0);
    const auto down = assemble::thermionic_electron_flux(1.0, 5.0, 0.0, 1.0, 1.0, 0.0);
    REQUIRE(close(-up.flux, std::exp(-5.0), 1e-15));
    REQUIRE(-down.flux == 1.0);
}

TEST_CASE("heterojunction: the 1D Poisson row across a permittivity and affinity step") {
    // Silicon | silicon with chi + 0.2 eV and eps_r + 1.3, n-type 1e17: the interface edge carries
    // the harmonic mean of the two permittivity ratios, the nodes their own charge with the band
    // shift s = 0.2 / V_T on the right.
    const auto x = graded_axis(8);
    const device::Device d =
        junction(*mesh::make_tensor_grid(x), physics::silicon_parameters, silicon_with(0.2, 1.3),
                 1e17, 1e17);
    const auto s = *assemble::make_scaling(d);
    const auto system = *EquilibriumPoisson::create(d, s);
    const std::size_t n = x.size();
    std::vector<double> psi(n);
    Noise noise{3};
    for (double& v : psi) v = 0.3 * noise.next();
    std::vector<double> f(n);
    system.residual(psi, f);
    const double nie = physics::intrinsic_density(physics::silicon(), 300.0) / s.Ns;
    const double et_r = (11.7 + 1.3) / 11.7;
    const double et_iface = 2.0 * et_r / (1.0 + et_r);
    const double shift = 0.2 / s.V_T;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        CAPTURE(i);
        const bool right = x[i] >= 1e-4;
        const double hl = (x[i] - x[i - 1]) / s.L_D, hr = (x[i + 1] - x[i]) / s.L_D;
        const auto et_of = [&](std::size_t a, std::size_t b) {
            const bool ra = x[a] >= 1e-4, rb = x[b] >= 1e-4;
            return ra != rb ? et_iface : (ra ? et_r : 1.0);
        };
        const double eta = psi[i] + (right ? shift : 0.0);
        const double rho = nie * std::exp(eta) - nie * std::exp(-eta) - 1e17 / s.Ns;
        const double expected = et_of(i, i + 1) * (psi[i + 1] - psi[i]) / hr -
                                et_of(i - 1, i) * (psi[i] - psi[i - 1]) / hl -
                                0.5 * (hl + hr) * rho;
        REQUIRE(std::abs(f[i] - expected) <= 1e-12 * (std::abs(expected) + 1.0));
    }
}

TEST_CASE("heterojunction: contacts and the neutral guess carry the band shift") {
    const auto x = graded_axis(10);
    const device::Device d = junction(*mesh::make_tensor_grid(x), physics::silicon_parameters,
                                      silicon_with(-0.25), -1e17, 1e17);
    const auto s = *assemble::make_scaling(d);
    const auto system = *EquilibriumPoisson::create(d, s);
    const auto homo = *EquilibriumPoisson::create(
        junction(*mesh::make_tensor_grid(x), physics::silicon_parameters,
                 physics::silicon_parameters, -1e17, 1e17),
        s);
    // The right side's intrinsic level is 0.25 eV shallower: its potential at the same carrier
    // density is 0.25 V higher. The left (node 0's material) is the reference.
    const std::size_t last = x.size() - 1;
    REQUIRE(system.contact_potential()[0] == homo.contact_potential()[0]);
    REQUIRE(close(system.contact_potential()[last] - homo.contact_potential()[last],
                  0.25 / s.V_T, 1e-12));
    const auto guess = system.charge_neutral_potential();
    REQUIRE(close(guess[last - 3] - homo.charge_neutral_potential()[last - 3], 0.25 / s.V_T,
                  1e-12));
    // The carriers at the guess are those of the homojunction: the shift moves psi, not n or p.
    std::vector<double> n1(x.size()), p1(x.size()), n2(x.size()), p2(x.size());
    system.carriers(guess, n1, p1);
    homo.carriers(homo.charge_neutral_potential(), n2, p2);
    for (std::size_t i = 0; i < x.size(); ++i) {
        CAPTURE(i);
        REQUIRE(close(n1[i], n2[i], 1e-11));
        REQUIRE(close(p1[i], p2[i], 1e-11));
    }
}

TEST_CASE("heterojunction: FD-Jacobian gates in 1D, 2D and 3D, with thermionic emission and FD") {
    // Si | Si with chi + 0.3 eV, eps_r + 1.3 and Eg0 + 0.1 eV, a p-n junction at the interface.
    physics::SemiconductorParameters right = silicon_with(0.3, 1.3);
    right.Eg0_eV += 0.1;
    const auto x = graded_axis(30);
    const auto xs = graded_axis(10);
    const auto y = uniform_axis(1e-4, 3);
    struct Case {
        assemble::PhysicsModels models;
        device::InterfaceTransport transport;
    };
    for (const Case& c : {Case{{}, drift_diffusion}, Case{{}, thermionic},
                          Case{{.fermi_dirac = true}, thermionic},
                          Case{{.field_mobility = true, .incomplete_ionization = true},
                               drift_diffusion}}) {
        const assemble::PhysicsModels& models = c.models;
        const auto make = [&](mesh::Mesh m) {
            return junction(std::move(m), physics::silicon_parameters, right, -1e17, 1e17,
                            c.transport);
        };
        CAPTURE(c.transport == thermionic, models.fermi_dirac, models.field_mobility);
        const auto [s1, x1] = probe_state(make(*mesh::make_tensor_grid(x)), 42, models);
        const auto [s2, x2] = probe_state(make(*mesh::make_tensor_grid(xs, y)), 43, models);
        const auto [s3, x3] = probe_state(make(*mesh::make_tensor_grid(xs, y, y)), 44, models);
        const double e1 = fd_jacobian_error(s1, x1);
        const double e2 = fd_jacobian_error(s2, x2);
        const double e3 = fd_jacobian_error(s3, x3);
        UNSCOPED_INFO("drift-diffusion worst column errors: 1D " << e1 << ", 2D " << e2 << ", 3D "
                                                                 << e3);
        REQUIRE(e1 <= 5e-5);
        REQUIRE(e2 <= 5e-5);
        REQUIRE(e3 <= 5e-5);
    }
    const device::Device d =
        junction(*mesh::make_tensor_grid(x), physics::silicon_parameters, right, -1e17, 1e17);
    const auto poisson = *EquilibriumPoisson::create(d, *assemble::make_scaling(d));
    std::vector<double> psi = poisson.charge_neutral_potential();
    Noise noise{9};
    for (double& v : psi) v += 0.02 * noise.next();
    const double e = fd_jacobian_error(poisson, psi);
    UNSCOPED_INFO("Poisson worst column error " << e);
    REQUIRE(e <= 5e-5);
}

TEST_CASE("heterojunction: the thermionic part of the Jacobian matches finite differences") {
    // J(TE) - J(SG) against finite differences of F(TE) - F(SG): only the interface edge differs,
    // so the difference isolates the thermionic fluxes' partials. Potential and majority-density
    // columns (minority densities are below what a relative step resolves, Unit 14).
    const device::Device d = junction(*mesh::make_tensor_grid(graded_axis(20)),
                                      physics::silicon_parameters, silicon_with(0.3), 1e17, 1e17);
    const device::Device dt = junction(*mesh::make_tensor_grid(graded_axis(20)),
                                       physics::silicon_parameters, silicon_with(0.3), 1e17, 1e17,
                                       thermionic);
    const auto [with, x0] = probe_state(dt, 5);
    const auto [without, unused] = probe_state(d, 5);
    std::vector<double> x = x0;
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
        if (c % 3 != 0 && x[c] < 1e-4) continue;
        const double base = x[c];
        const double step = c % 3 == 0 ? 1e-7 * std::max(std::abs(base), 1.0) : 1e-6 * base;
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
    UNSCOPED_INFO("worst thermionic column error " << worst);
    REQUIRE(largest > 0.0);
    REQUIRE(worst <= 1e-5);
}

TEST_CASE("heterojunction: the thermionic edge factor is hmean(v) times the scaled area") {
    // An n-n junction where electrons cross a 0.3 eV step: at a state with psi flat and n2 = 0, the
    // electron current through the interface edge is K n1 e^-u... evaluated through the public
    // edge currents: Jn = -K n1 g1 with g1 = e^(delta - ln(Nc2/Nc1)), same Nc, delta = s2 - s1.
    const auto x = graded_axis(10);
    const device::Device d = junction(*mesh::make_tensor_grid(x), physics::silicon_parameters,
                                      silicon_with(-0.3), 1e17, 1e17);
    const device::Device dt = junction(*mesh::make_tensor_grid(x), physics::silicon_parameters,
                                       silicon_with(-0.3), 1e17, 1e17, thermionic);
    const auto s = *assemble::make_scaling(d);
    const auto system = *DriftDiffusion::create(dt, s);
    std::vector<double> state(system.unknowns());
    for (std::size_t i = 0; i < x.size(); ++i) {
        state[3 * i] = 0.0;
        state[3 * i + 1] = 1.0;
        state[3 * i + 2] = 1e-6;
    }
    const auto currents = system.edge_currents(state);
    std::size_t iface = 0;
    while (!(x[iface] < 1e-4 && x[iface + 1] >= 1e-4)) ++iface;
    const double v = physics::emission_velocity_cm_s(
        physics::conduction_band_dos(physics::silicon(), 300.0), 300.0);
    const double K = v * s.L_D / s.D0;  // 1D: unit area
    const double delta = -0.3 / s.V_T;  // s_right - s_left
    // n1 = n2 = 1: Jn = K (g2 - g1), u = delta < 0: g1 = e^u, g2 = 1.
    REQUIRE(close(currents[iface].first, K * (1.0 - std::exp(delta)), 1e-13));
    // A homojunction edge is untouched.
    const auto plain = *DriftDiffusion::create(d, s);
    REQUIRE(plain.edge_currents(state)[iface + 3] == currents[iface + 3]);
}

TEST_CASE("heterojunction: field mobility across a material step is the mean of each side's") {
    // Silicon | silicon with chi + 0.1 and its own Canali parameters (v_sat 0.5x, beta 1): the
    // interface edge's mobility is hmean(mu_C,a(E), mu_C,b(E)), each end with its own material's
    // low-field mobility and saturation; a homojunction edge keeps Unit 13's mu_C(hmean).
    physics::SemiconductorParameters right = silicon_with(0.1);
    right.electron_saturation = {.v_sat_cm_s = 0.5e7, .beta = 1.0};
    right.electron_mobility.mu_max = 1000.0;
    const auto x = graded_axis(10);
    const device::Device d =
        junction(*mesh::make_tensor_grid(x), physics::silicon_parameters, right, 1e17, 1e17);
    const auto s = *assemble::make_scaling(d);
    const auto with = *DriftDiffusion::create(d, s, {.field_mobility = true});
    const auto without = *DriftDiffusion::create(d, s);
    std::vector<double> state(with.unknowns());
    Noise noise{17};
    for (std::size_t i = 0; i < x.size(); ++i) {
        state[3 * i] = 2.0 * noise.next();  // edge fields of order 1e4-1e6 V/cm
        state[3 * i + 1] = 1.0 + 0.1 * noise.next();
        state[3 * i + 2] = 1e-6;
    }
    const auto a = with.edge_currents(state), b = without.edge_currents(state);
    std::size_t iface = 0;
    while (!(x[iface] < 1e-4 && x[iface + 1] >= 1e-4)) ++iface;
    const physics::Semiconductor si = physics::silicon();
    const physics::Semiconductor r = *physics::Semiconductor::create(right);
    const auto electron = physics::Carrier::electron;
    const double mu_a = physics::caughey_thomas_mobility(si, electron, 1e17, 300.0);
    const double mu_b = physics::caughey_thomas_mobility(r, electron, 1e17, 300.0);
    const double E = s.V_T * std::abs(state[3 * (iface + 1)] - state[3 * iface]) /
                     (x[iface + 1] - x[iface]);
    const double ma =
        physics::canali_mobility(mu_a, E, si.parameters().electron_saturation).mobility;
    const double mb = physics::canali_mobility(mu_b, E, right.electron_saturation).mobility;
    const double factor = (2.0 * ma * mb / (ma + mb)) / (2.0 * mu_a * mu_b / (mu_a + mu_b));
    CAPTURE(E, factor);
    REQUIRE(factor < 0.9);  // the field is high enough to matter
    REQUIRE(close(a[iface].first, factor * b[iface].first, 1e-13));
    // And the Jacobian stays exact there.
    const auto [system, probe] = probe_state(d, 21, {.field_mobility = true});
    REQUIRE(fd_jacobian_error(system, probe) <= 5e-5);
}

// Unit 15 fixes: incomplete ionization, radiative recombination, band edges, a non-planar
// interface.

namespace {

// J(with) - J(without) against finite differences of F(with) - F(without), per column relative to
// that column's largest difference, on the columns `use` selects; step relative (1e-6) on
// densities, absolute (1e-7) on potentials.
template <class Use>
double part_error(const DriftDiffusion& with, const DriftDiffusion& without, std::vector<double> x,
                  Use use) {
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
        if (!use(c, x[c])) continue;
        const double base = x[c];
        const double step = c % 3 == 0 ? 1e-7 * std::max(std::abs(base), 1.0) : 1e-6 * base;
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
    REQUIRE(largest > 0.0);
    return worst;
}

}  // namespace

TEST_CASE("heterojunction: FD-Jacobian gates with incomplete ionization (Poisson and DD)") {
    // A 4H-SiC p-n junction (aluminium 220 meV: about 10% ionized) and a silicon one, Boltzmann
    // and Fermi-Dirac.
    const physics::SemiconductorParameters sic = physics::silicon_carbide_4h_parameters;
    for (const bool fd : {false, true}) {
        CAPTURE(fd);
        const assemble::PhysicsModels models{.fermi_dirac = fd, .incomplete_ionization = true};
        const device::Device d =
            junction(*mesh::make_tensor_grid(graded_axis(30)), sic, sic, -1e17, 1e17);
        const auto poisson = *EquilibriumPoisson::create(d, *assemble::make_scaling(d), models);
        std::vector<double> psi = poisson.charge_neutral_potential();
        Noise noise{31};
        for (double& v : psi) v += 0.02 * noise.next();
        REQUIRE(fd_jacobian_error(poisson, psi) <= 5e-5);
        const auto [system, x] = probe_state(d, 32, models);
        REQUIRE(fd_jacobian_error(system, x) <= 5e-5);
        const auto [s2, x2] = probe_state(junction(*mesh::make_tensor_grid(graded_axis(30)),
                                                   physics::silicon_parameters,
                                                   physics::silicon_parameters, -1e17, 1e17),
                                          33, models);
        REQUIRE(fd_jacobian_error(s2, x2) <= 5e-5);
    }
}

TEST_CASE("heterojunction: the ionization part of the Jacobian matches finite differences") {
    // The ionized doping depends on the majority densities (its own band's): their columns.
    const physics::SemiconductorParameters sic = physics::silicon_carbide_4h_parameters;
    const device::Device d =
        junction(*mesh::make_tensor_grid(graded_axis(20)), sic, sic, -1e17, 1e17);
    for (const bool fd : {false, true}) {
        CAPTURE(fd);
        const auto [with, x0] =
            probe_state(d, 41, {.fermi_dirac = fd, .incomplete_ionization = true});
        const auto [without, unused] = probe_state(d, 41, {.fermi_dirac = fd});
        const double e = part_error(with, without, x0, [](std::size_t c, double v) {
            return c % 3 != 0 && v > 1e-4;
        });
        UNSCOPED_INFO("worst ionization column error " << e);
        REQUIRE(e <= 1e-5);
    }
}

TEST_CASE("heterojunction: the radiative part of the Jacobian matches finite differences") {
    // A uniform GaAs state out of equilibrium (every flux zero): J(radiative) - J(none).
    mesh::Mesh m = *mesh::make_tensor_grid(uniform_axis(2e-4, 21));
    const std::size_t n_nodes = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"gaas", *physics::Semiconductor::create(
                                   physics::gallium_arsenide_parameters)}},
         .node_region = std::vector<device::RegionId>(n_nodes, 0),
         .donors = std::vector<double>(n_nodes, 1e17),
         .acceptors = std::vector<double>(n_nodes, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)}}});
    const auto s = *assemble::make_scaling(d);
    const auto with = *DriftDiffusion::create(d, s, {.srh = false, .auger = false});
    const auto without =
        *DriftDiffusion::create(d, s, {.srh = false, .auger = false, .radiative = false});
    std::vector<double> x(with.unknowns());
    for (std::size_t i = 0; i < n_nodes; ++i) {
        x[3 * i] = 0.1;
        x[3 * i + 1] = 1.0;
        x[3 * i + 2] = 1e-3;
    }
    const double e = part_error(with, without, x, [](std::size_t c, double) { return c >= 3; });
    UNSCOPED_INFO("worst radiative column error " << e);
    REQUIRE(e <= 1e-8);
    // The rate is B (n p - n_ie^2): the hole row of a bulk node, divided by its volume.
    std::vector<double> f(with.unknowns()), g(with.unknowns());
    with.residual(x, f);
    without.residual(x, g);
    const double nie = physics::intrinsic_density(*physics::Semiconductor::create(
                                                       physics::gallium_arsenide_parameters),
                                                   300.0);
    const double R = 7.2e-10 * (s.Ns * s.Ns * 1e-3 - nie * nie);
    const double V = (1e-5) / s.L_D;  // the uniform 1D control volume, scaled
    REQUIRE(close((f[3 * 5 + 2] - g[3 * 5 + 2]) / V, R / s.R0, 1e-12));
}

TEST_CASE("heterojunction: band edges are flat-Fermi at equilibrium and independent of node 0") {
    // GaAs | Al0.3Ga0.7As p-n and its mirror image (node 0 then in AlGaAs): at the neutral guess
    // the band edges of mirrored nodes agree, although the potentials differ by the band shift.
    // The axis straddles the interface (no node on it), so the mirror is exact.
    const auto graded = graded_axis(20);
    const auto x = *mesh::straddle_interface(graded, 1e-4, graded[21] - graded[20]);
    const physics::SemiconductorParameters gaas = physics::gallium_arsenide_parameters;
    const physics::SemiconductorParameters algaas = *physics::algaas_parameters(0.3);
    const device::Device d = junction(*mesh::make_tensor_grid(x), gaas, algaas, -1e17, 1e17);
    const device::Device m = junction(*mesh::make_tensor_grid(x), algaas, gaas, 1e17, -1e17);
    const auto pd = *EquilibriumPoisson::create(d, *assemble::make_scaling(d));
    const auto pm = *EquilibriumPoisson::create(m, *assemble::make_scaling(m));
    const auto bd = pd.band_edges(pd.charge_neutral_potential());
    const auto bm = pm.band_edges(pm.charge_neutral_potential());
    const std::size_t n = x.size();
    for (std::size_t i = 0; i < n; ++i) {
        CAPTURE(i);
        REQUIRE(std::abs(bd.electron_fermi[i]) <= 1e-12);
        REQUIRE(std::abs(bd.hole_fermi[i]) <= 1e-12);
        REQUIRE(std::abs(bd.conduction[i] - bm.conduction[n - 1 - i]) <= 1e-12);
        REQUIRE(std::abs(bd.valence[i] - bm.valence[n - 1 - i]) <= 1e-12);
    }
    REQUIRE(std::abs(pd.charge_neutral_potential()[0] - pm.charge_neutral_potential()[n - 1]) >
            1.0);  // the potentials themselves do differ
}

TEST_CASE("heterojunction: a non-planar interface in 2D (FD-Jacobian, thermionic)") {
    // GaAs on x < 1 um for y < 0.5 um and on x < 1.2 um above: the interface has a step, so its
    // edges run in both directions. Declared thermionic.
    const auto x = graded_axis(10);
    const auto y = uniform_axis(1e-4, 5);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& pt = m.points()[i];
        const double edge = pt[1] < 0.5e-4 ? 1e-4 : 1.2e-4;
        region[i] = pt[0] < edge ? 0 : 1;
        (pt[0] < 1e-4 ? acceptors[i] : donors[i]) = 1e17;
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"gaas",
                      *physics::Semiconductor::create(physics::gallium_arsenide_parameters)},
                     {"algaas", *physics::Semiconductor::create(*physics::algaas_parameters(0.3))}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .interfaces = {{"gaas", "algaas", thermionic}}});
    std::size_t vertical = 0;
    for (const mesh::Edge& e : d.mesh().edges()) {
        const bool cross = d.node_region()[static_cast<std::size_t>(e.first)] !=
                           d.node_region()[static_cast<std::size_t>(e.second)];
        if (cross && std::abs(d.mesh().points()[static_cast<std::size_t>(e.first)][0] -
                              d.mesh().points()[static_cast<std::size_t>(e.second)][0]) < 1e-12) {
            ++vertical;
        }
    }
    REQUIRE(vertical > 0);  // interface edges along y exist
    const auto [system, probe] = probe_state(d, 51, {.fermi_dirac = true});
    const double e = fd_jacobian_error(system, probe);
    UNSCOPED_INFO("worst column error " << e);
    REQUIRE(e <= 5e-5);
}
