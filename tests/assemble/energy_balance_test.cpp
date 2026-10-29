// Energy-balance transport in the drift-diffusion assembler (Unit 52; DECISIONS.md ED1-ED14): the
// finite-difference Jacobian gate of the five-unknown system (EC6: Boltzmann, Fermi-Dirac, impact
// ionization, Lombardi) in 1D, 2D and 3D, across a material step and into an oxide; a hot
// degenerate carrier at constant pressure (ED8 as ruled); at tau = 0 the continuity rows and
// currents of drift-diffusion; equilibrium tau = 0 exactly (EC5's assembly half); the uniform bulk
// state of Canali's velocity and the field-driven generation at 1e5 V/cm (EC3, as ruled); and the
// inputs the model refuses (ED11, EC15).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/fermi_dirac.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/insulator.hpp"
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

// Silicon with its bands, densities of states, mobility, saturation and relaxation times changed.
physics::Semiconductor other_material() {
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.electron_affinity_eV += 0.08;
    p.Eg0_eV += 0.05;
    p.Nc300 *= 1.6;
    p.Nv300 *= 0.7;
    p.electron_mobility.mu_max *= 0.8;
    p.hole_mobility.mu_max *= 1.1;
    p.electron_saturation = {.v_sat_cm_s = 0.9e7, .beta = 1.7};
    p.hole_saturation = {.v_sat_cm_s = 0.8e7, .beta = 1.0};
    p.energy_relaxation = {.tau_n_s = 0.3e-12, .tau_p_s = 0.25e-12};
    return *physics::Semiconductor::create(p);
}

struct Options {
    int D = 1;
    int nodes = 11;
    bool hetero = false;  // x >= 0.6 um in other_material()
    bool oxide = false;   // 2D/3D: y >= 0.3 um oxide (contacts on the silicon rows only)
    bool gate = false;    // 2D: a lumped gate on y_max
    bool schottky = false;
    bool thermionic = false;
    bool interface = false;
};

assemble::PhysicsModels models(bool energy = true) { return {.energy_balance = energy}; }

std::vector<double> biases(const Options& o) {
    std::vector<double> v{0.4, -0.3};
    if (o.gate) v.push_back(0.7);
    return v;
}

// A pn diode on [0, 1 um] along x (p 1e17 below 0.5 um, n 2e17 above), [0, 0.5 um] across.
device::Device diode(const Options& o) {
    const auto x = uniform(0.0, 1e-4, o.nodes);
    const auto y = uniform(0.0, 0.5e-4, o.D == 3 ? 4 : 6);
    mesh::Mesh m = o.D == 1   ? *mesh::make_tensor_grid(x)
                   : o.D == 2 ? *mesh::make_tensor_grid(x, y)
                              : *mesh::make_tensor_grid(x, y, y);
    const std::size_t n = m.node_count();
    const device::RegionId other = 1, oxide = o.hetero ? 2 : 1;
    std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    if (o.hetero) regions.push_back({"other", other_material()});
    if (o.oxide) regions.push_back({"oxide", physics::silicon_dioxide()});
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if (o.oxide && p[1] >= 0.3e-4 - 1e-12) {
            region[i] = oxide;
            continue;
        }
        if (o.hetero && p[0] >= 0.6e-4 - 1e-12) region[i] = other;
        (p[0] < 0.5e-4 ? acceptors : donors)[i] = p[0] < 0.5e-4 ? 1e17 : 2e17;
    }
    const auto silicon_nodes = [&](const char* patch) {
        std::vector<mesh::NodeId> v;
        for (const mesh::NodeId k : m.find_boundary(patch)->nodes) {
            if (!o.oxide || region[static_cast<std::size_t>(k)] != oxide) v.push_back(k);
        }
        return v;
    };
    std::vector<device::Contact> contacts{
        {"anode", device::ContactKind::ohmic, silicon_nodes("x_min")},
        {"cathode", device::ContactKind::ohmic, silicon_nodes("x_max")}};
    if (o.schottky) {
        contacts[1].kind = device::ContactKind::schottky;
        contacts[1].schottky = {.boundary = "x_max", .by_work_function = true, .work_function_eV = 4.6};
    }
    if (o.gate) {
        std::vector<mesh::NodeId> top;
        for (const mesh::NodeId k : m.find_boundary("y_max")->nodes) {
            const double xk = m.points()[static_cast<std::size_t>(k)][0];
            if (xk > 1e-12 && xk < 1e-4 - 1e-12) top.push_back(k);
        }
        contacts.push_back({"gate", device::ContactKind::gate, std::move(top),
                            {.boundary = "y_max", .oxide_thickness_cm = 5e-7,
                             .electrode = device::GateElectrode::n_poly}});
    }
    std::vector<device::Interface> interfaces;
    if (o.interface) {
        interfaces.push_back({.region_a = "silicon", .region_b = "oxide",
                              .fixed_charge_cm2 = 1e11});
    }
    if (o.thermionic) {
        interfaces.push_back(
            {"silicon", "other", device::InterfaceTransport::thermionic_emission});
    }
    auto device = device::Device::create({.mesh = std::move(m),
                                          .temperature_K = 300.0,
                                          .regions = std::move(regions),
                                          .node_region = std::move(region),
                                          .donors = std::move(donors),
                                          .acceptors = std::move(acceptors),
                                          .contacts = std::move(contacts),
                                          .interfaces = interfaces});
    if (!device) FAIL(device.error().message);
    return std::move(*device);
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

// A biased, heated state: psi falling 20 V_T across the device, n 1e16 and p 1e15 cm^-3, theta_n
// rising from 1 to 2.5 and theta_p to 1.6 along x, every unknown with relative noise. Insulator
// densities and rises 0.
std::vector<double> heated(const DriftDiffusion& dd, const device::Device& d,
                           const assemble::Scaling& s, double noise, std::uint64_t seed) {
    const std::size_t m = dd.stride();
    std::vector<double> x(dd.unknowns());
    Noise z{seed};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        const double u = d.mesh().points()[i][0] / 1e-4;
        x[m * i] = -20.0 * u * (1.0 + noise * z.next());
        const bool insulator = d.is_insulator(static_cast<mesh::NodeId>(i));
        x[m * i + 1] = insulator ? 0.0 : 1e16 / s.Ns * (1.0 + noise * z.next());
        x[m * i + 2] = insulator ? 0.0 : 1e15 / s.Ns * (1.0 + noise * z.next());
        if (m == 5 && !insulator) {
            x[m * i + 3] = (1.0 + 1.5 * u) * (1.0 + noise * z.next()) - 1.0;
            x[m * i + 4] = (1.0 + 0.6 * u) * (1.0 + noise * z.next()) - 1.0;
        }
    }
    return x;
}

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

// The finite-difference gate (as electrothermal_test.cpp): every column, central differences,
// the worst error over the column's largest entry, each column's the smallest at steps of 1e-5,
// 1e-6 and 1e-7 max(|x|, 1). With `relative`, steps of |x| itself (densities far below 1, which a
// step of max(|x|, 1) would drive negative); a row then cannot resolve a change below the rounding
// of its own terms, so each row's error beyond 64 eps S_r / (2 delta) counts, S_r the row's term
// size sum_c |J_rc x_c| + |f_r|.
double fd_jacobian_error(const DriftDiffusion& dd, std::vector<double> x, bool relative = false) {
    const std::size_t n = dd.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    auto jacobian = dd.make_jacobian();
    dd.evaluate(x, f, jacobian);
    const auto j = dense(jacobian);
    std::vector<double> terms(n, 0.0);
    for (std::size_t r = 0; r < n; ++r) {
        terms[r] = std::abs(f[r]);
        for (std::size_t c = 0; c < n; ++c) terms[r] += std::abs(j[r][c] * x[c]);
    }
    double worst = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        double scale = 1e-300;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][c]));
        double best = std::numeric_limits<double>::infinity();
        for (const double step : {1e-5, 1e-6, 1e-7}) {
            const double base = x[c];
            const double delta =
                step * (relative && base != 0.0 ? std::abs(base) : std::max(std::abs(base), 1.0));
            x[c] = base + delta;
            const double up = x[c];
            dd.residual(x, fp);
            x[c] = base - delta;
            const double down = x[c];
            dd.residual(x, fm);
            x[c] = base;
            double error = 0.0;
            for (std::size_t r = 0; r < n; ++r) {
                const double noise = relative ? 64.0 * std::numeric_limits<double>::epsilon() *
                                                    terms[r] / (up - down)
                                              : 0.0;
                error = std::max(error,
                                 std::abs((fp[r] - fm[r]) / (up - down) - j[r][c]) - noise);
            }
            best = std::min(best, error / scale);
        }
        worst = std::max(worst, best);
    }
    return worst;
}

DriftDiffusion make(const device::Device& d, const assemble::Scaling& s, const Options& o,
                    const assemble::PhysicsModels& m) {
    auto dd = DriftDiffusion::create(d, s, m);
    if (!dd) FAIL(dd.error().message);
    REQUIRE(dd->set_bias(biases(o)).has_value());
    return std::move(*dd);
}

DriftDiffusion make(const device::Device& d, const assemble::Scaling& s, const Options& o,
                    bool energy = true) {
    return make(d, s, o, models(energy));
}

// The heated state with degenerate and dilute carriers (slice 3, Fermi-Dirac): n from 1e16 to
// 1e20 cm^-3 and p from 1e19 down to 1e10 along x (both sides of N_c and of the 1e-6 N series
// switch), the rises as heated's.
std::vector<double> degenerate(const DriftDiffusion& dd, const device::Device& d,
                               const assemble::Scaling& s, double noise, std::uint64_t seed) {
    std::vector<double> x = heated(dd, d, s, noise, seed);
    Noise z{seed + 100};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        const double u = d.mesh().points()[i][0] / 1e-4;
        x[5 * i + 1] = std::pow(10.0, 16.0 + 4.0 * u) / s.Ns * (1.0 + noise * z.next());
        x[5 * i + 2] = std::pow(10.0, 19.0 - 9.0 * u) / s.Ns * (1.0 + noise * z.next());
    }
    return x;
}

// The heated state hot enough to ionize: theta_n to 20 and theta_p to 15 along x.
std::vector<double> ionizing(const DriftDiffusion& dd, const device::Device& d,
                             const assemble::Scaling& s, double noise, std::uint64_t seed) {
    std::vector<double> x = heated(dd, d, s, noise, seed);
    Noise z{seed + 200};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        const double u = d.mesh().points()[i][0] / 1e-4;
        x[5 * i + 3] = (1.0 + 19.0 * u) * (1.0 + noise * z.next()) - 1.0;
        x[5 * i + 4] = (1.0 + 14.0 * u) * (1.0 + noise * z.next()) - 1.0;
    }
    return x;
}

}  // namespace

TEST_CASE("energy balance: five unknowns per node, tau_n and tau_p last") {
    const device::Device d = diode({});
    const auto s = *assemble::make_scaling(d);
    const DriftDiffusion dd = make(d, s, {});
    REQUIRE(dd.energy_balance());
    REQUIRE(dd.stride() == 5u);
    REQUIRE(dd.dofs().slot(assemble::FieldId::electron_temperature) == 3u);
    REQUIRE(dd.dofs().slot(assemble::FieldId::hole_temperature) == 4u);
    const DriftDiffusion off = make(d, s, {}, false);
    REQUIRE_FALSE(off.energy_balance());
    REQUIRE(off.stride() == 3u);
}

TEST_CASE("energy balance: the Jacobian against finite differences (EC6, Boltzmann)") {
    struct Case {
        const char* name;
        Options o;
    };
    const Case cases[] = {
        {"1D", {.D = 1, .nodes = 13}},
        {"1D, material step", {.D = 1, .nodes = 13, .hetero = true}},
        {"2D", {.D = 2, .nodes = 7}},
        {"2D, oxide", {.D = 2, .nodes = 7, .oxide = true}},
        {"2D, gate", {.D = 2, .nodes = 7, .gate = true}},
        {"3D, oxide", {.D = 3, .nodes = 5, .oxide = true}},
        {"3D", {.D = 3, .nodes = 5}},
    };
    for (const Case& c : cases) {
        CAPTURE(c.name);
        const device::Device d = diode(c.o);
        const auto s = *assemble::make_scaling(d);
        const DriftDiffusion dd = make(d, s, c.o);
        for (const std::uint64_t seed : {1u, 2u}) {
            const double error = fd_jacobian_error(dd, heated(dd, d, s, 0.05, seed));
            UNSCOPED_INFO(c.name << " worst column error " << error);
            REQUIRE(error <= 1e-6);
        }
    }
}

TEST_CASE("energy balance: the Jacobian against finite differences (EC6, Fermi-Dirac)") {
    const Options cases[] = {{.D = 1, .nodes = 13},
                             {.D = 1, .nodes = 13, .hetero = true},
                             {.D = 2, .nodes = 7, .oxide = true},
                             {.D = 3, .nodes = 5}};
    for (const Options& o : cases) {
        CAPTURE(o.D, o.hetero, o.oxide);
        const device::Device d = diode(o);
        const auto s = *assemble::make_scaling(d);
        const DriftDiffusion dd = make(d, s, o, {.fermi_dirac = true, .energy_balance = true});
        for (const std::uint64_t seed : {1u, 2u}) {
            REQUIRE(fd_jacobian_error(dd, degenerate(dd, d, s, 0.05, seed), true) <= 1e-6);
        }
    }
}

TEST_CASE("energy balance: the Jacobian against finite differences (EC6, impact ionization)") {
    const Options cases[] = {{.D = 1, .nodes = 13},
                             {.D = 1, .nodes = 13, .hetero = true},
                             {.D = 2, .nodes = 7, .oxide = true},
                             {.D = 3, .nodes = 5}};
    for (const Options& o : cases) {
        CAPTURE(o.D, o.hetero, o.oxide);
        const device::Device d = diode(o);
        const auto s = *assemble::make_scaling(d);
        for (const bool fd : {false, true}) {
            CAPTURE(fd);
            DriftDiffusion dd = make(
                d, s, o, {.fermi_dirac = fd, .impact_ionization = true, .energy_balance = true});
            const DriftDiffusion off = make(d, s, o, {.fermi_dirac = fd, .energy_balance = true});
            const std::vector<double> x = ionizing(dd, d, s, 0.05, 3);
            // The generation is there: some continuity row moves by more than 1e-3 of its size.
            std::vector<double> on(dd.unknowns()), without(dd.unknowns());
            dd.residual(x, on);
            off.residual(x, without);
            double moved = 0.0;
            for (std::size_t i = 0; i < dd.node_count(); ++i) {
                const double r = std::abs(on[5 * i + 1] - without[5 * i + 1]);
                moved = std::max(moved, r / std::max(std::abs(without[5 * i + 1]), 1e-300));
            }
            UNSCOPED_INFO("largest relative change of an electron row " << moved);
            REQUIRE(moved > 1e-3);
            REQUIRE(fd_jacobian_error(dd, x, true) <= 1e-6);
            // Held (ED12's first stage): the continuity rows only.
            dd.hold_temperatures(x);
            REQUIRE(fd_jacobian_error(dd, x, true) <= 1e-6);
        }
    }
}

TEST_CASE("energy balance: the Jacobian against finite differences (EC6, Lombardi)") {
    // Surface edges along the oxide: the E_perp stencil reaches every transport row of both ends.
    const Options cases[] = {{.D = 2, .nodes = 7, .oxide = true},
                             {.D = 2, .nodes = 9, .oxide = true},
                             {.D = 3, .nodes = 5, .oxide = true}};
    for (const Options& o : cases) {
        CAPTURE(o.D, o.nodes);
        const device::Device d = diode(o);
        const auto s = *assemble::make_scaling(d);
        for (const bool fd : {false, true}) {
            CAPTURE(fd);
            DriftDiffusion dd = make(d, s, o,
                                     {.surface_mobility = assemble::SurfaceMobility::lombardi,
                                      .fermi_dirac = fd,
                                      .energy_balance = true});
            const std::vector<double> x = fd ? degenerate(dd, d, s, 0.05, 5)
                                             : heated(dd, d, s, 0.05, 5);
            REQUIRE(fd_jacobian_error(dd, x, true) <= 1e-6);
            dd.hold_temperatures(x);
            REQUIRE(fd_jacobian_error(dd, x, true) <= 1e-6);
        }
    }
}

TEST_CASE("energy balance: a hot degenerate carrier at constant pressure carries no flux (ED8)") {
    // One contact (x_min), no recombination: the last node's electron row is the last edge's
    // flux alone. A flat band, theta from 1.5 to 1.6 across that edge, and n_b at the pressure of
    // n_a: N theta^(5/2) F_{3/2}(eta) equal at both ends, n = N theta^(3/2) F_{1/2}(eta). The
    // continuum flux is 0; the edge's is 0 to the third order of its differences (and exactly 0
    // under Boltzmann statistics). A wrong sign or factor on the (5/2)(g_e - 1) Delta ln theta
    // term misses by about (g - 1) Delta ln theta.
    const int nodes = 5;
    mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 0.1e-4, nodes));
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(nodes, 0),
                                     .donors = std::vector<double>(nodes, 0.0),
                                     .acceptors = std::vector<double>(nodes, 0.0),
                                     .contacts = {{"a", device::ContactKind::ohmic, {0}}}});
    REQUIRE(d.has_value());
    const auto s = *assemble::make_scaling(*d);
    const double N = physics::conduction_band_dos(physics::silicon(), 300.0);
    for (const bool fd : {false, true}) {
        for (const double eta_a : {-3.0, 2.0, 6.0}) {
            CAPTURE(fd, eta_a);
            auto dd = DriftDiffusion::create(*d, s,
                                             {.srh = false,
                                              .auger = false,
                                              .bgn = false,
                                              .fermi_dirac = fd,
                                              .radiative = false,
                                              .energy_balance = true});
            REQUIRE(dd.has_value());
            REQUIRE(dd->set_bias(std::vector<double>{0.0}).has_value());
            const double ta = 1.5, tb = 1.6;
            // Boltzmann: F_j(eta) = e^eta.
            const auto F32 = [&](double eta) {
                return fd ? physics::fermi_three_halves(eta).value : std::exp(eta);
            };
            const auto F12 = [&](double eta) {
                return fd ? physics::fermi_half(eta).value : std::exp(eta);
            };
            const double target = std::pow(ta, 2.5) * F32(eta_a) / std::pow(tb, 2.5);
            double eta_b = eta_a;  // Newton on ln F_{3/2}(eta_b) = ln target, dlnF32/deta = F12/F32
            for (int k = 0; k < 60; ++k) {
                eta_b -= (std::log(F32(eta_b)) - std::log(target)) * F32(eta_b) / F12(eta_b);
            }
            std::vector<double> x(dd->unknowns(), 0.0);
            for (std::size_t i = 0; i < static_cast<std::size_t>(nodes); ++i) {
                const bool last = i == static_cast<std::size_t>(nodes - 1);
                const double theta = last ? tb : ta, eta = last ? eta_b : eta_a;
                x[5 * i + 1] = N * std::pow(theta, 1.5) * F12(eta) / s.Ns;
                x[5 * i + 2] = 1e3 / s.Ns;
                x[5 * i + 3] = i == 0 ? 0.0 : theta - 1.0;
                x[5 * i + 4] = 0.0;
            }
            std::vector<double> f(dd->unknowns());
            auto jacobian = dd->make_jacobian();
            dd->evaluate(x, f, jacobian);
            const auto j = dense(jacobian);
            const std::size_t row = 5 * (nodes - 1) + 1, col = row;  // the last node's n
            const double scale = std::abs(j[row][col] * x[col]);  // one of the opposing terms
            const double miss = std::abs(f[row]) / scale;
            UNSCOPED_INFO((fd ? "FD" : "Boltzmann") << " eta_a " << eta_a << ": flux / term "
                                                    << miss);
            REQUIRE(miss <= (fd ? 1e-4 : 1e-14));
        }
    }
}

TEST_CASE("energy balance: at tau = 0 the continuity rows and currents are drift-diffusion's") {
    // Boltzmann and Fermi-Dirac (Unit 14's flux, ED8 as ruled); Lombardi along an oxide.
    const Options cases[] = {{.D = 1, .nodes = 13, .hetero = true},
                             {.D = 2, .nodes = 7},
                             {.D = 2, .nodes = 7, .oxide = true}};
    for (const Options& o : cases) {
        for (const bool fd : {false, true}) {
            CAPTURE(o.D, o.oxide, fd);
            const device::Device d = diode(o);
            const auto s = *assemble::make_scaling(d);
            const assemble::SurfaceMobility surface =
                o.oxide ? assemble::SurfaceMobility::lombardi : assemble::SurfaceMobility::off;
            const DriftDiffusion eb = make(
                d, s, o, {.surface_mobility = surface, .fermi_dirac = fd, .energy_balance = true});
            const DriftDiffusion dd =
                make(d, s, o, {.surface_mobility = surface, .fermi_dirac = fd});
            std::vector<double> x5 = fd ? degenerate(eb, d, s, 0.05, 7) : heated(eb, d, s, 0.05, 7);
            std::vector<double> x3(dd.unknowns());
            for (std::size_t i = 0; i < dd.node_count(); ++i) {
                x5[5 * i + 3] = x5[5 * i + 4] = 0.0;
                for (std::size_t r = 0; r < 3; ++r) x3[3 * i + r] = x5[5 * i + r];
            }
            std::vector<double> f5(eb.unknowns()), f3(dd.unknowns());
            eb.residual(x5, f5);
            dd.residual(x3, f3);
            for (std::size_t i = 0; i < dd.node_count(); ++i) {
                for (std::size_t r = 0; r < 3; ++r) {
                    CAPTURE(i, r);
                    REQUIRE(std::abs(f5[5 * i + r] - f3[3 * i + r]) <=
                            1e-13 * std::max(1.0, std::abs(f3[3 * i + r])));
                }
            }
            const auto c5 = eb.terminal_currents(x5), c3 = dd.terminal_currents(x3);
            for (std::size_t k = 0; k < c3.size(); ++k) {
                REQUIRE(std::abs(c5[k] - c3[k]) <= 1e-13 * std::abs(c3[k]));
            }
        }
    }
}

TEST_CASE("energy balance: equilibrium has tau = 0 in every energy row") {
    // At zero bias the drift-diffusion equilibrium with tau = 0 leaves every energy row at 0:
    // no flux, no Joule heat, no relaxation, no net recombination (EC5's assembly half).
    const Options o{.D = 2, .nodes = 9, .oxide = true};
    const device::Device d = diode(o);
    const auto s = *assemble::make_scaling(d);
    for (const bool fd : {false, true}) {
        CAPTURE(fd);
        auto dd = DriftDiffusion::create(d, s, {.fermi_dirac = fd, .energy_balance = true});
        REQUIRE(dd.has_value());
        REQUIRE(dd->set_bias(std::vector<double>{0.0, 0.0}).has_value());
        // Any potential with the equilibrium densities (one Fermi level): the charge-neutral one.
        // The fluxes vanish analytically; numerically to the rounding of their opposing terms.
        const auto poisson = *assemble::EquilibriumPoisson::create(d, s);
        auto x = dd->state_from_potential(poisson.charge_neutral_potential());
        std::vector<double> f(dd->unknowns());
        dd->residual(x, f);
        const auto terms = dd->edge_flux_terms(x);
        const double scale = *std::ranges::max_element(terms);
        REQUIRE(scale > 0.0);
        for (std::size_t i = 0; i < dd->node_count(); ++i) {
            CAPTURE(i);
            REQUIRE(std::abs(f[5 * i + 3]) <= 1e-13 * scale);
            REQUIRE(std::abs(f[5 * i + 4]) <= 1e-13 * scale);
        }
    }
}

TEST_CASE("energy balance: the uniform bulk state of Canali's velocity is exact (EC3 at 1e5 V/cm)") {
    // EC3 at 1e5 V/cm, as ruled 2026-10-29: in saturation the cold carriers entering at a contact
    // are not screened within a resistor of practical length, so the uniform state is checked
    // where it is exact, in the discrete bulk. The neutral densities of each doping, a potential
    // falling 1e5 V/cm, and each carrier's rise at the local balance T_L + (2/3) q tau_w v E / k
    // with v = mu_C(mu0, E) E: every row of an interior node vanishes to its terms' rounding.
    constexpr double E = 1e5, k_over_q = base::k_B_J_per_K / base::q_C, T_L = 300.0;
    const physics::Semiconductor si = physics::silicon();
    const auto tau_w = physics::silicon_energy_relaxation;
    for (const double N : {1e16, 1e17, 1e18}) {
        for (const bool donors : {true, false}) {
            CAPTURE(N, donors);
            mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 0.4e-4, 41));
            const std::size_t n = m.node_count();
            auto d = device::Device::create(
                {.mesh = std::move(m),
                 .temperature_K = T_L,
                 .regions = {{"silicon", si}},
                 .node_region = std::vector<device::RegionId>(n, 0),
                 .donors = std::vector<double>(n, donors ? N : 0.0),
                 .acceptors = std::vector<double>(n, donors ? 0.0 : N),
                 .contacts = {{"left", device::ContactKind::ohmic, {0}},
                              {"right", device::ContactKind::ohmic,
                               {static_cast<mesh::NodeId>(n - 1)}}}});
            REQUIRE(d.has_value());
            const auto s = *assemble::make_scaling(*d);
            // Caughey-Thomas for both carriers, so mu0 is physics::caughey_thomas_mobility's.
            auto dd = DriftDiffusion::create(
                *d, s, {.hole_mobility = assemble::HoleMobility::caughey_thomas, .energy_balance = true});
            REQUIRE(dd.has_value());
            REQUIRE(dd->set_bias(std::vector<double>{0.0, E * 0.4e-4}).has_value());
            const auto poisson = *assemble::EquilibriumPoisson::create(*d, s);
            const std::vector<double> psi = poisson.charge_neutral_potential();
            std::vector<double> x = dd->state_from_potential(psi);  // the contacts at their bias
            // A carrier below 1e-10 of the majority relaxes with 1e-10 of the majority's density
            // (ED14): its Joule heat (its density) balances that relaxation, so its rise is the
            // local one times share = density / (1e-10 majority), solved by fixed point (share is
            // at most 1e-2 here, so the map contracts).
            const std::size_t mid = n / 2;
            const double n0 = x[5 * mid + 1], p0 = x[5 * mid + 2];
            const auto rise = [&](physics::Carrier c, double tau, double density, double other) {
                const double mu0 = physics::caughey_thomas_mobility(si, c, N, T_L);
                const physics::CanaliParameters& sat = physics::saturation(si, c);
                if (density >= 1e-10 * other) {
                    const double v = physics::canali_mobility(mu0, E, sat).mobility * E;
                    return 2.0 / 3.0 * tau * v * E / k_over_q / T_L;
                }
                const double share = density / (1e-10 * other);
                double t = 0.0;
                for (int k = 0; k < 200; ++k) {
                    const double mu =
                        physics::carrier_temperature_mobility(mu0, t * T_L, tau, sat).mobility;
                    t = 2.0 / 3.0 * tau * mu * E * E / k_over_q / T_L * share;
                }
                return t;
            };
            const double tn = rise(physics::Carrier::electron, tau_w.tau_n_s, n0, p0);
            const double tp = rise(physics::Carrier::hole, tau_w.tau_p_s, p0, n0);
            for (std::size_t i = 0; i < n; ++i) {
                x[5 * i] = psi[i] + E * d->mesh().points()[i][0] / s.V_T;
                x[5 * i + 3] = tn;
                x[5 * i + 4] = tp;
            }
            std::vector<double> f(dd->unknowns());
            auto jacobian = dd->make_jacobian();
            dd->evaluate(x, f, jacobian);
            const auto j = dense(jacobian);
            double worst = 0.0;
            for (std::size_t r = 5; r < 5 * (n - 1); ++r) {
                double scale = 0.0;  // the size of the row's terms: sum over c of |J_rc x_c|
                for (std::size_t c = 0; c < j[r].size(); ++c) scale += std::abs(j[r][c] * x[c]);
                CAPTURE(r);
                REQUIRE(scale > 0.0);
                worst = std::max(worst, std::abs(f[r]) / scale);
            }
            UNSCOPED_INFO("N " << N << (donors ? " n" : " p") << ": worst row " << worst);
            REQUIRE(worst <= 1e-12);
            // EC3's impact ionization (ED6) at the same state: the generation each interior
            // node's electron row gains equals the field-driven model's (Unit 19 with Canali's
            // mobility) at the state's psi, n, p, within 1e-3. E_eff(T) = E there.
            constexpr auto ct = assemble::HoleMobility::caughey_thomas;
            auto eb_ii = DriftDiffusion::create(
                *d, s, {.hole_mobility = ct, .impact_ionization = true, .energy_balance = true});
            auto dd_ii = DriftDiffusion::create(
                *d, s, {.hole_mobility = ct, .field_mobility = true, .impact_ionization = true});
            auto dd_off =
                DriftDiffusion::create(*d, s, {.hole_mobility = ct, .field_mobility = true});
            REQUIRE((eb_ii.has_value() && dd_ii.has_value() && dd_off.has_value()));
            for (DriftDiffusion* q : {&*eb_ii, &*dd_ii, &*dd_off}) {
                REQUIRE(q->set_bias(std::vector<double>{0.0, E * 0.4e-4}).has_value());
            }
            std::vector<double> x3(dd_ii->unknowns());
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t r = 0; r < 3; ++r) x3[3 * i + r] = x[5 * i + r];
            }
            std::vector<double> f_eb(eb_ii->unknowns()), f_dd(dd_ii->unknowns()),
                f_off(dd_off->unknowns());
            eb_ii->residual(x, f_eb);
            dd_ii->residual(x3, f_dd);
            dd_off->residual(x3, f_off);
            double generation = 0.0;
            for (std::size_t i = 1; i + 1 < n; ++i) {
                const double g_eb = f_eb[5 * i + 1] - f[5 * i + 1];
                const double g_dd = f_dd[3 * i + 1] - f_off[3 * i + 1];
                CAPTURE(i, g_eb, g_dd);
                REQUIRE(g_dd > 0.0);
                generation = std::max(generation, std::abs(g_eb - g_dd) / g_dd);
            }
            UNSCOPED_INFO("N " << N << ": worst generation " << generation);
            REQUIRE(generation <= 1e-3);
        }
    }
}

TEST_CASE("energy balance: the relaxation floor and held temperatures (ED14, ED12)") {
    const Options o{.D = 2, .nodes = 7};
    const device::Device d = diode(o);
    const auto s = *assemble::make_scaling(d);
    DriftDiffusion dd = make(d, s, o);
    // Holes 1e-12 of the electrons on the n side, electrons 1e-12 of the holes on the p side: the
    // floor is active there. The Jacobian gate holds across it. On odd nodes the majority is 1e18,
    // above the doping (1e17, 2e17), so the floor follows it; on even ones (1e16, 1e15) the doping.
    std::vector<double> x = heated(dd, d, s, 0.05, 11);
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        const bool n_side = d.mesh().points()[i][0] >= 0.5e-4;
        if (i % 2 == 1) x[5 * i + (n_side ? 1 : 2)] = 1e18 / s.Ns;
        (n_side ? x[5 * i + 2] : x[5 * i + 1]) = 1e-12 * (n_side ? x[5 * i + 1] : x[5 * i + 2]);
    }
    REQUIRE(fd_jacobian_error(dd, x, true) <= 1e-6);
    // The floor's own partial, d/d(majority) of the floored carrier's energy row, is far below the
    // majority column's Poisson entry, so it is checked entry by entry: a central difference in the
    // majority density against the entry, within 1e-6 of it. Without recombination, whose rate far
    // from n p = n_i^2 would swamp the row's rounding.
    {
        auto plain = DriftDiffusion::create(
            d, s, {.srh = false, .auger = false, .radiative = false, .energy_balance = true});
        REQUIRE(plain.has_value());
        REQUIRE(plain->set_bias(biases(o)).has_value());
        const DriftDiffusion& q = *plain;
        std::vector<double> f(q.unknowns()), fp(q.unknowns()), fm(q.unknowns());
        auto jacobian = q.make_jacobian();
        q.evaluate(x, f, jacobian);
        const auto j = dense(jacobian);
        std::size_t checked = 0, doping_set = 0;
        for (std::size_t i = 0; i < q.node_count(); ++i) {
            const double px = d.mesh().points()[i][0];
            if (px < 1e-12 || px > 1e-4 - 1e-12) continue;  // the contacts (x_min, x_max)
            // Next to the junction a neighbour's majority flux swamps the row's rounding.
            if (std::abs(px - 0.5e-4) < 0.2e-4) continue;
            const bool n_side = px >= 0.5e-4;
            const std::size_t row = 5 * i + (n_side ? 4 : 3), col = 5 * i + (n_side ? 1 : 2);
            const double base = x[col], h = 1e-6 * base;
            x[col] = base + h;
            q.residual(x, fp);
            x[col] = base - h;
            q.residual(x, fm);
            x[col] = base;
            const double fd = (fp[row] - fm[row]) / (2.0 * h);
            CAPTURE(i, fd, j[row][col], f[row], fp[row], fm[row], x[5 * i + 1], x[5 * i + 2],
                    x[5 * i + 3], x[5 * i + 4]);
            // Where the net doping exceeds the majority the floor follows the doping (ED14 as
            // amended) and the partial is 0.
            if (fd == 0.0) {
                REQUIRE(j[row][col] == 0.0);
                ++doping_set;
                continue;
            }
            REQUIRE(std::abs(fd - j[row][col]) <= 1e-6 * std::abs(j[row][col]));
            ++checked;
        }
        REQUIRE(checked > 0);
        REQUIRE(doping_set > 0);
    }
    // Held: the energy rows are tau - tau_held, the other rows those of the released system.
    std::vector<double> released(dd.unknowns()), held_f(dd.unknowns());
    dd.residual(x, released);
    std::vector<double> y = x;
    dd.hold_temperatures(x);
    REQUIRE(dd.temperatures_held());
    for (std::size_t i = 0; i < dd.node_count(); ++i) y[5 * i + 3] += 0.25;
    dd.residual(y, held_f);
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        CAPTURE(i);
        REQUIRE(held_f[5 * i + 3] == 0.25);
        REQUIRE(held_f[5 * i + 4] == 0.0);
    }
    dd.residual(x, held_f);
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        for (std::size_t r = 0; r < 3; ++r) REQUIRE(held_f[5 * i + r] == released[5 * i + r]);
    }
    REQUIRE(fd_jacobian_error(dd, y, true) <= 1e-6);
    dd.release_temperatures();
    REQUIRE_FALSE(dd.temperatures_held());
    dd.residual(x, held_f);
    REQUIRE(held_f == released);
}

TEST_CASE("energy balance: the inputs it refuses (ED11, EC15)") {
    const auto refused = [](const Options& o, assemble::PhysicsModels m, const std::string& what) {
        const device::Device d = diode(o);
        const auto s = *assemble::make_scaling(d);
        m.energy_balance = true;
        const auto dd = DriftDiffusion::create(d, s, m);
        REQUIRE_FALSE(dd.has_value());
        CAPTURE(dd.error().message);
        REQUIRE(dd.error().code == base::ErrorCode::invalid_input);
        REQUIRE(dd.error().message.find(what) != std::string::npos);
    };
    refused({}, {.electrothermal = true}, "electrothermal");
    refused({}, {.field_mobility = true}, "field_mobility");
    refused({}, {.btbt_local = true}, "band-to-band");
    refused({.schottky = true}, {}, "Schottky");
    refused({.hetero = true, .thermionic = true}, {}, "thermionic");
    refused({.D = 2, .oxide = true, .interface = true}, {}, "interface charge");
    // A semiconductor that is not silicon and gives no relaxation times.
    {
        physics::SemiconductorParameters p = physics::silicon_parameters;
        p.electron_mobility.mu_max *= 0.9;
        const auto x = uniform(0.0, 1e-4, 5);
        mesh::Mesh m = *mesh::make_tensor_grid(x);
        auto d = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"other", *physics::Semiconductor::create(p)}},
             .node_region = std::vector<device::RegionId>(5, 0),
             .donors = std::vector<double>(5, 1e17),
             .acceptors = std::vector<double>(5, 0.0),
             .contacts = {{"a", device::ContactKind::ohmic, {0}},
                          {"b", device::ContactKind::ohmic, {4}}}});
        REQUIRE(d.has_value());
        const auto dd =
            DriftDiffusion::create(*d, *assemble::make_scaling(*d), models());
        REQUIRE_FALSE(dd.has_value());
        REQUIRE(dd.error().message.find("relaxation") != std::string::npos);
    }
}
