// Electrothermal coupling in the drift-diffusion assembler (ARCHITECTURE.md section 11, Unit 23;
// DECISIONS.md T1-T14): the finite-difference Jacobian gate of the four-block system in 1D, 2D
// and 3D, across a material step and into an oxide; at theta = 1 the rows and currents of the
// isothermal system to round-off; at a uniform temperature away from T0 no current at
// equilibrium; the Joule part of every edge's heat is not negative; and the inputs the model
// refuses (T8, T11) or does not support yet.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

// A second semiconductor, silicon with its bands, densities of states, conductivity law and
// thermopower exponents changed: every per-material term of an edge between the two differs.
// With fermi, the thermopower exponents stay -1/2 (the Fermi-Dirac thermopower needs them).
physics::Semiconductor other_material(bool fermi = false) {
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.electron_affinity_eV += 0.08;
    p.Eg0_eV += 0.05;
    p.varshni_alpha_eV_per_K *= 1.2;
    p.Nc300 *= 1.6;
    p.Nv300 *= 0.7;
    p.electron_mobility.mu_max *= 0.8;
    p.thermal = {.conductivity_W_cmK = 0.9, .conductivity_exponent = 1.1,
                 .heat_capacity_J_cm3K = 1.7, .thermopower_exponent_n = -0.3,
                 .thermopower_exponent_p = 0.2};
    if (fermi) p.thermal.thermopower_exponent_n = p.thermal.thermopower_exponent_p = -0.5;
    return *physics::Semiconductor::create(p);
}

struct Options {
    int D = 1;
    int nodes = 11;
    bool hetero = false;  // x >= 0.6 um in other_material()
    bool oxide = false;   // 2D: y >= 0.3 um oxide (contacts on the silicon rows only)
    bool sinks = true;    // isothermal x_min, R_th x_max
    bool fermi = false;   // Fermi-Dirac statistics (the models, and other_material's exponents)
    bool field = false;   // field-dependent mobility (the models)
    bool gate = false;    // 2D: a p-polysilicon gate on y_max (inside the x ends)
    bool electrode = false;  // 2D with oxide: a p-polysilicon electrode on y_max
    bool thermionic = false;  // with hetero: the material step declared thermionic_emission
    bool interface = false;   // with oxide: silicon-oxide fixed charge and surface recombination
    bool traps = false;       // with interface: an interface trap level too (refused, T8)
    bool ionization = false;  // incomplete ionization (the models)
};

// The models of the fixture's options, with or without the electrothermal model.
assemble::PhysicsModels models(const Options& o, bool thermal = true) {
    return {.field_mobility = o.field,
            .fermi_dirac = o.fermi,
            .incomplete_ionization = o.ionization,
            .electrothermal = thermal};
}

// The biases of the fixture's contacts: anode, cathode, and the gate or electrode.
std::vector<double> biases(const Options& o) {
    std::vector<double> v{0.4, -0.3};
    if (o.gate || o.electrode) v.push_back(0.7);
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
    // Region ids: silicon 0, then "other" and "oxide" when used.
    const device::RegionId other = 1, oxide = o.hetero ? 2 : 1;
    std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    if (o.hetero) regions.push_back({"other", other_material(o.fermi)});
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
    if (o.gate || o.electrode) {
        std::vector<mesh::NodeId> top;
        for (const mesh::NodeId k : m.find_boundary("y_max")->nodes) {
            const double xk = m.points()[static_cast<std::size_t>(k)][0];
            if (o.electrode || (xk > 1e-12 && xk < 1e-4 - 1e-12)) top.push_back(k);
        }
        if (o.gate) {
            contacts.push_back({"gate", device::ContactKind::gate, std::move(top),
                                {.boundary = "y_max", .oxide_thickness_cm = 5e-7,
                                 .electrode = device::GateElectrode::p_poly}});
        } else {
            contacts.push_back({"electrode", device::ContactKind::electrode, std::move(top), {},
                                {.kind = device::GateElectrode::p_poly}});
        }
    }
    std::vector<device::Interface> interfaces;
    if (o.interface) {
        interfaces.push_back({.region_a = "silicon", .region_b = "oxide",
                              .fixed_charge_cm2 = 1e11, .recombination_velocity_n_cm_s = 1e4,
                              .recombination_velocity_p_cm_s = 2e4});
        if (o.traps) interfaces.back().traps.levels = {{.density_cm2 = 1e10, .energy_eV = 0.1}};
    }
    if (o.thermionic) {
        interfaces.push_back(
            {"silicon", "other", device::InterfaceTransport::thermionic_emission});
    }
    std::vector<device::ThermalContact> thermal;
    if (o.sinks) {
        thermal = {{"sink", "x_min", device::ThermalContactKind::isothermal, 310.0, 0.0},
                   {"package", "x_max", device::ThermalContactKind::resistance, 290.0, 1e-4}};
    }
    auto device = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = std::move(regions),
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = std::move(contacts),
         .interfaces = interfaces,
         .thermal_contacts = std::move(thermal)});
    if (!device) FAIL(device.error().message);
    REQUIRE(device->regions().size() == 1u + (o.hetero ? 1u : 0u) + (o.oxide ? 1u : 0u));
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

// A biased, heated state: psi falling `drop` (20) V_T across the device, n (1e16) and p 1e15
// cm^-3, theta rising from 1 to 1.3 along x, every unknown with relative noise. Insulator
// densities 0; the contact nodes are not stamped (their rows are Dirichlet rows).
std::vector<double> heated(const DriftDiffusion& dd, const device::Device& d,
                           const assemble::Scaling& s, double noise, std::uint64_t seed,
                           double theta_rise = 0.3, double n_cm3 = 1e16, double drop = 20.0) {
    const std::size_t m = dd.stride();
    std::vector<double> x(dd.unknowns());
    Noise z{seed};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        const double u = d.mesh().points()[i][0] / 1e-4;
        x[m * i] = -drop * u * (1.0 + noise * z.next());
        const bool insulator = d.is_insulator(static_cast<mesh::NodeId>(i));
        x[m * i + 1] = insulator ? 0.0 : n_cm3 / s.Ns * (1.0 + noise * z.next());
        x[m * i + 2] = insulator ? 0.0 : 1e15 / s.Ns * (1.0 + noise * z.next());
        if (m == 4) {  // tau
            x[4 * i + 3] = (1.0 + theta_rise * u) * (1.0 + noise * z.next()) - 1.0;
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

// The finite-difference gate (as impact_ionization_test.cpp): every column, central differences,
// the worst error over the column's largest entry, each column's error the smallest at steps of
// 1e-5, 1e-6 and 1e-7 max(|x|, 1).
// With a time step, of the step's rows.
double fd_jacobian_error(const DriftDiffusion& dd, std::vector<double> x,
                         const DriftDiffusion::TimeStep* step = nullptr) {
    const std::size_t n = dd.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    auto jacobian = dd.make_jacobian();
    if (step != nullptr) {
        dd.evaluate(x, *step, f, jacobian);
    } else {
        dd.evaluate(x, f, jacobian);
    }
    const auto residual = [&](std::span<double> out) {
        if (step != nullptr) {
            dd.residual(x, *step, out);
        } else {
            dd.residual(x, out);
        }
    };
    const auto j = dense(jacobian);
    double worst = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        double scale = 1e-300;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][c]));
        double best = std::numeric_limits<double>::infinity();
        for (const double relative : {1e-5, 1e-6, 1e-7}) {
            const double base = x[c];
            const double delta = relative * std::max(std::abs(base), 1.0);
            x[c] = base + delta;
            const double up = x[c];
            residual(fp);
            x[c] = base - delta;
            const double down = x[c];
            residual(fm);
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

}  // namespace

TEST_CASE("electrothermal assemble: the stride and the pattern") {
    const device::Device d = diode({});
    const auto s = *assemble::make_scaling(d);
    const auto off = *DriftDiffusion::create(d, s);
    const auto on = *DriftDiffusion::create(d, s, {.electrothermal = true});
    REQUIRE(off.stride() == 3);
    REQUIRE_FALSE(off.electrothermal());
    REQUIRE(on.stride() == 4);
    REQUIRE(on.electrothermal());
    REQUIRE(on.unknowns() == 4 * d.mesh().node_count());
    REQUIRE(on.heat_sinks_complete());
    REQUIRE_FALSE(off.heat_sinks_complete());
    // Full 4 x 4 blocks: 16 per node and 32 per edge.
    const auto j = on.make_jacobian();
    REQUIRE(j.nonzeros() == 16 * d.mesh().node_count() + 32 * d.mesh().edges().size());
    // Without thermal contacts the model has no sink.
    const device::Device bare = diode({.sinks = false});
    REQUIRE_FALSE(
        DriftDiffusion::create(bare, s, {.electrothermal = true})->heat_sinks_complete());
}

TEST_CASE("electrothermal assemble: the finite-difference Jacobian gate") {
    // Every block of the four-unknown system: fluxes with the thermopower, the edge-varying T,
    // the T-dependent mobility and n_ie, recombination, the ohmic contacts at their T, the heat
    // rows (conduction in both forms, carrier energies, Peltier at the contacts) and the sinks.
    struct Case {
        Options o;
        const char* name;
    };
    // Fermi-Dirac with degenerate electrons (1e20 cm^-3, about 3 kT above the band edge): the
    // degeneracy and thermal diffusion terms, and the contacts' implicit neutral equilibrium.
    const Case cases[] = {
        {{.D = 1}, "1D"},
        {{.D = 1, .hetero = true}, "1D two materials"},
        {{.D = 2, .nodes = 9}, "2D"},
        {{.D = 2, .nodes = 9, .hetero = true, .oxide = true}, "2D with oxide"},
        {{.D = 3, .nodes = 6, .hetero = true}, "3D"},
        {{.D = 1, .hetero = true, .fermi = true}, "1D Fermi-Dirac"},
        {{.D = 2, .nodes = 9, .hetero = true, .oxide = true, .fermi = true}, "2D Fermi-Dirac"},
        {{.D = 1, .hetero = true, .field = true}, "1D field mobility"},
        {{.D = 2, .nodes = 9, .hetero = true, .fermi = true, .field = true},
         "2D Fermi-Dirac and field mobility"},
        {{.D = 2, .nodes = 9, .gate = true}, "2D p-polysilicon gate"},
        {{.D = 2, .nodes = 9, .hetero = true, .oxide = true, .electrode = true},
         "2D p-polysilicon electrode"},
        {{.D = 1, .hetero = true, .thermionic = true}, "1D thermionic emission"},
        {{.D = 2, .nodes = 9, .hetero = true, .fermi = true, .thermionic = true},
         "2D thermionic emission, Fermi-Dirac"},
        {{.D = 2, .nodes = 9, .oxide = true, .interface = true}, "2D interface"},
        {{.D = 2, .nodes = 9, .oxide = true, .fermi = true, .interface = true},
         "2D interface, Fermi-Dirac"},
        {{.D = 1, .hetero = true, .ionization = true}, "1D incomplete ionization"},
        {{.D = 2, .nodes = 9, .hetero = true, .fermi = true, .ionization = true},
         "2D incomplete ionization, Fermi-Dirac"}};
    for (const Case& c : cases) {
        const device::Device d = diode(c.o);
        const auto s = *assemble::make_scaling(d);
        auto made = DriftDiffusion::create(d, s, models(c.o));
        if (!made) FAIL(c.name << ": " << made.error().message);
        auto& dd = *made;
        REQUIRE(dd.set_bias(biases(c.o)).has_value());
        const auto x = heated(dd, d, s, 0.05, 11, 0.3, c.o.fermi ? 1e20 : 1e16);
        const double error = fd_jacobian_error(dd, x);
        CAPTURE(c.name, error);
        REQUIRE(error <= 1e-6);
    }
}

TEST_CASE("electrothermal assemble: the finite-difference Jacobian gate of a time step") {
    // The storage terms: carriers in the continuity rows, rho c dT/dt and the carriers' stored
    // energy in the heat rows, from a history state that differs from the new one everywhere.
    for (const Options o : {Options{.D = 1, .hetero = true},
                            Options{.D = 2, .nodes = 9, .hetero = true, .oxide = true},
                            Options{.D = 1, .hetero = true, .fermi = true, .field = true},
                            Options{.D = 2, .nodes = 9, .gate = true},
                            Options{.D = 2, .nodes = 9, .oxide = true, .electrode = true},
          Options{.D = 1, .hetero = true, .thermionic = true},
          Options{.D = 1, .hetero = true, .fermi = true, .thermionic = true},
          Options{.D = 2, .nodes = 9, .oxide = true, .interface = true},
          Options{.D = 2, .nodes = 9, .oxide = true, .fermi = true, .interface = true},
          Options{.D = 1, .hetero = true, .ionization = true},
          Options{.D = 1, .hetero = true, .fermi = true, .ionization = true}}) {
        const device::Device d = diode(o);
        const auto s = *assemble::make_scaling(d);
        auto dd = *DriftDiffusion::create(d, s, models(o));
        REQUIRE(dd.set_bias(biases(o)).has_value());
        const double n_cm3 = o.fermi ? 1e20 : 1e16;
        const auto before = heated(dd, d, s, 0.05, 5, 0.25, n_cm3);
        const auto x = heated(dd, d, s, 0.05, 11, 0.3, n_cm3);
        const std::vector<double> storage = dd.storage(before);
        REQUIRE(storage.size() == 3 * dd.node_count());
        REQUIRE(dd.storage_width() == 3);
        const DriftDiffusion::TimeStep step{0.7, storage, {}};
        const double error = fd_jacobian_error(dd, x, &step);
        CAPTURE(o.D, o.fermi, error);
        REQUIRE(error <= 1e-6);
    }
}

TEST_CASE("electrothermal assemble: at theta = 1 the isothermal rows and currents") {
    for (const Options o :
         {Options{.D = 1}, Options{.D = 1, .hetero = true},
          Options{.D = 2, .nodes = 9, .hetero = true, .oxide = true},
          Options{.D = 1, .hetero = true, .fermi = true},
          Options{.D = 2, .nodes = 9, .hetero = true, .oxide = true, .fermi = true},
          Options{.D = 1, .hetero = true, .fermi = true, .field = true},
          Options{.D = 2, .nodes = 9, .gate = true},
          Options{.D = 2, .nodes = 9, .oxide = true, .electrode = true},
          Options{.D = 1, .hetero = true, .thermionic = true},
          Options{.D = 1, .hetero = true, .fermi = true, .thermionic = true},
          Options{.D = 2, .nodes = 9, .oxide = true, .interface = true},
          Options{.D = 2, .nodes = 9, .oxide = true, .fermi = true, .interface = true},
          Options{.D = 1, .hetero = true, .ionization = true},
          Options{.D = 1, .hetero = true, .fermi = true, .ionization = true}}) {
        const device::Device d = diode(o);
        const auto s = *assemble::make_scaling(d);
        auto iso = *DriftDiffusion::create(d, s, models(o, false));
        auto on = *DriftDiffusion::create(d, s, models(o));
        const std::vector<double> bias = biases(o);
        REQUIRE(iso.set_bias(bias).has_value());
        REQUIRE(on.set_bias(bias).has_value());
        const auto x3 = heated(iso, d, s, 0.05, 3, 0.3, o.fermi ? 1e20 : 1e16);
        std::vector<double> x4(on.unknowns());
        for (std::size_t i = 0; i < on.node_count(); ++i) {
            for (std::size_t r = 0; r < 3; ++r) x4[4 * i + r] = x3[3 * i + r];
            x4[4 * i + 3] = 0.0;  // tau: theta = 1
        }
        std::vector<double> f3(iso.unknowns()), f4(on.unknowns());
        iso.residual(x3, f3);
        on.residual(x4, f4);
        // Each row against the scale of its terms (the largest edge flux in it, or its density).
        const auto c3 = iso.edge_currents(x3), c4 = on.edge_currents(x4);
        double flux_scale = 0.0;
        for (const auto& [jn, jp] : c3) {
            flux_scale = std::max({flux_scale, std::abs(jn), std::abs(jp)});
        }
        double worst_row = 0.0, worst_current = 0.0;
        for (std::size_t i = 0; i < on.node_count(); ++i) {
            for (std::size_t r = 0; r < 3; ++r) {
                const double scale = std::max({flux_scale, std::abs(x3[3 * i + r]), 1.0});
                worst_row = std::max(worst_row, std::abs(f4[4 * i + r] - f3[3 * i + r]) / scale);
            }
        }
        for (std::size_t k = 0; k < c3.size(); ++k) {
            worst_current = std::max({worst_current, std::abs(c4[k].first - c3[k].first),
                                      std::abs(c4[k].second - c3[k].second)});
        }
        CAPTURE(o.D, o.hetero, o.fermi, worst_row, worst_current / flux_scale);
        REQUIRE(worst_row <= 1e-12);
        REQUIRE(worst_current <= 1e-12 * flux_scale);
    }
}

TEST_CASE("electrothermal assemble: no current at equilibrium at a uniform temperature") {
    // At theta = 1.25 with both quasi-Fermi levels at 0 (each density set from the band edges),
    // over a potential profile and across the material step: every flux vanishes to round-off of
    // its one-sided terms. With a temperature step instead, the fluxes do not vanish.
    // Under Fermi-Dirac statistics too, with the potential swinging the electrons degenerate.
    // And across a thermionic-emission step.
    for (const auto& [fermi, hetero, emission] :
         {std::tuple{false, false, false}, std::tuple{false, true, false},
          std::tuple{false, true, true}, std::tuple{true, false, false},
          std::tuple{true, true, false}, std::tuple{true, true, true}}) {
        const device::Device d =
            diode({.D = 1, .hetero = hetero, .fermi = fermi, .thermionic = emission});
        const auto s = *assemble::make_scaling(d);
        const auto dd =
            *DriftDiffusion::create(d, s, {.fermi_dirac = fermi, .electrothermal = true});
        for (const double step : {0.0, 0.2}) {
            std::vector<double> x(dd.unknowns());
            for (std::size_t i = 0; i < dd.node_count(); ++i) {
                const double u = d.mesh().points()[i][0] / 1e-4;
                x[4 * i] = 30.0 * std::sin(3.0 * u);
                x[4 * i + 3] = 0.25 + (u > 0.45 ? step : 0.0);  // tau: theta 1.25
            }
            // Each density by bisection on its logarithm until its quasi-Fermi level is 0
            // (E_Fn rises with n, E_Fp falls with p).
            std::vector<double> lo(2 * dd.node_count(), -300.0), hi(2 * dd.node_count(), 60.0);
            for (int it = 0; it < 200; ++it) {
                for (std::size_t i = 0; i < dd.node_count(); ++i) {
                    x[4 * i + 1] = std::exp(0.5 * (lo[2 * i] + hi[2 * i]));
                    x[4 * i + 2] = std::exp(0.5 * (lo[2 * i + 1] + hi[2 * i + 1]));
                }
                const auto bands = dd.band_edges(x);
                for (std::size_t i = 0; i < dd.node_count(); ++i) {
                    const double mn = 0.5 * (lo[2 * i] + hi[2 * i]);
                    const double mp = 0.5 * (lo[2 * i + 1] + hi[2 * i + 1]);
                    (bands.electron_fermi[i] > 0.0 ? hi[2 * i] : lo[2 * i]) = mn;
                    (bands.hole_fermi[i] < 0.0 ? hi[2 * i + 1] : lo[2 * i + 1]) = mp;
                }
            }
            const auto again = dd.band_edges(x);
            double largest = 0.0;
            for (std::size_t i = 0; i < dd.node_count(); ++i) {
                REQUIRE(std::abs(again.electron_fermi[i]) <= 1e-12);
                REQUIRE(std::abs(again.hole_fermi[i]) <= 1e-12);
                largest = std::max(largest, x[4 * i + 1] * s.Ns);
            }
            if (fermi) REQUIRE(largest > 1e20);  // degenerate somewhere
            const auto c = dd.edge_currents(x);
            double worst = 0.0;
            for (std::size_t k = 0; k < c.size(); ++k) {
                const auto& e = d.mesh().edges()[k];
                const auto a = static_cast<std::size_t>(e.first);
                const auto b = static_cast<std::size_t>(e.second);
                // The one-sided terms: of the size of the densities times the edge factor.
                const double size_n = std::max(x[4 * a + 1], x[4 * b + 1]);
                const double size_p = std::max(x[4 * a + 2], x[4 * b + 2]);
                worst = std::max({worst, std::abs(c[k].first) / size_n,
                                  std::abs(c[k].second) / size_p});
            }
            CAPTURE(fermi, hetero, emission, step, worst);
            if (step == 0.0) {
                REQUIRE(worst <= 1e-9);
            } else {
                REQUIRE(worst >= 1e-6);
            }
        }
    }
}

TEST_CASE("electrothermal assemble: impact ionization with gamma(T)") {
    // A field of about 3e5 V/cm (1160 V_T across 1 um): generation at every node, its temperature
    // dependence through the currents and gamma. The FD gate; and at theta = 1 the isothermal
    // rows, generation included.
    for (const int D : {1, 2}) {
        for (const bool fermi : {false, true}) {
            const Options o{.D = D, .nodes = D == 1 ? 11 : 9};
            const device::Device d = diode(o);
            const auto s = *assemble::make_scaling(d);
            const assemble::PhysicsModels on_models{
                .fermi_dirac = fermi, .impact_ionization = true, .electrothermal = true};
            auto dd = *DriftDiffusion::create(d, s, on_models);
            REQUIRE(dd.set_bias(biases(o)).has_value());
            const auto x = heated(dd, d, s, 0.05, 7, 0.3, 1e16, 1160.0);
            const double error = fd_jacobian_error(dd, x);
            CAPTURE(D, fermi, error);
            REQUIRE(error <= 1e-6);
            // The generation is there: the rows differ from the model without it.
            auto off_models = on_models;
            off_models.impact_ionization = false;
            auto plain = *DriftDiffusion::create(d, s, off_models);
            REQUIRE(plain.set_bias(biases(o)).has_value());
            std::vector<double> f_on(dd.unknowns()), f_off(dd.unknowns());
            dd.residual(x, f_on);
            plain.residual(x, f_off);
            double generation = 0.0;
            for (std::size_t k = 0; k < f_on.size(); ++k) {
                generation = std::max(generation, std::abs(f_on[k] - f_off[k]));
            }
            REQUIRE(generation > 0.0);
            // theta = 1: the isothermal system's rows.
            auto iso = *DriftDiffusion::create(
                d, s, {.fermi_dirac = fermi, .impact_ionization = true});
            REQUIRE(iso.set_bias(biases(o)).has_value());
            std::vector<double> x3(iso.unknowns()), x4 = x;
            for (std::size_t i = 0; i < iso.node_count(); ++i) {
                for (std::size_t r = 0; r < 3; ++r) x3[3 * i + r] = x[4 * i + r];
                x4[4 * i + 3] = 0.0;
            }
            std::vector<double> f3(iso.unknowns()), f4(dd.unknowns());
            iso.residual(x3, f3);
            dd.residual(x4, f4);
            double scale = 0.0, worst = 0.0;
            for (const auto& [jn, jp] : iso.edge_currents(x3)) {
                scale = std::max({scale, std::abs(jn), std::abs(jp)});
            }
            for (std::size_t i = 0; i < iso.node_count(); ++i) {
                for (std::size_t r = 1; r < 3; ++r) {
                    worst = std::max(worst, std::abs(f4[4 * i + r] - f3[3 * i + r]) / scale);
                }
            }
            CAPTURE(worst);
            REQUIRE(worst <= 1e-12);
        }
    }
}

TEST_CASE("electrothermal assemble: the Joule part of each edge's heat is not negative") {
    for (const Options o : {Options{.D = 1, .hetero = true}, Options{.D = 2, .nodes = 9}}) {
        const device::Device d = diode(o);
        const auto s = *assemble::make_scaling(d);
        auto dd = *DriftDiffusion::create(d, s, {.electrothermal = true});
        for (std::uint64_t seed = 1; seed <= 5; ++seed) {
            const auto x = heated(dd, d, s, 0.3, seed);
            for (const auto& [n, p] : dd.edge_joule_heat(x)) {
                REQUIRE(n >= 0.0);
                REQUIRE(p >= 0.0);
            }
        }
    }
}

TEST_CASE("electrothermal assemble: refused inputs") {
    const device::Device d = diode({});
    const auto s = *assemble::make_scaling(d);
    const auto refused = [&](const device::Device& dev, assemble::PhysicsModels models,
                             const std::string& text) {
        models.electrothermal = true;
        const auto dd = DriftDiffusion::create(dev, *assemble::make_scaling(dev), models);
        REQUIRE_FALSE(dd.has_value());
        REQUIRE(dd.error().code == base::ErrorCode::invalid_input);
        CAPTURE(dd.error().message);
        REQUIRE(dd.error().message.find(text) != std::string::npos);
    };
    // T8: band-to-band tunnelling and interface traps are refused, not frozen at T0.
    refused(d, {.btbt_local = true}, "band-to-band tunnelling (DECISIONS.md T8)");
    refused(diode({.D = 2, .nodes = 9, .oxide = true, .interface = true, .traps = true}), {},
            "interface traps (DECISIONS.md T8)");
    // T11: a material without thermal data.
    {
        physics::SemiconductorParameters p = physics::silicon_parameters;
        p.thermal = {};
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 1e-4, 5));
        auto left = m.find_boundary("x_min")->nodes;
        const auto bare = *device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"plain", *physics::Semiconductor::create(p)}},
             .node_region = std::vector<device::RegionId>(5, 0),
             .donors = std::vector<double>(5, 1e16),
             .acceptors = std::vector<double>(5, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)}}});
        refused(bare, {}, "region 'plain' has no thermal data");
    }
    // T3: Fermi-Dirac only with thermopower exponents of -1/2 (other_material's are not).
    refused(diode({.D = 1, .hetero = true}), {.fermi_dirac = true}, "exponents of -1/2");
    // Off, the same device and models are accepted.
    REQUIRE(DriftDiffusion::create(d, s, {.btbt_local = true}).has_value());
}
