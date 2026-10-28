// Schottky contacts in the drift-diffusion assembler (Unit 48; DECISIONS.md SD1-SD9, SC1, SC11,
// SC14): the finite-difference Jacobian gate on a 2D Delaunay and a 3D tensor mesh under both
// statistics, with and without Canali and band-gap narrowing; the pattern of an ohmic contact
// made Schottky unchanged; and the refusals of the device and of the assembler.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using assemble::EquilibriumPoisson;

namespace {

constexpr double L = 1e-4;  // [cm]

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

mesh::Mesh tensor(int D, int nodes) {
    const auto a = uniform(0.0, L, nodes);
    if (D == 1) return *mesh::make_tensor_grid(a);
    return D == 2 ? *mesh::make_tensor_grid(a, a) : *mesh::make_tensor_grid(a, a, a);
}

mesh::Mesh delaunay() {
    mesh::Geometry2D g{{{{{0, 0, 0}, {L, 0, 0}, {L, L, 0}, {0, L, 0}}, 0}},
                       {{"x_min", {0, 0, 0}, {0, L, 0}}, {"x_max", {L, 0, 0}, {L, L, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = L / 6.0;
    auto m = mesh::generate_mesh(g, o);
    REQUIRE(m.has_value());
    return std::move(m->mesh);
}

struct Diode {
    bool schottky = true;  // the anode on x_min; else ohmic
    device::SchottkyBarrier barrier{
        .boundary = "x_min", .barrier_eV = 0.7, .richardson_n = 110.0, .richardson_p = 30.0};
    bool sinks = false;  // the electrothermal model's thermal contacts
};

// An n-type silicon block (1e17 donors, 1e16 acceptors): the anode on x_min, an ohmic cathode on
// x_max.
std::expected<device::Device, base::Error> make_diode(mesh::Mesh m, const Diode& o = {}) {
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    std::vector<device::ThermalContact> thermal;
    if (o.sinks) {
        thermal = {{"sink", "x_max", device::ThermalContactKind::isothermal, 300.0, 0.0}};
    }
    device::Contact anode{"anode", o.schottky ? device::ContactKind::schottky
                                              : device::ContactKind::ohmic,
                          std::move(left)};
    anode.schottky = o.barrier;
    return device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(physics::silicon_parameters)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 1e17),
         .acceptors = std::vector<double>(n, 1e16),
         .contacts = {std::move(anode), {"cathode", device::ContactKind::ohmic, std::move(right)}},
         .thermal_contacts = std::move(thermal)});
}

device::Device diode(mesh::Mesh m, const Diode& o = {}) {
    auto d = make_diode(std::move(m), o);
    if (!d) FAIL(d.error().message);
    return std::move(*d);
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

// The probe state (as heterojunction_test.cpp): the charge-neutral guess with the anode at +0.3 V,
// psi + 0.02 noise, n and p times (1 + 0.01 noise), so n and p differ from n_B and p_B on the
// Schottky nodes.
std::vector<double> probe_state(DriftDiffusion& system, const device::Device& d,
                                const assemble::Scaling& s, const assemble::PhysicsModels& models,
                                std::uint64_t seed) {
    const std::vector<double> bias{0.3, 0.0};
    REQUIRE(system.set_bias(bias).has_value());
    const auto poisson = *EquilibriumPoisson::create(d, s, models);
    auto x = system.state_from_potential(poisson.charge_neutral_potential());
    Noise noise{seed};
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[3 * i] += 0.02 * noise.next();
        x[3 * i + 1] *= 1.0 + 0.01 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.01 * noise.next();
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

// The finite-difference gate (as surface_mobility_test.cpp): every column, central differences,
// the worst error over the column's largest entry, each column's the smallest at steps of 1e-5,
// 1e-6 and 1e-7 max(|x|, 1).
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

void require_invalid(const auto& made, const std::string& name) {
    CAPTURE(name);
    REQUIRE_FALSE(made.has_value());
    CAPTURE(made.error().message);
    REQUIRE(made.error().code == base::ErrorCode::invalid_input);
}

}  // namespace

TEST_CASE("Schottky assemble: the finite-difference Jacobian gate (SC11)") {
    for (const int D : {0, 3}) {  // 0: 2D Delaunay
        for (const bool fermi_dirac : {false, true}) {
            for (const bool canali : {false, true}) {
                for (const bool bgn : {false, true}) {
                    const device::Device d = diode(D == 0 ? delaunay() : tensor(3, 4));
                    const auto s = *assemble::make_scaling(d);
                    assemble::PhysicsModels models;
                    models.fermi_dirac = fermi_dirac;
                    models.field_mobility = canali;
                    models.bgn = bgn;
                    auto made = DriftDiffusion::create(d, s, models);
                    if (!made) FAIL(made.error().message);
                    const auto x = probe_state(*made, d, s, models, 7);
                    const double error = fd_jacobian_error(*made, x);
                    CAPTURE(D, fermi_dirac, canali, bgn, error);
                    REQUIRE(error <= 1e-6);
                }
            }
        }
    }
}

TEST_CASE("Schottky assemble: the rows a Schottky contact writes (SD3, SD4)") {
    // On a Schottky node the Poisson row is psi - psi_S alone and each continuity row keeps its
    // box terms plus the emission term; against the ohmic diode, whose rows on the contact node
    // are Dirichlet, the interior rows are the same.
    const device::Device ds = diode(tensor(1, 11));
    const device::Device dz = diode(tensor(1, 11), {.schottky = false});
    const auto s = *assemble::make_scaling(dz);
    const assemble::PhysicsModels models;
    auto sch = *DriftDiffusion::create(ds, s, models);
    auto ohm = *DriftDiffusion::create(dz, s, models);
    const auto x = probe_state(ohm, dz, s, models, 3);
    REQUIRE(sch.set_bias(std::vector<double>{0.3, 0.0}).has_value());
    std::vector<double> fs(sch.unknowns()), fz(ohm.unknowns());
    sch.residual(x, fs);
    ohm.residual(x, fz);
    for (std::size_t r = 3; r < fs.size(); ++r) {
        CAPTURE(r);
        REQUIRE(fs[r] == fz[r]);
    }
    // The node's Poisson row is psi - psi_S: a unit diagonal and nothing else; its continuity
    // rows read the neighbours (the box terms), unlike the ohmic node's Dirichlet rows.
    auto jacobian = sch.make_jacobian();
    std::vector<double> f(sch.unknowns());
    sch.evaluate(x, f, jacobian);
    const auto j = dense(jacobian);
    REQUIRE(j[0][0] == 1.0);
    for (std::size_t c = 1; c < j.size(); ++c) REQUIRE(j[0][c] == 0.0);
    REQUIRE(j[1][4] != 0.0);
    REQUIRE(j[2][5] != 0.0);
}

TEST_CASE("Schottky assemble: the pattern of an ohmic contact made Schottky (SC1)") {
    for (const int D : {1, 2, 3}) {
        const device::Device ds = diode(tensor(D, D == 3 ? 4 : 6));
        const device::Device dz = diode(tensor(D, D == 3 ? 4 : 6), {.schottky = false});
        const auto s = *assemble::make_scaling(dz);
        const auto sch = *DriftDiffusion::create(ds, s, assemble::PhysicsModels{});
        const auto ohm = *DriftDiffusion::create(dz, s, assemble::PhysicsModels{});
        CAPTURE(D);
        REQUIRE(sch.make_jacobian().has_same_pattern(ohm.make_jacobian()));
    }
}

TEST_CASE("Schottky assemble: a device whose only carrier contact is Schottky (SC13)") {
    // One contact: the Schottky anode fixes the potential and exchanges carriers (SD6), so the
    // topology check passes with no ohmic contact.
    mesh::Mesh m = tensor(1, 11);
    const std::size_t n = m.node_count();
    device::Contact anode{"anode", device::ContactKind::schottky, m.find_boundary("x_min")->nodes};
    anode.schottky = Diode{}.barrier;
    const auto made = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(physics::silicon_parameters)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 1e17),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {std::move(anode)}});
    if (!made) FAIL(made.error().message);
    const auto s = *assemble::make_scaling(*made);
    REQUIRE(DriftDiffusion::create(*made, s, {}).has_value());
}

TEST_CASE("Schottky assemble: refusals of the device and the assembler (SC14)") {
    // The device: a negative Richardson constant, a missing patch, a node off its patch, a
    // negative barrier, a work function not positive, an insulator node.
    {
        Diode o;
        o.barrier.richardson_n = -1.0;
        require_invalid(make_diode(tensor(1, 11), o), "negative Richardson (electrons)");
        o.barrier.richardson_n = 110.0;
        o.barrier.richardson_p = -1.0;
        require_invalid(make_diode(tensor(1, 11), o), "negative Richardson (holes)");
    }
    {
        Diode o;
        o.barrier.boundary = "nowhere";
        require_invalid(make_diode(tensor(1, 11), o), "missing patch");
        o.barrier.boundary = "x_max";
        require_invalid(make_diode(tensor(1, 11), o), "node off its patch");
    }
    {
        Diode o;
        o.barrier.barrier_eV = -0.1;
        require_invalid(make_diode(tensor(1, 11), o), "negative barrier");
        o.barrier.by_work_function = true;
        o.barrier.work_function_eV = 0.0;
        require_invalid(make_diode(tensor(1, 11), o), "work function not positive");
    }
    {
        // An oxide on x < 0.2 um: the Schottky node on x_min is an insulator node.
        mesh::Mesh m = tensor(1, 11);
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n, 1);
        std::vector<double> donors(n, 1e17);
        for (std::size_t i = 0; i < n; ++i) {
            if (m.points()[i][0] < 0.25 * L) {
                region[i] = 0;
                donors[i] = 0.0;
            }
        }
        device::Contact anode{"anode", device::ContactKind::schottky,
                              m.find_boundary("x_min")->nodes};
        anode.schottky = Diode{}.barrier;
        auto right = m.find_boundary("x_max")->nodes;
        const auto made = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
             .node_region = std::move(region),
             .donors = std::move(donors),
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {std::move(anode),
                          {"cathode", device::ContactKind::ohmic, std::move(right)}}});
        require_invalid(made, "insulator node");
    }
    // The assembler: Phi_Bn outside [0, Eg(T)] (by the barrier and by the work function), the
    // electrothermal model, a thermionic-emission interface.
    {
        Diode o;
        o.barrier.barrier_eV = 1.2;  // Eg(300 K) = 1.12 eV
        const device::Device d = diode(tensor(1, 11), o);
        require_invalid(DriftDiffusion::create(d, *assemble::make_scaling(d), {}), "Phi > Eg");
        o.barrier.by_work_function = true;
        o.barrier.work_function_eV = 3.9;  // chi = 4.05 eV: Phi_Bn = -0.15 eV
        const device::Device w = diode(tensor(1, 11), o);
        require_invalid(DriftDiffusion::create(w, *assemble::make_scaling(w), {}), "Phi < 0");
        o.barrier.work_function_eV = 4.75;
        const device::Device ok = diode(tensor(1, 11), o);
        REQUIRE(DriftDiffusion::create(ok, *assemble::make_scaling(ok), {}).has_value());
    }
    {
        const device::Device d = diode(tensor(1, 11), {.sinks = true});
        assemble::PhysicsModels models;
        models.electrothermal = true;
        require_invalid(DriftDiffusion::create(d, *assemble::make_scaling(d), models),
                        "electrothermal");
        // Isothermal, the same device assembles.
        REQUIRE(DriftDiffusion::create(d, *assemble::make_scaling(d), {}).has_value());
    }
    {
        // Silicon | silicon with chi + 0.2 eV at x = 0.5 um, the interface thermionic.
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, L, 21));
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n, 0);
        for (std::size_t i = 0; i < n; ++i) region[i] = m.points()[i][0] < 0.5 * L ? 0 : 1;
        physics::SemiconductorParameters shifted = physics::silicon_parameters;
        shifted.electron_affinity_eV += 0.2;
        device::Contact anode{"anode", device::ContactKind::schottky,
                              m.find_boundary("x_min")->nodes};
        anode.schottky = Diode{}.barrier;
        auto right = m.find_boundary("x_max")->nodes;
        const auto made = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"left", *physics::Semiconductor::create(physics::silicon_parameters)},
                         {"right", *physics::Semiconductor::create(shifted)}},
             .node_region = std::move(region),
             .donors = std::vector<double>(n, 1e17),
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {std::move(anode),
                          {"cathode", device::ContactKind::ohmic, std::move(right)}},
             .interfaces = {{"left", "right", device::InterfaceTransport::thermionic_emission}}});
        if (!made) FAIL(made.error().message);
        require_invalid(DriftDiffusion::create(*made, *assemble::make_scaling(*made), {}),
                        "thermionic interface");
    }
}

TEST_CASE("Schottky assemble: thermal equilibrium needs a Schottky contact at 0 V (SD6)") {
    const std::vector<device::ContactKind> kinds{device::ContactKind::schottky,
                                                 device::ContactKind::ohmic};
    REQUIRE(assemble::check_contact_bias(kinds, std::vector<double>{0.0, 0.0}, true).has_value());
    REQUIRE(assemble::check_contact_bias(kinds, std::vector<double>{0.2, 0.0}, false).has_value());
    const auto refused = assemble::check_contact_bias(kinds, std::vector<double>{0.2, 0.0}, true);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().code == base::ErrorCode::invalid_input);
    CHECK(refused.error().message.find("Schottky contact at 0 V") != std::string::npos);
}
