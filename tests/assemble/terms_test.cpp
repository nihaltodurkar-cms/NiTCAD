// The equation terms of the drift-diffusion rows on their own (Unit 40, DECISIONS.md FE2, FG3):
// each term's partials against central differences of its own residual, at non-equilibrium states
// where it is active, and the terms' declared order and domains (FE3).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "terms.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using assemble::EquilibriumPoisson;
using assemble::TermDomain;
using Terms = DriftDiffusion::Terms;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
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

// One term's finite-difference check, as the Unit 19 gate: every column, central differences of
// the term's own residual, the worst error over the largest entry of the column in the whole
// system's Jacobian (the existing gates' scale: alone, a term can hold a column whose only entry
// is far below what its residual resolves, such as the vanishing side of a Scharfetter-Gummel flux
// at a drop of hundreds of V_T, 1e-128 beside a flux of 300), each column's error the smallest at
// steps of 1e-5, 1e-6 and 1e-7 max(|x|, 1). And whether the term is active at x: a nonzero
// residual entry or partial.
struct Check {
    double error;
    bool active;
};

Check check(const DriftDiffusion& dd, const DriftDiffusion::Term& t, std::vector<double> x,
            const DriftDiffusion::TimeStep* step = nullptr) {
    const std::size_t n = dd.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    auto jacobian = dd.make_jacobian();
    Terms::evaluate(dd, t, x, step, f, jacobian);
    const bool active =
        std::ranges::any_of(f, [](double v) { return v != 0.0; }) ||
        std::ranges::any_of(jacobian.values(), [](double v) { return v != 0.0; });
    const auto j = dense(jacobian);
    auto scratch = dd.make_jacobian();
    if (step != nullptr) {
        dd.evaluate(x, *step, fp, scratch);
    } else {
        dd.evaluate(x, fp, scratch);
    }
    const auto whole = dense(scratch);
    double worst = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        double scale = 1e-300;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(whole[r][c]));
        double best = std::numeric_limits<double>::infinity();
        for (const double relative : {1e-5, 1e-6, 1e-7}) {
            const double base = x[c];
            const double h = relative * std::max(std::abs(base), 1.0);
            x[c] = base + h;
            const double up = x[c];
            Terms::evaluate(dd, t, x, step, fp, scratch);
            x[c] = base - h;
            const double down = x[c];
            Terms::evaluate(dd, t, x, step, fm, scratch);
            x[c] = base;
            double error = 0.0;
            for (std::size_t r = 0; r < n; ++r) {
                error = std::max(error, std::abs((fp[r] - fm[r]) / (up - down) - j[r][c]));
            }
            best = std::min(best, error / scale);
        }
        worst = std::max(worst, best);
    }
    return {worst, active};
}

// Every term of a system at x: its check within the tolerance; the names of the active ones are
// added to `active`.
void check_terms(const DriftDiffusion& dd, const std::vector<double>& x,
                 std::set<std::string>& active, const DriftDiffusion::TimeStep* step = nullptr,
                 double tolerance = 1e-6) {
    for (const auto& t : Terms::list(dd)) {
        const Check c = check(dd, *t, x, step);
        const std::string name(t->info().name);
        CAPTURE(name, c.error, c.active);
        CHECK(c.error <= tolerance);
        if (c.active) active.insert(name);
    }
}

// Every term of FE2 is active in some section (FR9): a term at zero everywhere checks nothing.
void require_active(const std::set<std::string>& active,
                    std::initializer_list<const char*> names) {
    for (const char* name : names) {
        CAPTURE(name);
        CHECK(active.contains(name));
    }
}

// Silicon with a radiative coefficient (silicon's set leaves it 0).
physics::Semiconductor radiative_silicon() {
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.radiative_cm3_s = 1e-14;
    return *physics::Semiconductor::create(p);
}

// A pn diode (1e18 / 1e18) on a graded 1D or 2D grid, ohmic contacts on x_min and x_max.
device::Device diode(int D) {
    const auto x = uniform(0.0, 2e-4, 41);
    const auto t = uniform(0.0, 1e-4, 4);
    mesh::Mesh m = D == 1 ? *mesh::make_tensor_grid(x) : *mesh::make_tensor_grid(x, t);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e18;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", radiative_silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

// The Unit 15b fixture (transient_test.cpp): oxide on [-5 nm, 0) with a metal electrode on
// x_min (in 2D on part of it), p-silicon with an n+ pocket (2D) and the substrate contact on
// x_max; a fixed charge, two trap levels, a trap band and surface recombination at the interface.
device::Device mos(int D) {
    const auto x = *mesh::straddle_interface(uniform(-5e-7, 1.2e-5, 30), 0.0, 1e-7);
    const auto t = uniform(0.0, 3e-6, 4);
    mesh::Mesh m = D == 1 ? *mesh::make_tensor_grid(x) : *mesh::make_tensor_grid(x, t);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        const bool oxide = p[0] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (oxide) continue;
        acceptors[i] = 1e17;
        if (D > 1 && p[1] > 2e-6 && p[0] < 3e-6) donors[i] = 1e19;
    }
    std::vector<mesh::NodeId> electrode;
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        if (D == 1 || m.points()[static_cast<std::size_t>(v)][1] < 1.5e-6) electrode.push_back(v);
    }
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = 3e11;
    f.traps.levels = {{.type = physics::TrapType::donor, .density_cm2 = 2e11, .energy_eV = -0.2,
                       .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 3e-16},
                      {.type = physics::TrapType::acceptor, .density_cm2 = 1e11,
                       .energy_eV = 0.25, .sigma_n_cm2 = 2e-16, .sigma_p_cm2 = 1e-15}};
    f.traps.bands = {{.type = physics::TrapType::acceptor, .density_cm2_eV = 5e11,
                      .energy_low_eV = -0.1, .energy_high_eV = 0.4}};
    f.recombination_velocity_n_cm_s = 1e3;
    f.recombination_velocity_p_cm_s = 4e2;
    device::Contact gate{"gate", device::ContactKind::electrode, std::move(electrode)};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .interfaces = {std::move(f)}});
}

// A p-type MOS-C with a lumped gate (gate_test.cpp): gate on x_min, substrate on x_max.
device::Device moscap() {
    mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 2e-5, 30));
    const std::size_t n = m.node_count();
    device::Contact g{"gate", device::ContactKind::gate, m.find_boundary("x_min")->nodes};
    g.gate = {.boundary = "x_min", .oxide_thickness_cm = 5e-7, .fixed_charge_cm2 = 2e11};
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::vector<double>(n, 1e17),
         .contacts = {std::move(g),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
}

// A uniformly doped block (1e16 donors) with its cells, ohmic contacts on x_min and x_max
// (impact_ionization_test.cpp, band_to_band_test.cpp).
device::Device block(int D, int nodes, double length_cm) {
    const auto a = uniform(0.0, length_cm, nodes);
    mesh::Mesh m = D == 1 ? *mesh::make_tensor_grid(a) : *mesh::make_tensor_grid(a, a);
    auto cells = D == 1 ? *mesh::TensorCells::create(m, a) : *mesh::TensorCells::create(m, a, a);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 1e16),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}},
         .cells = std::move(cells)});
}

// The charge-neutral state under bias with noise on every unknown (psi by 0.02 V_T, densities by
// 1%), the insulator nodes' densities left at 0 and, near an interface (0 < x < 1e-7 cm), densities
// away from equilibrium so the traps are resolved (transient_test.cpp).
std::vector<double> probe(const DriftDiffusion& dd, const device::Device& d, std::uint64_t seed) {
    const auto scaling = *assemble::make_scaling(d);
    const auto poisson = *EquilibriumPoisson::create(d, scaling);
    auto x = dd.state_from_potential(poisson.charge_neutral_potential());
    Noise noise{seed};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        x[3 * i] += 0.02 * noise.next();
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        x[3 * i + 1] *= 1.0 + 0.01 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.01 * noise.next();
        const double px = d.mesh().points()[i][0];
        if (!d.interfaces().empty() && px > 0.0 && px < 1e-7) {
            x[3 * i + 1] = 2e-3 * (1.0 + 0.01 * noise.next());
            x[3 * i + 2] = 5e-4 * (1.0 + 0.01 * noise.next());
        }
    }
    dd.stamp_contacts(x);
    return x;
}

// A field E0 [V/cm] at angle theta to x, n and p uniform, relative noise on every unknown
// (impact_ionization_test.cpp).
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

// A time step's history near x: storage with 5% noise on the densities, trap occupancies spread
// over [-0.05, 1.05] (transient_test.cpp).
struct History {
    std::vector<double> storage, traps;
};

History history(const DriftDiffusion& dd, std::vector<double> x, std::uint64_t seed) {
    Noise noise{seed};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        x[3 * i + 1] *= 1.0 + 0.05 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.05 * noise.next();
    }
    History h{dd.storage(x), std::vector<double>(dd.trap_slots())};
    for (double& f : h.traps) f = 0.5 + 0.55 * noise.next();
    return h;
}

}  // namespace

TEST_CASE("terms: each term's partials match finite differences of its own residual") {
    std::set<std::string> active;
    SECTION("a diode, every node model, steady and a time step") {
        for (const int D : {1, 2}) {
            const device::Device d = diode(D);
            const auto s = *assemble::make_scaling(d);
            for (const bool fermi_dirac : {false, true}) {
                const assemble::PhysicsModels models{.field_mobility = true,
                                                     .fermi_dirac = fermi_dirac,
                                                     .incomplete_ionization = true};
                auto dd = *DriftDiffusion::create(d, s, models);
                REQUIRE(dd.set_bias(std::vector<double>{0.3, 0.0}).has_value());
                const auto x = probe(dd, d, 5);
                CAPTURE(D, fermi_dirac);
                check_terms(dd, x, active);
                const History h = history(dd, x, 11);
                const DriftDiffusion::TimeStep step{1e-3, h.storage, h.traps};
                check_terms(dd, x, active, &step);
            }
        }
        require_active(active, {"ohmic contact rows", "Poisson charge", "incomplete ionization",
                                "SRH recombination", "Auger recombination",
                                "radiative recombination", "storage", "Poisson flux",
                                "carrier fluxes"});
    }
    SECTION("an oxide with an electrode and a trapped interface, steady and a time step") {
        for (const int D : {1, 2}) {
            const device::Device d = mos(D);
            const auto s = *assemble::make_scaling(d);
            auto dd = *DriftDiffusion::create(d, s);
            REQUIRE(dd.set_bias(std::vector<double>{-0.4, 0.1}).has_value());
            const auto x = probe(dd, d, 3);
            CAPTURE(D);
            check_terms(dd, x, active);
            const History h = history(dd, x, 13);
            const DriftDiffusion::TimeStep step{1e-3, h.storage, h.traps};
            check_terms(dd, x, active, &step);
        }
        require_active(active, {"insulator rows", "interfaces", "storage"});
    }
    SECTION("a lumped gate") {
        const device::Device d = moscap();
        const auto s = *assemble::make_scaling(d);
        auto dd = *DriftDiffusion::create(d, s);
        REQUIRE(dd.set_bias(std::vector<double>{0.5, 0.0}).has_value());
        check_terms(dd, probe(dd, d, 9), active);
        require_active(active, {"gate"});
    }
    SECTION("impact ionization at 3e5 V/cm") {
        for (const int D : {1, 2}) {
            const device::Device d = block(D, 11, 1e-4);
            const auto s = *assemble::make_scaling(d);
            const auto dd = *DriftDiffusion::create(d, s, {.impact_ionization = true});
            CAPTURE(D);
            check_terms(dd, sloped(dd, d, s, 3e5, D == 1 ? 0.0 : std::numbers::pi / 6.0, 1e16,
                                   1e15, 0.05, 7),
                        active);
        }
        require_active(active, {"impact ionization"});
    }
    SECTION("local and nonlocal band-to-band tunnelling at 4e6 V/cm") {
        for (const int D : {1, 2}) {
            const device::Device d = block(D, 13, 3e-6);
            const auto s = *assemble::make_scaling(d);
            for (const assemble::PhysicsModels models :
                 {assemble::PhysicsModels{.btbt_local = true},
                  assemble::PhysicsModels{.btbt_nonlocal = assemble::NonlocalTunnelling::kane}}) {
                auto dd = *DriftDiffusion::create(d, s, models);
                const auto x = sloped(dd, d, s, 4e6, D == 1 ? 0.0 : std::numbers::pi / 9.0,
                                      1e17, 1e3, 0.02, 3);
                if (dd.tunnelling()) {
                    dd.set_paths(dd.trace_paths(x));
                    REQUIRE_FALSE(dd.paths().paths.empty());
                }
                CAPTURE(D, models.btbt_local);
                check_terms(dd, x, active);
            }
        }
        require_active(active,
                       {"local band-to-band tunnelling", "nonlocal band-to-band tunnelling"});
    }
}

TEST_CASE("terms: the assembly order and domains of FE3") {
    // Every model that adds a term on, on a device with an interface: the full list in order.
    const device::Device d = mos(1);
    const auto s = *assemble::make_scaling(d);
    const auto dd = *DriftDiffusion::create(
        d, s, {.radiative = true, .incomplete_ionization = true, .impact_ionization = true,
               .btbt_local = true});
    std::vector<std::string> names;
    std::vector<TermDomain> domains;
    for (const auto& t : Terms::list(dd)) {
        names.emplace_back(t->info().name);
        domains.push_back(t->info().domain);
    }
    CHECK(names == std::vector<std::string>{"ohmic contact rows", "insulator rows",
                                            "Poisson charge", "incomplete ionization",
                                            "SRH recombination", "Auger recombination",
                                            "radiative recombination", "storage", "interfaces",
                                            "Poisson flux", "carrier fluxes", "impact ionization",
                                            "local band-to-band tunnelling"});
    // Node terms first, then the interface edges, the mesh edges, the global terms.
    CHECK(std::ranges::is_sorted(domains));
    // The assignments (DECISIONS.md FE3, FR3): the rows set rather than added to.
    for (const auto& t : Terms::list(dd)) {
        const std::string name(t->info().name);
        CAPTURE(name);
        CHECK(t->info().assigns == (name == "ohmic contact rows" || name == "insulator rows" ||
                                    name == "Poisson charge" || name == "incomplete ionization"));
    }
}
