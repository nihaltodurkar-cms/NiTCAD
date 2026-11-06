// The manufactured-solution hook (verification only; assemble/models.hpp ManufacturedSources;
// DECISIONS.md Unit 55 14.3, TVQ3, CV-Q4): absent, the composition and every residual and Jacobian
// value are those of the system without it; present, it adds V S / Ns to the Poisson rows,
// +V G_n / R0 to the electron rows and -V G_p / R0 to the hole rows of the nodes off the boundary
// rows, nothing to the Jacobian, and is refused where it is not defined.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "terms.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using assemble::ManufacturedSources;
using assemble::PhysicsModels;

namespace {

// 1D: oxide [-5e-7, 0] with an electrode at x_min, p-silicon [0, 2e-5] with an ohmic contact at
// x_max; 2D: the same extruded.
device::Device mos(int dimension) {
    std::vector<double> x;
    for (int i = 0; i <= 40; ++i) x.push_back(-5e-7 + (2e-5 + 5e-7) * i / 40.0);
    for (const double y : {-4e-7, -3e-7, -2e-7, -1e-7, 0.0}) x.push_back(y);
    std::ranges::sort(x);
    std::vector<double> t{0.0, 1e-6, 2e-6};
    mesh::Mesh m = dimension == 1 ? *mesh::make_tensor_grid(x) : *mesh::make_tensor_grid(x, t);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const bool oxide = m.points()[i][0] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (!oxide) acceptors[i] = 1e17;
    }
    device::Contact gate{"gate", device::ContactKind::electrode, m.find_boundary("x_min")->nodes};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .interfaces = {{"oxide", "silicon"}}});
}

// A state off equilibrium: psi a ramp, densities perturbed, the boundary rows stamped.
std::vector<double> probe(const DriftDiffusion& dd, std::size_t nodes) {
    std::vector<double> psi(nodes);
    for (std::size_t i = 0; i < nodes; ++i) psi[i] = -3.0 + 0.11 * static_cast<double>(i % 41);
    std::vector<double> x = dd.state_from_potential(psi);
    for (std::size_t i = 0; i < nodes; ++i) {
        x[3 * i + 1] *= 1.0 + 0.01 * static_cast<double>(i % 7);
        x[3 * i + 2] *= 1.0 - 0.005 * static_cast<double>(i % 5);
    }
    dd.stamp_contacts(x);
    return x;
}

struct Evaluated {
    std::vector<double> f;
    std::vector<double> values;
};

Evaluated evaluate(const DriftDiffusion& dd, const std::vector<double>& x) {
    Evaluated e{std::vector<double>(dd.unknowns()), {}};
    linalg::SparseMatrix j = dd.make_jacobian();
    dd.evaluate(x, e.f, j);
    e.values.assign(j.values().begin(), j.values().end());
    return e;
}

std::vector<std::string> names(const DriftDiffusion& dd) {
    std::vector<std::string> out;
    for (const auto& t : assemble::TermAccess::Terms::list(dd)) out.emplace_back(t->info().name);
    return out;
}

std::shared_ptr<const ManufacturedSources> sources(std::size_t nodes, double S, double Gn, double Gp) {
    auto ms = std::make_shared<ManufacturedSources>();
    for (std::size_t i = 0; i < nodes; ++i) {
        const double w = 1.0 + 0.25 * static_cast<double>(i % 3);
        ms->charge_cm3.push_back(S * w);
        ms->electron_generation.push_back(Gn * w);
        ms->hole_generation.push_back(Gp * w);
    }
    return ms;
}

// The term's own entry against V S / Ns or V G / R0 computed here: the two differ only by the
// rounding of the volume's and the scales' products (a few ulp).
bool close(double value, double expected) {
    return std::abs(value - expected) <= 1e-14 * std::abs(expected);
}

}  // namespace

TEST_CASE("manufactured sources: absent, nothing changes (TVQ3)", "[manufactured]") {
    for (const int D : {1, 2}) {
        const device::Device d = mos(D);
        const auto s = *assemble::make_scaling(d);
        const auto plain = *DriftDiffusion::create(d, s);
        PhysicsModels null_models;
        null_models.manufactured = nullptr;
        const auto with_null = *DriftDiffusion::create(d, s, null_models);
        CHECK(names(plain) == names(with_null));
        for (const auto& n : names(plain)) CHECK(n != "manufactured source");
        const auto x = probe(plain, d.mesh().node_count());
        const Evaluated a = evaluate(plain, x), b = evaluate(with_null, x);
        CHECK(a.f == b.f);
        CHECK(a.values == b.values);
    }
}

TEST_CASE("manufactured sources: the rows they enter, their scaling, no Jacobian (TVQ3)", "[manufactured]") {
    for (const int D : {1, 2}) {
        CAPTURE(D);
        const device::Device d = mos(D);
        const std::size_t nodes = d.mesh().node_count();
        const auto s = *assemble::make_scaling(d);
        const auto plain = *DriftDiffusion::create(d, s);
        PhysicsModels models;
        const double S = 3e15, Gn = 2e20, Gp = 5e19;
        const auto owned = sources(nodes, S, Gn, Gp);  // the system does not own its sources
        models.manufactured = owned.get();
        const auto dd = *DriftDiffusion::create(d, s, models);
        const auto list = names(dd);
        // One node term, after storage, before the edge terms.
        const auto at = std::ranges::find(list, "manufactured source");
        REQUIRE(at != list.end());
        REQUIRE(at != list.begin());
        CHECK(*(at - 1) == "storage");
        const auto x = probe(plain, nodes);
        const Evaluated a = evaluate(plain, x), b = evaluate(dd, x);
        CHECK(a.values == b.values);  // independent of the state: no partials
        // The term alone (as terms_test evaluates one term): its rows and nothing in the Jacobian.
        const auto& terms = assemble::TermAccess::Terms::list(dd);
        const auto term = std::ranges::find_if(
            terms, [](const auto& t) { return t->info().name == "manufactured source"; });
        REQUIRE(term != terms.end());
        std::vector<double> f(dd.unknowns());
        auto jacobian = dd.make_jacobian();
        assemble::TermAccess::Terms::evaluate(dd, **term, x, nullptr, f, jacobian);
        CHECK(std::ranges::all_of(jacobian.values(), [](double v) { return v == 0.0; }));
        const double LD = std::pow(s.L_D, D);
        std::vector<char> on_contact(nodes, 0);
        for (const auto& c : d.contacts()) {
            for (const mesh::NodeId v : c.nodes) on_contact[static_cast<std::size_t>(v)] = 1;
        }
        std::size_t poisson = 0, carriers = 0;
        for (std::size_t i = 0; i < nodes; ++i) {
            CAPTURE(i);
            const double V = d.mesh().volumes()[i] / LD;
            const ManufacturedSources& ms = *models.manufactured;
            const bool contact = on_contact[i] != 0;
            const bool insulator = d.is_insulator(static_cast<mesh::NodeId>(i));
            if (contact) {  // ohmic and electrode rows are boundary rows
                CHECK(f[3 * i] == 0.0);
                CHECK(f[3 * i + 1] == 0.0);
                CHECK(f[3 * i + 2] == 0.0);
                continue;
            }
            CHECK(close(f[3 * i], V * ms.charge_cm3[i] / s.Ns));
            ++poisson;
            if (insulator) {
                CHECK(f[3 * i + 1] == 0.0);
                CHECK(f[3 * i + 2] == 0.0);
                continue;
            }
            CHECK(close(f[3 * i + 1], V * ms.electron_generation[i] / s.R0));
            CHECK(close(f[3 * i + 2], -V * ms.hole_generation[i] / s.R0));
            ++carriers;
        }
        CHECK(poisson > carriers);  // the oxide's interior nodes took charge
        CHECK(carriers > 0);
    }
}

TEST_CASE("manufactured sources: an empty vector leaves its rows alone (TVQ3)", "[manufactured]") {
    const device::Device d = mos(1);
    const std::size_t nodes = d.mesh().node_count();
    const auto s = *assemble::make_scaling(d);
    const auto plain = *DriftDiffusion::create(d, s);
    PhysicsModels models;
    auto ms = std::make_shared<ManufacturedSources>();
    ms->electron_generation.assign(nodes, 1e20);
    models.manufactured = ms.get();
    const auto dd = *DriftDiffusion::create(d, s, models);
    const auto x = probe(plain, nodes);
    const Evaluated a = evaluate(plain, x), b = evaluate(dd, x);
    for (std::size_t i = 0; i < nodes; ++i) {
        CHECK(b.f[3 * i] == a.f[3 * i]);
        CHECK(b.f[3 * i + 2] == a.f[3 * i + 2]);
    }
}

TEST_CASE("manufactured sources: refused where not defined (TVQ3)", "[manufactured]") {
    const device::Device d = mos(1);
    const std::size_t nodes = d.mesh().node_count();
    const auto s = *assemble::make_scaling(d);
    const auto refused = [&](PhysicsModels m) { return !DriftDiffusion::create(d, s, m).has_value(); };
    PhysicsModels m;
    auto wrong = std::make_shared<ManufacturedSources>();
    wrong->charge_cm3.assign(nodes - 1, 0.0);
    m.manufactured = wrong.get();
    CHECK(refused(m));
    auto nan = std::make_shared<ManufacturedSources>();
    nan->hole_generation.assign(nodes, 0.0);
    nan->hole_generation[3] = std::numeric_limits<double>::quiet_NaN();
    m.manufactured = nan.get();
    CHECK(refused(m));
    const auto owned = sources(nodes, 1.0, 1.0, 1.0);
    m.manufactured = owned.get();
    m.electrothermal = true;
    CHECK(refused(m));
    m.electrothermal = false;
    m.energy_balance = true;
    CHECK(refused(m));
}
