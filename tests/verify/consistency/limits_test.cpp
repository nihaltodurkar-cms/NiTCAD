// C10, limit continuity (DECISIONS Unit 55 section 11.6 and the C10 ruling of 11.13): the difference to the limit
// model decays at the rate derived before measurement, on three 1D devices written here: energy balance -> drift-
// diffusion as tau_w -> 0 with v_sat -> infinity (EC4's limit; a terminal current, rate 1 in tau_w); Fermi-Dirac ->
// Boltzmann as the doping of a symmetric pn diode -> 0 (the potential, rate 1 in N; the limit model's device moves
// with N, so one Boltzmann result per value); electrothermal -> isothermal as kappa -> infinity with isothermal
// sinks at the device temperature (a terminal current, rate 1 in 1/kappa). A wrong stated rate and a corrupted
// point fail it; fewer than three values are NOT APPLICABLE; a difference at its floor is NOT EVALUATED.
#include "NiTCAD/verify/consistency/checks.hpp"

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace c = NiTCAD::verify::consistency;
namespace v = NiTCAD::verify;
namespace solve = NiTCAD::solve;
namespace device = NiTCAD::device;
namespace mesh = NiTCAD::mesh;
namespace physics = NiTCAD::physics;
namespace results = NiTCAD::results;
using v::Outcome;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> x(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) x[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return x;
}

// A 1D silicon pn diode on [0, L], acceptors N_A on the left half and donors N_D on the right; ohmic contacts
// left (x_min) and right (x_max); optional isothermal sinks at both ends.
device::Device diode(double L, int nodes, double N_A, double N_D, physics::SemiconductorParameters material,
                     bool sinks = false, double T = 300.0) {
    mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, L, nodes));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 0.5 * L ? acceptors[i] = N_A : donors[i] = N_D);
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    std::vector<device::ThermalContact> thermal;
    if (sinks) {
        thermal = {{"left", "x_min", device::ThermalContactKind::isothermal, T, 0.0},
                   {"right", "x_max", device::ThermalContactKind::isothermal, T, 0.0}};
    }
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = T,
                                     .regions = {{"silicon", *physics::Semiconductor::create(material)}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                  {"right", device::ContactKind::ohmic, std::move(right)}},
                                     .thermal_contacts = std::move(thermal)});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

results::BiasPoint solved(const device::Device& d, std::vector<double> bias, const solve::BiasOptions& o) {
    auto p = solve::solve_bias(d, bias, o);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

// The model and limit results of one demonstration, kept alive for the series' pointers.
struct Study {
    std::optional<device::Device> device;  // the model's device at the first value (for the mesh and data)
    std::vector<results::BiasPoint> model, limit;
    c::DeviceData data;
    c::LimitSeries series(std::vector<double> x, c::LimitQuantity q, std::size_t contact, double rate,
                          std::string parameter) const {
        c::LimitSeries s;
        s.parameter = std::move(parameter);
        s.x = std::move(x);
        for (const auto& p : model) s.model.push_back(&p);
        for (const auto& p : limit) s.limit.push_back(&p);
        s.quantity = q;
        s.contact = contact;
        s.stated_rate = rate;
        return s;
    }
};

const std::vector<double> tau_values{1e-14, 1e-15, 1e-16};

// (i) Energy balance with v_sat -> infinity against drift-diffusion, a forward diode (1e17 / 2e17, 1 um) at 0.6 V.
const Study& energy_balance() {
    static const std::optional<Study> s = []() -> std::optional<Study> {
        Study out;
        solve::BiasOptions eb;
        eb.models.energy_balance = true;
        eb.step_control.max_cuts = 20;
        for (const double tau : tau_values) {
            physics::SemiconductorParameters p = physics::silicon_parameters;
            p.energy_relaxation = {.tau_n_s = tau, .tau_p_s = tau};
            p.electron_saturation.v_sat_cm_s = 1e30;
            p.hole_saturation.v_sat_cm_s = 1e30;
            device::Device d = diode(1e-4, 101, 1e17, 2e17, p);
            out.model.push_back(solved(d, {0.6, 0.0}, eb));
            if (!out.device) out.device.emplace(std::move(d));
        }
        const device::Device dd = diode(1e-4, 101, 1e17, 2e17, physics::silicon_parameters);
        out.limit.push_back(solved(dd, {0.6, 0.0}, {}));
        out.data = NiTCAD::verify::benchmark::device_data(*out.device);
        return out;
    }();
    REQUIRE(s.has_value());
    return *s;
}

const std::vector<double> doping_values{1e15, 1e16, 1e17};

// (ii) Fermi-Dirac against Boltzmann at 0 V, a symmetric diode N / N, 10 um (wide enough for 1e15's depletion).
const Study& fermi_dirac() {
    static const std::optional<Study> s = []() -> std::optional<Study> {
        Study out;
        solve::BiasOptions fd;
        fd.models.fermi_dirac = true;
        for (const double N : doping_values) {
            device::Device d = diode(1e-3, 201, N, N, physics::silicon_parameters);
            out.model.push_back(solved(d, {0.0, 0.0}, fd));
            out.limit.push_back(solved(d, {0.0, 0.0}, {}));
            if (!out.device) out.device.emplace(std::move(d));
        }
        out.data = NiTCAD::verify::benchmark::device_data(*out.device);
        return out;
    }();
    REQUIRE(s.has_value());
    return *s;
}

const std::vector<double> kappa_scales{1.0, 10.0, 100.0};

// (iii) Electrothermal with kappa scaled by s against the isothermal device, a forward diode (1e17 / 1e17, 1 um,
// isothermal sinks at 300 K at both ends) at 0.7 V; x = 1 / s.
const Study& electrothermal() {
    static const std::optional<Study> s = []() -> std::optional<Study> {
        Study out;
        solve::BiasOptions et;
        et.models.electrothermal = true;
        et.newton.tol_update = 1e-12;
        for (const double scale : kappa_scales) {
            physics::SemiconductorParameters p = physics::silicon_parameters;
            p.thermal.conductivity_W_cmK *= scale;
            device::Device d = diode(1e-4, 41, 1e17, 1e17, p, true);
            out.model.push_back(solved(d, {0.7, 0.0}, et));
            if (!out.device) out.device.emplace(std::move(d));
        }
        solve::BiasOptions plain;
        plain.newton.tol_update = 1e-12;
        const device::Device iso = diode(1e-4, 41, 1e17, 1e17, physics::silicon_parameters);
        out.limit.push_back(solved(iso, {0.7, 0.0}, plain));
        out.data = NiTCAD::verify::benchmark::device_data(*out.device, et.models);
        return out;
    }();
    REQUIRE(s.has_value());
    return *s;
}

std::vector<double> reciprocal(const std::vector<double>& s) {
    std::vector<double> x;
    for (const double v : s) x.push_back(1.0 / v);
    return x;
}

void report(const char* name, const c::CheckResult& r) {
    for (const auto& d : r.diagnostics) {
        std::printf("C10 %s: %s %s = %.4e (floor %.3e)\n", name, d.name.c_str(), d.where.c_str(), d.value, d.scale);
    }
    for (const auto& p : r.parts) {
        std::printf("C10 %s: part %s worst %.4f (tol %.2f) at %s\n", name, p.name.c_str(), p.residual, p.tolerance,
                    p.where.c_str());
    }
}

}  // namespace

TEST_CASE("C10: energy balance -> drift-diffusion at rate 1 in tau_w", "[verify][consistency]") {
    const Study& s = energy_balance();
    const auto r = c::check_c10(s.series(tau_values, c::LimitQuantity::terminal_current, 0, 1.0, "tau_w [s]"),
                                s.device->mesh(), s.data);
    report("EB", r);
    CAPTURE(r.comparison.reason);
    REQUIRE(r.id == "C10");
    REQUIRE(r.comparison.outcome == Outcome::pass);
}

TEST_CASE("C10: Fermi-Dirac -> Boltzmann at rate 1 in N", "[verify][consistency]") {
    const Study& s = fermi_dirac();
    const auto r = c::check_c10(s.series(doping_values, c::LimitQuantity::potential, 0, 1.0, "N [cm^-3]"),
                                s.device->mesh(), s.data);
    report("FD", r);
    CAPTURE(r.comparison.reason);
    REQUIRE(r.comparison.outcome == Outcome::pass);
}

TEST_CASE("C10: electrothermal -> isothermal at rate 1 in 1/kappa", "[verify][consistency]") {
    const Study& s = electrothermal();
    const auto r = c::check_c10(s.series(reciprocal(kappa_scales), c::LimitQuantity::terminal_current, 0, 1.0,
                                         "1 / kappa scale"),
                                s.device->mesh(), s.data);
    report("ET", r);
    CAPTURE(r.comparison.reason);
    REQUIRE(r.comparison.outcome == Outcome::pass);
}

TEST_CASE("C10 fails on a wrong rate or a corrupted point", "[verify][consistency]") {
    const Study& s = energy_balance();
    const auto& m = s.device->mesh();
    const auto quantity = c::LimitQuantity::terminal_current;
    // The stated rate 0.75 and 1.25 (just outside 0.2 of the derived 1): fail.
    for (const double rate : {0.75, 1.25}) {
        const auto r = c::check_c10(s.series(tau_values, quantity, 0, rate, "tau_w [s]"), m, s.data);
        CAPTURE(rate);
        REQUIRE(r.comparison.outcome == Outcome::fail);
    }
    // The middle point's current moved so its difference doubles: the two rates become 0.7 and 1.3.
    Study bad = s;
    const double I_lim = s.limit[0].terminal_current[0];
    bad.model[1].terminal_current[0] = I_lim + 2.0 * (s.model[1].terminal_current[0] - I_lim);
    const auto r = c::check_c10(bad.series(tau_values, quantity, 0, 1.0, "tau_w [s]"), m, s.data);
    REQUIRE(r.comparison.outcome == Outcome::fail);
    REQUIRE(r.parts.front().where.find("1e-15") != std::string::npos);
    // The largest value's difference doubled: only the last rate (1e-15 .. 1e-14) moves, to about 1.23.
    Study last = s;
    last.model[0].terminal_current[0] = I_lim + 2.0 * (s.model[0].terminal_current[0] - I_lim);
    const auto rl = c::check_c10(last.series(tau_values, quantity, 0, 1.0, "tau_w [s]"), m, s.data);
    CAPTURE(rl.parts.front().residual, rl.parts.front().where);
    REQUIRE(rl.comparison.outcome == Outcome::fail);
    REQUIRE(rl.parts.front().where.find("1e-15 .. 1e-14") != std::string::npos);
    // The potential form on the FD study, the first model point's largest node difference doubled: the
    // L-infinity difference doubles, and the rate between it and its neighbour falls to 0.7.
    const Study& f = fermi_dirac();
    Study fbad = f;
    auto& psi = fbad.model[0].fields.potential_V;
    double worst = 0.0;
    std::size_t at = 0;
    for (std::size_t i = 0; i < psi.size(); ++i) {
        const double e = std::abs(psi[i] - f.limit[0].fields.potential_V[i]);
        if (e > worst) {
            worst = e;
            at = i;
        }
    }
    psi[at] += (psi[at] > f.limit[0].fields.potential_V[at] ? 1.0 : -1.0) * worst;  // the difference doubled
    REQUIRE(c::check_c10(fbad.series(doping_values, c::LimitQuantity::potential, 0, 1.0, "N"), f.device->mesh(),
                         f.data)
                .comparison.outcome == Outcome::fail);
}

TEST_CASE("C10 is NOT APPLICABLE or NOT EVALUATED without its inputs", "[verify][consistency]") {
    const Study& s = energy_balance();
    const auto& m = s.device->mesh();
    const auto q = c::LimitQuantity::terminal_current;

    auto two = s.series(tau_values, q, 0, 1.0, "tau_w");  // two values: NOT APPLICABLE
    two.x.pop_back();
    two.model.pop_back();
    const auto na = c::check_c10(two, m, s.data);
    REQUIRE(na.comparison.outcome == Outcome::not_applicable);
    REQUIRE(na.comparison.reason.find("limit parameter swept") != std::string::npos);

    auto at_floor = s.series(tau_values, q, 0, 1.0, "tau_w");  // a model point equal to the limit: no difference
    at_floor.model[2] = &s.limit[0];
    const auto ne = c::check_c10(at_floor, m, s.data);
    REQUIRE(ne.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(ne.comparison.reason.find("not above its floor") != std::string::npos);

    // A difference above one result's resolution but not above the two together: NOT EVALUATED.
    std::vector<results::BiasPoint> near(1, s.limit[0]);
    const double res = s.limit[0].terminal_current_resolution[0];
    near[0].terminal_current[0] += 1.8 * res;  // floor = res (its own) + res (the limit's)
    auto close = s.series(tau_values, q, 0, 1.0, "tau_w");
    close.model[2] = &near[0];
    const auto nc = c::check_c10(close, m, s.data);
    CAPTURE(res, nc.comparison.reason);
    REQUIRE(nc.comparison.outcome == Outcome::not_evaluated);
    REQUIRE(nc.comparison.reason.find("not above its floor") != std::string::npos);

    auto sizes = s.series(tau_values, q, 0, 1.0, "tau_w");  // model count not the values'
    sizes.model.pop_back();
    REQUIRE(c::check_c10(sizes, m, s.data).comparison.outcome == Outcome::not_evaluated);

    // Each guard names its input (the metrics library would refuse some of these too, with another reason).
    const auto because = [&](const c::LimitSeries& series, const char* why) {
        const auto r = c::check_c10(series, m, s.data);
        CAPTURE(r.comparison.reason, why);
        REQUIRE(r.comparison.outcome == Outcome::not_evaluated);
        REQUIRE(r.comparison.reason.find(why) != std::string::npos);
    };
    because(s.series({1e-14, 1e-14, 1e-16}, q, 0, 1.0, "tau_w"), "values not distinct");
    because(s.series({1e-14, -1e-15, 1e-16}, q, 0, 1.0, "tau_w"), "a value not finite and positive");
    because(s.series(tau_values, q, 0, 0.0, "tau_w"), "stated rate not finite and positive");
    because(s.series(tau_values, q, 5, 1.0, "tau_w"), "contact out of range");

    auto limits = s.series(tau_values, q, 0, 1.0, "tau_w");  // two limit results for three values
    limits.limit.push_back(&s.limit[0]);
    REQUIRE(c::check_c10(limits, m, s.data).comparison.outcome == Outcome::not_evaluated);

    auto null_model = s.series(tau_values, q, 0, 1.0, "tau_w");  // a result missing
    null_model.model[1] = nullptr;
    REQUIRE(c::check_c10(null_model, m, s.data).comparison.outcome == Outcome::not_evaluated);

    auto contact = s.series(tau_values, q, 5, 1.0, "tau_w");  // contact out of range
    REQUIRE(c::check_c10(contact, m, s.data).comparison.outcome == Outcome::not_evaluated);
}
