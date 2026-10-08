// Band-to-band tunnelling in the solves (ARCHITECTURE.md section 11, Unit 20): a silicon Zener
// diode (p+ and n+ 5e19 cm^-3) with the local and the calibrated nonlocal Kane rate: the reverse
// current follows the Kane form and does not plateau, the terminal current is the generated pairs'
// (pure generation), the equilibrium current (no D factor) measured, transverse-uniform 2D and 3D
// reduce to 1D, relocation settles and does not depend on the starting paths (a cycle and a set
// that never settles are reported), transient runs hold the steady current and re-trace under a
// ramp, the small-signal conductance is the DC derivative, a 2D corner junction tunnels before a
// planar one, the direct-gap WKB rate on GaAs with a test's masses, the cost, and the run record.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <expected>
#include <string>
#include <span>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/band_to_band.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "NiTCAD/solve/waveform.hpp"
#include "../../src/solve/relocation.hpp"

using namespace NiTCAD;
using assemble::NonlocalTunnelling;

namespace {

constexpr double length_cm = 1e-5;  // 100 nm
constexpr int nodes = 201;          // 0.5 nm spacing

std::vector<double> uniform(double a, double b, int count) {
    std::vector<double> v(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (count - 1);
    return v;
}

// An abrupt p+ n+ junction at the middle of x, with its cells; for D > 1 extruded over `width`
// in two nodes along y (and z). Anode (p side) on x_min, cathode on x_max.
device::Device zener(int D = 1, double N = 5e19,
                     physics::SemiconductorParameters material = physics::silicon_parameters,
                     double width = 2e-6) {
    const auto x = uniform(0.0, length_cm, nodes);
    const std::vector<double> t{0.0, width};
    mesh::Mesh m = D == 1   ? *mesh::make_tensor_grid(x)
                   : D == 2 ? *mesh::make_tensor_grid(x, t)
                            : *mesh::make_tensor_grid(x, t, t);
    auto cells = D == 1   ? *mesh::TensorCells::create(m, x)
                 : D == 2 ? *mesh::TensorCells::create(m, x, t)
                          : *mesh::TensorCells::create(m, x, t, t);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        (m.points()[i][0] < 0.5 * length_cm ? acceptors : donors)[i] = N;
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"semiconductor", *physics::Semiconductor::create(material)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .cells = std::move(cells)});
}

solve::BiasOptions with(bool local, NonlocalTunnelling nonlocal = NonlocalTunnelling::off) {
    solve::BiasOptions o;
    o.models.btbt_local = local;
    o.models.btbt_nonlocal = nonlocal;
    return o;
}

std::vector<std::vector<double>> reverse(double to, double step) {
    std::vector<std::vector<double>> points;
    for (double v = 0.0; v >= to - 1e-12; v -= step) points.push_back({v, 0.0});
    return points;
}

// The largest edge field of a 1D state [V/cm].
double peak_field(const results::BiasPoint& p) {
    const auto& psi = p.fields.potential_V;
    const double h = length_cm / (nodes - 1);
    double F = 0.0;
    for (std::size_t i = 0; i + 1 < psi.size(); ++i) {
        F = std::max(F, std::abs(psi[i + 1] - psi[i]) / h);
    }
    return F;
}

// q times the pairs the local Kane rate generates per second in a 1D state [A/cm^2]: each
// interior node at the mean of its two edges' fields (assemble's reconstruction on this grid).
double local_pairs(const results::BiasPoint& p, const device::Device& d) {
    const auto& psi = p.fields.potential_V;
    const double h = length_cm / (nodes - 1);
    double sum = 0.0;
    for (std::size_t i = 1; i + 1 < psi.size(); ++i) {
        const double F = std::abs(psi[i + 1] - psi[i - 1]) / (2.0 * h);
        sum += physics::kane_generation(3.5e21, 1.03e8, F).rate * d.mesh().volumes()[i];
    }
    return base::q_C * sum;
}

// q times the pairs of the nonlocal paths [A / cm^(3-D)].
double path_pairs(const results::BiasPoint& p, const device::Device& d) {
    double sum = 0.0;
    for (const results::TunnelPath& t : p.tunnel_paths) {
        sum += t.generation_cm3_s * d.mesh().volumes()[t.start_node];
    }
    return base::q_C * sum;
}

}  // namespace

TEST_CASE("btbt: the Zener diode with the local Kane rate") {
    // Reverse to -3 V in 0.25 V steps. The current is the generated pairs' within the current's
    // resolution (5e-4 A/cm^2 at this doping; SRH adds 4e-12); it grows by 2.1 times over the
    // last 0.25 V to -3 V (no plateau). For a linear field profile I ~ F^4 exp(-B / F) (the
    // generating width grows as F, and the integral of exp(-B / F(x)) over it as F^2 / B), so
    // ln(I / F^4) against 1 / F has slope -B: measured over -1.5 to -3 V.
    const device::Device d = zener();
    const auto s = solve::sweep_bias(d, reverse(-3.0, 0.25), with(true));
    REQUIRE(s.has_value());
    REQUIRE(!s->stopped);
    double previous = 0.0;
    for (const results::BiasPoint& p : s->points) {
        const double I = -p.terminal_current[0];
        CAPTURE(p.bias_V[0], I, local_pairs(p, d), p.terminal_current_resolution[0]);
        REQUIRE(std::abs(I - local_pairs(p, d)) <= p.terminal_current_resolution[0]);
        if (p.bias_V[0] <= -0.5) REQUIRE(I > previous);
        previous = I;
    }
    const auto at = [&](double V) {
        for (const auto& p : s->points) {
            if (std::abs(p.bias_V[0] - V) < 1e-9) return &p;
        }
        return &s->points.front();
    };
    const double grow = at(-3.0)->terminal_current[0] / at(-2.75)->terminal_current[0];
    // Least-squares slope of ln(I / F^4) against 1 / F.
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int m = 0;
    for (const auto& p : s->points) {
        if (p.bias_V[0] > -1.5 + 1e-9) continue;
        const double F = peak_field(p);
        const double X = 1.0 / F, Y = std::log(-p.terminal_current[0] / std::pow(F, 4));
        sx += X;
        sy += Y;
        sxx += X * X;
        sxy += X * Y;
        ++m;
    }
    const double slope = (m * sxy - sx * sy) / (m * sxx - sx * sx);
    CAPTURE(grow, slope, slope / -1.03e8);
    REQUIRE(grow > 2.0);
    REQUIRE(std::abs(slope / -1.03e8 - 1.0) < 0.05);
}

TEST_CASE("btbt: the Zener diode with the calibrated nonlocal rate") {
    // The same sweep with the nonlocal Kane rate: every point settles (at most two re-tracings),
    // its paths reach their crossings, and the current is the paths' pairs within the resolution.
    // Against the local rate (measured, ARCHITECTURE.md 6.2): 33 times smaller at -1 V, 1.9 at -3
    // V (the mean field over a path is below the peak field the local rate takes).
    const device::Device d = zener();
    const auto s = solve::sweep_bias(d, reverse(-3.0, 0.25), with(false, NonlocalTunnelling::kane));
    REQUIRE(s.has_value());
    REQUIRE(!s->stopped);
    for (const results::BiasPoint& p : s->points) {
        const double I = -p.terminal_current[0];
        CAPTURE(p.bias_V[0], I, path_pairs(p, d), p.path_relocations);
        REQUIRE(p.path_relocations <= 2);
        REQUIRE(std::abs(I - path_pairs(p, d)) <= p.terminal_current_resolution[0]);
        double largest = 0.0;
        for (const auto& t : p.tunnel_paths) largest = std::max(largest, t.generation_cm3_s);
        for (const auto& t : p.tunnel_paths) {
            if (t.generation_cm3_s >= 1e-10 * largest && largest > 0.0) REQUIRE(t.reached);
        }
    }
    const auto local = solve::sweep_bias(d, reverse(-3.0, 0.25), with(true));
    REQUIRE(local.has_value());
    const double r1 = local->points[4].terminal_current[0] / s->points[4].terminal_current[0];
    const double r3 = local->points[12].terminal_current[0] / s->points[12].terminal_current[0];
    CAPTURE(r1, r3);
    REQUIRE(s->points[12].terminal_current[0] / s->points[11].terminal_current[0] > 2.0);
    REQUIRE(r1 > r3);  // the two approach each other at high bias
}

TEST_CASE("btbt: the equilibrium current of pure generation") {
    // No Hurkx D factor (owner's decision): at 0 V the junction field generates pairs, and the
    // terminal current is theirs (within the resolution). The generation itself is exact, the
    // current resolved only to 5e-4 (5e19) and 1e-3 A/cm^2 (1e20). At 5e19 (peak field 2.8e6 V/cm
    // at 0 V) the pairs are 2.9e-9 A/cm^2, 2.9e-6 of the -1 V current; at 1e20 3.4e-4 A/cm^2,
    // 7.7e-5 of it (measured; the known limit of ARCHITECTURE.md 6.2).
    for (const double N : {5e19, 1e20}) {
        const device::Device d = zener(1, N);
        const std::vector<std::vector<double>> points{{0.0, 0.0}, {-1.0, 0.0}};
        const auto s = solve::sweep_bias(d, points, with(true));
        REQUIRE(s.has_value());
        REQUIRE(!s->stopped);
        const results::BiasPoint& zero = s->points[0];
        const double I0 = -zero.terminal_current[0], I1 = -s->points[1].terminal_current[0];
        CAPTURE(N, I0, local_pairs(zero, d), I1, I0 / I1, zero.terminal_current_resolution[0]);
        REQUIRE(std::abs(I0 - local_pairs(zero, d)) <= zero.terminal_current_resolution[0]);
        const double share = local_pairs(zero, d) / I1;
        REQUIRE(share > 0.0);  // pure generation: never 0 in a junction field
        if (N == 5e19) REQUIRE(share < 1e-5);
        if (N == 1e20) REQUIRE(share > 3e-5);  // grows steeply with the doping
    }
    // Off, the 0 V current is within its resolution of 0.
    const device::Device d = zener(1, 1e20);
    const auto off = solve::solve_bias(d, std::vector<double>{0.0, 0.0}, with(false));
    REQUIRE(off.has_value());
    REQUIRE(std::abs(off->terminal_current[0]) <= off->terminal_current_resolution[0]);
}

TEST_CASE("btbt: with generation alone the current is the generated pairs'") {
    // SRH, Auger and radiative recombination off: every carrier crossing a contact was generated
    // by tunnelling, at -2 V with both rates (currents of 0.3 and 0.08 A/cm^2).
    const device::Device d = zener();
    for (const bool local : {true, false}) {
        solve::BiasOptions o =
            with(local, local ? NonlocalTunnelling::off : NonlocalTunnelling::kane);
        o.models.srh = false;
        o.models.auger = false;
        o.models.radiative = false;
        const auto s = solve::sweep_bias(d, reverse(-2.0, 0.5), o);
        REQUIRE(s.has_value());
        REQUIRE(!s->stopped);
        const results::BiasPoint& p = s->points.back();
        const double I = -p.terminal_current[0];
        const double pairs = local ? local_pairs(p, d) : path_pairs(p, d);
        CAPTURE(local, I, pairs, I / pairs - 1.0, p.terminal_current_resolution[0]);
        REQUIRE(std::abs(I - pairs) <= p.terminal_current_resolution[0]);
        REQUIRE(std::abs(I / pairs - 1.0) < 1e-2);
    }
}

TEST_CASE("btbt: transverse-uniform 2D and 3D reproduce 1D") {
    // At -2 V, the current per unit width (2e-6 cm) and area reproduces the 1D current density.
    for (const bool local : {true, false}) {
        const auto o = with(local, local ? NonlocalTunnelling::off : NonlocalTunnelling::kane);
        const auto points = reverse(-2.0, 0.5);
        const auto one = solve::sweep_bias(zener(1), points, o);
        const auto two = solve::sweep_bias(zener(2), points, o);
        const auto three = solve::sweep_bias(zener(3), points, o);
        REQUIRE(!one->stopped);
        REQUIRE(!two->stopped);
        REQUIRE(!three->stopped);
        const double J1 = one->points.back().terminal_current[0];
        const double J2 = two->points.back().terminal_current[0] / 2e-6;
        const double J3 = three->points.back().terminal_current[0] / 4e-12;
        CAPTURE(local, J1, J2 / J1 - 1.0, J3 / J1 - 1.0);
        REQUIRE(std::abs(J2 / J1 - 1.0) < 1e-9);
        REQUIRE(std::abs(J3 / J1 - 1.0) < 1e-9);
    }
}

TEST_CASE("btbt: relocation settles and does not depend on the starting paths") {
    // -2 V reached in one step from equilibrium (paths traced at the stamped state) and by a ramp
    // (paths carried from point to point): the same current.
    const device::Device d = zener();
    const auto o = with(false, NonlocalTunnelling::kane);
    const auto direct = solve::solve_bias(d, std::vector<double>{-2.0, 0.0}, o);
    const auto ramp = solve::sweep_bias(d, reverse(-2.0, 0.25), o);
    REQUIRE(direct.has_value());
    REQUIRE(!ramp->stopped);
    const double a = direct->terminal_current[0], b = ramp->points.back().terminal_current[0];
    CAPTURE(a, b, a / b - 1.0, direct->path_relocations);
    REQUIRE(std::abs(a - b) <= 1e-6 * std::abs(b));
}

namespace {

// A stand-in for the relocation loop: the traced set is chosen by the test.
struct Paths {
    std::vector<int> paths;
};
struct Relocating {
    std::vector<Paths> traces;  // returned in turn, the last one repeated
    std::size_t next = 0;
    Paths frozen;
    [[nodiscard]] bool tunnelling() const { return true; }
    [[nodiscard]] const Paths& paths() const { return frozen; }
    void set_paths(Paths p) { frozen = std::move(p); }
    Paths trace_paths(std::span<const double>) {
        return traces[std::min(next++, traces.size() - 1)];
    }
    [[nodiscard]] bool paths_agree(const Paths& a, const Paths& b, std::span<const double>) const {
        return a.paths == b.paths;
    }
};

}  // namespace

TEST_CASE("btbt: relocation reports a cycle and a set that never settles") {
    std::vector<double> x(3, 0.0);
    int newton_runs = 0;
    const auto newton = [&]() -> std::expected<void, base::Error> {
        ++newton_runs;
        return {};
    };
    int relocations = 0;
    {  // settles: A, B, B
        Relocating r{.traces = {Paths{{1}}, Paths{{2}}, Paths{{2}}}};
        REQUIRE(solve::detail::solve_relocating(r, x, newton, relocations).has_value());
        REQUIRE(relocations == 1);
    }
    {  // cycles: A, B, A
        Relocating r{.traces = {Paths{{1}}, Paths{{2}}, Paths{{1}}, Paths{{2}}}};
        const auto e = solve::detail::solve_relocating(r, x, newton, relocations);
        REQUIRE_FALSE(e.has_value());
        REQUIRE(e.error().code == base::ErrorCode::non_convergence);
        REQUIRE(e.error().message.find("cycle") != std::string::npos);
    }
    {  // never settles: a new set every time
        std::vector<Paths> t;
        for (int k = 0; k < 20; ++k) t.push_back({{k}});
        Relocating r{.traces = t};
        const auto e = solve::detail::solve_relocating(r, x, newton, relocations);
        REQUIRE_FALSE(e.has_value());
        REQUIRE(e.error().message.find("did not settle") != std::string::npos);
        REQUIRE(relocations == solve::detail::max_relocations);
    }
}

TEST_CASE("btbt: transient runs hold the steady current and re-trace under a ramp",
          "[.btbt]") {
    const device::Device d = zener();
    const auto o = with(false, NonlocalTunnelling::kane);
    const auto sweep = solve::sweep_bias(d, reverse(-2.5, 0.25), o);
    REQUIRE(!sweep->stopped);
    const results::BiasPoint& at2 = sweep->points[8];   // -2 V
    const results::BiasPoint& at25 = sweep->points[10];  // -2.5 V
    // Constant -2 V: the steady current throughout, no re-tracing.
    {
        const std::vector<solve::Waveform> w{solve::Waveform::constant(-2.0),
                                             solve::Waveform::constant(0.0)};
        const auto t = solve::solve_transient(d, w, {.steady = o, .t_end_s = 1e-9}, &at2.fields);
        REQUIRE(t.has_value());
        REQUIRE(!t->stopped);
        double drift = 0.0;
        for (const auto& p : t->points) {
            const double I = p.terminal_current[0];
            drift = std::max(drift, std::abs(I / at2.terminal_current[0] - 1.0));
        }
        CAPTURE(drift, t->retraced_steps);
        REQUIRE(drift < 1e-6);
        REQUIRE(t->retraced_steps == 0);
    }
    // A ramp to -2.5 V over 1 ns, held to 3 ns: the paths are re-traced on the way and the end
    // is the steady -2.5 V state; fixed steps of 20 and 10 ps agree.
    const auto ramp = [&](double dt) {
        const std::vector<solve::Waveform> w{*solve::Waveform::ramp(-2.0, -2.5, 0.0, 1e-9),
                                             solve::Waveform::constant(0.0)};
        return solve::solve_transient(d, w,
                                      {.steady = o, .t_end_s = 3e-9, .dt_initial_s = dt,
                                       .dt_min_s = dt, .dt_max_s = dt, .adaptive = false},
                                      &at2.fields);
    };
    const auto a = ramp(2e-11), b = ramp(1e-11);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    REQUIRE(!a->stopped);
    REQUIRE(!b->stopped);
    const double Ia = a->points.back().terminal_current[0];
    const double Ib = b->points.back().terminal_current[0];
    CAPTURE(a->retraced_steps, a->lagged_steps, Ia, Ib, at25.terminal_current[0]);
    REQUIRE(a->retraced_steps > 0);
    REQUIRE(std::abs(Ia / at25.terminal_current[0] - 1.0) < 1e-4);
    REQUIRE(std::abs(Ia / Ib - 1.0) < 1e-4);
}

TEST_CASE("btbt: the small-signal conductance is the DC derivative") {
    // At -2 V, f = 0, with both rates: Y against fourth-order differences of the DC current
    // (h = 0.05 V). Measured 3.7e-5 (local) and 9.5e-5 (nonlocal); the nonlocal differences
    // themselves scatter about Y by -3e-5 to 1.3e-4 between h = 0.1 and 0.0125 V (each DC current
    // carries the resolution of this doping), so the bound is 3e-4.
    const device::Device d = zener();
    for (const bool local : {true, false}) {
        const auto o = with(local, local ? NonlocalTunnelling::off : NonlocalTunnelling::kane);
        const auto sweep = solve::sweep_bias(d, reverse(-2.0, 0.25), o);
        REQUIRE(!sweep->stopped);
        const results::BiasPoint& dc = sweep->points.back();
        const double h = 0.05;
        const auto I = [&](double V) {
            return solve::solve_bias(d, std::vector<double>{V, 0.0}, o, &dc.fields)
                ->terminal_current[0];
        };
        const double G =
            (8.0 * (I(-2.0 + h) - I(-2.0 - h)) - (I(-2.0 + 2 * h) - I(-2.0 - 2 * h))) / (12.0 * h);
        const auto ac = solve::solve_small_signal(d, std::vector<std::vector<double>>{{-2.0, 0.0}},
                                                  {.steady = o, .frequencies_Hz = {0.0}},
                                                  &dc.fields);
        REQUIRE(ac.has_value());
        REQUIRE(!ac->stopped);
        const double Y = ac->conductance(0, 0, 0, 0);
        CAPTURE(local, G, Y, Y / G - 1.0);
        REQUIRE(std::abs(Y / G - 1.0) < 3e-4);
    }
}

namespace {

// 2D: p+ 5e19 with an n+ 5e19 region on x >= 50 nm, over the whole height (planar) or only for
// y <= 25 nm (curved: its corner at (50, 25) nm); 100 nm by 50 nm, 1 nm spacing. Anode on x_min,
// cathode on the n+ part of x_max.
device::Device junction_2d(bool curved) {
    const auto x = uniform(0.0, 1e-5, 101), y = uniform(0.0, 5e-6, 51);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    auto cells = *mesh::TensorCells::create(m, x, y);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        const bool n_side = p[0] >= 5e-6 - 1e-12 && (!curved || p[1] <= 2.5e-6 + 1e-12);
        (n_side ? donors : acceptors)[i] = 5e19;
    }
    std::vector<mesh::NodeId> cathode;
    for (const mesh::NodeId v : m.find_boundary("x_max")->nodes) {
        if (donors[static_cast<std::size_t>(v)] > 0.0) cathode.push_back(v);
    }
    auto anode = m.find_boundary("x_min")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .cells = std::move(cells)});
}

}  // namespace

TEST_CASE("btbt: a curved junction tunnels before a planar one", "[.btbt]") {
    // Nonlocal Kane. The field crowds at the n+ region's corner: the curved junction's tunnelling
    // current per unit of its junction length exceeds the planar one's at every bias, and it
    // reaches 1e-3 A/cm per cm of junction at a smaller reverse bias.
    const auto o = with(false, NonlocalTunnelling::kane);
    const auto points = reverse(-2.0, 0.25);
    const auto planar = solve::sweep_bias(junction_2d(false), points, o);
    const auto curved = solve::sweep_bias(junction_2d(true), points, o);
    REQUIRE(!planar->stopped);
    REQUIRE(!curved->stopped);
    // Junction lengths: planar 50 nm across; curved 25 nm across plus 50 nm along y = 25 nm.
    const double Lp = 5e-6, Lc = 7.5e-6;
    double onset_p = 0.0, onset_c = 0.0;
    for (std::size_t k = 0; k < points.size(); ++k) {
        const double jp = -planar->points[k].terminal_current[0] / Lp;
        const double jc = -curved->points[k].terminal_current[0] / Lc;
        if (onset_p == 0.0 && jp >= 1e-3) onset_p = points[k][0];
        if (onset_c == 0.0 && jc >= 1e-3) onset_c = points[k][0];
        if (points[k][0] <= -1.0) {
            CAPTURE(points[k][0], jp, jc);
            REQUIRE(jc > jp);
        }
    }
    CAPTURE(onset_p, onset_c);
    REQUIRE(onset_c != 0.0);
    REQUIRE(onset_c > onset_p - 1e-9);  // a smaller reverse bias (or the same step)
}

TEST_CASE("btbt: the direct-gap WKB rate on GaAs with a test's masses") {
    // A GaAs p+ n+ junction (5e19) with m_c 0.067 and m_v 0.5 m0 (the masses of a test; no
    // material set ships them): the sweep settles and the current is the paths' pairs, within the
    // resolution and the Newton tolerance (the current reaches 5e5 A/cm^2 at -2 V: a direct gap
    // tunnels far more readily than silicon's fit).
    physics::SemiconductorParameters gaas = physics::gallium_arsenide_parameters;
    gaas.band_to_band.electron_mass = 0.067;
    gaas.band_to_band.hole_mass = 0.5;
    const device::Device d = zener(1, 5e19, gaas);
    const auto s = solve::sweep_bias(d, reverse(-2.0, 0.25),
                                     with(false, NonlocalTunnelling::direct_wkb));
    REQUIRE(s.has_value());
    REQUIRE(!s->stopped);
    for (const results::BiasPoint& p : s->points) {
        const double I = -p.terminal_current[0];
        CAPTURE(p.bias_V[0], I, path_pairs(p, d), p.path_relocations);
        REQUIRE(std::abs(I - path_pairs(p, d)) <= p.terminal_current_resolution[0] + 1e-7 * I);
    }
    REQUIRE(-s->points.back().terminal_current[0] > -s->points[4].terminal_current[0]);
}

TEST_CASE("btbt: the cost of the nonlocal paths", "[.btbt]") {
    // The 2D curved junction at -2 V: the Jacobian's entries and the time of a Newton evaluation
    // and factorization with the paths, against the model off (measured, ARCHITECTURE.md 6.2).
    const device::Device d = junction_2d(true);
    const auto s = *assemble::make_scaling(d);
    const auto dc = solve::solve_bias(d, std::vector<double>{-2.0, 0.0},
                                      with(false, NonlocalTunnelling::kane));
    REQUIRE(dc.has_value());
    std::vector<double> x(3 * d.mesh().node_count());
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        x[3 * i] = dc->fields.potential_V[i] / s.V_T;
        x[3 * i + 1] = dc->fields.n_cm3[i] / s.Ns;
        x[3 * i + 2] = dc->fields.p_cm3[i] / s.Ns;
    }
    const auto time = [&](const assemble::DriftDiffusion& dd) {
        auto j = dd.make_jacobian();
        std::vector<double> f(dd.unknowns());
        auto solver = *linalg::LinearSolver::create({});
        const auto t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < 5; ++k) {
            dd.evaluate(x, f, j);
            REQUIRE(solver.factorize(j).has_value());
        }
        return std::pair{j.nonzeros(),
                         std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                                 .count() / 5.0};
    };
    auto on = *assemble::DriftDiffusion::create(d, s, {.btbt_nonlocal = NonlocalTunnelling::kane});
    (void)on.set_bias(std::vector<double>{-2.0, 0.0});
    on.set_paths(on.trace_paths(x));
    auto off = *assemble::DriftDiffusion::create(d, s);
    (void)off.set_bias(std::vector<double>{-2.0, 0.0});
    const auto [nz_on, t_on] = time(on);
    const auto [nz_off, t_off] = time(off);
    CAPTURE(on.paths().paths.size(), nz_on, nz_off, t_on, t_off);
    REQUIRE(nz_on < 1.2 * nz_off);   // measured 444 paths, 1.033 times
    REQUIRE(t_on < 2.0 * t_off + 0.05);  // measured 1.14 times
}

TEST_CASE("btbt: a trace carries the paths") {
    // trace_bias with the nonlocal rate to -2 V: the paths are re-traced after each accepted
    // point (lagging one step), and the end point is the swept solve's within the resolution and
    // the corrector's tolerance.
    const device::Device d = zener();
    solve::TraceOptions o;
    o.steady = with(false, NonlocalTunnelling::kane);
    o.start_V = {0.0, 0.0};
    o.end_V = -2.0;
    o.step_V = 0.1;
    o.max_step_V = 0.5;
    const auto t = solve::trace_bias(d, o);
    REQUIRE(t.has_value());
    REQUIRE(!t->stopped);
    const auto s = solve::solve_bias(d, std::vector<double>{-2.0, 0.0}, o.steady);
    REQUIRE(s.has_value());
    const double a = t->points.back().terminal_current[0], b = s->terminal_current[0];
    int retraced = 0;
    for (const auto& p : t->points) retraced += p.path_relocations;
    CAPTURE(a, b, a / b - 1.0, t->points.size(), retraced);
    REQUIRE(t->points.back().bias_V[0] == -2.0);
    REQUIRE(retraced > 0);
    // From -1 V, where the paths carry current from the start: every point is solved with them.
    solve::TraceOptions from = o;
    from.start_V = {-1.0, 0.0};
    const auto u = solve::trace_bias(d, from);
    REQUIRE(u.has_value());
    REQUIRE(!u->stopped);
    for (const auto& p : u->points) REQUIRE_FALSE(p.tunnel_paths.empty());
    REQUIRE(std::abs(a - b) <= s->terminal_current_resolution[0] + 1e-6 * std::abs(b));
}

TEST_CASE("btbt: the run record") {
    const device::Device d = zener();
    const std::vector<std::vector<double>> p{{-1.0, 0.0}};
    const auto names = [](const results::RunRecord& r) {
        std::vector<std::string> v;
        for (const auto& [name, value] : r.settings) v.push_back(name);
        return v;
    };
    const auto off = solve::make_run_record(d, with(false), p);
    const auto local = solve::make_run_record(d, with(true), p);
    const auto nonlocal = solve::make_run_record(d, with(false, NonlocalTunnelling::kane), p);
    const auto has = [&](const results::RunRecord& r, const char* name) {
        const std::vector<std::string> v = names(r);
        return std::ranges::find(v, name) != v.end();
    };
    REQUIRE_FALSE(has(off, "models.btbt_local"));
    REQUIRE_FALSE(has(off, "models.btbt_nonlocal"));
    REQUIRE(has(local, "models.btbt_local"));
    REQUIRE(has(nonlocal, "models.btbt_nonlocal"));
    REQUIRE(off.input_identity != local.input_identity);
    REQUIRE(local.input_identity != nonlocal.input_identity);
    // The coefficients are in the identity when on.
    physics::SemiconductorParameters other = physics::silicon_parameters;
    other.band_to_band.B_V_per_cm = 1.1e8;
    const device::Device d2 = zener(1, 5e19, other);
    REQUIRE(solve::make_run_record(d2, with(true), p).input_identity != local.input_identity);
    REQUIRE(solve::make_run_record(d2, with(false), p).input_identity == off.input_identity);
}
