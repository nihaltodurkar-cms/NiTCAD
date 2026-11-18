// The extensibility test (Unit 57 GG8, GQ5 a, review R5): a test-only equation adds one node field
// chi to the drift-diffusion composition through SystemData::create_extended (composition.hpp),
// without a change to System or to any existing term. chi is an unknown at every node (DofMap's
// layout, FD5), held at 0 on insulator and contact nodes; its row is the box-method form of
//     d2 chi / du2 = psi - psi_ref - S(u),   u = x / L
// (L the device's length along x, psi scaled by V_T, S a manufactured source, 0 by default). A
// second test-only term feeds chi back, two-way: k V chi into the electron row and -k V chi into
// the hole row of every semiconductor node (V its scaled box volume), the sign of a generation.
// The checks: (i) the DofMap and the pattern grow by one field; (ii) both terms' partials, and
// every existing term's at the wider block, match finite differences, the cross blocks (n, p rows
// in chi; chi row in psi) included; (iii) Newton converges quadratically in its last iterations;
// (iv) at k = 0, psi, n and p equal the system's without chi to 1e-12; (v) at k != 0 on a 1D bar
// with a manufactured solution, the error falls at second order under two refinements.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numbers>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/equilibrium.hpp"
#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/newton.hpp"
#include "composition.hpp"
#include "scaled_device.hpp"
#include "terms.hpp"

using namespace NiTCAD;
using assemble::FieldId;
using assemble::System;
using assemble::TermDomain;
using Data = assemble::SystemData;
using State = assemble::SystemState;
using Terms = assemble::SystemData::Terms;

namespace {

// No FieldId names chi: the test's own id, outside the enumerators (fields.hpp is not edited).
constexpr auto chi_id = static_cast<FieldId>(1000);
constexpr assemble::Field chi_field{chi_id, "chi (test)", "1", assemble::FieldLocation::node,
                                    assemble::FieldKind::unknown, assemble::FieldDomain::all};
constexpr FieldId chi_rows[] = {chi_id};
constexpr FieldId chi_columns[] = {FieldId::potential, chi_id};
constexpr FieldId carrier_rows[] = {FieldId::electrons, FieldId::holes};

bool held(const Data& s, std::size_t i) { return s.contact_[i] >= 0 || s.insulator_[i] != 0; }

// The chi rows (assigned): on a held node chi; elsewhere the edges' sum of (chi_b - chi_a) / du
// less the box length times (psi - psi_ref - S). 1D meshes (the box length is half of each
// incident edge's du).
class ChiRows final : public Data::Term {
public:
    ChiRows(std::vector<double> u, double psi_ref, std::vector<double> source)
        : u_(std::move(u)), psi_ref_(psi_ref), source_(std::move(source)) {}
    assemble::TermInfo info() const noexcept override {
        return {"chi rows (test)", TermDomain::global, chi_rows, chi_columns, true};
    }
    void pattern(const Data& s, const State&, std::vector<linalg::Triplet>& t) const override {
        const std::size_t m = s.m_, c = s.built_in_;
        for (const Data::EdgeTerm& e : s.edges_) {
            t.push_back({static_cast<linalg::Index>(m * e.a + c),
                         static_cast<linalg::Index>(m * e.b + c), 0.0});
            t.push_back({static_cast<linalg::Index>(m * e.b + c),
                         static_cast<linalg::Index>(m * e.a + c), 0.0});
        }
    }
    void add(const Data& s, const State& st, const Terms::Evaluation& ev) const override {
        const std::size_t m = s.m_, c = s.built_in_, n = s.node_count();
        std::vector<double> box(n, 0.0);
        for (const Data::EdgeTerm& e : s.edges_) {
            const double du = std::abs(u_[e.b] - u_[e.a]);
            box[e.a] += 0.5 * du;
            box[e.b] += 0.5 * du;
        }
        for (std::size_t i = 0; i < n; ++i) {
            const double chi = ev.x[m * i + c];
            if (held(s, i)) {
                ev.f[m * i + c] = chi;
                if (ev.jacobian()) Terms::at(s, ev, i, c, c) = 1.0;
                continue;
            }
            const double S = source_.empty() ? 0.0 : source_[i];
            ev.f[m * i + c] = -box[i] * (ev.x[m * i] - psi_ref_ - S);
            if (ev.jacobian()) Terms::at(s, ev, i, c, 0) = -box[i];
        }
        for (const Data::EdgeTerm& e : s.edges_) {
            const double w = 1.0 / std::abs(u_[e.b] - u_[e.a]);
            const double flux = w * (ev.x[m * e.b + c] - ev.x[m * e.a + c]);
            if (!held(s, e.a)) {
                ev.f[m * e.a + c] += flux;
                if (ev.jacobian()) {
                    Terms::at(s, ev, e.a, c, c) -= w;
                    ev.values[assemble::detail::position(st.pattern_, m * e.a + c, m * e.b + c)] +=
                        w;
                }
            }
            if (!held(s, e.b)) {
                ev.f[m * e.b + c] -= flux;
                if (ev.jacobian()) {
                    Terms::at(s, ev, e.b, c, c) -= w;
                    ev.values[assemble::detail::position(st.pattern_, m * e.b + c, m * e.a + c)] +=
                        w;
                }
            }
        }
    }

private:
    std::vector<double> u_;
    double psi_ref_;
    std::vector<double> source_;
};

// chi fed back into the continuity rows: +k V chi (electrons), -k V chi (holes).
class ChiSource final : public Data::Term {
public:
    explicit ChiSource(double k) : k_(k) {}
    assemble::TermInfo info() const noexcept override {
        return {"chi source (test)", TermDomain::nodes, carrier_rows, chi_rows, false};
    }
    void add(const Data& s, const State&, const Terms::Evaluation& ev) const override {
        const std::size_t m = s.m_, c = s.built_in_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!Terms::semiconductor(s, i)) continue;
            const double kV = k_ * s.volume_[i];
            ev.f[m * i + 1] += kV * ev.x[m * i + c];
            ev.f[m * i + 2] -= kV * ev.x[m * i + c];
            if (ev.jacobian()) {
                Terms::at(s, ev, i, 1, c) += kV;
                Terms::at(s, ev, i, 2, c) -= kV;
            }
        }
    }

private:
    double k_;
};

// u = (x - x_min) / L at every node of a 1D device.
std::vector<double> unit_coordinate(const device::Device& d) {
    const auto& points = d.mesh().points();
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (const auto& p : points) {
        lo = std::min(lo, p[0]);
        hi = std::max(hi, p[0]);
    }
    std::vector<double> u;
    for (const auto& p : points) u.push_back((p[0] - lo) / (hi - lo));
    return u;
}

Data::CompositionExtension chi_extension(const device::Device& d, double k, double psi_ref,
                                         std::vector<double> source = {}) {
    return {{chi_field},
            {std::make_shared<ChiRows>(unit_coordinate(d), psi_ref, std::move(source)),
             std::make_shared<ChiSource>(k)}};
}

// 1D: oxide [-5e-7, 0] with a metal electrode at x_min, p-silicon [0, 2e-5] (1e17) with an ohmic
// contact at x_max, the oxide-silicon interface declared.
device::Device mos() {
    std::vector<double> x;
    for (int i = 0; i <= 40; ++i) x.push_back(-5e-7 + (2e-5 + 5e-7) * i / 40.0);
    for (const double y : {-4e-7, -3e-7, -2e-7, -1e-7, 0.0}) x.push_back(y);
    std::ranges::sort(x);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
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

// 1D silicon bar [0, L] (donors uniform, 1e16 by default), elements equal, ohmic contacts at both ends.
constexpr double bar_length_cm = 1e-4;
device::Device bar(int elements, double donors = 1e16) {
    std::vector<double> x;
    for (int i = 0; i <= elements; ++i) x.push_back(bar_length_cm * i / elements);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, donors),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}}});
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

// Whether (row, col) is in a matrix's pattern.
bool has(const linalg::SparseMatrix& a, std::size_t row, std::size_t col) {
    for (auto q = a.row_offsets()[row]; q < a.row_offsets()[row + 1]; ++q) {
        if (static_cast<std::size_t>(a.col_indices()[q]) == col) return true;
    }
    return false;
}

// One term's partials against central differences of its own residual (terms_test.cpp's gate):
// the worst column error over the column's largest entry in the whole Jacobian, the smallest at
// steps of 1e-5, 1e-6 and 1e-7 max(|x|, 1).
double fd_error(const System& sys, const Data::Term& t, std::vector<double> x) {
    const std::size_t n = sys.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    auto jacobian = sys.make_jacobian();
    Terms::evaluate(sys, t, x, nullptr, f, jacobian);
    const auto j = dense(jacobian);
    auto scratch = sys.make_jacobian();
    sys.evaluate(x, fp, scratch);
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
            Terms::evaluate(sys, t, x, nullptr, fp, scratch);
            x[c] = base - h;
            const double down = x[c];
            Terms::evaluate(sys, t, x, nullptr, fm, scratch);
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

// The charge-neutral potential's state (the solve layer's usual start), chi 0.
std::vector<double> neutral_start(const device::Device& d, const assemble::Scaling& s,
                                  const System& sys) {
    const System eq = *assemble::create_equilibrium(d, s);
    return sys.state_from_potential(assemble::charge_neutral_potential(eq));
}

results::ConvergenceRecord newton(const System& sys, std::vector<double>& x) {
    auto solver = *linalg::LinearSolver::create({});
    results::ConvergenceRecord record;
    solve::NewtonOptions options;
    options.tol_update = 1e-11;
    REQUIRE(solve::newton_solve(sys, x, options, solver, record).has_value());
    return record;
}

}  // namespace

TEST_CASE("GG8: the DofMap and the pattern grow by one field", "[extensibility]") {
    const device::Device d = mos();
    const auto s = *assemble::make_scaling(d);
    const System plain = *System::create(d, s);
    const System ext = *Data::create_extended(d, s, {}, chi_extension(d, 0.5, 0.0));
    const std::size_t nodes = d.mesh().node_count();
    REQUIRE(plain.dofs().block() == 3);
    REQUIRE(ext.dofs().block() == 4);
    CHECK(ext.unknowns() == 4 * nodes);
    CHECK(ext.dofs().slot(chi_id) == 3);
    for (const FieldId f : {FieldId::potential, FieldId::electrons, FieldId::holes}) {
        CHECK(ext.dofs().slot(f) == plain.dofs().slot(f));
    }
    CHECK(ext.dofs().field_of(4 * 7 + 3).name == "chi (test)");
    // The pattern: every entry of the plain system at its place in the wider block, then the
    // node blocks' new row and column (7 entries per node) and chi's two per edge.
    const linalg::SparseMatrix pp = plain.make_jacobian(), ep = ext.make_jacobian();
    const std::size_t edges = Data::of(ext).edges_.size();
    CHECK(ep.nonzeros() == pp.nonzeros() + 7 * nodes + 2 * edges);
    // Every entry of the plain pattern, at its place in the wider block.
    std::size_t kept = 0;
    for (std::size_t row = 0; row < static_cast<std::size_t>(pp.rows()); ++row) {
        const std::size_t node = row / 3, r = row % 3;
        for (auto q = pp.row_offsets()[row]; q < pp.row_offsets()[row + 1]; ++q) {
            const auto col = static_cast<std::size_t>(pp.col_indices()[q]);
            kept += has(ep, 4 * node + r, 4 * (col / 3) + col % 3) ? 1 : 0;
        }
    }
    CHECK(kept == pp.nonzeros());
    // The table's terms, then the extension's, in order.
    const auto plain_terms = Terms::list(plain);
    const auto ext_terms = Terms::list(ext);
    REQUIRE(ext_terms.size() == plain_terms.size() + 2);
    for (std::size_t k = 0; k < plain_terms.size(); ++k) {
        CHECK(ext_terms[k]->info().name == plain_terms[k]->info().name);
    }
    CHECK(ext_terms[ext_terms.size() - 2]->info().name == "chi rows (test)");
    CHECK(ext_terms.back()->info().name == "chi source (test)");
}

TEST_CASE("GG8: every term's partials match finite differences, the cross blocks included",
          "[extensibility]") {
    const device::Device d = mos();
    const auto s = *assemble::make_scaling(d);
    const System ext = *Data::create_extended(d, s, {}, chi_extension(d, 0.5, -2.0));
    const std::size_t nodes = d.mesh().node_count();
    // Off equilibrium: psi a ramp, the densities perturbed, chi nonzero.
    std::vector<double> psi(nodes);
    for (std::size_t i = 0; i < nodes; ++i) psi[i] = -3.0 + 0.11 * static_cast<double>(i % 41);
    std::vector<double> x = ext.state_from_potential(psi);
    for (std::size_t i = 0; i < nodes; ++i) {
        x[4 * i + 1] *= 1.0 + 0.01 * static_cast<double>(i % 7);
        x[4 * i + 2] *= 1.0 - 0.005 * static_cast<double>(i % 5);
        x[4 * i + 3] = 0.3 * std::sin(0.4 * static_cast<double>(i));
    }
    ext.stamp_contacts(x);
    for (const auto& t : Terms::list(ext)) {
        const std::string name(t->info().name);
        const double error = fd_error(ext, *t, x);
        CAPTURE(name, error);
        CHECK(error <= 1e-6);
    }
    // The cross blocks are there and nonzero: the n and p rows in chi, the chi row in psi.
    std::vector<double> f(ext.unknowns());
    auto j = ext.make_jacobian();
    ext.evaluate(x, f, j);
    const auto J = dense(j);
    const auto& data = Data::of(ext);
    std::size_t semiconductor = 0, free_chi = 0;
    for (std::size_t i = 0; i < nodes; ++i) {
        if (Terms::semiconductor(data, i)) {
            ++semiconductor;
            CHECK(J[4 * i + 1][4 * i + 3] == 0.5 * data.volume_[i]);
            CHECK(J[4 * i + 2][4 * i + 3] == -0.5 * data.volume_[i]);
        }
        if (!held(data, i)) {
            ++free_chi;
            CHECK(J[4 * i + 3][4 * i] < 0.0);
        } else {  // held at 0: the row is chi alone
            CHECK(f[4 * i + 3] == x[4 * i + 3]);
            CHECK(J[4 * i + 3][4 * i] == 0.0);
        }
    }
    CHECK(semiconductor > 30);
    CHECK(free_chi == semiconductor);
}

TEST_CASE("GG8: Newton converges quadratically; at k = 0 psi, n, p are the plain system's",
          "[extensibility]") {
    const device::Device d = mos();
    const auto s = *assemble::make_scaling(d);
    const System plain = *System::create(d, s);
    const std::size_t nodes = d.mesh().node_count();

    std::vector<double> x0 = neutral_start(d, s, plain);
    const results::ConvergenceRecord r0 = newton(plain, x0);

    for (const double k : {0.0, 1e-11}) {
        CAPTURE(k);
        const System ext = *Data::create_extended(d, s, {}, chi_extension(d, k, -2.0));
        std::vector<double> x = neutral_start(d, s, ext);
        const results::ConvergenceRecord r = newton(ext, x);
        REQUIRE(r.converged);
        // (iii) quadratic: the last pair of corrections (d_j, d_j+1) with d_j < 0.1 and d_j+1
        // above rounding (1e-13), its observed order log(d_j+1) / log(d_j) at least 1.8, and
        // the corrections falling from d_j on. Measured (2026-11-18): k = 0, 1.7e-5 -> 7.3e-11
        // (order 2.12); k = 1e-11, 1.2e-6 -> 5.5e-13 (2.07).
        std::vector<double> d_k;
        for (const auto& rec : r.iterations) d_k.push_back(rec.update);
        CAPTURE(d_k);
        std::size_t j = d_k.size();
        for (std::size_t q = 0; q + 1 < d_k.size(); ++q) {
            if (d_k[q] < 0.1 && d_k[q + 1] >= 1e-13) j = q;
        }
        REQUIRE(j < d_k.size());
        const double order = std::log(d_k[j + 1]) / std::log(d_k[j]);
        CAPTURE(j, order);
        CHECK(order >= 1.8);
        for (std::size_t q = j; q + 1 < d_k.size(); ++q) CHECK(d_k[q + 1] < d_k[q]);
        if (k != 0.0) {
            // Two-way: chi moved the densities (measured: the minority electrons by a factor of
            // about 5.5e3).
            double moved = 0.0;
            for (std::size_t i = 0; i < nodes; ++i) {
                for (std::size_t q = 1; q < 3; ++q) {
                    const double a = x0[3 * i + q], b = x[4 * i + q];
                    if (a > 0.0) moved = std::max(moved, std::abs(b - a) / a);
                }
            }
            CAPTURE(moved);
            CHECK(moved > 1.0);
            continue;
        }
        // (iv) at k = 0: psi, n, p equal the plain system's to 1e-12 (relative; psi against
        // max(|psi|, 1)).
        double worst = 0.0;
        for (std::size_t i = 0; i < nodes; ++i) {
            worst = std::max(worst, std::abs(x[4 * i] - x0[3 * i]) / std::max(std::abs(x0[3 * i]), 1.0));
            for (std::size_t q = 1; q < 3; ++q) {
                const double a = x0[3 * i + q], b = x[4 * i + q];
                worst = std::max(worst, a == 0.0 ? std::abs(b) : std::abs(b - a) / a);
            }
        }
        CAPTURE(worst, r0.iterations.size(), r.iterations.size());
        CHECK(worst <= 1e-12);
        // chi itself is nonzero off the held nodes (its row reads psi).
        double chi = 0.0;
        for (std::size_t i = 0; i < nodes; ++i) chi = std::max(chi, std::abs(x[4 * i + 3]));
        CHECK(chi > 1e-3);
    }
}

TEST_CASE("GG8: the Newton hooks measure and apply chi's correction", "[extensibility]") {
    // On an intrinsic bar at equilibrium (psi, n, p exact from the start, k = 0) with
    // psi_ref = psi_c - 200, chi = -100 u (1 - u) (|chi| up to 25, well beyond the clip of
    // max_update = 5), so only chi has to move, over several clipped steps. update_size must
    // measure its correction, or the first iteration stops with chi clipped short; apply_update
    // must add it. At the end every chi row is at its rounding. Intrinsic, because on a doped
    // device the minority densities' corrections carry the rounding of a large chi correction
    // (about eps |dchi| against a density near 1e-12 of the scale): measured on the 1D MOS, they
    // keep update_size near 1e-4 while chi moves, so Newton would not stop early even without
    // chi in the measure.
    const device::Device d = bar(32, 0.0);
    const auto s = *assemble::make_scaling(d);
    const std::size_t nodes = d.mesh().node_count();
    const System plain = *System::create(d, s);
    const double psi_c = neutral_start(d, s, plain)[0];
    const System far = *Data::create_extended(d, s, {}, chi_extension(d, 0.0, psi_c - 200.0));
    std::vector<double> x = far.state_from_potential(std::vector<double>(nodes, psi_c));
    const results::ConvergenceRecord r = newton(far, x);
    REQUIRE(r.converged);
    CHECK(r.iterations.size() >= 3);
    std::vector<double> f(far.unknowns());
    far.residual(x, f);
    double chi = 0.0, row = 0.0;
    for (std::size_t i = 0; i < nodes; ++i) {
        chi = std::max(chi, std::abs(x[4 * i + 3]));
        row = std::max(row, std::abs(f[4 * i + 3]));
    }
    CAPTURE(chi, row, r.iterations.size());
    CHECK(chi > 10.0);
    CHECK(row <= 1e-10 * chi);
}

TEST_CASE("GG8: at k != 0 the manufactured solution converges at second order",
          "[extensibility]") {
    // chi* = A sin(pi u), psi* = psi_c, n* = n_c, p* = p_c (the bar's uniform equilibrium). chi's
    // source S = (psi* - psi_ref) - chi*'' = A pi^2 sin(pi u) with psi_ref = psi_c, and the
    // manufactured generation G = -k chi* R0 (the TVQ3 hook: + V G / R0 into the electron rows,
    // - V G / R0 into the hole rows) cancels the k V chi* of the source term, so the exact
    // fields solve the continuous equations and every error is the discretization's.
    // k is in scaled rate units (R0 = D0 Ns / L_D^2, about 1e33 cm^-3 s^-1 here): at k = 1e-11
    // the discrete chi's O(h^2) error moves the minority holes (p_c about 2e4 cm^-3) by about
    // 10% at 32 elements, still in the linear range; at 1e-3 by a factor of 1e7. The gates are
    // chi's and the holes' L-infinity orders, in [1.8, 2.3] (p_th = 2). psi's and n's errors are
    // then at rounding (measured 6e-14, 9e-14 at 32 elements) and are not gated; at k = 1e-3
    // they too fell at order 2.0. Measured (2026-11-18): chi 1.6e-3, 4.0e-4, 1.0e-4; p relative
    // 9.5e-2, 2.4e-2, 5.9e-3.
    constexpr double A = 2.0, k = 1e-11;
    std::vector<double> chi_error, p_error;
    for (const int elements : {32, 64, 128}) {
        const device::Device d = bar(elements);
        const auto s = *assemble::make_scaling(d);
        const System plain = *System::create(d, s);
        const std::size_t nodes = d.mesh().node_count();
        std::vector<double> exact = neutral_start(d, s, plain);  // uniform: the exact state
        const double psi_c = exact[0], p_c = exact[2];
        const std::vector<double> u = unit_coordinate(d);
        std::vector<double> S(nodes);
        auto ms = std::make_shared<assemble::ManufacturedSources>();
        for (std::size_t i = 0; i < nodes; ++i) {
            const double chi = A * std::sin(std::numbers::pi * u[i]);
            S[i] = std::numbers::pi * std::numbers::pi * chi;
            ms->electron_generation.push_back(-k * chi * s.R0);
            ms->hole_generation.push_back(-k * chi * s.R0);
        }
        assemble::PhysicsModels models;
        models.manufactured = ms.get();
        const System ext = *Data::create_extended(d, s, models, chi_extension(d, k, psi_c, S));
        std::vector<double> x = ext.state_from_potential(std::vector<double>(nodes, psi_c));
        const results::ConvergenceRecord r = newton(ext, x);
        REQUIRE(r.converged);
        double ec = 0.0, ep = 0.0;
        for (std::size_t i = 0; i < nodes; ++i) {
            ec = std::max(ec, std::abs(x[4 * i + 3] - A * std::sin(std::numbers::pi * u[i])));
            ep = std::max(ep, std::abs(x[4 * i + 2] - p_c) / p_c);
        }
        chi_error.push_back(ec);
        p_error.push_back(ep);
    }
    CAPTURE(chi_error, p_error);
    for (std::size_t l = 1; l < chi_error.size(); ++l) {
        const double order = std::log2(chi_error[l - 1] / chi_error[l]);
        CAPTURE(l, order);
        CHECK(order >= 1.8);
        CHECK(order <= 2.3);
        const double hole_order = std::log2(p_error[l - 1] / p_error[l]);
        CAPTURE(hole_order);
        CHECK(hole_order >= 1.8);
        CHECK(hole_order <= 2.3);
    }
    CHECK(chi_error.back() <= 1e-3 * A);
    CHECK(p_error.back() <= 1e-2);
}
