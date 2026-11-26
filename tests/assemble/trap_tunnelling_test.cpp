// Hurkx trap-assisted tunnelling in drift-diffusion (Unit 67, DECISIONS.md TAC6, TAC7, TAC9 and
// review focus 1, 2): the node rate is SRH with each lifetime over 1 + Gamma(|E|), ohmic nodes
// and zero-field nodes keep plain SRH, the Jacobian against central differences (psi of the
// neighbours included, with local band-to-band tunnelling and impact ionization on too), and the
// refused combinations.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/trap_tunnelling.hpp"

using namespace NiTCAD;
using assemble::System;
using assemble::TrapTunnelling;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

// A uniformly doped silicon block, ohmic contacts on x_min and x_max.
device::Device block(int D, int nodes, int transverse, double length_cm, double width_cm,
                     double donors = 1e17) {
    const auto x = uniform(0.0, length_cm, nodes);
    const auto t = uniform(0.0, width_cm, transverse);
    mesh::Mesh m = D == 1   ? *mesh::make_tensor_grid(x)
                   : D == 2 ? *mesh::make_tensor_grid(x, t)
                            : *mesh::make_tensor_grid(x, t, t);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"semiconductor", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, donors),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}}});
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

// A uniform field E0 [V/cm] along the unit direction u, n and p uniform, relative noise on every
// unknown.
std::vector<double> sloped(const System& dd, const device::Device& d, const assemble::Scaling& s,
                           double E0, mesh::Point u, double noise, std::uint64_t seed,
                           double n_cm3, double p_cm3) {
    std::vector<double> x(dd.unknowns());
    Noise z{seed};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        const auto& p = d.mesh().points()[i];
        const double along = p[0] * u[0] + p[1] * u[1] + p[2] * u[2];
        x[3 * i] = -E0 * along / s.V_T + noise * z.next();
        x[3 * i + 1] = n_cm3 / s.Ns * (1.0 + noise * z.next());
        x[3 * i + 2] = p_cm3 / s.Ns * (1.0 + noise * z.next());
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

// As the Unit 19/20 gate: every column, central differences, the worst error over the column's
// largest entry, each column's error the smallest at steps of 1e-5, 1e-6 and 1e-7 max(|x|, 1).
double fd_jacobian_error(const System& dd, std::vector<double> x) {
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

mesh::Point direction(double degrees) {
    const double a = degrees * std::numbers::pi / 180.0;
    return {std::cos(a), std::sin(a), 0.0};
}

bool on_contact(const device::Device& d, std::size_t i, double length) {
    const double x = d.mesh().points()[i][0];
    return x == 0.0 || x == length;
}

}  // namespace

TEST_CASE("tat assemble: the rate is SRH with each lifetime over 1 + Gamma") {
    // BGN off: n_ie is silicon's intrinsic density, so a_n = ln(Nc / n_i), a_p = ln(Nv / n_i).
    // Unequal masses and densities, so that electrons and holes cannot be swapped unseen. The
    // state carries no current (n ~ e^psi, p ~ e^-psi on a 4 nm block, quasi-Fermi levels flat):
    // the residual difference is the recombination alone, not a small difference of fluxes.
    const double E0 = 5e5, L = 4e-7, N = 1e17, T = 300.0;
    const physics::Semiconductor si = physics::silicon();
    const double ni = physics::intrinsic_density(si, T);
    const double an = std::log(physics::conduction_band_dos(si, T) / ni);
    const double ap = std::log(physics::valence_band_dos(si, T) / ni);
    const double gn = physics::hurkx_enhancement(an, 0.2, T, E0).gamma;
    const double gp = physics::hurkx_enhancement(ap, 0.6, T, E0).gamma;
    const double tn = physics::scharfetter_lifetime(si, physics::Carrier::electron, N);
    const double tp = physics::scharfetter_lifetime(si, physics::Carrier::hole, N);
    const double n0 = 1e12, p0 = 1e10;  // np > n_i^2: net recombination
    const auto ratio = [&](double n, double p) {
        return (tp * (n + ni) + tn * (p + ni)) /
               (tp / (1.0 + gp) * (n + ni) + tn / (1.0 + gn) * (p + ni));
    };
    CAPTURE(gn, gp);
    REQUIRE(ratio(n0, p0) > 1.5);  // the field matters at 5e5 V/cm
    for (const int D : {1, 2, 3}) {
        const device::Device d = block(D, D == 3 ? 6 : 11, D == 3 ? 6 : 11, L, L, N);
        const auto s = *assemble::make_scaling(d);
        const auto none = *System::create(d, s, {.srh = false, .bgn = false});
        const auto srh = *System::create(d, s, {.bgn = false});
        const auto tat = *System::create(
            d, s,
            {.bgn = false, .trap_assisted_tunnelling = TrapTunnelling::hurkx, .tat_mass_n = 0.2,
             .tat_mass_p = 0.6});
        for (const double degrees : {0.0, 30.0}) {
            if (D == 1 && degrees != 0.0) continue;
            auto x = sloped(tat, d, s, E0, direction(degrees), 0.0, 1, n0, p0);
            for (std::size_t i = 0; i < tat.node_count(); ++i) {
                const double shift = x[3 * i] + 0.5 * E0 * L / s.V_T;  // psi from mid-block
                x[3 * i + 1] *= std::exp(shift);
                x[3 * i + 2] *= std::exp(-shift);
            }
            std::vector<double> f0(tat.unknowns()), fs(tat.unknowns()), ft(tat.unknowns());
            none.residual(x, f0);
            srh.residual(x, fs);
            tat.residual(x, ft);
            double worst = 0.0;
            for (std::size_t i = 0; i < tat.node_count(); ++i) {
                const double Rs = f0[3 * i + 1] - fs[3 * i + 1];
                const double Rt = f0[3 * i + 1] - ft[3 * i + 1];
                const double Rt_p = ft[3 * i + 2] - f0[3 * i + 2];
                if (on_contact(d, i, L)) continue;  // their rows are the contact conditions
                REQUIRE(Rs > 0.0);
                const double q = ratio(x[3 * i + 1] * s.Ns, x[3 * i + 2] * s.Ns);
                worst = std::max({worst, std::abs(Rt / Rs / q - 1.0), std::abs(Rt_p / Rs / q - 1.0)});
            }
            CAPTURE(D, degrees, worst);
            REQUIRE(worst < 1e-9);
        }
    }
}

TEST_CASE("tat contact nodes plain srh") {
    // Contact rows are assigned by the contact term, so compare the term's own write: with the
    // contact rows of both systems identical, the whole residual at the contact nodes agrees.
    const double L = 1e-5;
    const device::Device d = block(2, 11, 5, L, L);
    const auto s = *assemble::make_scaling(d);
    const auto srh = *System::create(d, s);
    const auto tat = *System::create(d, s, {.trap_assisted_tunnelling = TrapTunnelling::hurkx});
    const auto x = sloped(tat, d, s, 5e5, direction(0.0), 0.0, 1, 1e12, 3e11);
    std::vector<double> fs(tat.unknowns()), ft(tat.unknowns());
    srh.residual(x, fs);
    tat.residual(x, ft);
    for (std::size_t i = 0; i < tat.node_count(); ++i) {
        if (!on_contact(d, i, L)) continue;
        for (std::size_t r = 0; r < 3; ++r) REQUIRE(ft[3 * i + r] == fs[3 * i + r]);
    }
}

TEST_CASE("tat zero field keeps SRH") {
    // psi uniform: F = 0 at every node, Gamma = 0, the residual and Jacobian values are SRH's
    // bit for bit (plain SRH still added at a zero-field node).
    const double L = 1e-5;
    for (const int D : {1, 2}) {
        const device::Device d = block(D, 11, 5, L, L);
        const auto s = *assemble::make_scaling(d);
        const auto srh = *System::create(d, s);
        const auto tat = *System::create(d, s, {.trap_assisted_tunnelling = TrapTunnelling::hurkx});
        const auto x = sloped(tat, d, s, 0.0, direction(0.0), 0.0, 1, 1e12, 3e11);
        std::vector<double> fs(tat.unknowns()), ft(tat.unknowns());
        auto js = srh.make_jacobian();
        auto jt = tat.make_jacobian();
        srh.evaluate(x, fs, js);
        tat.evaluate(x, ft, jt);
        REQUIRE(ft == fs);
        const auto a = dense(js), b = dense(jt);
        REQUIRE(a == b);
    }
}

TEST_CASE("tat assemble: the finite-difference Jacobian gate") {
    // A noisy uniform field (psi noise 0.02 V_T) at 6e5 V/cm: Gamma and its psi partials through
    // the reconstructed field of every node, also with local band-to-band tunnelling and impact
    // ionization in the same system.
    const auto hurkx = TrapTunnelling::hurkx;
    for (const assemble::PhysicsModels& m :
         {assemble::PhysicsModels{.trap_assisted_tunnelling = hurkx},
          assemble::PhysicsModels{.btbt_local = true, .trap_assisted_tunnelling = hurkx},
          assemble::PhysicsModels{.impact_ionization = true, .trap_assisted_tunnelling = hurkx}}) {
        for (const int D : {1, 2, 3}) {
            const double L = 3e-6;
            const device::Device d = block(D, D == 3 ? 6 : 13, D == 3 ? 4 : 5, L, L);
            const auto s = *assemble::make_scaling(d);
            const auto dd = *System::create(d, s, m);
            const auto x = sloped(dd, d, s, 6e5, direction(D == 1 ? 0.0 : 20.0), 0.02, 3, 1e12,
                                  3e11);
            const double error = fd_jacobian_error(dd, x);
            CAPTURE(m.btbt_local, m.impact_ionization, D, error);
            REQUIRE(error <= 1e-6);
        }
    }
}

TEST_CASE("tat assemble: refusals") {
    const device::Device d = block(1, 11, 1, 1e-5, 1e-5);
    const auto s = *assemble::make_scaling(d);
    const auto refused = [&](assemble::PhysicsModels m, const std::string& names) {
        const auto r = System::create(d, s, m);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == base::ErrorCode::invalid_input);
        CAPTURE(r.error().message);
        CHECK(r.error().message.find(names) != std::string::npos);
    };
    const auto hurkx = TrapTunnelling::hurkx;
    refused({.srh = false, .trap_assisted_tunnelling = hurkx}, "srh");
    refused({.trap_assisted_tunnelling = hurkx, .electrothermal = true}, "electrothermal");
    refused({.trap_assisted_tunnelling = hurkx, .energy_balance = true}, "energy balance");
    refused({.quantum_confinement = assemble::QuantumConfinement::density_gradient,
             .trap_assisted_tunnelling = hurkx},
            "density gradient");
    for (const double bad : {0.0, -0.25, std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity()}) {
        refused({.trap_assisted_tunnelling = hurkx, .tat_mass_n = bad}, "tat_mass_n");
        refused({.trap_assisted_tunnelling = hurkx, .tat_mass_p = bad}, "tat_mass_p");
    }
    // Off, the masses are not read.
    REQUIRE(System::create(d, s, {.tat_mass_n = 0.0}).has_value());
}
