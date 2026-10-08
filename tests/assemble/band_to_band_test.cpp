// Band-to-band tunnelling in the drift-diffusion assembler (ARCHITECTURE.md section 11, Unit 20):
// the uniform-field gates (local Kane at A F^2 exp(-B / F) at every node, the calibrated nonlocal
// rate equal to it at each path's start, the direct-gap WKB rate equal to eq. (8), each at 0, 30
// and 45 degrees and in 2D and 3D), the finite-difference Jacobian gate for the three rates in 1D,
// 2D and 3D, pair conservation, the transverse-uniform reduction to 1D, tracing that stays in the
// semiconductor, and the refusals (local with nonlocal, no cells, a heterointerface, direct_wkb on
// silicon).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/band_to_band.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using assemble::NonlocalTunnelling;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

// GaAs with the masses of a test (m_c 0.067, m_v 0.5 m0, the textbook Gamma-valley electron and
// heavy-hole masses), for the direct-gap WKB rate only; no material set ships them.
physics::SemiconductorParameters direct_gap() {
    physics::SemiconductorParameters p = physics::gallium_arsenide_parameters;
    p.band_to_band.electron_mass = 0.067;
    p.band_to_band.hole_mass = 0.5;
    return p;
}

// A uniformly doped block with its cells, ohmic contacts on x_min and x_max; `transverse` nodes
// along y and z (the block is length_cm along x and width_cm across).
device::Device block(int D, int nodes, int transverse, double length_cm, double width_cm,
                     physics::SemiconductorParameters material = physics::silicon_parameters,
                     double donors = 1e17) {
    const auto x = uniform(0.0, length_cm, nodes);
    const auto t = uniform(0.0, width_cm, transverse);
    mesh::Mesh m = D == 1   ? *mesh::make_tensor_grid(x)
                   : D == 2 ? *mesh::make_tensor_grid(x, t)
                            : *mesh::make_tensor_grid(x, t, t);
    auto cells = D == 1   ? *mesh::TensorCells::create(m, x)
                 : D == 2 ? *mesh::TensorCells::create(m, x, t)
                          : *mesh::TensorCells::create(m, x, t, t);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"semiconductor", *physics::Semiconductor::create(material)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, donors),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}},
         .cells = std::move(cells)});
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

// A uniform field E0 [V/cm] along the unit direction u (x, y, z), n and p uniform, relative noise
// on every unknown. Electrons tunnel towards higher potential, against the field: along -u.
std::vector<double> sloped(const DriftDiffusion& dd, const device::Device& d,
                           const assemble::Scaling& s, double E0, mesh::Point u, double noise,
                           std::uint64_t seed, double n_cm3 = 1e17, double p_cm3 = 1e3) {
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

// As the Unit 19 gate: every column, central differences, the worst error over the column's
// largest entry, each column's error the smallest at steps of 1e-5, 1e-6 and 1e-7 max(|x|, 1).
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

// The scaled volume of node i.
double volume(const device::Device& d, const assemble::Scaling& s, std::size_t i) {
    return d.mesh().volumes()[i] / std::pow(s.L_D, d.mesh().dimension());
}

mesh::Point direction(double degrees) {
    const double a = degrees * std::numbers::pi / 180.0;
    return {std::cos(a), std::sin(a), 0.0};
}

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

}  // namespace

TEST_CASE("band to band assemble: refusals") {
    const device::Device d = block(1, 11, 1, 1e-5, 1e-5);
    const auto s = *assemble::make_scaling(d);
    const auto refused = [&](const device::Device& dev, assemble::PhysicsModels m) {
        const auto r = DriftDiffusion::create(dev, *assemble::make_scaling(dev), m);
        REQUIRE_FALSE(r.has_value());
        REQUIRE(r.error().code == base::ErrorCode::invalid_input);
    };
    refused(d, {.btbt_local = true, .btbt_nonlocal = NonlocalTunnelling::kane});
    refused(d, {.btbt_nonlocal = NonlocalTunnelling::direct_wkb});  // silicon
    // A GaAs region without masses is refused too; with the test's masses it is accepted.
    refused(block(1, 11, 1, 1e-5, 1e-5, physics::gallium_arsenide_parameters),
            {.btbt_nonlocal = NonlocalTunnelling::direct_wkb});
    const device::Device g = block(1, 11, 1, 1e-5, 1e-5, direct_gap());
    REQUIRE(DriftDiffusion::create(g, *assemble::make_scaling(g),
                                   {.btbt_nonlocal = NonlocalTunnelling::direct_wkb})
                .has_value());
    // No cells.
    {
        const auto x = uniform(0.0, 1e-5, 11);
        mesh::Mesh m = *mesh::make_tensor_grid(x);
        auto left = m.find_boundary("x_min")->nodes;
        auto right = m.find_boundary("x_max")->nodes;
        const device::Device bare = *device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"si", physics::silicon()}},
             .node_region = std::vector<device::RegionId>(11, 0),
             .donors = std::vector<double>(11, 1e17),
             .acceptors = std::vector<double>(11, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                          {"right", device::ContactKind::ohmic, std::move(right)}}});
        refused(bare, {.btbt_nonlocal = NonlocalTunnelling::kane});
        // The local model needs no cells.
        REQUIRE(DriftDiffusion::create(bare, *assemble::make_scaling(bare), {.btbt_local = true})
                    .has_value());
    }
    // A heterointerface (Si | Ge).
    {
        const auto x = uniform(0.0, 1e-5, 11);
        mesh::Mesh m = *mesh::make_tensor_grid(x);
        auto cells = *mesh::TensorCells::create(m, x);
        auto left = m.find_boundary("x_min")->nodes;
        auto right = m.find_boundary("x_max")->nodes;
        std::vector<device::RegionId> region(11, 0);
        for (std::size_t i = 6; i < 11; ++i) region[i] = 1;
        const device::Device hetero = *device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"si", physics::silicon()},
                         {"ge", *physics::Semiconductor::create(physics::germanium_parameters)}},
             .node_region = region,
             .donors = std::vector<double>(11, 1e17),
             .acceptors = std::vector<double>(11, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                          {"right", device::ContactKind::ohmic, std::move(right)}},
             .cells = std::move(cells)});
        refused(hetero, {.btbt_nonlocal = NonlocalTunnelling::kane});
    }
    // Cells of another mesh.
    {
        const auto x = uniform(0.0, 1e-5, 11), other = uniform(0.0, 2e-5, 11);
        mesh::Mesh m = *mesh::make_tensor_grid(x);
        const mesh::Mesh m2 = *mesh::make_tensor_grid(other);
        auto left = m.find_boundary("x_min")->nodes;
        auto right = m.find_boundary("x_max")->nodes;
        const auto r = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"si", physics::silicon()}},
             .node_region = std::vector<device::RegionId>(11, 0),
             .donors = std::vector<double>(11, 1e17),
             .acceptors = std::vector<double>(11, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                          {"right", device::ContactKind::ohmic, std::move(right)}},
             .cells = *mesh::TensorCells::create(m2, other)});
        REQUIRE_FALSE(r.has_value());
    }
    // Off: no paths, and the pattern does not change.
    const auto off = *DriftDiffusion::create(d, s);
    REQUIRE_FALSE(off.tunnelling());
    REQUIRE(off.trace_paths(std::vector<double>(off.unknowns(), 1.0)).paths.empty());
}

TEST_CASE("band to band assemble: local Kane in a uniform field at any angle") {
    // 4e6 V/cm: G = 3.5e21 1.6e13 exp(-25.75) = 3.7e23 cm^-3 s^-1.
    const double E0 = 4e6;
    const physics::KaneRate expected = physics::kane_generation(3.5e21, 1.03e8, E0);
    for (const int D : {1, 2, 3}) {
        const device::Device d = block(D, D == 3 ? 6 : 11, D == 3 ? 6 : 11, 1e-5, 1e-5);
        const auto s = *assemble::make_scaling(d);
        const auto on = *DriftDiffusion::create(d, s, {.btbt_local = true});
        const auto off = *DriftDiffusion::create(d, s);
        for (const double degrees : {0.0, 30.0, 45.0}) {
            if (D == 1 && degrees != 0.0) continue;
            // Low densities: the generation is read as a difference of residuals, whose flux
            // terms at n = 1e17 would be 1e9 times larger.
            const auto x = sloped(on, d, s, E0, direction(degrees), 0.0, 1, 1e3, 1e3);
            std::vector<double> fa(on.unknowns()), fb(on.unknowns());
            on.residual(x, fa);
            off.residual(x, fb);
            double worst = 0.0;
            for (std::size_t i = 0; i < on.node_count(); ++i) {
                const bool contact = d.mesh().points()[i][0] == 0.0 ||
                                     d.mesh().points()[i][0] == 1e-5;
                const double G = (fa[3 * i + 1] - fb[3 * i + 1]) / volume(d, s, i) * s.R0;
                const double Gp = -(fa[3 * i + 2] - fb[3 * i + 2]) / volume(d, s, i) * s.R0;
                if (contact) {
                    REQUIRE(G == 0.0);  // ohmic nodes do not generate
                    continue;
                }
                worst = std::max({worst, std::abs(G / expected.rate - 1.0),
                                  std::abs(Gp / expected.rate - 1.0)});
            }
            CAPTURE(D, degrees, worst);
            REQUIRE(worst < 1e-9);  // the field is reconstructed exactly at any angle
        }
    }
}

TEST_CASE("band to band assemble: nonlocal rates in a uniform field equal their closed forms") {
    // Silicon, kane: each reached path from a start whose straight path stays inside generates
    // A F^2 exp(-B / F) at the field itself (F = (psi_f - psi_i) / l). GaAs with the test's
    // masses, direct_wkb: eq. (11) reduces to eq. (8) where the energy window is wide.
    struct Case {
        NonlocalTunnelling kind;
        physics::SemiconductorParameters material;
        double field;
    };
    for (const Case& c : {Case{NonlocalTunnelling::kane, physics::silicon_parameters, 4e6},
                          Case{NonlocalTunnelling::direct_wkb, direct_gap(), 2e6}}) {
        for (const int D : {1, 2, 3}) {
            const double L = 6e-6;  // 60 nm
            const device::Device d =
                block(D, D == 3 ? 13 : 31, D == 1 ? 1 : (D == 3 ? 13 : 31), L, L, c.material);
            const auto s = *assemble::make_scaling(d);
            auto dd = *DriftDiffusion::create(d, s, {.btbt_nonlocal = c.kind});
            for (const double degrees : {0.0, 30.0, 45.0}) {
                if (D == 1 && degrees != 0.0) continue;
                const mesh::Point u = direction(degrees);
                const auto x = sloped(dd, d, s, c.field, u, 0.0, 1);
                dd.set_paths(dd.trace_paths(x));
                const auto states = dd.path_states(x);
                REQUIRE_FALSE(states.empty());
                const double Eg = physics::band_gap_eV(
                    *physics::Semiconductor::create(c.material), 300.0);
                const double W = Eg / c.field;  // the tunnel length in a uniform field [cm]
                double expected = 0.0;
                if (c.kind == NonlocalTunnelling::kane) {
                    expected = physics::kane_generation(3.5e21, 1.03e8, c.field).rate;
                } else {
                    const double mr = physics::tunnelling_reduced_mass(c.material.band_to_band);
                    expected = physics::wkb_uniform_field_rate(
                                   c.field * 100.0, {Eg * base::q_C, mr * base::m0_kg}) * 1e-6;
                }
                std::size_t compared = 0;
                double worst = 0.0;
                for (const auto& st : states) {
                    // Only starts whose straight path, and (WKB) a window of a gap behind and
                    // beyond it, stay inside the block.
                    const mesh::Point& p = d.mesh().points()[st.start];
                    const double reach = c.kind == NonlocalTunnelling::kane ? 1.6 * W : 2.2 * W;
                    const double back = c.kind == NonlocalTunnelling::kane ? 0.0 : 1.2 * W;
                    bool inside = true;
                    for (int a = 0; a < D; ++a) {
                        const double end = p[a] - reach * u[a], behind = p[a] + back * u[a];
                        if (std::min(end, behind) < 0.0 || std::max(end, behind) > L) {
                            inside = false;
                        }
                    }
                    if (!inside) continue;
                    ++compared;
                    REQUIRE(st.reached);
                    REQUIRE(close(st.field_V_per_cm, c.field, 1e-10));
                    worst = std::max(worst, std::abs(st.rate_cm3_s / expected - 1.0));
                }
                CAPTURE(static_cast<int>(c.kind), D, degrees, compared, worst);
                REQUIRE(compared > 0);
                REQUIRE(worst < (c.kind == NonlocalTunnelling::kane ? 1e-9 : 1e-6));
            }
        }
    }
}

TEST_CASE("band to band assemble: the finite-difference Jacobian gate") {
    // Uniform fields with noise on every unknown (psi by 0.02 V_T), frozen paths traced at the
    // state: the local rate, and the nonlocal rates with their crossing, mean field, WKB
    // integrals, band extrema and deposit weights all moving with psi.
    struct Case {
        assemble::PhysicsModels models;
        physics::SemiconductorParameters material;
        double field;
    };
    for (const Case& c :
         {Case{{.btbt_local = true}, physics::silicon_parameters, 4e6},
          Case{{.btbt_nonlocal = NonlocalTunnelling::kane}, physics::silicon_parameters, 4e6},
          Case{{.btbt_nonlocal = NonlocalTunnelling::direct_wkb}, direct_gap(), 2e6}}) {
        for (const int D : {1, 2, 3}) {
            const double L = 3e-6;
            const device::Device d = block(D, D == 3 ? 6 : 13, D == 3 ? 4 : 5, L, L, c.material);
            const auto s = *assemble::make_scaling(d);
            auto dd = *DriftDiffusion::create(d, s, c.models);
            const auto x = sloped(dd, d, s, c.field, direction(D == 1 ? 0.0 : 20.0), 0.02, 3);
            if (dd.tunnelling()) {
                dd.set_paths(dd.trace_paths(x));
                REQUIRE_FALSE(dd.paths().paths.empty());
            }
            const double error = fd_jacobian_error(dd, x);
            CAPTURE(c.models.btbt_local, static_cast<int>(c.models.btbt_nonlocal), D, error);
            REQUIRE(error <= 1e-6);
        }
    }
}

TEST_CASE("band to band assemble: nonlocal pairs balance exactly") {
    for (const NonlocalTunnelling kind :
         {NonlocalTunnelling::kane, NonlocalTunnelling::direct_wkb}) {
        const auto material =
            kind == NonlocalTunnelling::kane ? physics::silicon_parameters : direct_gap();
        const device::Device d = block(2, 21, 9, 4e-6, 2e-6, material);
        const auto s = *assemble::make_scaling(d);
        auto on = *DriftDiffusion::create(d, s, {.btbt_nonlocal = kind});
        const auto off = *DriftDiffusion::create(d, s);
        const auto x = sloped(on, d, s, kind == NonlocalTunnelling::kane ? 4e6 : 2e6,
                              direction(25.0), 0.02, 5, 1e3, 1e3);  // small fluxes (above)
        on.set_paths(on.trace_paths(x));
        std::vector<double> fa(on.unknowns()), fb(on.unknowns());
        on.residual(x, fa);
        off.residual(x, fb);
        double electrons = 0.0, holes = 0.0, scale = 0.0;
        for (std::size_t i = 0; i < on.node_count(); ++i) {
            electrons += fa[3 * i + 1] - fb[3 * i + 1];
            holes -= fa[3 * i + 2] - fb[3 * i + 2];
            scale += std::abs(fa[3 * i + 2] - fb[3 * i + 2]);
        }
        CAPTURE(static_cast<int>(kind), electrons, holes);
        REQUIRE(holes > 0.0);
        REQUIRE(std::abs(electrons - holes) <= 1e-13 * scale);
    }
}

TEST_CASE("band to band assemble: a transverse-uniform device reduces to 1D") {
    // The same x axis in 1D, 2D and 3D, a field along x: every node's generation (local), and
    // each path's rate and every node's electron and hole terms (nonlocal), per volume, agree.
    const double L = 6e-6;
    const auto run = [&](int D, const assemble::PhysicsModels& models) {
        const device::Device d = block(D, 31, D == 1 ? 1 : 4, L, 2e-6);
        const auto s = *assemble::make_scaling(d);
        auto on = *DriftDiffusion::create(d, s, models);
        const auto off = *DriftDiffusion::create(d, s);
        const auto x = sloped(on, d, s, 4e6, {1.0, 0.0, 0.0}, 0.0, 1);
        if (on.tunnelling()) on.set_paths(on.trace_paths(x));
        std::vector<double> fa(on.unknowns()), fb(on.unknowns());
        on.residual(x, fa);
        off.residual(x, fb);
        // Per x position (node index along x): electron and hole generation per volume.
        std::vector<double> g(2 * 31, 0.0);
        for (std::size_t i = 0; i < on.node_count(); ++i) {
            if (i % 31 != i && D > 1 && (i / 31) % 4 != 1) continue;  // one interior row
            const std::size_t k = i % 31;
            g[2 * k] = (fa[3 * i + 1] - fb[3 * i + 1]) / volume(d, s, i) * s.R0;
            g[2 * k + 1] = (fa[3 * i + 2] - fb[3 * i + 2]) / volume(d, s, i) * s.R0;
        }
        return g;
    };
    for (const assemble::PhysicsModels& m :
         {assemble::PhysicsModels{.btbt_local = true},
          assemble::PhysicsModels{.btbt_nonlocal = NonlocalTunnelling::kane}}) {
        const auto one = run(1, m), two = run(2, m), three = run(3, m);
        double largest = 0.0, worst = 0.0;
        for (std::size_t k = 0; k < one.size(); ++k) {
            largest = std::max(largest, std::abs(one[k]));
            worst = std::max({worst, std::abs(two[k] - one[k]), std::abs(three[k] - one[k])});
        }
        CAPTURE(m.btbt_local, worst, largest);
        REQUIRE(largest > 0.0);
        REQUIRE(worst <= 1e-12 * largest);
    }
}

TEST_CASE("band to band assemble: paths stay in the semiconductor") {
    // A silicon block under an oxide (y above 2e-6 cm), the field pointing into the oxide at 40
    // degrees: no sample has weight on an insulator node, and paths still start beside it.
    const auto x = uniform(0.0, 6e-6, 25), y = uniform(0.0, 3e-6, 13);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    auto cells = *mesh::TensorCells::create(m, x, y);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> donors(n, 1e17);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][1] > 2e-6 + 1e-12) {
            region[i] = 1;
            donors[i] = 0.0;
        }
    }
    std::vector<mesh::NodeId> left, right;
    for (const mesh::NodeId i : m.find_boundary("x_min")->nodes) {
        if (region[static_cast<std::size_t>(i)] == 0) left.push_back(i);
    }
    for (const mesh::NodeId i : m.find_boundary("x_max")->nodes) {
        if (region[static_cast<std::size_t>(i)] == 0) right.push_back(i);
    }
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"si", physics::silicon()},
                     {"oxide", physics::silicon_dioxide()}},
         .node_region = region,
         .donors = donors,
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}},
         .cells = std::move(cells)});
    const auto s = *assemble::make_scaling(d);
    auto dd = *DriftDiffusion::create(d, s, {.btbt_nonlocal = NonlocalTunnelling::kane});
    // Electrons go against the field: a field at 220 degrees sends them up and to the right.
    const auto state = sloped(dd, d, s, 4e6, direction(220.0), 0.0, 1);
    const assemble::TunnelPaths paths = dd.trace_paths(state);
    REQUIRE_FALSE(paths.paths.empty());
    bool beside = false;
    for (const assemble::TunnelPath& p : paths.paths) {
        if (d.mesh().points()[p.start][1] == 2e-6) beside = true;
        for (const auto& sample : p.forward) {
            for (int k = 0; k < sample.count; ++k) {
                const auto q = static_cast<std::size_t>(k);
                if (sample.weights[q] == 0.0) continue;
                REQUIRE_FALSE(d.is_insulator(static_cast<mesh::NodeId>(sample.nodes[q])));
            }
        }
    }
    REQUIRE(beside);
}
