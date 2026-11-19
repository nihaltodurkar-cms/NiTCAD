// L0 tests of the 1D Schrodinger-Poisson reference (verify/reference/schrodinger_poisson.hpp; Unit 62 slice 0,
// QP14). Every part against a second method:
//   - the tridiagonal eigensolver against the closed form of tridiag(-1, 2, -1) and against jacobi_eigen (n = 12);
//   - the Schrodinger discretisation: the infinite well against the discrete closed form (to rounding) and the
//     continuum one (observed order 2, uniform and smoothly graded grids); a finite well with a mass step against its
//     transcendental equations solved by Brent; the triangular well against the Airy zeros of NIST DLMF Table 9.9.1
//     (read 2026-11-19, provenance PA: a_1 ... a_5 = -2.33810 74105, -4.08794 94441, -5.52055 98281,
//     -6.78670 80901, -7.94413 35871);
//   - the Poisson part (classical electrons) against mos_capacitor.hpp (Q_s and the band bending);
//   - the self-consistent quantum solution: converged residuals, Gauss's law, orthonormal states, and refinement of
//     the grid, the subband count and the Schrodinger depth.
// Bounds were fixed before the first run (Unit 62 slice 0).
#include "NiTCAD/verify/reference/schrodinger_poisson.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <limits>
#include <numbers>
#include <vector>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/jacobi_eigen.hpp"
#include "NiTCAD/verify/reference/mos_capacitor.hpp"

namespace ref = NiTCAD::verify::reference;

namespace {

constexpr double pi = std::numbers::pi;
// CODATA 2018 (exact SI values; hbar = h / 2 pi).
constexpr double q_e = 1.602176634e-19, h_planck = 6.62607015e-34, m0 = 9.1093837015e-31, k_eV = 8.617333262e-5,
                 eps0 = 8.8541878128e-14;
const double hbar = h_planck / (2.0 * pi);
const double a_eVcm2 = hbar * hbar / (2.0 * m0) / q_e * 1e4;  // hbar^2 / (2 m0) [eV cm^2]
const double kT = k_eV * 300.0;

double observed_order(double e1, double e2) { return std::log2(e1 / e2); }

std::vector<double> uniform(double z0, double z1, std::size_t cells) {
    std::vector<double> z(cells + 1);
    for (std::size_t i = 0; i <= cells; ++i) z[i] = z0 + (z1 - z0) * static_cast<double>(i) / static_cast<double>(cells);
    z.back() = z1;
    return z;
}

}  // namespace

TEST_CASE("tridiagonal_eigen: tridiag(-1, 2, -1) against its closed form", "[verify][schrodinger_poisson]") {
    const std::size_t n = 400;
    const std::vector<double> d(n, 2.0), e(n - 1, -1.0);
    const auto r = ref::tridiagonal_eigen(d, e, 10);
    REQUIRE(r.has_value());
    for (std::size_t k = 0; k < 10; ++k) {
        const double exact = 2.0 - 2.0 * std::cos(static_cast<double>(k + 1) * pi / static_cast<double>(n + 1));
        CAPTURE(k, r->values[k], exact);
        CHECK(std::abs(r->values[k] - exact) <= 1e-14);
        CHECK(r->residuals[k] <= 1e-13);
        for (std::size_t j = 0; j < k; ++j) {
            double ip = 0.0;
            for (std::size_t i = 0; i < n; ++i) ip += r->vectors[k][i] * r->vectors[j][i];
            CHECK(std::abs(ip) <= 1e-12);
        }
    }
    CHECK(ref::sturm_count(d, e, 0.0) == 0);
    CHECK(ref::sturm_count(d, e, 4.0) == n);
    CHECK(ref::sturm_count(d, e, 2.0) == n / 2);
}

TEST_CASE("tridiagonal_eigen: against jacobi_eigen on n = 12", "[verify][schrodinger_poisson]") {
    const std::size_t n = 12;
    std::vector<double> d(n), e(n - 1);
    for (std::size_t i = 0; i < n; ++i) d[i] = 3.0 * std::sin(static_cast<double>(i + 1));
    for (std::size_t i = 0; i + 1 < n; ++i) e[i] = std::cos(static_cast<double>(2 * i + 1));
    std::vector<std::complex<double>> a(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) a[i * n + i] = d[i];
    for (std::size_t i = 0; i + 1 < n; ++i) a[i * n + i + 1] = a[(i + 1) * n + i] = e[i];
    const auto j = ref::jacobi_eigen(a, n);
    REQUIRE(j.status == ref::JacobiStatus::ok);
    const auto t = ref::tridiagonal_eigen(d, e, n);
    REQUIRE(t.has_value());
    for (std::size_t k = 0; k < n; ++k) {
        CAPTURE(k, t->values[k], j.eigenvalues[k]);
        CHECK(std::abs(t->values[k] - j.eigenvalues[k]) <= 1e-13 + j.off_norm);
        CHECK(t->residuals[k] <= 1e-13);
    }
}

TEST_CASE("tridiagonal_eigen and schrodinger_1d refuse bad input", "[verify][schrodinger_poisson]") {
    const std::vector<double> d{1.0, 2.0}, e{0.5}, bad{0.5, 0.5};
    CHECK_FALSE(ref::tridiagonal_eigen(d, bad, 1).has_value());
    CHECK_FALSE(ref::tridiagonal_eigen(d, e, 0).has_value());
    CHECK_FALSE(ref::tridiagonal_eigen(d, e, 3).has_value());
    const std::vector<double> z{0.0, 1.0, 2.0, 3.0}, V(3, 0.0), m(3, 1.0), m0v{1.0, 0.0, 1.0}, zbad{0.0, 1.0, 1.0, 2.0};
    CHECK(ref::schrodinger_1d(z, V, m, 1.0, 2).has_value());
    CHECK_FALSE(ref::schrodinger_1d(z, V, m0v, 1.0, 1).has_value());
    CHECK_FALSE(ref::schrodinger_1d(zbad, V, m, 1.0, 1).has_value());
    CHECK_FALSE(ref::schrodinger_1d(z, V, m, 1.0, 3).has_value());
    CHECK_FALSE(ref::schrodinger_1d(z, V, m, -1.0, 1).has_value());
}

TEST_CASE("schrodinger_1d: the infinite well", "[verify][schrodinger_poisson]") {
    const double L = 10e-7, m = 0.5;
    double err_prev[3] = {0.0, 0.0, 0.0};
    for (const std::size_t cells : {100u, 200u, 400u, 800u}) {
        const auto z = uniform(0.0, L, cells);
        const std::vector<double> V(cells, 0.0), mass(cells, m);
        const auto s = ref::schrodinger_1d(z, V, mass, a_eVcm2, 3);
        REQUIRE(s.has_value());
        const double h = L / static_cast<double>(cells);
        for (int k = 1; k <= 3; ++k) {
            const double E = s->energies_eV[static_cast<std::size_t>(k - 1)];
            const double discrete = 2.0 * a_eVcm2 / (m * h * h) * (1.0 - std::cos(k * pi * h / L));
            const double exact = a_eVcm2 * pi * pi * k * k / (m * L * L);
            CAPTURE(cells, k, E, discrete, exact);
            CHECK(std::abs(E - discrete) <= 4.0 * std::numeric_limits<double>::epsilon() * 4.0 * a_eVcm2 / (m * h * h));  // eps ||A||
            const double err = std::abs(E - exact);
            if (err_prev[k - 1] > 0.0) {
                const double p = observed_order(err_prev[k - 1], err);
                CAPTURE(p);
                CHECK(p >= 1.9);
                CHECK(p <= 2.1);
            }
            err_prev[k - 1] = err;
        }
        // Normalisation: sum w psi^2 = 1.
        double norm = 0.0;
        for (std::size_t i = 1; i < cells; ++i) norm += 0.5 * (z[i + 1] - z[i - 1]) * s->psi[0][i] * s->psi[0][i];
        CHECK(std::abs(norm - 1.0) <= 1e-12);
    }
    // A smoothly graded grid, z = L (u - 0.3 sin(2 pi u) / (2 pi)), u uniform: still second order.
    double prev = 0.0;
    for (const std::size_t cells : {200u, 400u, 800u, 1600u}) {
        std::vector<double> z(cells + 1);
        for (std::size_t i = 0; i <= cells; ++i) {
            const double u = static_cast<double>(i) / static_cast<double>(cells);
            z[i] = L * (u - 0.3 * std::sin(2.0 * pi * u) / (2.0 * pi));
        }
        const std::vector<double> V(cells, 0.0), mass(cells, m);
        const auto s = ref::schrodinger_1d(z, V, mass, a_eVcm2, 1);
        REQUIRE(s.has_value());
        const double err = std::abs(s->energies_eV[0] - a_eVcm2 * pi * pi / (m * L * L));
        if (prev > 0.0) {
            const double p = observed_order(prev, err);
            CAPTURE(cells, p);
            CHECK(p >= 1.8);
            CHECK(p <= 2.2);
        }
        prev = err;
    }
}

TEST_CASE("schrodinger_1d: a finite well with a mass step, against Brent", "[verify][schrodinger_poisson]") {
    const double mw = 0.2, mb = 0.6, V0 = 0.25, W = 6e-7, B = 8e-7;  // well width W, barriers B each side
    // Even states: (k / m_w) tan(k W / 2) = kappa / m_b; odd: -(k / m_w) cot(k W / 2) = kappa / m_b.
    auto kk = [&](double E) { return std::sqrt(E * mw / a_eVcm2); };
    auto kap = [&](double E) { return std::sqrt((V0 - E) * mb / a_eVcm2); };
    auto even = [&](const ref::DD& E) {
        const double e = E.to_double(), k = kk(e);
        return ref::DD(k / mw * std::tan(k * W / 2.0) - kap(e) / mb);
    };
    auto odd = [&](const ref::DD& E) {
        const double e = E.to_double(), k = kk(e);
        return ref::DD(-k / mw / std::tan(k * W / 2.0) - kap(e) / mb);
    };
    // Energies where k W / 2 = theta: E = a (2 theta / W)^2 / m_w.
    auto at = [&](double theta) { return a_eVcm2 * std::pow(2.0 * theta / W, 2) / mw; };
    const auto e0 = ref::brent_zero(even, ref::DD(1e-9), ref::DD(at(pi / 2 - 1e-9)), 1e-30);
    const auto e1 = ref::brent_zero(odd, ref::DD(at(pi / 2 + 1e-9)), ref::DD(std::min(at(pi - 1e-9), V0 - 1e-9)),
                                    1e-30);
    REQUIRE(e0.status == ref::BrentStatus::ok);
    REQUIRE(e1.status == ref::BrentStatus::ok);
    const double exact[2] = {e0.root.to_double(), e1.root.to_double()};
    double prev[2] = {0.0, 0.0}, last[2] = {0.0, 0.0}, richardson[2] = {0.0, 0.0};
    for (const std::size_t cells : {220u, 440u, 880u, 1760u}) {
        const double total = W + 2.0 * B;
        const auto z = uniform(0.0, total, cells);
        std::vector<double> V(cells), mass(cells);
        for (std::size_t c = 0; c < cells; ++c) {
            const double mid = 0.5 * (z[c] + z[c + 1]);
            const bool well = mid > B && mid < B + W;
            V[c] = well ? 0.0 : V0;
            mass[c] = well ? mw : mb;
        }
        const auto s = ref::schrodinger_1d(z, V, mass, a_eVcm2, 2);
        REQUIRE(s.has_value());
        for (int k = 0; k < 2; ++k) {
            const double err = std::abs(s->energies_eV[static_cast<std::size_t>(k)] - exact[k]);
            CAPTURE(cells, k, s->energies_eV[static_cast<std::size_t>(k)], exact[k]);
            if (cells == 1760u) richardson[k] = (4.0 * s->energies_eV[static_cast<std::size_t>(k)] - last[k]) / 3.0;
            last[k] = s->energies_eV[static_cast<std::size_t>(k)];
            if (prev[k] > 0.0) {
                const double p = observed_order(prev[k], err);
                CAPTURE(p);
                CHECK(p >= 1.9);
                CHECK(p <= 2.1);
            }
            prev[k] = err;
        }
    }
}

TEST_CASE("schrodinger_1d: the triangular well against DLMF's Airy zeros", "[verify][schrodinger_poisson]") {
    const double airy[5] = {2.3381074105, 4.0879494441, 5.5205598281, 6.7867080901, 7.9441335871};  // -a_k, Table 9.9.1
    const double m = 0.916, F = 1e6, L = 20e-7;  // F in eV/cm
    const double E0 = std::cbrt(a_eVcm2 / m) * std::pow(F, 2.0 / 3.0);
    double prev[5] = {}, last[5] = {}, richardson[5] = {};
    for (const std::size_t cells : {500u, 1000u, 2000u, 4000u}) {
        const auto z = uniform(0.0, L, cells);
        std::vector<double> V(cells), mass(cells, m);
        for (std::size_t c = 0; c < cells; ++c) V[c] = F * 0.5 * (z[c] + z[c + 1]);
        const auto s = ref::schrodinger_1d(z, V, mass, a_eVcm2, 5);
        REQUIRE(s.has_value());
        for (std::size_t k = 0; k < 5; ++k) {
            const double exact = E0 * airy[k];
            const double err = std::abs(s->energies_eV[k] - exact);
            CAPTURE(cells, k, s->energies_eV[k], exact);
            if (cells == 4000u) richardson[k] = (4.0 * s->energies_eV[k] - last[k]) / 3.0;
            last[k] = s->energies_eV[k];
            if (prev[k] > 0.0) {
                const double p = observed_order(prev[k], err);
                CAPTURE(p);
                CHECK(p >= 1.9);
                CHECK(p <= 2.1);
            }
            prev[k] = err;
        }
    }
    // Richardson from the two finest grids (order 2 checked above) against the printed zeros (10 digits).
    for (std::size_t k = 0; k < 5; ++k) {
        CAPTURE(k, richardson[k], E0 * airy[k]);
        CHECK(std::abs(richardson[k] - E0 * airy[k]) <= 1e-8 * E0 * airy[k]);
    }
}

namespace {

ref::MosSpInput mos_input(ref::SpElectrons electrons) {
    ref::MosSpInput in;
    in.q_C = q_e;
    in.kT_eV = kT;
    in.eps0_F_per_cm = eps0;
    in.a_eV_cm2 = a_eVcm2;
    in.eps_si_r = 11.7;
    in.eps_ox_r = 3.9;
    in.band_gap_eV = 1.12;
    in.conduction_dos_cm3 = 2.8e19;
    in.valence_dos_cm3 = 1.04e19;
    in.oxide_barrier_eV = 3.0;
    in.oxide_mass = 0.5;
    in.valleys = {{.mass_z = 0.916, .mass_dos = 0.19, .degeneracy = 2.0},
                  {.mass_z = 0.19, .mass_dos = std::sqrt(0.916 * 0.19), .degeneracy = 4.0}};
    in.acceptors_cm3 = 1e17;
    in.oxide_thickness_cm = 5e-7;
    in.oxide_cells = 100;
    in.h0_cm = 2e-9;
    in.ratio = 1.03;
    in.h_max_cm = 1e-7;
    in.depth_cm = 1e-4;
    in.quantum_depth_cm = 40e-7;
    in.electrons = electrons;
    in.subbands_per_valley = 8;
    return in;
}

}  // namespace

TEST_CASE("mos_schrodinger_poisson: classical electrons against mos_capacitor", "[verify][schrodinger_poisson]") {
    const auto in = mos_input(ref::SpElectrons::classical);
    const double ni = std::sqrt(in.conduction_dos_cm3 * in.valence_dos_cm3) * std::exp(-in.band_gap_eV / (2.0 * kT));
    const ref::MosInput mi{.q_C = q_e,
                           .eps_F_per_cm = in.eps_si_r * eps0,
                           .thermal_voltage_V = kT,
                           .ni_cm3 = ni,
                           .net_doping_cm3 = -in.acceptors_cm3,
                           .oxide_capacitance_F_per_cm2 = in.eps_ox_r * eps0 / in.oxide_thickness_cm,
                           .flat_band_V = 0.0};
    auto fine = in;  // the grid halved
    fine.oxide_cells *= 2;
    fine.h0_cm /= 2.0;
    fine.ratio = std::sqrt(in.ratio);
    fine.h_max_cm /= 2.0;
    for (const double vg : {-1.0, 0.5, 1.5, 2.5}) {
        const auto r = ref::mos_schrodinger_poisson(in, vg);
        const auto rf = ref::mos_schrodinger_poisson(fine, vg);
        REQUIRE(r.has_value());
        REQUIRE(rf.has_value());
        const auto m = ref::mos_point(mi, vg);
        REQUIRE(m.has_value());
        const double qs = m->charge_C_per_cm2.to_double(), bend = m->band_bending_V.to_double();
        const double eq = std::abs(r->charge_C_per_cm2 - qs), ef = std::abs(rf->charge_C_per_cm2 - qs);
        const double qr = (4.0 * rf->charge_C_per_cm2 - r->charge_C_per_cm2) / 3.0;
        const double br = (4.0 * rf->surface_potential_V - r->surface_potential_V) / 3.0;
        CAPTURE(vg, r->charge_C_per_cm2, rf->charge_C_per_cm2, qr, qs, r->surface_potential_V, br, bend);
        // Second order in the grid toward the exact charge, and the extrapolated values against it.
        if (eq > 1e-9 * std::abs(qs)) {
            const double p = observed_order(eq, ef);
            CAPTURE(p);
            CHECK(p >= 1.8);
            CHECK(p <= 2.2);
        }
        CHECK(std::abs(qr - qs) <= 1e-7 * std::abs(qs));
        CHECK(std::abs(br - bend) <= 1e-7);
        for (const auto* x : {&*r, &*rf}) {
            CHECK(std::abs(x->oxide_displacement_C_per_cm2 + x->charge_C_per_cm2) <= 1e-9 * std::abs(qs));
            CHECK(x->poisson_residual <= 1e-10);
        }
    }
}

TEST_CASE("mos_schrodinger_poisson: the self-consistent inversion layer", "[verify][schrodinger_poisson]") {
    const double vg = 2.0;
    const auto base = mos_input(ref::SpElectrons::quantum);
    const auto r = ref::mos_schrodinger_poisson(base, vg);
    REQUIRE(r.has_value());
    const auto c = ref::mos_schrodinger_poisson(mos_input(ref::SpElectrons::classical), vg);
    REQUIRE(c.has_value());
    CAPTURE(r->outer_iterations, r->last_update_V, r->poisson_residual, r->orthonormality, r->max_eigen_residual_eV);
    CHECK(r->poisson_residual <= 1e-9);
    CHECK(std::abs(r->oxide_displacement_C_per_cm2 + r->charge_C_per_cm2) <= 1e-8 * std::abs(r->charge_C_per_cm2));
    CHECK(r->orthonormality <= 1e-10);
    CHECK(r->max_eigen_residual_eV <= 1e-9);
    // Confinement: less inversion charge and a centroid further from the interface than classically.
    CAPTURE(r->inversion_cm2, c->inversion_cm2, r->centroid_cm * 1e7, c->centroid_cm * 1e7);
    CHECK(r->inversion_cm2 < c->inversion_cm2);
    CHECK(r->centroid_cm > c->centroid_cm);
    // The ground subband lies above the silicon conduction band edge at the interface.
    const double ecb = base.band_gap_eV - kT * std::log(base.valence_dos_cm3 / base.acceptors_cm3);
    CHECK(r->subband_eV[0][0] > ecb - r->surface_potential_V);

    auto relchange = [&](const ref::MosSpInput& in) {
        const auto s = ref::mos_schrodinger_poisson(in, vg);
        REQUIRE(s.has_value());
        return std::pair{std::abs(s->inversion_cm2 / r->inversion_cm2 - 1.0),
                         std::abs(s->centroid_cm - r->centroid_cm)};
    };
    auto fine = base;  // the grid halved
    fine.oxide_cells *= 2;
    fine.h0_cm /= 2.0;
    fine.ratio = std::sqrt(base.ratio);
    fine.h_max_cm /= 2.0;
    const auto [dn_grid, dz_grid] = relchange(fine);
    auto more = base;
    more.subbands_per_valley = 12;
    const auto [dn_sub, dz_sub] = relchange(more);
    auto deeper = base;
    deeper.quantum_depth_cm = 60e-7;
    const auto [dn_depth, dz_depth] = relchange(deeper);
    CAPTURE(dn_grid, dz_grid * 1e7, dn_sub, dz_sub * 1e7, dn_depth, dz_depth * 1e7);
    CHECK(dn_grid <= 1e-3);
    CHECK(dz_grid <= 1e-9);  // 0.01 nm
    CHECK(dn_sub <= 1e-4);
    CHECK(dn_depth <= 1e-4);

    auto fd = base;
    fd.statistics = ref::SpStatistics::fermi_dirac;
    const auto f = ref::mos_schrodinger_poisson(fd, vg);
    REQUIRE(f.has_value());
    CHECK(f->poisson_residual <= 1e-9);
    CHECK(f->inversion_cm2 <= r->inversion_cm2);
    std::printf("SP reference, 1e17, 5 nm SiO2, V_G - V_FB = 2 V: N_inv %.6e (classical %.6e) cm^-2, centroid %.4f "
                "(classical %.4f) nm, %d outer iterations; FD N_inv %.6e\n",
                r->inversion_cm2, c->inversion_cm2, r->centroid_cm * 1e7, c->centroid_cm * 1e7, r->outer_iterations,
                f->inversion_cm2);
}
