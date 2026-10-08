// Band-to-band tunnelling functions (ARCHITECTURE.md section 11, Unit 20): the Kane rate with the
// legacy silicon pair, the direct-gap WKB functions of Esseni et al. (2017) (kappa vanishing at the
// band edges, the antiderivatives against quadrature, the exact segment integrals and their
// partials, eq. (11) reducing to eq. (8) in a uniform field), and the parameter validation.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/band_to_band.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;
namespace base = NiTCAD::base;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// A synthetic direct gap for the WKB functions: 0.8 eV, m_r = 0.05 m0 (the masses of a test, not
// of any material).
const WkbBand band{0.8 * base::q_C, 0.05 * base::m0_kg};

// Composite Simpson on [a, b] with m (even) intervals.
template <class F>
double simpson(F f, double a, double b, int m) {
    const double h = (b - a) / m;
    double s = f(a) + f(b);
    for (int k = 1; k < m; ++k) s += (k % 2 == 1 ? 4.0 : 2.0) * f(a + k * h);
    return s * h / 3.0;
}

}  // namespace

TEST_CASE("band to band: the material sets") {
    const BandToBandParameters& si = silicon_parameters.band_to_band;
    REQUIRE(si.A_per_cm3_s == 3.5e21);
    REQUIRE(si.B_V_per_cm == 1.03e8);
    REQUIRE_FALSE(si.direct_gap);
    REQUIRE(tunnelling_reduced_mass(si) == 0.0);  // silicon never takes the direct-gap rate
    for (const auto* p : {&gallium_arsenide_parameters, &indium_gallium_arsenide_parameters}) {
        REQUIRE(p->band_to_band.direct_gap);
        REQUIRE(p->band_to_band.A_per_cm3_s == 0.0);
        REQUIRE(tunnelling_reduced_mass(p->band_to_band) == 0.0);  // no verified masses
    }
    REQUIRE(algaas_parameters(0.3)->band_to_band.direct_gap);
    for (const auto* p : {&germanium_parameters, &silicon_carbide_4h_parameters}) {
        REQUIRE(p->band_to_band == BandToBandParameters{});
    }
}

TEST_CASE("band to band: the Kane rate and its derivative") {
    const double A = 3.5e21, B = 1.03e8;
    for (const double F : {1e5, 5e5, 1e6, 3e6}) {
        CAPTURE(F);
        const KaneRate r = kane_generation(A, B, F);
        REQUIRE(close(r.rate, A * F * F * std::exp(-B / F), 1e-15));
        const double h = 1e-6 * F;
        const double fd =
            (kane_generation(A, B, F + h).rate - kane_generation(A, B, F - h).rate) / (2.0 * h);
        REQUIRE(close(r.d_dF, fd, 1e-7));
    }
    // Zero at F = 0, for A = 0, and below the underflow of exp(-B / F).
    REQUIRE(kane_generation(A, B, 0.0).rate == 0.0);
    REQUIRE(kane_generation(0.0, B, 1e6).rate == 0.0);
    const KaneRate low = kane_generation(A, B, B / 701.0);
    REQUIRE(low.rate == 0.0);
    REQUIRE(low.d_dF == 0.0);
}

TEST_CASE("band to band: kappa vanishes at the band edges; antiderivatives against quadrature") {
    REQUIRE(wkb_u(band) == 10.0);
    REQUIRE(wkb_kappa(0.0, band) == 0.0);
    REQUIRE(wkb_kappa(1.0, band) == 0.0);
    REQUIRE(wkb_kappa(-0.5, band) == 0.0);  // outside the gap: an allowed band
    REQUIRE(wkb_kappa(1e-12, band) > 0.0);  // no cancellation near the edge
    REQUIRE(wkb_kappa(1.0 - 1e-12, band) > 0.0);
    const double s = std::sqrt(band.reduced_mass_kg * band.gap_J);
    const double pi = std::numbers::pi;
    // Over the whole gap.
    REQUIRE(close(wkb_kappa_antiderivative(1.0, band) - wkb_kappa_antiderivative(0.0, band),
                  pi * s / (4.0 * base::hbar_J_s), 1e-14));
    REQUIRE(close(wkb_inverse_kappa_antiderivative(1.0, band) -
                      wkb_inverse_kappa_antiderivative(0.0, band),
                  pi * base::hbar_J_s / (2.0 * s), 1e-14));
    // Inside, against Simpson.
    for (const auto& [a, b] : {std::pair{0.1, 0.9}, std::pair{0.02, 0.3}, std::pair{0.6, 0.97}}) {
        CAPTURE(a, b);
        const double k = simpson([](double d) { return wkb_kappa(d, band); }, a, b, 4000);
        const double i = simpson([](double d) { return 1.0 / wkb_kappa(d, band); }, a, b, 4000);
        REQUIRE(close(wkb_kappa_antiderivative(b, band) - wkb_kappa_antiderivative(a, band), k,
                      1e-10));
        REQUIRE(close(wkb_inverse_kappa_antiderivative(b, band) -
                          wkb_inverse_kappa_antiderivative(a, band),
                      i, 1e-9));
    }
}

TEST_CASE("band to band: segment integrals, their partials and continuity") {
    const double L = 3e-9;  // m
    const auto check = [&](double da, double db) {
        CAPTURE(da, db);
        const WkbSegment s = wkb_segment(da, db, L, band);
        const double h = 1e-7;
        const WkbSegment pa = wkb_segment(da + h, db, L, band);
        const WkbSegment ma = wkb_segment(da - h, db, L, band);
        const WkbSegment pb = wkb_segment(da, db + h, L, band);
        const WkbSegment mb = wkb_segment(da, db - h, L, band);
        const auto near = [](double d, double fd, double ref) {
            return std::abs(d - fd) <= 1e-6 * std::abs(ref) + 1e-30;
        };
        REQUIRE(near(s.dI_k_da, (pa.I_k - ma.I_k) / (2 * h), s.I_k / std::abs(db - da)));
        REQUIRE(near(s.dI_k_db, (pb.I_k - mb.I_k) / (2 * h), s.I_k / std::abs(db - da)));
        REQUIRE(near(s.dI_ik_da, (pa.I_ik - ma.I_ik) / (2 * h), s.I_ik / std::abs(db - da)));
        REQUIRE(near(s.dI_ik_db, (pb.I_ik - mb.I_ik) / (2 * h), s.I_ik / std::abs(db - da)));
    };
    check(0.2, 0.5);
    check(0.7, 0.3);   // falling
    check(-0.3, 0.4);  // entering the gap
    check(0.6, 1.4);   // leaving it
    // Entirely outside the gap: nothing.
    const WkbSegment out = wkb_segment(1.2, 1.7, L, band);
    REQUIRE(out.I_k == 0.0);
    REQUIRE(out.I_ik == 0.0);
    // Continuous as an end crosses the band edge, and as the segment flattens.
    const WkbSegment a = wkb_segment(0.4, 1.0 - 1e-9, L, band);
    const WkbSegment b = wkb_segment(0.4, 1.0 + 1e-9, L, band);
    REQUIRE(close(a.I_k, b.I_k, 1e-8));
    // 1 / kappa ~ (1 - delta)^(-1/2) at the edge: the part beyond 1 - 1e-9 is of order sqrt(1e-9).
    REQUIRE(close(a.I_ik, b.I_ik, 1e-4));
    const WkbSegment flat = wkb_segment(0.5, 0.5 + 0.5e-10, L, band);
    const WkbSegment steep = wkb_segment(0.5, 0.5 + 2e-10, L, band);
    // Just above the flat threshold the antiderivatives' difference carries eps / D ~ 5e-7.
    REQUIRE(close(flat.I_k, steep.I_k, 1e-5));
    REQUIRE(close(flat.I_k, L * wkb_kappa(0.5, band), 1e-8));
    // The flat branch's partials are those of L f(midpoint), each end taking half (the exact
    // branch's are rounding there, eps / D^2).
    const double m = 0.5 + 0.25e-10, h = 1e-6;
    const double dk = (wkb_kappa(m + h, band) - wkb_kappa(m - h, band)) / (2 * h);
    const double k = wkb_kappa(m, band);
    REQUIRE(close(flat.dI_k_da, 0.5 * L * dk, 1e-6));
    REQUIRE(close(flat.dI_k_db, 0.5 * L * dk, 1e-6));
    REQUIRE(close(flat.dI_ik_da, -0.5 * L * dk / (k * k), 1e-6));
}

TEST_CASE("band to band: eq. (11) on a uniform field is eq. (8) at any resolution") {
    const double pi = std::numbers::pi;
    for (const double F : {1e8, 3e8, 1e9}) {  // V/m
        const double W = band.gap_J / (base::q_C * F);  // tunnel length
        for (const int segments : {1, 3, 17, 200}) {
            CAPTURE(F, segments);
            // delta runs linearly from 0 to 1 over W, then on to 1.4 (outside the gap).
            double Ik = 0.0, Iik = 0.0;
            const int total = segments + segments / 3 + 1;
            const double h = 1.4 * W / total;
            for (int k = 0; k < total; ++k) {
                const WkbSegment s =
                    wkb_segment(k * h / W, (k + 1) * h / W, h, band);
                Ik += s.I_k;
                Iik += s.I_ik;
            }
            REQUIRE(close(2.0 * Ik,
                          pi * std::sqrt(band.reduced_mass_kg) * std::pow(band.gap_J, 1.5) /
                              (2.0 * base::q_C * F * base::hbar_J_s),
                          1e-12));
            // A wide energy window (k_m^2 I_ik >> 1): the bracket is 1.
            const WkbPathRate r = wkb_path_rate(base::q_C * F, Ik, Iik, 1e40);
            REQUIRE(close(r.rate, wkb_uniform_field_rate(F, band), 1e-11));
        }
    }
}

TEST_CASE("band to band: the path rate's partials") {
    const double slope = 2e-11, Ik = 3.0, Iik = 4e-17, km2 = 5e16;
    const WkbPathRate r = wkb_path_rate(slope, Ik, Iik, km2);
    const auto fd = [&](int which) {
        double v[4] = {slope, Ik, Iik, km2};
        const double h = 1e-6 * v[which];
        v[which] += h;
        const double up = wkb_path_rate(v[0], v[1], v[2], v[3]).rate;
        v[which] -= 2 * h;
        const double down = wkb_path_rate(v[0], v[1], v[2], v[3]).rate;
        return (up - down) / (2 * h);
    };
    REQUIRE(close(r.d_slope, fd(0), 1e-7));
    REQUIRE(close(r.d_I_k, fd(1), 1e-7));
    REQUIRE(close(r.d_I_ik, fd(2), 1e-7));
    REQUIRE(close(r.d_km2, fd(3), 1e-7));
    REQUIRE(wkb_path_rate(slope, Ik, 0.0, km2).rate == 0.0);  // no barrier
}

TEST_CASE("band to band: parameter validation") {
    const auto rejected = [](BandToBandParameters b, const char* name) {
        SemiconductorParameters p = gallium_arsenide_parameters;
        p.band_to_band = b;
        const auto m = Semiconductor::create(p);
        REQUIRE_FALSE(m.has_value());
        REQUIRE(m.error().code == base::ErrorCode::invalid_input);
        REQUIRE(m.error().message.find(name) != std::string::npos);
    };
    rejected({-1.0, 1e8, true, 0.0, 0.0}, "band_to_band.A_per_cm3_s");
    rejected({1e21, 0.0, true, 0.0, 0.0}, "band_to_band.B_V_per_cm");
    rejected({0.0, 0.0, true, 0.06, 0.0}, "band_to_band.hole_mass");
    rejected({0.0, 0.0, false, 0.06, 0.4}, "direct gap");
    rejected({0.0, 0.0, true, 1.5, 1.5}, "reduced mass");
    rejected({0.0, 0.0, true, -0.1, 0.4}, "band_to_band.electron_mass");
    // A direct gap with both masses (a user's cited values) is accepted.
    SemiconductorParameters p = gallium_arsenide_parameters;
    p.band_to_band.electron_mass = 0.067;
    p.band_to_band.hole_mass = 0.5;
    REQUIRE(Semiconductor::create(p).has_value());
    REQUIRE(close(tunnelling_reduced_mass(p.band_to_band), 0.067 * 0.5 / 0.567, 1e-15));
}
