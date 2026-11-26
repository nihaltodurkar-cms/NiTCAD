// Unit 68 slice 2 (MAD8; DECISIONS.md amendment 6): the binary nitrides and the Caughey-Thomas
// temperature exponents.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>

#include "NiTCAD/physics/ionization.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/nitrides.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermal.hpp"

using namespace NiTCAD::physics;

namespace {
Semiconductor make(const SemiconductorParameters& p) {
    auto m = Semiconductor::create(p);
    if (!m) FAIL(m.error().message);
    return *m;
}
double mu(const Semiconductor& m, Carrier c, double N, double T) { return caughey_thomas_mobility(m, c, N, T); }
}  // namespace

TEST_CASE("nitrides: the binaries build, with their gaps and Q_c alignment", "[nitrides][unit68]") {
    const Semiconductor gan = make(gan_parameters()), aln = make(aln_parameters()), inn = make(inn_parameters());
    CHECK(std::abs(band_gap_eV(aln, 1e-6) - 6.10) < 1e-9);
    CHECK(std::abs(band_gap_eV(aln, 300.0) - (6.10 - 2.63e-3 * 9e4 / 2382.0)) < 1e-12);
    CHECK(std::abs(band_gap_eV(inn, 300.0) - 0.608) < 0.0005);
    CHECK(std::abs(band_gap_eV(gan, 300.0) - (3.47 - 7.7e-4 * 9e4 / 900.0)) < 1e-12);
    const double dEg_al = band_gap_eV(aln, 300.0) - band_gap_eV(gan, 300.0);
    const double dEg_in = band_gap_eV(inn, 300.0) - band_gap_eV(gan, 300.0);
    // Conduction step Q_c dEg, valence step the rest (option A).
    CHECK(std::abs((gan_parameters().electron_affinity_eV - aln_parameters().electron_affinity_eV) - 0.70 * dEg_al) < 1e-12);
    CHECK(std::abs((gan_parameters().electron_affinity_eV - inn_parameters().electron_affinity_eV) - 0.58 * dEg_in) < 1e-12);
    std::printf("nitrides 300 K: Eg GaN %.4f AlN %.4f InN %.4f eV; chi AlN %.5f InN %.5f eV\n", band_gap_eV(gan, 300.0),
                band_gap_eV(aln, 300.0), band_gap_eV(inn, 300.0), aln_parameters().electron_affinity_eV,
                inn_parameters().electron_affinity_eV);
    CHECK(inn_parameters().ionization.donor_eV == 0.0);  // fully ionized fallback
}

TEST_CASE("nitrides: electron mobility falls with doping (AlN's alpha sign)", "[nitrides][unit68]") {
    for (const SemiconductorParameters& p : {gan_parameters(), aln_parameters(), inn_parameters()}) {
        const Semiconductor m = make(p);
        for (const double T : {250.0, 300.0, 400.0}) {
            double last = mu(m, Carrier::electron, 0.0, T);
            for (const double N : {1e15, 1e16, 1e17, 1e18, 1e19, 1e20}) {
                const double v = mu(m, Carrier::electron, N, T);
                CHECK(v < last);
                last = v;
            }
        }
        CHECK(mu(m, Carrier::hole, 1e18, 350.0) == p.hole_mobility.mu_max);  // constant
    }
}

TEST_CASE("Caughey-Thomas exponents: each moves only its own parameter, slope exact", "[nitrides][unit68]") {
    SemiconductorParameters base = silicon_parameters;
    const Semiconductor si = make(base);
    const double N = 3e17, T = 420.0;
    // Exponents 0: the legacy value (every existing set).
    const auto& ct = base.electron_mobility;
    const double legacy = ct.mu_min + (ct.mu_max * std::pow(T / 300.0, ct.T_exponent) - ct.mu_min) /
                                          (1.0 + std::pow(N / ct.N_ref, ct.alpha));
    CHECK(mu(si, Carrier::electron, N, T) == legacy);
    const auto with = [&](auto set) {
        SemiconductorParameters p = base;
        set(p.electron_mobility);
        return make(p);
    };
    const Semiconductor a = with([](CaugheyThomasParameters& c) { c.mu_min_T_exponent = -1.0; });
    const Semiconductor b = with([](CaugheyThomasParameters& c) { c.N_ref_T_exponent = 1.5; });
    const Semiconductor c = with([](CaugheyThomasParameters& c) { c.alpha_T_exponent = 0.5; });
    const double t = T / 300.0;
    const auto ct_at = [&](double mn, double Nr, double al) {
        return mn + (ct.mu_max * std::pow(t, ct.T_exponent) - mn) / (1.0 + std::pow(N / Nr, al));
    };
    CHECK(std::abs(mu(a, Carrier::electron, N, T) - ct_at(ct.mu_min / t, ct.N_ref, ct.alpha)) < 1e-9);
    CHECK(std::abs(mu(b, Carrier::electron, N, T) - ct_at(ct.mu_min, ct.N_ref * std::pow(t, 1.5), ct.alpha)) < 1e-9);
    CHECK(std::abs(mu(c, Carrier::electron, N, T) - ct_at(ct.mu_min, ct.N_ref, ct.alpha * std::sqrt(t))) < 1e-9);
    for (const Semiconductor* m : {&a, &b, &c}) {
        for (const double n : {0.0, 1e15, 3e17, 1e20}) {
            const double h = 1e-3;
            const double fd = (mu(*m, Carrier::electron, n, T + h) - mu(*m, Carrier::electron, n, T - h)) / (2 * h);
            CHECK(std::abs(caughey_thomas_mobility_slope(*m, Carrier::electron, n, T) - fd) < 1e-6 * std::abs(fd) + 1e-9);
        }
    }
    const Semiconductor gan = make(gan_parameters());
    for (const double n : {0.0, 1e16, 1e18}) {
        const double h = 1e-3;
        const double fd = (mu(gan, Carrier::electron, n, 350 + h) - mu(gan, Carrier::electron, n, 350 - h)) / (2 * h);
        CHECK(std::abs(caughey_thomas_mobility_slope(gan, Carrier::electron, n, 350.0) - fd) < 1e-6 * std::abs(fd));
    }
}

namespace {
SemiconductorParameters alloy_of(bool al, double x) {
    auto p = al ? algan_parameters(x) : ingan_parameters(x);
    if (!p) FAIL(p.error().message);
    return *p;
}
double gap300(const SemiconductorParameters& p) { return band_gap_eV(make(p), 300.0); }
}  // namespace

TEST_CASE("nitride alloys: the binaries at the end points, continuous next to them", "[nitrides][unit68]") {
    CHECK(alloy_of(true, 0.0) == gan_parameters());
    CHECK(alloy_of(true, 1.0) == aln_parameters());
    CHECK(alloy_of(false, 0.0) == gan_parameters());
    CHECK(alloy_of(false, 1.0) == inn_parameters());
    for (const bool al : {true, false}) {
        const SemiconductorParameters b = al ? aln_parameters() : inn_parameters();
        const SemiconductorParameters lo = alloy_of(al, 1e-9), hi = alloy_of(al, 1.0 - 1e-9);
        CHECK(std::abs(lo.electron_affinity_eV - gan_parameters().electron_affinity_eV) < 1e-7);
        CHECK(std::abs(hi.electron_affinity_eV - b.electron_affinity_eV) < 1e-7);
        for (const double T : {200.0, 300.0, 450.0}) {
            CHECK(std::abs(band_gap_eV(make(lo), T) - band_gap_eV(make(gan_parameters()), T)) < 1e-7);
            if (al) CHECK(std::abs(band_gap_eV(make(hi), T) - band_gap_eV(make(b), T)) < 1e-7);
        }
        CHECK(std::abs(hi.Nc300 / b.Nc300 - 1.0) < 1e-7);
        CHECK(std::abs(hi.Nv300 / b.Nv300 - 1.0) < 1e-7);
        // Mobility next to the ends: the binaries' (owner test 8).
        for (const double N : {1e16, 1e18}) {
            CHECK(std::abs(mu(make(lo), Carrier::electron, N, 300) / mu(make(gan_parameters()), Carrier::electron, N, 300) - 1) < 1e-6);
            CHECK(std::abs(mu(make(hi), Carrier::electron, N, 300) / mu(make(b), Carrier::electron, N, 300) - 1) < 1e-6);
            CHECK(std::abs(mu(make(hi), Carrier::hole, N, 300) / mu(make(b), Carrier::hole, N, 300) - 1) < 1e-6);
        }
    }
    CHECK_FALSE(algan_parameters(1.01).has_value());
    CHECK_FALSE(ingan_parameters(-0.1).has_value());
}

TEST_CASE("nitride alloys: the bowed gap, the Q_c offset partition at 300 K", "[nitrides][unit68]") {
    const double g = gap300(gan_parameters());
    for (const double x : {0.15, 0.3, 0.7}) {
        const SemiconductorParameters al = alloy_of(true, x), in = alloy_of(false, x);
        CHECK(std::abs(gap300(al) - ((1 - x) * g + x * gap300(aln_parameters()) - 0.8 * x * (1 - x))) < 1e-12);
        CHECK(std::abs(gap300(in) - ((1 - x) * g + x * gap300(inn_parameters()) - 1.4 * x * (1 - x))) < 1e-12);
        // dEc = Q_c dEg; dEv the rest.
        CHECK(std::abs((4.10 - al.electron_affinity_eV) - 0.70 * (gap300(al) - g)) < 1e-12);
        CHECK(std::abs((4.10 - in.electron_affinity_eV) - 0.58 * (gap300(in) - g)) < 1e-12);
        // The gap's change from 300 K is linear in x.
        const auto shift = [](const SemiconductorParameters& p) { return band_gap_eV(make(p), 400.0) - gap300(p); };
        CHECK(std::abs(shift(al) - ((1 - x) * shift(gan_parameters()) + x * shift(aln_parameters()))) < 1e-12);
        // The affinity stays its 300 K value (v1).
        CHECK(al.electron_affinity_eV == alloy_of(true, x).electron_affinity_eV);
    }
    const SemiconductorParameters a30 = alloy_of(true, 0.3);
    std::printf("Al0.3Ga0.7N: Eg %.4f eV, dEc %.4f eV, dEv %.4f eV\n", gap300(a30), 4.10 - a30.electron_affinity_eV,
                (gap300(a30) - gap300(gan_parameters())) - (4.10 - a30.electron_affinity_eV));
}

TEST_CASE("nitride alloys: mobility dips inside, slope exact, Nv from the DOS mass", "[nitrides][unit68]") {
    const Semiconductor gan = make(gan_parameters()), aln = make(aln_parameters());
    const Semiconductor mid = make(alloy_of(true, 0.5));
    const double N = 1e17;
    CHECK(mu(mid, Carrier::electron, N, 300) < std::min(mu(gan, Carrier::electron, N, 300), mu(aln, Carrier::electron, N, 300)));
    const double expect = 1.0 / (0.5 / mu(gan, Carrier::electron, N, 300) + 0.5 / mu(aln, Carrier::electron, N, 300) + 0.25 / 40.0);
    CHECK(std::abs(mu(mid, Carrier::electron, N, 300) - expect) < 1e-9 * expect);
    for (const bool al : {true, false}) {
        const Semiconductor m = make(alloy_of(al, 0.3));
        for (const Carrier c : {Carrier::electron, Carrier::hole}) {
            for (const double n : {0.0, 1e17, 1e19}) {
                const double h = 1e-3;
                const double fd = (mu(m, c, n, 330 + h) - mu(m, c, n, 330 - h)) / (2 * h);
                CHECK(std::abs(caughey_thomas_mobility_slope(m, c, n, 330.0) - fd) < 1e-6 * std::abs(fd) + 1e-9);
            }
        }
    }
    // Nv (and Nc) from the DOS mass with the engine's constant; the alloy's Nv is m^(3/2) linear.
    const auto dos = [](double m) { return 4.83e15 * std::pow(m, 1.5) * std::pow(300.0, 1.5); };
    CHECK(std::abs(gan_parameters().Nv300 / dos(1.5) - 1) < 1e-14);
    CHECK(std::abs(aln_parameters().Nv300 / dos(2.7) - 1) < 1e-14);
    CHECK(std::abs(inn_parameters().Nv300 / 5.3e19 - 1) < 0.01);
    const double mh = std::pow(0.7 * std::pow(1.5, 1.5) + 0.3 * std::pow(2.7, 1.5), 2.0 / 3.0);
    CHECK(std::abs(alloy_of(true, 0.3).Nv300 / dos(mh) - 1) < 1e-12);
    const double me = 0.7 * 0.20 + 0.3 * 0.40 - 4.8e-3 * 0.21;
    CHECK(std::abs(alloy_of(true, 0.3).Nc300 / dos(me) - 1) < 1e-12);
}

TEST_CASE("nitrides: InN's full-ionization fallback, a configured level, its temperature ceiling", "[nitrides][unit68]") {
    // Depth 0: every donor ionized at any Fermi level; a level set (a [material] override) is not.
    const double kT = 0.025852;
    CHECK(ionized_density(1e18, 2.0, inn_parameters().ionization.donor_eV / kT, 2.0).value == 1e18);
    SemiconductorParameters p = inn_parameters();
    p.ionization.donor_eV = 0.05;
    const Semiconductor m = make(p);
    CHECK(ionized_density(1e18, 2.0, m.parameters().ionization.donor_eV / kT, 2.0).value < 1e18);
    const Semiconductor inn = make(inn_parameters());
    CHECK(check_temperature(inn, 480.0).has_value());
    CHECK_FALSE(check_temperature(inn, 500.0).has_value());
    // The alloys' end-member checks: InGaN inherits InN's ceiling only as x grows.
    CHECK(check_temperature(make(alloy_of(false, 0.2)), 480.0).has_value());
}
