// Unit 68 slice 2 (MAC10; DECISIONS.md amendment 6): the AlGaN/GaN 2DEG against Ambacher's
// analytic sheet density, and AlGaN/GaN heterojunctions under both interface transports (owner
// test 11).
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <numbers>
#include <string>

#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/physics/nitrides.hpp"
#include "NiTCAD/physics/polarization.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;

namespace {

// Ambacher et al. (M3; J. Appl. Phys. 85, 3222 (1999), eqs (11)-(13)): the 2DEG of a pseudomorphic
// Al(x)Ga(1-x)N barrier of thickness d under a gate of barrier phi_b,
//     n_s = sigma / e - eps0 eps_b / (d e^2) [e phi_b + E_F - dEc],
//     E_F = E_0 + pi hbar^2 n_s / m*,  E_0 = (9 pi hbar e^2 n_s / (8 eps0 eps_w sqrt(8 m*)))^(2/3),
// with this model's sigma, permittivities, dEc and GaN mass [cm^-2].
double ambacher_ns(double x, double d_cm, double phi_b) {
    const double hbar = 1.054571817e-34, e = 1.602176634e-19, eps0 = 8.8541878128e-12, m0 = 9.1093837015e-31;
    const physics::SemiconductorParameters g = physics::gan_parameters(), b = *physics::algan_parameters(x);
    const double a = physics::nitride_lattice("gan", 0)->a0_A;
    const double sigma = physics::polarization_sheet_density_cm2(
        physics::nitride_polarization(*physics::nitride_lattice("gan", 0), a) -
        physics::nitride_polarization(*physics::nitride_lattice("algan", x), a));
    const double dEc = g.electron_affinity_eV - b.electron_affinity_eV;
    const double m = 0.20 * m0;  // GaN, amendment 6
    double ns = sigma;
    for (int k = 0; k < 200; ++k) {
        const double n = ns * 1e4;  // m^-2
        const double E0 = std::pow(9 * std::numbers::pi * hbar * e * e * n / (8 * eps0 * g.eps_r * std::sqrt(8 * m)), 2.0 / 3.0);
        const double EF = (E0 + std::numbers::pi * hbar * hbar * n / m) / e;  // eV
        ns = 0.5 * ns + 0.5 * (sigma - eps0 * b.eps_r / (d_cm * 1e-2 * e) * (phi_b + EF - dEc) * 1e-4);
    }
    return ns;
}

// A 1D Al0.25Ga0.75N (25 nm) / GaN (1 um) stack, the surface at x_max.
std::string stack(const std::string& top, const std::string& extra) {
    return R"([mesh]
type = tensor
x = linspace(0, 0.9um, 46) + linspace(0.9um, 1.025um, 251)

[region buffer]
material = gan

[region barrier]
material = algan 0.25
box = 1um 1.025um

[doping]
profile = donor constant 1e15

[contact body]
kind = ohmic
patch = x_min

[contact top]
patch = x_max
)" + top + R"(
[models]
polarization = on
fermi_dirac = on
)" + extra;
}

// The electrons in the GaN below the interface [cm^-2].
double sheet(const deck::Plan& p, const results::NodeFields& f) {
    const auto pts = p.device.mesh().points();
    const auto vol = p.device.mesh().volumes();
    double ns = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (pts[i][0] < 1e-4) ns += f.n_cm3[i] * vol[i];
    }
    return ns - 1e15 * 1e-4;  // less the background
}

}  // namespace

TEST_CASE("MAC10: the 1D gated Al0.25Ga0.75N/GaN 2DEG against Ambacher", "[deck][unit68][hemt]") {
    const deck::Plan p = must(stack("kind = schottky\nbarrier = 0.83\n", "\n[equilibrium]\n"));
    const auto r = deck::run(p);
    if (!r) FAIL(r.error().message);
    const double ns = sheet(p, r->equilibrium->fields), ref = ambacher_ns(0.25, 25e-7, 0.83);
    std::printf("MAC10 1D: n_s %.4e cm^-2, Ambacher %.4e, off by %+.1f %%\n", ns, ref, 100.0 * (ns / ref - 1.0));
    CHECK(std::abs(ns / ref - 1.0) <= 0.10);
}

TEST_CASE("owner test 11: AlGaN/GaN under continuous and thermionic transport", "[deck][unit68][hemt]") {
    for (const char* transport : {"drift_diffusion", "thermionic"}) {
        std::string text = stack("kind = ohmic\n", std::string("\n[interface buffer barrier]\ntransport = ") + transport +
                                                       "\n\n[sweep]\ncontact top = 0 : 0.4 : 0.1\n");
        text.replace(text.find("donor constant 1e15"), 19, "donor constant 1e17");
        const deck::Plan p = must(text);
        const auto r = deck::run(p);
        if (!r) FAIL(r.error().message);
        REQUIRE(r->sweep.has_value());
        REQUIRE_FALSE(r->sweep->stopped.has_value());
        const auto& pts = r->sweep->points;
        // Equilibrium: no current at zero bias; conservation: the two terminals' currents cancel.
        const double scale = std::abs(pts.back().terminal_current[0]);
        CHECK(scale > 0.0);
        CHECK(std::abs(pts.front().terminal_current[0]) < 1e-9 * scale + 1e-15);
        for (const auto& b : pts) {
            CHECK(std::abs(b.terminal_current[0] + b.terminal_current[1]) < 1e-8 * scale);
        }
        std::printf("test 11 %s: I(0.4 V) %.4e A/cm^2, sum %.1e\n", transport, pts.back().terminal_current[1],
                    pts.back().terminal_current[0] + pts.back().terminal_current[1]);
    }
}

TEST_CASE("MAC10: the library HEMT's 2DEG under the gate against Ambacher", "[deck][unit68][hemt]") {
    std::ifstream in(std::filesystem::path(NITCAD_EXAMPLES_DIR) / "devices" / "algan_gan_hemt_2d.deck");
    std::stringstream s;
    s << in.rdbuf();
    std::string text = s.str();
    REQUIRE(text.find("[sweep]") != std::string::npos);
    text = text.substr(0, text.find("[sweep]")) + "[equilibrium]\n";
    const deck::Plan p = must(text);
    const auto r = deck::run(p);
    if (!r) FAIL(r.error().message);
    // The electrons in the buffer's column at the gate's middle (x = 1.5 um), per area.
    const auto pts = p.device.mesh().points();
    const auto vol = p.device.mesh().volumes();
    const double dx = 3e-4 / 40;
    double ns = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (std::abs(pts[i][0] - 1.5e-4) < 1e-9 && pts[i][1] < 1e-4) ns += r->equilibrium->fields.n_cm3[i] * vol[i] / dx;
    }
    const double ref = ambacher_ns(0.25, 25e-7, 0.83);
    std::printf("MAC10 2D HEMT: n_s %.4e cm^-2, Ambacher %.4e, off by %+.1f %%\n", ns, ref, 100.0 * (ns / ref - 1.0));
    CHECK(std::abs(ns / ref - 1.0) <= 0.10);
}

// Unit 68 slice 2 (owner 2026-11-27, lifting Unit 48 ruling 4's refusal): a Schottky contact on a device
// with a thermionic heterointerface. Zero current at zero bias; the terminal currents cancel at every bias.
TEST_CASE("Schottky contact with a thermionic AlGaN/GaN interface", "[deck][unit68][hemt]") {
    const auto check = [](const std::string& text, const char* name) {
        INFO(name);
        const deck::Plan p = must(text);
        const auto r = deck::run(p);
        if (!r) FAIL(r.error().message);
        REQUIRE(r->sweep.has_value());
        REQUIRE_FALSE(r->sweep->stopped.has_value());
        double scale = 0.0;
        for (const auto& b : r->sweep->points) {
            for (const double i : b.terminal_current) scale = std::max(scale, std::abs(i));
        }
        REQUIRE(scale > 0.0);
        for (const auto& b : r->sweep->points) {
            double sum = 0.0;
            for (const double i : b.terminal_current) sum += i;
            CHECK(std::abs(sum) < 1e-8 * scale);
        }
        return r->sweep->points;
    };
    // 1D: the gated stack, the Schottky gate swept forward.
    const auto one = check(stack("kind = schottky\nbarrier = 0.83\n",
                                 "\n[interface buffer barrier]\ntransport = thermionic\n\n[sweep]\ncontact top = 0 : 0.6 : 0.1\n"),
                           "1d");
    CHECK(std::abs(one.front().terminal_current[0]) < 1e-9 * std::abs(one.back().terminal_current[0]) + 1e-15);
    std::printf("Schottky+thermionic 1D: I(0.6 V) %.4e A/cm^2\n", one.back().terminal_current[1]);

    // 2D: the library HEMT (its interface thermionic) against the same with the continuous interface.
    std::ifstream in(std::filesystem::path(NITCAD_EXAMPLES_DIR) / "devices" / "algan_gan_hemt_2d.deck");
    std::stringstream s;
    s << in.rdbuf();
    const std::string deck = s.str();
    const std::string sweep = "[sweep]\ncontact drain = 0.1\ncontact gate = -4 : 0 : 1\n";
    const std::string head = deck.substr(0, deck.find("[sweep]"));
    REQUIRE(head.find("transport = thermionic") != std::string::npos);
    std::string continuous = head;
    continuous.replace(continuous.find("transport = thermionic"), 22, "transport = drift_diffusion");
    const auto dd = check(continuous + sweep, "hemt_dd");
    const auto te = check(head + sweep, "hemt_te");
    const double id_dd = dd.back().terminal_current[1], id_te = te.back().terminal_current[1];
    std::printf("HEMT Id(Vg 0, Vd 0.1): continuous %.4e, thermionic %.4e A/cm (%+.2f %%)\n", id_dd, id_te,
                100.0 * (id_te / id_dd - 1.0));
    CHECK(std::abs(id_te / id_dd - 1.0) < 0.05);  // along the channel the interface carries little current
}
