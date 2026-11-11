// The equilibrium abrupt heterojunction, Boltzmann statistics: the exact solution of 1D Poisson-Boltzmann by the
// first integral on each side, matched at the interface (Unit 55, TVR-5a, TVB-08; DECISIONS Unit 55 17.3, 17.12,
// 17.13). Base-free, in double-double, with an independent cross-check by Chebyshev collocation (I3).
//
// Problem. Material 1 on x < 0 and material 2 on x > 0, each with its permittivity eps, electron affinity chi, gap
// E_g, densities of states N_c and N_v and a uniform net doping C (donors positive). The alignment is the electron
// affinity rule (Anderson; read as Schubert's eqs (17.1), (17.2): Delta E_c = chi_1 - chi_2, valid when interface
// dipoles are neglected, PB): the vacuum level E_vac = E_c + chi is continuous, and phi = -E_vac / q is the
// potential. With the Fermi level at 0, a neutral bulk has E_F - E_c = kT ln(n_b / N_c), so
//     phi_b = V_T ln(n_b / N_c) - chi,        V_bi = phi_b2 - phi_b1.
// Each side has n_i = sqrt(N_c N_v) e^(-E_g / (2 V_T)) (E_g in eV, V_T in V) and the exact neutral densities of
// first_integral.hpp (n_b - p_b = C, n_b p_b = n_i^2).
//
// First integral (each side, t = (phi - phi_b) / V_T): (1/2) phi'^2 = (q V_T / eps) Phi_b, Phi_b = n_b G(t) +
// p_b G(-t), G(t) = e^t - 1 - t, so the displacement is D^2 = 2 q V_T eps Phi_b. At the interface phi and D are
// continuous (Gauss's law with no interface charge, derived here, PD), so the interface level t_1 = (phi_0 -
// phi_b1) / V_T solves
//     eps_1 Phi_1(t_1) = eps_2 Phi_2(t_1 - V_bi / V_T)
// between 0 and V_bi / V_T, where the left side rises from 0 and the right side falls to 0: one root, by Brent's
// procedure in double-double. The profile on each side is first_integral.hpp's bulk_profile from the interface level.
//
// The infinite-domain solution is the reference. A device of finite length with Dirichlet contacts at the bulk
// levels differs from it by about the tail V_T t at the contacts; hetero_potential reports t per point, and the
// collocation, which solves the finite-domain problem, measures the difference directly.
//
// Cross-check (I3): hetero_collocation solves y'' = (q / (eps_k V_T)) (n - p - C_k) for y = (phi - phi_b1) / V_T on
// [x_left, 0] and [0, x_right] with y(x_left) = 0, y(x_right) = V_bi / V_T, y continuous and eps y' continuous at
// every piece boundary (Boyd's patching with the flux, chebyshev.hpp), by damped Newton with a dense LU in double.
// It shares no code with the first integral.
//
// n_i is rounded to a double before the neutral densities are formed (first_integral.hpp takes it so); every
// quantity here is exact for that n_i, and N_v enters only through it.
#pragma once

#include <expected>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct HeteroSide {
    double eps_F_per_cm = 0.0;
    double affinity_eV = 0.0;
    double band_gap_eV = 0.0;
    double Nc_cm3 = 0.0;
    double Nv_cm3 = 0.0;
    double net_doping_cm3 = 0.0;  // donors positive, acceptors negative; not zero
};

struct HeterojunctionInput {
    double q_C = 0.0;
    double thermal_voltage_V = 0.0;
    HeteroSide left, right;  // x < 0 and x > 0
};

struct Heterojunction {
    HeterojunctionInput in;
    double ni_left_cm3 = 0.0, ni_right_cm3 = 0.0;  // as rounded (see the header comment)
    DD built_in_V;             // phi_b2 - phi_b1
    DD interface_V;            // phi(0) - phi_b1
    DD displacement_C_per_cm2; // |D| at the interface
    DD field_left_V_per_cm, field_right_V_per_cm;  // |phi'| at x = 0-, 0+ (D / eps)
    DD bend_left_V, bend_right_V;                  // phi(0) - phi_b1 and phi_b2 - phi(0)
    DD ec_left_eV, ec_right_eV;  // E_c - E_F at x = 0-, 0+; their difference is chi_1 - chi_2
    DD debye_left_cm, debye_right_cm;
    long long brent_evaluations = 0;
};

// Errors: an input not finite, or not positive where it must be (the dopings: not zero); Brent not converging.
[[nodiscard]] std::expected<Heterojunction, std::string> heterojunction(const HeterojunctionInput& in);

struct HeteroPoint {
    DD phi_V;  // phi(x) - phi_b1
    DD tail;   // |phi - phi_bulk| / V_T of the point's own side
};

// phi at each x (any order, any sign; x = 0 is the interface, taken with the right side). Errors: a Newton or
// quadrature failure.
[[nodiscard]] std::expected<std::vector<HeteroPoint>, std::string> hetero_potential(const Heterojunction& j,
                                                                                    std::span<const double> x_cm);

struct HeteroCollocation {
    std::vector<double> phi_V;    // phi - phi_b1 at the query points
    double interface_V = 0.0;     // at x = 0
    double displacement_C_per_cm2 = 0.0;  // eps_1 |phi'(0-)|
    double flux_mismatch = 0.0;   // |eps_1 phi'(0-) - eps_2 phi'(0+)| / (eps_1 |phi'(0-)|), from the solution
    double tail_estimate_V = 0.0; // max over pieces of |a_{N-1}| + |a_N|, times V_T
    int newton_iterations = 0;
    double last_update = 0.0;     // max |dy| of the last Newton step
};

// The finite-domain problem on [x_left, x_right] (x_left < 0 < x_right), `subdomains` pieces per side of Chebyshev
// degree `degree`. Errors: Newton not converging in 80 iterations, a singular matrix, bad arguments.
[[nodiscard]] std::expected<HeteroCollocation, std::string> hetero_collocation(const Heterojunction& j,
                                                                               double x_left_cm, double x_right_cm,
                                                                               int subdomains, int degree,
                                                                               std::span<const double> query_cm);

}  // namespace NiTCAD::verify::reference
