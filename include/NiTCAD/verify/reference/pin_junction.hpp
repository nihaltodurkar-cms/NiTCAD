// The abrupt p-i-n junction in the quasi-equilibrium (QE) model: the first integral across three regions, the field at
// any level, and the level at any position (Unit 55, TVR-5a, TVB-24; DECISIONS Unit 55 17.3, 17.13 (e)). Base-free,
// in double-double.
//
// The device. One material (one permittivity eps): the p side x < 0 (acceptors N_A), the i region 0 <= x <= d (net
// doping C_i, donors positive, -N_A < C_i < N_D; 0 for an intrinsic region), the n side x > d (donors N_D); the anode
// bias V <= 0 on the p side. The QE model as qe_junction.hpp: electrons at the n side's quasi-Fermi level and holes at
// the p side's, m = n_i e^(V / (2 V_T)). The level s = u - u_p is the reduced potential measured from the p bulk, and
// B = u_n - u_p = ln n_n + ln p_p - 2 ln m the barrier (V_bi - V = V_T B). Each density is referred to its own bulk,
// n = n_n e^(s - B) and p = p_p e^-s, the minority bulk densities carried as logarithms as in multiplication.hpp, so
// any reverse bias stays in range.
//
// The first integral (17.13 (e); derived here, PD). With E = V_T ds/dx > 0 and scale = 2 q / (eps V_T),
// E = V_T sqrt(scale Phi):
//   p side, 0 < s <= s_0:      Phi_p(s) = p_p G(-s) + n_p G(s),
//   i region, s_0 <= s <= s_d: Phi_i(s) = Phi_p(s_0) + n_n e^(s_0 - B) expm1(s - s_0) + p_p e^-s_0 expm1(s_0 - s)
//                                         - C_i (s - s_0),
//   n side, s_d <= s < B:      Phi_n(r) = n_n G(-r) + p_n G(r), r = B - s,
// G(t) = e^t - 1 - t. E is continuous at both edges (one permittivity). The n edge closes in one step: d/ds (Phi_i(s)
// - Phi_n(B - s)) = (n - p - C_i) - (n - p - N_D) = N_D - C_i is a constant (Gauss's law across the i region), so
// s_d = s_0 + (Phi_n(B - s_0) - Phi_p(s_0)) / (N_D - C_i), then one Newton correction on the computed residual. s_0
// is the root of L(s_0) = V_T int_{s_0}^{s_d} ds / E_i(s) = d (Brent): L falls as s_0 rises, from infinity (where
// Phi_i first touches 0 inside the region, or s_d reaches B) to 0 at the two-region junction's level, where
// s_d = s_0. The quadrature is adaptive Kronrod-15 at relative 1e-26. Phi_i > 0 on [s_0, s_d] is checked at its one
// stationary point (n - p = C_i) and at both ends.
//
// Positions. x(s) = -V_T int_s^{s_0} ds / E on the p side, V_T int_{s_0}^s ds / E in the i region and
// d + V_T int_{s_d}^s ds / E on the n side, the two sides in the logarithm of the distance to their bulk level
// (t = ln s, t = ln (B - s)), where the integrand V_T e^t / E stays bounded (it tends to the bulk's Debye length). The
// inverse, x -> s, by Brent in the same variable. dE/dx = (q / eps)(n - p - C) at a level, for the first-order terms
// of a one-sided slope (TVB-24).
//
// Second method (pin_collocation): s'' = (q / (eps V_T))(n - p - C) on [x_left, x_right] with s = 0 and s = B at the
// two ends, by Chebyshev collocation on pieces of each region (chebyshev.hpp; Boyd's patching with one permittivity:
// s and s' continuous), damped Newton from the depletion approximation and a dense LU in double. It shares no code
// with the first integral, and it solves the finite domain: its ends differ from the bulk levels by the tails, which
// are exponentially small when the ends lie many Debye lengths outside the depletion layers.
#pragma once

#include <expected>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct PinInput {
    double q_C = 0.0;
    double eps_F_per_cm = 0.0;
    double thermal_voltage_V = 0.0;
    double ni_cm3 = 0.0;
    double acceptors_cm3 = 0.0;       // x < 0
    double intrinsic_net_cm3 = 0.0;   // 0 <= x <= d: donors minus acceptors, in (-N_A, N_D)
    double donors_cm3 = 0.0;          // x > d
    double i_length_cm = 0.0;         // d > 0
};

struct PinJunction {
    PinInput in;
    DD VT, scale;            // E = V_T sqrt(scale Phi)
    DD ln_m, p_p, n_n, ln_np, ln_pn;
    DD barrier;              // B (reduced); V_bi - V = V_T B
    DD s0, sd;               // the levels at x = 0 and x = d
    DD phi0;                 // Phi_p(s_0) = Phi_i(s_0)
    DD field_p_edge_V_per_cm, field_n_edge_V_per_cm;  // E at x = 0 and x = d
    DD field_min_i_V_per_cm; // the smallest E in the i region
    DD depletion_p_cm, depletion_n_cm;  // eps E_edge / (q N) of each side, reported
    double length_residual_cm = 0.0;    // |L(s_0) - d| at the root
    double matching_residual = 0.0;     // |Phi_i(s_d) - Phi_n(B - s_d)| / Phi_n(B - s_d)
    long long evaluations = 0;          // length evaluations of the root search
};

// Errors: an input not finite and positive where it must be (C_i finite, in (-N_A, N_D)); a bias not finite or forward
// (V > 0); no barrier; the root not bracketed or Brent not converging; a quadrature not converging.
[[nodiscard]] std::expected<PinJunction, std::string> pin_junction(const PinInput& in, double bias_V);

// E [V/cm] at level s in (0, B).
[[nodiscard]] DD pin_field(const PinJunction& j, const DD& s);

// E [V/cm] on the n side at r = B - s levels from the n bulk, 0 < r <= B - s_d, without forming s (B - (B - r)
// would lose r's digits once r << 1).
[[nodiscard]] DD pin_field_n_bulk(const PinJunction& j, const DD& r);

// dE/dx = (q / eps)(n - p - C) [V/cm^2] at level s in (0, B), C the region's net doping (s_0 and s_d belong to the
// i region).
[[nodiscard]] DD pin_field_slope(const PinJunction& j, const DD& s);

// x(s) [cm] for s in (0, B). Errors: s outside (0, B); the quadrature not converging.
[[nodiscard]] std::expected<DD, std::string> pin_position(const PinJunction& j, const DD& s);

// s at position x [cm]. Errors: x deeper in a bulk than the level e^-700 reaches; Brent or a quadrature failing.
[[nodiscard]] std::expected<DD, std::string> pin_level(const PinJunction& j, double x_cm);

struct PinCollocation {
    double s0 = 0.0, sd = 0.0;          // s at x = 0 and x = d
    double field_p_edge_V_per_cm = 0.0;  // E at x = 0 (from the p side's piece)
    double field_n_edge_V_per_cm = 0.0;  // E at x = d (from the n side's piece)
    double field_mid_V_per_cm = 0.0;     // E at x = d / 2
    std::vector<double> levels;          // s at the query points
    double tail_estimate = 0.0;          // max over pieces of |a_{N-1}| + |a_N| (reduced)
    int newton_iterations = 0;
    double last_update = 0.0;
};

// The finite-domain problem on [x_left, x_right] (x_left < 0, x_right > d): each region cut into pieces no wider than
// `piece_debye` Debye lengths of the heavier side, Chebyshev degree `degree` on each. Errors: bad arguments, a
// singular matrix, Newton not converging in 200 iterations.
[[nodiscard]] std::expected<PinCollocation, std::string> pin_collocation(const PinJunction& j, double x_left_cm,
                                                                         double x_right_cm, int degree,
                                                                         double piece_debye,
                                                                         std::span<const double> query_cm);

}  // namespace NiTCAD::verify::reference
