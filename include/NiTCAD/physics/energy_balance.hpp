// The carrier-temperature mobility and effective field of the energy-balance model (Unit 52;
// DECISIONS.md ED4-ED6, the owner's ruling on beta = 1).
//
// In a uniform field E a carrier's energy balance is q v(E) E = (3/2) k (T - T_L) / tau_w, v the
// Canali drift velocity mu_C(E) E. Read backwards, a carrier temperature T defines the effective
// field E_eff(T) that heats it to T, and the energy-balance model takes
//     mu(T) = mu_C(E_eff(T)),   impact ionization alpha(E_eff(T)),
// so a uniform field reproduces Canali's v(E) and the field-driven ionization exactly.
// With x = mu0 E / v_sat and m(x) = (1 + x^beta)^(-1/beta) (so mu_C = mu0 m), the balance is
//     x^2 m(x) = a,   a = mu0 P / v_sat^2,   P = (3/2) (k / q) (T - T_L) / tau_w   [V/s].
// - beta = 2: y = x^2 solves y^2 = a^2 (1 + y), y = a (a + sqrt(a^2 + 4)) / 2, analytic through
//   a = 0 (y ~ a), so the mobility m = (1 + y)^(-1/2) continues smoothly below T_L, where cooled
//   carriers move faster (as in Hansch's model).
// - beta = 1: x = (a + sqrt(a^2 + 4 a)) / 2; any other beta: Newton in ln x. Near T_L these grow
//   like sqrt(a), and for beta < 2 the mobility's slope in T is infinite there (owner ruling,
//   2026-10-29): below a0 = 1e-8 (about 5 uK of heating in silicon) x is the cubic
//   c2 a^2 + c3 a^3 matching the exact x and slope at a0, with x = 0 for cooled carriers (a <= 0),
//   so x is C1 everywhere and mu = mu0 below T_L. The effective field takes this x for every beta
//   (it only feeds impact ionization, which is zero to all orders at such fields).
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

// Silicon's energy relaxation times (owner ruling EQ1, 2026-10-29): 0.35 ps for electrons and
// 0.4 ps for holes, the typical constants of A. Triebl, dissertation, TU Wien (iue.tuwien.ac.at/
// phd/triebl). Applied by the model to regions with silicon's mobility that give none.
inline constexpr EnergyRelaxationParameters silicon_energy_relaxation{.tau_n_s = 0.35e-12,
                                                                      .tau_p_s = 0.4e-12};

// The reduced balance's regularization threshold a0 (owner ruling on beta = 1).
inline constexpr double energy_balance_regularization = 1e-8;

// x(a) and dx/da of the balance x^2 m(x) = a.
struct ReducedField {
    double x;
    double d_a;
};

// The regularized x(a) the model uses (above). Preconditions (NITCAD_EXPECTS): a finite, beta
// finite and >= 1.
[[nodiscard]] ReducedField reduced_field(double a, double beta);

// The exact root for a > 0 by Newton in ln x for any beta (the method used for beta other than 1
// and 2), without regularization; exposed for the closed forms' tests (EC2). Preconditions: a
// finite and positive, beta finite and >= 1.
[[nodiscard]] ReducedField reduced_field_newton(double a, double beta);

struct CarrierTemperatureMobility {
    double mobility;      // mu(T) [cm^2/(V s)]
    double d_mu0;         // d mu / d mu0 at fixed T
    double d_rise;        // d mu / d (T - T_L) [cm^2/(V s K)]
    double field;         // E_eff(T) [V/cm]
    double field_d_mu0;   // d E_eff / d mu0 [V s / cm^3]
    double field_d_rise;  // d E_eff / d (T - T_L) [V/(cm K)]
};

// mu(T) and E_eff(T) for a carrier of low-field mobility mu0 (doping and surface terms included)
// heated by rise_K = T - T_L, with energy relaxation time tau_w_s and the carrier's Canali
// parameters. Preconditions (NITCAD_EXPECTS): mu0 finite and positive; rise_K finite; tau_w_s
// finite and positive; the Canali parameters as Semiconductor::create validates them.
[[nodiscard]] CarrierTemperatureMobility carrier_temperature_mobility(double mu0, double rise_K,
                                                                      double tau_w_s,
                                                                      const CanaliParameters& c);

}  // namespace NiTCAD::physics
