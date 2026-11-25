// Philips unified (Klaassen) low-field mobility (Unit 66; DECISIONS.md, Unit 66): D. B. M.
// Klaassen, Solid-State Electron. 35, 953 and 961 (1992), as reproduced by the Sandia Charon User
// Manual v2.2 section 7.4, DEVSIM and COMSOL (the papers themselves were not read). For carrier i,
// the other carrier j:
//
//     1 / mu_i = 1 / mu_i,L + 1 / mu_i,DAP                      (Matthiessen)
//     mu_i,L   = mu_max (T / 300)^-theta
//     mu_i,DAP = mu_N (N_sc / N_sc,eff) (N_ref / N_sc)^alpha + mu_c (n + p) / N_sc,eff
//     mu_N     = mu_max^2 / (mu_max - mu_min) (T / 300)^(3 alpha - 1.5)
//     mu_c     = mu_max mu_min / (mu_max - mu_min) (300 / T)^0.5
//     N_D*     = N_D (1 + 1 / (c_D + (N_ref,D / N_D)^2)), N_A* likewise   (clustering)
//     N_e,sc   = N_D* + N_A* + p,  N_e,sc,eff = N_D* + G(P_e) N_A* + p / F(P_e)   (holes mirrored)
//     P_i      = [f_CW / P_CW + f_BH / P_BH]^-1,
//     P_CW     = 3.97e13 N_sc^(-2/3) (T / 300)^2,  P_BH = 1.36e20 / (n + p) (m_i / m0) (T / 300)^2
//     G(P)     = 1 - s1 / [s2 + ((m0 / m_i)(T / 300))^s4 P]^s3 + s5 / [((m_i / m0)(300 / T))^s7 P]^s6,
//                held at its minimum below the minimiser
//     F(P)     = (r1 P^r6 + r2 + r3 m_i / m_j) / (P^r6 + r4 + r5 m_i / m_j)
//
// The mobility depends on n and p, so it returns their partials (and the temperature's) for
// Newton (ARCHITECTURE.md section 5, item 2). It is not differentiable where N_D = N_A = 0 and
// the screening carrier density is 0 (1 / mu_DAP goes as p^alpha there); the partials returned
// at such a point are the finite one-sided limits where they exist and otherwise 0.
#pragma once

#include <array>

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

struct KlaassenCarrierParameters {
    double mu_max;      // [cm^2/(V s)]
    double mu_min;      // [cm^2/(V s)]
    double N_ref;       // [cm^-3]
    double alpha;
    double theta;       // lattice exponent: mu_L = mu_max (T / 300)^-theta
    double mass_ratio;  // m_i / m0 of the screening functions
};

struct KlaassenParameters {
    KlaassenCarrierParameters electron;
    KlaassenCarrierParameters hole;
    double c_D, c_A;            // clustering coefficients
    double N_ref_D, N_ref_A;    // clustering reference densities [cm^-3]
    double f_CW, f_BH;          // screening weights
    std::array<double, 7> s;    // G(P)
    std::array<double, 6> r;    // F(P)
};

// Silicon: phosphorus electrons and boron holes (Charon Table 7-8; COMSOL's defaults), the
// common parameters of Charon eqs (71)-(74) and Table 7-9 and COMSOL (r1, r5 as those two; the
// ATLAS manual's 0.7463 and -0.8552 are not taken; DECISIONS.md, Unit 66).
inline constexpr KlaassenParameters klaassen_silicon{
    .electron = {.mu_max = 1414.0, .mu_min = 68.5, .N_ref = 9.20e16, .alpha = 0.711,
                 .theta = 2.285, .mass_ratio = 1.0},
    .hole = {.mu_max = 470.5, .mu_min = 44.9, .N_ref = 2.23e17, .alpha = 0.719, .theta = 2.247,
             .mass_ratio = 1.258},
    .c_D = 0.21,
    .c_A = 0.50,
    .N_ref_D = 4.0e20,
    .N_ref_A = 7.2e20,
    .f_CW = 2.459,
    .f_BH = 3.828,
    .s = {0.89233, 0.41372, 0.19778, 0.28227, 0.005978, 1.80618, 0.72169},
    .r = {0.7643, 2.2999, 6.5502, 2.3670, -0.01552, 0.6478}};

struct KlaassenMobility {
    double mobility;  // [cm^2/(V s)]
    double d_dn;      // [cm^2/(V s) cm^3]
    double d_dp;
    double d_dT;      // [cm^2/(V s K)]
};

// The temperature-dependent parts of one carrier's mobility at one temperature, computed once
// (they include a root solve for G's minimiser) and reused for every node at that temperature.
struct KlaassenAtTemperature {
    const KlaassenParameters* parameters;
    Carrier carrier;
    double temperature_K;
    double P_min;  // G's minimiser [1]
};

// Preconditions (NITCAD_EXPECTS): temperature_K finite and positive.
[[nodiscard]] KlaassenAtTemperature klaassen_at(const KlaassenParameters& K, Carrier carrier,
                                                double temperature_K);

// Preconditions (NITCAD_EXPECTS): N_D, N_A, n, p finite and >= 0.
[[nodiscard]] KlaassenMobility klaassen_mobility(const KlaassenAtTemperature& at, double N_D,
                                                 double N_A, double n, double p);

// Convenience: klaassen_mobility(klaassen_at(K, carrier, T), N_D, N_A, n, p).
[[nodiscard]] KlaassenMobility klaassen_mobility(const KlaassenParameters& K, Carrier carrier,
                                                 double N_D, double N_A, double n, double p,
                                                 double temperature_K);

// The neutral equilibrium densities of a node's doping: with N = N_D - N_A,
// n0 = N / 2 + sqrt((N / 2)^2 + n_i^2), p0 = n_i^2 / n0 (mirrored for N < 0, without
// cancellation). Preconditions: N_D, N_A finite and >= 0; n_i finite and positive.
struct EquilibriumDensities {
    double n, p;          // [cm^-3]
    double dn_dni, dp_dni;
};
[[nodiscard]] EquilibriumDensities neutral_equilibrium(double N_D, double N_A, double n_i);

// The `klaassen_equilibrium` mobility (KD2): klaassen_mobility at the node's neutral equilibrium
// densities with the given n_i; d_dT the total derivative, with d n_i / dT = n_i * dlnni_dT.
// d_dn and d_dp are those of klaassen_mobility at (n0, p0).
[[nodiscard]] KlaassenMobility klaassen_equilibrium_mobility(const KlaassenAtTemperature& at,
                                                             double N_D, double N_A, double n_i,
                                                             double dlnni_dT);

}  // namespace NiTCAD::physics
