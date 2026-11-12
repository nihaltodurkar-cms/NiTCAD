// Lattice heat and thermoelectric references for TVB-27 and TVB-28 (Unit 55 TVR-4; DECISIONS 16.3, 16.18 (k),
// 16.43). Derived here (PD), base-free: every constant and material parameter is an input (R-7).
//
// Conduction: kappa(T) = kappa300 (T / 300)^-a (Unit 23 T11's form). The Kirchhoff transform
//     u(T) = int_300^T kappa dT = kappa300 300 / (1 - a) ((T / 300)^(1 - a) - 1)        (a != 1)
// and its inverse T(u) = 300 (1 + (1 - a) u / (300 kappa300))^(1 / (1 - a)).
// (i) The slab with isothermal ends and no source: u linear in x, T(x) = T(u(T_0) + (x / L)(u(T_1) - u(T_0))).
// (ii) The self-heated n-resistor (Unit 23 T2, T4 with n = N throughout, holes neglected, r any): with Stratton's
//     flux J = q mu(T) N (E_c' + (1 + r)(k/q) T') and the heat H = J e_n', e_n = E_c + (5/2 + r) k T,
//         (kappa(T) T')' = -J^2 / (q mu(T) N) - (3/2)(k/q) J T'      (the 3/2 is (5/2 + r) - (1 + r), any r),
//     T = T_amb at both ends, and J = -q N V / int_0^L dx / mu(T) (the ohmic ends at one temperature, so E_c's
//     drop is -V and the thermodiffusion term integrates to 0). mu(T) = Caughey-Thomas with mu_max (T / 300)^gamma.
//     Solved in double-double: Dormand-Prince on (u, u', int dx / mu) at fixed steps, shooting on u'(0) by Brent
//     and J by the secant method on J + q N V / int dx / mu. Second method: the Kirchhoff-form finite differences on
//     2^14 and 2^15 cells (Newton, double), Richardson-extrapolated.
// (iii) The Seebeck rise of an n-bar at n = N with Nc = Nc300 (T / 300)^(3/2): at J = 0 the electron quasi-Fermi
//     level E_Fn = E_c + kT ln(N / Nc) satisfies E_Fn' = -k T' ((5/2 + r) - ln(N / Nc(T))), so
//         E_Fn(T_1) - E_Fn(T_0) = -k [A (T_1 - T_0) + (3/2)(T ln(T / 300) - T)|_{T_0}^{T_1}],
//     A = 5/2 + r - ln(N / Nc300) [eV]. Second method: Gauss-Kronrod on the integrand.
// (iv) The Peltier coefficient Pi_n = T P_n = -(k T / q)((5/2 + r) - ln(N / Nc(T))) [V].
#pragma once

#include <expected>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct Conductivity {
    double kappa300_W_cmK = 0.0;
    double exponent = 0.0;  // a
};

[[nodiscard]] DD kirchhoff(const Conductivity& c, const DD& T);
[[nodiscard]] DD kirchhoff_inverse(const Conductivity& c, const DD& u);

// (i) T at the fractions s = x / L of the slab.
[[nodiscard]] std::expected<std::vector<DD>, std::string> slab_temperature(const Conductivity& c, double T0,
                                                                          double T1, const std::vector<double>& s);

struct ResistorInput {
    double q_C = 0.0, k_eV_per_K = 0.0;
    double length_cm = 0.0, density_cm3 = 0.0, voltage_V = 0.0, ambient_K = 300.0;
    Conductivity kappa;
    double mu_min = 0.0, mu_max = 0.0, N_ref = 0.0, alpha = 0.0, gamma = 0.0;  // Caughey-Thomas
    long long steps = 512;
    bool run_second_method = true;  // the finite differences (about 1 s in Release); off in L0 (CV-Q19's Debug budget)
};

struct ResistorProfile {
    DD current_A_cm2;              // J (signed: negative for V > 0)
    std::vector<DD> temperature;   // at the requested x
    double steps_difference = 0.0;  // max |T(steps) - T(steps / 2)| [K] at the requested x
    double second_method = 0.0;     // max |T - T_FD| [K], the finite differences, Richardson-extrapolated
    int iterations = 0;             // J's secant iterations
};

// (ii) T at the requested positions x [cm] (in [0, L]).
[[nodiscard]] std::expected<ResistorProfile, std::string> self_heated_resistor(const ResistorInput& in,
                                                                              const std::vector<double>& x);

struct SeebeckInput {
    double k_eV_per_K = 0.0;
    double density_cm3 = 0.0, Nc300_cm3 = 0.0, thermopower_exponent = -0.5;
    double T0 = 0.0, T1 = 0.0;
};

struct SeebeckResult {
    DD rise_eV;                  // E_Fn(T1) - E_Fn(T0)
    double second_method = 0.0;  // |GK - closed form| [eV]
};

// (iii)
[[nodiscard]] std::expected<SeebeckResult, std::string> seebeck_rise(const SeebeckInput& in);

// (iv) Pi_n [V] at T.
[[nodiscard]] DD peltier_coefficient(double k_eV_per_K, double density_cm3, double Nc300_cm3,
                                     double thermopower_exponent, double T);

}  // namespace NiTCAD::verify::reference
