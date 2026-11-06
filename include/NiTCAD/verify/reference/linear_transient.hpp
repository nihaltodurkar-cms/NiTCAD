// Linear transients of a uniform semiconductor bar (Unit 55, TVR-4, TVB-19; DECISIONS Unit 55 16.2 (m), 16.18 (h)),
// in double-double, base-free.
//
// (i) Dielectric relaxation. A uniform bar of majority density N and mobility mu has the conductivity sigma = q mu N.
// For a small excess majority charge rho, continuity and Gauss give d rho / dt = -rho / tau_d + D rho'' with
//     tau_d = eps / (q mu N),
// so the charge's integral over a bar whose ends carry no excess flux decays exactly as e^{-t / tau_d}; diffusion
// only spreads it. Second method: Dormand-Prince on dQ/dt = -Q / tau_d (tight tolerance), against the exponential.
//
// (ii) Haynes and Shockley (Phys. Rev. 81, 835 (1951), PA for the model; the moments derived here, PD). Minority
// carriers in a uniform field E at low injection obey d p / dt = D p'' - v p' - p / tau + g, v = mu E (the
// carrier's drift velocity, signed). A source g(x, t) = R exp(-(x - x0)^2 / (2 s^2)) on (t_on, t_off] is a
// superposition of Gaussians released at t'; one released at age a = t - t' has amount A0 e^{-a / tau},
// A0 = R s sqrt(2 pi) dt', centre x0 + v a and variance s^2 + 2 D a. So at t >= t_off, with ages a in
// [a1, a2] = [t - t_off, t - t_on] weighted by w(a) = e^{-a / tau}:
//     area   = R s sqrt(2 pi) int w da,
//     centre = x0 + v <a>,
//     width^2 = s^2 + 2 D <a> + v^2 (<a^2> - <a>^2),        <f> = int f w da / int w da,
// in closed form (the exponential's moments on [a1, a2]). Second method: the three integrals by Gauss-Kronrod in
// double-double.
#pragma once

#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

// tau_d = eps / (q mu N). Errors: an input not finite and positive.
[[nodiscard]] std::expected<DD, std::string> dielectric_relaxation_time(double eps_F_per_cm, double q_C,
                                                                       double mobility_cm2_Vs, double density_cm3);

// The second method for (i): Q(t) / Q(0) by Dormand-Prince on dQ/dt = -Q / tau_d, against e^{-t / tau_d}; returns
// the largest absolute difference over t in [0, t_end] at the integrator's steps.
[[nodiscard]] std::expected<double, std::string> relaxation_second_method(const DD& tau_d, double t_end_s);

struct PulseInput {
    double diffusivity_cm2_s = 0.0;  // D > 0
    double velocity_cm_s = 0.0;      // v = mu E, signed
    double lifetime_s = 0.0;         // tau > 0 (infinite allowed: no recombination)
    double centre_cm = 0.0;          // x0
    double sigma_cm = 0.0;           // s > 0
    double rate_cm3_s = 0.0;         // R > 0, the source's peak
    double t_on_s = 0.0, t_off_s = 0.0;  // t_on < t_off
};

struct PulseMoments {
    DD area_cm2;   // int p dx [cm^-2]
    DD centre_cm;
    DD width_cm;   // the standard deviation
    // The second method's differences (relative to area and width; centre relative to v (t - t_on)).
    double area_difference = 0.0, centre_difference = 0.0, width_difference = 0.0;
};

// At t >= t_off. Errors: bad inputs, t < t_off, a quadrature failure.
[[nodiscard]] std::expected<PulseMoments, std::string> haynes_shockley_moments(const PulseInput& in, double t_s);

}  // namespace NiTCAD::verify::reference
