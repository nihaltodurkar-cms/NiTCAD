// Gate tunnelling through insulators (Unit 51; DECISIONS.md GD1-GD5). Pure functions of one barrier:
// the WKB transmission of a piecewise-linear barrier and the Tsu-Esaki current through it.
//
// Energies in eV on one electron-energy scale, lengths in cm, masses in units of m0. A hole process
// is evaluated on the mirrored scale (every energy negated) by the caller, so one routine serves the
// electron (ECB, EVB) and hole (HVB) processes.
//
// WKB (Gehring, "Simulation of Tunneling in Semiconductor Devices", TU Wien 2003, section 3.5.1):
//     TC(E) = exp(-2 S(E)),   S(E) = sum over segments of int sqrt(2 m q (W(x) - E)+) / hbar dx
// with W linear along each segment from W0 to W1 over its length L, so with a = W0 - E, b = W1 - E
//     int_0^L sqrt((a + (b - a) x / L)+) dx = L g(a, b),   g = (2/3) (a+^(3/2) - b+^(3/2)) / (a - b)
// exact with or without a turning point inside the segment (the clipped part contributes nothing),
// evaluated as (2/3) (a + sqrt(a b) + b) / (sqrt a + sqrt b) when both ends are inside (no
// cancellation as b -> a, where g -> sqrt a).
//
// Tsu-Esaki (Gehring section 3.2, 3.3.1), with A* = 4 pi m q k^2 / h^3 the side's Richardson constant:
//     J = (A* / k^2) int TC(E) N(E) dE,   N(E) = kT [ln(1 + e^((F1 - E)/kT)) - ln(1 + e^((F2 - E)/kT))]
// positive from side 1 (Fermi level F1) to side 2 (F2), zero when F1 = F2. With TC = 1 above E_c and
// a Boltzmann supply it is A* T^2 e^(-(E_c - F1)/kT), thermionic emission.
#pragma once

#include <span>

namespace NiTCAD::physics {

// sqrt(2 m0 q) / hbar [1 / (cm eV^(1/2))]: S of a 1 cm, 1 eV rectangular barrier at m = m0.
[[nodiscard]] double wkb_action_scale() noexcept;

// One linear piece of a barrier: its height at the two ends [eV] (the insulator's band edge, on the
// process's scale), its length [cm] and the tunnelling mass [m0] of its insulator.
struct BarrierSegment {
    double W0_eV;
    double W1_eV;
    double length_cm;
    double mass_m0;
};

// L g(a, b) [cm eV^(1/2)] and its partials in a and b. Zero (with zero partials) when a, b <= 0.
struct SegmentAction {
    double value;
    double d_a;
    double d_b;
};
[[nodiscard]] SegmentAction segment_action(double a_eV, double b_eV, double length_cm) noexcept;

// S(E) = sum sqrt(m) L g / scale; dS/dE and, when d_W is not empty (size 2 n, (W0, W1) per
// segment), dS/dW added to d_W (scaled by `weight`).
[[nodiscard]] double wkb_exponent(std::span<const BarrierSegment> barrier, double E_eV,
                                  double* d_E = nullptr, std::span<double> d_W = {},
                                  double weight = 1.0);

// TC(E) = exp(-2 S(E)).
[[nodiscard]] double wkb_transmission(std::span<const BarrierSegment> barrier, double E_eV);

// N(E) [eV] and its partials in F1, F2 and E. Exactly 0 when F1 == F2.
struct Supply {
    double value;
    double d_F1;
    double d_F2;
    double d_E;
};
[[nodiscard]] Supply tsu_esaki_supply(double E_eV, double F1_eV, double F2_eV, double kT_eV) noexcept;

// The current density through `barrier` from side 1 to side 2 over energies [E_low, E_high], cut at
// the barrier's highest point (no over-barrier current, GD3) and at max(F1, F2, E_low) + 40 kT (the
// supply is below e^-40 of its peak there). A composite 4-point Gauss-Legendre rule with panels of
// at most kT / 2; the partials are the exact derivatives of that discrete sum (the cut's too).
// J [A/cm^2] with A* in A/(cm^2 K^2). Zero with zero partials when the range is empty.
// d_W (size 2 n) receives dJ/dW per segment end when not empty.
struct TunnelCurrent {
    double J;
    double d_F1;
    double d_F2;
    double d_E_low;
    double d_E_high;
};
struct TunnelRange {
    double E_low_eV;
    double E_high_eV;  // +infinity: the barrier top
};
[[nodiscard]] TunnelCurrent tsu_esaki_current(std::span<const BarrierSegment> barrier,
                                              TunnelRange range, double F1_eV, double F2_eV,
                                              double temperature_K, double richardson_A_cm2K2,
                                              std::span<double> d_W = {});

}  // namespace NiTCAD::physics
