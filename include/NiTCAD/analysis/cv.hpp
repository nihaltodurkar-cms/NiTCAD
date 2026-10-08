// Capacitance-voltage and interface extraction (ARCHITECTURE.md section 11, Unit 24): the
// quasi-static capacitance of a gate charge, the accumulation capacitance, the doping profile from
// d(1/C^2)/dV, the flat-band voltage, and interface traps by the conductance method.
//
// The doping profile and the conductance method take capacitances per unit area (F/cm^2): a 1D
// device's results already are; scale a 2D or 3D device's (curve.hpp scaled) by 1 / (its gate
// width in cm) or 1 / (its gate area in cm^2). analysis does not read the device, so the material
// numbers (permittivity, oxide capacitance, doping) are passed in.
#pragma once

#include <cstddef>
#include <expected>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/small_signal.hpp"

namespace NiTCAD::analysis {

// dQ/dV of a gate-charge curve (curve.hpp gate_charge_curve) at every point, as a curve over the
// same biases (no resolution). Errors: those of derivative.
[[nodiscard]] std::expected<Curve, base::Error> quasi_static_capacitance(const Curve& charge);

// The largest capacitance of the curve: in strong accumulation the oxide capacitance, approached
// from below (the accumulation layer's own capacitance in series). Errors: it is not positive.
[[nodiscard]] std::expected<Extraction, base::Error> accumulation_capacitance(
    const Curve& capacitance);

struct DopingProfile {
    std::vector<double> depth_cm;    // the depletion width eps / C_d at each point
    std::vector<double> doping_cm3;  // 2 / (q eps |d(1/C_d^2)/dV|)
};

// The depletion approximation's profile from a C-V curve in depletion (no inversion charge):
// C_d = C with no series capacitance, else 1 / (1/C - 1/series_capacitance) (the oxide of a MOS
// capacitor); w = eps / C_d and N(w) = 2 / (q eps |d(1/C_d^2)/dV|), the derivative by curve.hpp
// derivative (first order at the ends). Pass the depletion part of the curve (curve.hpp slice).
// Errors: permittivity not positive; series_capacitance negative; a capacitance not positive or
// not below series_capacitance; the derivative's errors or a zero derivative.
[[nodiscard]] std::expected<DopingProfile, base::Error> doping_profile(
    const Curve& capacitance, double permittivity_F_per_cm, double series_capacitance = 0.0);

// 1 / (1/C_ox + L_D / eps) with the extrinsic Debye length L_D = sqrt(eps V_T / (q N)): the
// capacitance at flat band. Preconditions: every argument positive.
[[nodiscard]] double flat_band_capacitance(double oxide_capacitance, double permittivity_F_per_cm,
                                           double doping_cm3, double temperature_K);

// The gate bias where the capacitance first reaches flat_band (crossing, linear). Errors: those of
// crossing.
[[nodiscard]] std::expected<Extraction, base::Error> flat_band_voltage(const Curve& capacitance,
                                                                       double flat_band);

struct ConductancePeak {
    Extraction conductance_over_omega;  // (G_p / omega) at the peak, in F cm^(D-3)
    Extraction frequency;               // where it peaks, Hz; window: frequency indices
};

// The conductance method (Nicollian and Goetzberger): at one operating point, Y = Y_contact,contact
// with the oxide removed, Y_s = 1 / (1/Y - 1 / (i omega C_ox)), and G_p / omega = Re Y_s / omega
// over the positive frequencies. Its largest sample and the samples either side, refined by the
// parabola in ln f through their ln (G_p / omega) (through G_p / omega itself if a neighbour is not
// positive). A single trap level gives C_it / 2 at omega tau = 1; a continuum
// gives about q D_it / 2.5 (interface_trap_density). Errors: point or contact out of range;
// oxide_capacitance not positive; fewer than 3 positive frequencies, not increasing; the largest
// at the first or last of them.
[[nodiscard]] std::expected<ConductancePeak, base::Error> conductance_peak(
    const results::SmallSignal& small_signal, std::size_t point, std::size_t contact,
    double oxide_capacitance);

// D_it = 2.5 (G_p / omega)_max / q in cm^-2 eV^-1, (G_p / omega)_max per unit area (F/cm^2).
[[nodiscard]] double interface_trap_density(double peak_conductance_over_omega);

}  // namespace NiTCAD::analysis
