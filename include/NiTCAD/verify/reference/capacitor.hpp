// Capacitor references for TVB-05 and TVB-06 (DECISIONS Unit 55 16.3, 16.4; criteria 16.8).
//
// Exact electrostatics of linear dielectrics between equipotential electrodes:
// - a planar stack of layers t_i, eps_i in series:    C = eps0 / sum_i (t_i / eps_r,i)    [F/cm^2];
//   in layer i the field is E_i = Q / (eps0 eps_r,i), so E_i / E_j = eps_r,j / eps_r,i (D continuous);
// - a coaxial capacitor, radii a < b, per unit length: C = 2 pi eps0 eps_r / ln(b / a)      [F/cm].
// Both follow from Gauss's law with a uniform D (planar) or D = Q / (2 pi r) (coaxial); provenance PD
// (derived from G1 alone, no external source needed), with a second method in each function's L0 test:
// the series sum against a Gauss-Kronrod integral of dx / eps(x), the coaxial form against
// the integral of dr / (2 pi eps r).
// Computed in double-double (U_ref at rounding of the result); base only (TVD-19).
#pragma once

#include <span>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

// eps0 [F/cm] as the reference uses it: CODATA 2018's 8.8541878128e-12 F/m, divided by 100.
[[nodiscard]] DD vacuum_permittivity_F_cm();

// The series capacitance of planar layers [F/cm^2]. Preconditions (NITCAD_EXPECTS): the spans have the
// same, non-zero size; every thickness and permittivity finite and positive.
[[nodiscard]] DD series_capacitance(std::span<const double> thickness_cm, std::span<const double> eps_r);

// The coaxial capacitance per unit length [F/cm]. Preconditions: 0 < a < b, eps_r > 0, all finite.
[[nodiscard]] DD coaxial_capacitance(double inner_cm, double outer_cm, double eps_r);

}  // namespace NiTCAD::verify::reference
