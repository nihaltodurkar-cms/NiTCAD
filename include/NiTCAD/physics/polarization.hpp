// Wurtzite nitride polarization (Unit 68 slice 2, MAD9; DECISIONS.md amendment 6), Ga-face (metal
// face), the c-axis pointing from substrate to surface. The c-components [C/m^2]:
//     P = P_sp + P_pz,  P_pz = 2 (a - a0) / a0 (e31 - e33 c13 / c33),
// a the layer's in-plane lattice constant, a0 its relaxed one (O. Ambacher et al., J. Phys.:
// Condens. Matter 14, 3399 (2002), eqs (3), (4)). The bound charge is P.n on a surface and
// P_below - P_above on an interface between stacked layers.
// Sources: P_sp, e31, e33, c13, c33 and the P_sp bowing, the owner's values of 2026-11-27; a0, M2'
// (Levinshtein et al. 2001: GaN 3.189, AlN 3.112, InN 3.533 A). Alloys: P_sp linear with bowing,
// P_sp(x) = (1 - x) P_GaN + x P_B - b x (1 - x), b -0.021 (AlGaN) and -0.037 (InGaN) C/m^2; a0, e31,
// e33, c13 and c33 linear in x (Vegard), Ambacher's practice.
#pragma once

#include <optional>
#include <string_view>

namespace NiTCAD::physics {

struct NitrideLattice {
    double a0_A;   // relaxed in-plane lattice constant [Angstrom]
    double P_sp;   // spontaneous polarization [C/m^2]
    double e31, e33;  // piezoelectric constants [C/m^2]
    double c13, c33;  // elastic constants [GPa]
};

// The lattice of `gan`, `aln`, `inn` (x ignored), `algan x` or `ingan x` (0 <= x <= 1, the
// caller's checked fraction); nullopt for any other material name.
[[nodiscard]] std::optional<NitrideLattice> nitride_lattice(std::string_view material, double x);

// P_sp + P_pz [C/m^2] of a layer whose in-plane lattice constant is a_A [Angstrom].
[[nodiscard]] double nitride_polarization(const NitrideLattice& l, double a_A);

// A polarization [C/m^2] as a sheet number density [cm^-2] (q = 1.602176634e-19 C).
[[nodiscard]] double polarization_sheet_density_cm2(double P_C_m2);

}  // namespace NiTCAD::physics
