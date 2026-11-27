// beta-Ga2O3 (Unit 68 slice 4b, MAD11; DECISIONS.md amendment 9), `ga2o3_beta`.
//
// The owner's values and rulings of 2026-11-27 (unit68_slice4b_4c_values_2026-11-27.txt, outside
// the repo), each with its source; read by the engineer: K. Ghosh, U. Singisetti, J. Appl. Phys.
// 122, 035702 (2017), Table I (the [010] velocity-field law) and J. Appl. Phys. 124, 085707 (2018),
// Table I (the [010] impact ionization). Monoclinic; one value per parameter along [010] (Ghosh's
// y), the normal of a (010) wafer (owner ruling 2). n-type only: holes are minority-carrier
// placeholders (m_h 3.4, mu_p 10, no saturation law), and a deck refuses a node where acceptors
// exceed donors (owner ruling 7).
//
// The gap is the 300 K indirect gap, 4.43 eV, at every temperature (no bulk temperature law in
// hand; flagged in the citations). The electrons follow Barnes's transferred-electron law
// (field_mobility.hpp) with the Caughey-Thomas mobility of Jang et al. in Ghosh's mu0; Ghosh fitted
// it to about 0.5 MV/cm only.
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

[[nodiscard]] SemiconductorParameters ga2o3_beta_parameters();

}  // namespace NiTCAD::physics
