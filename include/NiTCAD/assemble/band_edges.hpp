// The band diagram of a state, as both assemblers report it (Unit 15): per node, in units of V_T,
// measured from the equilibrium Fermi level (that of a contact at 0 V). With eta = psi + s (the
// band shift, scaled_device):
//     E_c = g_n - eta,  E_v = -(g_p + eta)          (g = ln(N / n_ie), narrowing shared),
//     E_Fn = ln(n / n_ie) - ln gamma_n - eta,  E_Fp = -eta - ln(p / n_ie) + ln gamma_p,
// with ln gamma 0 under Boltzmann statistics. psi depends on which material node 0 is in; psi + s,
// and so these, do not.
#pragma once

#include <vector>

namespace NiTCAD::assemble {

struct BandEdges {
    std::vector<double> conduction;
    std::vector<double> valence;
    std::vector<double> electron_fermi;
    std::vector<double> hole_fermi;
};

}  // namespace NiTCAD::assemble
