// Private to the solve layer: the physical fields of an equilibrium-Poisson state, shared by
// solve_equilibrium and the quasi-static bias sweep.
#pragma once

#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::solve::detail {

// Potential in V and the Boltzmann carriers in cm^-3 of the scaled potential psi.
inline results::NodeFields equilibrium_fields(const assemble::EquilibriumPoisson& system,
                                              std::span<const double> psi,
                                              const assemble::Scaling& scaling) {
    const std::size_t n = system.unknowns();
    results::NodeFields f{std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    system.carriers(psi, f.n_cm3, f.p_cm3);
    for (std::size_t i = 0; i < n; ++i) {
        f.potential_V[i] = psi[i] * scaling.V_T;
        f.n_cm3[i] *= scaling.Ns;
        f.p_cm3[i] *= scaling.Ns;
    }
    return f;
}

// Scale from scaled gate charges to C / cm^(3-D): q Ns L_D^D.
inline double charge_scale(const assemble::Scaling& scaling, int dimension) {
    return base::q_C * scaling.Ns * std::pow(scaling.L_D, dimension);
}

}  // namespace NiTCAD::solve::detail
