// Thermionic emission at an abrupt heterointerface (ARCHITECTURE.md section 5, item 5, an
// interface-local model; Unit 15; legacy device.py emission_velocity, M33-S2).
//
// The emission velocity of a band with effective density of states N:
//     v = sqrt(k T / (2 pi m)),   m from N = 2 (2 pi m k T / h^2)^(3/2),
// which simplifies to v = (k T / h) (2 / N)^(1/3). It is the A* T^2 / (q N) of Richardson's law
// with A* = 4 pi q m k^2 / h^3 written through the mass that N itself implies, so it needs no new
// material constant and stays consistent with the Nc and Nv the equations use. The legacy records
// what that costs: for silicon (Nc = 2.86e19 cm^-3 at 300 K, a density-of-states mass of 1.09 m0)
// it gives 2.575e6 cm/s, where the tabulated Richardson constant (2.1 m0, six valleys) gives a
// factor 1.92 more; the velocity is good to about a factor of two on silicon.
#pragma once

namespace NiTCAD::physics {

// v [cm/s] for a band with effective density of states dos_cm3 [cm^-3] at temperature_K.
// Preconditions (NITCAD_EXPECTS): both finite and positive.
[[nodiscard]] double emission_velocity_cm_s(double dos_cm3, double temperature_K);

}  // namespace NiTCAD::physics
