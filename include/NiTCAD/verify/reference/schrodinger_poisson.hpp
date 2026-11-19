// A 1D Schrodinger-Poisson reference for the inversion layer of a MOS capacitor (Unit 62 slice 0, QP14;
// DECISIONS.md Unit 62). The reference the density-gradient model's QV6 gate compares against: class H,
// independence I2 (it calls no NiTCAD product code), provenance PD (derived here from the governing equations, with
// a second method for every part, below). Double precision: its U_ref is measured by grid refinement.
//
// 1. Symmetric tridiagonal eigenproblems (tridiagonal_eigen). The eigenvalue count below x is the number of negative
//    pivots of the LDL^T factorisation of T - x I (Sylvester's law of inertia): q_1 = d_1 - x,
//    q_i = d_i - x - e_{i-1}^2 / q_{i-1}; a zero pivot is replaced by -eps (|d| + |e|) scale (it counts as negative,
//    the limit from below). Bisection on that count, from the Gershgorin interval, to the last representable midpoint,
//    gives each requested eigenvalue (lowest first). Each eigenvector is found by inverse iteration: (T - lambda I) x = b
//    solved by Gaussian elimination with partial pivoting (a tridiagonal LU with two superdiagonals), three steps from
//    b = (1, ..., 1), each step orthogonalised against the vectors already found (Gram-Schmidt, twice) and normalised.
//    Reported per eigenpair: ||T x - lambda x||_2 (for a symmetric T, |lambda - lambda_true| <= that residual).
//    Second methods (tests): the closed form of tridiag(-1, 2, -1) and jacobi_eigen.hpp on n <= 16.
//
// 2. The effective-mass Schrodinger equation with a position-dependent mass (BenDaniel-Duke form),
//        -a d/dz ( (1 / m(z)) d psi / dz ) + V(z) psi = E psi,   a = hbar^2 / (2 m0) [eV cm^2], m in m0,
//    on nodes z_0 < ... < z_N with psi = 0 at z_0 and z_N, the mass and the potential constant per cell (a cell is
//    [z_i, z_{i+1}]; materials change at nodes). Box method at interior node i (cells h_{i-1}, h_i, box
//    w_i = (h_{i-1} + h_i) / 2):
//        -a [ (psi_{i+1} - psi_i) / (m_i h_i) - (psi_i - psi_{i-1}) / (m_{i-1} h_{i-1}) ]
//            + (V_{i-1} h_{i-1} + V_i h_i) / 2 psi_i = E w_i psi_i,
//    which keeps (1/m) psi' continuous across a cell boundary (the BenDaniel-Duke condition). With y = W^{1/2} psi
//    the generalised problem becomes the symmetric tridiagonal A = W^{-1/2} K W^{-1/2}; psi_i = y_i / sqrt(w_i), so
//    sum_i w_i psi_i^2 = 1 (the discrete normalisation integral of |psi|^2 dz = 1). Second methods (tests): the
//    infinite well's a pi^2 k^2 / (m L^2); the finite well with a mass step against its transcendental equation
//    solved by Brent; the triangular well a^{1/3} F^{2/3} |a_k| against the Airy zeros printed in NIST DLMF
//    Table 9.9.1 (read 2026-11-19, provenance PA).
//
// 3. The MOS capacitor at equilibrium (mos_schrodinger_poisson). E_F = 0 everywhere; energies in eV, lengths in cm.
//    Gate face at z = -t_ox, interface at z = 0, substrate to z = L. Potential phi [V], phi(L) = 0 (neutral bulk),
//    phi(-t_ox) = V_G - V_FB (so phi = 0 everywhere at flat band; no oxide or interface charge). Bands:
//        silicon  E_c(z) = E_cb - phi(z),  E_cb = E_g - kT ln(N_v / N_A)  (the Boltzmann neutral p-type bulk,
//                                          p_b = N_A, full ionisation; the bulk's minority electrons neglected),
//        oxide    E_c(z) = E_cb + dE_c - phi(z),
//    holes classical and Boltzmann, p = N_A exp(-phi / kT); acceptors N_A fully ionised; no charge in the oxide.
//    Electrons, by input:
//      - classical:  n = N_c exp(-E_c / kT) (Boltzmann; the mos_capacitor.hpp reference is its second method);
//      - quantum:    on the Schrodinger domain [-t_ox, z_q] (z_q = quantum_depth; psi = 0 at both ends; in the
//                    oxide every valley takes the oxide mass), for each valley v (quantisation mass m_z, density-of-
//                    states mass m_d, degeneracy g) and subband k,
//                        n(z) = sum_{v,k} g_v (m_d,v kT / (pi hbar^2)) F(-E_vk / kT) |psi_vk(z)|^2,
//                    F = exp (Boltzmann) or ln(1 + e^x) (Fermi-Dirac: the 2D density of states integrated); outside
//                    [0, z_q] in the silicon n is the classical value (negligible in a p-type bulk). m_d kT/(pi hbar^2)
//                    = m_d kT / (2 pi a).
//    Poisson, box method, the permittivity per cell, the charge per node over its silicon share s_i of the box:
//        -[ eps_i (phi_{i+1} - phi_i) / h_i - eps_{i-1} (phi_i - phi_{i-1}) / h_{i-1} ] = q (p_i - n_i - N_A) s_i.
//    Self-consistency: outer iterations. Given phi_old, the Schrodinger problems are solved on phi_old; then Poisson
//    is solved by Newton with the predicted density n_i(phi) = sum c_vk |psi_vk,i|^2 F((-E_vk + phi_i - phi_old,i) /
//    kT) (a rise of phi lowers the local band and fills the subbands). The outer loop stops when
//    max |phi - phi_old| <= outer_tolerance_V. At that fixed point the predicted and the quantum density coincide,
//    so the result does not depend on the predictor (it only speeds convergence); the reported residual is that of
//    Poisson with the quantum density of the final phi's own subbands. The start is the classical solution.
//    Outputs: the inversion density N_inv = sum n_i s_i over the silicon [cm^-2], its centroid
//    sum z_i n_i s_i / N_inv [cm], the semiconductor charge Q_s = q sum (p_i - n_i - N_A) s_i [C/cm^2], the oxide
//    displacement D_ox = eps_ox (phi(-t_ox) - phi(0)) / t_ox, whose sum with Q_s Gauss's law makes zero, and the
//    subband energies (relative to E_F) per valley.
//
// Approximations named (the H class): parabolic isotropic-in-plane valleys; one oxide mass for every valley; hard
// walls at the gate face and at z_q; holes classical; equilibrium; no image or exchange-correlation potential; the
// electrons' tail inside the oxide carries no charge (it is in no node's silicon share: not in Poisson or N_inv).
#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <string>
#include <vector>

namespace NiTCAD::verify::reference {

struct TridiagonalEigen {
    std::vector<double> values;                // ascending, the lowest `count`
    std::vector<std::vector<double>> vectors;  // unit 2-norm
    std::vector<double> residuals;             // ||T x - lambda x||_2 per pair
};

// d: the diagonal (n), e: the off-diagonal (n - 1). Errors: n == 0, e.size() != n - 1, a non-finite entry, count == 0
// or count > n.
[[nodiscard]] std::expected<TridiagonalEigen, std::string> tridiagonal_eigen(std::span<const double> d,
                                                                             std::span<const double> e,
                                                                             std::size_t count);

// The number of eigenvalues of the tridiagonal (d, e) strictly below x (the Sturm count).
[[nodiscard]] std::size_t sturm_count(std::span<const double> d, std::span<const double> e, double x);

struct Subbands {
    std::vector<double> energies_eV;            // ascending
    std::vector<std::vector<double>> psi;       // per subband, at every node (0 at both ends) [cm^-1/2]
    std::vector<double> residuals;              // the eigen-residual of each, in eV
};

// nodes z (N + 1, strictly increasing) [cm]; per cell (N): potential V [eV] and mass [m0]; a = hbar^2/(2 m0)
// [eV cm^2]. Errors: fewer than 3 nodes, sizes, non-increasing z, a non-finite or non-positive mass, a <= 0, count
// out of range (at most N - 1).
[[nodiscard]] std::expected<Subbands, std::string> schrodinger_1d(std::span<const double> z_cm,
                                                                  std::span<const double> cell_potential_eV,
                                                                  std::span<const double> cell_mass,
                                                                  double a_eV_cm2, std::size_t count);

struct Valley {
    double mass_z = 0.0;    // quantisation mass [m0]
    double mass_dos = 0.0;  // 2D density-of-states mass [m0]
    double degeneracy = 0.0;
};

enum class SpElectrons { classical, quantum };
enum class SpStatistics { boltzmann, fermi_dirac };

struct MosSpInput {
    // Constants.
    double q_C = 0.0;
    double kT_eV = 0.0;
    double eps0_F_per_cm = 0.0;
    double a_eV_cm2 = 0.0;  // hbar^2 / (2 m0)
    // Materials.
    double eps_si_r = 0.0, eps_ox_r = 0.0;
    double band_gap_eV = 0.0;
    double conduction_dos_cm3 = 0.0, valence_dos_cm3 = 0.0;  // N_c (classical electrons), N_v (holes)
    double oxide_barrier_eV = 0.0;                            // dE_c
    double oxide_mass = 0.0;                                  // [m0]
    std::vector<Valley> valleys;
    // Device.
    double acceptors_cm3 = 0.0;  // N_A > 0, with N_A << N_v
    double oxide_thickness_cm = 0.0;
    // Grid: n_oxide uniform oxide cells; silicon cells from h0 at the interface growing by `ratio` up to h_max, to
    // depth L; the Schrodinger domain ends at the first node at or past quantum_depth.
    std::size_t oxide_cells = 0;
    double h0_cm = 0.0, ratio = 1.0, h_max_cm = 0.0, depth_cm = 0.0, quantum_depth_cm = 0.0;
    // Model.
    SpElectrons electrons = SpElectrons::quantum;
    SpStatistics statistics = SpStatistics::boltzmann;
    std::size_t subbands_per_valley = 0;
    double outer_tolerance_V = 1e-10;
    int max_outer = 300;
};

struct MosSpResult {
    std::vector<double> z_cm, phi_V, n_cm3, p_cm3;
    std::size_t interface_node = 0;              // z = 0
    std::vector<std::vector<double>> subband_eV;  // per valley (quantum)
    double inversion_cm2 = 0.0;                   // N_inv
    double centroid_cm = 0.0;
    double charge_C_per_cm2 = 0.0;                // Q_s
    double oxide_displacement_C_per_cm2 = 0.0;    // D_ox
    double surface_potential_V = 0.0;             // phi(0)
    double poisson_residual = 0.0;                // max_i |F_i| / max_i(|q rho_i s_i| + flux scale), quantum density of the final phi
    double last_update_V = 0.0;                   // the last outer max |phi - phi_old|
    double max_eigen_residual_eV = 0.0;
    double orthonormality = 0.0;                  // max |<psi_j, psi_k>_w - delta_jk|
    int outer_iterations = 0;
};

// gate_V: V_G - V_FB. Errors: bad inputs; Newton or the outer loop not converging (named).
[[nodiscard]] std::expected<MosSpResult, std::string> mos_schrodinger_poisson(const MosSpInput& in, double gate_V);

}  // namespace NiTCAD::verify::reference
