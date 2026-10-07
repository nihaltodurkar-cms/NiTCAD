// Test helper: a port of the legacy pytcad/moscap.py MOSCapacitor (Boltzmann statistics, no
// density gradient; the M14 interface-trap term with set_interface_traps, Unit 15b), so the Unit 12
// gates can compare NiTCAD's general gate
// contact with the legacy 1D MOS-C solve on the legacy mesh, and port the legacy C-V checks
// (tests/test_cv_physics_validation.py). Independent of NiTCAD's assemble and solve layers: its own
// residual, its own tridiagonal Newton. Only the material model functions are shared.
//
// Legacy algorithm: scaled Poisson on x = L expm1(6 s) / expm1(6), s uniform on [0, 1]; row 0 is
// the half box with the gate flux kappa (Vg - Vfb - (psi_0 - psi_b)) entering; the last node is
// held at psi_b; Newton with each correction clipped to +-3, stopping when max |d| < tol after the
// update. C-V: phi_s = (psi_0 - psi_b) V_T, Qg = C_ox (Vg - V_FB - phi_s), C = numpy.gradient(Qg).
// M14 interface traps: row 0 loses dit_coeff (psi_0 - psi_b), dit_coeff = q D_it L_D / eps_s (the
// trapped charge -q D_it phi_s, zero at flat band), and Qg loses q D_it phi_s.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

struct LegacyGate {
    enum Kind { n_poly, p_poly, metal } kind = n_poly;
    double work_function_eV = 0.0;  // metal only
};

inline constexpr double legacy_eps_ox_r = 3.9;  // moscap.EPS_OX_R

// numpy.gradient(f, x) with the default edge_order = 1. Like numpy, it needs at least two points.
inline std::vector<double> numpy_gradient(const std::vector<double>& f,
                                          const std::vector<double>& x) {
    const std::size_t n = f.size();
    if (n < 2 || x.size() != n) throw std::invalid_argument("numpy_gradient: needs two points");
    std::vector<double> g(n);
    g[0] = (f[1] - f[0]) / (x[1] - x[0]);
    g[n - 1] = (f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2]);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        const double hs = x[i] - x[i - 1], hd = x[i + 1] - x[i];
        g[i] = (hs * hs * f[i + 1] + (hd * hd - hs * hs) * f[i] - hd * hd * f[i - 1]) /
               (hs * hd * (hd + hs));
    }
    return g;
}

// moscap.flatband_voltage [V], for silicon.
inline double legacy_flatband_voltage(double Nsub, double tox_cm, LegacyGate gate, double Qf,
                                      double T) {
    using namespace NiTCAD;
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T);
    const double Cox = legacy_eps_ox_r * base::eps0_F_per_cm / tox_cm;
    const double ni = physics::intrinsic_density(si, T);
    const double nie = physics::effective_intrinsic_density(si, std::abs(Nsub), T);
    const double Ns = std::max(std::abs(Nsub), ni);
    const double psi_b = std::asinh((Nsub / Ns) / (2.0 * (nie / Ns)));
    const double chi = si.parameters().electron_affinity_eV;
    const double Eg = physics::band_gap_eV(si, T);
    const double phi_m = gate.kind == LegacyGate::n_poly   ? chi
                         : gate.kind == LegacyGate::p_poly ? chi + Eg
                                                           : gate.work_function_eV;
    const double phi_semi = chi + 0.5 * Eg - psi_b * VT;
    return (phi_m - phi_semi) - base::q_C * Qf / Cox;
}

struct LegacyCV {
    std::vector<double> phi_s, Qg, C;
};

class LegacyMOSCapacitor {
public:
    LegacyMOSCapacitor(double Nsub, double tox_cm, LegacyGate gate = {}, double Qf = 0.0,
                       double T = 300.0, double L_cm = 2e-4, std::size_t nx = 1200)
        : Nsub(Nsub), T(T) {
        using namespace NiTCAD;
        const physics::Semiconductor si = physics::silicon();
        VT = base::thermal_voltage(T);
        eps_s = si.parameters().eps_r * base::eps0_F_per_cm;
        const double eps_ox = legacy_eps_ox_r * base::eps0_F_per_cm;
        Cox = eps_ox / tox_cm;
        ni = physics::intrinsic_density(si, T);
        const double nie = physics::effective_intrinsic_density(si, std::abs(Nsub), T);
        Ns = std::max(std::abs(Nsub), ni);
        LD = std::sqrt(eps_s * VT / (base::q_C * Ns));
        x.resize(nx);
        xs.resize(nx);
        for (std::size_t i = 0; i < nx; ++i) {
            const double s = static_cast<double>(i) / static_cast<double>(nx - 1);
            x[i] = L_cm * (std::expm1(6.0 * s) / std::expm1(6.0));
            xs[i] = x[i] / LD;
        }
        C = Nsub / Ns;
        nie_s = nie / Ns;
        psi_b = std::asinh(C / (2.0 * nie_s));
        kappa = eps_ox * LD / (eps_s * tox_cm);
        Vfb = legacy_flatband_voltage(Nsub, tox_cm, gate, Qf, T);
    }

    // The legacy residual of the scaled potential at gate bias Vg (all rows).
    [[nodiscard]] std::vector<double> residual(const std::vector<double>& psi, double Vg) const {
        const std::size_t n = x.size();
        std::vector<double> F(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double h_left = i > 0 ? xs[i] - xs[i - 1] : 0.0;
            const double h_right = i + 1 < n ? xs[i + 1] - xs[i] : 0.0;
            const double dV = 0.5 * (h_left + h_right);
            const double e = std::clamp(psi[i], -700.0, 700.0);
            const double rho = nie_s * std::exp(e) - nie_s * std::exp(-e) - C;
            if (i == 0) {
                F[0] = (psi[1] - psi[0]) / h_right +
                       kappa * (Vg / VT - Vfb / VT - (psi[0] - psi_b)) - dV * rho;
                if (D_it != 0.0) F[0] -= dit_coeff * (psi[0] - psi_b);
            } else if (i + 1 == n) {
                F[i] = psi[i] - psi_b;
            } else {
                F[i] = (psi[i + 1] - psi[i]) / h_right - (psi[i] - psi[i - 1]) / h_left - dV * rho;
            }
        }
        return F;
    }

    // moscap.solve_psi: the scaled potential at Vg, from psi0 or the bulk value.
    [[nodiscard]] std::vector<double> solve_psi(double Vg,
                                                const std::vector<double>* psi0 = nullptr,
                                                int max_iter = 200, double tol = 1e-10) const {
        const std::size_t n = x.size();
        std::vector<double> psi = psi0 != nullptr ? *psi0 : std::vector<double>(n, psi_b);
        psi[n - 1] = psi_b;
        std::vector<double> lo(n), main(n), up(n), d(n);
        for (int it = 0; it < max_iter; ++it) {
            const std::vector<double> F = residual(psi, Vg);
            for (std::size_t i = 0; i < n; ++i) {
                const double h_left = i > 0 ? xs[i] - xs[i - 1] : 0.0;
                const double h_right = i + 1 < n ? xs[i + 1] - xs[i] : 0.0;
                const double dV = 0.5 * (h_left + h_right);
                const double e = std::clamp(psi[i], -700.0, 700.0);
                const double dnp = nie_s * std::exp(e) + nie_s * std::exp(-e);
                lo[i] = up[i] = 0.0;
                if (i == 0) {
                    main[0] = -1.0 / h_right - kappa - dV * dnp;
                    if (D_it != 0.0) main[0] -= dit_coeff;
                    up[0] = 1.0 / h_right;
                } else if (i + 1 == n) {
                    main[i] = 1.0;
                } else {
                    lo[i] = 1.0 / h_left;
                    main[i] = -1.0 / h_right - 1.0 / h_left - dV * dnp;
                    up[i] = 1.0 / h_right;
                }
                d[i] = -F[i];
            }
            // Thomas algorithm (the matrix is diagonally dominant).
            for (std::size_t i = 1; i < n; ++i) {
                const double w = lo[i] / main[i - 1];
                main[i] -= w * up[i - 1];
                d[i] -= w * d[i - 1];
            }
            d[n - 1] /= main[n - 1];
            for (std::size_t i = n - 1; i-- > 0;) d[i] = (d[i] - up[i] * d[i + 1]) / main[i];
            double largest = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                const double step = std::clamp(d[i], -3.0, 3.0);
                psi[i] += step;
                largest = std::max(largest, std::abs(step));
            }
            if (largest < tol) break;
        }
        return psi;
    }

    // moscap.cv_sweep: each point starts from the previous one.
    [[nodiscard]] LegacyCV cv_sweep(const std::vector<double>& Vg) const {
        LegacyCV r;
        std::vector<double> guess;
        for (const double v : Vg) {
            const std::vector<double> psi = solve_psi(v, guess.empty() ? nullptr : &guess);
            guess = psi;
            const double ps = (psi[0] - psi_b) * VT;
            r.phi_s.push_back(ps);
            r.Qg.push_back(D_it != 0.0 ? Cox * (v - Vfb - ps) - NiTCAD::base::q_C * D_it * ps
                                       : Cox * (v - Vfb - ps));
        }
        r.C = numpy_gradient(r.Qg, Vg);
        return r;
    }

    // moscap.MOSCapacitor(D_it=...): the M14 interface-trap density [cm^-2 eV^-1].
    void set_interface_traps(double density_cm2_eV) {
        D_it = density_cm2_eV;
        dit_coeff = NiTCAD::base::q_C * D_it * LD / eps_s;
    }

    struct Landmarks {
        double phi_F, W_max, C_ox, C_min, V_th, V_FB;
    };
    // moscap.analytic_landmarks: depletion-approximation landmarks.
    [[nodiscard]] Landmarks analytic_landmarks() const {
        const double N = std::abs(Nsub);
        const double phiF = VT * std::log(N / ni);
        const double Wmax = std::sqrt(4.0 * eps_s * phiF / (NiTCAD::base::q_C * N));
        const double Cmin = 1.0 / (1.0 / Cox + Wmax / eps_s);
        const double sign = Nsub < 0 ? 1.0 : -1.0;
        const double Vth =
            Vfb + sign * (2 * phiF + std::sqrt(4.0 * eps_s * NiTCAD::base::q_C * N * phiF) / Cox);
        return {phiF, Wmax, Cox, Cmin, Vth, Vfb};
    }

    double Nsub, T, VT = 0, eps_s = 0, Cox = 0, ni = 0, Ns = 0, LD = 0;
    std::vector<double> x, xs;
    double C = 0, nie_s = 0, psi_b = 0, kappa = 0, Vfb = 0;
    double D_it = 0.0, dit_coeff = 0.0;
};
