// Test helper: a port of the legacy 1D drift-diffusion diode (core/src/device1d/device1d.cpp,
// residual_jacobian and contact_value: Boltzmann statistics, Scharfetter-Gummel fluxes,
// Caughey-Thomas mobility, Scharfetter SRH, no band-gap narrowing, no Auger) and of the legacy
// backward-Euler transient (pytcad/transient.py solve_transient with theta = 1 and its
// _newton_step), so the Unit 21 turn-off gate compares NiTCAD with the legacy algorithm computed
// here rather than with stored output. Independent of NiTCAD's assemble and solve layers: its own
// residual, Jacobian, Newton loop and step loop; only the material model functions and the sparse
// LU of linalg are shared.
//
// Legacy algorithm, as ported:
// - Scaling: Ns = max |doping| (above n_i), L_D = sqrt(eps V_T / (q Ns)), J0 = q D0 Ns / L_D,
//   R0 = D0 Ns / L_D^2, D0 = 1 cm^2/s; t0 = Ns / R0 seconds per scaled time. Boxes dV: half cells
//   at the ends.
// - Interior rows: Poisson (psi_{i+1} - psi_i)/h_i - (psi_i - psi_{i-1})/h_{i-1} - dV (n - p - C);
//   electrons Jn_i - Jn_{i-1} - R dV; holes Jp_i - Jp_{i-1} + R dV, with
//   Jn = a_n (n_{k+1} B(d) - n_k B(-d)), Jp = -a_p (p_{k+1} B(-d) - p_k B(d)),
//   d = psi_{k+1} - psi_k, a = hmean(mu) V_T / (D0 h); B(x) = x / expm1(x) (x clipped to +-700; 1 - x/2 + x^2/12 below
//   1e-4) and B' = B (1/x - 1/(1 - e^-x)) (-1/2 + x/6 - x^3/180 below 1e-4).
//   R = (n p - n_i^2) / (tau_p (n + n_i) + tau_n (p + n_i)) on physical densities, / R0.
// - Contacts: Dirichlet rows psi - psi0, n - n0, p - p0, with n0 p0 = n_i^2, n0 - p0 = C (the
//   majority from the square root) and psi0 = V / V_T + ln(n0 / n_i).
// - Transient step (theta = 1): the interior electron rows lose dV (n - n_old) / dt, the hole rows
//   gain dV (p - p_old) / dt (dt scaled by t0).
// - Newton (_newton_step): du = J^-1 (-F); dpsi clipped to +-max_dpsi (5); merit 0.5 |F|^2,
//   lambda halved up to 40 times until merit <= base (1 - 1e-4 lambda), else lambda = 0; then
//   psi += lambda dpsi, n = clip(n + lambda dn, 0.1 n, 10 n) (and p). Converged when
//   max(|dpsi| (the clipped full correction), max |n / n_before - 1|, max |p / p_before - 1|) <
//   tol_update; at most max_iter iterations.
// - Terminal current of the left contact: (Jn_0 + Jp_0) J0 [A/cm^2].
// The steady starting state is found by the same Newton without the storage term, from the
// charge-neutral guess, the bias raised in steps of 0.1 V (the converged state does not depend
// on the path; the legacy solve_bias converges to the same state to its tolerance).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

struct LegacyTurnoffRun {
    double forward_current;          // at the starting bias [A/cm^2]
    double transit_time;             // (L - x_j)^2 / (2 D_p) at the cathode [s]
    double dt;                       // the fixed step [s]
    std::vector<double> time;        // [s], from 0
    std::vector<double> current;     // left-contact conduction current [A/cm^2]
    std::vector<int> iterations;     // Newton iterations of each step (0 for the start)
    std::vector<double> step_lambda; // the last iteration's lambda of each step
};

class LegacyDiode1D {
public:
    // doping: N_D - N_A per node [cm^-3]; x ascending [cm].
    LegacyDiode1D(std::vector<double> x, std::vector<double> doping, double T = 300.0)
        : x_(std::move(x)), doping_(std::move(doping)) {
        const NiTCAD::physics::Semiconductor si = NiTCAD::physics::silicon();
        N_ = x_.size();
        VT_ = NiTCAD::base::thermal_voltage(T);
        eps_ = si.parameters().eps_r * NiTCAD::base::eps0_F_per_cm;
        ni_ = NiTCAD::physics::intrinsic_density(si, T);
        double dmax = 0.0;
        for (const double d : doping_) dmax = std::max(dmax, std::abs(d));
        Ns_ = std::max(dmax, ni_);
        LD_ = std::sqrt(eps_ * VT_ / (NiTCAD::base::q_C * Ns_));
        J0_ = NiTCAD::base::q_C * 1.0 * Ns_ / LD_;
        R0_ = 1.0 * Ns_ / (LD_ * LD_);
        h_.resize(N_ - 1);
        for (std::size_t k = 0; k + 1 < N_; ++k) h_[k] = (x_[k + 1] - x_[k]) / LD_;
        dV_.assign(N_, 0.0);
        dV_[0] = 0.5 * h_[0];
        dV_[N_ - 1] = 0.5 * h_[N_ - 2];
        for (std::size_t i = 1; i + 1 < N_; ++i) dV_[i] = 0.5 * (h_[i - 1] + h_[i]);
        C_.resize(N_);
        mu_n_.resize(N_);
        mu_p_.resize(N_);
        tau_n_.resize(N_);
        tau_p_.resize(N_);
        for (std::size_t i = 0; i < N_; ++i) {
            C_[i] = doping_[i] / Ns_;
            const double N = std::abs(doping_[i]);
            using NiTCAD::physics::Carrier;
            mu_n_[i] = NiTCAD::physics::caughey_thomas_mobility(si, Carrier::electron, N, T);
            mu_p_[i] = NiTCAD::physics::caughey_thomas_mobility(si, Carrier::hole, N, T);
            tau_n_[i] = NiTCAD::physics::scharfetter_lifetime(si, Carrier::electron, N);
            tau_p_[i] = NiTCAD::physics::scharfetter_lifetime(si, Carrier::hole, N);
        }
        an_.resize(N_ - 1);
        ap_.resize(N_ - 1);
        for (std::size_t k = 0; k + 1 < N_; ++k) {
            const auto hmean = [](double a, double b) { return 2.0 * a * b / (a + b); };
            an_[k] = hmean(mu_n_[k], mu_n_[k + 1]) * VT_ / 1.0 / h_[k];
            ap_[k] = hmean(mu_p_[k], mu_p_[k + 1]) * VT_ / 1.0 / h_[k];
        }
        nie_s_ = ni_ / Ns_;
        // The Jacobian pattern: 3x3 blocks on the diagonal and to each neighbour.
        std::vector<NiTCAD::linalg::Triplet> t;
        for (std::size_t i = 0; i < N_; ++i) {
            for (std::size_t j = (i == 0 ? 0 : i - 1); j <= std::min(i + 1, N_ - 1); ++j) {
                for (int r = 0; r < 3; ++r) {
                    for (int c = 0; c < 3; ++c) {
                        t.push_back({static_cast<NiTCAD::linalg::Index>(3 * i + r),
                                     static_cast<NiTCAD::linalg::Index>(3 * j + c), 0.0});
                    }
                }
            }
        }
        const auto size = static_cast<NiTCAD::linalg::Index>(3 * N_);
        pattern_ = *NiTCAD::linalg::SparseMatrix::from_triplets(size, size, t);
    }

    [[nodiscard]] double time_scale() const { return Ns_ / R0_; }
    [[nodiscard]] double cathode_hole_diffusivity() const { return VT_ * mu_p_[N_ - 1]; }

    // The legacy fixture: backward Euler at fixed dt = transit_time / steps_per_transit from the
    // steady state at v_forward, the anode switched to 0 V at t = 0, for `transits` transit
    // times.
    [[nodiscard]] LegacyTurnoffRun turnoff(double v_forward, double length, double junction,
                                           int steps_per_transit, int transits,
                                           int max_iter = 50, double tol_update = 1e-8) {
        State s = neutral_state();
        for (double v = 0.0;; v = std::min(v + 0.1, v_forward)) {
            const Contacts bc = contacts(v);
            stamp(s, bc);
            if (!newton(s, bc, nullptr, 0.0, 200, tol_update).converged) {
                throw std::runtime_error("legacy port: the steady state did not converge");
            }
            if (v == v_forward) break;
        }
        LegacyTurnoffRun run{};
        const Contacts forward = contacts(v_forward);
        run.forward_current = left_current(s, forward);
        const double Dp = cathode_hole_diffusivity();
        run.transit_time = (length - junction) * (length - junction) / (2.0 * Dp);
        run.dt = run.transit_time / steps_per_transit;
        run.time.push_back(0.0);
        run.current.push_back(run.forward_current);
        run.iterations.push_back(0);
        run.step_lambda.push_back(0.0);
        const Contacts off = contacts(0.0);
        const double dt_s = run.dt / time_scale();
        const int steps = steps_per_transit * transits;
        for (int k = 1; k <= steps; ++k) {
            const State old = s;
            const Outcome o = newton(s, off, &old, dt_s, max_iter, tol_update);
            if (!o.converged) throw std::runtime_error("legacy port: a step did not converge");
            run.time.push_back(k * run.dt);
            run.current.push_back(left_current(s, off));
            run.iterations.push_back(o.iterations);
            run.step_lambda.push_back(o.lambda);
        }
        return run;
    }

private:
    struct State {
        std::vector<double> psi, n, p;
    };
    struct Contacts {
        double psi0[2], n0[2], p0[2];
    };
    struct Outcome {
        bool converged;
        int iterations;
        double lambda;
    };

    static double clip700(double x) { return std::clamp(x, -700.0, 700.0); }
    static double bernoulli(double x_in) {
        const double x = clip700(x_in);
        if (std::abs(x) < 1e-4) return 1.0 - x / 2.0 + x * x / 12.0;
        return x / std::expm1(x);
    }
    static double dbernoulli(double x_in) {
        const double x = clip700(x_in);
        if (std::abs(x) < 1e-4) return -0.5 + x / 6.0 - x * x * x / 180.0;
        const double B = x / std::expm1(x);
        double em = -std::expm1(-x);
        if (std::abs(em) < 1e-300) em = 1e-300;
        return B * (1.0 / x - 1.0 / em);
    }

    [[nodiscard]] Contacts contacts(double v_left) const {
        Contacts c{};
        for (int side = 0; side < 2; ++side) {
            const std::size_t i = side == 0 ? 0 : N_ - 1;
            const double Cc = C_[i], nie = nie_s_;
            const double root = std::sqrt(Cc * Cc + 4.0 * nie * nie);
            double n0, p0;
            if (Cc >= 0.0) {
                n0 = 0.5 * (Cc + root);
                p0 = nie * nie / n0;
            } else {
                p0 = 0.5 * (-Cc + root);
                n0 = nie * nie / p0;
            }
            const double V = side == 0 ? v_left : 0.0;
            c.psi0[side] = V / VT_ + std::log(n0 / nie);
            c.n0[side] = n0;
            c.p0[side] = p0;
        }
        return c;
    }

    [[nodiscard]] State neutral_state() const {
        State s{std::vector<double>(N_), std::vector<double>(N_), std::vector<double>(N_)};
        for (std::size_t i = 0; i < N_; ++i) {
            const double Cc = C_[i], nie = nie_s_;
            const double root = std::sqrt(Cc * Cc + 4.0 * nie * nie);
            const double n0 = Cc >= 0.0 ? 0.5 * (Cc + root) : nie * nie / (0.5 * (-Cc + root));
            s.n[i] = n0;
            s.p[i] = nie * nie / n0;
            s.psi[i] = std::log(n0 / nie);
        }
        return s;
    }

    void stamp(State& s, const Contacts& c) const {
        for (int side = 0; side < 2; ++side) {
            const std::size_t i = side == 0 ? 0 : N_ - 1;
            s.psi[i] = c.psi0[side];
            s.n[i] = c.n0[side];
            s.p[i] = c.p0[side];
        }
    }

    // Edge fluxes of edge k (scaled).
    void fluxes(const State& s, std::size_t k, double& Jn, double& Jp) const {
        const double d = s.psi[k + 1] - s.psi[k];
        Jn = an_[k] * (s.n[k + 1] * bernoulli(d) - s.n[k] * bernoulli(-d));
        Jp = -ap_[k] * (s.p[k + 1] * bernoulli(-d) - s.p[k] * bernoulli(d));
    }

    [[nodiscard]] double left_current(const State& s, const Contacts&) const {
        double Jn, Jp;
        fluxes(s, 0, Jn, Jp);
        return (Jn + Jp) * J0_;
    }

    // F and (when jac) the Jacobian values at s; with old, the backward-Euler storage of step dt.
    void residual(const State& s, const Contacts& c, const State* old, double dt,
                  std::vector<double>& F, NiTCAD::linalg::SparseMatrix* jac) const {
        const std::size_t N = N_;
        F.assign(3 * N, 0.0);
        std::vector<double> Jn(N - 1), Jp(N - 1), dJn(N - 1), dJp(N - 1);
        std::vector<double> Bp(N - 1), Bm(N - 1);
        for (std::size_t k = 0; k + 1 < N; ++k) {
            const double d = s.psi[k + 1] - s.psi[k];
            Bp[k] = bernoulli(d);
            Bm[k] = bernoulli(-d);
            const double dBp = dbernoulli(d), dBm = dbernoulli(-d);
            Jn[k] = an_[k] * (s.n[k + 1] * Bp[k] - s.n[k] * Bm[k]);
            Jp[k] = -ap_[k] * (s.p[k + 1] * Bm[k] - s.p[k] * Bp[k]);
            dJn[k] = an_[k] * (s.n[k + 1] * dBp + s.n[k] * dBm);  // d Jn / d psi_{k+1}
            dJp[k] = ap_[k] * (s.p[k + 1] * dBm + s.p[k] * dBp);
        }
        std::span<double> v;
        if (jac != nullptr) {
            v = jac->values();
            std::fill(v.begin(), v.end(), 0.0);
        }
        const auto add = [&](std::size_t r, std::size_t col, double value) {
            if (jac == nullptr) return;
            v[position(r, col)] += value;
        };
        for (std::size_t i = 1; i + 1 < N; ++i) {
            const std::size_t eR = i, eL = i - 1;
            F[3 * i] = (s.psi[i + 1] - s.psi[i]) / h_[eR] - (s.psi[i] - s.psi[i - 1]) / h_[eL] -
                       dV_[i] * (s.n[i] - s.p[i] - C_[i]);
            add(3 * i, 3 * i, -1.0 / h_[eR] - 1.0 / h_[eL]);
            add(3 * i, 3 * (i + 1), 1.0 / h_[eR]);
            add(3 * i, 3 * (i - 1), 1.0 / h_[eL]);
            add(3 * i, 3 * i + 1, -dV_[i]);
            add(3 * i, 3 * i + 2, dV_[i]);
            const NiTCAD::physics::RecombinationRate rr = srh(s.n[i] * Ns_, s.p[i] * Ns_, i);
            const double Rs = rr.rate / R0_, dRn = rr.d_dn * Ns_ / R0_, dRp = rr.d_dp * Ns_ / R0_;
            F[3 * i + 1] = Jn[eR] - Jn[eL] - Rs * dV_[i];
            add(3 * i + 1, 3 * i + 1, -an_[eR] * Bm[eR] - an_[eL] * Bp[eL] - dRn * dV_[i]);
            add(3 * i + 1, 3 * (i + 1) + 1, an_[eR] * Bp[eR]);
            add(3 * i + 1, 3 * (i - 1) + 1, an_[eL] * Bm[eL]);
            add(3 * i + 1, 3 * i + 2, -dRp * dV_[i]);
            add(3 * i + 1, 3 * i, -dJn[eR] - dJn[eL]);
            add(3 * i + 1, 3 * (i + 1), dJn[eR]);
            add(3 * i + 1, 3 * (i - 1), dJn[eL]);
            F[3 * i + 2] = Jp[eR] - Jp[eL] + Rs * dV_[i];
            add(3 * i + 2, 3 * i + 2, ap_[eR] * Bp[eR] + ap_[eL] * Bm[eL] + dRp * dV_[i]);
            add(3 * i + 2, 3 * (i + 1) + 2, -ap_[eR] * Bm[eR]);
            add(3 * i + 2, 3 * (i - 1) + 2, -ap_[eL] * Bp[eL]);
            add(3 * i + 2, 3 * i + 1, dRn * dV_[i]);
            add(3 * i + 2, 3 * i, -dJp[eR] - dJp[eL]);
            add(3 * i + 2, 3 * (i + 1), dJp[eR]);
            add(3 * i + 2, 3 * (i - 1), dJp[eL]);
            if (old != nullptr) {
                F[3 * i + 1] -= dV_[i] / dt * (s.n[i] - old->n[i]);
                F[3 * i + 2] += dV_[i] / dt * (s.p[i] - old->p[i]);
                add(3 * i + 1, 3 * i + 1, -dV_[i] / dt);
                add(3 * i + 2, 3 * i + 2, dV_[i] / dt);
            }
        }
        for (int side = 0; side < 2; ++side) {
            const std::size_t i = side == 0 ? 0 : N - 1;
            F[3 * i] = s.psi[i] - c.psi0[side];
            F[3 * i + 1] = s.n[i] - c.n0[side];
            F[3 * i + 2] = s.p[i] - c.p0[side];
            for (int r = 0; r < 3; ++r) add(3 * i + r, 3 * i + r, 1.0);
        }
    }

    [[nodiscard]] NiTCAD::physics::RecombinationRate srh(double n, double p,
                                                         std::size_t i) const {
        const double excess = n * p - ni_ * ni_;
        const double den = tau_p_[i] * (n + ni_) + tau_n_[i] * (p + ni_);
        return {excess / den, (p * den - excess * tau_p_[i]) / (den * den),
                (n * den - excess * tau_n_[i]) / (den * den)};
    }

    [[nodiscard]] std::size_t position(std::size_t r, std::size_t c) const {
        const auto& off = pattern_.row_offsets();
        const auto& col = pattern_.col_indices();
        const auto b = col.begin() + off[r], e = col.begin() + off[r + 1];
        const auto it = std::lower_bound(b, e, static_cast<NiTCAD::linalg::Index>(c));
        return static_cast<std::size_t>(it - col.begin());
    }

    [[nodiscard]] double merit(const State& s, const Contacts& c, const State* old,
                               double dt) const {
        std::vector<double> F;
        residual(s, c, old, dt, F, nullptr);
        double m = 0.0;
        for (const double f : F) m += f * f;
        return 0.5 * m;
    }

    // The legacy _newton_step.
    Outcome newton(State& s, const Contacts& c, const State* old, double dt, int max_iter,
                   double tol) const {
        NiTCAD::linalg::SparseMatrix J = pattern_;
        auto solver = *NiTCAD::linalg::LinearSolver::create({});
        std::vector<double> F, rhs(3 * N_), du(3 * N_);
        double lambda = 0.0;
        for (int it = 0; it < max_iter; ++it) {
            residual(s, c, old, dt, F, &J);
            for (std::size_t k = 0; k < F.size(); ++k) rhs[k] = -F[k];
            if (!solver.factorize(J) || !solver.solve(rhs, du)) {
                throw std::runtime_error("legacy port: singular Jacobian");
            }
            std::vector<double> dpsi(N_), dn(N_), dp(N_);
            for (std::size_t i = 0; i < N_; ++i) {
                dpsi[i] = std::clamp(du[3 * i], -5.0, 5.0);
                dn[i] = du[3 * i + 1];
                dp[i] = du[3 * i + 2];
            }
            double base = 0.0;
            for (const double f : F) base += f * f;
            base *= 0.5;
            const auto trial = [&](double lam) {
                State t = s;
                for (std::size_t i = 0; i < N_; ++i) {
                    t.psi[i] = s.psi[i] + lam * dpsi[i];
                    t.n[i] = std::clamp(s.n[i] + lam * dn[i], 0.1 * s.n[i], 10.0 * s.n[i]);
                    t.p[i] = std::clamp(s.p[i] + lam * dp[i], 0.1 * s.p[i], 10.0 * s.p[i]);
                }
                return t;
            };
            lambda = 1.0;
            bool accepted = false;
            for (int k = 0; k < 40; ++k) {
                const double m = merit(trial(lambda), c, old, dt);
                if (std::isfinite(m) && m <= base * (1.0 - 1e-4 * lambda)) {
                    accepted = true;
                    break;
                }
                lambda *= 0.5;
            }
            if (!accepted) lambda = 0.0;
            const State before = s;
            s = trial(lambda);
            double err = 0.0;
            for (std::size_t i = 0; i < N_; ++i) {
                err = std::max({err, std::abs(dpsi[i]),
                                std::abs(s.n[i] / std::max(before.n[i], 1e-300) - 1.0),
                                std::abs(s.p[i] / std::max(before.p[i], 1e-300) - 1.0)});
            }
            if (err < tol) return {true, it + 1, lambda};
        }
        return {false, max_iter, lambda};
    }

    std::vector<double> x_, doping_;
    std::size_t N_ = 0;
    double VT_ = 0, eps_ = 0, ni_ = 0, Ns_ = 0, LD_ = 0, J0_ = 0, R0_ = 0, nie_s_ = 0;
    std::vector<double> h_, dV_, C_, mu_n_, mu_p_, tau_n_, tau_p_, an_, ap_;
    NiTCAD::linalg::SparseMatrix pattern_;
};
