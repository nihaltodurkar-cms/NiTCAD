// The 1D Schrodinger-Poisson reference (schrodinger_poisson.hpp; Unit 62 slice 0, QP14).
#include "NiTCAD/verify/reference/schrodinger_poisson.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace NiTCAD::verify::reference {

namespace {

constexpr double eps = std::numeric_limits<double>::epsilon();

bool all_finite(std::span<const double> v) {
    return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
}

// The scale of the matrix: max |d_i| + 2 max |e_i| (an upper bound of ||T||_2 by Gershgorin), at least the
// smallest normal number.
double matrix_scale(std::span<const double> d, std::span<const double> e) {
    double md = 0.0, me = 0.0;
    for (const double x : d) md = std::max(md, std::abs(x));
    for (const double x : e) me = std::max(me, std::abs(x));
    return std::max(md + 2.0 * me, std::numeric_limits<double>::min());
}

std::size_t sturm_count_scaled(std::span<const double> d, std::span<const double> e, double x, double pivmin) {
    std::size_t count = 0;
    double q = d[0] - x;
    if (std::abs(q) < pivmin) q = -pivmin;
    if (q < 0.0) ++count;
    for (std::size_t i = 1; i < d.size(); ++i) {
        q = d[i] - x - e[i - 1] * e[i - 1] / q;
        if (std::abs(q) < pivmin) q = -pivmin;
        if (q < 0.0) ++count;
    }
    return count;
}

// (T - shift I) x = b by Gaussian elimination with partial pivoting (the LAPACK dgttrf / dgtts2 scheme); a zero
// pivot is replaced by `pivmin` (inverse iteration at an eigenvalue: the solve is meant to be near-singular).
std::vector<double> shifted_solve(std::span<const double> d0, std::span<const double> e, double shift,
                                  std::vector<double> b, double pivmin) {
    const std::size_t n = d0.size();
    std::vector<double> d(n), dl(e.begin(), e.end()), du(e.begin(), e.end()), du2(n > 2 ? n - 2 : 0, 0.0);
    std::vector<char> swapped(n > 1 ? n - 1 : 0, 0);
    for (std::size_t i = 0; i < n; ++i) d[i] = d0[i] - shift;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        if (std::abs(d[i]) >= std::abs(dl[i])) {
            if (d[i] == 0.0) d[i] = pivmin;
            const double fact = dl[i] / d[i];
            dl[i] = fact;
            d[i + 1] -= fact * du[i];
        } else {
            const double fact = d[i] / dl[i];
            d[i] = dl[i];
            dl[i] = fact;
            const double temp = du[i];
            du[i] = d[i + 1];
            d[i + 1] = temp - fact * d[i + 1];
            if (i + 2 < n) {
                du2[i] = du[i + 1];
                du[i + 1] = -fact * du[i + 1];
            }
            swapped[i] = 1;
        }
    }
    if (d[n - 1] == 0.0) d[n - 1] = pivmin;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        if (swapped[i] == 0) {
            b[i + 1] -= dl[i] * b[i];
        } else {
            const double temp = b[i];
            b[i] = b[i + 1];
            b[i + 1] = temp - dl[i] * b[i];
        }
    }
    b[n - 1] /= d[n - 1];
    if (n > 1) b[n - 2] = (b[n - 2] - du[n - 2] * b[n - 1]) / d[n - 2];
    for (std::size_t k = n >= 3 ? n - 2 : 0; k-- > 0;) {
        b[k] = (b[k] - du[k] * b[k + 1] - du2[k] * b[k + 2]) / d[k];
    }
    return b;
}

void normalise(std::vector<double>& x) {
    double s = 0.0;
    for (const double v : x) s += v * v;
    s = std::sqrt(s);
    if (s > 0.0)
        for (double& v : x) v /= s;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

}  // namespace

std::size_t sturm_count(std::span<const double> d, std::span<const double> e, double x) {
    if (d.empty() || e.size() + 1 != d.size()) return 0;
    return sturm_count_scaled(d, e, x, eps * matrix_scale(d, e));
}

std::expected<TridiagonalEigen, std::string> tridiagonal_eigen(std::span<const double> d,
                                                               std::span<const double> e, std::size_t count) {
    const std::size_t n = d.size();
    if (n == 0) return std::unexpected("tridiagonal_eigen: empty matrix");
    if (e.size() + 1 != n) return std::unexpected("tridiagonal_eigen: off-diagonal size is not n - 1");
    if (!all_finite(d) || !all_finite(e)) return std::unexpected("tridiagonal_eigen: a non-finite entry");
    if (count == 0 || count > n) return std::unexpected("tridiagonal_eigen: count out of range");
    const double scale = matrix_scale(d, e);
    const double pivmin = eps * scale;
    double lo = d[0], hi = d[0];
    for (std::size_t i = 0; i < n; ++i) {
        const double r = (i > 0 ? std::abs(e[i - 1]) : 0.0) + (i + 1 < n ? std::abs(e[i]) : 0.0);
        lo = std::min(lo, d[i] - r);
        hi = std::max(hi, d[i] + r);
    }
    lo -= pivmin + 2.0 * eps * std::abs(lo);
    hi += pivmin + 2.0 * eps * std::abs(hi);

    TridiagonalEigen out;
    for (std::size_t k = 0; k < count; ++k) {
        double a = lo, b = hi;
        for (;;) {
            const double mid = 0.5 * (a + b);
            if (!(mid > a && mid < b) || b - a <= 2.0 * eps * std::max(std::abs(a), std::abs(b)) + pivmin) break;
            if (sturm_count_scaled(d, e, mid, pivmin) > k) {
                b = mid;
            } else {
                a = mid;
            }
        }
        const double lambda = 0.5 * (a + b);
        std::vector<double> x(n);
        for (std::size_t i = 0; i < n; ++i) x[i] = 1.0 + 0.13 * static_cast<double>(i % 7);
        for (int step = 0; step < 3; ++step) {
            x = shifted_solve(d, e, lambda, std::move(x), pivmin);
            for (int pass = 0; pass < 2; ++pass)
                for (const auto& v : out.vectors) {
                    const double c = dot(x, v);
                    for (std::size_t i = 0; i < n; ++i) x[i] -= c * v[i];
                }
            normalise(x);
        }
        double r2 = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double t = (d[i] - lambda) * x[i];
            if (i > 0) t += e[i - 1] * x[i - 1];
            if (i + 1 < n) t += e[i] * x[i + 1];
            r2 += t * t;
        }
        out.values.push_back(lambda);
        out.vectors.push_back(std::move(x));
        out.residuals.push_back(std::sqrt(r2));
    }
    return out;
}

std::expected<Subbands, std::string> schrodinger_1d(std::span<const double> z, std::span<const double> V,
                                                    std::span<const double> m, double a, std::size_t count) {
    if (z.size() < 3) return std::unexpected("schrodinger_1d: fewer than 3 nodes");
    const std::size_t cells = z.size() - 1;
    if (V.size() != cells || m.size() != cells) return std::unexpected("schrodinger_1d: per-cell sizes");
    if (!all_finite(z) || !all_finite(V) || !all_finite(m)) return std::unexpected("schrodinger_1d: non-finite input");
    for (std::size_t i = 0; i < cells; ++i) {
        if (!(z[i + 1] > z[i])) return std::unexpected("schrodinger_1d: nodes not strictly increasing");
        if (!(m[i] > 0.0)) return std::unexpected("schrodinger_1d: a mass not positive");
    }
    if (!(a > 0.0) || !std::isfinite(a)) return std::unexpected("schrodinger_1d: a not positive");
    const std::size_t n = cells - 1;  // interior nodes 1 .. cells - 1
    if (count == 0 || count > n) return std::unexpected("schrodinger_1d: count out of range");
    std::vector<double> w(n), d(n), e(n > 0 ? n - 1 : 0);
    for (std::size_t j = 0; j < n; ++j) {
        const std::size_t i = j + 1;
        const double hl = z[i] - z[i - 1], hr = z[i + 1] - z[i];
        w[j] = 0.5 * (hl + hr);
        d[j] = (a / (m[i - 1] * hl) + a / (m[i] * hr) + 0.5 * (V[i - 1] * hl + V[i] * hr)) / w[j];
    }
    for (std::size_t j = 0; j + 1 < n; ++j) {
        const std::size_t i = j + 1;
        const double hr = z[i + 1] - z[i];
        e[j] = -a / (m[i] * hr) / std::sqrt(w[j] * w[j + 1]);
    }
    auto eig = tridiagonal_eigen(d, e, count);
    if (!eig) return std::unexpected(eig.error());
    Subbands out;
    out.energies_eV = eig->values;
    out.residuals = eig->residuals;
    for (auto& y : eig->vectors) {
        double ymax = 0.0;
        for (const double v : y) ymax = std::max(ymax, std::abs(v));
        double sign = 1.0;
        for (const double v : y)
            if (std::abs(v) > 0.01 * ymax) {
                sign = v > 0.0 ? 1.0 : -1.0;
                break;
            }
        std::vector<double> psi(z.size(), 0.0);
        for (std::size_t j = 0; j < n; ++j) psi[j + 1] = sign * y[j] / std::sqrt(w[j]);
        out.psi.push_back(std::move(psi));
    }
    return out;
}

namespace {

// ln(1 + e^x) and its derivative, without overflow.
double log1p_exp(double x) { return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x)); }
double logistic(double x) { return x > 0.0 ? 1.0 / (1.0 + std::exp(-x)) : std::exp(x) / (1.0 + std::exp(x)); }

struct Grid {
    std::vector<double> z, h, eps;  // nodes; per cell: width, permittivity
    std::vector<char> silicon;      // per cell
    std::vector<double> share;      // per node: silicon share of its box
    std::size_t interface = 0, quantum_end = 0;
};

Grid make_grid(const MosSpInput& in) {
    Grid g;
    const double tox = in.oxide_thickness_cm;
    for (std::size_t i = 0; i <= in.oxide_cells; ++i)
        g.z.push_back(-tox + tox * static_cast<double>(i) / static_cast<double>(in.oxide_cells));
    g.z.back() = 0.0;
    g.interface = in.oxide_cells;
    double h = in.h0_cm, z = 0.0;
    while (z < in.depth_cm) {
        const double next = std::min(z + h, in.depth_cm);
        g.z.push_back(next);
        z = next;
        h = std::min(h * in.ratio, in.h_max_cm);
    }
    const std::size_t cells = g.z.size() - 1;
    g.h.resize(cells);
    g.eps.resize(cells);
    g.silicon.resize(cells);
    for (std::size_t c = 0; c < cells; ++c) {
        g.h[c] = g.z[c + 1] - g.z[c];
        g.silicon[c] = c >= in.oxide_cells ? 1 : 0;
        g.eps[c] = (g.silicon[c] != 0 ? in.eps_si_r : in.eps_ox_r) * in.eps0_F_per_cm;
    }
    g.share.assign(g.z.size(), 0.0);
    for (std::size_t i = 0; i < g.z.size(); ++i) {
        if (i > 0 && g.silicon[i - 1] != 0) g.share[i] += 0.5 * g.h[i - 1];
        if (i < cells && g.silicon[i] != 0) g.share[i] += 0.5 * g.h[i];
    }
    g.quantum_end = g.z.size() - 1;
    for (std::size_t i = g.interface; i < g.z.size(); ++i)
        if (g.z[i] >= in.quantum_depth_cm) {
            g.quantum_end = i;
            break;
        }
    return g;
}

// The subbands of every valley on the potential phi: the Schrodinger domain is nodes 0 .. quantum_end.
struct Levels {
    std::vector<Subbands> valley;
};

std::expected<Levels, std::string> solve_levels(const MosSpInput& in, const Grid& g, double ecb,
                                                const std::vector<double>& phi) {
    const std::size_t qn = g.quantum_end;
    std::vector<double> zq(g.z.begin(), g.z.begin() + static_cast<std::ptrdiff_t>(qn) + 1);
    std::vector<double> V(qn), m(qn);
    Levels out;
    for (const Valley& v : in.valleys) {
        for (std::size_t c = 0; c < qn; ++c) {
            const double offset = g.silicon[c] != 0 ? ecb : ecb + in.oxide_barrier_eV;
            V[c] = offset - 0.5 * (phi[c] + phi[c + 1]);
            m[c] = g.silicon[c] != 0 ? v.mass_z : in.oxide_mass;
        }
        auto s = schrodinger_1d(zq, V, m, in.a_eV_cm2, in.subbands_per_valley);
        if (!s) return std::unexpected(s.error());
        out.valley.push_back(std::move(*s));
    }
    return out;
}

struct Density {
    double value = 0.0, derivative = 0.0;  // n and dn/dphi
};

// The electron density at node i: quantum (predicted with phi - phi_ref; phi_ref = phi gives the exact quantum
// density of `levels`) inside the Schrodinger domain's silicon, classical elsewhere in the silicon.
Density electrons(const MosSpInput& in, const Grid& g, double ecb, const Levels* levels, std::size_t i,
                  double phi, double phi_ref) {
    const double kT = in.kT_eV;
    if (in.electrons == SpElectrons::classical || levels == nullptr || i < g.interface || i >= g.quantum_end) {
        const double nv = in.conduction_dos_cm3 * std::exp(-(ecb - phi) / kT);
        return {nv, nv / kT};
    }
    Density out;
    for (std::size_t v = 0; v < in.valleys.size(); ++v) {
        const Valley& val = in.valleys[v];
        const double c = val.degeneracy * val.mass_dos * kT / (2.0 * std::numbers::pi * in.a_eV_cm2);
        const Subbands& s = levels->valley[v];
        for (std::size_t k = 0; k < s.energies_eV.size(); ++k) {
            const double x = (-s.energies_eV[k] + phi - phi_ref) / kT;
            const double w = c * s.psi[k][i] * s.psi[k][i];
            if (in.statistics == SpStatistics::boltzmann) {
                const double f = std::exp(x);
                out.value += w * f;
                out.derivative += w * f / kT;
            } else {
                out.value += w * log1p_exp(x);
                out.derivative += w * logistic(x) / kT;
            }
        }
    }
    return out;
}

struct Poisson {
    std::vector<double> residual;  // per node (0 at Dirichlet nodes)
    double scale = 0.0;
};

Poisson poisson_residual(const MosSpInput& in, const Grid& g, double ecb, const Levels* levels,
                         const std::vector<double>& phi, const std::vector<double>& phi_ref) {
    const std::size_t N = g.z.size() - 1;
    Poisson out;
    out.residual.assign(N + 1, 0.0);
    for (std::size_t i = 1; i < N; ++i) {
        const double fl = g.eps[i - 1] * (phi[i] - phi[i - 1]) / g.h[i - 1];
        const double fr = g.eps[i] * (phi[i + 1] - phi[i]) / g.h[i];
        double rho = 0.0;
        if (g.share[i] > 0.0) {
            const double p = in.acceptors_cm3 * std::exp(-phi[i] / in.kT_eV);
            const double n = electrons(in, g, ecb, levels, i, phi[i], phi_ref[i]).value;
            rho = in.q_C * (p - n - in.acceptors_cm3) * g.share[i];
        }
        out.residual[i] = -(fr - fl) - rho;
        out.scale = std::max(out.scale, std::abs(fl) + std::abs(fr) + std::abs(rho));
    }
    return out;
}

// Newton on Poisson with the density of `levels` predicted from phi_ref. phi holds the Dirichlet values and the
// start; returns false when it does not converge.
bool solve_poisson(const MosSpInput& in, const Grid& g, double ecb, const Levels* levels,
                   std::vector<double>& phi, const std::vector<double>& phi_ref) {
    const std::size_t N = g.z.size() - 1;
    const std::size_t n = N - 1;
    std::vector<double> lo(n), di(n), up(n), rhs(n);
    for (int it = 0; it < 400; ++it) {
        const Poisson r = poisson_residual(in, g, ecb, levels, phi, phi_ref);
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t i = j + 1;
            const double cl = g.eps[i - 1] / g.h[i - 1], cr = g.eps[i] / g.h[i];
            double drho = 0.0;
            if (g.share[i] > 0.0) {
                const double p = in.acceptors_cm3 * std::exp(-phi[i] / in.kT_eV);
                const double dn = electrons(in, g, ecb, levels, i, phi[i], phi_ref[i]).derivative;
                drho = in.q_C * (-p / in.kT_eV - dn) * g.share[i];
            }
            lo[j] = -cl;
            up[j] = -cr;
            di[j] = cl + cr - drho;
            rhs[j] = -r.residual[i];
        }
        // Thomas (the matrix is diagonally dominant: drho <= 0).
        for (std::size_t j = 1; j < n; ++j) {
            const double f = lo[j] / di[j - 1];
            di[j] -= f * up[j - 1];
            rhs[j] -= f * rhs[j - 1];
        }
        rhs[n - 1] /= di[n - 1];
        for (std::size_t j = n - 1; j-- > 0;) rhs[j] = (rhs[j] - up[j] * rhs[j + 1]) / di[j];
        double step = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            const double dphi = std::clamp(rhs[j], -0.1, 0.1);
            phi[j + 1] += dphi;
            step = std::max(step, std::abs(dphi));
        }
        if (step <= 1e-13) return true;
    }
    return false;
}

}  // namespace

std::expected<MosSpResult, std::string> mos_schrodinger_poisson(const MosSpInput& in, double gate_V) {
    const double positive[] = {in.q_C, in.kT_eV, in.eps0_F_per_cm, in.a_eV_cm2, in.eps_si_r, in.eps_ox_r,
                               in.band_gap_eV, in.conduction_dos_cm3, in.valence_dos_cm3, in.acceptors_cm3,
                               in.oxide_thickness_cm, in.h0_cm, in.ratio, in.h_max_cm, in.depth_cm,
                               in.outer_tolerance_V};
    for (const double x : positive)
        if (!(x > 0.0) || !std::isfinite(x)) return std::unexpected("mos_schrodinger_poisson: an input not positive");
    if (!std::isfinite(gate_V)) return std::unexpected("mos_schrodinger_poisson: gate bias not finite");
    if (in.oxide_cells == 0) return std::unexpected("mos_schrodinger_poisson: no oxide cells");
    if (in.ratio < 1.0 || in.h_max_cm < in.h0_cm) return std::unexpected("mos_schrodinger_poisson: grading");
    if (!(in.acceptors_cm3 < in.valence_dos_cm3))  // Unit 62 slice 4: the Boltzmann bulk, E_F in the gap
        return std::unexpected("mos_schrodinger_poisson: N_A not below N_v (Boltzmann bulk)");
    if (in.electrons == SpElectrons::quantum) {
        if (in.valleys.empty() || in.subbands_per_valley == 0)
            return std::unexpected("mos_schrodinger_poisson: no valleys or subbands");
        for (const Valley& v : in.valleys)
            if (!(v.mass_z > 0.0) || !(v.mass_dos > 0.0) || !(v.degeneracy > 0.0))
                return std::unexpected("mos_schrodinger_poisson: a valley parameter not positive");
        if (!(in.oxide_mass > 0.0) || !(in.oxide_barrier_eV > 0.0))
            return std::unexpected("mos_schrodinger_poisson: oxide mass or barrier not positive");
        if (!(in.quantum_depth_cm > 0.0) || in.quantum_depth_cm >= in.depth_cm)
            return std::unexpected("mos_schrodinger_poisson: quantum depth outside the silicon");
    }

    const Grid g = make_grid(in);
    const std::size_t N = g.z.size() - 1;
    const double ecb = in.band_gap_eV - in.kT_eV * std::log(in.valence_dos_cm3 / in.acceptors_cm3);

    // Start: the classical solution, from the oxide's linear drop and phi = 0 in the silicon.
    std::vector<double> phi(N + 1, 0.0);
    for (std::size_t i = 0; i <= g.interface; ++i)
        phi[i] = gate_V * static_cast<double>(g.interface - i) / static_cast<double>(g.interface);
    MosSpInput classical = in;
    classical.electrons = SpElectrons::classical;
    if (!solve_poisson(classical, g, ecb, nullptr, phi, phi))
        return std::unexpected("mos_schrodinger_poisson: classical Newton did not converge");

    MosSpResult out;
    Levels levels;
    if (in.electrons == SpElectrons::quantum) {
        int it = 0;
        for (;; ++it) {
            if (it >= in.max_outer) return std::unexpected("mos_schrodinger_poisson: outer loop did not converge");
            auto lv = solve_levels(in, g, ecb, phi);
            if (!lv) return std::unexpected(lv.error());
            levels = std::move(*lv);
            const std::vector<double> phi_old = phi;
            if (!solve_poisson(in, g, ecb, &levels, phi, phi_old))
                return std::unexpected("mos_schrodinger_poisson: predictor Newton did not converge");
            double update = 0.0;
            for (std::size_t i = 0; i <= N; ++i) update = std::max(update, std::abs(phi[i] - phi_old[i]));
            out.last_update_V = update;
            if (update <= in.outer_tolerance_V) break;
        }
        out.outer_iterations = it + 1;
        auto lv = solve_levels(in, g, ecb, phi);  // the final phi's own subbands
        if (!lv) return std::unexpected(lv.error());
        levels = std::move(*lv);
        for (const Subbands& s : levels.valley) {
            out.subband_eV.push_back(s.energies_eV);
            for (const double r : s.residuals) out.max_eigen_residual_eV = std::max(out.max_eigen_residual_eV, r);
            for (std::size_t j = 0; j < s.psi.size(); ++j)
                for (std::size_t k = 0; k <= j; ++k) {
                    double ip = 0.0;
                    for (std::size_t i = 1; i < g.quantum_end; ++i) {
                        const double w = 0.5 * (g.z[i + 1] - g.z[i - 1]);
                        ip += w * s.psi[j][i] * s.psi[k][i];
                    }
                    out.orthonormality = std::max(out.orthonormality, std::abs(ip - (j == k ? 1.0 : 0.0)));
                }
        }
    }
    const Levels* lp = in.electrons == SpElectrons::quantum ? &levels : nullptr;
    const Poisson r = poisson_residual(in, g, ecb, lp, phi, phi);
    double rmax = 0.0;
    for (const double x : r.residual) rmax = std::max(rmax, std::abs(x));
    out.poisson_residual = rmax / r.scale;

    out.z_cm = g.z;
    out.phi_V = phi;
    out.interface_node = g.interface;
    out.n_cm3.assign(N + 1, 0.0);
    out.p_cm3.assign(N + 1, 0.0);
    double ninv = 0.0, zsum = 0.0, qs = 0.0;
    for (std::size_t i = g.interface; i <= N; ++i) {
        const double nn = electrons(in, g, ecb, lp, i, phi[i], phi[i]).value;
        const double pp = in.acceptors_cm3 * std::exp(-phi[i] / in.kT_eV);
        out.n_cm3[i] = nn;
        out.p_cm3[i] = pp;
        if (i == N) continue;  // the Dirichlet node's box is outside the solved equations
        ninv += nn * g.share[i];
        zsum += g.z[i] * nn * g.share[i];
        qs += in.q_C * (pp - nn - in.acceptors_cm3) * g.share[i];
    }
    out.inversion_cm2 = ninv;
    out.centroid_cm = ninv > 0.0 ? zsum / ninv : 0.0;
    out.charge_C_per_cm2 = qs;
    out.oxide_displacement_C_per_cm2 =
        in.eps_ox_r * in.eps0_F_per_cm * (phi[0] - phi[g.interface]) / in.oxide_thickness_cm;
    out.surface_potential_V = phi[g.interface];
    return out;
}

}  // namespace NiTCAD::verify::reference
