#include "NiTCAD/verify/reference/tunnel_generation.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/model_generation.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

namespace NiTCAD::verify::reference {

namespace {

namespace mg = model_generation;

constexpr double outer_rtol = 1e-13;   // the pairs: 1e-13 relative, far below any comparison's tolerance
constexpr double inner_rtol = 1e-24;   // a path's integrals
constexpr double kane_cut = 700.0;     // G = 0 where B / F > 700, as physics::kane_generation

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

// [a, b] cut at the levels of `at` that lie strictly inside, ascending.
std::vector<DD> pieces(const DD& a, const DD& b, std::initializer_list<DD> at) {
    std::vector<DD> mid;
    for (const DD& c : at) {
        if (c > a && c < b) mid.push_back(c);
    }
    std::sort(mid.begin(), mid.end(), [](const DD& x, const DD& y) { return x < y; });
    std::vector<DD> out{a};
    for (const DD& c : mid) {
        if (!(c == out.back())) out.push_back(c);
    }
    out.push_back(b);
    return out;
}

struct Sum {
    DD value = DD(0.0);
    double abs_error = 0.0;
    bool ok = true;
};

template <class F>
Sum kronrod(F&& f, const std::vector<DD>& at, double rtol, double abs_floor = 0.0) {
    Sum s;
    for (std::size_t k = 0; k + 1 < at.size(); ++k) {
        const AdaptiveResult a = qag15(f, at[k], at[k + 1], abs_floor, rtol);
        // qag's tests are comparisons, which a NaN fails silently: a non-finite piece is not converged (a NaN I_ik
        // passed as converged at TVB-24's rerun, 17.27).
        const bool finite = std::isfinite(a.value.hi) && std::isfinite(a.abs_error);
        s.ok = s.ok && finite && a.status == QuadratureStatus::ok;
        s.value += a.value;
        s.abs_error += a.abs_error;
    }
    return s;
}

// Every caller's integrand is normalised so that the whole integral is of order one; rtol is then also the absolute
// floor of each piece, so that a piece far smaller than the whole (the slope term's i region, where E' vanishes but
// for the carriers) is not held to rtol of itself.
template <class F>
Sum tanh_sinh_sum(F&& f, const std::vector<DD>& at, double rtol) {
    Sum s;
    for (std::size_t k = 0; k + 1 < at.size(); ++k) {
        const TanhSinhResult a = tanh_sinh(f, at[k], at[k + 1], rtol, rtol);
        s.ok = s.ok && a.status == TanhSinhStatus::ok;
        s.value += a.value;
        s.abs_error += a.abs_error;
    }
    return s;
}

DD kane_rate(const KaneInput& k, const DD& F) {
    if (!(F.hi * kane_cut > k.B_V_per_cm)) return DD(0.0);
    return mg::kane<DD>({DD(k.A_per_cm3_s), DD(k.B_V_per_cm), F})[0];
}

bool valid(const KaneInput& k) { return positive(k.A_per_cm3_s) && positive(k.B_V_per_cm); }

// The second method's length: V_T int ds / E over [a, b] by tanh-sinh, split at s_0 and s_d, in units of the first
// method's length `unit` (tanh-sinh's estimate is absolute, so the integrand is kept of order one). In the i region
// in s; on a side in the logarithm of the distance c to that side's bulk level (ds = c d ln c), since a path may end
// e^-30 or less from a bulk level, where E vanishes with c. The first method differs in rule (Kronrod) and integrand
// (1 / sqrt(scale Phi / c^2), from the first integral rather than from E).
std::expected<DD, std::string> length_second(const PinJunction& j, const DD& a, const DD& b, const DD& unit) {
    const std::vector<DD> at = pieces(a, b, {j.s0, j.sd});
    // In ln c the integrand is flat (the bulk's Debye length) for c << 1 and turns over near c = 1: cut there.
    const auto in_logs = [&](const DD& t_lo, const DD& t_hi, const auto& g) {
        return tanh_sinh_sum(g, pieces(t_lo, t_hi, {DD(-8.0), DD(-4.0), DD(-2.0), DD(-1.0), DD(0.0), DD(1.0), DD(2.0)}),
                             inner_rtol);
    };
    Sum s;
    for (std::size_t k = 0; k + 1 < at.size(); ++k) {
        const DD lo = at[k], hi = at[k + 1];
        Sum r;
        if (!(hi > j.s0)) {  // p side: c = s
            r = in_logs(log(lo), log(hi),
                        [&](const DD& t) { const DD c = exp(t); return j.VT * c / (pin_field(j, c) * unit); });
        } else if (!(lo < j.sd)) {  // n side: c = B - s, descending in s
            r = in_logs(log(j.barrier - hi), log(j.barrier - lo), [&](const DD& t) {
                const DD c = exp(t);
                return j.VT * c / (pin_field_n_bulk(j, c) * unit);
            });
        } else {
            r = tanh_sinh_sum([&](const DD& z) { return j.VT / (pin_field(j, z) * unit); }, {lo, hi}, inner_rtol);
        }
        s.ok = s.ok && r.ok;
        s.value += r.value;
    }
    if (!s.ok) return std::unexpected("tunnel_generation: a length's tanh-sinh quadrature did not converge");
    return s.value * unit;
}

std::expected<DD, std::string> length_of(const PinJunction& j, const DD& a, const DD& b) {
    const auto xa = pin_position(j, a);
    if (!xa) return std::unexpected(xa.error());
    const auto xb = pin_position(j, b);
    if (!xb) return std::unexpected(xb.error());
    return *xb - *xa;
}

// The WKB band and windows.
struct Wkb {
    mg::Band band;
    DD kv, kc;  // 2 m q V_T / hbar^2 [1/m^2 per level]
    DD g;       // the gap in levels
    WkbLevels levels;
};

std::expected<Wkb, std::string> wkb_setup(const PinJunction& j, const WkbInput& w) {
    if (!positive(w.gap_eV) || !positive(w.electron_mass) || !positive(w.hole_mass) || !positive(w.flat_field_V_per_cm)) {
        return std::unexpected("tunnel_generation: the WKB gap, masses and flat field must be finite and positive");
    }
    Wkb b;
    const DD m0(base::m0_kg), q(base::q_C);
    const DD mc = DD(w.electron_mass) * m0, mv = DD(w.hole_mass) * m0;
    b.band = {DD(w.gap_eV) * q, mc * mv / (mc + mv)};
    if (!(mg::u_of(b.band).hi > 1.0)) return std::unexpected("tunnel_generation: u = m0 / (2 m_r) must exceed 1");
    const DD h2 = mg::hbar() * mg::hbar();
    b.kv = DD(2.0) * mv * q * j.VT / h2;
    b.kc = DD(2.0) * mc * q * j.VT / h2;
    b.g = DD(w.gap_eV) / j.VT;
    const auto lv = wkb_levels(j, w);
    if (!lv) return std::unexpected(lv.error());
    b.levels = *lv;
    return b;
}

// kappa(delta) [1/m] with delta in [0, 1] (model_generation clips).
DD kappa_at(const Wkb& b, const DD& delta) { return mg::kappa(delta, b.band); }

struct WkbIntegrals {
    DD I_k, I_ik;  // [-], [m^2]
    bool ok = true;
};

// The method: tanh-sinh in delta over [0, 1] (s = s_a + g delta), split at the region edges.
WkbIntegrals wkb_integrals(const PinJunction& j, const Wkb& b, const DD& s_a) {
    const DD m_per_cm(1e-2);
    const auto E = [&](const DD& delta) { return pin_field(j, s_a + b.g * delta); };
    const std::vector<DD> at = pieces(DD(0.0), DD(1.0), {(j.s0 - s_a) / b.g, (j.sd - s_a) / b.g});
    // Normalized to order one: I_ik by kappa_0 / L_0, kappa_0 = sqrt(m_r E_g) / hbar, L_0 = g V_T / E(s_a).
    const DD kappa0 = mg::scale_of(b.band) / mg::hbar();
    const DD L0 = b.g * j.VT / pin_field(j, s_a) * m_per_cm;
    const Sum k = tanh_sinh_sum([&](const DD& dl) { return kappa_at(b, dl) * b.g * j.VT / E(dl) * m_per_cm; }, at,
                                inner_rtol);
    const Sum i = tanh_sinh_sum(
        [&](const DD& dl) { return b.g * j.VT / E(dl) * m_per_cm / kappa_at(b, dl) * kappa0 / L0; }, at, inner_rtol);
    return {k.value, i.value * L0 / kappa0, k.ok && i.ok};
}

// The second method: Kronrod-15 in t with delta = t^2 (3 - 2t), which removes 1 / kappa's band-edge singularities.
WkbIntegrals wkb_integrals_second(const PinJunction& j, const Wkb& b, const DD& s_a) {
    const DD m_per_cm(1e-2);
    const auto delta_of = [](const DD& t) { return t * t * (DD(3.0) - ldexp(t, 1)); };
    const auto t_of = [&](const DD& delta) {
        const BrentResult r = brent_zero([&](const DD& t) { return delta_of(t) - delta; }, DD(0.0), DD(1.0), 1e-30);
        return r.root;
    };
    std::vector<DD> at{DD(0.0)};
    for (const DD& c : {(j.s0 - s_a) / b.g, (j.sd - s_a) / b.g}) {
        if (c.hi > 0.0 && c.hi < 1.0) at.push_back(t_of(c));
    }
    std::sort(at.begin(), at.end(), [](const DD& x, const DD& y) { return x < y; });
    at.push_back(DD(1.0));
    const auto jac = [](const DD& t) { return DD(6.0) * t * (DD(1.0) - t); };  // d delta / dt
    const auto E = [&](const DD& t) { return pin_field(j, s_a + b.g * delta_of(t)); };
    // kappa from t with 1 - delta = (1 - t)^2 (1 + 2t), exact where delta rounds near 1 (1 - delta formed from delta
    // there is 0 within u_dd, and 1 / kappa infinite: TVB-24's rerun at -0.5 V, 17.27).
    const auto kappa_t = [&](const DD& t) {
        const DD r = DD(1.0) - t;
        return mg::kappa(delta_of(t), r * r * (DD(1.0) + ldexp(t, 1)), b.band);
    };
    const auto fk = [&](const DD& t) { return kappa_t(t) * b.g * j.VT / E(t) * m_per_cm * jac(t); };
    const auto fi = [&](const DD& t) { return b.g * j.VT / E(t) * m_per_cm / kappa_t(t) * jac(t); };
    // The tolerance is relative to the whole integral, not to each piece: where a region edge falls within ~1e-8 of
    // t = 1 its sliver carries ~1e-6 of the total, and a relative tolerance on the sliver alone would drive qag to
    // nodes where delta rounds to 1 and 1 / kappa is infinite.
    const auto floor_of = [&](const auto& f) {
        double whole = 0.0;
        for (std::size_t q = 0; q + 1 < at.size(); ++q) whole += qk15(f, at[q], at[q + 1]).resabs;
        return inner_rtol * whole;
    };
    const Sum k = kronrod(fk, at, inner_rtol, floor_of(fk));
    const Sum i = kronrod(fi, at, inner_rtol, floor_of(fi));
    return {k.value, i.value, k.ok && i.ok};
}

struct WkbRate {
    DD rate_cm3_s, slope, km2;
};

WkbRate wkb_rate(const PinJunction& j, const Wkb& b, const DD& s_a, const WkbIntegrals& in) {
    const DD q(base::q_C);
    const DD s_b = s_a + b.g;
    const DD slope = q * pin_field(j, s_a) * DD(100.0);  // |dE_v/dx| [J/m]
    const DD s_v = s_a < b.levels.s_lo ? s_a : b.levels.s_lo;
    const DD s_c = s_b > b.levels.s_hi ? s_b : b.levels.s_hi;
    const DD win_v = b.kv * (s_a - s_v), win_c = b.kc * (s_c - s_b);  // both >= 0 by the choice of s_v and s_c
    const DD km2 = win_v < win_c ? win_v : win_c;
    const DD G = mg::path_rate({slope, in.I_k, in.I_ik, km2})[0] * DD(1e-6);
    return {G, slope, km2};
}

}  // namespace

std::expected<TunnelPairs, std::string> local_kane_pairs(const PinJunction& j, const KaneInput& k) {
    if (!valid(k)) return std::unexpected("local_kane_pairs: A and B must be finite and positive");
    const auto f = [&](const DD& s) {
        const DD E = pin_field(j, s);
        const DD G = kane_rate(k, E);  // 0 where E underflows deep in a bulk: no 0 / 0
        return G.hi == 0.0 ? G : j.VT * G / E;  // G dx / ds
    };
    const std::vector<DD> at = pieces(DD(0.0), j.barrier, {j.s0, j.sd});
    const Sum a = kronrod(f, at, outer_rtol);
    if (!a.ok) return std::unexpected("local_kane_pairs: the quadrature did not converge");
    if (!(a.value.hi > 0.0)) return TunnelPairs{a.value, a.value, a.abs_error};
    const Sum b = tanh_sinh_sum([&](const DD& s) { return f(s) / a.value; }, at, 1e-20);
    if (!b.ok) return std::unexpected("local_kane_pairs: the second method did not converge");
    return TunnelPairs{a.value, b.value * a.value, a.abs_error};
}

std::expected<PathReference, std::string> kane_path(const PinJunction& j, double gap_eV, const KaneInput& k,
                                                    const DD& s_start) {
    if (!valid(k) || !positive(gap_eV)) return std::unexpected("kane_path: A, B and the gap must be finite and positive");
    if (!(s_start.hi > 0.0 && s_start < j.barrier)) return std::unexpected("kane_path: the start is outside (0, B)");
    PathReference p;
    const DD s_b = s_start + DD(gap_eV) / j.VT;
    if (!(s_b < j.barrier)) return p;
    const auto l = length_of(j, s_start, s_b);
    if (!l) return std::unexpected(l.error());
    const auto l2 = length_second(j, s_start, s_b, *l);
    if (!l2) return std::unexpected(l2.error());
    p.exists = true;
    p.length_cm = *l;
    p.field_V_per_cm = DD(gap_eV) / *l;
    p.rate_cm3_s = kane_rate(k, p.field_V_per_cm);
    p.second_length_cm = *l2;
    p.second_rate_cm3_s = kane_rate(k, DD(gap_eV) / *l2);
    return p;
}

std::expected<TunnelPairs, std::string> nonlocal_kane_pairs(const PinJunction& j, double gap_eV, const KaneInput& k) {
    if (!valid(k) || !positive(gap_eV)) {
        return std::unexpected("nonlocal_kane_pairs: A, B and the gap must be finite and positive");
    }
    const DD g = DD(gap_eV) / j.VT;
    const DD top = j.barrier - g;
    if (!(top.hi > 0.0)) return TunnelPairs{};
    std::string error;
    const auto f = [&](const DD& s_a, bool second) {
        const DD s_b = s_a + g;
        if (!(s_a.hi > 0.0) || !(s_b < j.barrier)) return DD(0.0);
        auto l = length_of(j, s_a, s_b);
        // Where the first method's mean field is below the cut (B / F > 700, deep in a bulk) G is 0 by definition;
        // the second length is not formed there, its integrand being exponentially large.
        if (l && second && !(kane_rate(k, DD(gap_eV) / *l).hi == 0.0)) l = length_second(j, s_a, s_b, *l);
        if (!l) {
            error = l.error();
            return DD(0.0);
        }
        const DD G = kane_rate(k, DD(gap_eV) / *l);
        return G.hi == 0.0 ? G : j.VT * G / pin_field(j, s_a);
    };
    const std::vector<DD> at = pieces(DD(0.0), top, {j.s0, j.sd, j.s0 - g, j.sd - g});
    const Sum a = kronrod([&](const DD& s) { return f(s, false); }, at, outer_rtol);
    if (!error.empty()) return std::unexpected(error);
    if (!a.ok) return std::unexpected("nonlocal_kane_pairs: the quadrature did not converge");
    if (!(a.value.hi > 0.0)) return TunnelPairs{a.value, a.value, a.abs_error};
    const Sum b = tanh_sinh_sum([&](const DD& s) { return f(s, true) / a.value; }, at, 1e-20);
    if (!error.empty()) return std::unexpected(error);
    if (!b.ok) return std::unexpected("nonlocal_kane_pairs: the second method did not converge");
    return TunnelPairs{a.value, b.value * a.value, a.abs_error};
}

std::expected<WkbLevels, std::string> wkb_levels(const PinJunction& j, const WkbInput& w) {
    if (!positive(w.flat_field_V_per_cm)) return std::unexpected("wkb_levels: the flat field must be positive");
    const DD flat(w.flat_field_V_per_cm);
    if (!(pin_field(j, j.s0) > flat) || !(pin_field(j, j.sd) > flat)) {
        return std::unexpected("wkb_levels: the field at an edge of the i region is below the flat field");
    }
    WkbLevels out;
    // p side: E rises with s on (0, s_0]; n side: E rises as r = B - s grows on (0, B - s_d].
    const DD tiny_p = j.s0 * DD(1e-30), rd = j.barrier - j.sd, tiny_n = rd * DD(1e-30);
    const BrentResult p = brent_zero([&](const DD& s) { return pin_field(j, s) - flat; }, tiny_p, j.s0, 4.0 * u_dd * j.s0.hi);
    const BrentResult n = brent_zero([&](const DD& r) { return pin_field(j, j.barrier - r) - flat; }, tiny_n, rd,
                                     4.0 * u_dd * rd.hi);
    if (p.status != BrentStatus::ok || n.status != BrentStatus::ok) {
        return std::unexpected("wkb_levels: Brent did not converge for a flat-field level");
    }
    out.s_lo = p.root;
    out.s_hi = j.barrier - n.root;
    return out;
}

std::expected<PathReference, std::string> wkb_path(const PinJunction& j, const WkbInput& w, const DD& s_start) {
    const auto b = wkb_setup(j, w);
    if (!b) return std::unexpected(b.error());
    if (!(s_start.hi > 0.0 && s_start < j.barrier)) return std::unexpected("wkb_path: the start is outside (0, B)");
    PathReference p;
    const DD s_b = s_start + b->g;
    if (!(s_b < j.barrier)) return p;
    const auto l = length_of(j, s_start, s_b);
    if (!l) return std::unexpected(l.error());
    const auto l2 = length_second(j, s_start, s_b, *l);
    if (!l2) return std::unexpected(l2.error());
    const WkbIntegrals a = wkb_integrals(j, *b, s_start);
    const WkbIntegrals c = wkb_integrals_second(j, *b, s_start);
    if (!a.ok || !c.ok) return std::unexpected("wkb_path: a path integral did not converge");
    const WkbRate r = wkb_rate(j, *b, s_start, a);
    p.exists = true;
    p.length_cm = *l;
    p.field_V_per_cm = DD(w.gap_eV) / *l;
    p.rate_cm3_s = r.rate_cm3_s;
    p.second_length_cm = *l2;
    p.second_rate_cm3_s = wkb_rate(j, *b, s_start, c).rate_cm3_s;
    p.I_k = a.I_k;
    p.I_ik_m2 = a.I_ik;
    p.km2_per_m2 = r.km2;
    p.slope_J_per_m = r.slope;
    return p;
}

namespace {

// The WKB pairs' integrand family over the starts: G V_T / E (weight 1) or G (E'/E) V_T / E / 2 (the slope term).
template <class W>
std::expected<TunnelPairs, std::string> wkb_integral(const PinJunction& j, const WkbInput& w, W&& weight) {
    const auto b = wkb_setup(j, w);
    if (!b) return std::unexpected(b.error());
    const DD lo = b->levels.s_lo, hi = b->levels.s_hi - b->g;
    if (!(hi > lo)) return TunnelPairs{};
    bool failed = false;
    const auto f = [&](const DD& s_a, bool second) {
        const WkbIntegrals in = second ? wkb_integrals_second(j, *b, s_a) : wkb_integrals(j, *b, s_a);
        if (!in.ok) failed = true;
        const DD G = wkb_rate(j, *b, s_a, in).rate_cm3_s;
        return G.hi == 0.0 ? G : G * j.VT / pin_field(j, s_a) * weight(s_a);
    };
    const DD mv(w.hole_mass), mc(w.electron_mass);
    const DD meet = (mv * b->levels.s_lo + mc * (b->levels.s_hi - b->g)) / (mv + mc);  // win_v = win_c
    const std::vector<DD> at = pieces(lo, hi, {j.s0, j.sd, j.s0 - b->g, j.sd - b->g, meet});
    const Sum a = kronrod([&](const DD& s) { return f(s, false); }, at, outer_rtol);
    if (failed || !a.ok) return std::unexpected("wkb_pairs: a quadrature did not converge");
    if (a.value.hi == 0.0) return TunnelPairs{a.value, a.value, a.abs_error};
    const Sum c = tanh_sinh_sum([&](const DD& s) { return f(s, true) / a.value; }, at, 1e-20);
    if (failed || !c.ok) return std::unexpected("wkb_pairs: the second method did not converge");
    return TunnelPairs{a.value, c.value * a.value, a.abs_error};
}

}  // namespace

std::expected<TunnelPairs, std::string> wkb_pairs(const PinJunction& j, const WkbInput& w) {
    return wkb_integral(j, w, [](const DD&) { return DD(1.0); });
}

std::expected<DD, std::string> wkb_slope_coefficient(const PinJunction& j, const WkbInput& w) {
    const auto P = wkb_pairs(j, w);
    if (!P) return std::unexpected(P.error());
    if (!(P->pairs_cm2_s.hi > 0.0)) return std::unexpected("wkb_slope_coefficient: no pairs");
    const auto S = wkb_integral(j, w, [&](const DD& s) { return pin_field_slope(j, s) / (ldexp(pin_field(j, s), 1)); });
    if (!S) return std::unexpected(S.error());
    return S->pairs_cm2_s / P->pairs_cm2_s;
}

}  // namespace NiTCAD::verify::reference
