// Fermi-Dirac integrals in double-double (Unit 55, TVR-1, TVB-02); the derivation is in the header.
#include "NiTCAD/verify/reference/fermi_dirac.hpp"

#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

#include <cmath>
#include <stdexcept>

namespace NiTCAD::verify::reference {

namespace {

constexpr double cut_L = 100.0;
constexpr double piece_epsrel = 1e-27;

int power(FermiOrder j) {  // p = 2j + 1
    switch (j) {
        case FermiOrder::minus_half: return 0;
        case FermiOrder::half: return 2;
        case FermiOrder::three_halves: return 4;
    }
    return 0;
}

void check_domain(const DD& x) {
    if (!(x.hi > -cut_L && x.hi <= 600.0)) throw std::invalid_argument("fermi_dirac reference: x outside (-100, 600]");
}

// 2 u^p / (1 + e^(u^2 - x)).
struct Integrand {
    int p;
    DD x;
    DD operator()(const DD& u) const {
        const DD y = u * u - x;
        const DD logistic = y.hi > 0.0 ? exp(-y) / (DD(1.0) + exp(-y)) : DD(1.0) / (DD(1.0) + exp(y));
        DD up(1.0);
        for (int i = 0; i < p; ++i) up = up * u;
        return ldexp(up * logistic, 1);
    }
};

// The pieces [0, sqrt(x)] (x > 0) and [max(0, sqrt(x)), U], U = sqrt(x + L).
struct Pieces {
    DD bounds[3];
    int count = 0;  // number of pieces
};
Pieces pieces(const DD& x) {
    Pieces p;
    const DD U = sqrt(x + DD(cut_L));
    if (x.hi > 0.0) {
        p.bounds[0] = DD(0.0);
        p.bounds[1] = sqrt(x);
        p.bounds[2] = U;
        p.count = 2;
    } else {
        p.bounds[0] = DD(0.0);
        p.bounds[1] = U;
        p.count = 1;
    }
    return p;
}

double tail(FermiOrder j, const DD& x) {
    const int p = power(j);
    const double U2 = x.hi + cut_L;
    const double U = std::sqrt(U2);
    // p = 0: I(0) = e^(-U^2) / (2U) - (1/2) I(-2) <= e^(-U^2) / (2U); p >= 1: the recurrence bound (header).
    const double denominator = p == 0 ? 2.0 : 2.0 - (p - 1) / U2;
    const double bound = std::pow(U, p - 1) * std::exp(-cut_L) / denominator;
    return 2.0 * bound * fermi_normalisation(j).to_double();
}

template <class Rule>
FermiQuadrature integrate(FermiOrder j, const DD& x, Rule&& rule) {
    check_domain(x);
    const Integrand f{power(j), x};
    const Pieces ps = pieces(x);
    FermiQuadrature out;
    out.ok = true;
    DD sum(0.0);
    double err = 0.0;
    for (int i = 0; i < ps.count; ++i) {
        const auto r = rule(f, ps.bounds[i], ps.bounds[i + 1]);
        sum += r.value;
        err += r.abs_error;
        out.evaluations += r.evaluations;
        out.ok = out.ok && r.ok;
    }
    const DD n = fermi_normalisation(j);
    out.value = sum * n;
    out.abs_error = err * n.to_double();
    out.tail_bound = tail(j, x);
    return out;
}

struct PieceResult {
    DD value;
    double abs_error;
    int evaluations;
    bool ok;
};

}  // namespace

DD fermi_normalisation(FermiOrder j) {
    const DD sp = sqrt(pi());
    switch (j) {
        case FermiOrder::minus_half: return DD(1.0) / sp;
        case FermiOrder::half: return DD(2.0) / sp;
        case FermiOrder::three_halves: return DD(4.0) / (DD(3.0) * sp);
    }
    return DD(0.0);
}

DD fermi_eta_factor(FermiOrder j) {
    switch (j) {
        case FermiOrder::minus_half: return DD(1.0) - sqrt2();                 // 1 - 2^(1/2)
        case FermiOrder::half: return DD(1.0) - DD(1.0) / sqrt2();             // 1 - 2^(-1/2)
        case FermiOrder::three_halves: return DD(1.0) - DD(0.5) / sqrt2();     // 1 - 2^(-3/2)
    }
    return DD(0.0);
}

FermiQuadrature fermi_tanh_sinh(FermiOrder j, const DD& x) {
    return integrate(j, x, [](const Integrand& f, const DD& a, const DD& b) {
        const TanhSinhResult r = tanh_sinh(f, a, b, 0.0, piece_epsrel);
        return PieceResult{r.value, r.abs_error, r.evaluations, r.status == TanhSinhStatus::ok};
    });
}

FermiQuadrature fermi_gauss_kronrod(FermiOrder j, const DD& x) {
    return integrate(j, x, [](const Integrand& f, const DD& a, const DD& b) {
        const AdaptiveResult r = qag15(f, a, b, 0.0, piece_epsrel);
        return PieceResult{r.value, r.abs_error, r.evaluations, r.status == QuadratureStatus::ok};
    });
}

FermiSeries fermi_series(FermiOrder j, const DD& x) {
    if (!(x.hi <= -2.0)) throw std::invalid_argument("fermi_dirac reference: the series is used for x <= -2 only");
    // k^(j+1) = k^(p/2 + 1/2): sqrt(k) times k^(p/2).
    const int half_p = power(j) / 2;
    const DD z = exp(x);
    FermiSeries out;
    DD zk = z;
    DD sum(0.0);
    for (int k = 1;; ++k) {
        DD kp = sqrt(DD(static_cast<double>(k)));
        for (int i = 0; i < half_p; ++i) kp = kp * DD(static_cast<double>(k));
        const DD term = zk / kp;
        sum = (k % 2 == 1) ? sum + term : sum - term;
        out.terms = k;
        zk = zk * z;
        if (zk.hi < 0x1p-120 * z.hi) {
            // Alternating with decreasing terms: the remainder is below the first omitted term.
            out.truncation_bound = (zk / sqrt(DD(static_cast<double>(k + 1)))).to_double();
            break;
        }
    }
    out.value = sum;
    return out;
}

DD fermi_integrand(FermiOrder j, const DD& x, const DD& u) {
    return Integrand{power(j), x}(u) * fermi_normalisation(j);
}

DD fermi_cut(const DD& x) { return sqrt(x + DD(cut_L)); }

}  // namespace NiTCAD::verify::reference
