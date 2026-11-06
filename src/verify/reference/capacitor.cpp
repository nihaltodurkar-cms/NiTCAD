#include "NiTCAD/verify/reference/capacitor.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "NiTCAD/verify/reference/codata2018.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

}  // namespace

DD vacuum_permittivity_F_cm() {
    // CODATA 2018 prints eps0 = 8.854 187 8128 e-12 F/m; one more decade per 1/100 m.
    static const DD value = codata2018::to_dd(codata2018::vacuum_permittivity, -2);
    return value;
}

DD series_capacitance(std::span<const double> thickness_cm, std::span<const double> eps_r) {
    if (thickness_cm.empty() || thickness_cm.size() != eps_r.size()) {
        throw std::invalid_argument("series_capacitance: one permittivity per layer, at least one layer");
    }
    DD inverse{0.0, 0.0};  // sum t_i / eps_r,i [cm]
    for (std::size_t i = 0; i < thickness_cm.size(); ++i) {
        if (!positive(thickness_cm[i]) || !positive(eps_r[i])) {
            throw std::invalid_argument("series_capacitance: thickness and permittivity finite and positive");
        }
        inverse += DD{thickness_cm[i], 0.0} / DD{eps_r[i], 0.0};
    }
    return vacuum_permittivity_F_cm() / inverse;
}

DD coaxial_capacitance(double inner_cm, double outer_cm, double eps_r) {
    if (!positive(inner_cm) || !positive(outer_cm) || !(inner_cm < outer_cm) || !positive(eps_r)) {
        throw std::invalid_argument("coaxial_capacitance: 0 < a < b and eps_r > 0, all finite");
    }
    const DD ratio = DD{outer_cm, 0.0} / DD{inner_cm, 0.0};
    return ldexp(pi(), 1) * vacuum_permittivity_F_cm() * DD{eps_r, 0.0} / log(ratio);
}

}  // namespace NiTCAD::verify::reference
