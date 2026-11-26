#include "NiTCAD/physics/polarization.hpp"

namespace NiTCAD::physics {
namespace {

constexpr NitrideLattice gan{3.189, -0.034, -0.49, 0.73, 106.0, 398.0};
constexpr NitrideLattice aln{3.112, -0.090, -0.60, 1.46, 108.0, 373.0};
constexpr NitrideLattice inn{3.533, -0.042, -0.57, 0.97, 92.0, 224.0};

NitrideLattice mix(const NitrideLattice& b, double x, double bowing) {
    const auto lin = [x](double u, double v) { return (1.0 - x) * u + x * v; };
    return {lin(gan.a0_A, b.a0_A),          lin(gan.P_sp, b.P_sp) - bowing * x * (1.0 - x),
            lin(gan.e31, b.e31),            lin(gan.e33, b.e33),
            lin(gan.c13, b.c13),            lin(gan.c33, b.c33)};
}

}  // namespace

std::optional<NitrideLattice> nitride_lattice(std::string_view material, double x) {
    if (material == "gan") return gan;
    if (material == "aln") return aln;
    if (material == "inn") return inn;
    if (material == "algan") return mix(aln, x, -0.021);
    if (material == "ingan") return mix(inn, x, -0.037);
    return std::nullopt;
}

double nitride_polarization(const NitrideLattice& l, double a_A) {
    return l.P_sp + 2.0 * (a_A - l.a0_A) / l.a0_A * (l.e31 - l.e33 * l.c13 / l.c33);
}

double polarization_sheet_density_cm2(double P_C_m2) { return P_C_m2 / 1.602176634e-19 * 1e-4; }

}  // namespace NiTCAD::physics
