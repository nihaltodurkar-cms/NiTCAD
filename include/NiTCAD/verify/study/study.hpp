// Convergence studies against exact solutions (Unit 55, TVR-2, DECISIONS Unit 55 14.5-14.7: CV-Q6-Q8).
//
// A study is a ladder of mesh levels, finest first, each with its spacing h, the error of one quantity or norm
// against the exact solution, and the iterative-error estimate U_iter of that error (the final Newton
// correction propagated to the metric, CV-Q8). From it:
//   - the observed orders p_k = ln(e_{k+1} / e_k) / ln(h_{k+1} / h_k) (verify/reference/metrics), any ratio;
//   - p-hat, the order of the finest pair (14.7);
//   - the iterative-error rule: U_iter <= U_disc / 10 at every level, U_disc the level's error itself (the
//     solution is exact), else the order is not evaluated;
//   - the asymptotic-range check from 4 levels (14.7; TVR-1's metrics): successive orders within 0.1. TVR-1's
//     GCI-ratio test does not apply to errors against an exact solution (they converge to 0, rho_k = r^-p);
//     with 3 levels it is not evaluated and p-hat is judged on the finest pair.
// The judgement against an order band [p_th - 0.2, p_th + 0.3] (§4.4 E (a)) is judge_order.
//
// Unstructured spacing (CV-Q7): h = (V / N)^(1/D), V the mesh's total control volume and N its node count.
//
// The library judges numbers only: the benchmarks build the meshes and run the solves (TVD-19: it depends on
// verify_reference and mesh, never on physics, assemble or solve).
#pragma once

#include <span>
#include <string>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::verify::study {

inline constexpr double iterative_fraction = 0.1;  // U_iter <= U_disc / 10 (CV-Q8, §5.2)
inline constexpr int asymptotic_levels = 4;        // the asymptotic-range check needs 4 levels (14.7)

// (V / N)^(1/D) for a mesh of dimension D, V the sum of its control volumes.
[[nodiscard]] double mean_spacing(const mesh::Mesh& m);

enum class StudyStatus { ok, invalid_input, not_evaluated };

struct OrderStudy {
    StudyStatus status = StudyStatus::ok;
    std::string reason;                 // empty when ok
    std::vector<double> h, error;       // finest first
    std::vector<double> orders;         // p_k from levels k, k + 1
    double observed = 0.0;              // p-hat: orders[0]
    double iterative_ratio = 0.0;       // max_k U_iter_k / e_k
    bool asymptotic_evaluated = false;  // >= 4 levels
    bool in_asymptotic_range = false;
    double order_spread = 0.0;          // max_k |p_{k+1} - p_k|
};

// h strictly increasing (finest first), e > 0 and finite, u_iter >= 0 and finite, all of one size >= 2.
// invalid_input otherwise; not_evaluated (with the reason) when an order is not finite.
[[nodiscard]] OrderStudy order_study(std::span<const double> h, std::span<const double> e,
                                     std::span<const double> u_iter);

struct OrderBand {
    double p_th;          // the scheme's order for the norm (§5.3)
    double below = 0.2;   // p-hat >= p_th - below
    double above = 0.3;   // p-hat <= p_th + above
};

struct OrderVerdict {
    bool judged = false;  // false: not evaluated (reason says why)
    bool pass = false;
    std::string reason;
};

// Not judged when the study is not ok, when the iterative error is not below a tenth of the discretisation
// error, or (4 levels or more) when the ladder is not in the asymptotic range; otherwise pass iff p-hat is in
// the band.
[[nodiscard]] OrderVerdict judge_order(const OrderStudy& s, const OrderBand& band);

}  // namespace NiTCAD::verify::study
