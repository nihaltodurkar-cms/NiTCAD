// Curves read from results, and the tools every extractor shares (ARCHITECTURE.md section 11,
// Unit 24; legacy gui/services/characterization.py and sweep_derived.py).
//
// The analysis layer reads results and never solves: every extractor is a pure function of data
// already computed. An extraction returns its value with the points it came from, and says when
// it was extrapolated past the data or used a current at or below the state's current resolution
// (BiasPoint::terminal_current_resolution: there the current is rounding noise). A curve that
// cannot give the quantity is an error with the reason (invalid_input), never a NaN.
// OLD / NEW / REASON: OLD, the legacy extractors returned NaN or None without a reason and sorted
// their input by x, dropping non-finite points. NEW, input is checked and kept in its own order.
// REASON, a curve traced through a fold (trace_bias) is not a function of its bias, and sorting it
// scrambles the branch; a reason makes a failed extraction actionable.
//
// Units follow the results: bias in V, a contact's current in A cm^(D-3) (A/cm^2 in 1D, A/cm in
// 2D, A in 3D), charge in C cm^(D-3), capacitance and conductance per V likewise.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string_view>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/small_signal.hpp"

namespace NiTCAD::analysis {

// Point indices into the curve (or the frequency list) an extraction used, inclusive.
struct Window {
    std::size_t first = 0;
    std::size_t last = 0;
};

struct Extraction {
    double value = 0.0;
    std::string_view unit;
    Window window;
    bool extrapolated = false;      // the value lies outside the data
    bool below_resolution = false;  // a point used has |y| at or below its resolution
};

// y against x in the order of the run. resolution, when not empty, holds one bound per point.
struct Curve {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> resolution;
};

// Errors (invalid_input): x and y of different sizes or fewer than 2 points; a value not finite
// (the error's index the point); resolution neither empty nor one finite non-negative value per
// point.
[[nodiscard]] std::expected<Curve, base::Error> make_curve(std::vector<double> x,
                                                           std::vector<double> y,
                                                           std::vector<double> resolution = {});

// The swept contact's bias against the sense contact's terminal current, with its resolution
// (left empty if any point lacks it). Errors: as make_curve; a contact out of range for a point.
[[nodiscard]] std::expected<Curve, base::Error> current_curve(const results::Sweep& sweep,
                                                              std::size_t swept,
                                                              std::size_t sense);

// The swept contact's bias against a gate's charge (BiasPoint::gate_charge). Errors as above.
[[nodiscard]] std::expected<Curve, base::Error> gate_charge_curve(const results::Sweep& sweep,
                                                                  std::size_t swept,
                                                                  std::size_t gate);

// The swept contact's bias against Im Y_row,column / omega at one frequency over the operating
// points. Errors as above; frequency out of range or zero.
[[nodiscard]] std::expected<Curve, base::Error> capacitance_curve(
    const results::SmallSignal& small_signal, std::size_t frequency, std::size_t row,
    std::size_t column, std::size_t swept);

// y and the resolution times factor (the resolution by |factor|): a sign flip (a p-channel or a
// reverse current) or a size (cm^(3-D)) to absolute units.
[[nodiscard]] Curve scaled(Curve curve, double factor);

// The points with x between the two bounds (either order; widened by 1e-9 of the distance between
// them, so a bias summed from steps, 0.05 * 14 = 0.7000000000000001, is not lost), in their order.
// Errors: fewer than 2.
[[nodiscard]] std::expected<Curve, base::Error> slice(const Curve& curve, double from_x,
                                                      double to_x);

[[nodiscard]] bool strictly_monotone(std::span<const double> x) noexcept;

// dy/dx at every point: second order on uneven spacing inside, first order one-sided at the ends
// (numpy.gradient with edge_order 1). Errors: x not strictly monotone.
[[nodiscard]] std::expected<std::vector<double>, base::Error> derivative(const Curve& curve);

enum class Scale : std::uint8_t {
    linear,
    logarithmic,  // in ln |y|: exact for an exponential between the points
};

// y at x, interpolated between the bracketing points (the window those two), or y of the point at
// x exactly (the window that point). Errors: x not strictly monotone; x outside
// the curve; logarithmic with a bracketing y zero or the two of opposite signs.
[[nodiscard]] std::expected<Extraction, base::Error> value_at(const Curve& curve, double x,
                                                              Scale scale = Scale::linear);

// The x where y first reaches target in the curve's order (from either side), interpolated
// between the two points around it (x itself if the first point equals it); the window those
// points. Errors: y never reaches it; logarithmic with target zero or of the other sign than the
// bracketing points.
[[nodiscard]] std::expected<Extraction, base::Error> crossing(const Curve& curve, double target,
                                                              Scale scale = Scale::linear);

}  // namespace NiTCAD::analysis
