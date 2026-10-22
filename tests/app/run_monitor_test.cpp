// The run monitor's content (Unit 26, VA3): the rows' cells, the convergence plot, the progress.
#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "NiTCAD/app/run_monitor.hpp"

using namespace NiTCAD;
using app::PointRow;

TEST_CASE("app run monitor: cells, the convergence plot and the progress", "[app]") {
    CHECK(app::monitor_columns().size() == 8);
    const PointRow eq{.phase = solve::Phase::equilibrium, .iterations = 8, .update = 1e-9, .residual = 2e-12,
                      .converged = true, .seconds = 0.25};
    const auto c = app::monitor_cells(eq, 0);
    REQUIRE(c.size() == 8);
    CHECK(c[0] == "1");
    CHECK(c[1] == "equilibrium");
    CHECK(c[3] == "8");
    CHECK(c[4] == "1e-09");
    CHECK(c[6] == "yes");
    CHECK(c[7] == "0.25");
    const PointRow bias{.phase = solve::Phase::bias, .point = 2, .point_count = 11, .iterations = 4};
    CHECK(app::monitor_cells(bias, 3)[2] == "3 of 11");
    CHECK(app::monitor_cells(bias, 3)[6].empty());
    const PointRow tr{.phase = solve::Phase::transient, .point = 7, .time_s = 2.5e-9};
    CHECK(app::monitor_cells(tr, 0)[2] == "t = 2.5e-09 s");
    const PointRow ac{.phase = solve::Phase::small_signal, .point = 0, .point_count = 3, .frequency_Hz = 1e6};
    CHECK(app::monitor_cells(ac, 0)[2] == "1 of 3, 1e+06 Hz");

    // The plot: one series of the positive updates on a log axis.
    const std::vector<app::Iteration> its{{1, 1.0, 0}, {2, 1e-3, 0}, {3, 0.0, 0}, {4, 1e-9, 1}};
    auto plot = app::convergence_plot(its);
    REQUIRE(plot);
    REQUIRE(plot->series.size() == 1);
    CHECK(plot->series[0].x == std::vector<double>{1, 2, 4});
    CHECK(plot->y.scale == render::AxisScale::log10);
    CHECK_FALSE(app::convergence_plot(std::vector<app::Iteration>{{1, 1.0, 0}}));

    // Progress: converged bias points over the planned ones; unknown without a plan.
    std::vector<PointRow> rows{eq, {.phase = solve::Phase::bias, .point = 0, .converged = true},
                               {.phase = solve::Phase::bias, .point = 1, .converged = true},
                               {.phase = solve::Phase::bias, .point = 2, .converged = false}};
    CHECK(app::progress_fraction(rows, 4) == 0.5);
    CHECK(app::progress_fraction(rows, 1) == 1.0);
    CHECK(app::progress_fraction(rows, 0) == -1.0);
}
