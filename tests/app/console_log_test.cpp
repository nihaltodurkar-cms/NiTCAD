// The console's text (Unit 26j): every Newton iteration once, each point summarised once, the totals and
// the build errors at the end, and the oldest lines dropped past the cap.
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "NiTCAD/app/console_log.hpp"

using namespace NiTCAD;
using app::ConsoleLog;
using app::Diagnostic;
using app::Iteration;
using app::PointRow;

namespace {

std::vector<std::string> text(const ConsoleLog& log) { return {log.lines().begin(), log.lines().end()}; }

}  // namespace

TEST_CASE("app console: the iterations and points of a run, followed as they arrive", "[app]") {
    ConsoleLog log;
    log.add("Run sweep.deck");
    std::vector<PointRow> rows{{.phase = solve::Phase::equilibrium, .iterations = 1, .update = 0.5}};
    std::vector<Iteration> its{{.number = 1, .update = 0.5, .row = 0}};
    log.follow(rows, its);  // not following before begin_run
    CHECK(log.total() == 1);

    log.begin_run();
    log.follow(rows, its);
    CHECK(text(log) == std::vector<std::string>{"Run sweep.deck", "  equilibrium", "    Newton   1  update 5.000e-01"});
    log.follow(rows, its);  // nothing new: nothing added
    CHECK(log.total() == 3);

    // The equilibrium converges and the first bias point starts in the same batch.
    rows[0] = {.phase = solve::Phase::equilibrium, .iterations = 2, .update = 1e-9, .residual = 1e-12, .converged = true,
               .seconds = 0.5};
    rows.push_back({.phase = solve::Phase::bias, .point = 0, .point_count = 2, .iterations = 1, .update = 0.1});
    its.push_back({.number = 2, .update = 1e-9, .row = 0});
    its.push_back({.number = 3, .update = 0.1, .row = 1});
    log.follow(rows, its);
    CHECK(text(log) == std::vector<std::string>{
                           "Run sweep.deck", "  equilibrium", "    Newton   1  update 5.000e-01",
                           "    Newton   2  update 1.000e-09",
                           "    converged in 2 iterations, update 1e-09, residual 1e-12, 0.50 s", "  bias point 1 of 2",
                           "    Newton   1  update 1.000e-01"});

    // The run stops in the bias point: finish summarises it unconverged, then the totals and the errors.
    rows[1].seconds = 0.25;
    const std::vector<Diagnostic> errors{{.line = 7, .message = "bad key"}, {.line = 0, .message = "no contact"}};
    log.finish(rows, its, 1.5, errors);
    const auto t = text(log);
    REQUIRE(t.size() == 11);
    CHECK(t[7] == "    not converged after 1 iteration, update 0.1, residual 0, 0.25 s");
    CHECK(t[8] == "  2 points, 3 Newton iterations, 1.50 s");
    CHECK(t[9] == "  line 7: bad key");
    CHECK(t[10] == "  no contact");
    log.finish(rows, its, 1.5, errors);  // finished: nothing more
    CHECK(log.total() == 11);
}

TEST_CASE("app console: a check with no points adds only its errors; the oldest lines are dropped past the cap",
          "[app]") {
    ConsoleLog log;
    log.begin_run();
    const std::vector<Diagnostic> errors{{.line = 3, .message = "unknown section"}};
    log.finish({}, {}, 0.1, errors);
    CHECK(text(log) == std::vector<std::string>{"  line 3: unknown section"});

    for (std::size_t k = 1; k <= ConsoleLog::max_lines; ++k) log.add(std::to_string(k));
    CHECK(log.total() == ConsoleLog::max_lines + 1);
    CHECK(log.lines().size() == ConsoleLog::keep_lines);
    CHECK(log.dropped() == ConsoleLog::max_lines + 1 - ConsoleLog::keep_lines);
    CHECK(log.lines().back() == std::to_string(ConsoleLog::max_lines));
    log.clear();
    CHECK(log.lines().empty());
    CHECK(log.total() == ConsoleLog::max_lines + 1);
}
