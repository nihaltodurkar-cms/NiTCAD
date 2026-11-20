// The console's text (Unit 26j): the banner once the plan is built, every Newton iteration once in a table,
// each point summarised once (green when it converged, amber when not), the summary and the build errors at
// the end; the filters; the rich text; and the oldest lines dropped past the cap.
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "NiTCAD/app/console_log.hpp"

using namespace NiTCAD;
using app::ConsoleLog;
using app::Diagnostic;
using app::Iteration;
using app::LogFilter;
using app::LogLevel;
using app::LogLine;
using app::PointRow;
using app::RunView;

namespace {

std::vector<LogLine> lines(const ConsoleLog& log) { return {log.lines().begin(), log.lines().end()}; }

}  // namespace

TEST_CASE("app console: the banner, the Newton table and the points of a run, followed as they arrive", "[app]") {
    ConsoleLog log;
    log.add(LogLevel::header, "Run sweep.deck");
    const app::PlanSummary summary{.description = {"Analysis   bias sweep"},
                                   .contacts = {"anode", "cathode"},
                                   .points = {{0.0, 0.0}, {0.5, 0.0}}};
    std::vector<PointRow> rows{{.phase = solve::Phase::equilibrium, .iterations = 1, .update = 0.5, .residual = 2.0}};
    std::vector<Iteration> its{{.number = 1, .update = 0.5, .row = 0, .residual = 2.0, .in_row = 1}};
    log.follow(RunView{.rows = rows, .iterations = its, .summary = &summary});  // not following before begin_run
    CHECK(log.total() == 1);

    log.begin_run();
    log.follow(RunView{.rows = rows, .iterations = its, .summary = &summary});
    const std::string rule = "  ------------------------------------------------------------";
    CHECK(lines(log) == std::vector<LogLine>{{LogLevel::header, "Run sweep.deck"},
                                             {LogLevel::header, rule},
                                             {LogLevel::info, "   Analysis   bias sweep"},
                                             {LogLevel::header, rule},
                                             {LogLevel::info, "  Equilibrium"},
                                             {LogLevel::detail, "     Iter      |update|    |residual|"},
                                             {LogLevel::detail, "        1    5.000e-01    2.000e+00"}});
    log.follow(RunView{.rows = rows, .iterations = its, .summary = &summary});  // nothing new: nothing added
    CHECK(log.total() == 7);

    // The equilibrium converges and the second bias point starts in the same batch.
    rows[0] = {.phase = solve::Phase::equilibrium, .iterations = 2, .update = 1e-9, .residual = 1e-12, .converged = true,
               .seconds = 0.5};
    rows.push_back({.phase = solve::Phase::bias, .point = 1, .point_count = 2, .iterations = 1, .update = 0.1});
    its.push_back({.number = 2, .update = 1e-9, .row = 0, .residual = 1e-12, .in_row = 2});
    its.push_back({.number = 3, .update = 0.1, .row = 1, .residual = 3.0, .in_row = 1});
    log.follow(RunView{.rows = rows, .iterations = its, .summary = &summary});
    auto l = lines(log);
    REQUIRE(l.size() == 12);
    CHECK(l[7] == LogLine{LogLevel::detail, "        2    1.000e-09    1.000e-12"});
    CHECK(l[8] == LogLine{LogLevel::success, "     Converged in 2 iterations: |update| 1.00e-09, |residual| 1.00e-12, 0.50 s"});
    CHECK(l[9] == LogLine{LogLevel::info, "  Bias point 2/2   V(anode) = 0.5 V, V(cathode) = 0 V"});
    CHECK(l[11] == LogLine{LogLevel::detail, "        1    1.000e-01    3.000e+00"});

    // The run stops in the bias point: finish summarises it unconverged, then the summary, why, the errors.
    rows[1].seconds = 0.25;
    const std::vector<Diagnostic> errors{{.line = 7, .message = "bad key"}, {.line = 0, .message = "no contact"}};
    const base::Error stopped{base::ErrorCode::cancelled, "cancelled", std::nullopt};
    const RunView end{.rows = rows, .iterations = its, .summary = &summary, .seconds = 1.5, .diagnostics = errors,
                      .stopped = &stopped};
    log.finish(end);
    l = lines(log);
    REQUIRE(l.size() == 18);
    CHECK(l[12] == LogLine{LogLevel::warning, "     Not converged after 1 iteration: |update| 1.00e-01, |residual| 0.00e+00, 0.25 s"});
    CHECK(l[13] == LogLine{LogLevel::header, "  Summary"});
    CHECK(l[14] == LogLine{LogLevel::info, "  2 points (1 converged), 3 Newton iterations, 1.5 per point, 1.50 s"});
    CHECK(l[15] == LogLine{LogLevel::warning, "  Stopped early: cancelled"});
    CHECK(l[16] == LogLine{LogLevel::error, "  line 7: bad key"});
    CHECK(l[17] == LogLine{LogLevel::error, "  no contact"});
    log.finish(end);  // finished: nothing more
    CHECK(log.total() == 18);

    // The filters, and the text Copy and Save log take.
    CHECK(app::shows(LogFilter::all, LogLevel::detail));
    CHECK_FALSE(app::shows(LogFilter::summary, LogLevel::detail));
    CHECK(app::shows(LogFilter::summary, LogLevel::success));
    CHECK(app::shows(LogFilter::warnings, LogLevel::warning));
    CHECK_FALSE(app::shows(LogFilter::warnings, LogLevel::info));
    CHECK_FALSE(app::shows(LogFilter::errors, LogLevel::warning));
    CHECK(app::filter_names().size() == 4);
    CHECK(log.text(LogFilter::errors) == "  line 7: bad key\n  no contact\n");
    CHECK(log.text(LogFilter::warnings).starts_with("     Not converged after 1 iteration"));
}

TEST_CASE("app console: the rich text escapes, colours by level and continues; the cap drops the oldest lines", "[app]") {
    ConsoleLog log;
    log.add(LogLevel::detail, "a{b}\\c");
    log.add(LogLevel::error, "\xCE\xBC" "m");  // "μm"
    const std::string all = app::console_rtf(log.lines(), 0, LogFilter::all, false);
    CHECK(all.starts_with("{\\rtf1"));
    CHECK(all.ends_with("}"));
    CHECK(all.find("\\cf1\\b0 a\\{b\\}\\\\c\\par\n\\cf6\\b \\u956?m}") != std::string::npos);
    const std::string errors = app::console_rtf(log.lines(), 0, LogFilter::errors, true);
    CHECK(errors.find("a\\{b") == std::string::npos);
    CHECK(errors.find("\\fs18 \\par\n\\cf6") != std::string::npos);  // continues: a break before the first line

    log.begin_run();
    const std::vector<Diagnostic> errs{{.line = 3, .message = "unknown section"}};
    log.finish(RunView{.diagnostics = errs});
    CHECK(log.lines().back() == LogLine{LogLevel::error, "  line 3: unknown section"});

    for (std::size_t k = 1; k <= ConsoleLog::max_lines; ++k) log.add(LogLevel::info, std::to_string(k));
    CHECK(log.total() == ConsoleLog::max_lines + 3);
    // Dropped down to keep_lines when the 20001st line came (3 + 19998 added), then two more.
    CHECK(log.lines().size() == ConsoleLog::keep_lines + 2);
    CHECK(log.dropped() == ConsoleLog::max_lines + 1 - ConsoleLog::keep_lines);
    CHECK(log.lines().back().text == std::to_string(ConsoleLog::max_lines));
    log.clear();
    CHECK(log.lines().empty());
    CHECK(log.total() == ConsoleLog::max_lines + 3);
}

TEST_CASE("app console: the Points table, its live and pending rows, and its text (design D)", "[app]") {
    const app::PlanSummary summary{.description = {}, .contacts = {"drain", "gate"},
                                   .points = {{0.05, 0.0}, {0.05, 0.5}, {0.05, 1.0}}};
    const std::vector<PointRow> rows{
        {.phase = solve::Phase::equilibrium, .iterations = 3, .update = 1e-9, .residual = 1e-12, .converged = true, .seconds = 0.4},
        {.phase = solve::Phase::bias, .point = 0, .point_count = 3, .iterations = 2, .update = 2e-9, .residual = 3e-12,
         .converged = true, .seconds = 0.2},
        {.phase = solve::Phase::bias, .point = 1, .point_count = 3, .iterations = 1, .update = 0.25, .residual = 7.0,
         .seconds = 0.1}};
    const std::vector<Iteration> its(6);
    auto t = app::points_table(RunView{.rows = rows, .iterations = its, .summary = &summary, .seconds = 61.0, .running = true});
    CHECK(t.columns == std::vector<std::string>{"Point", "V(gate)", "It", "|update|", "|residual|", "Time (s)", "Status"});
    CHECK(t.numeric == std::vector<bool>{false, true, true, true, true, true, false});
    REQUIRE(t.rows.size() == 4);
    CHECK(t.rows[0] == std::vector<std::string>{"Eq", "–", "3", "1.00e-09", "1.00e-12", "0.40", "converged"});
    CHECK(t.rows[1] == std::vector<std::string>{"1", "0.000", "2", "2.00e-09", "3.00e-12", "0.20", "converged"});
    CHECK(t.rows[2] == std::vector<std::string>{"2", "0.500", "1", "2.50e-01", "7.00e+00", "0.10", "solving"});
    CHECK(t.rows[3] == std::vector<std::string>{"3", "1.000", "", "", "", "", ""});
    CHECK(t.states == std::vector<app::PointState>{app::PointState::done, app::PointState::done, app::PointState::live,
                                                   app::PointState::pending});
    CHECK(t.progress == "Bias point 2 of 3 · 6 Newton iterations · 1:01");

    // Stopped there: the point that did not converge is a warning, the rest "not run".
    t = app::points_table(RunView{.rows = rows, .iterations = its, .summary = &summary, .seconds = 1.0});
    CHECK(t.states[2] == app::PointState::warning);
    CHECK(t.rows[2].back() == "not converged");
    CHECK(t.rows[3].back() == "not run");
    // Planned but not run (a check): every point pending, blank status.
    t = app::points_table(RunView{.summary = &summary});
    CHECK(t.rows.size() == 3);
    CHECK(t.rows[0].back().empty());
    CHECK(t.progress.empty());

    const app::PointsTable small{.columns = {"Point", "Note"}, .numeric = {false, false},
                                 .rows = {{"1", "a,b"}, {"2", "say \"hi\""}}, .states = {}, .progress = {}};
    CHECK(app::points_text(small, ',') == "Point,Note\n1,\"a,b\"\n2,\"say \"\"hi\"\"\"\n");
    CHECK(app::points_text(small, '\t') == "Point\tNote\n1\ta,b\n2\t\"say \"\"hi\"\"\"\n");
}