// The MOSFET guide's lessons (Unit 26p, GU2, GU5, DECISIONS.md).
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "NiTCAD/app/lesson.hpp"

using namespace NiTCAD::app;
namespace fs = std::filesystem;

namespace {

// A deck the lessons ride on: two parameters, so set and sweep have names to refer to.
constexpr const char* body = "\n[parameters]\nNa = 1e16  # body acceptors\ntox = 5nm  # oxide\n\n"
                             "[models]\nfield_mobility = on  # x\n";

std::string lesson_text(const std::string& head) { return head + body; }

LessonError error_of(const std::string& text) {
    const auto l = parse_lesson(text);
    REQUIRE_FALSE(l);
    return l.error();
}

}  // namespace

TEST_CASE("lesson: a full lesson parses", "[lesson]") {
    const auto l = parse_lesson(lesson_text("# lesson: Body doping\n"
                                            "# track: advanced\n"
                                            "# order: 10\n"
                                            "# step: First\n"
                                            "# A\n"
                                            "# b.\n"
                                            "#\n"
                                            "# C d.\n"
                                            "# step: Second\n"
                                            "# Sweep it.\n"
                                            "# do: sweep Na = 1e16, 1e17, 1e18\n"
                                            "# show: results\n"
                                            "# report: Threshold (max gm), Subthreshold swing\n"));
    REQUIRE(l);
    CHECK(l->title == "Body doping");
    CHECK(l->track == LessonTrack::advanced);
    CHECK(l->order == 10);
    REQUIRE(l->steps.size() == 2);
    CHECK(l->steps[0].heading == "First");
    CHECK(l->steps[0].prose == "A b.\n\nC d.");
    CHECK(l->steps[0].action.kind == LessonAction::Kind::none);
    CHECK(l->steps[0].line == 4);
    const LessonStep& s = l->steps[1];
    CHECK(s.action.kind == LessonAction::Kind::sweep);
    CHECK(s.action.section.empty());
    CHECK(s.action.name == "Na");
    CHECK(s.action.value == "1e16, 1e17, 1e18");
    CHECK(s.show == "results");
    CHECK(s.report == std::vector<std::string>{"Threshold (max gm)", "Subthreshold swing"});
}

TEST_CASE("lesson: each action form", "[lesson]") {
    const auto one = [](const std::string& action) {
        const auto l = parse_lesson(lesson_text("# lesson: L\n# track: basics\n# order: 1\n# step: S\n# do: " + action + "\n"));
        REQUIRE(l);
        REQUIRE(l->steps.size() == 1);
        return l->steps[0].action;
    };
    const auto set = one("set Na = 1e17");
    CHECK(set.kind == LessonAction::Kind::set);
    CHECK(set.name == "Na");
    CHECK(set.value == "1e17");
    CHECK(one("run").kind == LessonAction::Kind::run);
    const auto edit = one("edit [models] field_mobility = off");
    CHECK(edit.kind == LessonAction::Kind::edit);
    CHECK(edit.section == "models");
    CHECK(edit.name == "field_mobility");
    CHECK(edit.value == "off");
}

TEST_CASE("lesson: errors name their line", "[lesson]") {
    const std::string head = "# lesson: L\n# track: basics\n# order: 1\n# step: S\n";
    const auto unknown = error_of(lesson_text(head + "# do: dance\n"));
    CHECK(unknown.line == 5);
    CHECK(unknown.message.find("dance") != std::string::npos);
    CHECK(error_of(lesson_text("# lesson: L\n# track: basics\n# do: run\n")).line == 3);  // before the first step
    CHECK(error_of(lesson_text(head + "# do: run\n# do: run\n")).line == 6);              // a second action
    const auto lx = error_of(lesson_text(head + "# do: sweep Lx = 1, 2\n"));
    CHECK(lx.line == 5);
    CHECK(lx.message.find("Lx") != std::string::npos);
    CHECK(error_of(lesson_text(head + "# do: set Lx = 1\n")).line == 5);
    CHECK(error_of(lesson_text(head + "# do: sweep tox = 1e16, 1e17\n")).line == 5);  // numbers for a length
    CHECK(error_of(lesson_text("# lesson: L\n# track: expert\n# order: 1\n# step: S\n")).line == 2);
    CHECK(error_of(lesson_text(head + "# show: elsewhere\n")).line == 5);
    CHECK(error_of(lesson_text("# lesson: L\n# track: basics\n# order: 1\n")).line == 1);  // no step
    CHECK(error_of(std::string("# A plain deck\n") + body).line == 0);                                // not a lesson
}

TEST_CASE("lesson: lessons listed by track, then order", "[lesson]") {
    const fs::path dir = fs::temp_directory_path() / "nitcad_lesson_list_test";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const auto write = [&](const char* file, const std::string& text) { std::ofstream(dir / file, std::ios::binary) << text; };
    write("a.deck", lesson_text("# lesson: Third\n# track: advanced\n# order: 1\n# step: S\n"));
    write("b.deck", lesson_text("# lesson: Second\n# track: basics\n# order: 2\n# step: S\n"));
    write("c.deck", lesson_text("# lesson: First\n# track: basics\n# order: 1\n# step: S\n"));
    write("d.deck", std::string("# Not a lesson\n") + body);
    const auto ls = list_lessons(dir);
    REQUIRE(ls.size() == 3);
    CHECK(ls[0].title == "First");
    CHECK(ls[1].title == "Second");
    CHECK(ls[2].title == "Third");
    CHECK(ls[0].path == dir / "c.deck");
    CHECK(to_string(LessonTrack::intermediate) == "Intermediate");
    fs::remove_all(dir);
}

TEST_CASE("lesson: set_entry changes one line", "[lesson]") {
    const auto out = set_entry(body, "models", "field_mobility", "off");
    REQUIRE(out);
    CHECK(out->find("field_mobility = off  # x") != std::string::npos);
    CHECK(out->find("field_mobility = on") == std::string::npos);
    CHECK(out->find("Na = 1e16  # body acceptors") != std::string::npos);
    CHECK_FALSE(set_entry(body, "models", "auger", "on"));
    CHECK_FALSE(set_entry(body, "solver", "max_cuts", "6"));
}
