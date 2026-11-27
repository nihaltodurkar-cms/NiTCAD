// The guide pane (Unit 26p, GU3, DECISIONS.md): View > Guide shows it at the main window's right edge, beside every
// tab. A lesson picker, the step's heading and prose, its action as a button, what the action gave ("What happened"),
// Back and Next. The pane runs nothing itself: the action button hands the step to the host (the main window), which
// edits, sweeps or runs through the paths it already has and reports back with set_report or set_error. The lessons
// are app/lesson.hpp's.
#pragma once

#include <windows.h>

#include <cstddef>
#include <string>
#include <vector>

#include "NiTCAD/app/lesson.hpp"

namespace NiTCAD::app {

class GuidePane {
public:
    struct Host {
        virtual void guide_open(std::size_t lesson) = 0;     // a lesson picked
        virtual void guide_act(const LessonStep& step) = 0;  // the step's action button
    protected:
        ~Host() = default;
    };

    GuidePane(HWND parent, Host& host, int id);
    ~GuidePane();
    GuidePane(const GuidePane&) = delete;
    GuidePane& operator=(const GuidePane&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return hwnd_; }
    void set_lessons(std::vector<Lesson> lessons);
    [[nodiscard]] const std::vector<Lesson>& lessons() const noexcept { return lessons_; }
    void show_lesson(std::size_t lesson);  // its first step; the picker follows
    void go(int delta);                    // Back (-1), Next (+1), clamped to the lesson's steps
    void set_report(const std::string& text);
    void set_error(const std::string& text);  // shown in the report's place; the step stays
    void set_busy(bool busy);                 // the action button off while a run is busy
    void set_dpi(int dpi, HFONT font, HFONT bold);
    void act();  // as the action button

    [[nodiscard]] std::size_t lesson() const noexcept { return lesson_; }
    [[nodiscard]] std::size_t step() const noexcept { return step_; }
    [[nodiscard]] const LessonStep* current() const noexcept;
    // "Sweep Na = 1e16, 1e17, 1e18", "Set tox = 2nm", "Set field_mobility = off", "Run"; empty without an action.
    [[nodiscard]] std::string action_label() const;
    // What the controls show (for the tests).
    [[nodiscard]] std::string heading_text() const;
    [[nodiscard]] std::string prose_text() const;
    [[nodiscard]] std::string report_text() const;
    [[nodiscard]] std::string position_text() const;  // "Step 2 of 5"
    [[nodiscard]] HWND action_button() const noexcept { return action_; }
    [[nodiscard]] HWND prose_edit() const noexcept { return prose_; }

private:
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    void layout();
    void show_step();

    Host& host_;
    HWND hwnd_ = nullptr;
    HWND picker_ = nullptr, heading_ = nullptr, prose_ = nullptr, action_ = nullptr, report_label_ = nullptr,
         report_ = nullptr, back_ = nullptr, next_ = nullptr, position_ = nullptr;
    std::vector<Lesson> lessons_;
    std::size_t lesson_ = 0, step_ = 0;
    bool busy_ = false;
    bool error_ = false;
    int dpi_ = 96;
    HFONT font_ = nullptr, bold_ = nullptr;
};

}  // namespace NiTCAD::app
