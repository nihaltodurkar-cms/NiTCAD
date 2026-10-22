// The per-user state (Unit 26, VA1): recent files and the layout, kept across sessions.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>

#include "NiTCAD/app/app_state.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

TEST_CASE("app state: recent files, most recent first, no duplicates, at most eight", "[app]") {
    app::AppState s;
    for (int k = 0; k < 10; ++k) s.add_recent("f" + std::to_string(k) + ".deck");
    REQUIRE(s.recent.size() == app::AppState::max_recent);
    CHECK(s.recent.front() == "f9.deck");
    CHECK(s.recent.back() == "f2.deck");
    s.add_recent("f5.deck");
    CHECK(s.recent.front() == "f5.deck");
    CHECK(s.recent.size() == app::AppState::max_recent);
    s.remove_recent("f5.deck");
    CHECK(s.recent.size() == app::AppState::max_recent - 1);
    CHECK(s.recent.front() == "f9.deck");
}

TEST_CASE("app state: written and read back; damaged lines ignored", "[app]") {
    app::AppState s;
    s.add_recent(fs::path(u8"C:/decks/b\u00e4nder.deck"));
    s.add_recent("C:/decks/a.deck");
    s.layout = {10, 20, 1280, 800, true, 230, 270, 0};
    const app::AppState r = app::parse_state(app::serialize(s));
    CHECK(r.recent == s.recent);
    CHECK(r.layout.x == 10);
    CHECK(r.layout.y == 20);
    CHECK(r.layout.width == 1280);
    CHECK(r.layout.height == 800);
    CHECK(r.layout.maximized);
    CHECK(r.layout.tree_width == 230);
    CHECK(r.layout.inspector_width == 270);
    // Unknown keys, bad numbers and lines without '=' change nothing.
    const app::AppState d = app::parse_state("window = 1 2 x 4\nunknown = 3\ngarbage\ntree_width = -5\nrecent = a.deck\n"
                                             "recent = a.deck\nmaximized = 1\n");
    CHECK(d.layout.width == 0);
    CHECK(d.layout.tree_width == 0);
    CHECK(d.layout.maximized);
    CHECK(d.recent == std::vector<fs::path>{"a.deck"});

    const fs::path dir = app_test::scratch("state");
    REQUIRE(app::save_state(s, dir / "sub" / "app.ini"));
    CHECK(app::load_state(dir / "sub" / "app.ini").recent == s.recent);
    CHECK(app::load_state(dir / "none.ini").recent.empty());
}
