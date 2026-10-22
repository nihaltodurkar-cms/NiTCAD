// The command list (Unit 26, WG5): every command in a menu with a unique mnemonic, accelerators unique,
// the menu bar in the reference's order (VQ13).
#include <catch2/catch_test_macros.hpp>

#include <cctype>
#include <map>
#include <set>
#include <string>
#include <tuple>

#include "NiTCAD/app/commands.hpp"

using namespace NiTCAD;

namespace {

// The letter after '&', lower-cased; 0 without one.
char mnemonic(std::string_view text) {
    const auto amp = text.find('&');
    if (amp == std::string_view::npos || amp + 1 >= text.size()) return 0;
    return static_cast<char>(std::tolower(static_cast<unsigned char>(text[amp + 1])));
}

}  // namespace

TEST_CASE("app commands: menus, mnemonics and accelerators (WG5)", "[app]") {
    const auto bar = app::menu_bar();
    REQUIRE(bar.size() == 7);
    CHECK(bar[0] == "&File");
    CHECK(bar[2] == "&Device");
    CHECK(bar[3] == "&Mesh");
    CHECK(bar[4] == "&Physics");
    CHECK(bar[5] == "&Simulation");
    CHECK(bar[6] == "&View");
    std::set<char> top;
    for (const auto& m : bar) CHECK(top.insert(mnemonic(m)).second);

    std::set<int> ids;
    std::map<std::string, std::set<char>> letters;
    std::set<std::tuple<int, bool, bool, bool>> keys;
    for (const app::CommandInfo& c : app::commands()) {
        CAPTURE(c.text);
        CHECK(ids.insert(static_cast<int>(c.id)).second);
        CHECK(std::ranges::find(bar, c.menu) != bar.end());
        CHECK(app::find_command(c.id) == &c);
        // Every item has a mnemonic, unique in its menu, except the pan moves, which have keys.
        const char m = mnemonic(c.text);
        if (m != 0) CHECK(letters[std::string(c.menu)].insert(m).second);
        else CHECK(std::string_view(c.text).starts_with("Pan"));
        // Accelerators unique; every command with keys shows them.
        if (c.accelerator.key != 0) {
            CHECK(keys.insert({c.accelerator.key, c.accelerator.ctrl, c.accelerator.shift, c.accelerator.alt}).second);
            CHECK_FALSE(c.accelerator_text.empty());
        }
    }
    CHECK(app::find_command(app::Command::section_first) == nullptr);
    // The commands the plan names are all there.
    for (const app::Command c : {app::Command::open_deck, app::Command::save, app::Command::check, app::Command::run,
                                 app::Command::stop, app::Command::fit, app::Command::probe_mode,
                                 app::Command::measure_mode, app::Command::perspective, app::Command::next_pane}) {
        CHECK(app::find_command(c) != nullptr);
    }
}
