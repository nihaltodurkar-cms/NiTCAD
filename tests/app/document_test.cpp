// The deck document (Unit 26, WG3): load, save, the modified flag, and errors at their lines.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "NiTCAD/app/document.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

std::string text_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

}  // namespace

TEST_CASE("app document: the line an error names", "[app]") {
    const auto e = [](std::string m) { return base::Error{base::ErrorCode::invalid_input, std::move(m), std::nullopt}; };
    CHECK(app::error_line(e("d.deck:12: unknown key 'x' in [mesh]"), "d.deck") == 12);
    CHECK(app::error_line(e("d.deck:3: [mesh] generate2d: something"), "d.deck") == 3);
    CHECK(app::error_line(e("other.deck:12: x"), "d.deck") == 0);
    CHECK(app::error_line(e("x.deck:12: another deck of the same name length"), "d.deck") == 0);
    CHECK(app::error_line(e("d.deck: no line"), "d.deck") == 0);
    CHECK(app::error_line(e("d.deck:12 no colon"), "d.deck") == 0);
    CHECK(app::error_line(e("d.deck:0: zero"), "d.deck") == 0);
    CHECK(app::error_line(e("no source at all"), "d.deck") == 0);
}

TEST_CASE("app document: a check reports what nitcad check reports, at its line", "[app]") {
    const fs::path dir = app_test::scratch("document_check");
    const fs::path deck = dir / "d.deck";
    {
        const auto r = app::check_deck(app_test::diode_sweep(), deck);
        CHECK(r.diagnostics.empty());
        REQUIRE(r.plan);
        CHECK(r.plan->points.size() == 5);
        // Relative outputs are relative to the deck's directory.
        CHECK(r.plan->outputs.sweep == dir / "diode.nitsweep");
    }
    {
        // A meaning error: an unknown key on line 4 of this text.
        std::string text = app_test::diode_device;
        text.insert(text.find("x = linspace"), "max_sise = 1um\n");
        const int line = 4;  // "# A 1D...", "[mesh]", "type = tensor", then the inserted key
        const auto r = app::check_deck(text + "[sweep]\ncontact anode = 0\n", deck);
        CHECK_FALSE(r.plan);
        REQUIRE(r.diagnostics.size() == 1);
        CHECK(r.diagnostics[0].line == line);
        CHECK(r.diagnostics[0].message.starts_with("d.deck:4: "));
        CHECK(r.diagnostics[0].message.find("max_sise") != std::string::npos);
    }
    {
        // A syntax error: a missing '=' on line 2.
        const auto r = app::check_deck("[mesh]\ntype tensor\n", deck);
        REQUIRE(r.diagnostics.size() == 1);
        CHECK(r.diagnostics[0].line == 2);
    }
    {
        // An untitled deck is "<untitled>".
        const auto r = app::check_deck("[mesh]\ntype tensor\n", {});
        REQUIRE(r.diagnostics.size() == 1);
        CHECK(r.diagnostics[0].message.starts_with("<untitled>:2: "));
        CHECK(r.diagnostics[0].line == 2);
    }
}

TEST_CASE("app document: open, edit, save, save as", "[app]") {
    const fs::path dir = app_test::scratch("document_files");
    const fs::path a = app_test::write(dir, "a.deck", app_test::diode_sweep());
    auto d = app::DeckDocument::open(a);
    REQUIRE(d);
    CHECK(d->text() == app_test::diode_sweep());
    CHECK_FALSE(d->modified());
    CHECK(d->title() == "a.deck");
    d->set_text(d->text());  // the same text: not modified
    CHECK_FALSE(d->modified());
    d->set_text("# changed\n" + d->text());
    CHECK(d->modified());
    CHECK(d->title() == "a.deck *");
    REQUIRE(d->save());
    CHECK_FALSE(d->modified());
    CHECK(text_of(a) == d->text());
    CHECK_FALSE(fs::exists(dir / "a.deck.saving"));
    // Save as: the new path becomes the document's; the old file keeps what it had.
    d->set_text("# again\n");
    REQUIRE(d->save_as(dir / "b.deck"));
    CHECK(d->path() == dir / "b.deck");
    CHECK(text_of(dir / "b.deck") == "# again\n");
    CHECK(text_of(a).starts_with("# changed\n"));
    // A failed save keeps the path, the text and the modified flag.
    d->set_text("# third\n");
    CHECK_FALSE(d->save_as(dir / "missing" / "c.deck"));
    CHECK(d->path() == dir / "b.deck");
    CHECK(d->modified());
    // An untitled deck cannot be saved without a name; a missing file does not open.
    auto u = app::DeckDocument::untitled("[mesh]\n");
    CHECK(u.is_untitled());
    CHECK(u.title() == "Untitled");
    CHECK_FALSE(u.save());
    CHECK_FALSE(app::DeckDocument::open(dir / "none.deck"));
}
