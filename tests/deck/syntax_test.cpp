// The deck's syntax (Unit 31, deck/syntax.hpp) and its value readers (src/deck/values.hpp).
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "NiTCAD/deck/syntax.hpp"
#include "values.hpp"

using namespace NiTCAD;
using base::ErrorCode;

namespace {

// The error of parsing `text`, which must fail.
base::Error parse_error(const std::string& text) {
    const auto d = deck::parse(text, "t.deck");
    REQUIRE_FALSE(d.has_value());
    CHECK(d.error().code == ErrorCode::invalid_input);
    return d.error();
}

}  // namespace

TEST_CASE("deck syntax: sections, entries, comments, quotes and continuations", "[deck]") {
    const std::string text =
        "\xEF\xBB\xBF# a deck\r\n"
        "[mesh]\r\n"
        "type = tensor   # trailing comment\r\n"
        "x = linspace(0, 2um, 3) + \\\r\n"
        "    3um, 4um\r\n"
        "\r\n"
        "[interface   silicon    oxide]\n"
        "trap_level = acceptor 1e10 0.1\n"
        "[contact anode]\n"
        "patch   =   x_min\n"
        "[output]\n"
        "sweep = \"out dir/a#b.nitsweep\"   # a quoted '#' is kept\n"
        "[sweep]\n"
        "contact    anode = 0 : 0.8 : 0.05\n"
        "contact anode = 1\n";
    const auto d = deck::parse(text, "t.deck", "base");
    REQUIRE(d.has_value());
    CHECK(d->source == "t.deck");
    CHECK(d->directory == std::filesystem::path("base"));
    REQUIRE(d->sections.size() == 5);
    const auto& mesh = d->sections[0];
    CHECK(mesh.kind == "mesh");
    CHECK(mesh.name.empty());
    CHECK(mesh.line == 2);
    REQUIRE(mesh.entries.size() == 2);
    CHECK(mesh.entries[0].key == "type");
    CHECK(mesh.entries[0].value == "tensor");
    CHECK(mesh.entries[0].line == 3);
    CHECK(mesh.entries[1].key == "x");
    CHECK(mesh.entries[1].value == "linspace(0, 2um, 3) +  3um, 4um");
    CHECK(mesh.entries[1].line == 4);
    CHECK(d->sections[1].kind == "interface");
    CHECK(d->sections[1].name == "silicon oxide");
    CHECK(d->sections[1].line == 7);
    CHECK(d->sections[2].entries[0].value == "x_min");
    CHECK(d->sections[3].entries[0].value == "\"out dir/a#b.nitsweep\"");
    CHECK(deck::unquote(d->sections[3].entries[0].value) == "out dir/a#b.nitsweep");
    // Keys collapse blanks; repeated keys are kept in order.
    REQUIRE(d->sections[4].entries.size() == 2);
    CHECK(d->sections[4].entries[0].key == "contact anode");
    CHECK(d->sections[4].entries[1].key == "contact anode");
    CHECK(d->sections[4].entries[1].line == 15);
}

TEST_CASE("deck syntax: errors name the line", "[deck]") {
    CHECK(parse_error("type = tensor\n").message == "t.deck:1: an entry before the first [section] header");
    CHECK(parse_error("[mesh]\n\n[1mesh]\n").message.starts_with("t.deck:3: malformed header"));
    CHECK(parse_error("[mesh\n").message == "t.deck:1: a header must end with ']'");
    CHECK(parse_error("[mesh]\ntype tensor\n").message ==
          "t.deck:2: expected 'key = value' or a [section] header");
    CHECK(parse_error("[mesh]\ntype =   # nothing\n").message == "t.deck:2: 'type' has no value");
    CHECK(parse_error("[mesh]\n = tensor\n").message.starts_with("t.deck:2: malformed key"));
    CHECK(parse_error("[mesh]\nty/pe = tensor\n").message.starts_with("t.deck:2: malformed key"));
    CHECK(parse_error("[output]\nsweep = \"a.nitsweep\n").message == "t.deck:2: unterminated quote");
    CHECK(parse_error("[mesh]\nx = 1, \\").message == "t.deck:2: continuation at the end of the file");
    const auto nul = parse_error(std::string("[mesh]\ntype = a\0b\n", 18));
    CHECK(nul.message == "t.deck:2: NUL byte");
    REQUIRE(nul.context.has_value());
    CHECK(nul.context->index == 2u);
}

TEST_CASE("deck syntax: reading a file", "[deck]") {
    const auto missing = deck::read("no/such/file.deck");
    REQUIRE_FALSE(missing.has_value());
    CHECK(missing.error().message.find("cannot open") != std::string::npos);
}

TEST_CASE("deck values: numbers, lengths, integers, switches, series", "[deck]") {
    using namespace deck::detail;
    CHECK(*number("1e17") == 1e17);
    CHECK(*number(" -2.5 ") == -2.5);
    CHECK(*number("+3") == 3.0);
    CHECK_FALSE(number("1e17cm").has_value());
    CHECK_FALSE(number("").has_value());
    CHECK_FALSE(number("inf").has_value());
    CHECK_FALSE(number("nan").has_value());
    CHECK_FALSE(number("1,5").has_value());
    CHECK(*length("2um") == 2e-4);
    CHECK(*length("50nm") == 50 * 1e-7);
    CHECK(*length("3e-4cm") == 3e-4);
    CHECK(*length("3e-4") == 3e-4);
    CHECK_FALSE(length("2 mm").has_value());
    CHECK(*integer("201") == 201);
    CHECK_FALSE(integer("2.0").has_value());
    CHECK(*boolean("on"));
    CHECK_FALSE(*boolean("off"));
    CHECK_FALSE(boolean("1").has_value());

    const auto num = [](std::string_view s) { return number(s); };
    const auto r = series("0 : 0.8 : 0.05", num);
    REQUIRE(r.has_value());
    REQUIRE(r->size() == 17);
    CHECK(r->front() == 0.0);
    CHECK(r->back() == 0.8);  // the last value is `last` exactly
    CHECK((*r)[3] == 3 * 0.05);
    const auto down = series("0:-1:-0.25", num);
    REQUIRE(down.has_value());
    CHECK(*down == std::vector<double>{0, -0.25, -0.5, -0.75, -1});
    const auto partial = series("0:1:0.3", num);  // the last value is not on the grid
    REQUIRE(partial.has_value());
    CHECK(partial->size() == 4);
    CHECK(partial->back() < 1.0);
    CHECK(*series("0.1, 0.5,2", num) == std::vector<double>{0.1, 0.5, 2});
    CHECK(*series("7", num) == std::vector<double>{7});
    CHECK_FALSE(series("0:1", num).has_value());
    CHECK_FALSE(series("0:1:0", num).has_value());
    CHECK_FALSE(series("0:1:-0.1", num).has_value());
    CHECK_FALSE(series("1,,2", num).has_value());
}
