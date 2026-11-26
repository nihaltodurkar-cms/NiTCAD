// libFuzzer target: the deck's syntax layer on arbitrary text (the release-fuzz preset).
//
// Any input must give a Deck or an error; a crash, an ASan report or a violated contract
// (NITCAD_EXPECTS terminates) is a bug. Covered: parse, parameter resolution and substitution,
// value splitting, unquoting, and the doping profile parser on every value.
#include <cstddef>
#include <cstdint>
#include <string_view>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    using namespace NiTCAD::deck;
    const std::string_view text(reinterpret_cast<const char*>(data), size);
    const auto deck = parse(text);
    if (!deck) return 0;
    (void)parameters(*deck);
    const auto sub = substituted(*deck);
    const Deck& d = sub ? *sub : *deck;
    for (const Section& s : d.sections) {
        for (const Entry& e : s.entries) {
            (void)value_words(e.value);
            (void)unquote(e.value);
            (void)parse_profile(e.value);
        }
    }
    return 0;
}
