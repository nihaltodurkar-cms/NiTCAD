// The io layer's record encoding (Unit 29, src/io/record.hpp): exact round trips, little-endian
// layout, unknown keys skipped, absent keys defaulted, malformed bytes refused.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "record.hpp"

using namespace NiTCAD;
using io::detail::Reader;
using io::detail::Writer;

TEST_CASE("record: every value type reads back bit for bit") {
    const std::vector<double> reals{0.0,
                                    -0.0,
                                    1.0 / 3.0,
                                    -1e-308,
                                    std::numeric_limits<double>::denorm_min(),
                                    std::numeric_limits<double>::max(),
                                    std::numeric_limits<double>::infinity(),
                                    std::numeric_limits<double>::quiet_NaN()};
    const std::vector<std::int64_t> integers{0, -1, std::numeric_limits<std::int64_t>::min(),
                                             std::numeric_limits<std::int64_t>::max()};
    Writer nested;
    nested.text(1, "inner");
    Writer w;
    w.u64(1, 0xfedcba9876543210ULL);
    w.i64(2, -42);
    w.f64(3, -0.0);
    w.text(4, std::string_view("anode\0x", 7));
    w.f64s(5, reals);
    w.i64s(6, integers);
    w.record(7, nested);
    w.f64s(8, {});
    auto r = Reader::parse(w.bytes());
    REQUIRE(r.has_value());
    CHECK(*r->u64(1) == 0xfedcba9876543210ULL);
    CHECK(*r->i64(2) == -42);
    CHECK(std::bit_cast<std::uint64_t>(*r->f64(3)) == std::bit_cast<std::uint64_t>(-0.0));
    CHECK(*r->text(4) == std::string("anode\0x", 7));
    const auto back = *r->f64s(5);
    REQUIRE(back.size() == reals.size());
    for (std::size_t i = 0; i < reals.size(); ++i) {
        CHECK(std::bit_cast<std::uint64_t>(back[i]) == std::bit_cast<std::uint64_t>(reals[i]));
    }
    CHECK(*r->i64s(6) == integers);
    CHECK(*r->record(7)->text(1) == "inner");
    CHECK(r->f64s(8)->empty());
}

TEST_CASE("record: the layout is little-endian whatever the host") {
    Writer w;
    w.u64(0x0102, 0x0807060504030201ULL);
    const auto& b = w.bytes();
    REQUIRE(b.size() == 2 + 8 + 8);
    CHECK(std::to_integer<int>(b[0]) == 0x02);
    CHECK(std::to_integer<int>(b[1]) == 0x01);
    CHECK(std::to_integer<int>(b[2]) == 8);  // the length, low byte first
    for (int i = 0; i < 8; ++i) CHECK(std::to_integer<int>(b[10 + i]) == i + 1);
}

TEST_CASE("record: unknown keys are skipped, absent keys take defaults, lists keep order") {
    Writer w;
    w.u64(1, 7);
    w.text(999, "from a newer writer");
    for (int k = 0; k < 3; ++k) {
        Writer item;
        item.i64(1, k);
        w.record(2, item);
    }
    auto r = Reader::parse(w.bytes());
    REQUIRE(r.has_value());
    CHECK(*r->u64(1) == 7);
    CHECK(*r->u64_or(3, 11) == 11);
    CHECK(*r->i64_or(4, -3) == -3);
    CHECK(r->f64s_or_empty(5)->empty());
    const auto items = *r->records(2);
    REQUIRE(items.size() == 3);
    for (int k = 0; k < 3; ++k) CHECK(*items[k].i64(1) == k);
    // An absent field without a default is an error, not a silent zero.
    CHECK(r->u64(3).error().code == base::ErrorCode::invalid_input);
}

TEST_CASE("record: malformed bytes are refused") {
    Writer w;
    w.f64s(1, std::vector<double>{1.0, 2.0});
    const auto& b = w.bytes();
    // Truncated in the header and in the payload.
    CHECK_FALSE(Reader::parse(std::span(b).first(5)).has_value());
    CHECK_FALSE(Reader::parse(std::span(b).first(b.size() - 1)).has_value());
    // A field of the wrong size for its type.
    Writer odd;
    odd.text(1, "abc");
    auto r = Reader::parse(odd.bytes());
    REQUIRE(r.has_value());
    CHECK_FALSE(r->u64(1).has_value());
    CHECK_FALSE(r->f64s(1).has_value());
}
