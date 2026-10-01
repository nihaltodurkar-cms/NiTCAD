// Unit 1 scaffold test: proves the toolchain, the pinned dependencies and the test
// harness work together. It is replaced in role by real layer tests as units land.
#include <catch2/catch_test_macros.hpp>

#include <expected>
#include <stop_token>

static_assert(_MSVC_LANG == 202302L, "NiTCAD requires /std:c++23preview");

TEST_CASE("scaffold: C++23 library features are available") {
    const std::expected<int, int> ok{42};
    REQUIRE(ok.has_value());
    REQUIRE(*ok == 42);

    const std::expected<int, int> bad = std::unexpected(7);
    REQUIRE_FALSE(bad.has_value());
    REQUIRE(bad.error() == 7);

    std::stop_source source;
    REQUIRE_FALSE(source.stop_requested());
    source.request_stop();
    REQUIRE(source.stop_requested());
}
