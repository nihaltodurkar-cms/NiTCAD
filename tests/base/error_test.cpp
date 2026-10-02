#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <expected>

#include "NiTCAD/base/error.hpp"

using namespace NiTCAD::base;

TEST_CASE("error: every code has a distinct name") {
    static_assert(sizeof(ErrorCode) == 1);
    REQUIRE(to_string(ErrorCode::invalid_input) == "invalid_input");
    REQUIRE(to_string(ErrorCode::degenerate_mesh) == "degenerate_mesh");
    REQUIRE(to_string(ErrorCode::singular_system) == "singular_system");
    REQUIRE(to_string(ErrorCode::inaccurate_solve) == "inaccurate_solve");
    REQUIRE(to_string(ErrorCode::non_convergence) == "non_convergence");
    REQUIRE(to_string(ErrorCode::cancelled) == "cancelled");
    REQUIRE(to_string(ErrorCode::resource_exhausted) == "resource_exhausted");
}

TEST_CASE("error: an out-of-range code is reported as unknown, not undefined behaviour") {
    REQUIRE(to_string(static_cast<ErrorCode>(std::uint8_t{200})) == "unknown");
}

TEST_CASE("error: carried as std::expected with optional numeric context") {
    const std::expected<double, Error> failed = std::unexpected(Error{
        ErrorCode::non_convergence, "Newton did not converge",
        ErrorContext{.index = 12, .value = 3.5e-4}});

    REQUIRE_FALSE(failed.has_value());
    REQUIRE(failed.error().code == ErrorCode::non_convergence);
    REQUIRE(failed.error().message == "Newton did not converge");
    REQUIRE(failed.error().context.has_value());
    REQUIRE(failed.error().context->index == 12);
    REQUIRE(failed.error().context->value == 3.5e-4);

    const std::expected<double, Error> ok = 1.5;
    REQUIRE(ok.has_value());
}

TEST_CASE("error: context is optional and its fields are independently optional") {
    const Error bare{ErrorCode::invalid_input, "negative length", std::nullopt};
    REQUIRE_FALSE(bare.context.has_value());

    const Error partial{ErrorCode::singular_system, "zero pivot", ErrorContext{.index = 3, .value = std::nullopt}};
    REQUIRE(partial.context->index == 3);
    REQUIRE_FALSE(partial.context->value.has_value());
}
