// Error value for recoverable failures, carried as std::expected<T, Error>.
//
// Policy (ARCHITECTURE.md 6.7): recoverable errors are values; exceptions are reserved for
// exceptional conditions; a violated precondition is a bug and goes through NITCAD_EXPECTS
// (contract.hpp), not through Error. The std::string allocates only on the error path.
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace NiTCAD::base {

enum class ErrorCode : std::uint8_t {
    invalid_input,
    degenerate_mesh,
    singular_system,
    inaccurate_solve,
    non_convergence,
    cancelled,
    resource_exhausted,
};

// Optional numeric detail, so diagnostics need no string formatting inside kernels.
struct ErrorContext {
    std::optional<std::size_t> index;  // node, edge, row or bias step
    std::optional<double> value;       // residual or tolerance
};

struct Error {
    ErrorCode code;
    std::string message;
    std::optional<ErrorContext> context;
};

[[nodiscard]] constexpr std::string_view to_string(ErrorCode code) noexcept {
    switch (code) {
        case ErrorCode::invalid_input: return "invalid_input";
        case ErrorCode::degenerate_mesh: return "degenerate_mesh";
        case ErrorCode::singular_system: return "singular_system";
        case ErrorCode::inaccurate_solve: return "inaccurate_solve";
        case ErrorCode::non_convergence: return "non_convergence";
        case ErrorCode::cancelled: return "cancelled";
        case ErrorCode::resource_exhausted: return "resource_exhausted";
    }
    return "unknown";
}

}  // namespace NiTCAD::base
