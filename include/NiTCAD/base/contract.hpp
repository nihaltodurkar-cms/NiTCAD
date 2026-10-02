// Precondition checking (ARCHITECTURE.md 6.7, Q2).
//
// NITCAD_EXPECTS(cond) is always on, in debug and release. A violated precondition is a bug,
// not a state to recover from, so on violation it reports file, line and condition to stderr
// and terminates the process immediately. It never throws.
#pragma once

#include <string>
#include <string_view>

namespace NiTCAD::base {

// The text written to stderr on violation. Separate so it can be tested without terminating.
[[nodiscard]] std::string format_contract_violation(std::string_view condition,
                                                    std::string_view file, int line);

// Reports and terminates. Not for direct use; call through NITCAD_EXPECTS.
[[noreturn]] void contract_violation(const char* condition, const char* file, int line) noexcept;

}  // namespace NiTCAD::base

#define NITCAD_EXPECTS(cond)                                                          \
    do {                                                                              \
        if (!(cond)) [[unlikely]] {                                                   \
            ::NiTCAD::base::contract_violation(#cond, __FILE__, __LINE__);            \
        }                                                                             \
    } while (false)
