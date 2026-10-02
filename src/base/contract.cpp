#include "NiTCAD/base/contract.hpp"

#include <cstdio>
#include <intrin.h>

namespace NiTCAD::base {

namespace {
// FAST_FAIL_FATAL_APP_EXIT from winnt.h (Windows SDK 10.0.26100.0, value 7); defined here so
// this file does not need windows.h.
constexpr unsigned int fast_fail_fatal_app_exit = 7;
}  // namespace

std::string format_contract_violation(std::string_view condition, std::string_view file, int line) {
    std::string text = "NiTCAD contract violation: ";
    text += condition;
    text += " (";
    text += file;
    text += ':';
    text += std::to_string(line);
    text += ")\n";
    return text;
}

void contract_violation(const char* condition, const char* file, int line) noexcept {
    // Best effort: building the message can in principle fail on allocation, and this function
    // must not throw. Fall back to the raw pieces if it does.
    try {
        const std::string text = format_contract_violation(condition, file, line);
        std::fputs(text.c_str(), stderr);
    } catch (...) {
        std::fputs("NiTCAD contract violation: ", stderr);
        std::fputs(condition, stderr);
        std::fputc('\n', stderr);
    }
    std::fflush(stderr);
    __fastfail(fast_fail_fatal_app_exit);
}

}  // namespace NiTCAD::base
