// The process's memory counters for the Unit 44 benchmark (solver_timing.hpp), in their own
// translation unit so windows.h stays out of the test files (NOMINMAX and WIN32_LEAN_AND_MEAN come
// from the build).
#include <windows.h>
#include <psapi.h>

#include "solver_timing.hpp"

ProcessMemory process_memory() {
    PROCESS_MEMORY_COUNTERS_EX c{};
    c.cb = sizeof c;
    GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&c),
                         sizeof c);
    constexpr double MB = 1024.0 * 1024.0;
    return {static_cast<double>(c.PrivateUsage) / MB, static_cast<double>(c.PeakWorkingSetSize) / MB};
}
