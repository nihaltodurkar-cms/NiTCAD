// The provenance of a run (Unit 42; DECISIONS.md PV1-PV3). The only file that includes the
// generated build_info.hpp, so a new revision recompiles this file alone.
#include "NiTCAD/solve/provenance.hpp"

#include <windows.h>

#include <intrin.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <format>
#include <string>

#include "NiTCAD/build_info.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"

namespace NiTCAD::solve {

namespace {

// RtlGetVersion reports the real version; GetVersionEx reports what the manifest declares.
std::string os_version() {
    using RtlGetVersion = LONG(WINAPI*)(OSVERSIONINFOW*);
    const HMODULE ntdll = GetModuleHandleW(L"ntdll.dll");
    if (ntdll == nullptr) return {};
    const auto get = reinterpret_cast<RtlGetVersion>(GetProcAddress(ntdll, "RtlGetVersion"));
    OSVERSIONINFOW v{};
    v.dwOSVersionInfoSize = sizeof v;
    if (get == nullptr || get(&v) != 0) return {};
    return std::format("Windows {}.{}.{}", v.dwMajorVersion, v.dwMinorVersion, v.dwBuildNumber);
}

// CPUID leaves 0x80000002-4: the processor's brand string, trailing blanks removed.
std::string cpu_brand() {
    std::array<int, 4> r{};
    __cpuid(r.data(), static_cast<int>(0x80000000u));
    if (static_cast<unsigned>(r[0]) < 0x80000004u) return {};
    char brand[49]{};
    for (unsigned leaf = 0; leaf < 3; ++leaf) {
        __cpuid(r.data(), static_cast<int>(0x80000002u + leaf));
        std::memcpy(brand + 16 * leaf, r.data(), 16);
    }
    std::string s(brand);
    while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) s.pop_back();
    return s;
}

}  // namespace

results::Provenance provenance() {
    results::Provenance p;
    p.version = build_info::version;
    p.revision = build_info::revision;
    p.dirty = build_info::dirty;
    p.build_type = build_info::build_type;
    p.compiler = build_info::compiler;
    const linalg::EigenBuild eigen = linalg::eigen_build();
    p.eigen_version = eigen.version;
    p.eigen_cache_bytes.assign(eigen.cache_bytes.begin(), eigen.cache_bytes.end());
    if (const auto mkl = linalg::loaded_mkl_runtime()) p.mkl_version = mkl->version;
    p.os = os_version();
    p.cpu = cpu_brand();
    p.logical_cpus = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
    const auto now = std::chrono::floor<std::chrono::milliseconds>(std::chrono::system_clock::now());
    p.started_utc = std::format("{:%FT%TZ}", now);
    return p;
}

}  // namespace NiTCAD::solve
