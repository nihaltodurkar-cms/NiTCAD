#include "NiTCAD/linalg/mkl_runtime.hpp"

#include <windows.h>

#include <atomic>
#include <mutex>
#include <string>
#include <system_error>
#include <utility>

#include "mkl_api.hpp"

namespace NiTCAD::linalg {

namespace {

constexpr int interface_lp64 = 0;    // MKL_INTERFACE_LP64
constexpr int threading_intel = 0;   // MKL_THREADING_INTEL

struct State {
    std::mutex mutex;
    HMODULE module = nullptr;
    detail::MklApi api{};
    std::optional<MklRuntime> runtime;
    std::atomic<const detail::MklApi*> published{nullptr};
};

State& state() {
    static State s;
    return s;
}

// The path as UTF-8 text for a message (path::string() throws for characters outside the code page).
std::string display(const std::filesystem::path& p) {
    const std::u8string u = p.u8string();
    return {reinterpret_cast<const char*>(u.data()), u.size()};
}

base::Error invalid(std::string message) {
    return {base::ErrorCode::invalid_input, std::move(message), std::nullopt};
}

template <class Function>
bool resolve(HMODULE module, const char* name, Function& out) {
    out = reinterpret_cast<Function>(reinterpret_cast<void*>(GetProcAddress(module, name)));
    return out != nullptr;
}

// A 1 x 1 system through every phase: the runtime and its threading layer really work.
bool check(const detail::MklApi& api) {
    void* pt[64] = {};
    detail::MklInt iparm[64] = {};
    const detail::MklInt mtype = 11, one = 1, quiet = 0, analyze_to_solve = 13, release = -1;
    api.pardisoinit(pt, &mtype, iparm);
    iparm[34] = 1;  // zero-based indices
    const detail::MklInt ia[2] = {0, 1}, ja[1] = {0};
    double a = 2.0, b = 4.0, x = 0.0;
    detail::MklInt perm = 0, error = 0;
    api.pardiso(pt, &one, &one, &mtype, &analyze_to_solve, &one, &a, ia, ja, &perm, &one, iparm,
                &quiet, &b, &x, &error);
    const bool solved = error == 0 && x == 2.0;
    detail::MklInt released = 0;
    api.pardiso(pt, &one, &one, &mtype, &release, &one, &a, ia, ja, &perm, &one, iparm, &quiet,
                &b, &x, &released);
    return solved;
}

}  // namespace

std::expected<MklRuntime, base::Error> load_mkl_runtime(const std::filesystem::path& path) {
    State& s = state();
    const std::lock_guard lock(s.mutex);
    if (s.runtime) return *s.runtime;

    std::error_code ec;
    const std::filesystem::path full = std::filesystem::absolute(path, ec);
    if (ec) return std::unexpected(invalid("the MKL runtime path is not valid"));
    // Its own libraries (mkl_core, the threading layer, libiomp5md) are found next to it.
    HMODULE module = LoadLibraryExW(full.c_str(), nullptr, LOAD_WITH_ALTERED_SEARCH_PATH);
    if (module == nullptr) {
        return std::unexpected(invalid("the MKL runtime could not be loaded: " + display(full)));
    }
    detail::MklApi api{};
    int (*set_interface)(int) = nullptr;
    int (*set_threading)(int) = nullptr;
    void (*version)(char*, int) = nullptr;
    const bool resolved = resolve(module, "pardisoinit", api.pardisoinit) &&
                          resolve(module, "pardiso", api.pardiso) &&
                          resolve(module, "MKL_Set_Num_Threads_Local", api.set_num_threads_local) &&
                          resolve(module, "MKL_Set_Interface_Layer", set_interface) &&
                          resolve(module, "MKL_Set_Threading_Layer", set_threading) &&
                          resolve(module, "MKL_Get_Version_String", version);
    if (!resolved) {
        FreeLibrary(module);
        return std::unexpected(invalid("an MKL entry point is missing: " + display(full)));
    }
    if (set_interface(interface_lp64) != interface_lp64 ||
        set_threading(threading_intel) != threading_intel) {
        FreeLibrary(module);
        return std::unexpected(invalid("the MKL interface or threading layer could not be set"));
    }
    if (!check(api)) {
        FreeLibrary(module);
        return std::unexpected(invalid("the MKL runtime failed its check solve"));
    }
    char text[256] = {};
    version(text, static_cast<int>(sizeof(text)) - 1);
    s.module = module;
    s.api = api;
    s.runtime = MklRuntime{.path = full, .version = text};
    s.published.store(&s.api, std::memory_order_release);
    return *s.runtime;
}

std::optional<MklRuntime> loaded_mkl_runtime() {
    State& s = state();
    const std::lock_guard lock(s.mutex);
    return s.runtime;
}

const detail::MklApi* detail::mkl_api() noexcept {
    return state().published.load(std::memory_order_acquire);
}

}  // namespace NiTCAD::linalg
