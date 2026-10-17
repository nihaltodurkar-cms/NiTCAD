// The Eigen backend's results do not depend on the core a process runs on (Unit 41, DECISIONS.md).
// Eigen detects the cache sizes it blocks by once per process, so the check needs one process per
// core: run without arguments, the probe starts itself once per logical CPU of its affinity mask,
// each child pinned to its CPU before it runs. Each child solves a 22 500-unknown nonsymmetric
// five-point system with the default (Eigen SparseLU) LinearSolver and prints a 64-bit hash of the
// solution's bit patterns. Exit 0 when every child succeeds with the same hash. On a hybrid CPU
// without the fixed cache sizes (eigen_sparse_lu.cpp), the P- and E-core hashes differ.
#include <windows.h>

#include <bit>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace {

using namespace NiTCAD;

// Deterministic values in [-1, 1).
struct Noise {
    std::uint64_t s;
    double next() {
        std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        z ^= z >> 31;
        return static_cast<double>(z >> 11) * 0x1p-53 * 2.0 - 1.0;
    }
};

int child() {
    const int g = 150, n = g * g;
    Noise z{1};
    std::vector<linalg::Triplet> t;
    const auto entry = [&](int r, int c, double v) {
        t.push_back({static_cast<linalg::Index>(r), static_cast<linalg::Index>(c), v});
    };
    for (int y = 0; y < g; ++y) {
        for (int x = 0; x < g; ++x) {
            const int k = y * g + x;
            entry(k, k, 4.0 + 0.1 * z.next());
            if (x > 0) entry(k, k - 1, -1.0 + 0.3 * z.next());
            if (x + 1 < g) entry(k, k + 1, -1.0 + 0.3 * z.next());
            if (y > 0) entry(k, k - g, -1.0 + 0.3 * z.next());
            if (y + 1 < g) entry(k, k + g, -1.0 + 0.3 * z.next());
        }
    }
    const auto a = linalg::SparseMatrix::from_triplets(n, n, t);
    if (!a) return 2;
    std::vector<double> b(static_cast<std::size_t>(n)), x(b.size());
    for (double& v : b) v = z.next();
    auto solver = linalg::LinearSolver::create({});
    if (!solver || !solver->factorize(*a) || !solver->solve(b, x)) return 3;
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const double v : x) {
        const auto w = std::bit_cast<std::uint64_t>(v);
        for (int k = 0; k < 8; ++k) {
            h ^= (w >> (8 * k)) & 0xffu;
            h *= 0x100000001b3ull;
        }
    }
    std::printf("%016llx\n", static_cast<unsigned long long>(h));
    return 0;
}

// Runs this executable as a child pinned to `cpu`; its output, or "" if it failed.
std::string run_pinned(const std::wstring& self, int cpu) {
    SECURITY_ATTRIBUTES inherit{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
    HANDLE read = nullptr, write = nullptr;
    if (!CreatePipe(&read, &write, &inherit, 0)) return {};
    SetHandleInformation(read, HANDLE_FLAG_INHERIT, 0);
    STARTUPINFOW si{};
    si.cb = sizeof si;
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdOutput = write;
    si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    std::wstring command = L"\"" + self + L"\" child";
    PROCESS_INFORMATION pi{};
    const BOOL started = CreateProcessW(nullptr, command.data(), nullptr, nullptr, TRUE,
                                        CREATE_SUSPENDED, nullptr, nullptr, &si, &pi);
    CloseHandle(write);
    if (!started) {
        CloseHandle(read);
        return {};
    }
    std::string out;
    if (SetProcessAffinityMask(pi.hProcess, DWORD_PTR{1} << cpu) &&
        ResumeThread(pi.hThread) != static_cast<DWORD>(-1)) {
        char buffer[256];
        DWORD got = 0;
        while (ReadFile(read, buffer, sizeof buffer, &got, nullptr) && got > 0) {
            out.append(buffer, got);
        }
        WaitForSingleObject(pi.hProcess, INFINITE);
        DWORD code = 1;
        if (!GetExitCodeProcess(pi.hProcess, &code) || code != 0) out.clear();
    } else {
        TerminateProcess(pi.hProcess, 1);
    }
    CloseHandle(read);
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
    return out;
}

}  // namespace

int wmain(int argc, wchar_t** argv) {
    if (argc > 1 && std::wstring(argv[1]) == L"child") return child();
    wchar_t path[MAX_PATH];
    if (GetModuleFileNameW(nullptr, path, MAX_PATH) == 0) return 4;
    DWORD_PTR process = 0, system = 0;
    if (!GetProcessAffinityMask(GetCurrentProcess(), &process, &system)) return 4;
    std::string first;
    int cpus = 0, differ = 0;
    for (int cpu = 0; cpu < static_cast<int>(8 * sizeof(DWORD_PTR)); ++cpu) {
        if ((process >> cpu & 1) == 0) continue;
        const std::string hash = run_pinned(path, cpu);
        std::printf("cpu %d: %s\n", cpu, hash.empty() ? "failed" : hash.c_str());
        if (hash.empty()) return 5;
        if (cpus++ == 0) first = hash;
        if (hash != first) ++differ;
    }
    std::printf("%d CPUs, %d differ from the first one\n", cpus, differ);
    return differ == 0 ? 0 : 1;
}
