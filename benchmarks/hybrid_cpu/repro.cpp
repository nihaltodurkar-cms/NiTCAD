// Unit 41 (DECISIONS.md), minimal reproducer: on a hybrid CPU, Eigen's results depend in the last
// bits on the core that first queries the cache sizes. Prints the executing core's CPUID leaf-4
// data-cache sizes, the sizes Eigen adopted, and FNV-1a hashes of a dense product and a SparseLU
// solve. Standalone (Eigen only; not part of the CMake build): see README.md.
//   repro auto                      Eigen queries the cache sizes itself (as NiTCAD does today)
//   repro fixed [l1 l2 l3 [n]]      Eigen::setCpuCacheSizes(l1, l2, l3) first (bytes; default
//                                   48 KiB, 2 MiB, 8 MiB); n the dense product size (300)
#include <intrin.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

std::uint64_t fnv(const double* v, std::size_t n) {
    std::uint64_t h = 1469598103934665603ull;
    for (std::size_t i = 0; i < n; ++i) {
        const auto b = std::bit_cast<std::uint64_t>(v[i]);
        for (int k = 0; k < 8; ++k) {
            h ^= (b >> (8 * k)) & 0xFF;
            h *= 1099511628211ull;
        }
    }
    return h;
}

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

// The executing core's data/unified cache sizes by level, from CPUID leaf 4.
void core_caches(long long out[4]) {
    for (int i = 0; i < 4; ++i) out[i] = 0;
    for (int id = 0; id < 16; ++id) {
        int r[4];
        __cpuidex(r, 4, id);
        const int type = r[0] & 0x1F;
        if (type == 0) break;
        if (type != 1 && type != 3) continue;
        const int level = (r[0] >> 5) & 0x7;
        const long long ways = ((r[1] >> 22) & 0x3FF) + 1, parts = ((r[1] >> 12) & 0x3FF) + 1,
                        line = (r[1] & 0xFFF) + 1, sets = static_cast<long long>(r[2]) + 1;
        if (level < 4) out[level] = ways * parts * line * sets;
    }
}

}  // namespace

int main(int argc, char** argv) {
    const bool fixed = argc > 1 && std::strcmp(argv[1], "fixed") == 0;
    long long c[4];
    core_caches(c);
    if (fixed) {
        const std::ptrdiff_t l1 = argc > 2 ? std::atoll(argv[2]) : 48 * 1024;
        const std::ptrdiff_t l2 = argc > 3 ? std::atoll(argv[3]) : 2 * 1024 * 1024;
        const std::ptrdiff_t l3 = argc > 4 ? std::atoll(argv[4]) : 8 * 1024 * 1024;
        Eigen::setCpuCacheSizes(l1, l2, l3);
    }
    // Dense: C = A B, n x n (blocked by the cache sizes).
    const int n = argc > 5 ? std::atoi(argv[5]) : 300;
    Eigen::MatrixXd A(n, n), B(n, n);
    Noise z{1};
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            A(i, j) = z.next();
            B(i, j) = z.next();
        }
    }
    const Eigen::MatrixXd C = A * B;
    // Sparse: a nonsymmetric 2D five-point operator on 150 x 150 nodes, SparseLU.
    const int g = 150, N = g * g;
    std::vector<Eigen::Triplet<double>> t;
    for (int y = 0; y < g; ++y) {
        for (int x = 0; x < g; ++x) {
            const int k = y * g + x;
            t.emplace_back(k, k, 4.0 + 0.1 * z.next());
            if (x > 0) t.emplace_back(k, k - 1, -1.0 + 0.3 * z.next());
            if (x + 1 < g) t.emplace_back(k, k + 1, -1.0 + 0.3 * z.next());
            if (y > 0) t.emplace_back(k, k - g, -1.0 + 0.3 * z.next());
            if (y + 1 < g) t.emplace_back(k, k + g, -1.0 + 0.3 * z.next());
        }
    }
    Eigen::SparseMatrix<double> M(N, N);
    M.setFromTriplets(t.begin(), t.end());
    M.makeCompressed();
    Eigen::VectorXd b(N);
    for (int k = 0; k < N; ++k) b(k) = z.next();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.compute(M);
    const Eigen::VectorXd x = lu.solve(b);
    std::printf("core L1d %lld L2 %lld L3 %lld | eigen l1 %td l2 %td l3 %td | dense %016llx"
                " | sparse %016llx\n",
                c[1], c[2], c[3], Eigen::l1CacheSize(), Eigen::l2CacheSize(),
                Eigen::l3CacheSize(), static_cast<unsigned long long>(fnv(C.data(), C.size())),
                static_cast<unsigned long long>(fnv(x.data(), x.size())));
    return 0;
}
