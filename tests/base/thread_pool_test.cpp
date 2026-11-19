// base::ThreadPool (Unit 58 SD2): every index in exactly one call, the fixed chunks whatever the
// thread count, reuse across calls, work on more than one thread, and for_chunks without a pool.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <mutex>
#include <set>
#include <thread>
#include <utility>
#include <vector>

#include "NiTCAD/base/thread_pool.hpp"

using NiTCAD::base::ThreadPool;

TEST_CASE("thread pool: each index once, in the fixed chunks, for every thread count") {
    for (const int threads : {1, 2, 3, 6, 8}) {
        ThreadPool pool(threads);
        REQUIRE(pool.threads() == threads);
        for (const std::size_t chunk : {std::size_t{1}, std::size_t{64}, std::size_t{8192}}) {
            for (const std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{7},
                                        std::size_t{8192}, std::size_t{100003}}) {
                CAPTURE(threads, chunk, n);
                for (int repeat = 0; repeat < 2; ++repeat) {  // the pool is reused
                    std::vector<std::atomic<int>> hits(n);
                    std::mutex m;
                    std::set<std::pair<std::size_t, std::size_t>> calls;
                    pool.parallel_for(n, chunk, [&](std::size_t begin, std::size_t end) {
                        for (std::size_t i = begin; i < end; ++i) hits[i].fetch_add(1);
                        const std::scoped_lock lock(m);
                        calls.emplace(begin, end);
                    });
                    std::size_t wrong = 0;
                    for (const auto& h : hits) wrong += h.load() != 1;
                    REQUIRE(wrong == 0);
                    std::set<std::pair<std::size_t, std::size_t>> expected;
                    for (std::size_t b = 0; b < n; b += chunk) expected.emplace(b, std::min(n, b + chunk));
                    REQUIRE(calls == expected);
                }
            }
        }
    }
}

TEST_CASE("thread pool: the chunks run on more than one thread") {
    ThreadPool pool(4);
    std::mutex m;
    std::set<std::thread::id> ids;
    pool.parallel_for(64, 1, [&](std::size_t, std::size_t) {
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
        const std::scoped_lock lock(m);
        ids.insert(std::this_thread::get_id());
    });
    REQUIRE(ids.size() > 1);
}

TEST_CASE("thread pool: for_chunks without a pool is one call over the whole range") {
    std::vector<std::pair<std::size_t, std::size_t>> calls;
    NiTCAD::base::for_chunks(nullptr, 100, 8, [&](std::size_t b, std::size_t e) { calls.emplace_back(b, e); });
    REQUIRE(calls == std::vector<std::pair<std::size_t, std::size_t>>{{0, 100}});
    calls.clear();
    NiTCAD::base::for_chunks(nullptr, 0, 8, [&](std::size_t b, std::size_t e) { calls.emplace_back(b, e); });
    REQUIRE(calls.empty());
}
