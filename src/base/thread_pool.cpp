// ThreadPool (thread_pool.hpp).
#include "NiTCAD/base/thread_pool.hpp"

#include <chrono>
#include <immintrin.h>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::base {

namespace {

// How long a worker spins for the next job before it sleeps (Unit 58 slice 1: GMRES's kernels come
// a fraction of a millisecond apart; a sleeping worker takes tens of microseconds to wake).
constexpr auto spin_window = std::chrono::milliseconds(2);

}  // namespace

ThreadPool::ThreadPool(int threads) {
    NITCAD_EXPECTS(threads >= 1);
    workers_.reserve(static_cast<std::size_t>(threads - 1));
    for (int i = 1; i < threads; ++i) {
        workers_.emplace_back([this](const std::stop_token& stop) { work(stop); });
    }
}

ThreadPool::~ThreadPool() {
    for (std::jthread& w : workers_) w.request_stop();
    {
        const std::scoped_lock lock(mutex_);
    }
    start_.notify_all();
    workers_.clear();  // joins
}

void ThreadPool::run_chunks() noexcept {
    for (;;) {
        const std::size_t k = next_.fetch_add(1, std::memory_order_relaxed);
        if (k >= chunks_) return;
        const std::size_t begin = k * chunk_;
        const std::size_t end = begin + chunk_ < n_ ? begin + chunk_ : n_;
        call_(body_, begin, end);
    }
}

void ThreadPool::run(std::size_t n, std::size_t chunk, void* body, Trampoline call) {
    NITCAD_EXPECTS(chunk >= 1);
    if (n == 0) return;
    const std::size_t chunks = (n + chunk - 1) / chunk;
    if (workers_.empty()) {
        for (std::size_t k = 0; k < chunks; ++k) {
            const std::size_t begin = k * chunk;
            call(body, begin, begin + chunk < n ? begin + chunk : n);
        }
        return;
    }
    // busy_ is 0 here: the previous job's workers have all left it.
    n_ = n;
    chunk_ = chunk;
    chunks_ = chunks;
    body_ = body;
    call_ = call;
    next_.store(0, std::memory_order_relaxed);
    busy_.store(workers_.size(), std::memory_order_relaxed);
    generation_.fetch_add(1);  // seq_cst: publishes the job, ordered before sleeping_'s load
    if (sleeping_.load() > 0) {
        { const std::scoped_lock lock(mutex_); }  // a worker between its check and its wait
        start_.notify_all();
    }
    run_chunks();
    while (busy_.load(std::memory_order_acquire) != 0) _mm_pause();
}

void ThreadPool::work(const std::stop_token& stop) {
    std::uint64_t seen = 0;
    for (;;) {
        // Spin for the next job, then sleep.
        const auto until = std::chrono::steady_clock::now() + spin_window;
        bool started = false;
        for (unsigned k = 1; !stop.stop_requested(); ++k) {
            if (generation_.load(std::memory_order_acquire) != seen) {
                started = true;
                break;
            }
            _mm_pause();
            if ((k & 1023u) == 0 && std::chrono::steady_clock::now() > until) break;
        }
        if (!started) {
            std::unique_lock lock(mutex_);
            sleeping_.fetch_add(1);  // seq_cst, before the predicate's load
            const bool woke = start_.wait(lock, stop, [&] { return generation_.load() != seen; });
            sleeping_.fetch_sub(1);
            if (!woke) return;  // stopped
        }
        seen = generation_.load(std::memory_order_acquire);
        run_chunks();
        busy_.fetch_sub(1, std::memory_order_release);
    }
}

}  // namespace NiTCAD::base
