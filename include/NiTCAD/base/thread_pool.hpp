// A persistent thread pool with a deterministic parallel_for (Unit 58 SD2; ARCHITECTURE.md 6.8).
//
// parallel_for(n, chunk, body) splits [0, n) into chunks of a fixed size, chunk, that the caller
// chooses independently of the thread count: chunk k is [k chunk, min(n, (k + 1) chunk)). body is
// called once per chunk, on the calling thread or a worker, and may write only outputs that belong
// to its own chunk's indices. Which thread runs a chunk is not fixed, but nothing a kernel computes
// depends on it: the chunks and the order of operations within each are the same for every thread
// count, so a kernel whose outputs are per index gives bit-identical results with 1 to N threads.
// A reduction over chunks would not be; the pool offers none (DECISIONS.md Unit 58, slice 1).
//
// A pool of one thread has no workers and runs every chunk on the caller, in chunk order. A pool is
// used by one caller at a time (it is not reentrant: body must not call parallel_for on the same
// pool). body must not throw; an exception escaping it terminates the process.
#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <vector>

namespace NiTCAD::base {

class ThreadPool {
public:
    // threads >= 1 (NITCAD_EXPECTS): threads - 1 workers, the caller being the last.
    explicit ThreadPool(int threads);
    ~ThreadPool();
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    [[nodiscard]] int threads() const noexcept { return static_cast<int>(workers_.size()) + 1; }

    // chunk >= 1 (NITCAD_EXPECTS). Returns once every chunk has run. n == 0 calls nothing.
    template <class Body>
    void parallel_for(std::size_t n, std::size_t chunk, Body&& body) {
        using B = std::remove_reference_t<Body>;
        run(n, chunk, const_cast<std::remove_const_t<B>*>(std::addressof(body)),
            [](void* b, std::size_t begin, std::size_t end) noexcept {
                (*static_cast<B*>(b))(begin, end);
            });
    }

private:
    using Trampoline = void (*)(void*, std::size_t, std::size_t) noexcept;
    void run(std::size_t n, std::size_t chunk, void* body, Trampoline call);
    void work(const std::stop_token& stop);
    void run_chunks() noexcept;

    std::vector<std::jthread> workers_;
    // A worker spins on generation_ for a while after each job (kernels come in quick runs, as in
    // GMRES's Gram-Schmidt), then sleeps on start_; the caller spins until busy_ is 0.
    std::mutex mutex_;
    std::condition_variable_any start_;
    std::atomic<std::uint64_t> generation_{0};  // one per parallel_for
    std::atomic<std::size_t> busy_{0};          // workers still in the current generation
    std::atomic<std::size_t> sleeping_{0};      // workers on start_ (or about to be)
    // The current job, written before generation_ advances; read by workers only after they
    // see the new generation, and not rewritten until busy_ is 0.
    std::size_t n_ = 0, chunk_ = 1, chunks_ = 0;
    void* body_ = nullptr;
    Trampoline call_ = nullptr;
    std::atomic<std::size_t> next_{0};  // the next unclaimed chunk
};

// body over [0, n) in chunks: on the pool, or (no pool) on the caller as one call over [0, n). The
// same results either way for a per-index kernel.
template <class Body>
void for_chunks(ThreadPool* pool, std::size_t n, std::size_t chunk, Body&& body) {
    if (pool == nullptr || pool->threads() == 1 || n <= chunk) {
        if (n > 0) body(std::size_t{0}, n);
        return;
    }
    pool->parallel_for(n, chunk, body);
}

}  // namespace NiTCAD::base
