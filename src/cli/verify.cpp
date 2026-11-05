// `nitcad verify` (verify.hpp).
#include "verify.hpp"

#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/io/verify_report.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/solve/provenance.hpp"

namespace NiTCAD::cli {

namespace {

namespace bench = verify::benchmark;

int usage(std::span<const bench::Benchmark> registry, const std::string& why) {
    std::fprintf(stderr, "nitcad verify: %s\nusage: nitcad verify [--tier fast|full|nightly] [--id <id>] "
                         "[--report <file>]\n                     [--mkl <mkl_rt.dll>] [--threads <n>]\nknown ids:",
                 why.c_str());
    for (const auto& b : registry) std::fprintf(stderr, " %.*s", static_cast<int>(b.id.size()), b.id.data());
    std::fputs("\n", stderr);
    return 64;
}

void print(const bench::BenchmarkResult& r) {
    std::map<std::string, int> counts;
    for (const auto& c : r.comparisons) ++counts[std::string(verify::to_string(c.comparison.outcome))];
    const verify::Outcome o = bench::overall(r);
    std::string line = r.id + " " +
                       (o == verify::Outcome::not_applicable ? std::string("unverified")
                                                             : std::string(verify::to_string(o))) +
                       ": " + r.title + " (" + std::to_string(r.comparisons.size()) + " comparisons";
    for (const auto& [name, n] : counts) line += "; " + std::to_string(n) + " " + name;
    std::printf("%s)\n", line.c_str());
    if (r.run_error) std::printf("  run error: %s\n", r.run_error->c_str());
    for (const auto& n : r.notes) std::printf("  note: %s\n", n.c_str());
}

}  // namespace

int verify(int argc, char** argv, std::span<const bench::Benchmark> registry) {
    std::string tier = "fast", id, report, mkl, threads_text = "1";
    for (int i = 2; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (i + 1 >= argc || (a != "--tier" && a != "--id" && a != "--report" && a != "--mkl" &&
                            a != "--threads")) {
            return usage(registry, "unknown or incomplete option '" + std::string(a) + "'");
        }
        std::string& target = a == "--tier"     ? tier
                               : a == "--id"     ? id
                               : a == "--report" ? report
                               : a == "--mkl"    ? mkl
                                                 : threads_text;
        target = argv[++i];
    }
    if (tier != "fast" && tier != "full" && tier != "nightly") return usage(registry, "unknown tier '" + tier + "'");
    int threads = 0;
    const auto [end, ec] = std::from_chars(threads_text.data(), threads_text.data() + threads_text.size(), threads);
    if (ec != std::errc{} || end != threads_text.data() + threads_text.size() || threads < 1 || threads > 1024) {
        return usage(registry, "--threads must be 1 to 1024, not '" + threads_text + "'");
    }
    bench::set_pardiso_threads(threads);
    if (!mkl.empty()) {  // PARDISO for the solves that use it (TVB-04's 3D levels)
        if (const auto loaded = linalg::load_mkl_runtime(mkl); !loaded) {
            return usage(registry, "--mkl: " + loaded.error().message);
        }
    }
    const bench::Tier level =
        tier == "fast" ? bench::Tier::fast : tier == "full" ? bench::Tier::full : bench::Tier::nightly;
    std::vector<const bench::Benchmark*> selected;
    for (const auto& b : registry) {
        if (!id.empty() && b.id != id) continue;
        if (b.tier > level) {  // a benchmark runs in its own tier and those above it (CV-Q9)
            if (!id.empty()) return usage(registry, "'" + id + "' is not in tier " + tier);
            continue;
        }
        selected.push_back(&b);
    }
    if (selected.empty()) return usage(registry, "unknown id '" + id + "'");

    results::Provenance provenance = solve::provenance();
    const auto start = std::chrono::steady_clock::now();
    std::vector<bench::BenchmarkResult> results;
    for (const bench::Benchmark* b : selected) {
        bench::BenchmarkResult r;
        try {
            r = b->run(level);
        } catch (const std::exception& e) {
            r.id = std::string(b->id);
            r.title = std::string(b->title);
            r.run_error = std::string("exception: ") + e.what();
        }
        print(r);
        results.push_back(std::move(r));
    }
    provenance.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const bench::Summary s = bench::summarize(results);
    std::printf("%s\n", bench::summary_line(s, results.size()).c_str());

    int status = s.fail > 0 ? 1 : 0;
    for (const auto& r : results) {
        if (r.run_error) std::fprintf(stderr, "nitcad verify: %s could not run: %s\n", r.id.c_str(), r.run_error->c_str());
        if (bench::overall(r) == verify::Outcome::fail && !r.run_error) {
            std::fprintf(stderr, "nitcad verify: %s failed\n", r.id.c_str());
        }
    }
    if (!report.empty()) {
        // The deck solves use the default Eigen backend on one thread (RC-Q4); with --mkl, the PARDISO solves (TVB-04's
        // 3D levels) use --threads.
        const std::uint32_t used = linalg::loaded_mkl_runtime() ? static_cast<std::uint32_t>(threads) : 1;
        const io::VerifyReport rep = bench::make_report(results, provenance, used, tier);
        if (const auto saved = io::save_verify_report(report, rep); !saved) {
            std::fprintf(stderr, "nitcad verify: %s\n", saved.error().message.c_str());
            return 1;
        }
        std::ofstream csv(report + ".csv", std::ios::binary);
        csv << bench::report_csv(rep);
        csv.close();
        if (!csv) {
            std::fprintf(stderr, "nitcad verify: cannot write '%s.csv'\n", report.c_str());
            return 1;
        }
        std::printf("report: %s (and %s.csv)\n", report.c_str(), report.c_str());
    }
    return status;
}

}  // namespace NiTCAD::cli
