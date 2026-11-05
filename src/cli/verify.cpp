// `nitcad verify` (verify.hpp).
#include "verify.hpp"

#include <chrono>
#include <cstdio>
#include <exception>
#include <fstream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/io/verify_report.hpp"
#include "NiTCAD/solve/provenance.hpp"

namespace NiTCAD::cli {

namespace {

namespace bench = verify::benchmark;

int usage(std::span<const bench::Benchmark> registry, const std::string& why) {
    std::fprintf(stderr, "nitcad verify: %s\nusage: nitcad verify [--tier fast|full|nightly] [--id <id>] "
                         "[--report <file>]\nknown ids:", why.c_str());
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
    std::string tier = "fast", id, report;
    for (int i = 2; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (i + 1 >= argc || (a != "--tier" && a != "--id" && a != "--report")) {
            return usage(registry, "unknown or incomplete option '" + std::string(a) + "'");
        }
        std::string& target = a == "--tier" ? tier : a == "--id" ? id : report;
        target = argv[++i];
    }
    if (tier != "fast" && tier != "full" && tier != "nightly") return usage(registry, "unknown tier '" + tier + "'");
    std::vector<const bench::Benchmark*> selected;
    for (const auto& b : registry) {
        if (id.empty() || b.id == id) selected.push_back(&b);  // every TVR-1 benchmark is in every tier
    }
    if (selected.empty()) return usage(registry, "unknown id '" + id + "'");

    results::Provenance provenance = solve::provenance();
    const auto start = std::chrono::steady_clock::now();
    std::vector<bench::BenchmarkResult> results;
    for (const bench::Benchmark* b : selected) {
        bench::BenchmarkResult r;
        try {
            r = b->run();
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
        // The benchmarks run on one thread: the deck solves use the default Eigen backend, not PARDISO (RC-Q4).
        const io::VerifyReport rep = bench::make_report(results, provenance, 1, tier);
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
