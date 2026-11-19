// `nitcad compare` (compare.hpp; DECISIONS.md Unit 60, XC3).
#include "compare.hpp"

#include <charconv>
#include <cmath>
#include <cstdio>
#include <string>
#include <string_view>

#include "NiTCAD/analysis/compare.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/results/run.hpp"

namespace NiTCAD::cli {

namespace {

int usage(const char* why) {
    std::fprintf(stderr,
                 "nitcad compare: %s\nusage: nitcad compare <a.nitsweep> <b.nitsweep> [--rtol <x>] "
                 "[--potential-tol <V>]\n",
                 why);
    return 64;
}

bool parse_positive(std::string_view text, double& out) {
    double v = 0.0;
    const auto [end, ec] = std::from_chars(text.data(), text.data() + text.size(), v);
    if (ec != std::errc{} || end != text.data() + text.size() || !std::isfinite(v) || v <= 0.0) return false;
    out = v;
    return true;
}

// The machine and build a sweep came from, one line.
void print_origin(const char* label, const results::Sweep& s) {
    const results::Provenance& p = s.run.provenance;
    std::printf("%s: %s, revision %s%s, %s, %s\n", label, p.version.empty() ? "unknown" : p.version.c_str(),
                p.revision.empty() ? "unknown" : p.revision.c_str(), p.dirty ? " (modified)" : "",
                p.build_type.empty() ? "unknown" : p.build_type.c_str(), p.cpu.empty() ? "unknown CPU" : p.cpu.c_str());
}

}  // namespace

int compare(int argc, char** argv) {
    std::string a_path, b_path;
    analysis::CompareTolerances t;
    for (int i = 2; i < argc; ++i) {
        const std::string_view arg = argv[i];
        if (arg == "--rtol" || arg == "--potential-tol") {
            if (i + 1 >= argc) return usage("a value is missing");
            if (!parse_positive(argv[++i], arg == "--rtol" ? t.relative : t.potential_V))
                return usage("a tolerance must be a positive number");
        } else if (arg.starts_with("--")) {
            return usage("unknown option");
        } else if (a_path.empty()) {
            a_path = arg;
        } else if (b_path.empty()) {
            b_path = arg;
        } else {
            return usage("too many files");
        }
    }
    if (b_path.empty()) return usage("two sweep files are needed");

    const auto a = io::load_sweep(a_path);
    if (!a) {
        std::fprintf(stderr, "nitcad compare: %s: %s\n", a_path.c_str(), a.error().message.c_str());
        return 1;
    }
    const auto b = io::load_sweep(b_path);
    if (!b) {
        std::fprintf(stderr, "nitcad compare: %s: %s\n", b_path.c_str(), b.error().message.c_str());
        return 1;
    }

    print_origin("a", *a);
    print_origin("b", *b);
    const analysis::SweepComparison c = analysis::compare_sweeps(*a, *b, t);
    std::printf("input identity %s, mesh digest %s (reported, not compared)\n",
                c.identity_equal ? "equal" : "different", c.mesh_digest_equal ? "equal" : "different");
    if (c.outcome == analysis::CompareOutcome::not_comparable) {
        std::fprintf(stderr, "nitcad compare: not comparable: %s\n", c.reason.c_str());
        return 1;
    }
    for (const auto& q : c.quantities) {
        std::printf("%-24s max difference %.3e  ratio %.3e  at point %zu, index %zu  (%zu values)  %s\n",
                    q.name.c_str(), q.max_difference, q.max_ratio, q.point, q.index, q.compared,
                    q.pass() ? "ok" : "DIFFERS");
    }
    std::printf("%s\n", analysis::to_string(c.outcome).data());
    return c.outcome == analysis::CompareOutcome::match ? 0 : 4;
}

}  // namespace NiTCAD::cli
