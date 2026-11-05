// The VR8 test build of `nitcad verify` (DECISIONS Unit 55 11.13 RC-Q7): the shipped registry plus the test-only
// entry TVB-TEST-FAIL, whose one comparison fails, so the exit status of a fail can be shown; it is compiled only into
// this executable, never into nitcad. Also `verify_probe compare <a> <b>`: 0 when two saved reports are equal with
// provenance.started_utc and wall_s cleared (VR7, RC-Q4), 1 otherwise or when one cannot be loaded.
#include <cstdio>
#include <string_view>
#include <vector>

#include "NiTCAD/io/verify_report.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "verify.hpp"

namespace bench = NiTCAD::verify::benchmark;
namespace v = NiTCAD::verify;

namespace {

bench::BenchmarkResult always_fails() {
    bench::BenchmarkResult r;
    r.id = "TVB-TEST-FAIL";
    r.title = "test-only: one comparison that fails";
    r.info.id = "TVB-TEST-FAIL";
    r.info.cls = v::ReferenceClass::E;
    r.info.independence = v::Independence::I3;
    r.info.provenance = v::Provenance::PA;
    r.info.citation = "test";
    r.info.version = "test";
    r.info.equations = {1};
    r.info.form = "1 = 2";
    r.info.derivation = "none";
    r.info.verification = "none";
    bench::ComparisonRecord c{"x", "1", 0.0, 1.0, 2.0, 0.0, 1.0, "q_sim == q_ref", {}};
    c.comparison = v::judge(r.info, {}, false, c.criterion);
    r.comparisons.push_back(c);
    return r;
}

int compare(const char* a, const char* b) {
    auto ra = NiTCAD::io::load_verify_report(a);
    auto rb = NiTCAD::io::load_verify_report(b);
    if (!ra || !rb) {
        std::fputs("verify_probe: a report could not be loaded\n", stderr);
        return 1;
    }
    for (auto* r : {&*ra, &*rb}) {
        r->provenance.started_utc.clear();
        r->provenance.wall_s = 0.0;
    }
    if (*ra == *rb) return 0;
    std::fputs("verify_probe: the reports differ\n", stderr);
    return 1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string_view(argv[1]) == "compare") return compare(argv[2], argv[3]);
    std::vector<bench::Benchmark> registry(bench::benchmarks().begin(), bench::benchmarks().end());
    registry.push_back({"TVB-TEST-FAIL", "test-only: one comparison that fails", bench::Tier::fast, &always_fails});
    if (argc >= 2 && std::string_view(argv[1]) == "verify") return NiTCAD::cli::verify(argc, argv, registry);
    std::fputs("usage: verify_probe verify ... | verify_probe compare <a> <b>\n", stderr);
    return 64;
}
