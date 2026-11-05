// The verification report on disk (Unit 55, TVR-1; io/verify_report.hpp). Keys are local to this content kind.
#include "NiTCAD/io/verify_report.hpp"

#include <string>
#include <utility>

#include "container.hpp"
#include "record.hpp"

namespace NiTCAD::io {

namespace {

using detail::Reader;
using detail::Writer;

namespace key {
// Report.
constexpr std::uint16_t provenance = 1, threads = 2, tier = 3, evidence = 4, benchmark = 5;
// Provenance.
constexpr std::uint16_t version = 1, revision = 2, dirty = 3, build_type = 4, compiler = 5, eigen = 6,
                        cache = 7, mkl = 8, os = 9, cpu = 10, cpus = 11, started = 12, wall = 13;
// Benchmark.
constexpr std::uint16_t id = 1, title = 2, reference = 3, outcome = 4, run_error = 5, note = 6, comparison = 7;
// Comparison.
constexpr std::uint16_t quantity = 1, unit = 2, parameter = 3, q_sim = 4, q_ref = 5, u_ref = 6, u_num = 7,
                        metric = 8, criterion = 9, c_outcome = 10, reason = 11;
}  // namespace key

// A list of texts, each a one-field nested record (the record format lists records, not texts).
void write_texts(Writer& w, std::uint16_t k, const std::vector<std::string>& texts) {
    for (const std::string& t : texts) {
        Writer one;
        one.text(1, t);
        w.record(k, one);
    }
}

std::expected<std::vector<std::string>, base::Error> read_texts(const Reader& r, std::uint16_t k) {
    NITCAD_IO_TRY(list, r.records(k));
    std::vector<std::string> out;
    for (const Reader& one : list) {
        NITCAD_IO_TRY(t, one.text(1));
        out.push_back(std::move(t));
    }
    return out;
}

Writer write_provenance(const results::Provenance& p) {
    Writer w;
    w.text(key::version, p.version);
    w.text(key::revision, p.revision);
    w.u64(key::dirty, p.dirty ? 1 : 0);
    w.text(key::build_type, p.build_type);
    w.text(key::compiler, p.compiler);
    w.text(key::eigen, p.eigen_version);
    w.i64s(key::cache, p.eigen_cache_bytes);
    w.text(key::mkl, p.mkl_version);
    w.text(key::os, p.os);
    w.text(key::cpu, p.cpu);
    w.u64(key::cpus, p.logical_cpus);
    w.text(key::started, p.started_utc);
    w.f64(key::wall, p.wall_s);
    return w;
}

std::expected<results::Provenance, base::Error> read_provenance(const Reader& r) {
    results::Provenance p;
    NITCAD_IO_TRY(version, r.text(key::version));
    NITCAD_IO_TRY(revision, r.text(key::revision));
    NITCAD_IO_TRY(dirty, r.u64(key::dirty));
    NITCAD_IO_TRY(build_type, r.text(key::build_type));
    NITCAD_IO_TRY(compiler, r.text(key::compiler));
    NITCAD_IO_TRY(eigen, r.text(key::eigen));
    NITCAD_IO_TRY(cache, r.i64s(key::cache));
    NITCAD_IO_TRY(mkl, r.text(key::mkl));
    NITCAD_IO_TRY(os, r.text(key::os));
    NITCAD_IO_TRY(cpu, r.text(key::cpu));
    NITCAD_IO_TRY(cpus, r.u64(key::cpus));
    NITCAD_IO_TRY(started, r.text(key::started));
    NITCAD_IO_TRY(wall, r.f64(key::wall));
    p.version = std::move(version);
    p.revision = std::move(revision);
    p.dirty = dirty != 0;
    p.build_type = std::move(build_type);
    p.compiler = std::move(compiler);
    p.eigen_version = std::move(eigen);
    p.eigen_cache_bytes = std::move(cache);
    p.mkl_version = std::move(mkl);
    p.os = std::move(os);
    p.cpu = std::move(cpu);
    p.logical_cpus = static_cast<std::uint32_t>(cpus);
    p.started_utc = std::move(started);
    p.wall_s = wall;
    return p;
}

Writer write_comparison(const VerifyComparison& c) {
    Writer w;
    w.text(key::quantity, c.quantity);
    w.text(key::unit, c.unit);
    w.f64(key::parameter, c.parameter);
    w.f64(key::q_sim, c.q_sim);
    w.f64(key::q_ref, c.q_ref);
    w.f64(key::u_ref, c.u_ref);
    w.text(key::u_num, c.u_num);
    w.f64(key::metric, c.metric);
    w.text(key::criterion, c.criterion);
    w.text(key::c_outcome, c.outcome);
    w.text(key::reason, c.reason);
    return w;
}

std::expected<VerifyComparison, base::Error> read_comparison(const Reader& r) {
    VerifyComparison c;
    NITCAD_IO_TRY(quantity, r.text(key::quantity));
    NITCAD_IO_TRY(unit, r.text(key::unit));
    NITCAD_IO_TRY(parameter, r.f64(key::parameter));
    NITCAD_IO_TRY(q_sim, r.f64(key::q_sim));
    NITCAD_IO_TRY(q_ref, r.f64(key::q_ref));
    NITCAD_IO_TRY(u_ref, r.f64(key::u_ref));
    NITCAD_IO_TRY(u_num, r.text(key::u_num));
    NITCAD_IO_TRY(metric, r.f64(key::metric));
    NITCAD_IO_TRY(criterion, r.text(key::criterion));
    NITCAD_IO_TRY(outcome, r.text(key::c_outcome));
    NITCAD_IO_TRY(reason, r.text(key::reason));
    c.quantity = std::move(quantity);
    c.unit = std::move(unit);
    c.parameter = parameter;
    c.q_sim = q_sim;
    c.q_ref = q_ref;
    c.u_ref = u_ref;
    c.u_num = std::move(u_num);
    c.metric = metric;
    c.criterion = std::move(criterion);
    c.outcome = std::move(outcome);
    c.reason = std::move(reason);
    return c;
}

Writer write_benchmark(const VerifyBenchmark& b) {
    Writer w;
    w.text(key::id, b.id);
    w.text(key::title, b.title);
    w.text(key::reference, b.reference);
    w.text(key::outcome, b.outcome);
    w.text(key::run_error, b.run_error);
    write_texts(w, key::note, b.notes);
    for (const VerifyComparison& c : b.comparisons) w.record(key::comparison, write_comparison(c));
    return w;
}

std::expected<VerifyBenchmark, base::Error> read_benchmark(const Reader& r) {
    VerifyBenchmark b;
    NITCAD_IO_TRY(id, r.text(key::id));
    NITCAD_IO_TRY(title, r.text(key::title));
    NITCAD_IO_TRY(reference, r.text(key::reference));
    NITCAD_IO_TRY(outcome, r.text(key::outcome));
    NITCAD_IO_TRY(run_error, r.text(key::run_error));
    NITCAD_IO_TRY(comparisons, r.records(key::comparison));
    b.id = std::move(id);
    b.title = std::move(title);
    b.reference = std::move(reference);
    b.outcome = std::move(outcome);
    b.run_error = std::move(run_error);
    NITCAD_IO_TRY(note_texts, read_texts(r, key::note));
    b.notes = std::move(note_texts);
    for (const Reader& c : comparisons) {
        NITCAD_IO_TRY(comparison, read_comparison(c));
        b.comparisons.push_back(std::move(comparison));
    }
    return b;
}

}  // namespace

std::expected<void, base::Error> save_verify_report(const std::filesystem::path& path, const VerifyReport& report) {
    Writer w;
    w.record(key::provenance, write_provenance(report.provenance));
    w.u64(key::threads, report.threads);
    w.text(key::tier, report.tier);
    write_texts(w, key::evidence, report.reference_evidence);
    for (const VerifyBenchmark& b : report.benchmarks) w.record(key::benchmark, write_benchmark(b));
    return detail::write_file(path, Content::verify_report, w);
}

std::expected<VerifyReport, base::Error> load_verify_report(const std::filesystem::path& path) {
    NITCAD_IO_TRY(bytes, detail::read_file(path, Content::verify_report));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    VerifyReport report;
    NITCAD_IO_TRY(p, r.record(key::provenance));
    NITCAD_IO_TRY(provenance, read_provenance(p));
    report.provenance = std::move(provenance);
    NITCAD_IO_TRY(threads, r.u64(key::threads));
    report.threads = static_cast<std::uint32_t>(threads);
    NITCAD_IO_TRY(tier, r.text(key::tier));
    report.tier = std::move(tier);
    NITCAD_IO_TRY(evidence, read_texts(r, key::evidence));
    report.reference_evidence = std::move(evidence);
    NITCAD_IO_TRY(benchmarks, r.records(key::benchmark));
    for (const Reader& b : benchmarks) {
        NITCAD_IO_TRY(benchmark, read_benchmark(b));
        report.benchmarks.push_back(std::move(benchmark));
    }
    return report;
}

}  // namespace NiTCAD::io
