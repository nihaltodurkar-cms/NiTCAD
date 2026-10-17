// The assembly trace (Unit 40, DECISIONS.md FE7, FG1, FG2): with the environment variable
// NITCAD_TRACE_FILE set to a path, a test run writes there one line per value the assembler's
// public evaluations hand to the observer (assemble/observer.hpp), in call order under each test
// case and section: its index in the section, its tag, its length and a 64-bit FNV-1a hash of its
// bit patterns (-0.0 differs from 0.0; every NaN payload counts). Two runs of the same tests
// whose traces are equal saw bit-identical residuals, Jacobians, states and derived quantities;
// the first differing line names the test, the section and the call. Without the variable the
// listener does nothing. Linked into the test executables that assemble drift-diffusion systems.
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/catch_test_case_info.hpp>

#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <span>
#include <string>
#include <string_view>

#include "NiTCAD/assemble/observer.hpp"

namespace {

std::ofstream* trace = nullptr;
std::uint64_t call = 0;

std::uint64_t fnv1a(std::uint64_t h, std::uint64_t word) {
    for (int b = 0; b < 8; ++b) {
        h ^= (word >> (8 * b)) & 0xffu;
        h *= 0x100000001b3ull;
    }
    return h;
}

void record(std::string_view tag, std::span<const double> data) {
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const char c : tag) h = fnv1a(h, static_cast<unsigned char>(c));
    h = fnv1a(h, data.size());
    for (const double v : data) h = fnv1a(h, std::bit_cast<std::uint64_t>(v));
    char hex[17];
    std::snprintf(hex, sizeof hex, "%016llx", static_cast<unsigned long long>(h));
    *trace << call++ << '\t' << tag << '\t' << data.size() << '\t' << hex << '\n';
}

std::string environment(const char* name) {
    char* value = nullptr;
    std::size_t length = 0;
    if (_dupenv_s(&value, &length, name) != 0 || value == nullptr) return {};
    const std::unique_ptr<char, decltype(&std::free)> owned(value, &std::free);
    return value;
}

class TraceListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;

    void testRunStarting(const Catch::TestRunInfo&) override {
        const std::string path = environment("NITCAD_TRACE_FILE");
        if (path.empty()) return;
        file_.open(path, std::ios::binary | std::ios::trunc);
        if (!file_) {
            std::fprintf(stderr, "trace: cannot open %s\n", path.c_str());
            std::exit(2);
        }
        trace = &file_;
        NiTCAD::assemble::set_observer(&record);
    }
    void testCaseStarting(const Catch::TestCaseInfo& info) override {
        if (trace != nullptr) *trace << "case\t" << info.name << '\n';
    }
    void sectionStarting(const Catch::SectionInfo& info) override {
        if (trace == nullptr) return;
        *trace << "section\t" << info.name << '\n';
        call = 0;
    }
    void testRunEnded(const Catch::TestRunStats&) override {
        if (trace == nullptr) return;
        NiTCAD::assemble::set_observer(nullptr);
        file_.flush();
        trace = nullptr;
    }

private:
    std::ofstream file_;
};

}  // namespace

CATCH_REGISTER_LISTENER(TraceListener)
