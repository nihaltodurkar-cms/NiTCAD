#include <catch2/catch_test_macros.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "NiTCAD/base/contract.hpp"

namespace {
// Runs the probe executable (path injected by CMake) with stdout discarded and stderr captured.
struct ProbeResult {
    int exit_code;
    std::string stderr_text;
};

ProbeResult run_probe(const std::string& mode) {
    const std::filesystem::path capture =
        std::filesystem::temp_directory_path() / ("nitcad_contract_probe_" + mode + ".txt");
    // The outer quotes are required by cmd.exe when the command itself contains quotes.
    const std::string command = "\"\"" NITCAD_CONTRACT_PROBE_PATH "\" " + mode + " >NUL 2>\"" +
                                capture.string() + "\"\"";
    const int code = std::system(command.c_str());
    std::string text;
    {
        std::ifstream in(capture, std::ios::binary);
        text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }  // closed before the file is removed
    std::filesystem::remove(capture);
    return {code, text};
}
}  // namespace

TEST_CASE("contract: the violation text names the condition, file and line") {
    const std::string text = NiTCAD::base::format_contract_violation("x > 0", "a/b.cpp", 42);
    REQUIRE(text == "NiTCAD contract violation: x > 0 (a/b.cpp:42)\n");
}

TEST_CASE("contract: a satisfied precondition has no effect and evaluates its condition once") {
    int calls = 0;
    NITCAD_EXPECTS(++calls == 1);
    REQUIRE(calls == 1);
    NITCAD_EXPECTS(true);
}

TEST_CASE("contract: a satisfied precondition lets the probe process exit normally") {
    const ProbeResult result = run_probe("pass");
    REQUIRE(result.exit_code == 0);
}

TEST_CASE("contract: a violated precondition terminates the process and says what failed") {
    const ProbeResult result = run_probe("violate");
    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("NiTCAD contract violation: violate ? argc == 99 : argc >= 1") !=
            std::string::npos);
    REQUIRE(result.stderr_text.find("contract_probe.cpp") != std::string::npos);
}
