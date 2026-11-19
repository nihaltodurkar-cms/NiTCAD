// [transient] restart and checkpoint (Unit 61, RD4; gates RG4 deck part, RG7, RG8): a transient deck
// stopped part-way writes a checkpoint with its resume state; the same deck with restart continues
// it, and the joined run equals the uninterrupted run bit for bit. The checkpoint round-trips bit
// for bit; a damaged file, another run's state, the final step, a steady checkpoint as a restart
// and a transient checkpoint as initial or as a sweep's restart are refused; `nitcad run` completes
// a restart (exit 0); a stopped run writes its checkpoint (the cancelled runs above).
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stop_token>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/io/files.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_device;
using deck_test::must;
using deck_test::scratch;
namespace fs = std::filesystem;

namespace {

bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
bool same(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!same(a[i], b[i])) return false;
    return true;
}

void same_point(const results::TimePoint& a, const results::TimePoint& b) {
    CHECK(same(a.time_s, b.time_s));
    CHECK(same(a.step_s, b.step_s));
    CHECK(a.order == b.order);
    CHECK(same(a.error_ratio, b.error_ratio));
    CHECK(same(a.terminal_current, b.terminal_current));
    CHECK(same(a.displacement_current, b.displacement_current));
    CHECK(same(a.contact_charge, b.contact_charge));
}

const std::string body = R"([transient]
contact anode = pulse(0, 0.5, 100ps, 200ps, 50ps, 50ps)
t_end = 600ps
dt_initial = 0.1ps
output_times = 200ps, 400ps
)";

// Runs `text` in `dir`; with stop_after > 0 the run is cancelled right after that many accepted steps
// (the stop is asked on the accepting attempt's converged iteration, so the next loop top sees it
// with the step accepted; a step rejected after converging is not counted).
deck::RunReport run_deck(const std::string& text, const fs::path& dir, std::size_t stop_after = 0) {
    const deck::Plan plan = must(text, dir);
    std::stop_source source;
    solve::RunControl control{.stop = source.get_token()};
    std::size_t accepted_before = 0;
    double last_time = 0.0;
    if (stop_after > 0) {
        control.progress = [&](const solve::Progress& p) {
            if (p.phase != solve::Phase::transient || !p.converged) return;
            if (p.time_s > last_time) {
                last_time = p.time_s;
                if (++accepted_before == stop_after) source.request_stop();
            }
        };
    }
    auto r = deck::run(plan, control);
    if (!r) FAIL(r.error().message);
    REQUIRE(r->transient.has_value());
    return std::move(*r);
}

std::string text_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

struct Cli {
    int status;
    std::string output;
};
Cli cli(const std::string& arguments) {
    const std::string command = "\"\"" NITCAD_CLI_PATH "\" " + arguments + " 2>&1\"";
    FILE* pipe = _popen(command.c_str(), "r");
    REQUIRE(pipe != nullptr);
    std::string out;
    char buffer[512];
    while (std::fgets(buffer, sizeof buffer, pipe)) out += buffer;
    return {_pclose(pipe), out};
}

}  // namespace

TEST_CASE("deck transient restart: a stopped deck continued equals the uninterrupted run (RG4, RG7)",
          "[deck][restart]") {
    const fs::path dir = scratch("transient_restart");
    const auto full = run_deck(diode_device + body + "[output]\ntransient = full.nittran\n", dir);
    REQUIRE_FALSE(full.transient->stopped.has_value());
    const auto& f = *full.transient;

    for (const std::size_t stop_after : {std::size_t{3}, std::size_t{20}, f.points.size() / 2}) {
        INFO("stopped after " << stop_after << " steps of " << f.points.size() - 1);
        const auto first = run_deck(diode_device + body + "[output]\ncheckpoint = part.nitckpt\ntransient = part.nittran\n",
                                    dir, stop_after);
        REQUIRE(first.transient->stopped.has_value());
        const std::size_t k = first.transient->points.size() - 1;

        // RG7: the checkpoint holds the resume state bit for bit.
        const auto c = io::load_checkpoint(dir / "part.nitckpt", must(diode_device + body, dir).device.mesh());
        REQUIRE(c.has_value());
        CHECK(c->kind == io::CheckpointKind::transient);
        REQUIRE(c->transient.has_value());
        REQUIRE(first.transient->resume.has_value());
        const auto& a = *first.transient->resume;
        const auto& b = *c->transient;
        REQUIRE(a.history.size() == b.history.size());
        for (std::size_t i = 0; i < a.history.size(); ++i) {
            CHECK(same(a.history[i].t, b.history[i].t));
            CHECK(same(a.history[i].x, b.history[i].x));
            CHECK(same(a.history[i].storage, b.history[i].storage));
            CHECK(same(a.history[i].traps, b.history[i].traps));
            CHECK(same(a.history[i].charges, b.history[i].charges));
        }
        same_point(a.point, b.point);
        CHECK(same(a.time_s, b.time_s));
        CHECK(same(a.next_step_s, b.next_step_s));
        CHECK(a.landing == b.landing);
        CHECK(a.attempt == b.attempt);
        CHECK(a.identity == b.identity);
        CHECK(same(a.fields.potential_V, b.fields.potential_V));
        CHECK(c->bias_V == a.point.bias_V);

        const auto rest = run_deck(diode_device + body + "restart = part.nitckpt\n[output]\ntransient = rest.nittran\n", dir);
        REQUIRE_FALSE(rest.transient->stopped.has_value());
        const auto& r = *rest.transient;
        CHECK(r.run.input_identity == f.run.input_identity);
        REQUIRE(r.points.size() == f.points.size() - k);
        for (std::size_t i = 0; i < r.points.size(); ++i) {
            INFO("point " << k + i);
            same_point(r.points[i], f.points[k + i]);
        }
        // The file written holds the same (RQ8: its first point the checkpoint's step).
        const auto file = io::load_transient(dir / "rest.nittran");
        REQUIRE(file.has_value());
        REQUIRE(file->points.size() == r.points.size());
        same_point(file->points.front(), f.points[k]);
    }
}

TEST_CASE("deck transient restart: refusals (RG4, RG7)", "[deck][restart]") {
    const fs::path dir = scratch("transient_restart_refusals");
    const auto part = run_deck(diode_device + body + "[output]\ncheckpoint = part.nitckpt\n", dir, 5);
    REQUIRE(part.transient->stopped.has_value());
    const auto done = run_deck(diode_device + body + "[output]\ncheckpoint = done.nitckpt\n", dir);
    REQUIRE_FALSE(done.transient->stopped.has_value());
    REQUIRE(deck::run(must(diode_device + "[sweep]\ncontact anode = 0\n[output]\ncheckpoint = steady.nitckpt\n", dir)));

    const auto refused = [&](const std::string& text, const std::string& what) {
        INFO(text);
        const auto r = deck::run(must(text, dir));
        REQUIRE_FALSE(r.has_value());
        INFO(r.error().message);
        CHECK(r.error().message.find(what) != std::string::npos);
    };
    // Another run: a longer t_end (its step defaults differ), other waveforms.
    std::string longer = body;
    longer.replace(longer.find("t_end = 600ps"), 13, "t_end = 700ps");
    refused(diode_device + longer + "restart = part.nitckpt\n", "belongs to another run");
    std::string other = body;
    other.replace(other.find("pulse(0, 0.5"), 12, "pulse(0, 0.4");
    refused(diode_device + other + "restart = part.nitckpt\n", "belongs to another run");
    refused(diode_device + body + "restart = done.nitckpt\n", "final step");
    refused(diode_device + body + "restart = steady.nitckpt\n", "holds a steady state");
    refused(diode_device + body + "initial = part.nitckpt\n", "use restart to continue");
    refused(diode_device + "[sweep]\ncontact anode = 0 : 0.2 : 0.1\nrestart = part.nitckpt\n", "resume state");
    refused(diode_device + "[sweep]\ncontact anode = 0 : 0.2 : 0.1\ninitial = part.nitckpt\n", "use restart");

    // A damaged checkpoint (one byte of the resume record flipped) does not load.
    std::string bytes = text_of(dir / "part.nitckpt");
    bytes[bytes.size() - 20] = static_cast<char>(bytes[bytes.size() - 20] ^ 0x5a);
    {
        std::ofstream out(dir / "damaged.nitckpt", std::ios::binary);
        out << bytes;
    }
    const auto mesh = must(diode_device + body, dir).device.mesh();
    CHECK_FALSE(io::load_checkpoint(dir / "damaged.nitckpt", mesh).has_value());
    refused(diode_device + body + "restart = damaged.nitckpt\n", "");
}

TEST_CASE("deck transient restart: nitcad run completes a restart (RG8)",
          "[deck][restart]") {
    const fs::path dir = scratch("transient_restart_cli");
    const auto part = run_deck(diode_device + body + "[output]\ncheckpoint = part.nitckpt\n", dir, 7);
    REQUIRE(part.transient->stopped.has_value());
    {
        std::ofstream deck(dir / "rest.deck");
        deck << diode_device + body + "restart = part.nitckpt\n[output]\ntransient = rest.nittran\n";
    }
    const Cli rest = cli("run \"" + (dir / "rest.deck").string() + "\"");
    INFO(rest.output);
    CHECK(rest.status == 0);
    CHECK(io::load_transient(dir / "rest.nittran").has_value());

}
