// [transient] restart and checkpoint (Unit 61, RD4; gates RG4 deck part, RG7, RG8): a transient deck
// stopped part-way writes a checkpoint with its resume state; the same deck with restart continues
// it, and the joined run equals the uninterrupted run bit for bit. The checkpoint round-trips bit
// for bit; a damaged file, another run's state, the final step, a steady checkpoint as a restart
// and a transient checkpoint as initial or as a sweep's restart are refused; `nitcad run` completes
// a restart (exit 0); a stopped run writes its checkpoint (the cancelled runs above). [trace]
// restart and checkpoint likewise (RD5; RG5 deck part).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stop_token>
#include <string>
#include <vector>

#include <windows.h>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/io/files.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::diode_device;
using deck_test::compensated_device;
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

// [trace] restart (Unit 61, RD5; RG5 deck part): a stopped trace deck continued equals the
// uninterrupted one bit for bit; its checkpoint is a trace's, still a steady point for a sweep.
TEST_CASE("deck trace restart: a stopped trace continued equals the uninterrupted trace (RG5)", "[deck][restart]") {
    const fs::path dir = scratch("trace_restart");
    const std::string trace = "[trace]\ncontact = anode\nend = 0.9\nstep = 0.05\nmax_step = 0.1\n";
    const auto run_trace = [&](const std::string& text, std::size_t stop_at) {
        const deck::Plan plan = must(text, dir);
        std::stop_source source;
        std::size_t count = 0;
        solve::RunControl control{.stop = source.get_token()};
        if (stop_at > 0) {
            control.progress = [&](const solve::Progress& p) {
                if (p.phase == solve::Phase::bias && p.point_count == 0 && p.converged && ++count == stop_at)
                    source.request_stop();
            };
        }
        auto r = deck::run(plan, control);
        if (!r) FAIL(r.error().message);
        REQUIRE(r->sweep.has_value());
        return std::move(*r->sweep);
    };
    const auto full = run_trace(diode_device + trace, 0);
    REQUIRE_FALSE(full.stopped.has_value());
    const auto mesh = must(diode_device + trace, dir).device.mesh();
    for (const std::size_t stop_at : {std::size_t{1}, std::size_t{5}}) {
        INFO("stopped at converged attempt " << stop_at);
        const auto first = run_trace(diode_device + trace + "[output]\ncheckpoint = part.nitckpt\n", stop_at);
        REQUIRE(first.stopped.has_value());
        const std::size_t k = first.points.size() - 1;
        const auto c = io::load_checkpoint(dir / "part.nitckpt", mesh);
        REQUIRE(c.has_value());
        CHECK(c->kind == io::CheckpointKind::trace);
        REQUIRE(c->trace.has_value());
        CHECK(same(c->trace->state, first.resume->state));
        CHECK(same(c->trace->tangent, first.resume->tangent));
        CHECK(same(c->trace->step, first.resume->step));
        CHECK(c->trace->entry == first.resume->entry);
        CHECK(c->trace->points == first.resume->points);
        CHECK(same(c->trace->point.terminal_current, first.resume->point.terminal_current));

        const auto rest = run_trace(diode_device + trace + "restart = part.nitckpt\n[output]\nsweep = rest.nitsweep\n", 0);
        REQUIRE_FALSE(rest.stopped.has_value());
        REQUIRE(rest.points.size() == full.points.size() - k);
        for (std::size_t i = 0; i < rest.points.size(); ++i) {
            INFO("point " << k + i);
            CHECK(same(rest.points[i].bias_V, full.points[k + i].bias_V));
            CHECK(same(rest.points[i].terminal_current, full.points[k + i].terminal_current));
            CHECK(same(rest.points[i].fields.potential_V, full.points[k + i].fields.potential_V));
        }
        CHECK(io::load_sweep(dir / "rest.nitsweep")->points.size() == rest.points.size());
    }
    // Refusals: a steady checkpoint, another trace. The trace checkpoint is still a steady point:
    // initial takes its fields.
    REQUIRE(deck::run(must(diode_device + "[sweep]\ncontact anode = 0.2\n[output]\ncheckpoint = steady.nitckpt\n", dir)));
    const auto refused = [&](const std::string& text, const std::string& what) {
        INFO(text);
        const auto r = deck::run(must(text, dir));
        REQUIRE_FALSE(r.has_value());
        INFO(r.error().message);
        CHECK(r.error().message.find(what) != std::string::npos);
    };
    refused(diode_device + trace + "restart = steady.nitckpt\n", "holds no trace's resume state");
    std::string other = trace;
    other.replace(other.find("end = 0.9"), 9, "end = 0.8");
    refused(diode_device + other + "restart = part.nitckpt\n", "belongs to another run");
    CHECK(deck::run(must(diode_device + "[sweep]\ncontact anode = 0.3\ninitial = part.nitckpt\n", dir)));
}

// [small_signal] restart (Unit 61, RD6; RG6): after the operating point a checkpoint holds (a
// sweep's or a small-signal run's), the remaining points equal the uninterrupted run's to 1e-9
// relative at every frequency; a checkpoint at the last point, a bias that is not a point and a
// transient's checkpoint are refused.
TEST_CASE("deck small-signal restart: the remaining points equal the uninterrupted run's (RG6)", "[deck][restart]") {
    const fs::path dir = scratch("small_signal_restart");
    const std::string ac = "[small_signal]\ncontact anode = 0 : 0.4 : 0.1\nfrequencies = 1MHz, 1GHz\n";
    const auto full = deck::run(must(compensated_device + ac, dir));
    if (!full) FAIL(full.error().message);
    const auto& f = *full->small_signal;
    REQUIRE(f.points.size() == 5);
    const auto close = [](double a, double b) { return std::abs(a - b) <= 1e-9 * std::max(std::abs(a), std::abs(b)) + 1e-30; };
    REQUIRE(deck::run(must(compensated_device + "[sweep]\ncontact anode = 0 : 0.2 : 0.1\n[output]\ncheckpoint = sweep.nitckpt\n", dir)));
    REQUIRE(deck::run(must(compensated_device + "[small_signal]\ncontact anode = 0 : 0.2 : 0.1\nfrequencies = 1MHz\n"
                           "[output]\ncheckpoint = ac.nitckpt\n", dir)));
    for (const std::string checkpoint : {"sweep.nitckpt", "ac.nitckpt"}) {
        INFO(checkpoint);
        const auto rest = deck::run(must(compensated_device + ac + "restart = " + checkpoint + "\n", dir));
        if (!rest) FAIL(rest.error().message);
        CHECK(rest->skipped_points == 3);
        const auto& r = *rest->small_signal;
        REQUIRE(r.points.size() == 2);
        for (std::size_t i = 0; i < 2; ++i) {
            const auto& a = r.points[i];
            const auto& b = f.points[3 + i];
            CHECK(a.dc.bias_V == b.dc.bias_V);
            for (std::size_t c = 0; c < a.dc.terminal_current.size(); ++c)
                CHECK(close(a.dc.terminal_current[c], b.dc.terminal_current[c]));
            REQUIRE(a.admittance.size() == b.admittance.size());
            for (std::size_t q = 0; q < a.admittance.size(); ++q) {
                REQUIRE(a.admittance[q].size() == b.admittance[q].size());
                for (std::size_t e = 0; e < a.admittance[q].size(); ++e) {
                    INFO("point " << 3 + i << ", frequency " << q << ", entry " << e);
                    CHECK(close(a.admittance[q][e].real(), b.admittance[q][e].real()));
                    CHECK(close(a.admittance[q][e].imag(), b.admittance[q][e].imag()));
                }
            }
        }
    }
    REQUIRE(deck::run(must(compensated_device + "[sweep]\ncontact anode = 0 : 0.4 : 0.1\n[output]\ncheckpoint = last.nitckpt\n", dir)));
    REQUIRE(deck::run(must(compensated_device + "[sweep]\ncontact anode = 0.05\n[output]\ncheckpoint = off.nitckpt\n", dir)));
    REQUIRE(deck::run(must(compensated_device + body + "[output]\ncheckpoint = tran.nitckpt\n", dir)));
    const auto refused = [&](const std::string& file, const std::string& what) {
        const auto r = deck::run(must(compensated_device + ac + "restart = " + file + "\n", dir));
        REQUIRE_FALSE(r.has_value());
        INFO(r.error().message);
        CHECK(r.error().message.find(what) != std::string::npos);
    };
    refused("last.nitckpt", "nothing is left to solve");
    refused("off.nitckpt", "not an operating point");
    refused("tran.nitckpt", "not an operating point");
}

namespace {

// Runs `nitcad run <deck>` in its own process group and sends it Ctrl-Break once its output holds
// `lines` progress lines with `marker`: the run stops as on a user's Ctrl-C. The exit status, or -1
// when the break could not be sent (no console to send it on) or the run ended first.
int run_and_break(const fs::path& deck, std::size_t lines, const std::string& marker) {
    SECURITY_ATTRIBUTES sa{sizeof sa, nullptr, TRUE};
    HANDLE read = nullptr, write = nullptr;
    REQUIRE(CreatePipe(&read, &write, &sa, 0));
    SetHandleInformation(read, HANDLE_FLAG_INHERIT, 0);
    STARTUPINFOW si{};
    si.cb = sizeof si;
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdOutput = write;
    si.hStdError = write;
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    std::wstring command = L"\"" + fs::path(NITCAD_CLI_PATH).wstring() + L"\" run \"" + deck.wstring() + L"\"";
    PROCESS_INFORMATION pi{};
    REQUIRE(CreateProcessW(nullptr, command.data(), nullptr, nullptr, TRUE, CREATE_NEW_PROCESS_GROUP, nullptr,
                           nullptr, &si, &pi));
    CloseHandle(write);
    std::string out;
    char buffer[256];
    DWORD got = 0;
    bool sent = false, delivered = false;
    while (ReadFile(read, buffer, sizeof buffer, &got, nullptr) && got > 0) {
        out.append(buffer, got);
        if (sent) continue;
        std::size_t count = 0;
        for (std::size_t at = out.find(marker); at != std::string::npos; at = out.find(marker, at + 1)) ++count;
        if (count >= lines) {
            sent = true;
            delivered = GenerateConsoleCtrlEvent(CTRL_BREAK_EVENT, pi.dwProcessId) != 0;
        }
    }
    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD status = 0;
    GetExitCodeProcess(pi.hProcess, &status);
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    CloseHandle(read);
    INFO(out);
    return delivered ? static_cast<int>(status) : -1;
}

void write_deck(const fs::path& p, const std::string& text) {
    std::ofstream f(p);
    f << text;
}

}  // namespace

// RG8: `nitcad run` on a transient and on a trace stopped by Ctrl-Break (exit 2) writes a
// checkpoint that a restart completes (exit 0).
TEST_CASE("deck restart: nitcad run stopped by Ctrl-Break writes a checkpoint a restart completes (RG8)",
          "[deck][restart]") {
    const fs::path dir = scratch("restart_break");
    const auto mesh = must(diode_device + body, dir).device.mesh();
    const std::string trace = "[trace]\ncontact = anode\nend = 0.9\nstep = 0.01\nmax_step = 0.01\n";
    struct Run {
        std::string name, analysis, marker;
        io::CheckpointKind kind;
    };
    for (const Run& r : {Run{"transient", body, "t = ", io::CheckpointKind::transient},
                         Run{"trace", trace, "traced point", io::CheckpointKind::trace}}) {
        INFO(r.name);
        write_deck(dir / (r.name + ".deck"), diode_device + r.analysis + "[output]\ncheckpoint = " + r.name + ".nitckpt\n");
        const int status = run_and_break(dir / (r.name + ".deck"), 3, r.marker);
        if (status == -1) SKIP("no console to send Ctrl-Break on");
        CHECK(status == 2);
        const auto c = io::load_checkpoint(dir / (r.name + ".nitckpt"), mesh);
        REQUIRE(c.has_value());
        CHECK(c->kind == r.kind);
        write_deck(dir / (r.name + "_rest.deck"),
                   diode_device + r.analysis + "restart = " + r.name + ".nitckpt\n");
        const Cli rest = cli("run \"" + (dir / (r.name + "_rest.deck")).string() + "\"");
        INFO(rest.output);
        CHECK(rest.status == 0);
    }
}

TEST_CASE("deck transient restart: with density gradient, bit for bit (Unit 62, QV12)", "[deck][restart]") {
    const fs::path dir = scratch("transient_restart_dg");
    const std::string device = diode_device + "[models]\nquantum_confinement = density_gradient\n";
    const auto full = run_deck(device + body, dir);
    REQUIRE_FALSE(full.transient->stopped.has_value());
    const auto& f = *full.transient;
    REQUIRE_FALSE(f.snapshots.empty());
    CHECK(f.snapshots.back().fields.quantum_potential_n_eV.size() == 201);
    for (const std::size_t stop_after : {std::size_t{3}, f.points.size() / 2}) {
        INFO("stopped after " << stop_after << " steps of " << f.points.size() - 1);
        const auto first = run_deck(device + body + "[output]\ncheckpoint = part.nitckpt\n", dir, stop_after);
        REQUIRE(first.transient->stopped.has_value());
        const std::size_t k = first.transient->points.size() - 1;
        const auto c = io::load_checkpoint(dir / "part.nitckpt", must(device + body, dir).device.mesh());
        REQUIRE(c.has_value());
        REQUIRE(c->transient.has_value());
        CHECK(same(c->transient->fields.quantum_potential_n_eV, first.transient->resume->fields.quantum_potential_n_eV));
        const auto rest = run_deck(device + body + "restart = part.nitckpt\n", dir);
        REQUIRE_FALSE(rest.transient->stopped.has_value());
        const auto& r = *rest.transient;
        REQUIRE(r.points.size() == f.points.size() - k);
        for (std::size_t i = 0; i < r.points.size(); ++i) {
            INFO("point " << k + i);
            same_point(r.points[i], f.points[k + i]);
        }
        const auto& a = r.snapshots.back().fields;
        const auto& b = f.snapshots.back().fields;
        CHECK(same(a.potential_V, b.potential_V));
        CHECK(same(a.quantum_potential_n_eV, b.quantum_potential_n_eV));
        CHECK(same(a.quantum_potential_p_eV, b.quantum_potential_p_eV));
        // The model changed: another run's identity, refused.
        const auto plain = deck::run(must(diode_device + body + "restart = part.nitckpt\n", dir));
        REQUIRE_FALSE(plain.has_value());
    }
}
