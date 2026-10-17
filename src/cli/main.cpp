// nitcad: the batch runner (ARCHITECTURE.md 11.3, Unit 31; DECISIONS.md DK10).
//
//   nitcad check <deck>   parse and build the deck, print what it describes, solve nothing
//   nitcad run <deck>     build the deck, run its analysis, write its outputs, and print the run's
//                         provenance (Unit 42)
//   nitcad --version      print the build: version, git revision, build type, compiler, Eigen
//
// Exit status: 0 done; 1 the deck or the run failed (the message on stderr); 2 the sweep stopped
// early (its completed points and outputs are written; the reason on stderr); 64 usage.
// Ctrl+C or Ctrl+Break asks the run to stop: the solver stops at its next check and the completed
// points are written, as for any early stop.
#include <windows.h>

#include <cstdio>
#include <stop_token>
#include <string>
#include <string_view>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/solve/provenance.hpp"

namespace {

std::stop_source stop;

BOOL WINAPI on_console_event(DWORD event) {
    if (event == CTRL_C_EVENT || event == CTRL_BREAK_EVENT) {
        stop.request_stop();
        return TRUE;
    }
    return FALSE;
}

int usage() {
    std::fputs("usage: nitcad check <deck>\n       nitcad run <deck>\n       nitcad --version\n",
               stderr);
    return 64;
}

// The build, one line: "nitcad 0.0.0 (revision <commit>[, modified], Release, MSVC ..., Eigen ...)".
void print_build(const NiTCAD::results::Provenance& p) {
    std::printf("nitcad %s (revision %s%s, %s, %s, Eigen %s)\n", p.version.c_str(),
                p.revision.c_str(), p.dirty ? ", modified" : "", p.build_type.c_str(),
                p.compiler.c_str(), p.eigen_version.c_str());
}

// A run's provenance (Unit 42): the build, the platform, the start and the wall time.
void print_provenance(const NiTCAD::results::Provenance& p) {
    std::fputs("build: ", stdout);
    print_build(p);
    std::printf("platform: %s, %s, %u logical CPUs%s%s\n", p.os.c_str(), p.cpu.c_str(),
                p.logical_cpus, p.mkl_version.empty() ? "" : ", MKL ", p.mkl_version.c_str());
    std::printf("started: %s, wall time %.3f s\n", p.started_utc.c_str(), p.wall_s);
}

int fail(const NiTCAD::base::Error& e) {
    std::fprintf(stderr, "nitcad: %s\n", e.message.c_str());
    return 1;
}

const char* kind_name(NiTCAD::device::ContactKind k) {
    switch (k) {
    case NiTCAD::device::ContactKind::ohmic: return "ohmic";
    case NiTCAD::device::ContactKind::gate: return "gate";
    case NiTCAD::device::ContactKind::electrode: return "electrode";
    }
    return "?";
}

void describe(const NiTCAD::deck::Plan& plan) {
    const auto& d = plan.device;
    std::printf("mesh: %dD, %zu nodes, %zu edges\n", d.mesh().dimension(), d.mesh().node_count(),
                d.mesh().edges().size());
    for (const auto& r : d.regions()) {
        std::printf("region %s (%s)\n", r.name.c_str(),
                    NiTCAD::device::is_insulator(r) ? "insulator" : "semiconductor");
    }
    for (const auto& c : d.contacts()) {
        std::printf("contact %s: %s, %zu nodes\n", c.name.c_str(), kind_name(c.kind), c.nodes.size());
    }
    if (plan.analysis == NiTCAD::deck::Analysis::equilibrium) {
        std::puts("analysis: equilibrium");
    } else {
        std::printf("analysis: sweep of %zu point(s)%s\n", plan.points.size(),
                    plan.restart.empty() ? "" : ", restarted from a checkpoint");
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--version") {
        print_build(NiTCAD::solve::provenance());
        return 0;
    }
    if (argc != 3) return usage();
    const std::string_view command = argv[1];
    if (command != "check" && command != "run") return usage();
    const auto deck = NiTCAD::deck::read(argv[2]);
    if (!deck) return fail(deck.error());
    const auto plan = NiTCAD::deck::build(*deck);
    if (!plan) return fail(plan.error());
    describe(*plan);
    if (command == "check") return 0;

    SetConsoleCtrlHandler(on_console_event, TRUE);
    NiTCAD::solve::RunControl control{.stop = stop.get_token(), .progress = {}};
    control.progress = [](const NiTCAD::solve::Progress& p) {
        if (!p.converged) return;
        if (p.phase == NiTCAD::solve::Phase::equilibrium) {
            std::printf("equilibrium: converged in %d iterations\n", p.iteration);
        } else {
            std::printf("point %zu of %zu: converged (%d iterations)\n", p.point + 1, p.point_count,
                        p.iteration);
        }
        std::fflush(stdout);
    };
    const auto report = NiTCAD::deck::run(*plan, control);
    if (!report) return fail(report.error());
    if (report->skipped_points > 0) {
        std::printf("restart: %zu point(s) already done\n", report->skipped_points);
    }
    for (const auto& p : report->written) std::printf("wrote %s\n", p.string().c_str());
    // A sweep's own record; an equilibrium run has none, so the provenance of now, without a
    // wall time.
    print_provenance(report->sweep ? report->sweep->run.provenance : NiTCAD::solve::provenance());
    if (report->sweep && report->sweep->stopped) {
        std::fprintf(stderr, "nitcad: the sweep stopped after %zu of %zu point(s): %s\n",
                     report->sweep->points.size(), plan->points.size() - report->skipped_points,
                     report->sweep->stopped->message.c_str());
        return 2;
    }
    return 0;
}
