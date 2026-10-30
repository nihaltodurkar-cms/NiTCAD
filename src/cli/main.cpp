// nitcad: the batch runner (ARCHITECTURE.md 11.3, Unit 31; DECISIONS.md DK10).
//
//   nitcad check <deck>   parse and build the deck, print what it describes, solve nothing
//   nitcad run <deck>     build the deck, run its analysis, write its outputs, and print the run's
//                         provenance (Unit 42)
//   nitcad plot <deck or result file> ... --out <image.png>
//                         draw a view or a curve of a deck's results (Unit 25a; plot.hpp)
//   nitcad --version      print the build: version, git revision, build type, compiler, Eigen
//
// Exit status: 0 done; 1 the deck, the run or the plot failed (the message on stderr); 2 a sweep,
// trace, transient or small-signal run stopped early (its completed part and outputs are written;
// the reason on stderr; Unit 43 for the last three); 64 usage.
// Ctrl+C or Ctrl+Break asks the run to stop: the solver stops at its next check and the completed
// points are written, as for any early stop.
#include <windows.h>

#include <algorithm>
#include <cstdio>
#include <stop_token>
#include <string>
#include <string_view>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/solve/provenance.hpp"
#include "plot.hpp"

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
    std::fputs("usage: nitcad check <deck>\n       nitcad run <deck>\n"
               "       nitcad plot <deck or result file> [options] --out <image.png>\n"
               "       nitcad --version\n",
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
    if (e.context && e.context->index) {
        std::fprintf(stderr, "  (index %zu", *e.context->index);
        if (e.context->value) std::fprintf(stderr, ", value %g", *e.context->value);
        std::fputs(")\n", stderr);
    }
    return 1;
}

// A waveform's shape, from its kind and the number of its parameters (two per corner).
const char* waveform_name(NiTCAD::solve::Waveform::Kind kind, std::size_t parameters) {
    if (kind == NiTCAD::solve::Waveform::Kind::sine) return "sine";
    return parameters == 2 ? "constant" : "piecewise linear";
}

const char* kind_name(NiTCAD::device::ContactKind k) {
    switch (k) {
    case NiTCAD::device::ContactKind::ohmic: return "ohmic";
    case NiTCAD::device::ContactKind::gate: return "gate";
    case NiTCAD::device::ContactKind::electrode: return "electrode";
    case NiTCAD::device::ContactKind::schottky: return "schottky";
    }
    return "?";
}

void describe(const NiTCAD::deck::Plan& plan) {
    const auto& d = plan.device;
    std::printf("mesh: %dD, %zu nodes, %zu edges\n", d.mesh().dimension(), d.mesh().node_count(),
                d.mesh().edges().size());
    if (plan.generated) {
        const auto& r = plan.generated->report;
        if (d.mesh().dimension() == 3) {
            std::printf("generated: %zu tetrahedra, %zu slivers, %zu exempt, smallest dihedral %.2f deg, "
                        "largest radius-edge %.3f\n",
                        r.tetrahedra, r.slivers, r.exempt_tetrahedra, r.min_dihedral_deg, r.max_radius_edge);
        } else {
            std::printf("generated: %zu triangles, %zu exempt, smallest angle %.2f deg\n", r.triangles,
                        r.exempt_triangles, r.min_angle_deg);
        }
        if (plan.generated->generation_digest != 0) {
            std::printf("generation digest: %016llx (%s)\n",
                        static_cast<unsigned long long>(plan.generated->generation_digest), plan.mesh_type.c_str());
        }
        if (const auto* c = d.interface_certificate()) {
            std::printf("certificate: P4 %zu cross faces (%zu in corner zones)", c->report.cross_faces,
                        c->report.corner_faces);
            if (c->patches) {
                std::printf("; P6 %zu patch-boundary pieces (%zu in zones)", c->patches->cross_pieces,
                            c->patches->corner_pieces);
            }
            std::putchar('\n');
        }
    }
    for (const auto& r : d.regions()) {
        std::printf("region %s (%s)\n", r.name.c_str(),
                    NiTCAD::device::is_insulator(r) ? "insulator" : "semiconductor");
    }
    for (const auto& c : d.contacts()) {
        std::printf("contact %s: %s, %zu nodes\n", c.name.c_str(), kind_name(c.kind), c.nodes.size());
    }
    using NiTCAD::deck::Analysis;
    switch (plan.analysis) {
    case Analysis::equilibrium: std::puts("analysis: equilibrium"); break;
    case Analysis::sweep:
        std::printf("analysis: sweep of %zu point(s)%s\n", plan.points.size(),
                    plan.restart.empty() ? "" : ", restarted from a checkpoint");
        break;
    case Analysis::transient: {
        std::printf("analysis: transient to %g s\n", plan.transient.t_end_s);
        const auto contacts = d.contacts();
        for (std::size_t k = 0; k < contacts.size(); ++k) {
            std::printf("  contact %s: %s\n", contacts[k].name.c_str(),
                        waveform_name(plan.waveforms[k].kind(), plan.waveforms[k].parameters().size()));
        }
        break;
    }
    case Analysis::small_signal: {
        const auto& f = plan.small_signal.frequencies_Hz;
        std::printf("analysis: small signal at %zu point(s), %zu frequenc%s from %g to %g Hz\n",
                    plan.points.size(), f.size(), f.size() == 1 ? "y" : "ies",
                    *std::min_element(f.begin(), f.end()), *std::max_element(f.begin(), f.end()));
        break;
    }
    case Analysis::trace:
        std::printf("analysis: trace of contact %s from %g to %g V\n",
                    d.contacts()[plan.trace.contact].name.c_str(), plan.trace.start_V[plan.trace.contact],
                    plan.trace.end_V);
        break;
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--version") {
        print_build(NiTCAD::solve::provenance());
        return 0;
    }
    if (argc >= 2 && std::string_view(argv[1]) == "plot") return NiTCAD::cli::plot(argc, argv);
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
    const bool tracing = plan->analysis == NiTCAD::deck::Analysis::trace;
    control.progress = [tracing](const NiTCAD::solve::Progress& p) {
        if (!p.converged) return;
        using NiTCAD::solve::Phase;
        if (p.phase == Phase::equilibrium) {
            std::printf("equilibrium: converged in %d iterations\n", p.iteration);
        } else if (p.phase == Phase::transient) {
            std::printf("t = %.6g s: converged (%d iterations)\n", p.time_s, p.iteration);
        } else if (p.phase == Phase::small_signal) {
            std::printf("point %zu of %zu, %.6g Hz: solved\n", p.point + 1, p.point_count, p.frequency_Hz);
        } else if (tracing) {
            std::printf("traced point %zu: converged (%d iterations)\n", p.point + 1, p.iteration);
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
        if (plan->analysis == NiTCAD::deck::Analysis::trace) {
            std::fprintf(stderr, "nitcad: the trace stopped after %zu point(s): %s\n",
                         report->sweep->points.size(), report->sweep->stopped->message.c_str());
        } else {
            std::fprintf(stderr, "nitcad: the sweep stopped after %zu of %zu point(s): %s\n",
                         report->sweep->points.size(), plan->points.size() - report->skipped_points,
                         report->sweep->stopped->message.c_str());
        }
        return 2;
    }
    if (report->transient && report->transient->stopped) {
        const auto& t = *report->transient;
        std::fprintf(stderr, "nitcad: the transient stopped at t = %g s after %zu step(s): %s\n",
                     t.points.empty() ? 0.0 : t.points.back().time_s,
                     t.points.empty() ? std::size_t{0} : t.points.size() - 1, t.stopped->message.c_str());
        return 2;
    }
    if (report->small_signal && report->small_signal->stopped) {
        std::fprintf(stderr, "nitcad: the small-signal run stopped after %zu of %zu point(s): %s\n",
                     report->small_signal->points.size(), plan->points.size(),
                     report->small_signal->stopped->message.c_str());
        return 2;
    }
    return 0;
}
