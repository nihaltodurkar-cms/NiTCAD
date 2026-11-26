// nitcad: the batch runner (ARCHITECTURE.md 11.3, Unit 31; DECISIONS.md DK10).
//
//   nitcad check <deck>   parse and build the deck, print what it describes, solve nothing
//   nitcad run <deck>     build the deck, run its analysis, write its outputs, and print the run's
//                         provenance (Unit 42)
//   check and run take [--set <name>=<value>]...: a [parameters] value overridden, as if written in the deck
//                         (Unit 26n, DT6; deck/dimension_sweep.hpp); the deck file is not edited
//   nitcad sweep <deck> <name>=<values> [--set <name>=<value>]...
//                         the deck run once per value (`20nm, 30nm` or `20nm : 45nm : 5nm`), each run's
//                         outputs in <deck stem>_sweep_<name>/<name>=<value>/ (Unit 26n, DT6)
//   nitcad plot <deck or result file> ... --out <image.png>
//                         draw a view or a curve of a deck's results (Unit 25a; plot.hpp)
//   nitcad verify [--tier fast|full|nightly] [--id <id>] [--report <file>]
//                 [--mkl <mkl_rt.dll>] [--threads <n>]   (PARDISO for TVB-04's 3D levels)
//                         run the verification benchmarks of TVR-1 (Unit 55; verify.hpp, its own exit status)
//   nitcad learn manifest <manifest> <history>...
//                         the manifest of adaptation histories for cross-run learning (Unit 17b; learn.hpp)
//   nitcad compare <a.nitsweep> <b.nitsweep> [--rtol <x>] [--potential-tol <V>]
//                         compare two sweeps by tolerances (Unit 60; compare.hpp, its own exit
//                         status: 0 match, 4 differ, 1 unreadable or not comparable, 64 usage)
//   nitcad materials [<name> [<composition>]...]
//                         list the material registry, or one material's parameters and sources (Unit 68)
//   nitcad --version      print the build: version, git revision, build type, compiler, Eigen
//
// Exit status: 0 done; 1 the deck, the run or the plot failed (the message on stderr); 2 a sweep,
// trace, transient or small-signal run stopped early (its completed part and outputs are written;
// the reason on stderr; Unit 43 for the last three); 3 an [adapt] run with certify = required
// whose result is not certified (Unit 17, AD13: its outputs are written, the certificates printed);
// 64 usage.
// Ctrl+C or Ctrl+Break asks the run to stop: the solver stops at its next check and the completed
// points are written, as for any early stop.
#include <windows.h>

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <expected>
#include <filesystem>
#include <format>
#include <stop_token>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/dimension_sweep.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/physics/materials.hpp"
#include "NiTCAD/solve/provenance.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "plot.hpp"
#include "compare.hpp"
#include "learn.hpp"
#include "verify.hpp"

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
    std::fputs("usage: nitcad check <deck> [--set <name>=<value>]...\n"
               "       nitcad run <deck> [--set <name>=<value>]...\n"
               "       nitcad sweep <deck> <name>=<values> [--set <name>=<value>]...\n"
               "       nitcad plot <deck or result file> [options] --out <image.png>\n"
               "       nitcad verify [--tier fast|full|nightly] [--id <id>] [--report <file>]\n"
               "                     [--mkl <mkl_rt.dll>] [--threads <n>]\n"
               "       nitcad learn manifest <manifest> <history.nithist>...\n"
               "       nitcad learn train <manifest> <model.nitlearn> [--holdout <history.nithist>]...\n"
               "       nitcad compare <a.nitsweep> <b.nitsweep> [--rtol <x>] [--potential-tol <V>]\n"
               "       nitcad materials [<name> [<composition>]...]\n"
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
        if (r.prisms != 0) {
            // A swept mesh (Unit 54, UD10): the cross-section's quality, the product's counts.
            std::printf("swept: %zu axis lines, %zu prisms, %zu tetrahedra; %zu layers, %zu corner blocks; "
                        "cross-section smallest angle %.2f deg (%zu exempt)\n",
                        r.axis_lines, r.prisms, r.tetrahedra, r.layers, r.corner_blocks, r.min_angle_deg,
                        r.exempt_triangles);
        } else if (d.mesh().dimension() == 3) {
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
    if (plan.adapt) {
        const auto& o = plan.adapt->options;
        std::string goals;
        bool derived = false, vt = false;
        for (const auto& g : o.goals) {
            derived = derived || g.kind != NiTCAD::solve::GoalKind::current;
            vt = vt || g.kind == NiTCAD::solve::GoalKind::vt;
        }
        // Unit 17c (DG9): with a derived goal each goal is "<kind> of <contact>"; else as before.
        for (const auto& g : o.goals) {
            goals += (goals.empty() ? "" : ", ") +
                     (derived ? std::string(NiTCAD::solve::to_string(g.kind)) + " of " : std::string()) +
                     d.contacts()[g.contact].name;
        }
        if (vt) goals += std::format(" (vt to {:g} V)", o.vt_tolerance);
        std::printf("adapt: %s%s to %g (relative), policy %s, weighting %s, verify %s, certify %s, "
                    "at most %d passes\n",
                    derived ? "" : "current of ", goals.c_str(), o.tolerance, NiTCAD::solve::to_string(o.policy),
                    o.weighting == NiTCAD::solve::Weighting::adjoint ? "adjoint" : "none", o.verify ? "richardson" : "off",
                    plan.adapt->certify_required ? "required" : "report", o.max_passes);
        if (const auto& m = plan.adapt->predictor) {
            std::printf("adapt: predictor %s (%016llx): format %llu, schema %llu, %zu trees, trained on %zu devices; "
                        "held-out rank correlation %.3f (ridge %.3f)\n",
                        plan.adapt->predictor_path.generic_string().c_str(),
                        static_cast<unsigned long long>(plan.adapt->predictor_digest),
                        static_cast<unsigned long long>(m->format), static_cast<unsigned long long>(m->schema),
                        m->trees.size(), m->training.size(), m->metrics.rank_correlation, m->metrics.ridge_rank_correlation);
        }
    }
}

// An adaptive run's passes and certificates (Unit 17, AD8, AD13).
void print_adapt(const NiTCAD::results::AdaptReport& a) {
    // Unit 17b (LB8): the learned predictor's use, or why it was not used.
    if (a.predictor_status == "used") {
        std::printf("adapt: predictor %s (%016llx) predicted the start mesh\n", a.predictor.c_str(),
                    static_cast<unsigned long long>(a.predictor_digest));
    } else if (a.predictor_status == "out_of_range") {
        std::printf("adapt: out of the predictor's training range (%.1f%% of the start mesh's edges): "
                    "the a priori prediction used\n",
                    100.0 * a.predictor_out_of_range);
    }
    for (std::size_t k = 0; k < a.passes.size(); ++k) {
        const auto& p = a.passes[k];
        std::printf("pass %zu (%s): %llu nodes, %llu DOFs, relative estimate %.3g, %llu marked, %.3f s\n", k,
                    p.source.c_str(), static_cast<unsigned long long>(p.nodes), static_cast<unsigned long long>(p.dofs),
                    p.relative_estimate, static_cast<unsigned long long>(p.marked),
                    p.seconds_mesh + p.seconds_transfer + p.seconds_solve + p.seconds_estimate);
    }
    std::printf("adapt: stopped (%s), policy %s, fallback %s; verification on %llu nodes, %.3f s\n", a.stop.c_str(),
                a.policy.c_str(), a.fallback.c_str(), static_cast<unsigned long long>(a.verification_nodes),
                a.verification_seconds);
    for (const auto& c : a.certificates) {
        const std::string kind = NiTCAD::results::goal_kind(a, c.goal);
        if (kind == "current") {
            std::printf("certificate: current of %s at point %llu: %.10g A, estimate %.3g, verified error %.3g "
                        "(order %.2f), allowed %.3g: %s\n",
                        a.goal_contacts[c.goal].c_str(), static_cast<unsigned long long>(c.point), c.value, c.estimate,
                        c.verified_error, c.order, c.allowed, c.status.c_str());
            continue;
        }
        // Unit 17c (DG6, DG8): a derived goal's certificate covers the mesh's error at this gate step,
        // not the step's own; a moved window is named.
        const std::string unit = NiTCAD::results::goal_unit(a, c.goal);
        std::printf("certificate: %s of %s (window at point %llu): %.10g %s, estimate %.3g, verified error %.3g "
                    "(order %.2f), allowed %.3g: %s%s (the mesh's error at this sweep's gate step, not the step's)\n",
                    kind.c_str(), a.goal_contacts[c.goal].c_str(), static_cast<unsigned long long>(c.point), c.value,
                    unit.c_str(), c.estimate, c.verified_error, c.order, c.allowed, c.status.c_str(),
                    c.window_moved ? ", window moved on M'" : "");
    }
    std::printf("adapt: %s (an estimate and an independent Richardson check, not a proof)\n", a.status.c_str());
}

}  // namespace

namespace {

// Build, run and report one plan: the exit status of `nitcad run`.
int run_plan(const NiTCAD::deck::Plan& plan) {
    NiTCAD::solve::RunControl control{.stop = stop.get_token(), .progress = {}};
    const bool tracing = plan.analysis == NiTCAD::deck::Analysis::trace;
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
    const auto report = NiTCAD::deck::run(plan, control);
    if (!report) return fail(report.error());
    if (report->skipped_points > 0) {
        std::printf("restart: %zu point(s) already done\n", report->skipped_points);
    }
    if (report->adapt) print_adapt(*report->adapt);
    for (const auto& p : report->written) std::printf("wrote %s\n", p.string().c_str());
    // The run record's warnings (Unit 62).
    const NiTCAD::results::RunRecord* record = report->sweep          ? &report->sweep->run
                                               : report->transient    ? &report->transient->run
                                               : report->small_signal ? &report->small_signal->run
                                                                      : nullptr;
    if (record != nullptr) {
        for (const std::string& w : record->warnings) std::printf("warning: %s\n", w.c_str());
    }
    // A sweep's own record; an equilibrium run has none, so the provenance of now, without a
    // wall time.
    print_provenance(report->sweep ? report->sweep->run.provenance : NiTCAD::solve::provenance());
    if (report->sweep && report->sweep->stopped) {
        if (plan.analysis == NiTCAD::deck::Analysis::trace) {
            std::fprintf(stderr, "nitcad: the trace stopped after %zu point(s): %s\n",
                         report->sweep->points.size(), report->sweep->stopped->message.c_str());
        } else {
            std::fprintf(stderr, "nitcad: the sweep stopped after %zu of %zu point(s): %s\n",
                         report->sweep->points.size(), plan.points.size() - report->skipped_points,
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
                     report->small_signal->points.size(), plan.points.size() - report->skipped_points,
                     report->small_signal->stopped->message.c_str());
        return 2;
    }
    if (report->adapt && report->adapt->certify_required && report->adapt->status != "certified") {
        std::fprintf(stderr, "nitcad: the adaptive result is %s; the deck requires it certified (certify = required)\n",
                     report->adapt->status.c_str());
        return 3;
    }
    return 0;
}

// `--set name=value` pairs from argv[from..]. Errors: anything else, a `--set` without its pair.
std::expected<std::vector<NiTCAD::deck::Override>, std::string> overrides_of(int argc, char** argv, int from) {
    std::vector<NiTCAD::deck::Override> out;
    for (int i = from; i < argc; ++i) {
        if (std::string_view(argv[i]) != "--set" || i + 1 == argc) return std::unexpected("expected --set <name>=<value>");
        auto o = NiTCAD::deck::parse_override(argv[++i]);
        if (!o) return std::unexpected(o.error());
        out.push_back(std::move(*o));
    }
    return out;
}

int usage_error(const std::string& what) {
    std::fprintf(stderr, "nitcad: %s\n", what.c_str());
    return usage();
}

// nitcad sweep <deck> <name>=<values> [--set <name>=<value>]... (Unit 26n, DT6): the deck run once per value of one
// parameter, each run's outputs in its own folder beside the deck. Every value is checked before the first run.
// Exit status: 0 every run done; 1 a value's build or run failed (the runs before it kept); 2 stopped (Ctrl+C after
// the current run, or a run stopped early); 3 a run's [adapt] result not certified where required (the rest still
// run); 64 usage, an unknown parameter or a value of the wrong kind.
int sweep(int argc, char** argv) {
    if (argc < 4) return usage();
    const std::filesystem::path deck_path = argv[2];
    const auto deck = NiTCAD::deck::read(deck_path);
    if (!deck) return fail(deck.error());
    const auto params = NiTCAD::deck::parameters(*deck);
    if (!params) return fail(params.error());
    const auto swept = NiTCAD::deck::parse_override(argv[3]);
    if (!swept) return usage_error(swept.error());
    const auto held = overrides_of(argc, argv, 4);
    if (!held) return usage_error(held.error());
    const auto p = std::ranges::find(*params, swept->name, &NiTCAD::deck::Parameter::name);
    if (p == params->end()) return usage_error(std::format("{}: no parameter '{}' in [parameters]", deck->source, swept->name));
    const auto values = NiTCAD::deck::sweep_values(swept->value, p->length);
    if (!values) return usage_error(std::format("{}: {}", swept->name, values.error()));
    std::vector<NiTCAD::deck::Deck> decks;
    for (const std::string& v : *values) {
        std::vector<NiTCAD::deck::Override> o = *held;
        o.push_back({.name = swept->name, .value = v});
        auto d = NiTCAD::deck::with_overrides(*deck, o);
        if (!d) return usage_error(d.error().message);
        decks.push_back(std::move(*d));
    }
    const std::filesystem::path folder = NiTCAD::deck::sweep_directory(deck_path, swept->name);
    int status = 0;
    for (std::size_t k = 0; k < decks.size(); ++k) {
        if (stop.stop_requested()) {
            std::fprintf(stderr, "nitcad: the sweep of %s stopped after %zu of %zu value(s)\n", swept->name.c_str(), k,
                         decks.size());
            return 2;
        }
        const std::string& v = (*values)[k];
        const std::filesystem::path point = NiTCAD::deck::sweep_point_directory(folder, swept->name, v);
        std::printf("sweep %s = %s (%zu of %zu): %s\n", swept->name.c_str(), v.c_str(), k + 1, decks.size(),
                    point.string().c_str());
        std::fflush(stdout);
        auto plan = NiTCAD::deck::build(decks[k]);
        if (!plan) return fail(plan.error());
        plan->outputs = NiTCAD::deck::with_outputs_in(std::move(plan->outputs), decks[k].directory, point);
        std::error_code ec;
        std::filesystem::create_directories(point, ec);
        if (ec) {
            std::fprintf(stderr, "nitcad: cannot create %s: %s\n", point.string().c_str(), ec.message().c_str());
            return 1;
        }
        describe(*plan);
        const int s = run_plan(*plan);
        if (s == 1 || s == 2) return s;
        if (s == 3) status = 3;
    }
    return status;
}

// `nitcad materials [<name> [<composition>]...]` (Unit 68, MAD12): the registry, one line per
// material, or one material's set, every parameter a deck can set and the set's sources.
int materials(int argc, char** argv) {
    using namespace NiTCAD::physics;
    const auto kind = [](MaterialKind k) {
        return k == MaterialKind::semiconductor ? "semiconductor" : k == MaterialKind::insulator ? "insulator" : "metal";
    };
    if (argc == 2) {
        for (const MaterialEntry& e : material_registry()) {
            std::printf("%-16s %-14s %s\n", std::string(e.name).c_str(), kind(e.kind), std::string(e.summary).c_str());
        }
        std::puts("\nmetal presets (a contact's metal = <name>):");
        for (const MetalPreset& m : metal_presets()) {
            const std::string value = m.kind == MetalPresetKind::metal      ? std::format("{:.12g} eV", m.work_function_eV)
                                      : m.kind == MetalPresetKind::silicide ? std::format("barrier {:.12g} eV on n-Si (Schottky, silicon)", m.barrier_eV)
                                                                            : std::string(m.source.where);
            std::printf("%-16s %s\n", std::string(m.name).c_str(), value.c_str());
        }
        return 0;
    }
    if (const MetalPreset* m = find_metal_preset(argv[2]); m && argc == 3) {
        const std::string value = m->kind == MetalPresetKind::metal ? std::format(", work function {:.12g} eV", m->work_function_eV)
                                  : m->kind == MetalPresetKind::silicide
                                      ? std::format(", silicide barrier {:.12g} eV on n-Si", m->barrier_eV)
                                      : std::string{};
        std::printf("%s: metal preset%s\n  source: %s (%s)\n", argv[2], value.c_str(), std::string(m->source.source).c_str(),
                    std::string(m->source.where).c_str());
        return 0;
    }
    const MaterialEntry* e = find_material(argv[2]);
    if (!e) {
        std::vector<std::string_view> names;
        for (const MaterialEntry& m : material_registry()) names.push_back(m.name);
        for (const MetalPreset& m : metal_presets()) names.push_back(m.name);
        const std::string_view nearest = nearest_name(argv[2], names);
        std::fprintf(stderr, "nitcad materials: unknown material '%s'%s (nitcad materials lists them)\n", argv[2],
                     nearest.empty() ? "" : std::format(" (did you mean '{}'?)", nearest).c_str());
        return 1;
    }
    std::vector<double> arguments;
    for (int i = 3; i < argc; ++i) {
        char* end = nullptr;
        const double v = std::strtod(argv[i], &end);
        if (end == argv[i] || *end != '\0') return usage_error(std::format("'{}' is not a number", argv[i]));
        arguments.push_back(v);
    }
    if (arguments.size() < e->arguments_min || arguments.size() > e->arguments_max) {
        return usage_error(e->arguments_max == 0 ? std::format("'{}' takes no arguments", e->name) : std::string(e->usage));
    }
    const auto set = e->make(arguments);
    if (!set) return fail(set.error());
    std::printf("%s: %s, %s\n", std::string(e->name).c_str(), kind(e->kind), std::string(e->summary).c_str());
    for (std::size_t i = 0; i < e->sources.size(); ++i) {
        const Citation& c = e->sources[i];
        std::printf("  source [%zu]: %s (%s)\n", i + 1, std::string(c.source).c_str(), std::string(c.where).c_str());
    }
    if (e->parameter_sources.empty()) std::puts("  (the sources above cover the whole set)");
    for (const auto& [path, value] : parameter_values(*set)) {
        std::string cited;
        for (const ParameterCitation& c : e->parameter_sources) {
            if (c.path == path) cited = std::format("   [{}] {}", c.source + 1, c.basis);
        }
        std::printf("  %s = %s%s\n", std::string(path).c_str(), std::format("{:.12g}", value).c_str(), cited.c_str());
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--version") {
        print_build(NiTCAD::solve::provenance());
        return 0;
    }
    if (argc >= 2 && std::string_view(argv[1]) == "plot") return NiTCAD::cli::plot(argc, argv);
    if (argc >= 2 && std::string_view(argv[1]) == "learn") return NiTCAD::cli::learn(argc, argv);
    if (argc >= 2 && std::string_view(argv[1]) == "compare") return NiTCAD::cli::compare(argc, argv);
    if (argc >= 2 && std::string_view(argv[1]) == "materials") return materials(argc, argv);
    if (argc >= 2 && std::string_view(argv[1]) == "verify") {
        return NiTCAD::cli::verify(argc, argv, NiTCAD::verify::benchmark::benchmarks());
    }
    SetConsoleCtrlHandler(on_console_event, TRUE);
    if (argc >= 2 && std::string_view(argv[1]) == "sweep") return sweep(argc, argv);
    if (argc < 3) return usage();
    const std::string_view command = argv[1];
    if (command != "check" && command != "run") return usage();
    const auto overrides = overrides_of(argc, argv, 3);
    if (!overrides) return usage_error(overrides.error());
    const auto read = NiTCAD::deck::read(argv[2]);
    if (!read) return fail(read.error());
    if (const auto params = NiTCAD::deck::parameters(*read); !params) return fail(params.error());
    const auto deck = NiTCAD::deck::with_overrides(*read, *overrides);
    if (!deck) return usage_error(deck.error().message);
    const auto plan = NiTCAD::deck::build(*deck);
    if (!plan) return fail(plan.error());
    describe(*plan);
    if (command == "check") return 0;
    return run_plan(*plan);
}
