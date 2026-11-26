#include "NiTCAD/deck/run.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>
#include <system_error>
#include <utility>

#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/io/adapt_report.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/io/learn_files.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/solve/adapt_learn.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/solve/transient.hpp"

namespace NiTCAD::deck {

namespace {

base::Error invalid(std::string message) {
    return {base::ErrorCode::invalid_input, std::move(message), std::nullopt};
}

// The same bias to 1e-12 relative (1e-12 V near zero): a series written as 0 : 0.3 : 0.1 ends at
// 0.3 exactly, the fourth value of 0 : 0.6 : 0.1 is 3 x 0.1 = 0.30000000000000004 (DK9).
bool same_bias(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > 1e-12 * std::max({1.0, std::abs(a[i]), std::abs(b[i])})) return false;
    }
    return true;
}

// A restart's checkpoint and the run agree on density gradient (Unit 62, QV12): the checkpoint
// carries the quantum potentials exactly when the run solves the model.
std::optional<base::Error> restart_confinement(const std::filesystem::path& restart, const io::Checkpoint& c,
                                               const assemble::PhysicsModels& models) {
    const bool run = models.quantum_confinement == assemble::QuantumConfinement::density_gradient;
    const bool saved = !c.fields.quantum_potential_n_eV.empty();
    if (run == saved) return std::nullopt;
    return invalid("the restart checkpoint '" + restart.string() + "' was solved " +
                   (saved ? "with" : "without") + " quantum_confinement density_gradient and this run " +
                   (run ? "solves" : "does not solve") + " it: start from it with initial instead");
}

std::string format(double v) {
    char buffer[32];
    const int n = std::snprintf(buffer, sizeof buffer, "%.17g", v);
    return std::string(buffer, static_cast<std::size_t>(n));
}

std::expected<void, base::Error> write_text(const std::filesystem::path& path, const std::string& text);

}  // namespace

std::expected<void, base::Error> write_table(const std::filesystem::path& path,
                                             const device::Device& device,
                                             const results::Sweep& sweep,
                                             const results::AdaptReport* adapt) {
    // With the weighted current (Unit 32) every point has it: two more columns per contact.
    const bool weighted = !sweep.points.empty() && !sweep.points.front().weighted_current.empty();
    std::string text = "point";
    for (const device::Contact& c : device.contacts()) text += "," + c.name + " V";
    for (const device::Contact& c : device.contacts()) text += "," + c.name + " A";
    if (weighted) {
        for (const device::Contact& c : device.contacts()) text += "," + c.name + " weighted A";
        for (const device::Contact& c : device.contacts()) {
            text += "," + c.name + " weighted resolution A";
        }
    }
    // Gate tunnelling (Unit 51, GD10): each tunnelling contact's current by process.
    const bool tunnelling = !sweep.points.empty() && !sweep.points.front().tunnelling_ecb.empty();
    std::vector<std::size_t> tunnels;
    for (std::size_t c = 0; c < device.contacts().size(); ++c) {
        if (tunnelling && device.contacts()[c].tunnelling != 0) tunnels.push_back(c);
    }
    for (const std::size_t c : tunnels) {
        for (const char* what : {" ECB A", " HVB A", " EVB A"}) {
            text += "," + device.contacts()[c].name + what;
        }
    }
    text += ",substeps";
    // The certificates (Unit 17, AD13): per goal, at the adaptation points; empty at the others.
    // Unit 17c (DG9): a derived goal's columns are "<kind> <contact> ...", in its unit, filled at its
    // certificate's point.
    if (adapt != nullptr) {
        for (std::size_t g = 0; g < adapt->goal_contacts.size(); ++g) {
            const std::string name = results::goal_name(*adapt, g);
            const std::string unit = results::goal_unit(*adapt, g);
            text += "," + name + " certificate," + name + " estimate " + unit + "," + name + " verified error " + unit +
                    "," + name + " allowed " + unit + "," + name + " order";
        }
    }
    text += "\n";
    for (std::size_t k = 0; k < sweep.points.size(); ++k) {
        const results::BiasPoint& p = sweep.points[k];
        text += std::to_string(k);
        for (const double v : p.bias_V) text += "," + format(v);
        for (const double i : p.terminal_current) text += "," + format(i);
        if (weighted) {
            for (const double i : p.weighted_current) text += "," + format(i);
            for (const double i : p.weighted_current_resolution) text += "," + format(i);
        }
        for (const std::size_t c : tunnels) {
            text += "," + format(p.tunnelling_ecb[c]) + "," + format(p.tunnelling_hvb[c]) + "," +
                    format(p.tunnelling_evb[c]);
        }
        text += "," + std::to_string(p.bias_substeps);
        if (adapt != nullptr) {
            for (std::size_t g = 0; g < adapt->goal_contacts.size(); ++g) {
                const auto c = std::ranges::find_if(adapt->certificates, [&](const results::AdaptCertificateRow& r) {
                    return r.goal == g && r.point == k;
                });
                if (c == adapt->certificates.end()) {
                    text += ",,,,,";
                } else {
                    text += "," + c->status + "," + format(c->estimate) + "," + format(c->verified_error) + "," +
                            format(c->allowed) + "," + format(c->order);
                }
            }
        }
        text += "\n";
    }
    return write_text(path, text);
}

namespace {

std::string pass_name(std::size_t pass) {
    char name[32];
    std::snprintf(name, sizeof name, "pass_%02zu", pass);
    return name;
}

}  // namespace

std::filesystem::path pass_device_file(const std::filesystem::path& directory, std::size_t pass) {
    return directory / (pass_name(pass) + ".nitdev");
}

std::filesystem::path pass_mesh_file(const std::filesystem::path& directory, std::size_t pass) {
    return directory / (pass_name(pass) + ".nitmesh");
}

std::expected<void, base::Error> write_transient_table(const std::filesystem::path& path,
                                                       const device::Device& device,
                                                       const results::Transient& run) {
    std::string text = "time s,step s,order";
    for (const char* what : {" V", " A", " conduction A", " displacement A"}) {
        for (const device::Contact& c : device.contacts()) text += "," + c.name + what;
    }
    // Heavy-ion tracks (Unit 65, HD9): each contact's collected charge and the tracks' deposit so far.
    const bool ions = !run.points.empty() && run.points.front().heavy_ion_charge.has_value();
    std::vector<std::vector<double>> collected;
    if (ions) {
        for (std::size_t c = 0; c < device.contacts().size(); ++c) {
            auto q = analysis::collected_charge(run, c);
            if (!q) return std::unexpected(std::move(q.error()));
            collected.push_back(std::move(*q));
            text += "," + device.contacts()[c].name + " collected C";
        }
        text += ",deposited C";
    }
    text += "\n";
    for (std::size_t k = 0; k < run.points.size(); ++k) {
        const results::TimePoint& p = run.points[k];
        text += format(p.time_s) + "," + format(p.step_s) + "," + std::to_string(p.order);
        for (const auto* v : {&p.bias_V, &p.terminal_current, &p.conduction_current, &p.displacement_current}) {
            for (const double x : *v) text += "," + format(x);
        }
        if (ions) {
            for (const auto& q : collected) text += "," + format(q[k]);
            text += "," + format(p.heavy_ion_charge.value_or(0.0));
        }
        text += "\n";
    }
    return write_text(path, text);
}

std::expected<void, base::Error> write_small_signal_table(const std::filesystem::path& path,
                                                          const device::Device& device,
                                                          const results::SmallSignal& run) {
    const auto contacts = device.contacts();
    std::string text = "point,frequency Hz";
    for (const device::Contact& c : contacts) text += "," + c.name + " V";
    for (const char* what : {"G ", "C "}) {
        for (const device::Contact& i : contacts) {
            for (const device::Contact& j : contacts) {
                text += std::string(",") + what + i.name + "," + j.name + (what[0] == 'G' ? " S" : " F");
            }
        }
    }
    text += "\n";
    for (std::size_t k = 0; k < run.points.size(); ++k) {
        for (std::size_t f = 0; f < run.frequency_Hz.size(); ++f) {
            text += std::to_string(k) + "," + format(run.frequency_Hz[f]);
            for (const double v : run.points[k].dc.bias_V) text += "," + format(v);
            for (std::size_t i = 0; i < run.contacts; ++i) {
                for (std::size_t j = 0; j < run.contacts; ++j) text += "," + format(run.conductance(k, f, i, j));
            }
            for (std::size_t i = 0; i < run.contacts; ++i) {
                for (std::size_t j = 0; j < run.contacts; ++j) text += "," + format(run.capacitance(k, f, i, j));
            }
            text += "\n";
        }
    }
    return write_text(path, text);
}

namespace {

// Atomic like the binary files: a temporary beside the target, renamed over it.
std::expected<void, base::Error> write_text(const std::filesystem::path& path, const std::string& text) {
    std::filesystem::path temporary = path;
    temporary += ".tmp";
    {
        std::ofstream f(temporary, std::ios::binary | std::ios::trunc);
        f.write(text.data(), static_cast<std::streamsize>(text.size()));
        f.flush();
        if (!f) {
            f.close();
            std::error_code ignored;
            std::filesystem::remove(temporary, ignored);
            return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                               "cannot write '" + path.string() + "'", std::nullopt});
        }
    }
    std::error_code ec;
    std::filesystem::rename(temporary, path, ec);
    if (ec) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                           "cannot write '" + path.string() + "': " + ec.message(),
                                           std::nullopt});
    }
    return {};
}

}  // namespace

std::expected<RunReport, base::Error> run(const Plan& plan, const solve::RunControl& control) {
    RunReport report;
    if (plan.options.linear.backend == linalg::SolverBackend::mkl_pardiso) {
        auto mkl = linalg::load_mkl_runtime(plan.mkl_runtime);
        if (!mkl) return std::unexpected(std::move(mkl.error()));
    }
    const auto written = [&report](const std::filesystem::path& p,
                                   std::expected<void, base::Error> r) -> std::expected<void, base::Error> {
        if (!r) return r;
        report.written.push_back(p);
        return {};
    };
    // FD14: a run on a mesh generated from a complex names the mesh type and its generation digest.
    const auto stamp = [&plan](results::RunRecord& run) {
        if (plan.generated && plan.generated->generation_digest != 0) {
            run.mesh_type = plan.mesh_type;
            run.mesh_generation_digest = plan.generated->generation_digest;
        }
    };
    const Outputs& out = plan.outputs;
    if (plan.adapt) {
        // Unit 17 (AD8): the loop rebuilds the deck per pass; every output is the final mesh's, the
        // passes' devices and meshes go to adapt_meshes as they are solved.
        const AdaptPlan& a = *plan.adapt;
        solve::AdaptOptions o = a.options;
        if (!out.adapt_meshes.empty()) {
            std::error_code ec;
            std::filesystem::create_directories(out.adapt_meshes, ec);
            if (ec) {
                return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                                   "cannot create '" + out.adapt_meshes.string() + "': " + ec.message(),
                                                   std::nullopt});
            }
            o.on_mesh = [&](std::size_t pass, const solve::AdaptableMesh& m) -> std::expected<void, base::Error> {
                const auto dev = pass_device_file(out.adapt_meshes, pass);
                if (auto w = written(dev, io::save_device(dev, m.device)); !w) return w;
                const auto msh = pass_mesh_file(out.adapt_meshes, pass);
                return written(msh, m.generated ? io::save_generated_mesh(msh, *m.generated) : io::save_mesh(msh, m.device.mesh()));
            };
        }
        // Unit 17c (DG10): feature 12, the run's goal kind: 0 when every goal is a current, else its
        // first derived goal's (1 vt, 2 gm, 3 swing).
        double goal_kind = 0.0;
        for (const solve::AdaptGoal& g : o.goals) {
            if (g.kind != solve::GoalKind::current) {
                goal_kind = static_cast<double>(static_cast<unsigned>(g.kind));
                break;
            }
        }
        // Unit 17b (LB11): the learning history, each pass's features as it is solved.
        std::optional<solve::HistoryRecorder> recorder;
        if (!out.adapt_history.empty()) {
            recorder.emplace(plan.points, plan.options, goal_kind);
            o.on_mesh = [&, saving = std::move(o.on_mesh)](std::size_t pass,
                                                           const solve::AdaptableMesh& m) -> std::expected<void, base::Error> {
                if (saving) {
                    if (auto w = saving(pass, m); !w) return w;
                }
                return recorder->record(pass, m);
            };
        }
        // Unit 17b (LB4, LB8): the learned predictor replaces ID2's prediction of the start mesh (the
        // deck's own mesh, AD10) unless too many of its edges are outside the model's training range.
        std::string predictor_status;
        double predictor_out = 0.0;
        if (a.predictor) {
            auto f = solve::edge_features(plan.device, plan.points, plan.options, a.predictor->materials, goal_kind);
            if (!f) return std::unexpected(std::move(f.error()));
            predictor_out = solve::out_of_range_fraction(*a.predictor, *f);
            if (predictor_out > solve::out_of_range_limit) {
                predictor_status = "out_of_range";
            } else {
                predictor_status = "used";
                o.intelligent.first_prediction = [model = a.predictor, &plan, goal_kind](const device::Device& d) {
                    auto t = solve::learned_targets(*model, d, plan.points, plan.options, goal_kind);
                    return t ? std::move(*t) : solve::SizeTargets{};
                };
            }
        }
        auto r = solve::adapt(device_builder(plan), plan.points, plan.options, o, control);
        if (!r) return std::unexpected(std::move(r.error()));
        solve::stamp_adapt(r->sweep.run, o, a.certify_required);
        if (a.predictor) solve::stamp_predictor(r->sweep.run, a.predictor_digest, *a.predictor, predictor_status == "used");
        const solve::AdaptableMesh& fin = *r->mesh;
        if (fin.generated && fin.generated->generation_digest != 0) {
            r->sweep.run.mesh_type = plan.mesh_type;
            r->sweep.run.mesh_generation_digest = fin.generated->generation_digest;
        }
        if (recorder) {
            auto h = recorder->finish(fin);
            if (!h) return std::unexpected(std::move(h.error()));
            h->deck_text = a.deck ? a.deck->text : std::string();
            h->run = r->sweep.run;
            if (auto w = written(out.adapt_history, io::save_adapt_history(out.adapt_history, *h)); !w) {
                return std::unexpected(w.error());
            }
        }
        report.adapt = solve::adapt_report(*r, o, fin.device, a.certify_required);
        if (a.predictor) {
            report.adapt->predictor = a.predictor_path.generic_string();
            report.adapt->predictor_digest = a.predictor_digest;
            report.adapt->predictor_status = predictor_status;
            report.adapt->predictor_out_of_range = predictor_out;
        }
        if (!out.device.empty()) {
            if (auto w = written(out.device, io::save_device(out.device, fin.device)); !w) return std::unexpected(w.error());
        }
        if (!out.mesh.empty()) {
            auto m = fin.generated ? io::save_generated_mesh(out.mesh, *fin.generated) : io::save_mesh(out.mesh, fin.device.mesh());
            if (auto w = written(out.mesh, std::move(m)); !w) return std::unexpected(w.error());
        }
        if (!out.sweep.empty()) {
            if (auto w = written(out.sweep, io::save_sweep(out.sweep, r->sweep)); !w) return std::unexpected(w.error());
        }
        if (!out.table.empty()) {
            if (auto w = written(out.table, write_table(out.table, fin.device, r->sweep, &*report.adapt)); !w) {
                return std::unexpected(w.error());
            }
        }
        if (!out.checkpoint.empty() && !r->sweep.points.empty()) {
            auto c = io::save_checkpoint(out.checkpoint, fin.device.mesh(), r->sweep.run, r->sweep.points.back());
            if (auto w = written(out.checkpoint, std::move(c)); !w) return std::unexpected(w.error());
        }
        if (!out.adapt.empty()) {
            if (auto w = written(out.adapt, io::save_adapt_report(out.adapt, *report.adapt)); !w) return std::unexpected(w.error());
        }
        report.sweep = std::move(r->sweep);
        report.adapted = std::move(r->mesh);
        return report;
    }
    if (!out.device.empty()) {
        if (auto r = written(out.device, io::save_device(out.device, plan.device)); !r) return std::unexpected(r.error());
    }
    if (!out.mesh.empty()) {
        auto r = plan.generated ? io::save_generated_mesh(out.mesh, *plan.generated)
                                : io::save_mesh(out.mesh, plan.device.mesh());
        if (auto w = written(out.mesh, std::move(r)); !w) return std::unexpected(w.error());
    }

    if (plan.analysis == Analysis::equilibrium) {
        const solve::EquilibriumOptions options{.newton = plan.options.newton,
                                                .linear = plan.options.linear,
                                                .Ns_override = plan.options.Ns_override,
                                                .models = plan.options.models};
        auto eq = solve::solve_equilibrium(plan.device, options, control);
        if (!eq) return std::unexpected(std::move(eq.error()));
        if (!out.checkpoint.empty()) {
            const std::vector<std::vector<double>> zero{std::vector<double>(plan.device.contacts().size(), 0.0)};
            results::BiasPoint point;
            point.bias_V = zero.front();
            point.fields = eq->fields;
            auto record = solve::make_run_record(plan.device, plan.options, zero);
            stamp(record);
            auto r = io::save_checkpoint(out.checkpoint, plan.device.mesh(), record, point);
            if (auto w = written(out.checkpoint, std::move(r)); !w) return std::unexpected(w.error());
        }
        report.equilibrium = std::move(*eq);
        return report;
    }

    // `initial` (Unit 61, RD2): the checkpoint's fields are the first state's starting guess. Its bias
    // need not be the run's (RQ3); it must belong to this mesh (load_checkpoint).
    std::optional<io::Checkpoint> initial;
    if (!plan.initial.empty()) {
        auto c = io::load_checkpoint(plan.initial, plan.device.mesh());
        if (!c) return std::unexpected(std::move(c.error()));
        if (c->kind == io::CheckpointKind::transient) {
            return std::unexpected(invalid("the initial checkpoint '" + plan.initial.string() +
                                           "' holds a transient's resume state: use restart to continue that run"));
        }
        initial = std::move(*c);
    }
    const results::NodeFields* guess = initial ? &initial->fields : nullptr;

    // Unit 43: the analysis's own options, with the plan's models and solver as their `steady`.
    if (plan.analysis == Analysis::transient) {
        solve::TransientOptions o = plan.transient;
        o.steady = plan.options;
        // A restart (Unit 61, RD4): the run continues from the checkpoint's resume state, which
        // solve_transient checks belongs to this run (its identity).
        std::optional<io::Checkpoint> resume;
        if (!plan.restart.empty()) {
            auto c = io::load_checkpoint(plan.restart, plan.device.mesh());
            if (!c) return std::unexpected(std::move(c.error()));
            if (c->kind != io::CheckpointKind::transient || !c->transient) {
                return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                               "' holds a steady state, not a transient's resume state: use initial "
                                               "to start a transient from it"));
            }
            if (auto e = restart_confinement(plan.restart, *c, o.steady.models)) return std::unexpected(std::move(*e));
            resume = std::move(*c);
        }
        auto t = solve::solve_transient(plan.device, plan.waveforms, o, guess, control,
                                        resume ? &*resume->transient : nullptr);
        if (!t) return std::unexpected(std::move(t.error()));
        stamp(t->run);
        if (!out.checkpoint.empty() && t->resume) {
            auto r = io::save_transient_checkpoint(out.checkpoint, plan.device.mesh(), *t);
            if (auto w = written(out.checkpoint, std::move(r)); !w) return std::unexpected(w.error());
        }
        if (!out.transient.empty()) {
            if (auto w = written(out.transient, io::save_transient(out.transient, *t)); !w) return std::unexpected(w.error());
        }
        if (!out.table.empty()) {
            if (auto w = written(out.table, write_transient_table(out.table, plan.device, *t)); !w) return std::unexpected(w.error());
        }
        report.transient = std::move(*t);
        return report;
    }
    if (plan.analysis == Analysis::small_signal) {
        solve::SmallSignalOptions o = plan.small_signal;
        o.steady = plan.options;
        // A restart (Unit 61, RD6): as the sweep's (DK9), after the operating point the checkpoint
        // holds, from its state, at every frequency of the remaining points.
        std::vector<std::vector<double>> points = plan.points;
        std::optional<io::Checkpoint> restart;
        if (!plan.restart.empty()) {
            auto c = io::load_checkpoint(plan.restart, plan.device.mesh());
            if (!c) return std::unexpected(std::move(c.error()));
            if (c->kind == io::CheckpointKind::transient) {
                return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                               "' holds a transient's resume state, not an operating point"));
            }
            if (auto e = restart_confinement(plan.restart, *c, o.steady.models)) return std::unexpected(std::move(*e));
            std::size_t k = 0;
            while (k < points.size() && !same_bias(points[k], c->bias_V)) ++k;
            if (k == points.size()) {
                return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                               "' holds a bias that is not an operating point of the run"));
            }
            if (k + 1 == points.size()) {
                return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                               "' holds the run's last operating point: nothing is left to solve"));
            }
            report.skipped_points = k + 1;
            points.erase(points.begin(), points.begin() + static_cast<std::ptrdiff_t>(k + 1));
            if (!o.steady.thermal_bias_K.empty()) {
                o.steady.thermal_bias_K.erase(o.steady.thermal_bias_K.begin(),
                                              o.steady.thermal_bias_K.begin() + static_cast<std::ptrdiff_t>(k + 1));
            }
            o.steady.step_control.initial_bias_V = c->bias_V;
            restart = std::move(*c);
        }
        auto s = solve::solve_small_signal(plan.device, points, o, restart ? &restart->fields : guess, control);
        if (!s) return std::unexpected(std::move(s.error()));
        stamp(s->run);
        if (!out.small_signal.empty()) {
            if (auto w = written(out.small_signal, io::save_small_signal(out.small_signal, *s)); !w) return std::unexpected(w.error());
        }
        if (!out.table.empty()) {
            if (auto w = written(out.table, write_small_signal_table(out.table, plan.device, *s)); !w) return std::unexpected(w.error());
        }
        if (!out.checkpoint.empty() && !s->points.empty()) {
            auto r = io::save_checkpoint(out.checkpoint, plan.device.mesh(), s->run, s->points.back().dc);
            if (auto w = written(out.checkpoint, std::move(r)); !w) return std::unexpected(w.error());
        }
        report.small_signal = std::move(*s);
        return report;
    }
    if (plan.analysis == Analysis::trace) {
        solve::TraceOptions o = plan.trace;
        o.steady = plan.options;
        // A restart (Unit 61, RD5): the trace continues from the checkpoint's resume state.
        std::optional<io::Checkpoint> resume;
        if (!plan.restart.empty()) {
            auto c = io::load_checkpoint(plan.restart, plan.device.mesh());
            if (!c) return std::unexpected(std::move(c.error()));
            if (c->kind != io::CheckpointKind::trace || !c->trace) {
                return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                               "' holds no trace's resume state: use initial to start a trace "
                                               "from its fields"));
            }
            if (auto e = restart_confinement(plan.restart, *c, o.steady.models)) return std::unexpected(std::move(*e));
            resume = std::move(*c);
        }
        auto s = solve::trace_bias(plan.device, o, guess, control, resume ? &*resume->trace : nullptr);
        if (!s) return std::unexpected(std::move(s.error()));
        stamp(s->run);
        if (!out.sweep.empty()) {
            if (auto w = written(out.sweep, io::save_sweep(out.sweep, *s)); !w) return std::unexpected(w.error());
        }
        if (!out.table.empty()) {
            if (auto w = written(out.table, write_table(out.table, plan.device, *s)); !w) return std::unexpected(w.error());
        }
        if (!out.checkpoint.empty() && !s->points.empty()) {
            // With the resume state (Unit 61, RD5) when the trace got past its first tangent.
            auto r = s->resume ? io::save_trace_checkpoint(out.checkpoint, plan.device.mesh(), *s)
                               : io::save_checkpoint(out.checkpoint, plan.device.mesh(), s->run, s->points.back());
            if (auto w = written(out.checkpoint, std::move(r)); !w) return std::unexpected(w.error());
        }
        report.sweep = std::move(*s);
        return report;
    }

    solve::BiasOptions options = plan.options;
    std::vector<std::vector<double>> points = plan.points;
    std::optional<io::Checkpoint> restart;
    if (!plan.restart.empty()) {
        auto c = io::load_checkpoint(plan.restart, plan.device.mesh());
        if (!c) return std::unexpected(std::move(c.error()));
        if (c->kind == io::CheckpointKind::transient) {
            return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                           "' holds a transient's resume state, not a steady point of the sweep"));
        }
        if (auto e = restart_confinement(plan.restart, *c, options.models)) return std::unexpected(std::move(*e));
        std::size_t k = 0;
        while (k < points.size() && !same_bias(points[k], c->bias_V)) ++k;
        if (k == points.size()) {
            return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                           "' holds a bias that is not a point of the sweep"));
        }
        if (k + 1 == points.size()) {
            return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                           "' holds the sweep's last point: nothing is left to solve"));
        }
        report.skipped_points = k + 1;
        points.erase(points.begin(), points.begin() + static_cast<std::ptrdiff_t>(k + 1));
        if (!options.thermal_bias_K.empty()) {
            options.thermal_bias_K.erase(options.thermal_bias_K.begin(),
                                         options.thermal_bias_K.begin() + static_cast<std::ptrdiff_t>(k + 1));
        }
        options.step_control.initial_bias_V = c->bias_V;
        restart = std::move(*c);
    }
    auto sweep = solve::sweep_bias(plan.device, points, options, restart ? &restart->fields : guess, control);
    if (!sweep) return std::unexpected(std::move(sweep.error()));
    stamp(sweep->run);
    if (!out.sweep.empty()) {
        if (auto w = written(out.sweep, io::save_sweep(out.sweep, *sweep)); !w) return std::unexpected(w.error());
    }
    if (!out.table.empty()) {
        if (auto w = written(out.table, write_table(out.table, plan.device, *sweep)); !w) return std::unexpected(w.error());
    }
    if (!out.checkpoint.empty() && !sweep->points.empty()) {
        auto r = io::save_checkpoint(out.checkpoint, plan.device.mesh(), sweep->run, sweep->points.back());
        if (auto w = written(out.checkpoint, std::move(r)); !w) return std::unexpected(w.error());
    }
    report.sweep = std::move(*sweep);
    return report;
}

}  // namespace NiTCAD::deck
