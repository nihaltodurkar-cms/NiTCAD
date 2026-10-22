// The results workspace (workspace.hpp).
#include "NiTCAD/app/workspace.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <format>
#include <utility>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::app {

namespace {

namespace fs = std::filesystem;

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

// The drawing geometry of the device: its tensor cells, a generated mesh's simplices, or the graph.
std::expected<void, base::Error> make_geometry(Workspace& w) {
    const mesh::Mesh& m = w.device->mesh();
    if (w.device->cells()) {
        auto g = render::Geometry::tensor(m, *w.device->cells());
        if (!g) return std::unexpected(std::move(g.error()));
        w.geometry = std::make_unique<render::Geometry>(std::move(*g));
    } else if (w.generated) {
        auto g = render::Geometry::simplex(m, w.generated->simplices);
        if (!g) return std::unexpected(std::move(g.error()));
        w.geometry = std::make_unique<render::Geometry>(std::move(*g));
    } else {
        w.geometry = std::make_unique<render::Geometry>(render::Geometry::graph(m));
    }
    return {};
}

// The bias points the views read, as nitcad plot lists them (src/cli/plot.cpp points_of).
void collect_points(Workspace& w) {
    w.points.clear();
    if (w.sweep) w.points = w.sweep->points;
    if (w.transient) {
        for (const auto& s : w.transient->snapshots) {
            results::BiasPoint p;
            p.fields = s.fields;
            p.bands = s.bands;
            p.edge_current_n = s.edge_current_n;
            p.edge_current_p = s.edge_current_p;
            p.rates = s.rates;
            if (s.point < w.transient->points.size()) p.bias_V = w.transient->points[s.point].bias_V;
            w.points.push_back(std::move(p));
        }
    }
    if (w.small_signal) {
        for (const auto& p : w.small_signal->points) w.points.push_back(p.dc);
    }
    if (w.equilibrium) w.points.push_back(*w.equilibrium);
}

// The results belong to the device's mesh: by the mesh digest when they carry one (Unit 45), and every
// point's fields one per node.
std::expected<void, base::Error> check_pairing(const Workspace& w) {
    const mesh::Mesh& m = w.device->mesh();
    const results::RunRecord* r = w.record();
    if (r && r->mesh_digest != 0 && r->mesh_digest != mesh::graph_digest(m)) {
        return std::unexpected(invalid("the result was computed on another mesh: its mesh digest is not the device's"));
    }
    for (const auto& p : w.points) {
        if (p.fields.potential_V.size() != m.node_count()) {
            return std::unexpected(invalid("the result has " + std::to_string(p.fields.potential_V.size()) +
                                           " values per field and the device's mesh " + std::to_string(m.node_count()) +
                                           " nodes: they do not belong together"));
        }
    }
    return {};
}

std::expected<void, base::Error> load_result_file(const fs::path& path, Workspace& w) {
    const std::string ext = path.extension().string();
    if (ext == ".nitsweep") {
        auto s = io::load_sweep(path);
        if (!s) return std::unexpected(std::move(s.error()));
        w.sweep = std::move(*s);
    } else if (ext == ".nittran") {
        auto t = io::load_transient(path);
        if (!t) return std::unexpected(std::move(t.error()));
        w.transient = std::move(*t);
    } else if (ext == ".nitac") {
        auto a = io::load_small_signal(path);
        if (!a) return std::unexpected(std::move(a.error()));
        w.small_signal = std::move(*a);
    } else {
        return std::unexpected(invalid("a result file is .nitsweep, .nittran or .nitac: " + path.string()));
    }
    return {};
}

results::BiasPoint equilibrium_point(const results::EquilibriumResult& e) {
    results::BiasPoint p;
    p.fields = e.fields;
    p.bands = e.bands;
    p.gate_charge = e.gate_charge;
    p.interface_trap_charge = e.interface_trap_charge;
    p.convergence = e.convergence;
    return p;
}

std::string hex(std::uint64_t v) { return std::format("{:016x}", v); }

std::string number(double v) { return std::format("{:.6g}", v); }

}  // namespace

const results::RunRecord* Workspace::record() const noexcept {
    if (sweep) return &sweep->run;
    if (transient) return &transient->run;
    if (small_signal) return &small_signal->run;
    return nullptr;
}

std::expected<Workspace, base::Error> workspace_from_run(deck::Plan plan, deck::RunReport report) {
    Workspace w;
    w.analysis = plan.analysis;
    w.device = std::make_unique<device::Device>(std::move(plan.device));
    if (plan.generated) w.generated = std::make_unique<mesh::GeneratedMesh>(std::move(*plan.generated));
    w.sweep = std::move(report.sweep);
    w.transient = std::move(report.transient);
    w.small_signal = std::move(report.small_signal);
    if (report.equilibrium && !w.sweep && !w.transient && !w.small_signal) {
        w.equilibrium = equilibrium_point(*report.equilibrium);
    }
    if (auto g = make_geometry(w); !g) return std::unexpected(std::move(g.error()));
    collect_points(w);
    if (auto c = check_pairing(w); !c) return std::unexpected(std::move(c.error()));
    return w;
}

std::expected<Workspace, base::Error> workspace_from_deck(const fs::path& deck_path) {
    auto parsed = deck::read(deck_path);
    if (!parsed) return std::unexpected(std::move(parsed.error()));
    auto plan = deck::build(*parsed);
    if (!plan) return std::unexpected(std::move(plan.error()));
    return workspace_from_plan(std::move(*plan), deck_path);
}

std::expected<Workspace, base::Error> workspace_from_plan(deck::Plan plan, const fs::path& deck_path) {
    Workspace w;
    w.source = deck_path;
    w.analysis = plan.analysis;
    const deck::Outputs out = plan.outputs;
    const deck::Analysis analysis = plan.analysis;
    w.device = std::make_unique<device::Device>(std::move(plan.device));
    if (plan.generated) w.generated = std::make_unique<mesh::GeneratedMesh>(std::move(*plan.generated));
    if (auto g = make_geometry(w); !g) return std::unexpected(std::move(g.error()));
    // The results the deck's outputs hold, if it has run.
    using deck::Analysis;
    const fs::path file = analysis == Analysis::transient      ? out.transient
                          : analysis == Analysis::small_signal ? out.small_signal
                          : analysis == Analysis::equilibrium  ? out.checkpoint
                                                               : out.sweep;
    if (!file.empty() && fs::exists(file)) {
        if (analysis == Analysis::equilibrium) {
            auto c = io::load_checkpoint(file, w.device->mesh());
            if (!c) return std::unexpected(std::move(c.error()));
            results::BiasPoint p;
            p.fields = c->fields;
            w.equilibrium = std::move(p);
        } else if (auto r = load_result_file(file, w); !r) {
            return std::unexpected(std::move(r.error()));
        }
    }
    collect_points(w);
    if (auto c = check_pairing(w); !c) return std::unexpected(std::move(c.error()));
    return w;
}

std::expected<Workspace, base::Error> workspace_from_result(const fs::path& result, const fs::path& device,
                                                           const fs::path& mesh_file) {
    Workspace w;
    w.source = result;
    auto d = io::load_device(device);
    if (!d) return std::unexpected(std::move(d.error()));
    w.device = std::make_unique<device::Device>(std::move(*d));
    if (!mesh_file.empty()) {
        auto g = io::load_generated_mesh(mesh_file);
        if (!g) return std::unexpected(std::move(g.error()));
        if (mesh::graph_digest(g->mesh) != mesh::graph_digest(w.device->mesh())) {
            return std::unexpected(invalid("the mesh file " + mesh_file.string() + " is not the device's mesh"));
        }
        w.generated = std::make_unique<mesh::GeneratedMesh>(std::move(*g));
    }
    if (auto r = load_result_file(result, w); !r) return std::unexpected(std::move(r.error()));
    if (auto g = make_geometry(w); !g) return std::unexpected(std::move(g.error()));
    collect_points(w);
    if (auto c = check_pairing(w); !c) return std::unexpected(std::move(c.error()));
    return w;
}

std::vector<Line> run_record(const Workspace& w) {
    std::vector<Line> out;
    const results::RunRecord* r = w.record();
    const mesh::Mesh& m = w.device->mesh();
    out.push_back({"Mesh", std::format("{}D, {} nodes, {} edges", m.dimension(), m.node_count(), m.edges().size())});
    if (w.sweep) out.push_back({"Points", std::to_string(w.sweep->points.size())});
    if (w.transient) {
        out.push_back({"Steps", std::to_string(w.transient->points.empty() ? 0 : w.transient->points.size() - 1)});
        out.push_back({"Snapshots", std::to_string(w.transient->snapshots.size())});
        out.push_back({"Rejected steps", std::to_string(w.transient->rejected_steps)});
    }
    if (w.small_signal) {
        out.push_back({"Operating points", std::to_string(w.small_signal->points.size())});
        out.push_back({"Frequencies", std::to_string(w.small_signal->frequency_Hz.size())});
    }
    if (w.equilibrium) {
        out.push_back({"Equilibrium", w.equilibrium->convergence.converged
                                          ? std::format("converged in {} iterations", w.equilibrium->convergence.iterations.size())
                                          : std::string("not converged")});
    }
    if (!r) return out;
    const std::optional<base::Error>& stopped = w.sweep          ? w.sweep->stopped
                                                : w.transient    ? w.transient->stopped
                                                                 : w.small_signal->stopped;
    out.push_back({"Ended", stopped ? "stopped early: " + stopped->message : std::string("completed")});
    out.push_back({"Input identity", hex(r->input_identity)});
    out.push_back({"Mesh digest", r->mesh_digest != 0 ? hex(r->mesh_digest) : std::string("not recorded")});
    out.push_back({"Settings", std::to_string(r->settings.size())});
    const results::Provenance& p = r->provenance;
    out.push_back({"Version", p.version});
    out.push_back({"Revision", p.revision + (p.dirty ? " (modified)" : "")});
    out.push_back({"Build", p.build_type + ", " + p.compiler + ", Eigen " + p.eigen_version});
    if (!p.mkl_version.empty()) out.push_back({"MKL", p.mkl_version});
    out.push_back({"Platform", p.os + ", " + p.cpu + ", " + std::to_string(p.logical_cpus) + " logical CPUs"});
    out.push_back({"Started", p.started_utc});
    out.push_back({"Wall time", number(p.wall_s) + " s"});
    for (const auto& [name, value] : r->settings) out.push_back({"  " + name, number(value)});
    return out;
}

std::vector<PointHistory> convergence_history(const Workspace& w) {
    std::vector<PointHistory> out;
    if (w.sweep) {
        for (std::size_t k = 0; k < w.sweep->points.size(); ++k) {
            const auto& c = w.sweep->points[k].convergence;
            out.push_back({k, c.iterations, c.converged});
        }
        if (w.sweep->unfinished) out.push_back({w.sweep->points.size(), w.sweep->unfinished->iterations, false});
    } else if (w.small_signal) {
        for (std::size_t k = 0; k < w.small_signal->points.size(); ++k) {
            const auto& c = w.small_signal->points[k].dc.convergence;
            out.push_back({k, c.iterations, c.converged});
        }
    } else if (w.transient) {
        for (std::size_t k = 0; k < w.transient->points.size(); ++k) {
            const auto& c = w.transient->points[k].convergence;
            out.push_back({k, c.iterations, c.converged});
        }
    } else if (w.equilibrium) {
        out.push_back({0, w.equilibrium->convergence.iterations, w.equilibrium->convergence.converged});
    }
    return out;
}

std::size_t swept_contact(const std::vector<std::vector<double>>& biases) {
    std::size_t best = 0;
    double spread = -1.0;
    for (std::size_t k = 0; !biases.empty() && k < biases.front().size(); ++k) {
        double lo = 1e300, hi = -1e300;
        for (const auto& b : biases) {
            lo = std::min(lo, b[k]);
            hi = std::max(hi, b[k]);
        }
        if (hi - lo > spread) {
            spread = hi - lo;
            best = k;
        }
    }
    return best;
}

std::vector<Figure> extracted_figures(const Workspace& w) {
    std::vector<Figure> out;
    if (!w.sweep || w.sweep->points.size() < 2) return out;
    const auto contacts = w.device->contacts();
    std::vector<std::vector<double>> biases;
    for (const auto& p : w.sweep->points) biases.push_back(p.bias_V);
    const std::size_t swept = swept_contact(biases);
    std::string curve_name;
    const auto add = [&](std::string name, const std::expected<analysis::Extraction, base::Error>& e, std::string unit) {
        if (!e) return;
        out.push_back({std::move(name), curve_name, e->value, unit.empty() ? std::string(e->unit) : std::move(unit),
                       e->extrapolated, e->below_resolution});
    };
    const int dim = w.device->mesh().dimension();
    const std::string current_unit = dim == 1 ? "A/cm\xC2\xB2" : dim == 2 ? "A/cm" : "A";
    if (contacts[swept].kind != device::ContactKind::ohmic) {
        // A transfer curve: the contact carrying the most current at the last point is sensed (the swept
        // gate or electrode carries none in a steady state).
        std::size_t sense = swept;
        double best = -1.0;
        const auto& last = w.sweep->points.back().terminal_current;
        for (std::size_t k = 0; k < contacts.size() && k < last.size(); ++k) {
            if (std::abs(last[k]) > best) {
                best = std::abs(last[k]);
                sense = k;
            }
        }
        if (sense == swept) return out;
        auto curve = analysis::current_curve(*w.sweep, swept, sense);
        if (!curve) return out;
        // An n-channel current rises with the gate; a p-channel curve is read negated.
        const bool negative = curve->y.back() < 0.0;
        const analysis::Curve transfer = negative ? analysis::scaled(*curve, -1.0) : *curve;
        const double drain_V = w.sweep->points.back().bias_V[sense];
        curve_name = "I(" + contacts[sense].name + ") vs V(" + contacts[swept].name + ")";
        add("Threshold (max gm)", analysis::threshold_max_gm(transfer, drain_V), "V");
        add("Subthreshold swing", analysis::subthreshold_swing(transfer), "mV/dec");
        add("Max transconductance", analysis::max_transconductance(transfer), current_unit + "/V");
        add("On/off ratio", analysis::on_off_ratio(transfer, transfer.x.back(), transfer.x.front()), "");
    } else {
        auto curve = analysis::current_curve(*w.sweep, swept, swept);
        if (!curve) return out;
        curve_name = "I(" + contacts[swept].name + ") vs V(" + contacts[swept].name + ")";
        add("Current at the last bias",
            analysis::value_at(*curve, curve->x.back(), analysis::Scale::linear), current_unit);
    }
    return out;
}

}  // namespace NiTCAD::app
