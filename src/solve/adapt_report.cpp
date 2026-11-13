// The adaptation's settings, run identity and report (Unit 17, slice 6; DECISIONS.md AD8, AD13, AR8).
#include <bit>
#include <span>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "NiTCAD/solve/adapt.hpp"

namespace NiTCAD::solve {

namespace {

// 64-bit FNV-1a, as the run records' digest (run_record.cpp).
class Digest {
public:
    explicit Digest(std::uint64_t h) : h_(h) {}
    void bytes(const void* data, std::size_t size) {
        const auto* p = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < size; ++i) {
            h_ ^= p[i];
            h_ *= 0x100000001b3ull;
        }
    }
    void integer(std::uint64_t v) { bytes(&v, sizeof v); }
    void text(std::string_view s) {
        integer(s.size());
        bytes(s.data(), s.size());
    }
    [[nodiscard]] std::uint64_t value() const noexcept { return h_; }

private:
    std::uint64_t h_;
};


std::vector<std::uint64_t> counts(std::span<const std::size_t> v) {
    return {v.begin(), v.end()};
}

}  // namespace

std::vector<std::pair<std::string, double>> adapt_settings(const AdaptOptions& o, bool certify_required) {
    std::vector<std::pair<std::string, double>> s;
    const auto put = [&s](std::string name, double v) { s.emplace_back(std::move(name), v); };
    put("adapt.tolerance", o.tolerance);
    for (const AdaptGoal& g : o.goals) put("adapt.goal", static_cast<double>(g.contact));
    for (const std::size_t k : o.at) put("adapt.at", static_cast<double>(k));
    put("adapt.weighting", o.weighting == Weighting::adjoint ? 1.0 : 0.0);
    put("adapt.policy", o.policy == Policy::intelligent ? 1.0 : 0.0);
    put("adapt.theta", o.theta);
    put("adapt.grading", o.grading);
    put("adapt.min_size", o.min_size);
    put("adapt.max_passes", static_cast<double>(o.max_passes));
    put("adapt.max_nodes", static_cast<double>(o.max_nodes));
    put("adapt.verify", o.verify ? 1.0 : 0.0);
    put("adapt.certify_required", certify_required ? 1.0 : 0.0);
    for (std::size_t k = 0; k < indicator_count; ++k) {
        put(std::string("adapt.indicators.") + to_string(static_cast<Indicator>(k)), o.indicators.enabled[k] ? 1.0 : 0.0);
    }
    put("adapt.indicators.potential_vt", o.indicators.potential_vt);
    put("adapt.indicators.density_decades", o.indicators.density_decades);
    put("adapt.indicators.doping_decades", o.indicators.doping_decades);
    put("adapt.indicators.generation_share", o.indicators.generation_share);
    put("adapt.max_points", static_cast<double>(o.intelligent.max_points));
    put("adapt.miss_margin", o.intelligent.miss_margin);
    put("adapt.budget_dofs", static_cast<double>(o.intelligent.dof_budget));
    put("adapt.budget_seconds", o.intelligent.seconds_budget);
    return s;
}

void stamp_adapt(results::RunRecord& run, const AdaptOptions& o, bool certify_required) {
    Digest d(run.input_identity);
    d.text("adapt");
    for (auto& [name, value] : adapt_settings(o, certify_required)) {
        d.text(name);
        d.integer(std::bit_cast<std::uint64_t>(value));
        run.settings.emplace_back(std::move(name), value);
    }
    d.integer(run.mesh_digest);
    run.input_identity = d.value();
}

results::AdaptReport adapt_report(const AdaptResult& r, const AdaptOptions& o, const device::Device& device,
                                  bool certify_required) {
    results::AdaptReport out;
    out.run = r.sweep.run;
    for (const AdaptGoal& g : o.goals) {
        out.goal_contacts.push_back(g.contact < device.contacts().size() ? device.contacts()[g.contact].name : "?");
    }
    out.tolerance = o.tolerance;
    out.points.assign(r.points.begin(), r.points.end());
    out.policy = to_string(r.policy);
    out.weighting = o.weighting == Weighting::adjoint ? "adjoint" : "none";
    out.verify = o.verify;
    out.certify_required = certify_required;
    for (std::size_t k = 0; k < indicator_count; ++k) out.indicators.emplace_back(to_string(static_cast<Indicator>(k)));
    out.stop = to_string(r.stop);
    out.rejected_stops = r.rejected_stops;
    out.fallback = to_string(r.fallback);
    out.fallback_pass = r.fallback_pass;
    out.misses = r.misses;
    out.junction_points = r.junction_points;
    out.interface_points = r.interface_points;
    out.predicted_nodes = r.predicted_nodes;
    out.prediction_rounds = r.prediction_rounds;
    out.coarse_nodes = r.coarse_nodes;
    out.coarse_seconds = r.coarse_seconds;
    out.verification_nodes = r.verification_nodes;
    out.verification_seconds = r.verification_seconds;
    out.verification_goal = r.verification_goal;
    for (const AdaptPass& p : r.passes) {
        results::AdaptPassRow row;
        row.nodes = p.nodes;
        row.edges = p.edges;
        row.dofs = p.dofs;
        row.mesh_digest = p.mesh_digest;
        row.source = to_string(p.source);
        row.goal = p.goal;
        row.resolution = p.resolution;
        row.change = p.change;
        row.estimate = p.estimate;
        row.relative_estimate = p.relative_estimate;
        row.estimate_nodes = p.estimate_nodes;
        row.marked = p.marking.count;
        row.marked_by_indicator = counts(p.marking.by_kind);
        row.marked_by_region = counts(p.marking.by_region);
        row.nodes_with_excess = counts(p.nodes_with_excess);
        row.newton_iterations = p.newton_iterations;
        row.matrix_nonzeros = p.matrix_nonzeros;
        row.factor_fill = p.factor_fill;
        row.rejected = p.rejected;
        row.miss = p.miss;
        row.predicted_seconds = p.predicted_seconds;
        row.predicted_seconds_nodes = p.predicted_seconds_nodes;
        row.candidate_count = p.candidate_count;
        row.candidates_taken = p.candidates_taken;
        row.predicted_relative = p.predicted_relative;
        row.predicted_dofs = p.predicted_dofs;
        for (const AdaptCandidate& c : p.candidates) {
            row.candidates.push_back({.axis = c.axis,
                                      .index = c.index,
                                      .x = c.location[0],
                                      .y = c.location[1],
                                      .z = c.location[2],
                                      .region = c.region,
                                      .indicator = to_string(c.indicator),
                                      .order = c.order,
                                      .levels = c.levels,
                                      .benefit = c.benefit,
                                      .cost = c.cost,
                                      .work = c.work,
                                      .added_nodes = c.added_nodes,
                                      .added_dofs = c.added_dofs,
                                      .added_nonzeros = c.added_nonzeros,
                                      .realized_ratio = c.realized_ratio,
                                      .realized = c.realized});
        }
        for (const OrderFit& f : p.orders) {
            row.orders.push_back({.indicator = to_string(f.kind), .region = f.region, .order = f.order, .samples = f.samples});
        }
        const CostModel& m = p.cost_model;
        row.cost_a = m.a;
        row.cost_b = m.b;
        row.cost_c = m.c;
        row.cost_d = m.d;
        row.cost_mesh_per_node = m.mesh_per_node;
        row.cost_transfer_per_node = m.transfer_per_node;
        row.cost_residual = m.residual;
        row.cost_alpha = m.alpha;
        row.cost_beta = m.beta;
        row.cost_residual_nodes = m.residual_nodes;
        row.cost_samples = m.samples;
        row.seconds_mesh = p.seconds_mesh;
        row.seconds_transfer = p.seconds_transfer;
        row.seconds_solve = p.seconds_solve;
        row.seconds_estimate = p.seconds_estimate;
        row.seconds_estimate_factorize = p.seconds_estimate_factorize;
        row.seconds_transposed = p.seconds_transposed;
        out.passes.push_back(std::move(row));
    }
    bool any_not = false, any_unverified = false;
    for (const GoalCertificate& c : r.certificates) {
        out.certificates.push_back({.goal = c.goal,
                                    .point = c.point,
                                    .value = c.value,
                                    .estimate = c.estimate,
                                    .verified_error = c.verified_error,
                                    .order = c.order,
                                    .allowed = c.allowed,
                                    .status = to_string(c.status)});
        any_not = any_not || c.status == CertificateStatus::not_certified;
        any_unverified = any_unverified || c.status == CertificateStatus::unverified;
    }
    out.status = any_not ? "not_certified" : any_unverified || r.certificates.empty() ? "unverified" : "certified";
    out.final_mesh_digest = r.passes.empty() ? 0 : r.passes.back().mesh_digest;
    return out;
}

}  // namespace NiTCAD::solve
