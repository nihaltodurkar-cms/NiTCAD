// The adaptation report on disk (Unit 17, slice 6; io/adapt_report.hpp). Keys are local to this content kind.
#include "NiTCAD/io/adapt_report.hpp"

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "container.hpp"
#include "record.hpp"

namespace NiTCAD::io {

namespace {

using detail::corrupt;
using detail::Reader;
using detail::Writer;

namespace key {
// Report.
constexpr std::uint16_t run = 1, goal_contact = 2, tolerance = 3, points = 4, policy = 5, weighting = 6, verify = 7,
                        certify = 8, indicator = 9, stop = 10, rejected_stops = 11, fallback = 12, fallback_pass = 13,
                        misses = 14, prediction = 15, coarse_nodes = 16, coarse_seconds = 17, verification_nodes = 18,
                        verification_seconds = 19, verification_goal = 20, pass = 21, certificate = 22, status = 23,
                        final_mesh = 24;
// Pass.
constexpr std::uint16_t sizes = 1, mesh_digest = 2, source = 3, goal = 4, resolution = 5, change = 6, estimate = 7,
                        relative_estimate = 8, estimate_nodes = 9, marked = 10, by_indicator = 11, by_region = 12,
                        with_excess = 13, newton = 14, flags = 15, predicted = 16, allocation = 17, candidate = 18,
                        order = 19, cost = 20, seconds = 21;
// Candidate.
constexpr std::uint16_t c_ints = 1, c_reals = 2, c_indicator = 3;
// Order.
constexpr std::uint16_t o_indicator = 1, o_ints = 2, o_order = 3;
// Certificate.
constexpr std::uint16_t g_ints = 1, g_reals = 2, g_status = 3;
}  // namespace key

std::vector<std::int64_t> ints(std::initializer_list<std::uint64_t> v) {
    std::vector<std::int64_t> out;
    for (const std::uint64_t x : v) out.push_back(static_cast<std::int64_t>(x));
    return out;
}

std::vector<std::int64_t> ints(const std::vector<std::uint64_t>& v) {
    std::vector<std::int64_t> out;
    for (const std::uint64_t x : v) out.push_back(static_cast<std::int64_t>(x));
    return out;
}

std::vector<std::uint64_t> unsigned_of(const std::vector<std::int64_t>& v) {
    std::vector<std::uint64_t> out;
    for (const std::int64_t x : v) out.push_back(static_cast<std::uint64_t>(x));
    return out;
}

// A list of texts, each a one-field nested record.
void write_texts(Writer& w, std::uint16_t k, const std::vector<std::string>& texts) {
    for (const std::string& t : texts) {
        Writer one;
        one.text(1, t);
        w.record(k, one);
    }
}

std::expected<std::vector<std::string>, base::Error> read_texts(const Reader& r, std::uint16_t k) {
    NITCAD_IO_TRY(list, r.records(k));
    std::vector<std::string> out;
    for (const Reader& one : list) {
        NITCAD_IO_TRY(t, one.text(1));
        out.push_back(std::move(t));
    }
    return out;
}

// A fixed-length array, or corrupt.
std::expected<std::vector<std::int64_t>, base::Error> ints_of(const Reader& r, std::uint16_t k, std::size_t n,
                                                              const char* what) {
    NITCAD_IO_TRY(v, r.i64s(k));
    if (v.size() != n) return std::unexpected(corrupt(std::string("adaptation report ") + what));
    return v;
}

std::expected<std::vector<double>, base::Error> reals_of(const Reader& r, std::uint16_t k, std::size_t n,
                                                         const char* what) {
    NITCAD_IO_TRY(v, r.f64s(k));
    if (v.size() != n) return std::unexpected(corrupt(std::string("adaptation report ") + what));
    return v;
}

Writer write_pass(const results::AdaptPassRow& p) {
    namespace k = key;
    Writer w;
    w.i64s(k::sizes, ints({p.nodes, p.edges, p.dofs}));
    w.u64(k::mesh_digest, p.mesh_digest);
    w.text(k::source, p.source);
    w.f64s(k::goal, p.goal);
    w.f64s(k::resolution, p.resolution);
    w.f64s(k::change, p.change);
    w.f64s(k::estimate, p.estimate);
    w.f64(k::relative_estimate, p.relative_estimate);
    w.u64(k::estimate_nodes, p.estimate_nodes);
    w.u64(k::marked, p.marked);
    w.i64s(k::by_indicator, ints(p.marked_by_indicator));
    w.i64s(k::by_region, ints(p.marked_by_region));
    w.i64s(k::with_excess, ints(p.nodes_with_excess));
    w.i64s(k::newton, std::vector<std::int64_t>{p.newton_iterations, static_cast<std::int64_t>(p.matrix_nonzeros),
                                                static_cast<std::int64_t>(p.factor_fill)});
    w.i64s(k::flags, std::vector<std::int64_t>{p.rejected ? 1 : 0, p.miss ? 1 : 0});
    w.f64s(k::predicted, std::vector<double>{p.predicted_seconds, p.predicted_seconds_nodes, p.predicted_relative});
    w.i64s(k::allocation, ints({p.candidate_count, p.candidates_taken, p.predicted_dofs}));
    for (const results::AdaptCandidateRow& c : p.candidates) {
        Writer cw;
        cw.i64s(key::c_ints, std::vector<std::int64_t>{c.axis, static_cast<std::int64_t>(c.index),
                                                       static_cast<std::int64_t>(c.region), c.levels,
                                                       static_cast<std::int64_t>(c.added_nodes),
                                                       static_cast<std::int64_t>(c.added_dofs),
                                                       static_cast<std::int64_t>(c.added_nonzeros)});
        cw.f64s(key::c_reals, std::vector<double>{c.x, c.y, c.z, c.order, c.benefit, c.cost, c.work, c.realized_ratio,
                                                  c.realized});
        cw.text(key::c_indicator, c.indicator);
        w.record(k::candidate, cw);
    }
    for (const results::AdaptOrderRow& o : p.orders) {
        Writer ow;
        ow.text(key::o_indicator, o.indicator);
        ow.i64s(key::o_ints, ints({o.region, o.samples}));
        ow.f64(key::o_order, o.order);
        w.record(k::order, ow);
    }
    w.f64s(k::cost, std::vector<double>{p.cost_a, p.cost_b, p.cost_c, p.cost_d, p.cost_mesh_per_node,
                                        p.cost_transfer_per_node, p.cost_residual, p.cost_alpha, p.cost_beta,
                                        p.cost_residual_nodes, static_cast<double>(p.cost_samples)});
    w.f64s(k::seconds, std::vector<double>{p.seconds_mesh, p.seconds_transfer, p.seconds_solve, p.seconds_estimate});
    return w;
}

std::expected<results::AdaptPassRow, base::Error> read_pass(const Reader& r, std::size_t goals) {
    namespace k = key;
    results::AdaptPassRow p;
    NITCAD_IO_TRY(sz, ints_of(r, k::sizes, 3, "pass sizes"));
    p.nodes = static_cast<std::uint64_t>(sz[0]);
    p.edges = static_cast<std::uint64_t>(sz[1]);
    p.dofs = static_cast<std::uint64_t>(sz[2]);
    NITCAD_IO_TRY(digest, r.u64(k::mesh_digest));
    p.mesh_digest = digest;
    NITCAD_IO_TRY(source, r.text(k::source));
    p.source = std::move(source);
    NITCAD_IO_TRY(goal, r.f64s_or_empty(k::goal));
    NITCAD_IO_TRY(resolution, r.f64s_or_empty(k::resolution));
    NITCAD_IO_TRY(change, r.f64s_or_empty(k::change));
    NITCAD_IO_TRY(estimate, r.f64s_or_empty(k::estimate));
    if (goal.size() != goals || resolution.size() != goals || (!change.empty() && change.size() != goals) ||
        (!estimate.empty() && estimate.size() != goals)) {
        return std::unexpected(corrupt("adaptation report pass goals"));
    }
    p.goal = std::move(goal);
    p.resolution = std::move(resolution);
    p.change = std::move(change);
    p.estimate = std::move(estimate);
    NITCAD_IO_TRY(rel, r.f64(k::relative_estimate));
    p.relative_estimate = rel;
    NITCAD_IO_TRY(en, r.u64(k::estimate_nodes));
    p.estimate_nodes = en;
    NITCAD_IO_TRY(marked, r.u64(k::marked));
    p.marked = marked;
    NITCAD_IO_TRY(bi, r.i64s(k::by_indicator));
    NITCAD_IO_TRY(br, r.i64s(k::by_region));
    NITCAD_IO_TRY(we, r.i64s(k::with_excess));
    p.marked_by_indicator = unsigned_of(bi);
    p.marked_by_region = unsigned_of(br);
    p.nodes_with_excess = unsigned_of(we);
    NITCAD_IO_TRY(nw, ints_of(r, k::newton, 3, "pass Newton counts"));
    p.newton_iterations = nw[0];
    p.matrix_nonzeros = static_cast<std::uint64_t>(nw[1]);
    p.factor_fill = static_cast<std::uint64_t>(nw[2]);
    NITCAD_IO_TRY(fl, ints_of(r, k::flags, 2, "pass flags"));
    p.rejected = fl[0] != 0;
    p.miss = fl[1] != 0;
    NITCAD_IO_TRY(pr, reals_of(r, k::predicted, 3, "pass predictions"));
    p.predicted_seconds = pr[0];
    p.predicted_seconds_nodes = pr[1];
    p.predicted_relative = pr[2];
    NITCAD_IO_TRY(al, ints_of(r, k::allocation, 3, "pass allocation"));
    p.candidate_count = static_cast<std::uint64_t>(al[0]);
    p.candidates_taken = static_cast<std::uint64_t>(al[1]);
    p.predicted_dofs = static_cast<std::uint64_t>(al[2]);
    NITCAD_IO_TRY(cands, r.records(k::candidate));
    for (const Reader& cr : cands) {
        results::AdaptCandidateRow c;
        NITCAD_IO_TRY(ci, ints_of(cr, key::c_ints, 7, "candidate"));
        NITCAD_IO_TRY(cf, reals_of(cr, key::c_reals, 9, "candidate"));
        NITCAD_IO_TRY(ind, cr.text(key::c_indicator));
        c.axis = ci[0];
        c.index = static_cast<std::uint64_t>(ci[1]);
        c.region = static_cast<std::uint64_t>(ci[2]);
        c.levels = ci[3];
        c.added_nodes = static_cast<std::uint64_t>(ci[4]);
        c.added_dofs = static_cast<std::uint64_t>(ci[5]);
        c.added_nonzeros = static_cast<std::uint64_t>(ci[6]);
        c.x = cf[0];
        c.y = cf[1];
        c.z = cf[2];
        c.order = cf[3];
        c.benefit = cf[4];
        c.cost = cf[5];
        c.work = cf[6];
        c.realized_ratio = cf[7];
        c.realized = cf[8];
        c.indicator = std::move(ind);
        p.candidates.push_back(std::move(c));
    }
    NITCAD_IO_TRY(orders, r.records(k::order));
    for (const Reader& orr : orders) {
        results::AdaptOrderRow o;
        NITCAD_IO_TRY(ind, orr.text(key::o_indicator));
        NITCAD_IO_TRY(oi, ints_of(orr, key::o_ints, 2, "order"));
        NITCAD_IO_TRY(ov, orr.f64(key::o_order));
        o.indicator = std::move(ind);
        o.region = static_cast<std::uint64_t>(oi[0]);
        o.samples = static_cast<std::uint64_t>(oi[1]);
        o.order = ov;
        p.orders.push_back(std::move(o));
    }
    NITCAD_IO_TRY(cost, reals_of(r, k::cost, 11, "pass cost model"));
    p.cost_a = cost[0];
    p.cost_b = cost[1];
    p.cost_c = cost[2];
    p.cost_d = cost[3];
    p.cost_mesh_per_node = cost[4];
    p.cost_transfer_per_node = cost[5];
    p.cost_residual = cost[6];
    p.cost_alpha = cost[7];
    p.cost_beta = cost[8];
    p.cost_residual_nodes = cost[9];
    p.cost_samples = static_cast<std::uint64_t>(cost[10]);
    NITCAD_IO_TRY(sec, reals_of(r, k::seconds, 4, "pass times"));
    p.seconds_mesh = sec[0];
    p.seconds_transfer = sec[1];
    p.seconds_solve = sec[2];
    p.seconds_estimate = sec[3];
    return p;
}

}  // namespace

std::expected<void, base::Error> save_adapt_report(const std::filesystem::path& path, const results::AdaptReport& a) {
    namespace k = key;
    Writer w;
    w.record(k::run, detail::write_run_record(a.run));
    write_texts(w, k::goal_contact, a.goal_contacts);
    w.f64(k::tolerance, a.tolerance);
    w.i64s(k::points, ints(a.points));
    w.text(k::policy, a.policy);
    w.text(k::weighting, a.weighting);
    w.u64(k::verify, a.verify ? 1 : 0);
    w.u64(k::certify, a.certify_required ? 1 : 0);
    write_texts(w, k::indicator, a.indicators);
    w.text(k::stop, a.stop);
    w.u64(k::rejected_stops, a.rejected_stops);
    w.text(k::fallback, a.fallback);
    w.i64(k::fallback_pass, a.fallback_pass);
    w.u64(k::misses, a.misses);
    w.i64s(k::prediction, std::vector<std::int64_t>{static_cast<std::int64_t>(a.junction_points),
                                                    static_cast<std::int64_t>(a.interface_points),
                                                    static_cast<std::int64_t>(a.predicted_nodes), a.prediction_rounds});
    w.u64(k::coarse_nodes, a.coarse_nodes);
    w.f64(k::coarse_seconds, a.coarse_seconds);
    w.u64(k::verification_nodes, a.verification_nodes);
    w.f64(k::verification_seconds, a.verification_seconds);
    w.f64s(k::verification_goal, a.verification_goal);
    for (const results::AdaptPassRow& p : a.passes) w.record(k::pass, write_pass(p));
    for (const results::AdaptCertificateRow& c : a.certificates) {
        Writer cw;
        cw.i64s(key::g_ints, ints({c.goal, c.point}));
        cw.f64s(key::g_reals, std::vector<double>{c.value, c.estimate, c.verified_error, c.order, c.allowed});
        cw.text(key::g_status, c.status);
        w.record(k::certificate, cw);
    }
    w.text(k::status, a.status);
    w.u64(k::final_mesh, a.final_mesh_digest);
    return detail::write_file(path, Content::adapt_report, w);
}

std::expected<results::AdaptReport, base::Error> load_adapt_report(const std::filesystem::path& path) {
    namespace k = key;
    NITCAD_IO_TRY(bytes, detail::read_file(path, Content::adapt_report));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    results::AdaptReport a;
    NITCAD_IO_TRY(run_record, r.record(k::run));
    NITCAD_IO_TRY(run, detail::read_run_record(run_record));
    a.run = std::move(run);
    NITCAD_IO_TRY(goals, read_texts(r, k::goal_contact));
    a.goal_contacts = std::move(goals);
    NITCAD_IO_TRY(tol, r.f64(k::tolerance));
    a.tolerance = tol;
    NITCAD_IO_TRY(points, r.i64s(k::points));
    a.points = unsigned_of(points);
    NITCAD_IO_TRY(policy, r.text(k::policy));
    a.policy = std::move(policy);
    NITCAD_IO_TRY(weighting, r.text(k::weighting));
    a.weighting = std::move(weighting);
    NITCAD_IO_TRY(verify, r.u64(k::verify));
    a.verify = verify != 0;
    NITCAD_IO_TRY(certify, r.u64(k::certify));
    a.certify_required = certify != 0;
    NITCAD_IO_TRY(indicators, read_texts(r, k::indicator));
    a.indicators = std::move(indicators);
    NITCAD_IO_TRY(stop, r.text(k::stop));
    a.stop = std::move(stop);
    NITCAD_IO_TRY(rejected, r.u64(k::rejected_stops));
    a.rejected_stops = rejected;
    NITCAD_IO_TRY(fallback, r.text(k::fallback));
    a.fallback = std::move(fallback);
    NITCAD_IO_TRY(fallback_pass, r.i64(k::fallback_pass));
    a.fallback_pass = fallback_pass;
    NITCAD_IO_TRY(misses, r.u64(k::misses));
    a.misses = misses;
    NITCAD_IO_TRY(pred, ints_of(r, k::prediction, 4, "prediction"));
    a.junction_points = static_cast<std::uint64_t>(pred[0]);
    a.interface_points = static_cast<std::uint64_t>(pred[1]);
    a.predicted_nodes = static_cast<std::uint64_t>(pred[2]);
    a.prediction_rounds = pred[3];
    NITCAD_IO_TRY(cn, r.u64(k::coarse_nodes));
    a.coarse_nodes = cn;
    NITCAD_IO_TRY(cs, r.f64(k::coarse_seconds));
    a.coarse_seconds = cs;
    NITCAD_IO_TRY(vn, r.u64(k::verification_nodes));
    a.verification_nodes = vn;
    NITCAD_IO_TRY(vs, r.f64(k::verification_seconds));
    a.verification_seconds = vs;
    NITCAD_IO_TRY(vg, r.f64s_or_empty(k::verification_goal));
    const std::size_t goals_by_points = a.goal_contacts.size() * a.points.size();
    if (!vg.empty() && vg.size() != goals_by_points) return std::unexpected(corrupt("adaptation report verification goals"));
    a.verification_goal = std::move(vg);
    NITCAD_IO_TRY(passes, r.records(k::pass));
    for (const Reader& p : passes) {
        NITCAD_IO_TRY(row, read_pass(p, goals_by_points));
        a.passes.push_back(std::move(row));
    }
    NITCAD_IO_TRY(certs, r.records(k::certificate));
    for (const Reader& cr : certs) {
        results::AdaptCertificateRow c;
        NITCAD_IO_TRY(ci, ints_of(cr, key::g_ints, 2, "certificate"));
        NITCAD_IO_TRY(cf, reals_of(cr, key::g_reals, 5, "certificate"));
        NITCAD_IO_TRY(status, cr.text(key::g_status));
        c.goal = static_cast<std::uint64_t>(ci[0]);
        c.point = static_cast<std::uint64_t>(ci[1]);
        if (c.goal >= a.goal_contacts.size()) return std::unexpected(corrupt("adaptation report certificate goal"));
        c.value = cf[0];
        c.estimate = cf[1];
        c.verified_error = cf[2];
        c.order = cf[3];
        c.allowed = cf[4];
        c.status = std::move(status);
        a.certificates.push_back(std::move(c));
    }
    NITCAD_IO_TRY(status, r.text(k::status));
    a.status = std::move(status);
    NITCAD_IO_TRY(final_mesh, r.u64(k::final_mesh));
    a.final_mesh_digest = final_mesh;
    return a;
}

}  // namespace NiTCAD::io
