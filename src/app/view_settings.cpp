// View settings into scenes (view_settings.hpp). Each view is built with the calls and choices of
// `nitcad plot` (src/cli/plot.cpp) for the same options, so the window and the command draw the same
// picture (WG4); the app does not include the command's code (AQ4), and the WG4 test compares the two
// images.
#include "NiTCAD/app/view_settings.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <utility>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::app {

namespace {

using namespace NiTCAD::render;

constexpr std::array<std::pair<ViewKind, std::string_view>, 10> view_names{{{ViewKind::solid, "solid"},
                                                                          {ViewKind::mesh, "mesh"},
                                                                          {ViewKind::field, "field"},
                                                                          {ViewKind::slice, "slice"},
                                                                          {ViewKind::iso, "iso"},
                                                                          {ViewKind::bands, "bands"},
                                                                          {ViewKind::profile, "profile"},
                                                                          {ViewKind::iv, "iv"},
                                                                          {ViewKind::cv, "cv"},
                                                                          {ViewKind::current, "current"}}};

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

std::string current_unit(int dimension) { return dimension == 1 ? "A/cm\xC2\xB2" : dimension == 2 ? "A/cm" : "A"; }

std::size_t contact_index(const device::Device& d, const std::string& name, std::size_t fallback) {
    for (std::size_t k = 0; k < d.contacts().size(); ++k) {
        if (d.contacts()[k].name == name) return k;
    }
    return fallback;
}

}  // namespace

std::string_view to_string(ViewKind v) noexcept {
    for (const auto& [k, name] : view_names) {
        if (k == v) return name;
    }
    return "?";
}

std::optional<ViewKind> find_view(std::string_view name) noexcept {
    for (const auto& [k, n] : view_names) {
        if (n == name) return k;
    }
    return std::nullopt;
}

bool is_plot(ViewKind v) noexcept {
    return v == ViewKind::bands || v == ViewKind::profile || v == ViewKind::iv || v == ViewKind::cv ||
           v == ViewKind::current;
}

std::expected<Scene, base::Error> build_scene(const Workspace& w, const ViewSettings& o) {
    if (!w.device || !w.geometry) return std::unexpected(invalid("nothing is loaded"));
    const device::Device& device = *w.device;
    const Geometry& geometry = *w.geometry;
    const mesh::Mesh& m = device.mesh();
    const auto& points = w.points;

    const auto pick = [&]() -> std::expected<const results::BiasPoint*, base::Error> {
        if (points.empty()) return std::unexpected(invalid("the results hold no field to plot"));
        const std::size_t k = o.point.value_or(points.size() - 1);
        if (k >= points.size()) {
            return std::unexpected(invalid("--point " + std::to_string(k) + " is out of range (" +
                                           std::to_string(points.size()) + " points)"));
        }
        const auto& p = points[k];
        if (p.fields.potential_V.size() != m.node_count()) {
            return std::unexpected(invalid("the result has " + std::to_string(p.fields.potential_V.size()) +
                                           " values per field and the device's mesh " + std::to_string(m.node_count()) +
                                           " nodes: they do not belong together"));
        }
        return &p;
    };

    ViewOptions view;
    view.contacts = o.contacts;
    view.mesh_edges = o.edges;
    if (o.clip) {
        if (m.dimension() != 3) return std::unexpected(invalid("a clipping box needs a 3D device"));
        if (o.kind != ViewKind::solid && o.kind != ViewKind::field) {
            return std::unexpected(invalid("a clipping box applies to the solid and field views"));
        }
        Box b;
        for (std::size_t k = 0; k < 3; ++k) {
            b.lo[k] = o.clip->lo_um[k] * 1e-4;
            b.hi[k] = o.clip->hi_um[k] * 1e-4;
        }
        view.clip = b;
    }
    view.transparent_insulators = !o.opaque;
    view.projection = o.perspective ? Projection::perspective : Projection::orthographic;
    if (o.azimuth) view.azimuth_deg = *o.azimuth;
    if (o.elevation) view.elevation_deg = *o.elevation;

    const auto field = find_field(o.field);
    if (!field) return std::unexpected(invalid("unknown field " + o.field));
    const auto scale_kind = [&]() -> std::optional<ScaleKind> {
        if (o.scale.empty()) return field_info(*field).default_scale;
        if (o.scale == "linear") return ScaleKind::linear;
        if (o.scale == "log") return ScaleKind::log10;
        if (o.scale == "signed-log") return ScaleKind::signed_log10;
        return std::nullopt;
    }();
    if (!scale_kind) return std::unexpected(invalid("unknown scale " + o.scale));
    const auto colormap = find_colormap(o.colormap);
    if (!colormap) return std::unexpected(invalid("unknown colormap " + o.colormap));

    const auto values_for = [&]() -> std::expected<std::vector<double>, base::Error> {
        if (*field == Field::net_doping || *field == Field::donors || *field == Field::acceptors) {
            return field_values(device, results::BiasPoint{}, *field);
        }
        auto p = pick();
        if (!p) return std::unexpected(std::move(p.error()));
        auto v = field_values(device, geometry, **p, *field);
        if (!v) return std::unexpected(std::move(v.error()));
        // A vector magnitude on a log scale: without a range, the ten decades below its largest value.
        if (is_vector_magnitude(*field) && *scale_kind == ScaleKind::log10 && !o.range) {
            double top = 0.0;
            for (const double x : *v) if (std::isfinite(x)) top = std::max(top, x);
            if (top > 0.0) {
                const double floor = std::pow(10.0, std::ceil(std::log10(top)) - 10.0);
                for (double& x : *v) if (std::isfinite(x) && x < floor) x = floor;
            }
        }
        return *v;
    };
    const auto field_options = [&](const std::vector<double>& values) -> std::expected<FieldOptions, base::Error> {
        FieldOptions f;
        f.view = view;
        f.map = *colormap;
        f.contours = o.contours;
        f.title = std::string(field_info(*field).title);
        const bool all_zero = is_vector_magnitude(*field) &&
                              std::ranges::none_of(values, [](double x) { return std::isfinite(x) && x != 0.0; });
        auto s = o.range    ? ColorScale::create(*scale_kind, o.range->first, o.range->second, 1e10)
                 : all_zero ? ColorScale::create(ScaleKind::linear, 0.0, 1.0, 1e10)
                            : ColorScale::fit(*scale_kind, values, 1e10);
        if (!s) return std::unexpected(std::move(s.error()));
        f.scale = *s;
        return f;
    };
    const auto wrap = [](auto&& s) -> std::expected<Scene, base::Error> {
        if (!s) return std::unexpected(std::move(s.error()));
        return Scene(std::move(*s));
    };

    switch (o.kind) {
    case ViewKind::solid:
    case ViewKind::mesh:
        return wrap(o.kind == ViewKind::solid ? solid_view(device, geometry, view) : mesh_view(device, geometry, view));
    case ViewKind::field:
    case ViewKind::slice: {
        auto v = values_for();
        if (!v) return std::unexpected(std::move(v.error()));
        auto f = field_options(*v);
        if (!f) return std::unexpected(std::move(f.error()));
        if (m.dimension() == 1) {
            return wrap(profile_plot(m, *v, std::string(field_info(*field).title),
                                     *scale_kind == ScaleKind::log10 ? AxisScale::log10 : AxisScale::linear));
        }
        if (o.kind == ViewKind::slice) {
            if (o.slices.empty()) return std::unexpected(invalid("--view slice needs --slice"));
            std::vector<Plane> planes;
            for (const auto& s : o.slices) {
                Plane p;
                p.normal[static_cast<std::size_t>(s.axis)] = 1.0;
                p.offset = s.at_um * 1e-4;
                planes.push_back(p);
            }
            return wrap(slice_view(device, geometry, *v, planes, *f));
        }
        return wrap(field_view(device, geometry, *v, *f));
    }
    case ViewKind::iso: {
        if (o.levels.empty()) return std::unexpected(invalid("--view iso needs --levels"));
        if (m.dimension() == 1) return std::unexpected(invalid("isosurfaces need a 2D or 3D device (a 1D field is a profile)"));
        auto v = values_for();
        if (!v) return std::unexpected(std::move(v.error()));
        auto f = field_options(*v);
        if (!f) return std::unexpected(std::move(f.error()));
        return wrap(iso_view(device, geometry, *v, o.levels, *f));
    }
    case ViewKind::profile:
    case ViewKind::bands: {
        auto p = pick();
        if (!p) return std::unexpected(std::move(p.error()));
        if (o.kind == ViewKind::bands) return wrap(band_plot(m, (*p)->bands));
        auto v = values_for();
        if (!v) return std::unexpected(std::move(v.error()));
        return wrap(profile_plot(m, *v, std::string(field_info(*field).title), o.log ? AxisScale::log10 : AxisScale::linear));
    }
    case ViewKind::iv:
    case ViewKind::cv:
    case ViewKind::current: {
        const auto contacts = device.contacts();
        std::vector<SeriesInput> series;
        std::string x_label, y_label;
        if (o.kind == ViewKind::iv) {
            if (!w.sweep) return std::unexpected(invalid("an I-V plot needs a sweep's results"));
            std::vector<std::vector<double>> biases;
            for (const auto& p : w.sweep->points) biases.push_back(p.bias_V);
            const std::size_t swept = contact_index(device, o.contact, swept_contact(biases));
            const std::size_t sense = contact_index(device, o.sense, swept);
            auto c = analysis::current_curve(*w.sweep, swept, sense);
            if (!c) return std::unexpected(std::move(c.error()));
            series.push_back({"I(" + contacts[sense].name + ")", std::move(*c)});
            x_label = "V(" + contacts[swept].name + ") (V)";
            y_label = "I (" + current_unit(m.dimension()) + ")";
        } else if (o.kind == ViewKind::cv) {
            if (!w.small_signal) return std::unexpected(invalid("a C-V plot needs a small-signal run's results"));
            std::vector<std::vector<double>> biases;
            for (const auto& p : w.small_signal->points) biases.push_back(p.dc.bias_V);
            const std::size_t swept = contact_index(device, o.contact, swept_contact(biases));
            const std::size_t sense = contact_index(device, o.sense, swept);
            const std::size_t f = o.frequency.value_or(0);
            if (f >= w.small_signal->frequency_Hz.size()) return std::unexpected(invalid("the frequency is out of range"));
            auto c = analysis::capacitance_curve(*w.small_signal, f, sense, swept, swept);
            if (!c) return std::unexpected(std::move(c.error()));
            series.push_back({"C at " + std::to_string(w.small_signal->frequency_Hz[f]) + " Hz", std::move(*c)});
            x_label = "V(" + contacts[swept].name + ") (V)";
            y_label = m.dimension() == 1 ? "C (F/cm\xC2\xB2)" : m.dimension() == 2 ? "C (F/cm)" : "C (F)";
        } else {
            if (!w.transient) return std::unexpected(invalid("a current plot needs a transient run's results"));
            for (std::size_t k = 0; k < contacts.size(); ++k) {
                if (!o.contact.empty() && contacts[k].name != o.contact) continue;
                std::vector<double> t, i;
                for (const auto& p : w.transient->points) {
                    t.push_back(p.time_s);
                    i.push_back(p.terminal_current[k]);
                }
                auto c = analysis::make_curve(std::move(t), std::move(i));
                if (!c) return std::unexpected(std::move(c.error()));
                series.push_back({contacts[k].name, std::move(*c)});
            }
            x_label = "t (s)";
            y_label = "I (" + current_unit(m.dimension()) + ")";
        }
        return wrap(curve_plot(series, x_label, y_label, AxisScale::linear, o.log ? AxisScale::log10 : AxisScale::linear));
    }
    }
    return std::unexpected(invalid("unknown view"));
}

std::vector<ViewKind> available_views(const Workspace& w) {
    std::vector<ViewKind> out;
    if (!w.device) return out;
    const int dim = w.device->mesh().dimension();
    const bool elements = w.geometry && w.geometry->element_count() > 0;
    if (dim >= 2 && elements) {
        out.push_back(ViewKind::solid);
        out.push_back(ViewKind::mesh);
    }
    out.push_back(dim == 1 ? ViewKind::profile : ViewKind::field);
    if (dim == 3 && elements) out.push_back(ViewKind::slice);
    if (dim >= 2 && elements && w.has_results()) out.push_back(ViewKind::iso);
    if (dim == 1 && w.has_results()) out.push_back(ViewKind::bands);
    if (w.sweep && w.sweep->points.size() >= 2) out.push_back(ViewKind::iv);
    if (w.small_signal && w.small_signal->points.size() >= 2) out.push_back(ViewKind::cv);
    if (w.transient && w.transient->points.size() >= 2) out.push_back(ViewKind::current);
    return out;
}

std::vector<std::string> available_fields(const Workspace& w) {
    std::vector<std::string> out;
    if (!w.device) return out;
    const bool results = w.has_results();
    const results::BiasPoint* last = results ? &w.points.back() : nullptr;
    const bool elements = w.geometry && w.geometry->element_count() > 0;
    for (int k = 0; k <= static_cast<int>(Field::last); ++k) {
        const Field f = static_cast<Field>(k);
        const bool doping = f == Field::net_doping || f == Field::donors || f == Field::acceptors;
        if (!doping && !results) continue;
        if (is_vector_magnitude(f) && !elements) continue;
        if (!doping && last) {
            // The rates only when the result stores them: ask for the values and keep what works.
            if (f >= Field::recombination && f <= Field::btbt_generation &&
                !field_values(*w.device, *last, f).has_value()) {
                continue;
            }
        }
        out.emplace_back(field_info(f).key);
    }
    return out;
}

std::vector<std::string> plot_arguments(const ViewSettings& s) {
    std::vector<std::string> a{"--view", std::string(to_string(s.kind))};
    const bool fieldish = s.kind == ViewKind::field || s.kind == ViewKind::slice || s.kind == ViewKind::iso ||
                          s.kind == ViewKind::profile;
    if (fieldish) a.insert(a.end(), {"--field", s.field});
    if (s.point) a.insert(a.end(), {"--point", std::to_string(*s.point)});
    if (!s.scale.empty() && fieldish) a.insert(a.end(), {"--scale", s.scale});
    if (s.colormap != "turbo") a.insert(a.end(), {"--colormap", s.colormap});
    if (s.range) a.insert(a.end(), {"--range", std::format("{},{}", s.range->first, s.range->second)});
    if (s.contours) a.emplace_back("--contours");
    if (s.edges) a.emplace_back("--edges");
    if (s.opaque) a.emplace_back("--opaque");
    if (s.perspective) a.emplace_back("--perspective");
    if (s.azimuth) a.insert(a.end(), {"--azimuth", std::format("{}", *s.azimuth)});
    if (s.elevation) a.insert(a.end(), {"--elevation", std::format("{}", *s.elevation)});
    for (const auto& p : s.slices) {
        a.insert(a.end(), {"--slice", std::string(1, static_cast<char>('x' + p.axis)) + "=" + std::format("{}", p.at_um)});
    }
    if (s.clip) {
        const auto& c = *s.clip;
        a.insert(a.end(), {"--clip", std::format("{},{},{}:{},{},{}", c.lo_um[0], c.lo_um[1], c.lo_um[2], c.hi_um[0],
                                                 c.hi_um[1], c.hi_um[2])});
    }
    if (!s.levels.empty()) {
        std::string l;
        for (const double v : s.levels) l += (l.empty() ? "" : ",") + std::format("{}", v);
        a.insert(a.end(), {"--levels", l});
    }
    if (!s.contact.empty()) a.insert(a.end(), {"--contact", s.contact});
    if (!s.sense.empty()) a.insert(a.end(), {"--sense", s.sense});
    if (s.frequency) a.insert(a.end(), {"--frequency", std::to_string(*s.frequency)});
    if (s.log) a.emplace_back("--log");
    return a;
}

}  // namespace NiTCAD::app
